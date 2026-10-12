MODULE workproc
  ! This module contains basic work subroutines
  USE matform
  USE misc
  USE linalg
  USE globvars
  USE iso_fortran_env, ONLY: int64
  USE qrlinalg, ONLY: qr_real_state, QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, &
                      QR_ERR_ALLOCATION, QR_ERR_INVALID_STATE, &
                      QR_ERR_DIMENSION_MISMATCH, QR_ERR_NO_CONVERGENCE, &
                      QR_ERR_NONPOSITIVE_OVERLAP

  IMPLICIT NONE

  ! State shared between the inverse-iteration drivers (BasisEnlI, OptCycleI,
  ! FullOpt1I) and the energy routines EnergyIA/EnergyIAM/EnergyIB:
  !   WrkP_LastINVITEnergy   best energy seen so far during FullOpt1I; EnergyIB
  !                          re-anchors the shift on it while WrkP_RefreshShiftInIB
  !                          is set (whole-basis optimization only, nfru==0).
  !   WrkP_WrongStateCount   trial points refused because inverse iteration
  !                          converged on a level other than Glob_WhichEigenvalue
  !                          (IsRequestedEigenstate); not a solver failure.
  REAL(wp) :: WrkP_LastINVITEnergy = 0.0_wp
  LOGICAL  :: WrkP_RefreshShiftInIB = .FALSE.
  INTEGER  :: WrkP_WrongStateCount = 0

  !==================================================================
  ! State of the QR solution method ('Q')
  !==================================================================
  ! The Q method keeps the physical basis in its canonical order. Unlike
  ! the G method it must not move the functions being optimized to the end
  ! of the basis, because the factors owned by qrlinalg describe one
  ! particular row and column order. ActiveFunction maps an optimizer
  ! block to that fixed basis order; ActivePosition is the inverse map.
  ! Physical H, S, nonlinear parameters, raw diagonal overlaps,
  ! derivatives and linear coefficients stay in their Glob_ arrays. The
  ! factors live on rank 0 only (qrlinalg is serial); the physical
  ! matrices and the relationship flags are replicated on every rank.
  ! Each numerical Q entry point constructs or updates this state through
  ! the transaction helpers (AssembleQTrial, ApplyQTrial, AppendQTrial,
  ! TrimQFactors) before requesting an eigenpair (SolveQ).
  !------------------------------------------------------------------
  INTEGER, PARAMETER :: Q_METHOD_SUCCESS = QR_SUCCESS
  INTEGER, PARAMETER :: Q_METHOD_INVALID_ARGUMENT = QR_ERR_INVALID_ARGUMENT
  INTEGER, PARAMETER :: Q_METHOD_ALLOCATION_ERROR = QR_ERR_ALLOCATION
  INTEGER, PARAMETER :: Q_METHOD_INVALID_STATE = QR_ERR_INVALID_STATE
  INTEGER, PARAMETER :: Q_METHOD_DIMENSION_MISMATCH = QR_ERR_DIMENSION_MISMATCH

  TYPE :: QMethodWorkspace
    ! The QR state, meaningful on rank 0 only
    TYPE(qr_real_state) :: Factors
    ! MatrixOrder is the active order represented by the physical matrices,
    ! Capacity the leading dimension allocated for Glob_H and Glob_S,
    ! MaxActive the largest simultaneous optimization block of the step
    INTEGER :: MatrixOrder = 0
    INTEGER :: Capacity = 0
    INTEGER :: MaxActive = 0
    INTEGER :: NumActive = 0
    ! Relationships owned by workproc: qrlinalg cannot know whether a
    ! caller has changed Glob_H or Glob_S since the last update
    LOGICAL :: MatricesAreCanonical = .FALSE.
    LOGICAL :: FactorsMatchMatrices = .FALSE.
    LOGICAL :: MatrixParametersAreStored = .FALSE.
    LOGICAL :: TrialIsReady = .FALSE.
    LOGICAL :: TrialHasDerivatives = .FALSE.
    ! ActiveFunction(a) is the canonical basis index represented by
    ! optimizer block a; ActivePosition(i) is a, or zero when function i
    ! is inactive
    INTEGER, ALLOCATABLE :: ActiveFunction(:)
    INTEGER, ALLOCATABLE :: ActivePosition(:)
    ! A trial is assembled completely before any physical column or QR
    ! factor is changed (full conceptual symmetric columns, though only the
    ! lower triangles of Glob_H and Glob_S are canonical); the previous
    ! columns let a failed multi-column transaction be restored
    REAL(wp), ALLOCATABLE :: PreviousH(:, :), PreviousS(:, :)
    REAL(wp), ALLOCATABLE :: TrialH(:, :), TrialS(:, :)
    REAL(wp), ALLOCATABLE :: PreviousDiagS(:), TrialDiagS(:)
    ! MatrixParam records the nonlinear parameters represented by the
    ! current physical columns (independent of Glob_NonlinParam, where
    ! DRMNG writes its next requested point before the assembly starts);
    ! PreviousParam belongs to PreviousH/S; AcceptedParam is the
    ! optimizer's best point
    REAL(wp), ALLOCATABLE :: MatrixParam(:, :), PreviousParam(:, :), AcceptedParam(:, :)
    REAL(wp) :: AcceptedEnergy = ZERO
    LOGICAL  :: AcceptedPointIsStored = .FALSE.
    ! Diagnostics of the last solve: relative action residuals of the
    ! eigenpair and of the factorization
    REAL(wp) :: LastEigenpairResidual = HUGE(ONE)
    REAL(wp) :: LastFactorResidual = HUGE(ONE)
    INTEGER  :: FreshFactorizations = 0
    ! qrlinalg needs distinct input and output vectors; DeltaH and DeltaS
    ! hold one replacement column and serve the residual checks
    REAL(wp), ALLOCATABLE :: InitialVector(:), SolvedVector(:)
    REAL(wp), ALLOCATABLE :: DeltaH(:), DeltaS(:)
  END TYPE QMethodWorkspace

  TYPE(QMethodWorkspace), SAVE :: Q_Workspace


CONTAINS


  SUBROUTINE ReadIOFile()
    !==================================================================
    ! Subroutine ReadIOFile
    !==================================================================
    ! Reads the input/output data file Glob_DataFileName on rank 0 and
    ! broadcasts everything. The file is also the RESTART file written by
    ! SaveResults, so the record order here and there must agree.
    !
    ! Record layout: BASIS_TYPE (optional), PARTICLES, FIXED_INDEX
    ! (optional), MASSES, CHARGES, REPULSION_SCALING_PARAM[_PLUS|_MINUS]
    ! (optional, up to 3), ATTRACTION_SCALING_PARAM (optional), SYMMETRY,
    ! BASIS_SIZE, CURRENT_ENERGY, WHICH_EIGENVALUE, EIGVAL_TOLERANCE,
    ! INV_IT_PARAM, LAST/BEST/WORST_EIGVAL_TOL, GENERATOR_PARAM, then the
    ! BBOP script, the optimization history and the nonlinear parameters,
    ! each block between separator lines.
    !
    ! Optional records are read into a character buffer and parsed with an
    ! internal READ, so BACKSPACE reliably pushes back a non-matching line.
    ! Every record is validated: EC0100..EC0104 are fatal at once,
    ! EC0210..EC0238 are collected per block and then abort every rank
    ! (MPI_Abort), WC0001..WC0003 are warnings.
    !==================================================================

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: OpenFileErr  ! IOSTAT returned by OPEN
    INTEGER  :: ReadErr      ! IOSTAT of an internal READ
    INTEGER  :: ReadLineErr  ! IOSTAT of a whole-record READ
    INTEGER  :: ReadInt      ! scratch integer field
    REAL(wp) :: ReadRealA    ! scratch real field
    REAL(wp) :: ReadRealB    ! declared, currently unused

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: ReadRealArr  ! currently unused

    ! Packing buffer used to broadcast CHARACTER data as INTEGER.
    ! Sized for the longest item that passes through it.
    INTEGER :: WorkInt(MAX(MAX(Glob_YOperatorStringLength, 20), &
                           Glob_FileNameLength))

    ! Packing buffers for broadcasting the derived-type history array
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: WorkBuffReal
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: WorkBuffInt

    INTEGER :: i, j            ! loop counters
    INTEGER :: j1, j2, j3, j4  ! trimmed file-name lengths
    INTEGER :: Line            ! current line number, for messages

    CHARACTER(70)  :: ReadChar       ! field label of the current record
    CHARACTER(5)   :: ReadBasisType  ! tag from the optional BASIS_TYPE line
    CHARACTER(256) :: ReadLine       ! whole-record buffer

    LOGICAL :: ErrorInDataFile  ! .TRUE. if a fatal error was found
    LOGICAL :: IsBBOPStep       ! loop control for counting BBOP steps

    ! -- data-file validation ---------------------------------------
    INTEGER :: k               ! loop over nonlinear parameters
    INTEGER :: NumOfWarnings   ! non-fatal problems found
    INTEGER :: NumOfValues     ! values counted on a record
    INTEGER :: ValuesExpected  ! values a record should carry
    LOGICAL :: InValue         ! token-scanner state

    ! Whole-record buffer for the nonlinear-parameter lines. Each line
    ! carries 2 integers plus Glob_npt reals written 29 columns wide, so
    ! the length is derived from the compile-time maximum rather than
    ! guessed. NOTE: this is deliberately NOT the 256-byte ReadLine used
    ! for the optional header records, which are short.
    CHARACTER(128+32*Glob_npt_MaxAllowed) :: LongLine


    ErrorInDataFile = .FALSE.
    NumOfWarnings = 0


    !==================================================================
    ! Open the data file (rank 0 only)
    !==================================================================
    ! Only rank 0 touches the file. If it cannot be opened we set the
    ! error flag, broadcast it, and abort on every rank together.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      OPEN(1, FILE=Glob_DataFileName, STATUS='old', IOSTAT=OpenFileErr)

      IF (OpenFileErr /= 0) THEN
        WRITE (*, *) 'Error EC0100 in DataFileInit: data file not found - ', TRIM(ADJUSTL(Glob_DataFileName))
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop


    !==================================================================
    ! Optional BASIS_TYPE line, then the mandatory PARTICLES line
    !==================================================================
    ! BASIS_TYPE may precede PARTICLES. When present it is checked
    ! against Glob_BasisType, so that a data file written for a
    ! different basis is refused instead of being silently misread.
    !
    ! PARTICLES gives the TOTAL particle count N. Glob_n counts
    ! pseudoparticles (Jacobi coordinates), which is one fewer; it is fixed
    ! at compile time, so the record is only checked against it.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) Line = 0

    IF (Glob_ProcID == 0) THEN

      IF (Verbose >= 2) WRITE(*, *) 'Reading initial conditions from data file ', TRIM(ADJUSTL(Glob_DataFileName))

      ! -- optional BASIS_TYPE ---------------------------------------
      READ(1, '(A)', IOSTAT=ReadLineErr) ReadLine
      IF (ReadLineErr == 0) READ(ReadLine, *, IOSTAT=ReadErr) ReadChar(1:10), ReadBasisType

      IF ((ReadLineErr == 0) .AND. (ReadErr == 0) .AND. (ReadChar(1:10) == 'BASIS_TYPE')) THEN

        IF (ReadBasisType /= Glob_BasisType) THEN
          WRITE(*, *) 'Error EC0101 in data file: basis type specifier "', ReadBasisType, &
            '" does not match the basis type of this code "', Glob_BasisType, '"'
          ErrorInDataFile = .TRUE.
        ELSE
          IF (Verbose >= 2) WRITE(*, '(1x,a10,1x,a5)') ReadChar(1:10), ReadBasisType
          Glob_BasisTypeSupplied = .TRUE.
          Line = Line + 1
        ENDIF

      ELSE
        ! Not a BASIS_TYPE line - push the record back for the next READ.
        IF (ReadLineErr == 0) BACKSPACE 1
      ENDIF

      ! -- mandatory PARTICLES ---------------------------------------
      READ(1, *) ReadChar(1:9), ReadInt
      IF (Verbose >= 2) WRITE(*, '(1x,a9,1x,i6)') ReadChar(1:9), ReadInt
      Line = Line + 1

      IF ((ReadInt < 2) .OR. (ReadChar(1:9) /= 'PARTICLES')) THEN
        WRITE(*, *) 'Error EC0102 in data file, line ', Line
        ErrorInDataFile = .TRUE.
      ENDIF

      ! The particle count must match the compiled-in count EXACTLY:
      ! array sizes throughout the code are fixed at compile time, so a
      ! file describing a different system cannot be run through this
      ! executable.
      IF (ReadInt-1 /= Glob_n) THEN
        IF (Verbose >= 2) WRITE (*, *) 'The version of the code you are running was compiled for the case'
        IF (Verbose >= 2) WRITE (*, *) 'when the number of particles in the system is equal to', &
          Glob_AllowedNumOfParticles
        IF (Verbose >= 2) WRITE (*, *) 'while the number of particles specified in the input file is', ReadInt
        WRITE (*, *) 'Please make appropriate changes. Program will now stop.'
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_BasisTypeSupplied, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop


    !==================================================================
    ! Optional FIXED_INDEX line
    !==================================================================
    ! May appear between PARTICLES and MASSES. When present, every basis
    ! function keeps the same z-index (r-premultiplier power) instead of
    ! having one generated and optimized per function.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      READ(1, '(A)', IOSTAT=ReadLineErr) ReadLine
      IF (ReadLineErr == 0) READ(ReadLine, *, IOSTAT=ReadErr) ReadChar(1:11), ReadInt

      IF ((ReadLineErr == 0) .AND. (ReadErr == 0) .AND. (ReadChar(1:11) == 'FIXED_INDEX')) THEN

        Glob_IsIndexFixed = .TRUE.
        Glob_IndexFixedValue = ReadInt

        IF (Verbose >= 2) WRITE(*, '(1x,a11,1x,i6)') ReadChar(1:11), Glob_IndexFixedValue
        Line = Line + 1

        IF ((Glob_IndexFixedValue < 1) .OR. (Glob_IndexFixedValue > Glob_n)) THEN
          WRITE(*, *) 'Error EC0103 in data file, line ', Line
          IF (Verbose >= 2) WRITE(*, *) 'FIXED_INDEX value must be in the range from 1 to', Glob_n
          ErrorInDataFile = .TRUE.
        ENDIF

      ELSE
        ! Not a FIXED_INDEX line - push the record back.
        IF (ReadLineErr == 0) BACKSPACE 1
      ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

    CALL MPI_BCAST(Glob_IsIndexFixed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_IndexFixedValue, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! MASSES
    !==================================================================
    ! Masses of the PARTICLES, not of the pseudoparticles: m1 is the
    ! reference particle and m2..m_{n+1} the rest, so there are
    ! Glob_n+1 values on the line.
    !------------------------------------------------------------------
    ALLOCATE(Glob_Mass(Glob_n+1))

    IF (Glob_ProcID == 0) THEN

      READ(1, *, IOSTAT=ReadErr) ReadChar(1:6), Glob_Mass(1:Glob_n+1)
      Line = Line + 1

      IF (ReadErr /= 0) THEN
        WRITE(*, *) 'Error EC0210 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Could not read the MASSES line. Expected the label MASSES'
        IF (Verbose >= 2) WRITE(*, *) 'followed by', Glob_n+1, 'mass values (one per particle).'
        ErrorInDataFile = .TRUE.
      ELSE

        IF (ReadChar(1:6) /= 'MASSES') THEN
          WRITE(*, *) 'Error EC0211 in data file, line ', Line
          IF (Verbose >= 2) WRITE(*, *) 'Expected the label MASSES but found - ', ReadChar(1:6)
          ErrorInDataFile = .TRUE.
        ENDIF

        ! A zero or negative mass makes the mass matrix singular, which
        ! shows up much later as a meaningless energy rather than an error.
        DO i = 1, Glob_n+1
          IF (Glob_Mass(i) <= ZERO) THEN
            WRITE(*, *) 'Error EC0212 in data file, line ', Line
            IF (Verbose >= 2) WRITE(*, *) 'Mass of particle', i, 'is not positive:', Glob_Mass(i)
            ErrorInDataFile = .TRUE.
          ENDIF
        ENDDO

        IF (Verbose >= 2) WRITE(*, '(1x,a6)', ADVANCE='no') ReadChar(1:6)
        IF (Verbose >= 2) CALL writerealarradv(6, Glob_Mass, Glob_n+1)

      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_Mass, Glob_n+1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! CHARGES
    !==================================================================
    ! q0 is the charge of the reference particle; q1..qn are the
    ! pseudoparticle charges.
    !------------------------------------------------------------------
    ALLOCATE(Glob_PseudoCharge(Glob_n))

    IF (Glob_ProcID == 0) THEN

      READ(1, *, IOSTAT=ReadErr) ReadChar(1:7), Glob_PseudoCharge0, Glob_PseudoCharge(1:Glob_n)
      Line = Line + 1

      IF (ReadErr /= 0) THEN
        WRITE(*, *) 'Error EC0213 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Could not read the CHARGES line. Expected the label CHARGES'
        IF (Verbose >= 2) WRITE(*, *) 'followed by the reference charge and', Glob_n, 'pseudoparticle charges.'
        ErrorInDataFile = .TRUE.
      ELSE

        IF (ReadChar(1:7) /= 'CHARGES') THEN
          WRITE(*, *) 'Error EC0214 in data file, line ', Line
          IF (Verbose >= 2) WRITE(*, *) 'Expected the label CHARGES but found - ', ReadChar(1:7)
          ErrorInDataFile = .TRUE.
        ENDIF

        ! Every charge zero means there is no interaction at all. That is
        ! a legal input but almost always a mistake, so warn rather than stop.
        IF ((Glob_PseudoCharge0 == ZERO) .AND. (ALL(Glob_PseudoCharge(1:Glob_n) == ZERO))) THEN
          IF (Verbose >= 1) WRITE(*, *) 'Warning WC0001 in data file, line ', Line
          IF (Verbose >= 2) WRITE(*, *) 'All particle charges are zero - there is no Coulomb interaction.'
          NumOfWarnings = NumOfWarnings + 1
        ENDIF

        IF (Verbose >= 2) WRITE(*, '(1x,a7)', ADVANCE='no') ReadChar(1:7)
        IF (Verbose >= 2) CALL writereal(6, Glob_PseudoCharge0)
        IF (Verbose >= 2) CALL writerealarradv(6, Glob_PseudoCharge, Glob_n)

      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_PseudoCharge0, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_PseudoCharge, Glob_n, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Optional REPULSION_SCALING_PARAM lines (up to three)
    !==================================================================
    ! REPULSION_SCALING_PARAM, _PLUS (positive charges) and _MINUS (negative
    ! charges) may appear in any order; the 23-character prefix is tested
    ! first and the longer names are told apart by their length. The loop
    ! runs three times; a non-matching record is pushed back each time,
    ! which leaves the file positioned in front of it.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      Glob_RepulsionScalingParam = 1.0_wp
      Glob_RepScalParamSupplied = .FALSE.
      Glob_RepulsionScalingParamPlus = 1.0_wp
      Glob_RepScalParamPlusSupplied = .FALSE.
      Glob_RepulsionScalingParamMinus = 1.0_wp
      Glob_RepScalParamMinusSupplied = .FALSE.

      DO i = 1, 3

        READ(1, '(A)', IOSTAT=ReadLineErr) ReadLine
        IF (ReadLineErr /= 0) EXIT

        READ(ReadLine, *, IOSTAT=ReadErr) ReadChar(1:29), ReadRealA

        IF ((ReadErr /= 0) .OR. (ReadChar(1:23) /= 'REPULSION_SCALING_PARAM')) THEN

          BACKSPACE 1

        ELSE

          IF (ReadChar(1:28) == 'REPULSION_SCALING_PARAM_PLUS') THEN

            Glob_RepulsionScalingParamPlus = ReadRealA
            Glob_RepScalParamPlusSupplied = .TRUE.
            IF (Verbose >= 2) WRITE(*, '(1x,a28)', ADVANCE='no') ReadChar(1:28)
            IF (Verbose >= 2) CALL writerealadv(6, Glob_RepulsionScalingParamPlus)

          ELSEIF (ReadChar(1:29) == 'REPULSION_SCALING_PARAM_MINUS') THEN

            Glob_RepulsionScalingParamMinus = ReadRealA
            Glob_RepScalParamMinusSupplied = .TRUE.
            ! NOTE: a28 against a 29-character argument, so the echoed
            ! label loses its trailing 'S'. Screen output only - the
            ! value that is read and stored is unaffected.
            IF (Verbose >= 2) WRITE(*, '(1x,a28)', ADVANCE='no') ReadChar(1:29)
            IF (Verbose >= 2) CALL writerealadv(6, Glob_RepulsionScalingParamMinus)

          ELSE

            Glob_RepulsionScalingParam = ReadRealA
            Glob_RepScalParamSupplied = .TRUE.
            IF (Verbose >= 2) WRITE(*, '(1x,a23)', ADVANCE='no') ReadChar(1:23)
            IF (Verbose >= 2) CALL writerealadv(6, Glob_RepulsionScalingParam)

          ENDIF

          Line = Line + 1

        ENDIF

      ENDDO

    ENDIF

    CALL MPI_BCAST(Glob_RepScalParamSupplied, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RepulsionScalingParam, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RepScalParamPlusSupplied, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RepulsionScalingParamPlus, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RepScalParamMinusSupplied, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RepulsionScalingParamMinus, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Optional ATTRACTION_SCALING_PARAM line
    !==================================================================
    ! Same buffered-read / BACKSPACE scheme as the repulsion block above.
    ! When the line is absent the parameter falls back to 1.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      READ(1, '(A)', IOSTAT=ReadLineErr) ReadLine
      IF (ReadLineErr == 0) READ(ReadLine, *, IOSTAT=ReadErr) ReadChar(1:24), Glob_AttractionScalingParam

      IF ((ReadLineErr /= 0) .OR. (ReadErr /= 0) .OR. (ReadChar(1:24) /= 'ATTRACTION_SCALING_PARAM')) THEN

        Glob_AttractionScalingParam = 1.0_wp
        Glob_AttrScalParamSupplied = .FALSE.
        IF (ReadLineErr == 0) BACKSPACE 1

      ELSE

        Glob_AttrScalParamSupplied = .TRUE.
        IF (Verbose >= 2) WRITE(*, '(1x,a24)', ADVANCE='no') ReadChar(1:24)
        IF (Verbose >= 2) CALL writerealadv(6, Glob_AttractionScalingParam)
        Line = Line + 1

      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_AttrScalParamSupplied, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_AttractionScalingParam, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! SYMMETRY - the Young operator string
    !==================================================================
    ! The string encodes the Young operator that fixes the permutational
    ! symmetry of the wave function.
    !
    ! CHARACTER data is not broadcast directly here: each character is
    ! converted to its ICHAR code, the integer array is broadcast, and
    ! CHAR converts it back on every rank.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:8), Glob_YOperatorString
      j = LEN_TRIM(Glob_YOperatorString)
      IF (Verbose >= 2) WRITE(*, '(1x,a8)', ADVANCE='no') ReadChar(1:8)
      IF (Verbose >= 2) CALL writestringadv(6, Glob_YOperatorString, j)
      Line = Line + 1

      IF (ReadChar(1:8) /= 'SYMMETRY') THEN
        WRITE(*, *) 'Error EC0215 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label SYMMETRY but found - ', ReadChar(1:8)
        ErrorInDataFile = .TRUE.
      ENDIF

      ! An empty Young operator leaves the symmetry undefined.
      IF (j == 0) THEN
        WRITE(*, *) 'Error EC0216 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'The SYMMETRY line carries an empty Young operator string.'
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    DO i = 1, Glob_YOperatorStringLength
      WorkInt(i) = ICHAR(Glob_YOperatorString(i:i))
    ENDDO

    CALL MPI_BCAST(WorkInt, Glob_YOperatorStringLength, MPI_INTEGER, 0, &
                   MPI_COMM_WORLD, Glob_MPIErrCode)

    DO i = 1, Glob_YOperatorStringLength
      Glob_YOperatorString(i:i) = CHAR(WorkInt(i))
    ENDDO


    !==================================================================
    ! BASIS_SIZE - current number of basis functions
    !==================================================================
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:10), Glob_CurrBasisSize
      IF (Verbose >= 2) WRITE(*, '(1x,a10,1x,i6)') ReadChar(1:10), Glob_CurrBasisSize
      Line = Line + 1

      IF (ReadChar(1:10) /= 'BASIS_SIZE') THEN
        WRITE(*, *) 'Error EC0217 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label BASIS_SIZE but found - ', ReadChar(1:10)
        ErrorInDataFile = .TRUE.
      ENDIF

      IF (Glob_CurrBasisSize < 0) THEN
        WRITE(*, *) 'Error EC0218 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'BASIS_SIZE is negative:', Glob_CurrBasisSize
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_CurrBasisSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! CURRENT_ENERGY - energy of the basis as it stands in the file
    !==================================================================
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:14), Glob_CurrEnergy
      IF (Verbose >= 2) WRITE(*, '(1x,a14)', ADVANCE='no') ReadChar(1:14)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_CurrEnergy)
      Line = Line + 1

      IF (ReadChar(1:14) /= 'CURRENT_ENERGY') THEN
        WRITE(*, *) 'Error EC0219 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label CURRENT_ENERGY but found - ', ReadChar(1:14)
        ErrorInDataFile = .TRUE.
      ENDIF

      ! A bound state has negative energy. A positive value usually means
      ! the file was truncated or the columns are misaligned.
      IF (Glob_CurrEnergy >= ZERO) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Warning WC0002 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'CURRENT_ENERGY is not negative:', Glob_CurrEnergy
        NumOfWarnings = NumOfWarnings + 1
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_CurrEnergy, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! WHICH_EIGENVALUE - which eigenvalue to target
    !==================================================================
    ! Counted from 1 for the lowest. Used when the GSEP solution method
    ! is 'G'.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:16), Glob_WhichEigenvalue
      IF (Verbose >= 2) WRITE(*, '(1x,a16,1x,i6)') ReadChar(1:16), Glob_WhichEigenvalue
      Line = Line + 1

      IF (ReadChar(1:16) /= 'WHICH_EIGENVALUE') THEN
        WRITE(*, *) 'Error EC0220 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label WHICH_EIGENVALUE but found - ', ReadChar(1:16)
        ErrorInDataFile = .TRUE.
      ENDIF

      IF (Glob_WhichEigenvalue < 1) THEN
        WRITE(*, *) 'Error EC0221 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'WHICH_EIGENVALUE must be 1 or greater, found:', Glob_WhichEigenvalue
        ErrorInDataFile = .TRUE.
      ENDIF

      ! A basis of N functions has only N eigenvalues.
      IF ((Glob_CurrBasisSize > 0) .AND. (Glob_WhichEigenvalue > Glob_CurrBasisSize)) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Warning WC0003 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'WHICH_EIGENVALUE is', Glob_WhichEigenvalue, 'but the basis holds only', &
          Glob_CurrBasisSize, 'functions.'
        NumOfWarnings = NumOfWarnings + 1
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_WhichEigenvalue, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! EIGVAL_TOLERANCE - requested eigenvalue accuracy
    !==================================================================
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:16), Glob_EigvalTol
      IF (Verbose >= 2) WRITE(*, '(1x,a16)', ADVANCE='no') ReadChar(1:16)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_EigvalTol)
      Line = Line + 1

      IF (ReadChar(1:16) /= 'EIGVAL_TOLERANCE') THEN
        WRITE(*, *) 'Error EC0222 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label EIGVAL_TOLERANCE but found - ', ReadChar(1:16)
        ErrorInDataFile = .TRUE.
      ENDIF

      IF (Glob_EigvalTol <= ZERO) THEN
        WRITE(*, *) 'Error EC0223 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'EIGVAL_TOLERANCE must be positive, found:', Glob_EigvalTol
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_EigvalTol, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! INV_IT_PARAM - inverse-iteration shift factor
    !==================================================================
    ! The approximate eigenvalue used by the 'I' solver is
    !   Glob_ApproxEnergy = Glob_CurrEnergy * Glob_InvItParameter
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:14), Glob_InvItParameter
      IF (Verbose >= 2) WRITE(*, '(1x,a14)', ADVANCE='no') ReadChar(1:14)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_InvItParameter)
      Line = Line + 1

      IF (ReadChar(1:14) /= 'INVITPARAMETER') THEN
        WRITE(*, *) 'Error EC0224 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label INVITPARAMETER but found - ', ReadChar(1:14)
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_InvItParameter, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Eigenvalue accuracy trackers
    !==================================================================
    !   LAST_EIGVAL_TOL   accuracy achieved at the last solver call
    !   BEST_EIGVAL_TOL   best accuracy seen so far
    !   WORST_EIGVAL_TOL  worst accuracy seen so far
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:15), Glob_LastEigvalTol
      IF (Verbose >= 2) WRITE(*, '(1x,a15)', ADVANCE='no') ReadChar(1:15)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_LastEigvalTol)
      Line = Line + 1

      IF (ReadChar(1:15) /= 'LAST_EIGVAL_TOL') THEN
        WRITE(*, *) 'Error EC0225 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label LAST_EIGVAL_TOL but found - ', ReadChar(1:15)
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_LastEigvalTol, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:15), Glob_BestEigvalTol
      IF (Verbose >= 2) WRITE(*, '(1x,a15)', ADVANCE='no') ReadChar(1:15)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_BestEigvalTol)
      Line = Line + 1

      IF (ReadChar(1:15) /= 'BEST_EIGVAL_TOL') THEN
        WRITE(*, *) 'Error EC0226 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label BEST_EIGVAL_TOL but found - ', ReadChar(1:15)
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_BestEigvalTol, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:16), Glob_WorstEigvalTol
      IF (Verbose >= 2) WRITE(*, '(1x,a16)', ADVANCE='no') ReadChar(1:16)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_WorstEigvalTol)
      Line = Line + 1

      IF (ReadChar(1:16) /= 'WORST_EIGVAL_TOL') THEN
        WRITE(*, *) 'Error EC0227 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label WORST_EIGVAL_TOL but found - ', ReadChar(1:16)
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_WorstEigvalTol, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! GENERATOR_PARAM - random trial-function generator parameters
    !==================================================================
    ! Control how GenerateTrialParam draws candidate basis functions:
    !   Glob_RG_p1  shape of the distribution
    !   Glob_RG_s1  scale parameter 1
    !   Glob_RG_s2  scale parameter 2
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, *) ReadChar(1:15), Glob_RG_p1, Glob_RG_s1, Glob_RG_s2
      IF (Verbose >= 2) WRITE(*, '(1x,a15)', ADVANCE='no') ReadChar(1:15)
      IF (Verbose >= 2) CALL writereal(6, Glob_RG_p1)
      IF (Verbose >= 2) CALL writereal(6, Glob_RG_s1)
      IF (Verbose >= 2) CALL writerealadv(6, Glob_RG_s2)
      Line = Line + 1

      IF (ReadChar(1:15) /= 'GENERATOR_PARAM') THEN
        WRITE(*, *) 'Error EC0228 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'Expected the label GENERATOR_PARAM but found - ', ReadChar(1:15)
        ErrorInDataFile = .TRUE.
      ENDIF

      ! s1 and s2 are the scale parameters of the trial-function
      ! distribution; a non-positive scale makes GenerateTrialParam
      ! produce degenerate candidates for the whole run.
      IF ((Glob_RG_s1 <= ZERO) .OR. (Glob_RG_s2 <= ZERO)) THEN
        WRITE(*, *) 'Error EC0229 in data file, line ', Line
        IF (Verbose >= 2) WRITE(*, *) 'GENERATOR_PARAM scales must be positive, found:', Glob_RG_s1, Glob_RG_s2
        ErrorInDataFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(Glob_RG_p1, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RG_s1, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_RG_s2, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    ! Stop here if any header record above was malformed, rather than
    ! carrying on into the BBOP block with the file out of step.
    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)

    ! Every WC warning above is non-fatal, so report the tally once here
    ! instead of letting single lines scroll past unnoticed.
    IF ((Glob_ProcID == 0) .AND. (NumOfWarnings > 0)) THEN
      IF (Verbose >= 2) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) '*** The data file was accepted with', NumOfWarnings, 'warning(s).'
      IF (Verbose >= 2) WRITE(*, *) '*** Check the WC codes above before trusting the results.'
      IF (Verbose >= 2) WRITE(*, *)
    ENDIF


    !==================================================================
    ! Three separator / heading lines
    !==================================================================
    ! Copied straight through to the screen; they carry no data.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      READ(1, '(a70)') ReadChar(1:70)
      IF (Verbose >= 2) WRITE(*, '(a70)') ReadChar(1:70)
      READ(1, '(a70)') ReadChar(1:70)
      IF (Verbose >= 2) WRITE(*, '(a70)') ReadChar(1:70)
      READ(1, '(a70)') ReadChar(1:70)
      IF (Verbose >= 2) WRITE(*, '(a70)') ReadChar(1:70)
      Line = Line + 3
    ENDIF


    !==================================================================
    ! BBOP, pass 1 of 3: count the steps
    !==================================================================
    ! The BBOP block has unknown length, so it is read three times: pass 1
    ! counts the lines starting with a known action keyword and backspaces
    ! over them, pass 2 reads Action and GSEPSolutionMethod, pass 3 re-reads
    ! each line with the right parameter list. Every keyword listed here
    ! must also have a CASE in pass 3, or the read pointer desynchronizes.
    !------------------------------------------------------------------
    ReadChar(1:70) = ' '

    IF (Glob_ProcID == 0) THEN

      Glob_NumOfBBOPSteps = 0
      IsBBOPStep = .TRUE.

      DO WHILE (IsBBOPStep)

        READ(1, *) ReadChar(1:9)

        IF ((ReadChar(1:9) == 'BASIS_ENL') .OR. (ReadChar(1:9) == 'OPT_CYCLE') .OR. &
            (ReadChar(1:9) == 'FULL_OPT1') .OR. (ReadChar(1:9) == 'EXPC_VALS') .OR. &
            (ReadChar(1:9) == 'OVERLAP_D') .OR. (ReadChar(1:9) == 'ELIM_LCFN') .OR. &
            (ReadChar(1:9) == 'ELIM_LND1') .OR. (ReadChar(1:9) == 'SEPR_LND1') .OR. &
            (ReadChar(1:9) == 'SEPR_FLCF') .OR. (ReadChar(1:9) == 'SAVE_FILE') .OR. &
            (ReadChar(1:9) == 'SAVE_HSWF') .OR. (ReadChar(1:9) == 'SAVE_HS_R') .OR. &
            (ReadChar(1:9) == 'DENSITIES')) THEN
          Glob_NumOfBBOPSteps = Glob_NumOfBBOPSteps + 1
        ELSE
          IsBBOPStep = .FALSE.
        ENDIF

      ENDDO

      ! Rewind over the counted lines PLUS the one that ended the loop.
      DO i = 1, Glob_NumOfBBOPSteps+1
        BACKSPACE 1
      ENDDO

    ENDIF

    CALL MPI_BCAST(Glob_NumOfBBOPSteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    ! The counting loop above stops at the first line that is not a recognized
    ! action, so a BBOP block that is empty or whose first line is misspelled
    ! yields zero steps. Without this check the run would allocate a zero-length
    ! Glob_BBOP and carry on reading from the wrong line.
    IF (Glob_NumOfBBOPSteps <= 0) THEN
      IF (Glob_ProcID == 0) WRITE (*, *) 'Error EC0104 in data file: no BBOP steps found'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    ENDIF

    ALLOCATE(Glob_BBOP(Glob_NumOfBBOPSteps))


    !==================================================================
    ! BBOP, pass 2 of 3: read the action name and solution method
    !==================================================================
    IF (Glob_ProcID == 0) THEN

      DO i = 1, Glob_NumOfBBOPSteps
        READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod
      ENDDO

      DO i = 1, Glob_NumOfBBOPSteps
        BACKSPACE 1
      ENDDO

      !================================================================
      ! Read BBOP (third pass: parse full parameter lists)
      !================================================================
      ! Each action type has a different parameter count:
      !   BASIS_ENL:  Method A B C D E Q R
      !   OPT_CYCLE:  Method A B C D E F G Q R H
      !   FULL_OPT1:  Method A B C D Q R E F FileName1
      !   EXPC_VALS:  Method A
      !   DENSITIES:  Method A FileName1..4
      !   OVERLAP_D:  Method A [FileName1]   (FileName1 optional, default overlap.txt)
      !   SAVE_HSWF:  Method A FileName1..4
      !   SAVE_HS_R:  Method A FileName1 FileName2
      !   ELIM_LCFN, ELIM_LND1:  Method A Q FileName1
      !   SEPR_LND1, SEPR_FLCF:  Method A Q R FileName1
      !   SAVE_FILE:  A FileName1  (no Method field)
      !--------------------------------------------------------------
      Glob_IsOptCycleScripted = .FALSE.


      DO i = 1, Glob_NumOfBBOPSteps

        SELECT CASE (Glob_BBOP(i)%Action(1:9))

        CASE ('BASIS_ENL')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, &
            Glob_BBOP(i)%D, Glob_BBOP(i)%E, Glob_BBOP(i)%Q, Glob_BBOP(i)%R
        ! write(*,'(1x,a9,1x,a1,5(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
        !                   Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writerealadv(6,Glob_BBOP(i)%R)

        CASE ('OPT_CYCLE')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, &
            Glob_BBOP(i)%D, Glob_BBOP(i)%E, Glob_BBOP(i)%F, Glob_BBOP(i)%G, &
            Glob_BBOP(i)%Q, Glob_BBOP(i)%R, Glob_BBOP(i)%H
          ! write(*,'(1x,a9,1x,a1,7(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
          !                 Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E, &
          !                 Glob_BBOP(i)%F,Glob_BBOP(i)%G
          ! call writereal(6,Glob_BBOP(i)%Q)
          ! call writereal(6,Glob_BBOP(i)%R)
          ! write(*,'(1x,i6)') Glob_BBOP(i)%H
          Glob_IsOptCycleScripted = .TRUE.

        CASE ('FULL_OPT1')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, &
            Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
            Glob_BBOP(i)%E, Glob_BBOP(i)%F, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,a1,4(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
        !                 Glob_BBOP(i)%C,Glob_BBOP(i)%D
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writereal(6,Glob_BBOP(i)%R)
        ! write(*,'(2(1x,i6),1x)',advance='no') Glob_BBOP(i)%E,Glob_BBOP(i)%F
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('EXPC_VALS')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
        ! write(*,'(1x,a9,1x,a1,1x,i6)') Glob_BBOP(i)%Action(1:9),  &
        !        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A

        ! DENSITIES: the same four file names as in the reference codes
        ! (CF grid, CF output, density grid, density output); only the
        ! nucleus-nucleus correlation function is available here.
        CASE ('DENSITIES')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName3(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3 = LEN_TRIM(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4 = LEN_TRIM(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))

        ! OVERLAP_D: Method A, plus an OPTIONAL 4th field naming the file that
        ! receives the full spectrum (default: Glob_OverlapFileName). The record
        ! is read whole; the mandatory fields go through a list-directed READ
        ! and the file name is copied verbatim as the 4th blank-delimited
        ! token, so a name containing '/' (a path) is not cut short the way a
        ! list-directed READ would cut it.
        CASE ('OVERLAP_D')
          READ(1, '(a)') ReadLine
          READ(ReadLine, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          Glob_BBOP(i)%FileName1 = Glob_OverlapFileName
          NumOfValues = 0
          InValue = .FALSE.
          DO j = 1, LEN(ReadLine)
            IF (ReadLine(j:j) /= ' ') THEN
              IF (.NOT. InValue) THEN
                InValue = .TRUE.
                NumOfValues = NumOfValues+1
                IF (NumOfValues == 4) THEN
                  k = INDEX(ReadLine(j:), ' ')
                  IF (k == 0) k = LEN(ReadLine)-j+2
                  Glob_BBOP(i)%FileName1 = ReadLine(j:j+k-2)
                  EXIT
                ENDIF
              ENDIF
            ELSE
              InValue = .FALSE.
            ENDIF
          ENDDO
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))

        CASE ('ELIM_LCFN')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('ELIM_LND1')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('SEPR_LND1')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,a1,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writereal(6,Glob_BBOP(i)%R)
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('SEPR_FLCF')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
        !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        ! call writereal(6,Glob_BBOP(i)%Q)
        ! call writereal(6,Glob_BBOP(i)%R)
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('SAVE_FILE')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%A, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
        ! write(*,'(1x,a9,1x,i6,1x)',advance='no') Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A
        ! call writestringadv(6,Glob_BBOP(i)%FileName1,j)

        CASE ('SAVE_HSWF')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName3(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3 = LEN_TRIM(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4 = LEN_TRIM(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))

        ! SAVE_HS_R writes out H and S only. There is no eigenvector or wave
        ! function to save, so it takes two file names instead of the four of
        ! SAVE_HSWF.
        CASE ('SAVE_HS_R')
          READ(1, *) Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A, Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength)
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))

        ENDSELECT

      ENDDO

      ! Separator line that closes the BBOP block
      READ(1, '(a70)') ReadChar(1:70)
      ! write(*,'(a70)')  ReadChar(1:70)

    ENDIF


    !==================================================================
    ! Broadcast the BBOP script to every rank
    !==================================================================
    ! Glob_BBOP is a derived type, so it cannot be broadcast in one go.
    ! Each component is sent separately, and the CHARACTER components
    ! (Action, GSEPSolutionMethod, FileName1..4) go through the same
    ! ICHAR / MPI_BCAST / CHAR packing used for the Young operator.
    !------------------------------------------------------------------
    DO i = 1, Glob_NumOfBBOPSteps

      ! -- Action: 9 characters --------------------------------------
      DO j = 1, 9
        WorkInt(j) = ICHAR(Glob_BBOP(i)%Action(j:j))
      ENDDO

      CALL MPI_BCAST(WorkInt, 9, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      DO j = 1, 9
        Glob_BBOP(i)%Action(j:j) = CHAR(WorkInt(j))
      ENDDO

      ! -- GSEPSolutionMethod: a single character --------------------
      j = ICHAR(Glob_BBOP(i)%GSEPSolutionMethod)
      CALL MPI_BCAST(j, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      Glob_BBOP(i)%GSEPSolutionMethod = CHAR(j)

      ! -- Integer parameters A..H -----------------------------------
      CALL MPI_BCAST(Glob_BBOP(i)%A, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%B, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%C, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%D, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%E, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%F, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%G, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%H, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      ! -- Real parameters Q, R --------------------------------------
      CALL MPI_BCAST(Glob_BBOP(i)%Q, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_BBOP(i)%R, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      ! -- FileName1 -------------------------------------------------
      DO j = 1, Glob_FileNameLength
        WorkInt(j) = ICHAR(Glob_BBOP(i)%FileName1(j:j))
      ENDDO

      CALL MPI_BCAST(WorkInt, Glob_FileNameLength, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      DO j = 1, Glob_FileNameLength
        Glob_BBOP(i)%FileName1(j:j) = CHAR(WorkInt(j))
      ENDDO

      ! -- FileName2 -------------------------------------------------
      DO j = 1, Glob_FileNameLength
        WorkInt(j) = ICHAR(Glob_BBOP(i)%FileName2(j:j))
      ENDDO

      CALL MPI_BCAST(WorkInt, Glob_FileNameLength, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      DO j = 1, Glob_FileNameLength
        Glob_BBOP(i)%FileName2(j:j) = CHAR(WorkInt(j))
      ENDDO

      ! -- FileName3 -------------------------------------------------
      DO j = 1, Glob_FileNameLength
        WorkInt(j) = ICHAR(Glob_BBOP(i)%FileName3(j:j))
      ENDDO

      CALL MPI_BCAST(WorkInt, Glob_FileNameLength, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      DO j = 1, Glob_FileNameLength
        Glob_BBOP(i)%FileName3(j:j) = CHAR(WorkInt(j))
      ENDDO

      ! -- FileName4 -------------------------------------------------
      DO j = 1, Glob_FileNameLength
        WorkInt(j) = ICHAR(Glob_BBOP(i)%FileName4(j:j))
      ENDDO

      CALL MPI_BCAST(WorkInt, Glob_FileNameLength, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      DO j = 1, Glob_FileNameLength
        Glob_BBOP(i)%FileName4(j:j) = CHAR(WorkInt(j))
      ENDDO

    ENDDO

    CALL MPI_BCAST(Glob_IsOptCycleScripted, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Allocate the per-basis-function arrays
    !==================================================================
    !   Glob_History     optimization history, one record per function
    !   Glob_FuncNum     function numbering
    !   Glob_PWR         r-premultiplier power (z-index) per function
    !   Glob_NonlinParam Cholesky elements of L_k, shape (npt x N)
    !------------------------------------------------------------------
    ALLOCATE(Glob_History(Glob_CurrBasisSize))
    ALLOCATE(Glob_FuncNum(Glob_CurrBasisSize))
    ALLOCATE(Glob_PWR(Glob_CurrBasisSize))
    ALLOCATE(Glob_NonlinParam(Glob_npt, Glob_CurrBasisSize))


    !==================================================================
    ! Nothing further to read when the basis is empty
    !==================================================================
    ! A zero-size basis is legitimate: it is how a run that builds its
    ! basis from scratch starts. The arrays above are allocated with
    ! size zero and there is no history or parameter block to read.
    !------------------------------------------------------------------
    IF (Glob_CurrBasisSize == 0) THEN
      ! Stop reading if basis size is zero
      CLOSE(1)
      RETURN
    ENDIF


    ALLOCATE(WorkBuffReal(Glob_CurrBasisSize))
    ALLOCATE(WorkBuffInt(Glob_CurrBasisSize))


    !==================================================================
    ! Optimization history
    !==================================================================
    ! One record per basis function:
    !   <index> <Energy> <CyclesDone> <InitFuncAtLastStep>
    !           <NumOfEnergyEvalDuringFullOpt>
    !
    ! The leading index is read into j and discarded - the records are
    ! taken to be in order. Glob_History is a derived type, so each
    ! component is packed into a flat buffer, broadcast, and unpacked.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      DO i = 1, Glob_CurrBasisSize

        READ(1, *, IOSTAT=ReadErr) j, Glob_History(i)%Energy, Glob_History(i)%CyclesDone, &
          Glob_History(i)%InitFuncAtLastStep, &
          Glob_History(i)%NumOfEnergyEvalDuringFullOpt

        IF (ReadErr /= 0) THEN
          WRITE(*, *) 'Error EC0230 in data file: optimization history block'
          WRITE(*, *) 'Failed to read history record', i, 'of', Glob_CurrBasisSize
          IF (Verbose >= 2) WRITE(*, *) 'Each record needs 5 values: index, energy, cycles done,'
          IF (Verbose >= 2) WRITE(*, *) 'initial function at last step, number of energy evaluations.'
          ErrorInDataFile = .TRUE.
          EXIT
        ENDIF

        ! The leading index is informational, but if it does not run 1..N
        ! then a record is missing or the block is out of order.
        IF (j /= i) THEN
          WRITE(*, *) 'Error EC0231 in data file: optimization history block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'carries index', j
          IF (Verbose >= 2) WRITE(*, *) 'The records are out of order or one of them is missing.'
          ErrorInDataFile = .TRUE.
        ENDIF

        ! None of these three counters can be negative in a file this code wrote.
        IF ((Glob_History(i)%CyclesDone < 0) .OR. (Glob_History(i)%InitFuncAtLastStep < 0) .OR. &
            (Glob_History(i)%NumOfEnergyEvalDuringFullOpt < 0)) THEN
          WRITE(*, *) 'Error EC0232 in data file: optimization history block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'has a negative counter:', &
            Glob_History(i)%CyclesDone, Glob_History(i)%InitFuncAtLastStep, &
            Glob_History(i)%NumOfEnergyEvalDuringFullOpt
          ErrorInDataFile = .TRUE.
        ENDIF

      ENDDO

      READ(1, *) ReadChar(1:70)
      IF (Verbose >= 2) WRITE(*, *)

    ENDIF

    ! Abort before the packing loops below, which would otherwise pack
    ! history entries that were never filled in.
    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)

    ! -- Energy ------------------------------------------------------
    DO i = 1, Glob_CurrBasisSize
      WorkBuffReal(i) = Glob_History(i)%Energy
    ENDDO

    CALL MPI_BCAST(WorkBuffReal, Glob_CurrBasisSize, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    DO i = 1, Glob_CurrBasisSize
      Glob_History(i)%Energy = WorkBuffReal(i)
    ENDDO

    ! -- CyclesDone --------------------------------------------------
    DO i = 1, Glob_CurrBasisSize
      WorkBuffInt(i) = Glob_History(i)%CyclesDone
    ENDDO

    CALL MPI_BCAST(WorkBuffInt, Glob_CurrBasisSize, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    DO i = 1, Glob_CurrBasisSize
      Glob_History(i)%CyclesDone = WorkBuffInt(i)
    ENDDO

    ! -- InitFuncAtLastStep ------------------------------------------
    DO i = 1, Glob_CurrBasisSize
      WorkBuffInt(i) = Glob_History(i)%InitFuncAtLastStep
    ENDDO

    CALL MPI_BCAST(WorkBuffInt, Glob_CurrBasisSize, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    DO i = 1, Glob_CurrBasisSize
      Glob_History(i)%InitFuncAtLastStep = WorkBuffInt(i)
    ENDDO

    ! -- NumOfEnergyEvalDuringFullOpt --------------------------------
    DO i = 1, Glob_CurrBasisSize
      WorkBuffInt(i) = Glob_History(i)%NumOfEnergyEvalDuringFullOpt
    ENDDO

    CALL MPI_BCAST(WorkBuffInt, Glob_CurrBasisSize, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    DO i = 1, Glob_CurrBasisSize
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = WorkBuffInt(i)
    ENDDO


    DEALLOCATE(WorkBuffReal)
    DEALLOCATE(WorkBuffInt)


    !==================================================================
    ! Nonlinear parameters of the basis functions
    !==================================================================
    ! One record per basis function:
    !   <FuncNum> <PWR> <param_1> ... <param_npt>
    !
    ! Read list-directed, so no format string has to change when
    ! Glob_npt changes. The leading function number is read into j and
    ! discarded; the numbering is regenerated below.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      ValuesExpected = 2 + Glob_npt

      DO i = 1, Glob_CurrBasisSize

        ! Read the whole record first so that the values on it can be
        ! counted. A list-directed READ straight from unit 1 silently
        ! ignores surplus values, and quietly runs on into the NEXT record
        ! when values are missing, so a file written for a different
        ! Glob_npt would be accepted and misread without a word.
        ! This assumes one record per basis function, which is how
        ! SaveResults writes them.
        READ(1, '(A)', IOSTAT=ReadLineErr) LongLine

        IF (ReadLineErr /= 0) THEN
          WRITE(*, *) 'Error EC0233 in data file: nonlinear parameter block'
          IF (Verbose >= 2) WRITE(*, *) 'The file ends after', i-1, 'of', Glob_CurrBasisSize, 'records.'
          ErrorInDataFile = .TRUE.
          EXIT
        ENDIF

        ! Count the blank-, tab- or comma-separated values on the record.
        NumOfValues = 0
        InValue = .FALSE.

        DO k = 1, LEN(LongLine)
          IF ((LongLine(k:k) == ' ') .OR. (LongLine(k:k) == ',') .OR. &
              (LongLine(k:k) == CHAR(9))) THEN
            InValue = .FALSE.
          ELSE
            IF (.NOT. InValue) NumOfValues = NumOfValues + 1
            InValue = .TRUE.
          ENDIF
        ENDDO

        IF (NumOfValues /= ValuesExpected) THEN
          WRITE(*, *) 'Error EC0234 in data file: nonlinear parameter block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'carries', NumOfValues, 'values but', &
            ValuesExpected, 'are expected:'
          IF (Verbose >= 2) WRITE(*, *) 'function number, z-index, and', Glob_npt, 'nonlinear parameters.'
          IF (Verbose >= 2) WRITE(*, *) 'This file was most likely written for a different particle number.'
          ErrorInDataFile = .TRUE.
          EXIT
        ENDIF

        READ(LongLine, *, IOSTAT=ReadErr) j, Glob_PWR(i), Glob_NonlinParam(1:Glob_npt, i)

        IF (ReadErr /= 0) THEN
          WRITE(*, *) 'Error EC0235 in data file: nonlinear parameter block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'holds a field that is not a number.'
          ErrorInDataFile = .TRUE.
          EXIT
        ENDIF

        IF (j /= i) THEN
          WRITE(*, *) 'Error EC0236 in data file: nonlinear parameter block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'carries function number', j
          IF (Verbose >= 2) WRITE(*, *) 'The records are out of order or one of them is missing.'
          ErrorInDataFile = .TRUE.
        ENDIF

        ! The z-index indexes the precomputed gamma tables in globvars.f90
        ! (0-based: a power of 0 is a plain Gaussian).
        ! Out of range it reads past the end of those tables, which surfaces
        ! far away inside MatrixElementsOpt instead of here.
        IF ((Glob_PWR(i) < 0) .OR. (Glob_PWR(i) > Glob_MaxPowerAllowed)) THEN
          WRITE(*, *) 'Error EC0237 in data file: nonlinear parameter block'
          IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'has z-index', Glob_PWR(i), &
            'outside the allowed range 0 ..', Glob_MaxPowerAllowed
          ErrorInDataFile = .TRUE.
        ENDIF

        ! A NaN is the one value that compares unequal to itself. It reaches
        ! a data file when a diverged run saves its state, and would poison
        ! every matrix element computed from this basis function.
        DO k = 1, Glob_npt
          IF (Glob_NonlinParam(k, i) /= Glob_NonlinParam(k, i)) THEN
            WRITE(*, *) 'Error EC0238 in data file: nonlinear parameter block'
            IF (Verbose >= 2) WRITE(*, *) 'Record', i, 'parameter', k, 'is not a number (NaN).'
            ErrorInDataFile = .TRUE.
          ENDIF
        ENDDO

      ENDDO

    ENDIF

    CALL MPI_BCAST(ErrorInDataFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (ErrorInDataFile) CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)

    CALL MPI_BCAST(Glob_NonlinParam, Glob_npt*Glob_CurrBasisSize, &
                   MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_PWR, Glob_CurrBasisSize, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Function numbering
    !==================================================================
    ! Numbered sequentially in the order read, so the numbers stored in
    ! the file are informational only.
    !------------------------------------------------------------------
    DO i = 1, Glob_CurrBasisSize
      Glob_FuncNum(i) = i
    ENDDO


    CLOSE(1)

  END SUBROUTINE ReadIOFile


  SUBROUTINE SaveResults(FileName, Sort)
    !==================================================================
    ! Subroutine SaveResults
    !==================================================================
    ! Writes the current state in the layout ReadIOFile reads back (the
    ! file is the RESTART file; keep the record order of the two routines
    ! in step). Rank 0 writes, every other rank returns at once.
    !
    ! Arguments (both optional):
    !   FileName - target file; default Glob_DataFileName.
    !   Sort     - 'yes'/'YES' writes the functions in order of their
    !              function numbers and needs Glob_IntWorkArrForSaveResults
    !              (>= Glob_CurrBasisSize entries) allocated by the caller.
    !
    ! Optional records (BASIS_TYPE, FIXED_INDEX, the scaling parameters) are
    ! written only when the file that was read carried them. Guards: the
    ! CURRENT_ENERGY failure sentinel is never written, Glob_History(0) is
    ! not indexed for an empty basis, and an unknown BBOP action stops the
    ! run (CASE DEFAULT) instead of being dropped silently.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(*) :: FileName  ! target file, defaults to Glob_DataFileName
    CHARACTER(*) :: Sort      ! 'yes'/'YES' to sort by function number
    OPTIONAL :: FileName, Sort

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j            ! loop counter / trimmed name length
    INTEGER :: j1, j2, j3, j4  ! trimmed file-name lengths
    LOGICAL :: SortNeeded      ! .TRUE. when Sort requested sorting


    IF (Glob_ProcID == 0) THEN


      !==================================================================
      ! Open the output file
      !==================================================================
      ! STATUS='replace' because this file IS the restart file - it is
      ! rewritten in place on every save.
      !------------------------------------------------------------------
      IF (PRESENT(FileName)) THEN
        OPEN(1, FILE=FileName, STATUS='replace')
      ELSE
        OPEN(1, FILE=Glob_DataFileName, STATUS='replace')
      ENDIF


      !==================================================================
      ! Physical system: BASIS_TYPE, PARTICLES, FIXED_INDEX
      !==================================================================
      ! The two optional records are written only if the file that was
      ! read carried them, so that a round trip does not invent records.
      !------------------------------------------------------------------
      IF (Glob_BasisTypeSupplied) WRITE(1, '(1x,a10,1x,a5)') 'BASIS_TYPE', Glob_BasisType

      WRITE(1, '(1x,a9,1x,i6)') 'PARTICLES', Glob_n+1

      IF (Glob_IsIndexFixed) WRITE(1, '(1x,a11,1x,i6)') 'FIXED_INDEX', Glob_IndexFixedValue


      !==================================================================
      ! MASSES and CHARGES
      !==================================================================
      ! MASSES carries one value per PARTICLE (Glob_n+1 of them).
      ! CHARGES carries the reference charge q0 followed by the Glob_n
      ! pseudoparticle charges.
      !------------------------------------------------------------------
      WRITE(1, '(1x,a6)', ADVANCE='no') 'MASSES'
      CALL writerealarradv(1, Glob_Mass, Glob_n+1)

      WRITE(1, '(1x,a7)', ADVANCE='no') 'CHARGES'
      CALL writereal(1, Glob_PseudoCharge0)
      CALL writerealarradv(1, Glob_PseudoCharge, Glob_n)


      !==================================================================
      ! Optional interaction-scaling records
      !==================================================================
      ! Each is written only when the corresponding *Supplied flag was
      ! set while the file was read.
      !------------------------------------------------------------------
      IF (Glob_RepScalParamSupplied) THEN
        WRITE(1, '(1x,a23)', ADVANCE='no') 'REPULSION_SCALING_PARAM'
        CALL writerealadv(1, Glob_RepulsionScalingParam)
      ENDIF

      IF (Glob_RepScalParamPlusSupplied) THEN
        WRITE(1, '(1x,a28)', ADVANCE='no') 'REPULSION_SCALING_PARAM_PLUS'
        CALL writerealadv(1, Glob_RepulsionScalingParamPlus)
      ENDIF

      IF (Glob_RepScalParamMinusSupplied) THEN
        WRITE(1, '(1x,a29)', ADVANCE='no') 'REPULSION_SCALING_PARAM_MINUS'
        CALL writerealadv(1, Glob_RepulsionScalingParamMinus)
      ENDIF

      IF (Glob_AttrScalParamSupplied) THEN
        WRITE(1, '(1x,a24)', ADVANCE='no') 'ATTRACTION_SCALING_PARAM'
        CALL writerealadv(1, Glob_AttractionScalingParam)
      ENDIF


      !==================================================================
      ! SYMMETRY and BASIS_SIZE
      !==================================================================
      i = LEN_TRIM(Glob_YOperatorString)
      WRITE(1, '(1x,a8)', ADVANCE='no') 'SYMMETRY'
      CALL writestringadv(1, Glob_YOperatorString, i)

      WRITE(1, '(1x,a10,1x,i6)') 'BASIS_SIZE', Glob_CurrBasisSize


      !==================================================================
      ! CURRENT_ENERGY sentinel guard
      !==================================================================
      ! Never write the rejection sentinel into the data file. A failed
      ! energy evaluation returns 1e31/1e33 so that the optimizer discards
      ! that trial point; if such a value reaches CURRENT_ENERGY the file
      ! stops being a usable restart, because the next run turns it
      ! straight into the inverse-iteration shift.
      !
      ! Keep the last physical energy instead, and say so loudly.
      !------------------------------------------------------------------
      IF (ABS(Glob_CurrEnergy) > 1.0E10_wp) THEN

        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) '*** WARNING: CURRENT_ENERGY holds the failure sentinel ***'
        WRITE(*, '(a)', ADVANCE='no') ' value = '
        CALL writerealadv(6, Glob_CurrEnergy)
        WRITE(*, *) 'Writing the last physical energy to the file instead, so that'
        WRITE(*, *) 'it remains restartable. The results of this step are suspect.'
        WRITE(*, *)

        ! Scan the history BACKWARDS for the most recent physical energy.
        ! Looking only at the current basis size is not enough: when a step
        ! fails, that entry usually holds the same bad value, and the file
        ! then gets written with it anyway - which is exactly how
        ! CURRENT_ENERGY = 2**39 reached the data file and made every
        ! restart from it die on the first solve.
        Glob_CurrEnergy = ZERO

        DO i = Glob_CurrBasisSize, 1, -1
          IF (ABS(Glob_History(i)%Energy) < 1.0E10_wp) THEN
            IF (Glob_History(i)%Energy /= ZERO) THEN
              Glob_CurrEnergy = Glob_History(i)%Energy
              WRITE(*, *) 'Using the energy recorded at basis size ', i
              EXIT
            ENDIF
          ENDIF
        ENDDO

        IF (Glob_CurrEnergy == ZERO) THEN
          WRITE(*, *) 'No physical energy anywhere in the history. The data file'
          WRITE(*, *) 'cannot be made restartable automatically - CURRENT_ENERGY'
          WRITE(*, *) 'must be set by hand to an eigenvalue of this basis.'
        ENDIF

      ENDIF


      !==================================================================
      ! CURRENT_ENERGY and the eigenvalue-solver records
      !==================================================================
      WRITE(1, '(1x,a14)', ADVANCE='no') 'CURRENT_ENERGY'
      CALL writerealadv(1, Glob_CurrEnergy)

      WRITE(1, '(1x,a16,1x,i6)') 'WHICH_EIGENVALUE', Glob_WhichEigenvalue

      WRITE(1, '(1x,a16)', ADVANCE='no') 'EIGVAL_TOLERANCE'
      CALL writerealadv(1, Glob_EigvalTol)

      WRITE(1, '(1x,a14)', ADVANCE='no') 'INVITPARAMETER'
      CALL writerealadv(1, Glob_InvItParameter)

      WRITE(1, '(1x,a15)', ADVANCE='no') 'LAST_EIGVAL_TOL'
      CALL writerealadv(1, Glob_LastEigvalTol)

      WRITE(1, '(1x,a15)', ADVANCE='no') 'BEST_EIGVAL_TOL'
      CALL writerealadv(1, Glob_BestEigvalTol)

      WRITE(1, '(1x,a16)', ADVANCE='no') 'WORST_EIGVAL_TOL'
      CALL writerealadv(1, Glob_WorstEigvalTol)


      !==================================================================
      ! GENERATOR_PARAM - random trial-function generator parameters
      !==================================================================
      WRITE(1, '(1x,a15)', ADVANCE='no') 'GENERATOR_PARAM'
      CALL writereal(1, Glob_RG_p1)
      CALL writereal(1, Glob_RG_s1)
      CALL writerealadv(1, Glob_RG_s2)


      !==================================================================
      ! Separator, summary line, separator
      !==================================================================
      ! ReadIOFile reads these three records as plain text and discards
      ! them, so the summary line is for the reader, not for the code.
      !
      ! ZERO-BASIS GUARD: with an empty basis Glob_History has size zero,
      ! so Glob_History(Glob_CurrBasisSize) would be Glob_History(0) - an
      ! out of bounds access on an unallocated element. Write a row of
      ! zeros in that case instead.
      !------------------------------------------------------------------
      WRITE(1, *) '=============================='

      IF (Glob_CurrBasisSize > 0) THEN

        i = Glob_CurrBasisSize
        WRITE(1, '(1x,i6)', ADVANCE='no') i
        CALL writereal(1, Glob_History(i)%Energy)
        WRITE(1, '(3(1x,i6))') Glob_History(i)%CyclesDone, &
          Glob_History(i)%InitFuncAtLastStep, Glob_History(i)%NumOfEnergyEvalDuringFullOpt

      ELSE

        WRITE(1, '(1x,i6)', ADVANCE='no') 0
        CALL writereal(1, ZERO)
        WRITE(1, '(3(1x,i6))') 0, 0, 0

      ENDIF

      WRITE(1, *) '=============================='


      !==================================================================
      ! BBOP (Basis Building and Optimization Program)
      !==================================================================
      ! One record per step, with the parameter list of its action; the
      ! layout must match the CASE ReadIOFile uses to read it back. CASE
      ! DEFAULT stops the run on an unknown action rather than dropping the
      ! step from the restart file.
      !------------------------------------------------------------------
      DO i = 1, Glob_NumOfBBOPSteps

        SELECT CASE (Glob_BBOP(i)%Action(1:9))

        CASE ('BASIS_ENL')
          WRITE(1, '(1x,a9,1x,a1,5(1x,i6))', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A, Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C, Glob_BBOP(i)%D, Glob_BBOP(i)%E
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writerealadv(1, Glob_BBOP(i)%R)

        CASE ('OPT_CYCLE')
          WRITE(1, '(1x,a9,1x,a1,7(1x,i6))', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A, Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C, Glob_BBOP(i)%D, Glob_BBOP(i)%E, &
            Glob_BBOP(i)%F, Glob_BBOP(i)%G
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writereal(1, Glob_BBOP(i)%R)
          WRITE(1, '(1x,i6)') Glob_BBOP(i)%H

        CASE ('FULL_OPT1')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,4(1x,i6))', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A, Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C, Glob_BBOP(i)%D
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writereal(1, Glob_BBOP(i)%R)
          WRITE(1, '(2(1x,i6),1x)', ADVANCE='no') Glob_BBOP(i)%E, Glob_BBOP(i)%F
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('EXPC_VALS')
          WRITE(1, '(1x,a9,1x,a1,1x,i6)') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A

        CASE ('DENSITIES')
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3 = LEN_TRIM(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4 = LEN_TRIM(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') &
            Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writestring(1, Glob_BBOP(i)%FileName1, j1)
          CALL writestring(1, Glob_BBOP(i)%FileName2, j2)
          CALL writestring(1, Glob_BBOP(i)%FileName3, j3)
          CALL writestringadv(1, Glob_BBOP(i)%FileName4, j4)

        ! OVERLAP_D carries the spectrum file name as its 4th field, so the
        ! line reads back the way it was given (the default name is written
        ! when the line had none).
        CASE ('OVERLAP_D')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          IF (j == 0) THEN
            Glob_BBOP(i)%FileName1 = Glob_OverlapFileName
            j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          ENDIF
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('ELIM_LCFN')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('ELIM_LND1')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('SEPR_LND1')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,i6)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writereal(1, Glob_BBOP(i)%R)
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('SEPR_FLCF')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writereal(1, Glob_BBOP(i)%Q)
          CALL writereal(1, Glob_BBOP(i)%R)
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        ! SAVE_FILE carries no solution-method field.
        CASE ('SAVE_FILE')
          j = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,i6,1x)', ADVANCE='no') Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%A
          CALL writestringadv(1, Glob_BBOP(i)%FileName1, j)

        CASE ('SAVE_HSWF')
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3 = LEN_TRIM(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4 = LEN_TRIM(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') &
            Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writestring(1, Glob_BBOP(i)%FileName1, j1)
          CALL writestring(1, Glob_BBOP(i)%FileName2, j2)
          CALL writestring(1, Glob_BBOP(i)%FileName3, j3)
          CALL writestringadv(1, Glob_BBOP(i)%FileName4, j4)

        ! SAVE_HS_R writes H and S only, so it carries two file names.
        CASE ('SAVE_HS_R')
          j1 = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2 = LEN_TRIM(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          WRITE(1, '(1x,a9,1x,a1,1x,i6)', ADVANCE='no') &
            Glob_BBOP(i)%Action(1:9), Glob_BBOP(i)%GSEPSolutionMethod, Glob_BBOP(i)%A
          CALL writestring(1, Glob_BBOP(i)%FileName1, j1)
          CALL writestringadv(1, Glob_BBOP(i)%FileName2, j2)

        CASE DEFAULT
          WRITE(*, *) 'Error in SaveResults: unsupported BBOP action: ', &
                     TRIM(Glob_BBOP(i)%Action)
          CLOSE(1)
          STOP

        ENDSELECT

      ENDDO

      WRITE(1, *) '=============================='


      !==================================================================
      ! Optimization history
      !==================================================================
      ! One record per basis function:
      !   <index> <Energy> <CyclesDone> <InitFuncAtLastStep>
      !           <NumOfEnergyEvalDuringFullOpt>
      !------------------------------------------------------------------
      DO i = 1, Glob_CurrBasisSize
        WRITE(1, '(1x,i6)', ADVANCE='no') i
        CALL writereal(1, Glob_History(i)%Energy)
        WRITE(1, '(3(1x,i6))') Glob_History(i)%CyclesDone, &
          Glob_History(i)%InitFuncAtLastStep, Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      ENDDO

      WRITE(1, *) '=============================='


      !==================================================================
      ! Nonlinear parameters of the basis functions
      !==================================================================
      ! One record per basis function:
      !   <FuncNum> <PWR> <param_1> ... <param_npt>
      !
      ! When sorting is requested, Glob_IntWorkArrForSaveResults is used
      ! to invert Glob_FuncNum, so that the functions come out ordered by
      ! function number. That array must already be allocated by the
      ! caller - see the header note.
      !------------------------------------------------------------------
      SortNeeded = .FALSE.

      IF (PRESENT(Sort)) THEN
        IF ((Sort == 'yes') .OR. (Sort == 'YES')) SortNeeded = .TRUE.
      ENDIF

      IF (SortNeeded) THEN

        ! Invert the numbering: work(FuncNum(i)) = i
        DO i = 1, Glob_CurrBasisSize
          Glob_IntWorkArrForSaveResults(Glob_FuncNum(i)) = i
        ENDDO

        DO i = 1, Glob_CurrBasisSize
          WRITE(1, '(1x,i6,1x,i6)', ADVANCE='no') i, Glob_PWR(Glob_IntWorkArrForSaveResults(i))
          CALL writerealarradv(1, Glob_NonlinParam(1:Glob_npt, Glob_IntWorkArrForSaveResults(i)), Glob_npt)
        ENDDO

      ELSE

        DO i = 1, Glob_CurrBasisSize
          WRITE(1, '(1x,i6,1x,i6)', ADVANCE='no') i, Glob_PWR(i)
          CALL writerealarradv(1, Glob_NonlinParam(1:Glob_npt, i), Glob_npt)
        ENDDO

      ENDIF


      CLOSE(1)

    ENDIF

  END SUBROUTINE SaveResults


  SUBROUTINE ReadBlackList()
    !==================================================================
    ! Subroutine ReadBlackList
    !==================================================================
    ! Reads Glob_BlackListFileName: one function number per line (any
    ! order, repetitions harmless) of the functions the cyclic optimization
    ! must NOT touch. Rank 0 reads, the result is broadcast.
    !
    ! Result: Glob_lbf is the largest number in the file (the size of
    ! Glob_Blacklisted) and Glob_Blacklisted(i) is .TRUE. for a listed
    ! function. Callers must test i <= Glob_lbf before indexing. The file
    ! is read twice: first to validate and find the maximum, then, after a
    ! REWIND, to mark the entries. A missing or empty file gives Glob_lbf=0
    ! and a quiet return (Glob_Blacklisted stays unallocated); a bad value
    ! is fatal (EC0110, MPI_Abort).
    !==================================================================

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: OpenFileErr  ! IOSTAT returned by OPEN
    INTEGER :: ReadErr      ! IOSTAT of the record READ
    INTEGER :: i            ! function number / loop counter

    LOGICAL :: ErrorInFile  ! file missing, or a bad value in it


    ErrorInFile = .FALSE.


    !==================================================================
    ! Open the black list file (rank 0 only)
    !==================================================================
    ! A missing file is NOT an error - it simply means nothing is
    ! blacklisted. ErrorInFile is reused here as "there is no usable
    ! file", and the RETURN below is the quiet exit, not a failure.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      OPEN(1, FILE=Glob_BlackListFileName, STATUS='old', IOSTAT=OpenFileErr)

      IF (OpenFileErr /= 0) THEN
        WRITE(*, *)
        IF (Verbose >= 2) WRITE(*, *) 'Black list file not found: ', TRIM(Glob_BlackListFileName)
        WRITE(*, *) 'No basis functions will be blacklisted during cyclic optimization'
        WRITE(*, *)
        ErrorInFile = .TRUE.
      ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    ! Every rank must take this exit together, hence the broadcast above.
    IF (ErrorInFile) THEN
      Glob_lbf = 0
      RETURN
    ENDIF


    !==================================================================
    ! First pass: validate the entries and find the largest one
    !==================================================================
    ! Glob_lbf ends up holding the largest function number in the file,
    ! which is the size Glob_Blacklisted has to be. A number outside
    ! 1..Glob_CurrBasisSize does not name a basis function that exists,
    ! so it is fatal rather than ignorable.
    !
    ! The loop ends on the first read failure, which is normally just
    ! end of file.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      Glob_lbf = 0
      ReadErr = 0

      DO WHILE ((ReadErr == 0) .AND. (.NOT. ErrorInFile))

        READ(1, *, IOSTAT=ReadErr) i

        IF (ReadErr == 0) THEN
          IF ((i <= Glob_CurrBasisSize) .AND. (i > 0)) THEN
            Glob_lbf = MAX(i, Glob_lbf)
          ELSE
            ErrorInFile = .TRUE.
          ENDIF
        ENDIF

      ENDDO

    ENDIF

    CALL MPI_BCAST(ErrorInFile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (ErrorInFile) THEN

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0110 in ReadBlackList: incorrect values in file ', TRIM(Glob_BlackListFileName)
        CLOSE(1)
      ENDIF

      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

    ENDIF


    !==================================================================
    ! An empty file blacklists nothing
    !==================================================================
    ! Glob_lbf is still 0, so there is nothing to allocate and nothing
    ! to broadcast. Note that Glob_Blacklisted stays UNALLOCATED on this
    ! path - see the caller contract in the header.
    !------------------------------------------------------------------
    CALL MPI_BCAST(Glob_lbf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (Glob_lbf == 0) THEN

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 2) WRITE(*, *) 'Black list file is empty - ', TRIM(Glob_BlackListFileName)
        WRITE(*, *) 'no basis functions will be blacklisted during cyclic optimization'
        WRITE(*, *)
      ENDIF

      RETURN

    ENDIF


    !==================================================================
    ! Second pass: mark the blacklisted functions
    !==================================================================
    ! The array is allocated on EVERY rank, filled on rank 0, and then
    ! broadcast. REWIND puts the file back to the first record before
    ! reading it again.
    !
    ! Indexing Glob_Blacklisted(i) is safe here without a further test:
    ! the first pass already rejected anything outside 1..Glob_lbf.
    !------------------------------------------------------------------
    ALLOCATE(Glob_Blacklisted(Glob_lbf))

    IF (Glob_ProcID == 0) THEN

      Glob_Blacklisted(1:Glob_lbf) = .FALSE.

      REWIND(1)
      ReadErr = 0

      DO WHILE (ReadErr == 0)
        READ(1, *, IOSTAT=ReadErr) i
        IF (ReadErr == 0) Glob_Blacklisted(i) = .TRUE.
      ENDDO

      CLOSE(1)

      ! Echo the list that was read, so the run log records exactly
      ! which functions were frozen.
      WRITE(*, *)
      IF (Verbose >= 2) WRITE(*, *) 'Black list file has been read - ', TRIM(Glob_BlackListFileName)
      WRITE(*, *) 'The following basis functions will be blacklisted during cyclic optimization:'

      DO i = 1, Glob_lbf
        IF (Glob_Blacklisted(i)) WRITE(*, '(1x,i6)', ADVANCE='no') i
      ENDDO

      WRITE(*, *)
      WRITE(*, *)

    ENDIF

    CALL MPI_BCAST(Glob_Blacklisted, Glob_lbf, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


  END SUBROUTINE ReadBlackList


  SUBROUTINE ProgramDataInit()
    !==================================================================
    ! Subroutine ProgramDataInit
    !==================================================================
    ! Builds, once, every quantity that depends only on the physical system
    ! and the symmetry; call it right after ReadIOFile.
    !
    ! In order: Glob_AbsTolForDSYGVX; the charge products
    ! Glob_PseudoChargeMatrix (bare), Glob_ScaledPseudoChargeMatrix (scaled)
    ! and Glob_ChargeMatrix (the 1..n form matelem reads); Glob_MassMatrix;
    ! Glob_bvc; the pair transpositions Glob_Transposit; the Young operator
    ! Y and Y^{+}Y (Glob_YMatr/Glob_YCoeff/Glob_NumYTerms and the YHY
    ! counterparts) parsed from Glob_YOperatorString, e.g. (1+P34)(1-P12);
    ! the identical-particle sets and equivalent pairs (Glob_IdentPartList,
    ! Glob_EqvPairList and their counts).
    !
    ! A product of pair-permutation OPERATORS corresponds to the REVERSED
    ! product of the matrices acting on the nonlinear parameters, so the
    ! factor loops run backwards. Errors: EC0115 (illegal symbol in the
    ! Young operator), EC0116 (unbalanced brackets); both abort all ranks.
    !==================================================================

    IMPLICIT NONE

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
    INTEGER      :: n, npart, indexh
    INTEGER      :: i, j, k, p, q, t, s, w, ii, jj, kk
    CHARACTER(1) :: c1, cc1
    INTEGER      :: StrLen, NumFactY
    INTEGER      :: TotNumOfYTerms, TotNumOfYHYTerms, CurrNumOfTerms
    INTEGER      :: L, R, FirstLPos, LastRPos, MaxNumTermsInFact
    INTEGER      :: Coeff, Cf3
    INTEGER      :: pi, pj, pt, ps
    LOGICAL      :: AreTermsIdentical

    ! Charge products for the scaling block below: qi and qj are the two
    ! charges, qq their bare product, and sqq that product after the
    ! repulsion / attraction scaling parameters have been applied.
    REAL(wp) :: qi, qj, qq, sqq

    ! Lightest and heaviest pseudoparticle masses, used to decide whether
    ! all pseudoparticle masses are the same.
    REAL(wp) :: ml, mh

    INTEGER, ALLOCATABLE, DIMENSION(:)       :: TempSymCoeff, TempSymCoeff1
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: NumTermsInYOpFact
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: TempSymMatr, TempSymMatr1
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: Matr1, Matr2, Matr3, Matr4
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: IdentParticleSet
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: IdentPseudoPartPairSet

    CHARACTER(Glob_YOperatorStringLength), ALLOCATABLE, DIMENSION(:) :: YOpStr, YHOpStr


    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Initializing program data'


    !==================================================================
    ! Machine-dependent constants
    !==================================================================
    ! DLAMCH('S') is the LAPACK safe-minimum. Twice that is the absolute
    ! tolerance handed to DSYGVX.
    !------------------------------------------------------------------
    Glob_AbsTolForDSYGVX = 2*DLAMCH('S')

    n = Glob_n
    npart = n+1


    !==================================================================
    ! Charge-product matrices, bare and scaled
    !==================================================================
    ! Indexed 0..n: index 0 is the reference particle. Glob_PseudoChargeMatrix
    ! holds qi*qj; Glob_ScaledPseudoChargeMatrix holds the same product times
    ! Glob_AttractionScalingParam (qq < 0), or Glob_RepulsionScalingParam
    ! times _Plus (both charges positive) or _Minus (both negative). Without
    ! scaling records every parameter is 1 and the matrices are identical.
    !------------------------------------------------------------------
    ALLOCATE(Glob_PseudoChargeMatrix(0:n, 0:n))
    ALLOCATE(Glob_ScaledPseudoChargeMatrix(0:n, 0:n))

    Glob_PseudoChargeMatrix(0:n, 0:n) = ZERO
    Glob_ScaledPseudoChargeMatrix(0:n, 0:n) = ZERO

    DO j = 0, n

      IF (j == 0) THEN
        qj = Glob_PseudoCharge0
      ELSE
        qj = Glob_PseudoCharge(j)
      ENDIF

      DO i = 0, n

        IF (i == 0) THEN
          qi = Glob_PseudoCharge0
        ELSE
          qi = Glob_PseudoCharge(i)
        ENDIF

        qq = qi*qj

        IF (qq < 0.0_wp) THEN
          sqq = qq*Glob_AttractionScalingParam
        ELSE
          IF ((qi > 0.0_wp) .AND. (qj > 0.0_wp)) THEN
            sqq = qq*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamPlus
          ELSE
            sqq = qq*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamMinus
          ENDIF
        ENDIF

        Glob_PseudoChargeMatrix(i, j) = qq
        Glob_ScaledPseudoChargeMatrix(i, j) = sqq

      ENDDO

    ENDDO


    !==================================================================
    ! Charge-product matrix Glob_ChargeMatrix, the form matelem reads
    !==================================================================
    ! Indexed 1..n: the diagonal holds q_i*q_0 and the strict lower
    ! triangle (j,i), j>i, holds q_i*q_j; the upper triangle is not read.
    ! Filled from Glob_ScaledPseudoChargeMatrix so the optional scaling
    ! records act on it; without them it holds the bare charge products.
    !------------------------------------------------------------------
    ALLOCATE(Glob_ChargeMatrix(n, n))

    Glob_ChargeMatrix(1:n, 1:n) = ZERO

    DO i = 1, n
      Glob_ChargeMatrix(i, i) = Glob_ScaledPseudoChargeMatrix(i, 0)
      DO j = i+1, n
        Glob_ChargeMatrix(j, i) = Glob_ScaledPseudoChargeMatrix(i, j)
      ENDDO
    ENDDO


    !==================================================================
    ! Mass matrix M
    !==================================================================
    ! Every element carries 1/(2*m1) from the reference particle; the
    ! diagonal picks up 1/(2*m_{i+1}) in addition.
    !------------------------------------------------------------------
    ALLOCATE(Glob_MassMatrix(n, n))

    Glob_MassMatrix(1:n, 1:n) = ONEHALF/Glob_Mass(1)

    DO i = 1, n
      Glob_MassMatrix(i, i) = Glob_MassMatrix(i, i)+ONEHALF/Glob_Mass(i+1)
    ENDDO


    !==================================================================
    ! Density-vector coefficients Glob_bvc
    !==================================================================
    ! Used when evaluating particle densities. Column i holds the
    ! coefficients for particle i.
    !------------------------------------------------------------------
    ALLOCATE(Glob_bvc(n, npart))

    Glob_MassTotal = SUM(Glob_Mass(1:npart))

    DO i = 1, npart
      Glob_bvc(1:n, i) = -Glob_Mass(2:n+1)/Glob_MassTotal
    ENDDO

    DO i = 2, npart
      Glob_bvc(i-1, i) = Glob_bvc(i-1, i)+ONE
    ENDDO


    !==================================================================
    ! Lightest and heaviest pseudoparticle
    !==================================================================
    ! The reference particle is excluded from the search. ml and mh end
    ! up holding the lightest and heaviest pseudoparticle mass, indexh
    ! the index of the heaviest.
    !
    ! When the two differ the pseudoparticle masses are not all equal,
    ! which switches off the simplifications used elsewhere in the code.
    !------------------------------------------------------------------
    indexh = 0
    mh = 0
    ml = 2*Glob_MassTotal

    DO i = 1, n

      IF (Glob_Mass(i+1) < ml) THEN
        ml = Glob_Mass(i+1)
      ENDIF

      IF (Glob_Mass(i+1) > mh) THEN
        mh = Glob_Mass(i+1)
        indexh = i
      ENDIF

    ENDDO

    IF (ABS(ml - mh) > 1.d-14) THEN
      Glob_ArePseudoParticleMassesTheSame = .FALSE.
    ENDIF


    !==================================================================
    ! Pair-transposition matrices Pij
    !==================================================================
    ! Glob_Transposit(1:n,1:n,i,j) is the matrix of the transposition that
    ! swaps particles i and j, acting on the matrix of nonlinear parameters.
    ! P1i (i/=1) is the unit matrix with column i-1 replaced by -1 in every
    ! row; Pij (i,j/=1) is the unit matrix with rows i-1 and j-1 exchanged.
    !------------------------------------------------------------------
    ALLOCATE(Glob_Transposit(n, n, npart, npart))

    ! First set all of them to be unit matrices
    Glob_Transposit(1:n, 1:n, 1:npart, 1:npart) = 0

    DO i = 1, npart
      DO j = 1, npart
        DO k = 1, n
          Glob_Transposit(k, k, i, j) = 1
        ENDDO
      ENDDO
    ENDDO

    ! Now continue depending on type of transposition (P1i or Pij)
    DO i = 2, npart
      Glob_Transposit(1:n, i-1, 1, i) = -1
    ENDDO

    DO i = 2, npart
      DO j = i+1, npart
        Glob_Transposit(i-1, i-1, i, j) = 0
        Glob_Transposit(j-1, j-1, i, j) = 0
        Glob_Transposit(j-1, i-1, i, j) = 1
        Glob_Transposit(i-1, j-1, i, j) = 1
        Glob_Transposit(i-1, i-1, j, i) = 0
        Glob_Transposit(j-1, j-1, j, i) = 0
        Glob_Transposit(j-1, i-1, j, i) = 1
        Glob_Transposit(i-1, j-1, j, i) = 1
      ENDDO
    ENDDO


    !==================================================================
    ! Young operator, stage 1: strip blanks and multiplication signs
    !==================================================================
    ! Every removal shifts the tail of the string left by one and blanks
    ! the vacated position, so the string stays contiguous.
    !------------------------------------------------------------------
    StrLen = LEN_TRIM(Glob_YOperatorString)

    DO i = 1, StrLen
      c1 = Glob_YOperatorString(i:i)
      IF ((c1 == ' ') .OR. (c1 == '*')) THEN
        DO j = i, StrLen-1
          Glob_YOperatorString(j:j) = Glob_YOperatorString(j+1:j+1)
        ENDDO
        Glob_YOperatorString(j:j) = ' '
      ENDIF
    ENDDO

    StrLen = LEN_TRIM(Glob_YOperatorString)


    !==================================================================
    ! Young operator, stage 2: reject illegal symbols
    !==================================================================
    ! Only digits, 'P', the signs, '*' and the two brackets are legal.
    ! Anything else means the expression cannot be parsed, so stop here
    ! rather than build a wrong operator.
    !------------------------------------------------------------------
    DO i = 1, StrLen

      c1 = Glob_YOperatorString(i:i)

      IF ((c1 /= '1') .AND. (c1 /= '2') .AND. (c1 /= '3') .AND. (c1 /= '4') .AND. (c1 /= '5') .AND. &
          (c1 /= '6') .AND. (c1 /= '7') .AND. (c1 /= '8') .AND. (c1 /= '9') .AND. (c1 /= '0') .AND. &
          (c1 /= 'P') .AND. (c1 /= '+') .AND. (c1 /= '-') .AND. (c1 /= '*') .AND. (c1 /= ')') .AND. &
          (c1 /= '(')) THEN
        WRITE(*, *) 'Error EC0115 in ProgramDataInit: the Young operator expression'
        WRITE(*, *) 'contains wrong symbols'
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF

    ENDDO


    !==================================================================
    ! Young operator, stage 3: check that the brackets balance
    !==================================================================
    ! R also ends up holding the number of bracketed factors, which the
    ! next stage needs.
    !------------------------------------------------------------------
    L = 0
    R = 0

    DO i = 1, StrLen
      IF (Glob_YOperatorString(i:i) == ')') R = R+1
      IF (Glob_YOperatorString(i:i) == '(') L = L+1
    ENDDO

    IF (R /= L) THEN
      WRITE(*, *) 'Error EC0116 in ProgramDataInit: the numer of left and right brackets in the'
      WRITE(*, *) 'Young operator is different'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF


    !==================================================================
    ! Young operator, stage 4: count the factors
    !==================================================================
    ! NumFactY is the number of factors in the Young operator,
    ! FirstLPos is the position of the first left bracket,
    ! LastRPos is the position of the last right bracket.
    !
    ! Each bracketed group is one factor. Anything standing OUTSIDE the
    ! brackets - before the first one, or between two of them - is a
    ! factor in its own right, which is what the inner scan below counts
    ! through k.
    !------------------------------------------------------------------
    IF (R /= 0) THEN

      FirstLPos = SCAN(Glob_YOperatorString(1:StrLen), '(')
      LastRPos = SCAN(Glob_YOperatorString(1:StrLen), ')', BACK=.TRUE.)
      NumFactY = R
      i = 0

      DO j = 1, R

        k = 0
        i = i+1
        c1 = Glob_YOperatorString(i:i)

        ! Walk up to the next '(' - a non-'*' character on the way
        ! means there is an unbracketed factor here.
        DO WHILE (c1 /= '(')
          IF (c1 /= '*') k = 1
          i = i+1
          c1 = Glob_YOperatorString(i:i)
        ENDDO

        IF (k == 1) NumFactY = NumFactY+1

        ! Skip to the matching ')'
        DO WHILE (c1 /= ')')
          i = i+1
          c1 = Glob_YOperatorString(i:i)
        ENDDO

      ENDDO

      ! A trailing unbracketed factor
      IF (Glob_YOperatorString(StrLen:StrLen) /= ')') NumFactY = NumFactY+1

    ELSE

      NumFactY = 1

    ENDIF


    !==================================================================
    ! Young operator, stage 5: split into one factor per YOpStr entry
    !==================================================================
    ! Each entry of YOpStr holds exactly one factor with the brackets
    ! removed. A '+' or '-' is prefixed to the first term of a factor
    ! when the source did not carry one, so that every term in every
    ! factor starts with an explicit sign. Multiplication signs are
    ! dropped.
    !------------------------------------------------------------------
    ALLOCATE(YOpStr(NumFactY))

    DO i = 1, NumFactY
      YOpStr(i) = ' '
    ENDDO

    IF (R == 0) THEN

      ! No brackets at all - the whole string is a single factor.
      c1 = Glob_YOperatorString(1:1)

      IF ((c1 /= '+') .OR. (c1 /= '-')) THEN
        YOpStr(1)(1:1) = '+'
        YOpStr(1)(2:StrLen+1) = Glob_YOperatorString(1:StrLen)
      ELSE
        YOpStr(1)(1:StrLen) = Glob_YOperatorString(1:StrLen)
      ENDIF

    ELSE

      i = 1
      k = 1
      p = i
      q = 0
      c1 = Glob_YOperatorString(i:i)

      ! Leading unbracketed factor, if any
      IF ((c1 /= '(') .AND. (c1 /= '+') .AND. (c1 /= '-')) THEN
        q = 1
        YOpStr(k)(1:1) = '+'
      ENDIF

      DO WHILE (Glob_YOperatorString(i:i) /= '(')
        i = i+1
      ENDDO

      IF (i > 1) THEN
        YOpStr(k)(p+q:i-1+q) = Glob_YOperatorString(p:i-1)
        ! if (YOpStr(k)(i-1+q:i-1+q)=='*') YOpStr(k)(i-1+q:i-1+q)=' '
        k = k+1
      ENDIF

      ! One pass per bracketed group
      DO j = 1, R

        i = i+1
        p = i
        q = 0
        c1 = Glob_YOperatorString(i:i)

        IF ((c1 /= ')') .AND. (c1 /= '+') .AND. (c1 /= '-')) THEN
          q = 1
          YOpStr(k)(1:1) = '+'
        ENDIF

        DO WHILE (Glob_YOperatorString(i:i) /= ')')
          i = i+1
        ENDDO

        YOpStr(k)(1+q:i+q-p) = Glob_YOperatorString(p:i-1)
        k = k+1
        i = i+1
        c1 = Glob_YOperatorString(i:i)

        IF (c1 == '*') THEN
          i = i+1
          c1 = Glob_YOperatorString(i:i)
        ENDIF

        ! An unbracketed factor sitting between this group and the next
        IF ((c1 /= '(') .AND. (c1 /= ' ')) THEN
          p = i
          YOpStr(k)(1:1) = '+'
          DO WHILE ((c1 /= '(') .AND. (c1 /= ' '))
            i = i+1
            c1 = Glob_YOperatorString(i:i)
          ENDDO
          YOpStr(k)(2:i+1-p) = Glob_YOperatorString(p:i-1)
          k = k+1
        ENDIF

      ENDDO

    ENDIF

    ! Print all factors in the Young operator
    ! j=StrLen+1
    ! do i=1,NumFactY
    !  write (*,'(1x,i3,1x,a3,a<j>)') i,':  ',YOpStr(i)(1:j)
    ! enddo


    !==================================================================
    ! Young operator, stage 6: the factors of Y^{+}
    !==================================================================
    ! Y^{+} is the reverse of Y: the order of the factors is reversed,
    ! and within each factor the permutation products come in reverse
    ! order too. Each permutation token 'Pij' is three characters wide,
    ! which is what the 3*(k-t) arithmetic below steps over.
    !------------------------------------------------------------------
    ALLOCATE(YHOpStr(NumFactY))

    DO i = 1, NumFactY
      YHOpStr(i) = ' '
    ENDDO

    DO i = NumFactY, 1, -1

      s = NumFactY-i+1
      j = 1
      c1 = YOpStr(s)(j:j)

      DO WHILE (c1 /= ' ')

        IF (c1 == 'P') THEN

          k = 0  ! k counts the number of Permutations in the current term
          t = 0

          DO WHILE ((c1 /= '+') .AND. (c1 /= '-') .AND. (c1 /= ' '))
            IF (c1 == 'P') k = k+1
            t = t+1
            c1 = YOpStr(s)(j+t:j+t)
          ENDDO

          ! Copy the k permutation tokens out in reverse order
          DO t = 1, k
            YHOpStr(i)(j+3*(k-t):j+3*(k-t)+2) = YOpStr(s)(j+3*(t-1):j+3*(t-1)+2)
          ENDDO

          j = j+3*k

        ELSE

          YHOpStr(i)(j:j) = c1
          j = j+1

        ENDIF

        c1 = YOpStr(s)(j:j)

      ENDDO

    ENDDO


    !==================================================================
    ! Young operator, stage 7: count the terms
    !==================================================================
    ! A term is anything introduced by a '+' or a '-'. The total number
    ! of terms in the non-simplified operator is the product over the
    ! factors, which is why it grows so quickly.
    !------------------------------------------------------------------
    ALLOCATE(NumTermsInYOpFact(NumFactY))

    TotNumOfYTerms = 1

    DO k = 1, NumFactY

      j = 0

      DO i = 1, Glob_YOperatorStringLength
        IF ((YOpStr(k)(i:i) == '+') .OR. (YOpStr(k)(i:i) == '-')) j = j+1
      ENDDO

      NumTermsInYOpFact(k) = j
      TotNumOfYTerms = TotNumOfYTerms*j
      TotNumOfYHYTerms = TotNumOfYTerms*TotNumOfYTerms

    ENDDO

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y operator:     ', TotNumOfYTerms
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y^{+}Y operator:', TotNumOfYHYTerms
    ENDIF

    ALLOCATE(Matr1(1:n, 1:n))
    ALLOCATE(Matr2(1:n, 1:n))
    ALLOCATE(Matr3(1:n, 1:n))
    ALLOCATE(Matr4(1:n, 1:n))


    !==================================================================
    ! Young operator, stage 8: expand Y and collect identical terms
    !==================================================================
    ! Multiplies out the factors of YOpStr into Glob_YMatr/Glob_YCoeff.
    ! The factor loop runs from NumFactY down to 1 (reversed operator
    ! product, see above). After every factor duplicate matrices have their
    ! coefficients summed into the first occurrence and are removed; t
    ! counts the survivors and the arrays are reallocated to that size.
    !------------------------------------------------------------------
    CurrNumOfTerms = NumTermsInYOpFact(NumFactY)
    ALLOCATE(TempSymCoeff(CurrNumOfTerms))
    ALLOCATE(TempSymMatr(n, n, CurrNumOfTerms))

    DO j = NumFactY, 1, -1

      ! reading the current factor
      k = 0
      i = 1
      c1 = YOpStr(j)(i:i)
      p = i

      DO WHILE (c1 /= ' ')

        i = i+1
        c1 = YOpStr(j)(i:i)

        ! Advance to the next token boundary
        DO WHILE ((c1 /= 'P') .AND. (c1 /= '+') .AND. (c1 /= '-') .AND. (i < Glob_YOperatorStringLength))
          i = i+1
          c1 = YOpStr(j)(i:i)
        ENDDO

        ! An explicit numeric coefficient, or just a bare sign
        IF (i-p > 1) THEN
          READ(YOpStr(j)(p:i-1), *) Coeff
        ELSE
          IF (YOpStr(j)(i-1:i-1) == '+') THEN
            Coeff = 1
          ELSE
            Coeff = -1
          ENDIF
        ENDIF

        ! Start from the identity and multiply in each Pij of this term
        Matr1 = Glob_Transposit(1:n, 1:n, 1, 1)

        DO WHILE (c1 == 'P')

          READ(YOpStr(j)(i+1:i+1), *) p
          READ(YOpStr(j)(i+2:i+2), *) q
          Matr2(1:n, 1:n) = Glob_Transposit(1:n, 1:n, p, q)
          Matr4(1:n, 1:n) = Matr1(1:n, 1:n)

          DO ii = 1, n
            DO jj = 1, n
              w = 0
              DO kk = 1, n
                w = w+Matr2(ii, kk)*Matr4(kk, jj)
              ENDDO
              Matr1(ii, jj) = w
            ENDDO
          ENDDO

          i = i+3
          c1 = YOpStr(j)(i:i)

        ENDDO

        k = k+1
        p = i

        IF (j /= NumFactY) THEN

          ! The first term of the factor is held back in Matr3/Cf3 and
          ! applied to the whole accumulator after the loop; the rest
          ! extend the accumulator in blocks of t.
          IF (k == 1) THEN
            Matr3(1:n, 1:n) = Matr1(1:n, 1:n)
            Cf3 = Coeff
          ELSE
            DO s = 1, t
              Matr2(1:n, 1:n) = TempSymMatr(1:n, 1:n, s)
              q = t*(k-1)+s
              DO ii = 1, n
                DO jj = 1, n
                  w = 0
                  DO kk = 1, n
                    w = w+Matr2(ii, kk)*Matr1(kk, jj)
                  ENDDO
                  TempSymMatr(ii, jj, q) = w
                ENDDO
              ENDDO
              TempSymCoeff(q) = Coeff*TempSymCoeff(s)
            ENDDO
          ENDIF

        ELSE

          ! The very first factor seeds the accumulator
          TempSymMatr(1:n, 1:n, k) = Matr1(1:n, 1:n)
          TempSymCoeff(k) = Coeff

        ENDIF

      ENDDO

      ! Apply the held-back first term of this factor
      IF (j /= NumFactY) THEN
        DO s = 1, t
          Matr2(1:n, 1:n) = TempSymMatr(1:n, 1:n, s)
          DO ii = 1, n
            DO jj = 1, n
              w = 0
              DO kk = 1, n
                w = w+Matr2(ii, kk)*Matr3(kk, jj)
              ENDDO
              TempSymMatr(ii, jj, s) = w
            ENDDO
          ENDDO
          TempSymCoeff(s) = Cf3*TempSymCoeff(s)
        ENDDO
      ENDIF

      ! mark the identical terms (adding their coefficients
      ! and setting all of them but one to zero)
      t = CurrNumOfTerms

      DO i = 1, CurrNumOfTerms

        IF (TempSymCoeff(i) == 0) CYCLE

        DO s = i+1, CurrNumOfTerms
          IF (TempSymCoeff(s) == 0) CYCLE
          IF (ALL(TempSymMatr(1:n, 1:n, i) == TempSymMatr(1:n, 1:n, s))) THEN
            TempSymCoeff(i) = TempSymCoeff(i)+TempSymCoeff(s)
            IF (TempSymCoeff(i) == 0) t = t-1
            TempSymCoeff(s) = 0
            t = t-1
          ENDIF
        ENDDO

      ENDDO

      ! reallocate arrays containing symmetry terms
      ! to allow for multiplication by the next factor
      IF (j /= 1) THEN

        ALLOCATE(TempSymCoeff1(t))
        ALLOCATE(TempSymMatr1(n, n, t))

        ! Compact the survivors into the temporaries
        s = 0
        DO i = 1, CurrNumOfTerms
          IF (TempSymCoeff(i) /= 0) THEN
            s = s+1
            TempSymCoeff1(s) = TempSymCoeff(i)
            TempSymMatr1(1:n, 1:n, s) = TempSymMatr(1:n, 1:n, i)
          ENDIF
        ENDDO

        CurrNumOfTerms = t*NumTermsInYOpFact(j-1)

        DEALLOCATE(TempSymCoeff)
        DEALLOCATE(TempSymMatr)
        ALLOCATE(TempSymCoeff(CurrNumOfTerms))
        ALLOCATE(TempSymMatr(n, n, CurrNumOfTerms))

        TempSymCoeff(1:t) = TempSymCoeff1(1:t)
        TempSymMatr(1:n, 1:n, 1:t) = TempSymMatr1(1:n, 1:n, 1:t)

        DEALLOCATE(TempSymCoeff1)
        DEALLOCATE(TempSymMatr1)

      ENDIF

    ENDDO

    ! Copy the survivors into the global Y arrays
    Glob_NumYTerms = t
    ALLOCATE(Glob_YCoeff(Glob_NumYTerms))
    ALLOCATE(Glob_YMatr(n, n, Glob_NumYTerms))

    s = 0
    DO i = 1, CurrNumOfTerms
      IF (TempSymCoeff(i) /= 0) THEN
        s = s+1
        Glob_YCoeff(s) = TempSymCoeff(i)
        Glob_YMatr(1:n, 1:n, s) = TempSymMatr(1:n, 1:n, i)
      ENDIF
    ENDDO

    DEALLOCATE(TempSymCoeff)
    DEALLOCATE(TempSymMatr)

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i9)') 'Total number of terms in simplified Y operator:        ', Glob_NumYTerms
    ENDIF


    !==================================================================
    ! Young operator, stage 9: expand Y^{+}Y the same way
    !==================================================================
    ! Multiply the factors of YHOpStr into the matrices and coefficients
    ! of Y that stage 8 just produced, so the accumulator starts out
    ! holding Y rather than the identity.
    !
    ! The same reversal rule applies: a product of permutation OPERATORS
    ! is the REVERSED product of the matrices, so the factor loop again
    ! runs backwards.
    !------------------------------------------------------------------
    CurrNumOfTerms = NumTermsInYOpFact(1)*Glob_NumYTerms
    ALLOCATE(TempSymCoeff(CurrNumOfTerms))
    ALLOCATE(TempSymMatr(n, n, CurrNumOfTerms))
    ! TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
    ! TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)
    TempSymCoeff(1:Glob_NumYTerms) = Glob_YCoeff(1:Glob_NumYTerms)
    TempSymMatr(1:n, 1:n, 1:Glob_NumYTerms) = Glob_YMatr(1:n, 1:n, 1:Glob_NumYTerms)
    t = Glob_NumYTerms

    DO j = NumFactY, 1, -1

      ! reading the current factor
      k = 0
      i = 1
      c1 = YHOpStr(j)(i:i)
      p = i

      DO WHILE (c1 /= ' ')

        i = i+1
        c1 = YHOpStr(j)(i:i)

        DO WHILE ((c1 /= 'P') .AND. (c1 /= '+') .AND. (c1 /= '-') .AND. (i < Glob_YOperatorStringLength))
          i = i+1
          c1 = YHOpStr(j)(i:i)
        ENDDO

        IF (i-p > 1) THEN
          READ(YHOpStr(j)(p:i-1), *) Coeff
        ELSE
          IF (YHOpStr(j)(i-1:i-1) == '+') THEN
            Coeff = 1
          ELSE
            Coeff = -1
          ENDIF
        ENDIF

        Matr1 = Glob_Transposit(1:n, 1:n, 1, 1)

        DO WHILE (c1 == 'P')

          READ(YHOpStr(j)(i+1:i+1), *) p
          READ(YHOpStr(j)(i+2:i+2), *) q
          Matr2(1:n, 1:n) = Glob_Transposit(1:n, 1:n, p, q)
          Matr4(1:n, 1:n) = Matr1(1:n, 1:n)

          DO ii = 1, n
            DO jj = 1, n
              w = 0
              DO kk = 1, n
                w = w+Matr2(ii, kk)*Matr4(kk, jj)
              ENDDO
              Matr1(ii, jj) = w
            ENDDO
          ENDDO

          i = i+3
          c1 = YHOpStr(j)(i:i)

        ENDDO

        k = k+1
        p = i

        ! Unlike stage 8 there is no special case for the first factor:
        ! the accumulator already holds Y, so every factor multiplies it.
        IF (k == 1) THEN
          Matr3(1:n, 1:n) = Matr1(1:n, 1:n)
          Cf3 = Coeff
        ELSE
          DO s = 1, t
            Matr2(1:n, 1:n) = TempSymMatr(1:n, 1:n, s)
            q = t*(k-1)+s
            DO ii = 1, n
              DO jj = 1, n
                w = 0
                DO kk = 1, n
                  w = w+Matr2(ii, kk)*Matr1(kk, jj)
                ENDDO
                TempSymMatr(ii, jj, q) = w
              ENDDO
            ENDDO
            TempSymCoeff(q) = Coeff*TempSymCoeff(s)
          ENDDO
        ENDIF

      ENDDO

      ! Apply the held-back first term of this factor
      DO s = 1, t
        Matr2(1:n, 1:n) = TempSymMatr(1:n, 1:n, s)
        DO ii = 1, n
          DO jj = 1, n
            w = 0
            DO kk = 1, n
              w = w+Matr2(ii, kk)*Matr3(kk, jj)
            ENDDO
            TempSymMatr(ii, jj, s) = w
          ENDDO
        ENDDO
        TempSymCoeff(s) = Cf3*TempSymCoeff(s)
      ENDDO

      ! mark the identical terms (adding their coefficients
      ! and setting all of them but one to zero)
      t = CurrNumOfTerms

      DO i = 1, CurrNumOfTerms

        IF (TempSymCoeff(i) == 0) CYCLE

        DO s = i+1, CurrNumOfTerms
          IF (TempSymCoeff(s) == 0) CYCLE
          IF (ALL(TempSymMatr(1:n, 1:n, i) == TempSymMatr(1:n, 1:n, s))) THEN
            TempSymCoeff(i) = TempSymCoeff(i)+TempSymCoeff(s)
            IF (TempSymCoeff(i) == 0) t = t-1
            TempSymCoeff(s) = 0
            t = t-1
          ENDIF
        ENDDO

      ENDDO

      ! reallocate arrays containing symmetry terms
      ! to allow for multiplication by the next factor
      IF (j /= 1) THEN

        ALLOCATE(TempSymCoeff1(t))
        ALLOCATE(TempSymMatr1(n, n, t))

        s = 0
        DO i = 1, CurrNumOfTerms
          IF (TempSymCoeff(i) /= 0) THEN
            s = s+1
            TempSymCoeff1(s) = TempSymCoeff(i)
            TempSymMatr1(1:n, 1:n, s) = TempSymMatr(1:n, 1:n, i)
          ENDIF
        ENDDO

        CurrNumOfTerms = t*NumTermsInYOpFact(NumFactY-j+2)

        DEALLOCATE(TempSymCoeff)
        DEALLOCATE(TempSymMatr)
        ALLOCATE(TempSymCoeff(CurrNumOfTerms))
        ALLOCATE(TempSymMatr(n, n, CurrNumOfTerms))

        TempSymCoeff(1:t) = TempSymCoeff1(1:t)
        TempSymMatr(1:n, 1:n, 1:t) = TempSymMatr1(1:n, 1:n, 1:t)

        DEALLOCATE(TempSymCoeff1)
        DEALLOCATE(TempSymMatr1)

      ENDIF

    ENDDO

    ! Copy the survivors into the global Y^{+}Y arrays
    Glob_NumYHYTerms = t
    ALLOCATE(Glob_YHYCoeff(Glob_NumYHYTerms))
    ALLOCATE(Glob_YHYMatr(n, n, Glob_NumYHYTerms))

    s = 0
    DO i = 1, CurrNumOfTerms
      IF (TempSymCoeff(i) /= 0) THEN
        s = s+1
        Glob_YHYCoeff(s) = TempSymCoeff(i)
        Glob_YHYMatr(1:n, 1:n, s) = TempSymMatr(1:n, 1:n, i)
      ENDIF
    ENDDO

    DEALLOCATE(TempSymCoeff)
    DEALLOCATE(TempSymMatr)

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i9)') 'Total number of terms in simplified Y^{+}Y operator:   ', Glob_NumYHYTerms
    ENDIF

    ! Debugging output of the Y and Y^{+}Y terms (file symterms_new.txt) was
    ! removed here; see NEW_workproc.f90.bak4_20260924 if it is needed again.

    ! Debugging dump of Glob_YHYCoeff/Glob_YHYMatr was removed here; see
    ! NEW_workproc.f90.bak4_20260924 if it is needed again.


    !==================================================================
    ! Release the Young-operator work arrays
    !==================================================================
    ! Deallocated in reverse order of allocation.
    !------------------------------------------------------------------
    DEALLOCATE(Matr4)
    DEALLOCATE(Matr3)
    DEALLOCATE(Matr2)
    DEALLOCATE(Matr1)
    DEALLOCATE(NumTermsInYOpFact)
    DEALLOCATE(YHOpStr)
    DEALLOCATE(YOpStr)


    !==================================================================
    ! Which particles are identical
    !==================================================================
    ! Determined from the masses and charges alone; needed to symmetrize
    ! two-particle expectation values. IdentParticleSet(i) labels the set
    ! of particle i (equal labels = identical particles); the largest label
    ! is the number of sets. Particle 1 is the reference particle, whose
    ! charge is Glob_PseudoCharge0 - hence the separate j==1 branch.
    !------------------------------------------------------------------
    ALLOCATE(IdentParticleSet(npart))

    IdentParticleSet(1) = 1
    k = 1

    DO i = 2, npart

      s = 0
      j = 0

      ! Look for an earlier particle with the same mass and charge
      DO WHILE ((j < i-1) .AND. (s == 0))

        j = j+1

        IF (j > 1) THEN
          IF ((Glob_Mass(j) == Glob_Mass(i)) .AND. (Glob_PseudoCharge(j-1) == Glob_PseudoCharge(i-1))) THEN
            IdentParticleSet(i) = IdentParticleSet(j)
            s = 1
          ENDIF
        ELSE
          ! j=1 case
          IF ((Glob_Mass(j) == Glob_Mass(i)) .AND. (Glob_PseudoCharge0 == Glob_PseudoCharge(i-1))) THEN
            IdentParticleSet(i) = IdentParticleSet(j)
            s = 1
          ENDIF
        ENDIF

      ENDDO

      ! No match - this particle opens a new set
      IF (s == 0) THEN
        k = k+1
        IdentParticleSet(i) = k
      ENDIF

    ENDDO

    Glob_NumOfIdentPartSets = MAXVAL(IdentParticleSet(1:npart))


    !==================================================================
    ! Which pairs of pseudoparticles are equivalent
    !==================================================================
    ! IdentPseudoPartPairSet(i,j) labels the equivalence class of the
    ! pair (i,j), so two pairs with the same label are equivalent. The
    ! largest label is the number of non-equivalent pairs.
    !
    ! A DIAGONAL element does not stand for a pair of pseudoparticles
    ! but for a single one, which itself corresponds to a pair of
    ! particles - which is why the i==j branch maps to particle 1.
    !------------------------------------------------------------------
    ALLOCATE(IdentPseudoPartPairSet(1:n, 1:n))

    IdentPseudoPartPairSet(1:n, 1:n) = 0
    k = 0

    DO i = 1, n
      DO j = i, n

        ! Particles behind the pseudoparticle pair (i,j)
        IF (i == j) THEN
          pi = 1; pj = j+1
        ELSE
          pi = i+1; pj = j+1
        ENDIF

        ! Scan the pairs already classified for an equivalent one
        w = 0

        DO s = 1, i

          IF (s == i) THEN
            q = j-1
          ELSE
            q = n
          ENDIF

          DO t = s, q

            IF (w == 1) CYCLE

            IF (s == t) THEN
              ps = 1; pt = t+1
            ELSE
              ps = s+1; pt = t+1
            ENDIF

            IF ((IdentParticleSet(ps) == IdentParticleSet(pi)) .AND. &
                (IdentParticleSet(pt) == IdentParticleSet(pj))) THEN
              w = 1
              IdentPseudoPartPairSet(i, j) = IdentPseudoPartPairSet(s, t)
            ENDIF

          ENDDO

        ENDDO

        ! No equivalent pair found - open a new class
        IF (w == 0) THEN
          k = k+1
          IdentPseudoPartPairSet(i, j) = k
        ENDIF

      ENDDO
    ENDDO

    Glob_NumOfNoneqvPairSets = MAXVAL(IdentPseudoPartPairSet(1:n, 1:n))


    !==================================================================
    ! Build Glob_NumOfPartInIdentPartSet and Glob_IdentPartList
    !==================================================================
    ! Invert IdentParticleSet into a per-set list of particle numbers.
    !------------------------------------------------------------------
    ALLOCATE(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
    ALLOCATE(Glob_IdentPartList(npart, Glob_NumOfIdentPartSets))

    Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets) = 0
    Glob_IdentPartList(1:npart, 1:Glob_NumOfIdentPartSets) = 0

    DO i = 1, npart
      k = IdentParticleSet(i)
      Glob_NumOfPartInIdentPartSet(k) = Glob_NumOfPartInIdentPartSet(k)+1
      Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k), k) = i
    ENDDO


    !==================================================================
    ! Build Glob_NumOfPairsInEqvPairSet and Glob_EqvPairList
    !==================================================================
    ! The same inversion for the pair classes. Glob_EqvPairList(1,m,k)
    ! and Glob_EqvPairList(2,m,k) hold the two pseudoparticles of the
    ! m-th pair in class k.
    !------------------------------------------------------------------
    ALLOCATE(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
    ALLOCATE(Glob_EqvPairList(2, n*(n+1)/2, Glob_NumOfNoneqvPairSets))

    Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets) = 0
    Glob_EqvPairList(1:2, 1:n*(n+1)/2, 1:Glob_NumOfNoneqvPairSets) = 0

    DO i = 1, n
      DO j = i, n
        k = IdentPseudoPartPairSet(i, j)
        Glob_NumOfPairsInEqvPairSet(k) = Glob_NumOfPairsInEqvPairSet(k)+1
        Glob_EqvPairList(1, Glob_NumOfPairsInEqvPairSet(k), k) = i
        Glob_EqvPairList(2, Glob_NumOfPairsInEqvPairSet(k), k) = j
      ENDDO
    ENDDO


    DEALLOCATE(IdentPseudoPartPairSet)
    DEALLOCATE(IdentParticleSet)


  END SUBROUTINE ProgramDataInit


  SUBROUTINE GenerateTrialParam(nfun, x, m, k, method_used)
    !==================================================================
    ! Subroutine GenerateTrialParam
    !==================================================================
    ! Generates the nonlinear parameters and premultiplier powers of nfun
    ! trial functions at once, driven by the GENERATOR_PARAM record:
    ! Glob_RG_p1 (probability of method 1), Glob_RG_s1 and Glob_RG_s2
    ! (spreads of methods 1 and 2). Both methods perturb a random existing
    ! function (nfun consecutive ones when nfun>1):
    !   Method 1 - each parameter is drawn from a normal distribution
    !              centred on the prototype's with deviation Glob_RG_s1*x(i).
    !   Method 2 - one normal draw r (0,Glob_RG_s2) scales every parameter
    !              by (1+r); redrawn while |1+r| lies in [0.8,1.2], which
    !              would give a nearly linearly dependent function.
    ! Powers are EVEN and lie in 2..PWRMax (PWRMax = largest even value
    ! <= Glob_MaxPowerAllowed). Method 1 keeps the prototype's power with
    ! probability PWRChangeProb, else draws a different even one; method 2
    ! always inherits it. Glob_IsIndexFixed overrides both with
    ! Glob_IndexFixedValue.
    ! Input:  nfun - number of trial functions.
    ! Output: x(1:Glob_npt,1:nfun), m(1:nfun); k - the (first) prototype
    !         function, 0 for an empty basis; method_used - 1, 2, or 0 for
    !         an empty basis.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: nfun               ! number of trial functions wanted
    REAL(wp) :: x(Glob_npt, nfun)  ! out: nonlinear parameters
    INTEGER  :: m(nfun)            ! out: premultiplier powers
    INTEGER  :: k                  ! out: prototype basis function
    INTEGER  :: method_used        ! out: 1, 2, or 0 for an empty basis

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j
    INTEGER :: p     ! declared, currently unused

    ! REAL(8) rather than REAL(wp) on purpose: DRNOR is declared
    ! real(8) in misc.f90 whatever PREC the rest of the build uses,
    ! so these match the generator they are fed from.
    REAL(8) :: r
    REAL(8) :: sumf  ! declared, currently unused

    ! Bounds of the uniform distribution used when the basis is empty
    REAL(8) :: Lmin = -0.5_8
    REAL(8) :: Lmax = 0.5_8

    ! Probability that a generated function KEEPS the power of its
    ! prototype instead of drawing a new one (method 1 only).
    REAL(8) :: PWRChangeProb = 0.5_8

    ! Largest EVEN premultiplier power this build allows. Powers index
    ! the gamma tables in globvars.f90, which are generated up to
    ! Glob_MaxPowerAllowed, so nothing above it may ever be produced.
    INTEGER, PARAMETER :: PWRMax = 2*(Glob_MaxPowerAllowed/2)


    IF (Glob_CurrBasisSize == 0) THEN

      !==================================================================
      ! Empty basis: draw everything uniformly
      !==================================================================
      ! There is no prototype to perturb, so the nonlinear parameters come
      ! from a flat distribution on [Lmin,Lmax] and the power from a flat
      ! distribution over the even values in 2..PWRMax.
      !------------------------------------------------------------------
      DO i = 1, nfun

        DO j = 1, Glob_np
          CALL RANDOM_NUMBER(r)
          x(j, i) = r*(Lmax-Lmin)+Lmin
        ENDDO

        CALL RANDOM_NUMBER(r)
        m(i) = 2*(1+INT(r*(PWRMax/2)))
        IF (m(i) > PWRMax) m(i) = PWRMax

      ENDDO

      method_used = 0
      k = 0

    ELSE

      IF (Glob_CurrBasisSize < nfun) THEN

        !==================================================================
        ! Basis shorter than the batch: every function picks its own prototype
        !==================================================================
        ! There are not enough existing functions to take a consecutive run
        ! of nfun, so each trial function draws its own prototype index k.
        ! On return k therefore holds the LAST prototype used.
        !------------------------------------------------------------------
        CALL RANDOM_NUMBER(r)

        IF (r < GLob_RG_p1) THEN

          ! method 1
          DO i = 1, nfun

            CALL RANDOM_NUMBER(r)
            k = INT(r*(Glob_CurrBasisSize))+1

            DO j = 1, Glob_npt
              x(j, i) = (Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j, k)
            ENDDO

            CALL RANDOM_NUMBER(r)

            IF (r > PWRChangeProb) THEN

              ! Keep the prototype's power
              m(i) = Glob_PWR(k)

            ELSE

              ! Draw a different even power
              j = Glob_PWR(k)
              m(i) = j
              DO WHILE (m(i) == j)
                CALL RANDOM_NUMBER(r)
                m(i) = 2*(1+INT(r*(PWRMax/2)))
                IF (m(i) > PWRMax) m(i) = PWRMax
              ENDDO

            ENDIF

          ENDDO

          method_used = 1

        ELSE

          ! method 2
          DO i = 1, nfun

            CALL RANDOM_NUMBER(r)
            k = INT(r*(Glob_CurrBasisSize))+1

            ! Redraw until the scale factor is far enough from 1
            r = Glob_RG_s2*drnor()+ONE
            DO WHILE ((ABS(r) > 0.8E0_wp) .AND. (ABS(r) < 1.2E0_wp))
              r = Glob_RG_s2*drnor()+ONE
            ENDDO

            DO j = 1, Glob_npt
              x(j, i) = r*Glob_NonlinParam(j, k)
            ENDDO

            m(i) = Glob_PWR(k)

          ENDDO

          method_used = 2

        ENDIF

      ELSE

        !==================================================================
        ! Normal case: take a consecutive run of nfun prototypes
        !==================================================================
        ! k is drawn once and indexes the FIRST function of the run, so the
        ! batch perturbs functions k, k+1, ... k+nfun-1.
        !------------------------------------------------------------------
        CALL RANDOM_NUMBER(r)
        k = INT(r*(Glob_CurrBasisSize-nfun+1))+1

        CALL RANDOM_NUMBER(r)

        IF (r < GLob_RG_p1) THEN

          ! method 1
          DO i = 1, nfun

            DO j = 1, Glob_npt
              x(j, i) = (Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j, k+i-1)
            ENDDO

            CALL RANDOM_NUMBER(r)

            ! The PWRMax<4 guard leaves fewer than two even powers to
            ! choose from, which would make the redraw loop below spin
            ! forever; in that case the prototype's power is kept.
            IF ((r > PWRChangeProb) .OR. (PWRMax < 4)) THEN

              m(i) = Glob_PWR(k+i-1)

            ELSE

              j = Glob_PWR(k+i-1)
              m(i) = j
              DO WHILE (m(i) == j)
                CALL RANDOM_NUMBER(r)
                m(i) = 2*(1+INT(r*(PWRMax/2)))
                IF (m(i) > PWRMax) m(i) = PWRMax
              ENDDO

            ENDIF

          ENDDO

          method_used = 1

        ELSE

          ! method 2
          ! One scale factor is drawn for the WHOLE batch here, unlike
          ! the short-basis branch above which draws one per function.
          r = Glob_RG_s2*drnor()+ONE
          DO WHILE ((ABS(r) > 0.8E0_wp) .AND. (ABS(r) < 1.2E0_wp))
            r = Glob_RG_s2*drnor()+ONE
          ENDDO

          DO i = 1, nfun

            DO j = 1, Glob_npt
              x(j, i) = r*Glob_NonlinParam(j, k+i-1)
            ENDDO

            m(i) = Glob_PWR(k+i-1)

          ENDDO

          method_used = 2

        ENDIF

      ENDIF

    ENDIF


    !==================================================================
    ! A fixed power overrides everything generated above
    !==================================================================
    ! If the power is fixed then we set it to the fixed value
    IF (Glob_IsIndexFixed) m(1:nfun) = Glob_IndexFixedValue


  END SUBROUTINE GenerateTrialParam


  SUBROUTINE ComputeOverlapPenalty(MaxPairOverlapPenalty, OverlapThreshold2, TotalPenalty)
    !==================================================================
    ! Subroutine ComputeOverlapPenalty
    !==================================================================
    ! Computes the pair-overlap penalty that full optimization may add to
    ! the energy to keep basis functions away from pair linear dependence:
    !   P = sum over pairs (i<j, j > Glob_nfru) of Pij,
    !   Pij = (Sij^2 - t^2)*b/(1-t^2) for Sij^2 > t^2, else 0,
    !   b = MaxPairOverlapPenalty, t^2 = OverlapThreshold2, Sij from Glob_S.
    ! P is zero below the threshold and rises towards b as |S| -> 1, so the
    ! objective stays finite and differentiable even at a point that already
    ! violates the threshold. Only pairs with at least one optimized function
    ! (Glob_nfru+1..Glob_nfa) are summed. The pair loops are split across
    ! ranks and combined with one MPI_ALLREDUCE; TotalPenalty is identical
    ! on every rank.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp) :: MaxPairOverlapPenalty  ! b, the penalty at |S| = 1
    REAL(wp) :: OverlapThreshold2      ! t^2, the squared threshold
    REAL(wp) :: TotalPenalty           ! out: P, same on every rank

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j, k   ! loop counters
    INTEGER :: nbands    ! whole bands of rows per rank
    INTEGER :: leftover  ! rows left after the whole bands
    LOGICAL :: oddband   ! alternates the sweep direction

    REAL(wp) :: pen_coeff  ! b/(1-t^2)
    REAL(wp) :: tp         ! this rank's partial sum


    !==================================================================
    ! Guards
    !==================================================================
    ! FullOpt1G/FullOpt1I never switch the penalty on with
    ! OverlapThreshold >= ONE, so ONE-OverlapThreshold2 cannot be zero; the
    ! tests repeat that check here and make a non-positive threshold or
    ! maximum penalty mean "no penalty". The early return is collective-safe
    ! because all three values are identical on every rank.
    !------------------------------------------------------------------
    TotalPenalty = ZERO
    IF (OverlapThreshold2 <= ZERO) RETURN
    IF (OverlapThreshold2 >= ONE) RETURN
    IF (MaxPairOverlapPenalty <= ZERO) RETURN

    tp = ZERO
    pen_coeff = MaxPairOverlapPenalty/(ONE-OverlapThreshold2)


    ! Local-work mode counts as a single process here: in that mode each rank
    ! holds its OWN trial function in Glob_S and needs the complete penalty for
    ! it, so the pair loops must NOT be split across ranks.
    IF ((Glob_NumOfProcs == 1) .OR. (Glob_LocalWorkMode)) THEN

      !==================================================================
      ! Serial path: one rank does every pair
      !==================================================================
      ! First the frozen-against-optimized block, then the triangle among
      ! the optimized functions themselves.
      !------------------------------------------------------------------
      ! In case of a single MPI process we do not split work
      DO i = 1, Glob_nfru
        DO j = Glob_nfru+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) &
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
        ENDDO
      ENDDO

      DO i = Glob_nfru+1, Glob_nfa
        DO j = i+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) &
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
        ENDDO
      ENDDO

    ELSE

      !==================================================================
      ! Parallel path: split the pairs across ranks
      !==================================================================
      ! The frozen-against-optimized rectangle is shared round robin over i.
      ! The triangle among the optimized functions (j=i+1..Glob_nfa) is dealt
      ! out in BANDS of Glob_NumOfProcs rows with the direction reversed on
      ! every other band (oddband), so each rank gets one long and one short
      ! row per pair of bands. nbands = whole bands, leftover = remaining rows.
      !------------------------------------------------------------------
      DO i = 1+Glob_ProcID, Glob_nfru, Glob_NumOfProcs
        DO j = Glob_nfru+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) &
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
        ENDDO
      ENDDO

      nbands = Glob_nfo/Glob_NumOfProcs
      leftover = MOD(Glob_nfo, Glob_NumOfProcs)
      oddband = .TRUE.

      DO k = 1, nbands

        ! Alternate the direction so long and short rows pair up
        IF (oddband) THEN
          i = Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
          oddband = .FALSE.
        ELSE
          i = Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
          oddband = .TRUE.
        ENDIF

        DO j = i+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) &
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
        ENDDO

      ENDDO  ! k

      ! The rows that did not fill a whole band. Only ranks whose index
      ! lands inside the remaining strip take one.
      IF (leftover > 1) THEN

        IF (oddband) THEN
          i = Glob_nfa-Glob_ProcID+1
        ELSE
          i = Glob_nfa-leftover+1+Glob_ProcID
        ENDIF

        IF ((i < Glob_nfa) .AND. (i > Glob_nfa-leftover)) THEN
          DO j = i+1, Glob_nfa
            IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) &
              tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
          ENDDO
        ENDIF

      ENDIF

    ENDIF


    !==================================================================
    ! Combine the partial sums
    !==================================================================
    ! Skipped in local-work mode: there the ranks no longer execute the
    ! same sequence of operations, so this collective would hang (see
    ! Glob_LocalWorkMode in globvars.f90). The new frame never sets that
    ! flag yet, so today this always takes the MPI_ALLREDUCE branch; the
    ! guard is here so the routine is already correct when BasisEnl brings
    ! the per-rank trial mode across.
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      TotalPenalty = tp
    ELSE
      CALL MPI_ALLREDUCE(tp, TotalPenalty, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF


  END SUBROUTINE ComputeOverlapPenalty


  SUBROUTINE ComputeOverlapPenaltyAndAddGradient(MaxPairOverlapPenalty, OverlapThreshold2, &
                                                 TotalPenalty, WkGR)
    !==================================================================
    ! Subroutine ComputeOverlapPenaltyAndAddGradient
    !==================================================================
    ! Same penalty as ComputeOverlapPenalty, plus its contribution to the
    ! energy gradient: dP/da = 2*pen_coeff*Sij*dSij/da with Sij the
    ! NORMALIZED overlap. Glob_D holds the derivative of the RAW overlap,
    ! without the derivative of the normalization factor (which cancels in
    ! dH - E*dS but not in a penalty on S alone), so that term is restored
    ! here: the -ONEHALF*...*Glob_S contributions, with
    ! temp2 = sqrt(Glob_diagS(i)/Glob_diagS(j)) multiplying the j-side and
    ! dividing the i-side correction. The gradient is NOT reduced here; the
    ! caller reduces WkGR once, after this returns.
    ! Result: TotalPenalty (combined, same on all ranks); WkGR (added to in
    ! place, still partial per rank).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp) :: MaxPairOverlapPenalty  ! b, the penalty at |S| = 1
    REAL(wp) :: OverlapThreshold2      ! t^2, the squared threshold
    REAL(wp) :: TotalPenalty           ! out: P, same on every rank
    REAL(wp) :: WkGR(:)                ! in/out: gradient, added to

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j, k   ! loop counters
    INTEGER :: m         ! nonlinear parameter index
    INTEGER :: nbands    ! whole bands of rows per rank
    INTEGER :: leftover  ! rows left after the whole bands
    LOGICAL :: oddband   ! alternates the sweep direction

    REAL(wp) :: pen_coeff  ! b/(1-t^2)
    REAL(wp) :: tp         ! this rank's partial sum
    REAL(wp) :: temp1      ! 2*pen_coeff*Sij
    REAL(wp) :: temp2      ! sqrt(diagS(i)/diagS(j))


    !==================================================================
    ! Guards
    !==================================================================
    ! FullOpt1G/FullOpt1I never switch the penalty on with
    ! OverlapThreshold >= ONE, so ONE-OverlapThreshold2 cannot be zero; the
    ! tests repeat that check here and make a non-positive threshold or
    ! maximum penalty mean "no penalty". The early return is collective-safe
    ! because all three values are identical on every rank; WkGR is left
    ! untouched, correctly.
    !------------------------------------------------------------------
    TotalPenalty = ZERO
    IF (OverlapThreshold2 <= ZERO) RETURN
    IF (OverlapThreshold2 >= ONE) RETURN
    IF (MaxPairOverlapPenalty <= ZERO) RETURN

    tp = ZERO
    pen_coeff = MaxPairOverlapPenalty/(ONE-OverlapThreshold2)


    ! Local-work mode counts as a single process here: in that mode each rank
    ! holds its OWN trial function in Glob_S and needs the complete penalty for
    ! it, so the pair loops must NOT be split across ranks.
    IF ((Glob_NumOfProcs == 1) .OR. (Glob_LocalWorkMode)) THEN

      !==================================================================
      ! Serial path: one rank does every pair
      !==================================================================
      ! Frozen-against-optimized first. Only the j side is accumulated
      ! there, because i is a frozen function and its parameters are not
      ! being varied.
      !------------------------------------------------------------------
      ! In case of a single MPI process we do not split work
      DO i = 1, Glob_nfru
        DO j = Glob_nfru+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) THEN
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
            temp1 = 2*pen_coeff*Glob_S(j, i)
            temp2 = SQRT(Glob_diagS(i)/Glob_diagS(j))
            DO m = 1, Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m) = WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, j-Glob_nfru, i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, j-Glob_nfru, j)*Glob_S(j, i)*temp2)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      ! Both functions of the pair are being optimized here, so BOTH
      ! sides of the gradient are accumulated.
      DO i = Glob_nfru+1, Glob_nfa
        DO j = i+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) THEN
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
            temp1 = 2*pen_coeff*Glob_S(j, i)
            temp2 = SQRT(Glob_diagS(i)/Glob_diagS(j))
            DO m = 1, Glob_npt
              WkGR((i-Glob_nfru-1)*Glob_npt+m) = WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, i-Glob_nfru, j) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, i-Glob_nfru, i)*Glob_S(j, i)/temp2)
            ENDDO
            DO m = 1, Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m) = WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, j-Glob_nfru, i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, j-Glob_nfru, j)*Glob_S(j, i)*temp2)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

    ELSE

      !==================================================================
      ! Parallel path: split the pairs across ranks
      !==================================================================
      ! Same balancing as ComputeOverlapPenalty: round robin over the
      ! rectangle, bands of Glob_NumOfProcs rows with alternating direction
      ! over the triangle (nbands whole bands, leftover remaining rows).
      !------------------------------------------------------------------
      DO i = 1+Glob_ProcID, Glob_nfru, Glob_NumOfProcs
        DO j = Glob_nfru+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) THEN
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
            temp1 = 2*pen_coeff*Glob_S(j, i)
            temp2 = SQRT(Glob_diagS(i)/Glob_diagS(j))
            DO m = 1, Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m) = WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, j-Glob_nfru, i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, j-Glob_nfru, j)*Glob_S(j, i)*temp2)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      nbands = Glob_nfo/Glob_NumOfProcs
      leftover = MOD(Glob_nfo, Glob_NumOfProcs)
      oddband = .TRUE.

      DO k = 1, nbands

        ! Alternate the direction so long and short rows pair up
        IF (oddband) THEN
          i = Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
          oddband = .FALSE.
        ELSE
          i = Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
          oddband = .TRUE.
        ENDIF

        DO j = i+1, Glob_nfa
          IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) THEN
            tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
            temp1 = 2*pen_coeff*Glob_S(j, i)
            temp2 = SQRT(Glob_diagS(i)/Glob_diagS(j))
            DO m = 1, Glob_npt
              WkGR((i-Glob_nfru-1)*Glob_npt+m) = WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, i-Glob_nfru, j) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, i-Glob_nfru, i)*Glob_S(j, i)/temp2)
            ENDDO
            DO m = 1, Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m) = WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m, j-Glob_nfru, i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m, j-Glob_nfru, j)*Glob_S(j, i)*temp2)
            ENDDO
          ENDIF
        ENDDO

      ENDDO  ! k

      ! The rows that did not fill a whole band. Only ranks whose index
      ! lands inside the remaining strip take one.
      IF (leftover > 1) THEN

        IF (oddband) THEN
          i = Glob_nfa-Glob_ProcID+1
        ELSE
          i = Glob_nfa-leftover+1+Glob_ProcID
        ENDIF

        IF ((i < Glob_nfa) .AND. (i > Glob_nfa-leftover)) THEN
          DO j = i+1, Glob_nfa
            IF (Glob_S(j, i)*Glob_S(j, i) > OverlapThreshold2) THEN
              tp = tp+pen_coeff*(Glob_S(j, i)*Glob_S(j, i)-OverlapThreshold2)
              temp1 = 2*pen_coeff*Glob_S(j, i)
              temp2 = SQRT(Glob_diagS(i)/Glob_diagS(j))
              DO m = 1, Glob_npt
                WkGR((i-Glob_nfru-1)*Glob_npt+m) = WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                  +temp1*(Glob_D(Glob_npt+m, i-Glob_nfru, j) &
                                                          -ONEHALF*Glob_D(Glob_npt+m, i-Glob_nfru, i)*Glob_S(j, i)/temp2)
              ENDDO
              DO m = 1, Glob_npt
                WkGR((j-Glob_nfru-1)*Glob_npt+m) = WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                  +temp1*(Glob_D(Glob_npt+m, j-Glob_nfru, i) &
                                                          -ONEHALF*Glob_D(Glob_npt+m, j-Glob_nfru, j)*Glob_S(j, i)*temp2)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

    ENDIF


    !==================================================================
    ! Combine the partial penalty sums
    !==================================================================
    ! Only the PENALTY is reduced; WkGR stays partial (see the header).
    ! Skipped in local-work mode, where the ranks do not execute the same
    ! sequence of collectives (never set in this frame yet).
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      TotalPenalty = tp
    ELSE
      CALL MPI_ALLREDUCE(tp, TotalPenalty, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF


  END SUBROUTINE ComputeOverlapPenaltyAndAddGradient


  FUNCTION EnergyGA(Nmin, Nmax, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Function EnergyGA
    !==================================================================
    ! Energy from DSYGVX (LAPACK) for the level Glob_WhichEigenvalue in a
    ! basis of Nmax functions. The nonlinear parameters are in
    ! Glob_NonlinParam; the matrix elements of the first Nmin-1 functions
    ! must be stored already (Nmin=0: none); AreMatElemNeeded=.FALSE. skips
    ! the matrix elements and only solves.
    ! ErrorCode: 0 on success; Nmax+i means the leading minor of S of size
    ! i is not positive definite; on failure the function returns 1e33.
    ! Rank 0 solves and broadcasts. In LOCAL-WORK MODE every rank holds a
    ! different trial function and must solve its own problem without
    ! collectives (Glob_LocalWorkMode is .FALSE. throughout this frame).
    !==================================================================

    IMPLICIT NONE

    REAL(wp) :: EnergyGA  ! function result

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    INTEGER :: Nmin, Nmax        ! basis range to work over
    LOGICAL :: AreMatElemNeeded  ! .FALSE. to solve only
    INTEGER :: ErrorCode         ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j               ! loop counters
    INTEGER :: NumOfEigvalsFound  ! M returned by DSYGVX

    REAL(wp) :: Evalue  ! the eigenvalue found
    REAL(wp) :: EVs(1)  ! W, holds the one eigenvalue
    REAL(wp) :: Z(1)    ! eigenvector placeholder

    ! LAPACK specifies IFAIL as dimension (N) for DSYGVX. Declaring it
    ! smaller lets the library write past the end of the array and corrupt
    ! the stack with any implementation that touches more than the first
    ! element: the Netlib reference writes only IFAIL(1:M), but that is not
    ! guaranteed and MKL differs. Z(1) is safe by contrast - JOBZ='N' here,
    ! so no eigenvector is ever referenced.
    INTEGER :: IFAIL(Glob_HSLeadDim)


    IF (AreMatElemNeeded) CALL ComputeMatElem(Nmin, Nmax)


    IF (Nmax == 1) THEN

      !==================================================================
      ! Single basis function: the energy is the one diagonal element
      !==================================================================
      EnergyGA = Glob_diagH(1)
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Mirror the lower triangles of H and S into the upper ones
      !==================================================================
      ! DSYGVX is called with UPLO='U', so it reads the upper triangle.
      ! The authoritative copies live in the lower triangles, and the
      ! diagonals are kept separately in Glob_diagH and are identically
      ! ONE for S because the basis functions are normalised.
      !------------------------------------------------------------------
      ! Copying H and S matrix elements from the lower triangles
      ! to the upper ones. The diagonals are copied from the global
      ! arrays where they are stored
      DO i = 1, Nmax
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO

      DO i = 1, Nmax
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO


      !==================================================================
      ! Solve the GSEP
      !==================================================================
      ! Rank 0 solves and broadcasts below. In local-work mode every rank
      ! solves its own problem instead - see the header note.
      !------------------------------------------------------------------
      IF (Glob_LocalWorkMode .OR. (Glob_ProcID == 0)) THEN

        CALL DSYGVX(1, 'N', 'I', 'U', Nmax, Glob_H, Glob_HSLeadDim, Glob_S, Glob_HSLeadDim, &
                    ZERO, ZERO, Glob_WhichEigenvalue, Glob_WhichEigenvalue, Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound, EVs, Z, Nmax, Glob_WorkForDSYGVX, &
                    Glob_LWorkForDSYGVX, Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
        !$      VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
        !$      LWORK, IWORK, IFAIL, INFO )
        Evalue = EVs(1)

      ENDIF

      ! ComputeOverlapPenalty carries its own local-work-mode guard.
      IF (Glob_OverlapPenaltyAllowed) CALL ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2, Glob_TotalOverlapPenalty)

      IF (.NOT. Glob_LocalWorkMode) &
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


      !==================================================================
      ! A failed solve is a rejected trial point, not a broken run
      !==================================================================
      ! DSYGVX failing during a stochastic search only means an ill-conditioned
      ! trial basis: return 1e33, which every caller rejects, instead of the
      ! meaningless Evalue. The sentinel is what SaveResults looks for before
      ! writing CURRENT_ENERGY. The caller counts failures and stops on EC0126
      ! when they exceed Glob_MaxFracOfTrialFailsAllowed. The early RETURN is
      ! collective-safe: ErrorCode was just broadcast.
      !------------------------------------------------------------------
      IF (ErrorCode /= 0) THEN

        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'Warning in EnergyGA: routine DSYGVX failed at this point'
          WRITE(*, *) 'ErrorCode = ', ErrorCode
          WRITE(*, *) 'Returning huge energy; the caller will reject this point'
        ENDIF

        EnergyGA = 1.0E33_wp
        Glob_EnergyGACounter = Glob_EnergyGACounter+1
        RETURN

      ENDIF


      !==================================================================
      ! Publish the eigenvalue and add the penalty
      !==================================================================
      IF (.NOT. Glob_LocalWorkMode) &
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      IF (Glob_OverlapPenaltyAllowed) THEN
        EnergyGA = Evalue+Glob_TotalOverlapPenalty
      ELSE
        EnergyGA = Evalue
      ENDIF

    ENDIF


    Glob_EnergyGACounter = Glob_EnergyGACounter+1


  END FUNCTION EnergyGA


  FUNCTION EnergyGAM(Nmin, Nmax, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Function EnergyGAM
    !==================================================================
    ! EnergyGA with JOBZ='V': also returns the linear coefficients in
    ! Glob_c, which the linear-coefficient rejection test needs. On failure
    ! it returns 1e33 AND zeroes Glob_c. Rank 0 solves and broadcasts the
    ! energy and Glob_c (the callers decide on the coefficients, so the
    ! ranks must agree); in LOCAL-WORK MODE no collectives run. DSYGVX
    ! destroys only the UPPER triangles of Glob_H and Glob_S, so they are
    ! rebuilt from the lower ones on every call.
    !==================================================================

    IMPLICIT NONE

    REAL(wp) :: EnergyGAM  ! function result

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    INTEGER :: Nmin, Nmax        ! basis range to work over
    LOGICAL :: AreMatElemNeeded  ! .FALSE. to solve only
    INTEGER :: ErrorCode         ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j               ! loop counters
    INTEGER :: NumOfEigvalsFound  ! M returned by DSYGVX

    REAL(wp) :: Evalue  ! the eigenvalue found
    REAL(wp) :: EVs(1)  ! W, holds the one eigenvalue

    ! LAPACK specifies IFAIL as dimension (N) for DSYGVX. Declaring it
    ! smaller lets the library write past the end of the array and corrupt
    ! the stack. It is not hypothetical here: JOBZ='V', so on failure LAPACK
    ! writes the indices of the eigenvectors that did not converge into
    ! IFAIL. See the same note in EnergyGA.
    INTEGER :: IFAIL(Glob_HSLeadDim)


    IF (AreMatElemNeeded) CALL ComputeMatElem(Nmin, Nmax)


    IF (Nmax == 1) THEN

      !==================================================================
      ! Single basis function: energy is the diagonal, coefficient is 1
      !==================================================================
      EnergyGAM = Glob_diagH(1)
      Glob_c(1) = ONE
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Mirror the lower triangles of H and S into the upper ones
      !==================================================================
      ! Copying H and S matrix elements from the lower triangles
      ! to the upper ones. The diagonals are copied from the global
      ! arrays where they are stored
      DO i = 1, Nmax
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO

      DO i = 1, Nmax
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO


      !==================================================================
      ! Solve the GSEP with JOBZ='V'
      !==================================================================
      ! The eigenvector is written straight into Glob_c.
      !------------------------------------------------------------------
      IF (Glob_LocalWorkMode .OR. (Glob_ProcID == 0)) THEN

        CALL DSYGVX(1, 'V', 'I', 'U', Nmax, Glob_H, Glob_HSLeadDim, Glob_S, Glob_HSLeadDim, &
                    ZERO, ZERO, Glob_WhichEigenvalue, Glob_WhichEigenvalue, Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound, EVs, Glob_c, Nmax, Glob_WorkForDSYGVX, &
                    Glob_LWorkForDSYGVX, Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
        !$      VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
        !$      LWORK, IWORK, IFAIL, INFO )
        Evalue = EVs(1)

      ENDIF

      ! ComputeOverlapPenalty carries its own local-work-mode guard.
      IF (Glob_OverlapPenaltyAllowed) CALL ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2, Glob_TotalOverlapPenalty)

      IF (.NOT. Glob_LocalWorkMode) &
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


      !==================================================================
      ! A failed solve is a rejected trial point, not a broken run
      !==================================================================
      ! Return 1e33 and ZERO Glob_c: zero coefficients cannot trip the
      ! linear-coefficient threshold, so a failed solve does not get a good
      ! function thrown away. The caller counts failures (EC0126); the early
      ! RETURN is collective-safe because ErrorCode was just broadcast.
      !------------------------------------------------------------------
      IF (ErrorCode /= 0) THEN

        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'Warning in EnergyGAM: routine DSYGVX failed at this point'
          WRITE(*, *) 'ErrorCode = ', ErrorCode
          WRITE(*, *) 'Returning huge energy; the caller will reject this point'
        ENDIF

        Glob_c(1:Nmax) = ZERO
        EnergyGAM = 1.0E33_wp
        RETURN

      ENDIF


      !==================================================================
      ! Publish the eigenvalue and coefficients, add the penalty
      !==================================================================
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Glob_c, Nmax, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF

      IF (Glob_OverlapPenaltyAllowed) THEN
        EnergyGAM = Evalue+Glob_TotalOverlapPenalty
      ELSE
        EnergyGAM = Evalue
      ENDIF

    ENDIF

    ! Deliberately does NOT bump Glob_EnergyGACounter: this is a
    ! diagnostic re-solve for the acceptance test, not an optimization
    ! step, and counting it would distort EnergyGA's call statistics.


  END FUNCTION EnergyGAM


  SUBROUTINE EnergyGB(Evalue, Gradient, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Subroutine EnergyGB
    !==================================================================
    ! Energy AND gradient with respect to the nonlinear parameters of the
    ! last Glob_nfo functions, with DSYGVX for the level
    ! Glob_WhichEigenvalue. Used by BasisEnlG, OptCycleG and FullOpt1G: the
    ! whole window of Glob_nfo functions is differentiated in one call.
    ! Matrix elements of the first Glob_nfru functions must be stored
    ! (else call ComputeMatElemAndDeriv first or pass
    ! AreMatElemNeeded=.TRUE.). Gradient = (dEdvechL_{nfru+1}, ...,
    ! dEdvechL_{nfa}).
    ! ErrorCode: 0 on success; <= Glob_nfa an eigenvector did not converge;
    ! Glob_nfa+i the leading minor of S of size i is not positive definite.
    ! On failure: huge energy and ZERO gradient. Rank 0 solves and
    ! broadcasts; the gradient loops are split across ranks and combined by
    ! one MPI_ALLREDUCE (all switched off in local-work mode).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp) :: Evalue                       ! out: the energy
    REAL(wp) :: Gradient(Glob_npt*Glob_nfo)  ! out: the gradient
    LOGICAL  :: AreMatElemNeeded             ! .FALSE. to solve only
    INTEGER  :: ErrorCode                    ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: nfo, nfa, nfru, npt  ! local copies of the window bounds
    INTEGER :: i, j, k, l, m        ! loop counters
    INTEGER :: nbands, leftover     ! (declared, currently unused)
    LOGICAL :: oddband              ! (declared, currently unused)
    INTEGER :: N                    ! (declared, currently unused)
    INTEGER :: NumOfEigvalsFound    ! M returned by DSYGVX

    ! First index and stride of the rank-split gradient loops. Set once
    ! below so that local-work mode can turn the splitting off without
    ! duplicating the loop bodies.
    INTEGER :: LoopFirst, LoopStride

    REAL(wp) :: EVs(1)                         ! W, holds the one eigenvalue
    REAL(wp) :: W(Glob_npt_MaxAllowed), t, t2  ! gradient accumulators
    REAL(wp) :: pen_coeff                      ! (declared, currently unused)

    ! LAPACK specifies IFAIL as dimension (N) for DSYGVX. Declaring it
    ! smaller lets the library write past the end of the array and corrupt
    ! the stack. It is not hypothetical here: JOBZ='V', so on failure LAPACK
    ! writes the indices of the eigenvectors that did not converge into
    ! IFAIL. See the same note in EnergyGA.
    INTEGER :: IFAIL(Glob_HSLeadDim)


    nfo = Glob_nfo
    nfa = Glob_nfa
    npt = Glob_npt
    nfru = Glob_nfru

    IF (AreMatElemNeeded) CALL ComputeMatElemAndDeriv(nfru+1, nfa)


    IF (nfa == 1) THEN

      !==================================================================
      ! Single basis function: energy is the Rayleigh quotient
      !==================================================================
      Evalue = Glob_diagH(1)/Glob_diagS(1)
      Glob_c(1) = ONE
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Mirror the lower triangles of H and S into the upper ones
      !==================================================================
      ! DSYGVX is called with UPLO='U' and destroys only the upper
      ! triangles; the matrix elements themselves live in the lower ones,
      ! which is why the upper triangles are rebuilt on every call.
      !------------------------------------------------------------------
      ! Copying H and S matrix elements from the lower triangles
      ! to the upper ones. The diagonals are copied from the global
      ! arrays where they are stored
      DO i = 1, nfa
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO

      DO i = 1, nfa
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO


      !==================================================================
      ! Solve the GSEP with JOBZ='V'
      !==================================================================
      ! The eigenvector is written straight into Glob_c, and the gradient
      ! below is built from it.
      !------------------------------------------------------------------
      IF (Glob_LocalWorkMode .OR. (Glob_ProcID == 0)) THEN

        CALL DSYGVX(1, 'V', 'I', 'U', nfa, Glob_H, Glob_HSLeadDim, Glob_S, Glob_HSLeadDim, &
                    ZERO, ZERO, Glob_WhichEigenvalue, Glob_WhichEigenvalue, Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound, EVs, Glob_c, nfa, Glob_WorkForDSYGVX, &
                    Glob_LWorkForDSYGVX, Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
        !$      VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
        !$      LWORK, IWORK, IFAIL, INFO )
        Evalue = EVs(1)

      ENDIF

      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Glob_c, nfa, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF

    ENDIF


    !==================================================================
    ! A failed solve is an unusable trial point, not a broken run
    !==================================================================
    ! On failure Glob_c holds garbage that the gradient loops would read at
    ! once, so hand the optimizer a huge energy (1e31, consumed by a line
    ! search) AND a zero gradient; the line search then backs out of the
    ! region where the basis went linearly dependent. Failures are counted
    ! by the caller (EC0126); the early RETURN is collective-safe.
    !------------------------------------------------------------------
    IF (ErrorCode /= 0) THEN

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Warning in EnergyGB: routine DSYGVX failed at this point'
        WRITE(*, *) 'ErrorCode = ', ErrorCode
        WRITE(*, *) 'Returning huge energy and zero gradient so the optimizer rejects this step'
      ENDIF

      Evalue = 1.0E31_wp
      Gradient(1:nfo*npt) = ZERO
      Glob_EnergyGBCounter = Glob_EnergyGBCounter+1
      RETURN

    ENDIF


    !==================================================================
    ! Computing gradient
    !==================================================================
    ! Two of the three loops are shared out across ranks and the
    ! MPI_ALLREDUCE at the end combines them: W(m) and TWO*t*W(m) are
    ! partial per rank, and the -t2*(...) correction is applied by each rank
    ! to its own slice of m only, so after the sum each m carries it once.
    ! In local-work mode LoopFirst/LoopStride are 1/1 and nothing is reduced.
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      LoopFirst = 1
      LoopStride = 1
    ELSE
      LoopFirst = 1+Glob_ProcID
      LoopStride = Glob_NumOfProcs
    ENDIF

    DO k = 1, nfo

      W(1:npt) = ZERO

      DO l = LoopFirst, nfa, LoopStride
        t = Glob_c(l)
        DO m = 1, npt
          W(m) = W(m)+t*(Glob_D(m, k, l)-Evalue*Glob_D(m+npt, k, l))
        ENDDO
      ENDDO

      t = Glob_c(k+nfru)
      t2 = t*t

      DO m = 1, npt
        Glob_WkGR((k-1)*npt+m) = TWO*t*W(m)
      ENDDO

      DO m = LoopFirst, npt, LoopStride
        Glob_WkGR((k-1)*npt+m) = Glob_WkGR((k-1)*npt+m)-t2*(Glob_D(m, k, k+nfru) &
                                                            -Evalue*Glob_D(m+npt, k, k+nfru))
      ENDDO

    ENDDO


    !==================================================================
    ! Overlap penalty, added to both the energy and the gradient
    !==================================================================
    ! ComputeOverlapPenaltyAndAddGradient adds into Glob_WkGR and does
    ! NOT reduce it - the reduction below does that once, which is why
    ! the penalty must be added before it.
    !------------------------------------------------------------------
    IF ((Glob_OverlapPenaltyAllowed) .AND. (nfa /= 1)) THEN
      CALL ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty, Glob_OverlapPenaltyThreshold2, &
                                               Glob_TotalOverlapPenalty, Glob_WkGR)
      Evalue = Evalue+Glob_TotalOverlapPenalty
    ENDIF


    !==================================================================
    ! Combine the partial gradients
    !==================================================================
    ! Skipped in local-work mode, where Glob_WkGR already holds this
    ! rank's complete gradient for its own trial function.
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      Gradient(1:nfo*npt) = Glob_WkGR(1:nfo*npt)
    ELSE
      CALL MPI_ALLREDUCE(Glob_WkGR, Gradient, nfo*npt, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF


    Glob_EnergyGBCounter = Glob_EnergyGBCounter+1


  END SUBROUTINE EnergyGB


  FUNCTION EnergyIA(Nmin, Nmax, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Function EnergyIA
    !==================================================================
    ! Energy from GSEPIIS (inverse iteration with a shift) in a basis of
    ! Nmax functions: the level returned is the one CLOSEST to
    ! Glob_ApproxEnergy, not the one chosen by index as in EnergyGA. The
    ! matrix elements of the first Nmin-1 functions must be stored (Nmin=0:
    ! none); AreMatElemNeeded=.FALSE. only solves. Glob_c is seeded with
    ! Glob_LastEigvector, so each solve starts from the previous answer.
    ! ErrorCode is passed through from GSEPIIS. Results are broadcast from
    ! rank 0 (not in local-work mode, never set in this frame).
    !==================================================================

    IMPLICIT NONE

    REAL(wp) :: EnergyIA  ! function result

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    INTEGER :: Nmin, Nmax        ! basis range to work over
    LOGICAL :: AreMatElemNeeded  ! .FALSE. to solve only
    INTEGER :: ErrorCode         ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    REAL(wp) :: Evalue           ! the eigenvalue found
    INTEGER  :: NumOfIterations  ! inverse iterations used
    LOGICAL  :: IsWrongState     ! converged on a level other than Glob_WhichEigenvalue


    IF (AreMatElemNeeded) CALL ComputeMatElem(Nmin, Nmax)

    IsWrongState = .FALSE.

    IF (Nmax == 1) THEN

      !==================================================================
      ! Single basis function: the energy is the one diagonal element
      !==================================================================
      EnergyIA = Glob_diagH(1)
      NumOfIterations = 1
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Solve by inverse iteration, seeded with the previous eigenvector
      !==================================================================
      Glob_c(1:Nmax) = Glob_LastEigvector(1:Nmax)

      CALL GSEPIIS(Nmin, Nmax, Glob_H, Glob_HSLeadDim, Glob_invD, Glob_S, Glob_HSLeadDim, &
                   Glob_ApproxEnergy, Glob_c, Glob_WorkForGSEPIIS, Glob_EigvalTol, &
                   Evalue, Glob_LastEigvector, Glob_LastEigvalTol, Glob_MaxIterForGSEPIIS, &
                   -1, NumOfIterations, ErrorCode)
      ! GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      ! lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)


      !==================================================================
      ! Track the accuracy actually achieved
      !==================================================================
      ! Glob_EigvalTol is what was ASKED for; Glob_LastEigvalTol is what
      ! the solver delivered. The best and worst seen over the run are
      ! carried in the data file so a restart keeps the history.
      !------------------------------------------------------------------
      IF (Glob_LastEigvalTol > Glob_WorstEigvalTol) Glob_WorstEigvalTol = Glob_LastEigvalTol
      IF (Glob_LastEigvalTol > Glob_BestEigvalTol) Glob_BestEigvalTol = Glob_LastEigvalTol

      ! ComputeOverlapPenalty carries its own local-work-mode guard.
      IF (Glob_OverlapPenaltyAllowed) CALL ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2, Glob_TotalOverlapPenalty)

      ! Skipped in local-work mode: each process is solving for its own trial
      ! function, so there is nothing to agree on and the collective would
      ! hang (see Glob_LocalWorkMode in globvars).
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF

      IF (Glob_OverlapPenaltyAllowed) THEN
        EnergyIA = Evalue+Glob_TotalOverlapPenalty
      ELSE
        EnergyIA = Evalue
      ENDIF


      !==================================================================
      ! Verdict on the solve
      !==================================================================
      ! IsEigenpairUsable accepts ErrorCode=2 when the residual is below
      ! Glob_EigvalTolUsable; IsRequestedEigenstate refuses a solve that
      ! converged on a level other than Glob_WhichEigenvalue (huge energy,
      ! ErrorCode=0, so it is not counted as a solver failure). Both verdicts
      ! are taken on rank 0 and broadcast.
      !------------------------------------------------------------------
      IF (IsEigenpairUsable(ErrorCode)) ErrorCode = 0
      IsWrongState = (ErrorCode == 0) .AND. (.NOT. IsRequestedEigenstate(Nmax))
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(IsWrongState, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
      IF (IsWrongState) THEN
        WrkP_WrongStateCount = WrkP_WrongStateCount+1
        EnergyIA = 1.0E31_wp
      ENDIF
      CALL ReportWrongStateOnce(Nmax, ErrorCode)

    ENDIF


    !==================================================================
    ! Call statistics
    !==================================================================
    ! Counter1 counts the calls and Counter2 accumulates the iterations,
    ! so their ratio is the average number of inverse iterations per
    ! solve - the figure reported at the end of an optimization step.
    !------------------------------------------------------------------
    Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIACounter = Glob_EnergyIACounter+1


  END FUNCTION EnergyIA


  FUNCTION EnergyIAM(Nmin, Nmax, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Function EnergyIAM
    !==================================================================
    ! EnergyIA with a properly NORMALIZED eigenvector, so that the linear
    ! coefficients in Glob_c are usable for the rejection test. Differences
    ! from EnergyIA: GSEPIIS is called with SpecifNorm=0 (normalization);
    ! Glob_c is refreshed from Glob_LastEigvector after the solve; for
    ! Nmax==1 the energy is Glob_H(1,1)+Glob_ApproxEnergy (Glob_H is
    ! SHIFTED on this path) and the single coefficient is ONE. ErrorCode is
    ! passed through from GSEPIIS; results are broadcast from rank 0.
    !==================================================================

    IMPLICIT NONE

    REAL(wp) :: EnergyIAM  ! function result

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    INTEGER :: Nmin, Nmax        ! basis range to work over
    LOGICAL :: AreMatElemNeeded  ! .FALSE. to solve only
    INTEGER :: ErrorCode         ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    REAL(wp) :: Evalue           ! the eigenvalue found
    INTEGER  :: NumOfIterations  ! inverse iterations used
    LOGICAL  :: IsWrongState     ! converged on a level other than Glob_WhichEigenvalue


    IF (AreMatElemNeeded) CALL ComputeMatElem(Nmin, Nmax)

    IsWrongState = .FALSE.

    IF (Nmax == 1) THEN

      !==================================================================
      ! Single basis function
      !==================================================================
      ! Glob_H holds the SHIFTED matrix here, so Glob_ApproxEnergy is
      ! added back to recover the actual energy. Note this differs from
      ! EnergyIA, which reads the unshifted Glob_diagH(1).
      !------------------------------------------------------------------
      EnergyIAM = Glob_H(1, 1)+Glob_ApproxEnergy
      Glob_c(1) = ONE
      NumOfIterations = 1
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Solve by inverse iteration, asking for a normalized eigenvector
      !==================================================================
      ! Glob_c is seeded with the previous eigenvector, which keeps the
      ! iteration count low because successive trial functions differ only
      ! slightly. SpecifNorm = 0 (rather than -1 as in EnergyIA) is what
      ! makes GSEPIIS normalize the result.
      !------------------------------------------------------------------
      Glob_c(1:Nmax) = Glob_LastEigvector(1:Nmax)

      CALL GSEPIIS(Nmin, Nmax, Glob_H, Glob_HSLeadDim, Glob_invD, Glob_S, Glob_HSLeadDim, &
                   Glob_ApproxEnergy, Glob_c, Glob_WorkForGSEPIIS, Glob_EigvalTol, &
                   Evalue, Glob_LastEigvector, Glob_LastEigvalTol, Glob_MaxIterForGSEPIIS, &
                   0, NumOfIterations, ErrorCode)
      ! GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      ! lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)

      ! Take the NORMALIZED vector, not the iterate Glob_c was seeded with
      Glob_c(1:Nmax) = Glob_LastEigvector(1:Nmax)


      !==================================================================
      ! Track the accuracy actually achieved
      !==================================================================
      ! Glob_EigvalTol is what was ASKED for; Glob_LastEigvalTol is what
      ! the solver delivered. The best and worst seen over the run are
      ! carried in the data file so a restart keeps the history.
      !------------------------------------------------------------------
      IF (Glob_LastEigvalTol > Glob_WorstEigvalTol) Glob_WorstEigvalTol = Glob_LastEigvalTol
      IF (Glob_LastEigvalTol > Glob_BestEigvalTol) Glob_BestEigvalTol = Glob_LastEigvalTol

      ! ComputeOverlapPenalty carries its own local-work-mode guard.
      IF (Glob_OverlapPenaltyAllowed) CALL ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2, Glob_TotalOverlapPenalty)

      ! Skipped in local-work mode: each process is solving for its own trial
      ! function, so there is nothing to agree on and the collective would
      ! hang (see Glob_LocalWorkMode in globvars).
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF

      IF (Glob_OverlapPenaltyAllowed) THEN
        EnergyIAM = Evalue+Glob_TotalOverlapPenalty
      ELSE
        EnergyIAM = Evalue
      ENDIF


      !==================================================================
      ! Verdict on the solve
      !==================================================================
      ! IsEigenpairUsable accepts ErrorCode=2 when the residual is below
      ! Glob_EigvalTolUsable; IsRequestedEigenstate refuses a solve that
      ! converged on a level other than Glob_WhichEigenvalue (huge energy,
      ! ErrorCode=0, so it is not counted as a solver failure). Both verdicts
      ! are taken on rank 0 and broadcast.
      !------------------------------------------------------------------
      IF (IsEigenpairUsable(ErrorCode)) ErrorCode = 0
      IsWrongState = (ErrorCode == 0) .AND. (.NOT. IsRequestedEigenstate(Nmax))
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(IsWrongState, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
      IF (IsWrongState) THEN
        WrkP_WrongStateCount = WrkP_WrongStateCount+1
        EnergyIAM = 1.0E31_wp
      ENDIF
      CALL ReportWrongStateOnce(Nmax, ErrorCode)

    ENDIF


    !==================================================================
    ! Call statistics
    !==================================================================
    ! Counter1 counts the calls and Counter2 accumulates the iterations,
    ! so their ratio is the average number of inverse iterations per
    ! solve. Note this shares EnergyIA's counters rather than keeping
    ! its own.
    !------------------------------------------------------------------
    Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIACounter = Glob_EnergyIACounter+1


  END FUNCTION EnergyIAM


  SUBROUTINE EnergyIB(Evalue, Gradient, AreMatElemNeeded, ErrorCode)
    !==================================================================
    ! Subroutine EnergyIB
    !==================================================================
    ! Energy AND gradient with respect to the nonlinear parameters of the
    ! last Glob_nfo functions, with GSEPIIS: the level returned is the one
    ! CLOSEST to Glob_ApproxEnergy (inverse-iteration twin of EnergyGB).
    ! Matrix elements of the first Glob_nfru functions must be stored (else
    ! call ComputeMatElemAndDeriv first or pass AreMatElemNeeded=.TRUE.).
    ! Gradient = (dEdvechL_{nfru+1}, ..., dEdvechL_{nfa}). ErrorCode is
    ! passed through from GSEPIIS; on failure: huge energy and ZERO
    ! gradient. GSEPIIS handles its own parallelism; the two broadcasts, the
    ! rank-split gradient loops and the reduction are skipped in local-work
    ! mode (never set in this frame).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp) :: Evalue                       ! out: the energy
    REAL(wp) :: Gradient(Glob_npt*Glob_nfo)  ! out: the gradient
    LOGICAL  :: AreMatElemNeeded             ! .FALSE. to solve only
    INTEGER  :: ErrorCode                    ! out: 0 on success

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: nfo, nfa, nfru, npt  ! local copies of the window bounds
    INTEGER :: i, j, k, l, m        ! loop counters
    INTEGER :: nbands, leftover     ! (declared, currently unused)
    LOGICAL :: oddband              ! (declared, currently unused)
    INTEGER :: NumOfIterations      ! inverse iterations used

    ! First index and stride of the rank-split gradient loops. Set once
    ! below so that local-work mode can turn the splitting off without
    ! duplicating the loop bodies.
    INTEGER :: LoopFirst, LoopStride

    REAL(wp) :: W(Glob_npt_MaxAllowed), t, t2  ! gradient accumulators
    REAL(wp) :: pen_coeff                      ! (declared, currently unused)
    LOGICAL  :: IsWrongState                   ! converged on a level other than Glob_WhichEigenvalue


    nfo = Glob_nfo
    nfa = Glob_nfa
    npt = Glob_npt
    nfru = Glob_nfru

    !==================================================================
    ! Refresh the inverse-iteration shift
    !==================================================================
    ! A shift frozen at the pre-optimization energy drifts away from the
    ! target eigenvalue as the energy descends and GSEPIIS convergence
    ! degrades. When the whole basis is optimized (nfru==0) every call
    ! recomputes and refactorizes everything, so moving the shift is free;
    ! with nfru>0 the reused leading block is valid for one shift only.
    !------------------------------------------------------------------
    IF (WrkP_RefreshShiftInIB .AND. (nfru == 0) .AND. AreMatElemNeeded) THEN
      IF (ABS(WrkP_LastINVITEnergy) < 1.0E10_wp) &
        Glob_ApproxEnergy = WrkP_LastINVITEnergy*Glob_InvItParameter
    ENDIF

    IF (AreMatElemNeeded) CALL ComputeMatElemAndDeriv(nfru+1, nfa)

    IsWrongState = .FALSE.

    IF (nfa == 1) THEN

      !==================================================================
      ! Single basis function
      !==================================================================
      ! Glob_H holds the SHIFTED matrix on this path, so Glob_ApproxEnergy
      ! is added back to recover the actual energy.
      !------------------------------------------------------------------
      Evalue = Glob_H(1, 1)+Glob_ApproxEnergy
      Glob_c(1) = ONE
      NumOfIterations = 1
      ErrorCode = 0

    ELSE

      !==================================================================
      ! Solve by inverse iteration, asking for a normalized eigenvector
      !==================================================================
      ! Glob_c is seeded with the previous eigenvector, which keeps the
      ! iteration count low because successive trial functions differ only
      ! slightly. SpecifNorm = 0 makes GSEPIIS normalize the result, which
      ! the gradient below needs.
      !------------------------------------------------------------------
      Glob_c(1:nfa) = Glob_LastEigvector(1:nfa)

      CALL GSEPIIS(nfru+1, nfa, Glob_H, Glob_HSLeadDim, Glob_invD, Glob_S, Glob_HSLeadDim, &
                   Glob_ApproxEnergy, Glob_c, Glob_WorkForGSEPIIS, Glob_EigvalTol, &
                   Evalue, Glob_LastEigvector, Glob_LastEigvalTol, Glob_MaxIterForGSEPIIS, &
                   0, NumOfIterations, ErrorCode)
      ! GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      ! lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)

      ! Take the NORMALIZED vector, not the iterate Glob_c was seeded with
      Glob_c(1:nfa) = Glob_LastEigvector(1:nfa)

      ! Glob_EigvalTol is what was ASKED for; Glob_LastEigvalTol is what
      ! the solver delivered. The best and worst seen over the run are
      ! carried in the data file so a restart keeps the history.
      IF (Glob_LastEigvalTol > Glob_WorstEigvalTol) Glob_WorstEigvalTol = Glob_LastEigvalTol
      IF (Glob_LastEigvalTol > Glob_BestEigvalTol) Glob_BestEigvalTol = Glob_LastEigvalTol

      ! Skipped in local-work mode: each process is solving for its own trial
      ! function, so there is nothing to agree on and the collective would
      ! hang (see Glob_LocalWorkMode in globvars).
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF

      !==================================================================
      ! Verdict on the solve
      !==================================================================
      ! See EnergyIA for the two tests. A point on the wrong level is
      ! handed back below with a huge energy and a zero gradient and
      ! ErrorCode=0, so the line search backs off without the caller
      ! counting a solver failure.
      !------------------------------------------------------------------
      IF (IsEigenpairUsable(ErrorCode)) ErrorCode = 0
      IsWrongState = (ErrorCode == 0) .AND. (.NOT. IsRequestedEigenstate(nfa))
      IF (.NOT. Glob_LocalWorkMode) THEN
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(IsWrongState, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
      CALL ReportWrongStateOnce(nfa, ErrorCode)

    ENDIF


    !==================================================================
    ! A failed solve is an unusable trial point, not a broken run
    !==================================================================
    ! Glob_c would feed garbage into the gradient loops, so hand the
    ! optimizer a huge energy (1e31, as EnergyGB) and a zero gradient; the
    ! line search backs out of the linearly dependent region. The early
    ! RETURN is collective-safe because ErrorCode was just broadcast.
    !------------------------------------------------------------------
    IF (ErrorCode /= 0) THEN

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Warning in EnergyIB: GSEP solution failed at this point'
        WRITE(*, *) 'ErrorCode = ', ErrorCode
        WRITE(*, *) 'Returning huge energy and zero gradient so the optimizer rejects this step'
      ENDIF

      Evalue = 1.0E31_wp
      Gradient(1:nfo*npt) = ZERO
      Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
      Glob_EnergyIBCounter = Glob_EnergyIBCounter+1
      RETURN

    ENDIF


    !==================================================================
    ! A point on the wrong level is refused, not counted as a failure
    !==================================================================
    ! See IsRequestedEigenstate for why this is not optional. The
    ! optimizer gets a huge energy and a zero gradient: the line search
    ! shortens the step and backs away from the region where a lower
    ! level crossed the shift.
    !------------------------------------------------------------------
    IF (IsWrongState) THEN
      WrkP_WrongStateCount = WrkP_WrongStateCount+1
      Evalue = 1.0E31_wp
      Gradient(1:nfo*npt) = ZERO
      Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
      Glob_EnergyIBCounter = Glob_EnergyIBCounter+1
      RETURN
    ENDIF

    ! Track the best energy for the shift refresh on subsequent calls
    IF (WrkP_RefreshShiftInIB .AND. (Evalue < WrkP_LastINVITEnergy)) WrkP_LastINVITEnergy = Evalue


    !==================================================================
    ! Computing gradient
    !==================================================================
    ! Two of the three loops are shared out across ranks and the reduction
    ! at the end combines them: W(m) and 2*t*W(m) are partial per rank, and
    ! the -t2*(...) correction is applied by each rank to its own slice of m
    ! only, so after the sum each m carries it once. In local-work mode
    ! LoopFirst/LoopStride are 1/1 and nothing is reduced.
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      LoopFirst = 1
      LoopStride = 1
    ELSE
      LoopFirst = 1+Glob_ProcID
      LoopStride = Glob_NumOfProcs
    ENDIF

    DO k = 1, nfo

      W(1:npt) = ZERO

      DO l = LoopFirst, nfa, LoopStride
        t = Glob_c(l)
        DO m = 1, npt
          W(m) = W(m)+t*(Glob_D(m, k, l)-Evalue*Glob_D(m+npt, k, l))
        ENDDO
      ENDDO

      t = Glob_c(k+nfru)
      t2 = t*t

      DO m = 1, npt
        Glob_WkGR((k-1)*npt+m) = 2*t*W(m)
      ENDDO

      DO m = LoopFirst, npt, LoopStride
        Glob_WkGR((k-1)*npt+m) = Glob_WkGR((k-1)*npt+m)-t2*(Glob_D(m, k, k+nfru) &
                                                            -Evalue*Glob_D(m+npt, k, k+nfru))
      ENDDO

    ENDDO


    !==================================================================
    ! Overlap penalty, added to both the energy and the gradient
    !==================================================================
    ! ComputeOverlapPenaltyAndAddGradient adds into Glob_WkGR and does
    ! NOT reduce it - the reduction below does that once, which is why
    ! the penalty must be added before it.
    !------------------------------------------------------------------
    IF ((Glob_OverlapPenaltyAllowed) .AND. (nfa /= 1)) THEN
      CALL ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty, Glob_OverlapPenaltyThreshold2, &
                                               Glob_TotalOverlapPenalty, Glob_WkGR)
      Evalue = Evalue+Glob_TotalOverlapPenalty
    ENDIF


    !==================================================================
    ! Combine the partial gradients
    !==================================================================
    ! Skipped in local-work mode, where Glob_WkGR already holds this
    ! rank's complete gradient for its own trial function.
    !------------------------------------------------------------------
    IF (Glob_LocalWorkMode) THEN
      Gradient(1:nfo*npt) = Glob_WkGR(1:nfo*npt)
    ELSE
      CALL MPI_ALLREDUCE(Glob_WkGR, Gradient, nfo*npt, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF


    !==================================================================
    ! Call statistics
    !==================================================================
    ! Counter1 counts the calls and Counter2 accumulates the iterations,
    ! so their ratio is the average number of inverse iterations per
    ! solve.
    !------------------------------------------------------------------
    Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIBCounter = Glob_EnergyIBCounter+1


  END SUBROUTINE EnergyIB


  SUBROUTINE ReadSwapFileAndDistributeData(IsSwapFileOK)
    !==================================================================
    ! Subroutine ReadSwapFileAndDistributeData
    !==================================================================
    ! Reads H and S back from the swap file (when it exists, Glob_UseSwapFile
    ! is set and the stored basis size equals Glob_CurrBasisSize) and
    ! broadcasts them, so a step resumes without recomputing the matrix
    ! elements; IsSwapFileOK reports success. ONE matrix is stored: lower
    ! part including the diagonal = H, upper part = S, followed by
    ! Glob_diagS (StoreMatricesInSwapFile writes exactly this; a swap file
    ! with any other layout is rejected and the elements are recomputed).
    ! Unpacking: 'G' copies S into the lower triangle of Glob_S and the
    ! diagonal of H into Glob_diagH; 'I' fills Glob_S symmetrically with a
    ! unit diagonal and converts Glob_H to the SHIFTED H - Glob_ApproxEnergy*S;
    ! 'Q' copies S into the lower triangle of Glob_S with a unit diagonal and
    ! leaves the unshifted H, diagonal included, in the lower triangle of Glob_H.
    ! The collectives are never guarded by Glob_LocalWorkMode: this runs at
    ! the start of a step and the matrices are identical on every rank.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    LOGICAL :: IsSwapFileOK  ! out: .TRUE. if the data was recovered

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j         ! loop counters
    INTEGER :: OpenFileErr  ! IOSTAT of the OPEN and of each READ


    IsSwapFileOK = .FALSE.


    !==================================================================
    ! Read the swap file (rank 0 only)
    !==================================================================
    ! Every READ carries IOSTAT, including the one that fetches the
    ! stored basis size: a truncated file, or one written by the other
    ! frame, must be REJECTED so the elements get recomputed, not raise
    ! an unhandled runtime error.
    !------------------------------------------------------------------
    IF (Glob_UseSwapFile) THEN

      IF (Glob_ProcID == 0) THEN

        OPEN(1, FILE=Glob_SwapFileName, FORM='unformatted', STATUS='old', IOSTAT=OpenFileErr)

        IF (OpenFileErr == 0) THEN

          READ(1, IOSTAT=OpenFileErr) i

          IF ((OpenFileErr == 0) .AND. (i == Glob_CurrBasisSize)) THEN

            ! Reading H and S from a matrix stored in file
            ! Remember that the lower part (including the diagonal)
            ! contains elements of H, while the upper part contains S
            IF (Verbose >= 2) WRITE(*, '(1x,a33)', ADVANCE='no') 'Reading H and S from swap file...'

            DO j = 1, Glob_CurrBasisSize
              IF (OpenFileErr == 0) THEN
                READ(1, IOSTAT=OpenFileErr) Glob_H(1:Glob_CurrBasisSize, j)
              ENDIF
            ENDDO

            ! reading diagS
            IF (OpenFileErr == 0) READ(1, IOSTAT=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize)

            IF (OpenFileErr == 0) THEN
              IsSwapFileOK = .TRUE.
              IF (Verbose >= 2) WRITE(*, *) 'completed'
              CLOSE(1)
            ELSE
              WRITE(*, *) 'failed'
            ENDIF

            ! erase information in swap file to free disc space
            OPEN(1, FILE=Glob_SwapFileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)
            WRITE(1) 'Swap file is empty'
            CLOSE(1)

          ENDIF

        ENDIF

      ENDIF

    ENDIF

    CALL MPI_BCAST(IsSwapFileOK, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)


    !==================================================================
    ! Distribute the matrices and unpack them
    !==================================================================
    ! If swap file is OK then send the data to all processes
    !------------------------------------------------------------------
    IF (IsSwapFileOK) THEN

      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a35)', ADVANCE='no') 'Sending H and S to all processes...'

      CALL MPI_BCAST(Glob_H, Glob_HSLeadDim*Glob_HSLeadDim, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_diagS, Glob_CurrBasisSize, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      ! Remember that the lower part (including the diagonal)
      ! contains elements of H, while the upper part contains S

      ! Restoring proper storage of data for Glob_GSEPSolutionMethod='G'
      IF (Glob_GSEPSolutionMethod == 'G') THEN
        DO i = 1, Glob_CurrBasisSize
          DO j = 1, i-1
            Glob_S(i, j) = Glob_H(j, i)
          ENDDO
          Glob_diagH(i) = Glob_H(i, i)
        ENDDO
      ENDIF

      ! Restoring proper storage of data for Glob_GSEPSolutionMethod='I'
      ! Inverse iteration works on the SHIFTED matrix, so H is converted
      ! to H - Glob_ApproxEnergy*S here. Glob_S must be filled first,
      ! because the shift of the off-diagonal H elements reads it.
      IF (Glob_GSEPSolutionMethod == 'I') THEN

        DO i = 1, Glob_CurrBasisSize
          DO j = 1, i-1
            Glob_S(i, j) = Glob_H(j, i)
            Glob_S(j, i) = Glob_H(j, i)
          ENDDO
          Glob_S(i, i) = ONE
        ENDDO

        DO i = 1, Glob_CurrBasisSize
          Glob_H(i, i) = Glob_H(i, i)-Glob_ApproxEnergy
          DO j = i+1, Glob_CurrBasisSize
            Glob_H(j, i) = Glob_H(j, i)-Glob_ApproxEnergy*Glob_S(j, i)
          ENDDO
        ENDDO

      ENDIF

      ! 'Q' uses the same compact layout as 'G' but keeps both diagonals in
      ! place: the lower triangles are the only authoritative storage.
      IF (Glob_GSEPSolutionMethod == 'Q') THEN
        DO i = 1, Glob_CurrBasisSize
          DO j = 1, i-1
            Glob_S(i, j) = Glob_H(j, i)
          ENDDO
          Glob_S(i, i) = ONE
        ENDDO
      ENDIF

      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'completed'

    ELSE

      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) 'Matrices H and S were not read from swap file'
        IF (Verbose >= 2) WRITE(*, *) 'All H and S matrix elements need be (re)computed'
      ENDIF

    ENDIF


  END SUBROUTINE ReadSwapFileAndDistributeData


  SUBROUTINE StoreMatricesInSwapFile()
    !==================================================================
    ! Subroutine StoreMatricesInSwapFile
    !==================================================================
    ! Writes H and S (with their diagonals) into the swap file when
    ! Glob_UseSwapFile is set, except on the last BBOP step. ONE matrix is
    ! written: lower part including the diagonal = H, upper part = S, then
    ! Glob_diagS; ReadSwapFileAndDistributeData reads exactly this layout
    ! (the two must change together).
    !*** SIDE EFFECT: Glob_H IS MODIFIED IN PLACE AND NOT RESTORED: S is
    ! copied into its upper triangle, for 'G' the diagonal is replaced by
    ! Glob_diagH, for 'I' the lower triangle and diagonal are UN-SHIFTED
    ! (Glob_ApproxEnergy added back); for 'Q' only the upper triangle is
    ! written, so its factors stay valid. Every call site deallocates the step
    ! workspace right afterwards, so this is safe by convention only; a call
    ! placed earlier in a step would corrupt Glob_H. Glob_S and Glob_diagS
    ! are only read.
    !==================================================================

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: i, j         ! loop counters
    INTEGER :: OpenFileErr  ! IOSTAT of the OPEN and of each WRITE


    IF (Glob_UseSwapFile) THEN

      ! If swap file is allowed to use then write H and S matrix
      ! elements into it
      IF (Glob_CurrBBOPStep /= Glob_NumOfBBOPSteps) THEN

        IF (Glob_ProcID == 0) THEN

          OPEN(1, FILE=Glob_SwapFileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)

          IF (OpenFileErr == 0) THEN

            WRITE(1) Glob_CurrBasisSize

            ! We store a matrix whose lower part (including the diagonal)
            ! contains elements of H, while the upper part contains S.
            ! Also, we store the diagonal of S.

            !==================================================================
            ! Pack S into the upper triangle of Glob_H - see the side-effect
            ! note in the header: Glob_H is left modified
            !==================================================================
            IF (Glob_GSEPSolutionMethod == 'G') THEN
              DO i = 1, Glob_CurrBasisSize
                DO j = 1, i-1
                  Glob_H(j, i) = Glob_S(i, j)
                ENDDO
                Glob_H(i, i) = Glob_diagH(i)
              ENDDO
            ENDIF

            ! On the 'I' path Glob_H holds the SHIFTED matrix, so the
            ! shift is added back here to store the true H. Each element
            ! is touched exactly once: column i takes S above the
            ! diagonal and the un-shifted H on and below it.
            IF (Glob_GSEPSolutionMethod == 'I') THEN
              DO i = 1, Glob_CurrBasisSize
                DO j = 1, i-1
                  Glob_H(j, i) = Glob_S(i, j)
                ENDDO
                Glob_H(i, i) = Glob_H(i, i)+Glob_ApproxEnergy
                DO j = i+1, Glob_CurrBasisSize
                  Glob_H(j, i) = Glob_H(j, i)+Glob_ApproxEnergy*Glob_S(j, i)
                ENDDO
              ENDDO
            ENDIF
            ! On the 'Q' path only the unused upper triangle of Glob_H is
            ! borrowed for S; the canonical lower triangles and the H
            ! diagonal stay as they are, so live factors still match.
            IF (Glob_GSEPSolutionMethod == 'Q') THEN
              DO i = 1, Glob_CurrBasisSize
                DO j = 1, i-1
                  Glob_H(j, i) = Glob_S(i, j)
                ENDDO
              ENDDO
            ENDIF


            !==================================================================
            ! Write it out, checking every operation
            !==================================================================
            ! A failed write is reported but not fatal: the file is simply not
            ! usable next time, and the matrix elements get recomputed.
            !------------------------------------------------------------------
            WRITE(*, *)
            WRITE(*, '(1x,a33)', ADVANCE='no') 'Writing H and S into swap file... '

            DO i = 1, Glob_CurrBasisSize
              IF (OpenFileErr == 0) THEN
                WRITE(1, IOSTAT=OpenFileErr) Glob_H(1:Glob_CurrBasisSize, i)
              ENDIF
            ENDDO

            IF (OpenFileErr == 0) THEN
              WRITE(1, IOSTAT=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize)
            ENDIF

            IF (OpenFileErr == 0) THEN
              IF (Verbose >= 2) WRITE(*, *) 'completed'
              IF (Verbose >= 2) WRITE(*, *)
            ELSE
              WRITE(*, *) 'failed to complete'
              WRITE(*, *)
            ENDIF

          ENDIF

          ! Outside the IF above on purpose: closing a unit that was
          ! never connected is legal and does nothing.
          CLOSE(1)

        ENDIF

      ENDIF

    ENDIF


  END SUBROUTINE StoreMatricesInSwapFile


  SUBROUTINE ReadHessianFile(V, IVLMAT, D, nvar, FileName, IsHessFileOK)
    !==================================================================
    ! Subroutine ReadHessianFile
    !==================================================================
    ! Reads a Hessian saved by SaveHessianFile into the DRMNG work array V
    ! at V(IVLMAT), packed as (nvar*(nvar+1))/2 elements, plus the scaling
    ! vector D when Glob_FullOptSaveD is set, so a full optimization resumes
    ! with its curvature information. IsHessFileOK is .TRUE. only if all
    ! was read. A file written for another nvar is refused with a WARNING
    ! (WC0100..WC0102), not an error: the caller falls back to the default
    ! initialization. Rank 0 only, by design: the caller runs this and
    ! DRMNG under IF (Glob_ProcID==0) and broadcasts IV and x afterwards.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp)                       :: V(*)          ! DRMNG work array; Hessian goes at V(IVLMAT)
    REAL(wp)                       :: D(*)          ! scaling vector, read only if Glob_FullOptSaveD
    INTEGER                        :: IVLMAT        ! offset of the Hessian inside V
    INTEGER                        :: nvar          ! number of variables in the problem
    CHARACTER(Glob_FileNameLength) :: FileName
    LOGICAL                        :: IsHessFileOK  ! out: .TRUE. if the Hessian was recovered

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i            ! (declared, currently unused)
    INTEGER :: j            ! nvar as recorded in the file
    INTEGER :: OpenFileErr  ! IOSTAT of the OPEN and of each READ


    IsHessFileOK = .FALSE.

    IF (Glob_ProcID == 0) THEN

      OPEN(1, FILE=FileName, FORM='unformatted', STATUS='old', IOSTAT=OpenFileErr)

      IF (OpenFileErr == 0) THEN

        !==================================================================
        ! The file opened - read it
        !==================================================================
        ! OpenFileErr is reused from here on to carry the IOSTAT of each
        ! READ, so it can no longer be used to decide whether the unit is
        ! connected. The unit is closed unconditionally at the end of this
        ! branch instead - see the note there.
        !------------------------------------------------------------------
        READ(1, IOSTAT=OpenFileErr) j

        IF ((j == nvar) .AND. (OpenFileErr == 0)) THEN

          IF (Verbose >= 2) WRITE(*, '(1x,a28)', ADVANCE='no') 'Reading Hessian from file...'

          READ(1, IOSTAT=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)

          IF ((OpenFileErr == 0) .AND. (Glob_FullOptSaveD)) READ(1, IOSTAT=OpenFileErr) D(1:nvar)

          IF (OpenFileErr == 0) THEN
            IsHessFileOK = .TRUE.
            IF (Verbose >= 2) WRITE(*, *) 'done'
          ELSE
            WRITE(*, *) 'failed'
            IF (Verbose >= 1) WRITE(*, *) 'Warning WC0100: Default Hessian initialization must be used'
          ENDIF

        ELSE

          ! Wrong size, or the header could not be read at all. The usual
          ! cause is simply that the basis has grown since the file was
          ! written.
          IF (Verbose >= 1) WRITE(*, *) 'Warning WC0101: Hessian file ', FileName
          WRITE(*, *) 'is inconsistent with the dimension of the current optimization problem'
          WRITE(*, *) 'and will not be read'
          IF (Verbose >= 2) WRITE(*, *) 'Default Hessian initialization must be used'

        ENDIF

        ! Close on EVERY path that got the file open, not just the one
        ! that read it successfully. Unit 1 is shared with the data file,
        ! the swap file, the black list and SaveResults, so a unit left
        ! connected here makes the NEXT unrelated OPEN(1,...) fail - and
        ! the most common way to land in this routine's failure paths is
        ! a merely stale Hessian file, which must not break anything.
        CLOSE(1)

      ELSE

        IF (Verbose >= 1) WRITE(*, *) 'Warning WC0102: Cannot open Hessian file ', FileName
        IF (Verbose >= 2) WRITE(*, *) 'Default Hessian initialization must be used'

      ENDIF

    ENDIF


  END SUBROUTINE ReadHessianFile


  SUBROUTINE SaveHessianFile(V, IVLMAT, D, nvar, FileName, IsSuccess)
    !==================================================================
    ! Subroutine SaveHessianFile
    !==================================================================
    ! Writes the DRMNG Hessian from V(IVLMAT) to FileName, packed as
    ! (nvar*(nvar+1))/2 elements, preceded by nvar so the reader can refuse
    ! a file of another problem size, and followed by the scaling vector D
    ! when Glob_FullOptSaveD is set (LAST, so a run that no longer wants D
    ! can still read the Hessian). ReadHessianFile reads exactly this
    ! layout. IsSuccess is .TRUE. only if every write succeeded. Rank 0
    ! only, like ReadHessianFile.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    REAL(wp)                       :: V(*)       ! DRMNG work array; Hessian starts at V(IVLMAT)
    REAL(wp)                       :: D(*)       ! scaling vector, written if Glob_FullOptSaveD
    INTEGER                        :: IVLMAT     ! offset of the Hessian inside V
    INTEGER                        :: nvar       ! number of variables in the problem
    CHARACTER(Glob_FileNameLength) :: FileName
    LOGICAL                        :: IsSuccess  ! out: .TRUE. if the file was written

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i            ! (declared, currently unused)
    INTEGER :: OpenFileErr  ! IOSTAT of the OPEN and of each WRITE


    IsSuccess = .FALSE.

    IF (Glob_ProcID == 0) THEN

      OPEN(1, FILE=FileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)

      IF (OpenFileErr == 0) THEN

        !==================================================================
        ! The file opened - write it
        !==================================================================
        ! OpenFileErr is reused from here on to carry the IOSTAT of each
        ! WRITE, so it can no longer say whether the unit is connected. The
        ! unit is closed unconditionally at the end of this branch instead.
        !------------------------------------------------------------------
        WRITE(1, IOSTAT=OpenFileErr) nvar

        IF (OpenFileErr == 0) THEN

          IF (Verbose >= 2) WRITE(*, '(1x,a17)', ADVANCE='no') 'Saving Hessian...'

          WRITE(1, IOSTAT=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)

          ! The OpenFileErr==0 test is what stops a FAILED Hessian write
          ! from being masked by a SUCCEEDING D write: without it a
          ! successful second write resets OpenFileErr and the routine
          ! reports success over a truncated file. ReadHessianFile guards
          ! the matching read the same way.
          IF ((OpenFileErr == 0) .AND. (Glob_FullOptSaveD)) WRITE(1, IOSTAT=OpenFileErr) D(1:nvar)

          IF (OpenFileErr == 0) THEN
            IsSuccess = .TRUE.
            IF (Verbose >= 2) WRITE(*, *) 'done'
          ELSE
            WRITE(*, *) 'failed'
          ENDIF

        ENDIF

        ! Close on EVERY path that got the file open, including the one
        ! where the nvar write failed and the block above was skipped.
        ! Unit 1 is shared with the data file, the swap file, the black
        ! list and SaveResults, so a unit left connected here makes the
        ! next unrelated OPEN(1,...) fail.
        CLOSE(1)

      ENDIF

      IF ((.NOT. IsSuccess) .AND. (Verbose >= 1)) WRITE(*, *) 'Warning WC0105: Hessian was not saved'

    ENDIF


  END SUBROUTINE SaveHessianFile


  SUBROUTINE PermuteFunctions(fb, fe, FuncNumTemp, NonlinParamTemp)
    !==================================================================
    ! Subroutine PermuteFunctions
    !==================================================================
    ! Moves basis functions fb..fe to the END of the basis and slides the
    ! functions that followed them down, so that the functions being
    ! optimized are the trailing block that the matrix element and gradient
    ! routines recompute. Glob_NonlinParam, Glob_PWR and Glob_FuncNum are
    ! permuted together; the matrix elements are handled by
    ! PermuteMatrixElements.
    ! Workspace (caller): FuncNumTemp (>= fe-fb+1), NonlinParamTemp
    ! (>= (fe-fb+1)*Glob_npt). nfco = fe-fb+1 is the block size and
    ! fbn = fb + Glob_CurrBasisSize - fe its destination. The shift loop
    ! reads fe+1..Glob_CurrBasisSize and writes fb..fbn-1 ascending, reading
    ! ahead of where it writes; run backwards it would corrupt the basis.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: fb, fe                        ! block to move to the end
    INTEGER  :: FuncNumTemp(*)                ! integer scratch
    REAL(wp) :: NonlinParamTemp(Glob_npt, *)  ! parameter scratch

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i     ! loop counter
    INTEGER :: fbn   ! where the moved block starts
    INTEGER :: nfco  ! size of the moved block


    ! The caller must pass a non-empty range lying inside the current
    ! basis. A bad range here would not fail loudly - it would quietly
    ! scramble Glob_NonlinParam and every later result with it - so it
    ! is checked rather than assumed. Every rank runs this routine, so
    ! all of them take this branch together.
    IF ((fb < 1) .OR. (fb > fe) .OR. (fe > Glob_CurrBasisSize)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0239 in PermuteFunctions: invalid function range'
        WRITE(*, *) 'fb =', fb, ' fe =', fe, ' Glob_CurrBasisSize =', Glob_CurrBasisSize
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    ENDIF


    nfco = fe-fb+1
    fbn = fb+Glob_CurrBasisSize-fe


    !==================================================================
    ! Nonlinear parameters
    !==================================================================
    NonlinParamTemp(1:Glob_npt, 1:nfco) = Glob_NonlinParam(1:Glob_npt, fb:fe)

    DO i = fb, fbn-1
      Glob_NonlinParam(1:Glob_npt, i) = Glob_NonlinParam(1:Glob_npt, nfco+i)
    ENDDO

    Glob_NonlinParam(1:Glob_npt, fbn:Glob_CurrBasisSize) = NonlinParamTemp(1:Glob_npt, 1:nfco)


    !==================================================================
    ! Premultiplier powers
    !==================================================================
    ! FuncNumTemp is reused as scratch here and again below; the two uses
    ! are sequential, so one array of size nfco serves both.
    !------------------------------------------------------------------
    FuncNumTemp(1:nfco) = Glob_PWR(fb:fe)

    DO i = fb, fbn-1
      Glob_PWR(i) = Glob_PWR(nfco+i)
    ENDDO

    Glob_PWR(fbn:Glob_CurrBasisSize) = FuncNumTemp(1:nfco)


    !==================================================================
    ! Function numbers
    !==================================================================
    FuncNumTemp(1:nfco) = Glob_FuncNum(fb:fe)

    DO i = fb, fbn-1
      Glob_FuncNum(i) = Glob_FuncNum(nfco+i)
    ENDDO

    Glob_FuncNum(fbn:Glob_CurrBasisSize) = FuncNumTemp(1:nfco)


  END SUBROUTINE PermuteFunctions


  SUBROUTINE PermuteFunctions2(fb1, fe1, fe2, FuncNumTemp, NonlinParamTemp)
    !==================================================================
    ! Subroutine PermuteFunctions2
    !==================================================================
    ! Exchanges the two ADJACENT blocks fb1..fe1 and fe1+1..fe2 of basis
    ! functions (Glob_NonlinParam, Glob_PWR, Glob_FuncNum together); the
    ! matrix elements are handled by PermuteMatrixElements2.
    ! Workspace (caller): FuncNumTemp (>= fe1-fb1+1, the FIRST block),
    ! NonlinParamTemp (>= (fe1-fb1+1)*Glob_npt). Only the first block is
    ! saved: k = fe1-fb1+1, fbn = fb1 + fe2 - fe1; the second block slides
    ! down into fb1..fbn-1 (ascending loop, reading ahead of the writes)
    ! and the saved block is dropped into fbn..fe2.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: fb1, fe1  ! first block
    INTEGER :: fe2       ! end of the second block
    ! (it begins at fe1+1)
    INTEGER  :: FuncNumTemp(*)                ! integer scratch
    REAL(wp) :: NonlinParamTemp(Glob_npt, *)  ! parameter scratch

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i    ! loop counter
    INTEGER :: fbn  ! where the first block lands
    INTEGER :: k    ! size of the first block


    ! Both blocks must be non-empty and lie inside the current basis:
    ! fb1 <= fe1 for the first, fe1 < fe2 for the second. A bad range
    ! would quietly scramble the basis rather than fail, so it is
    ! checked. Every rank runs this routine, so all of them take this
    ! branch together.
    IF ((fb1 < 1) .OR. (fb1 > fe1) .OR. (fe1 >= fe2) .OR. (fe2 > Glob_CurrBasisSize)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0240 in PermuteFunctions2: invalid function range'
        WRITE(*, *) 'fb1 =', fb1, ' fe1 =', fe1, ' fe2 =', fe2, &
                   ' Glob_CurrBasisSize =', Glob_CurrBasisSize
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    ENDIF


    k = fe1-fb1+1
    fbn = fb1+fe2-fe1


    !==================================================================
    ! Nonlinear parameters
    !==================================================================
    NonlinParamTemp(1:Glob_npt, 1:k) = Glob_NonlinParam(1:Glob_npt, fb1:fe1)

    DO i = fb1, fbn-1
      Glob_NonlinParam(1:Glob_npt, i) = Glob_NonlinParam(1:Glob_npt, k+i)
    ENDDO

    Glob_NonlinParam(1:Glob_npt, fbn:fe2) = NonlinParamTemp(1:Glob_npt, 1:k)


    !==================================================================
    ! Premultiplier powers
    !==================================================================
    FuncNumTemp(1:k) = Glob_PWR(fb1:fe1)

    DO i = fb1, fbn-1
      Glob_PWR(i) = Glob_PWR(k+i)
    ENDDO

    Glob_PWR(fbn:fe2) = FuncNumTemp(1:k)


    !==================================================================
    ! Function numbers
    !==================================================================
    FuncNumTemp(1:k) = Glob_FuncNum(fb1:fe1)

    DO i = fb1, fbn-1
      Glob_FuncNum(i) = Glob_FuncNum(k+i)
    ENDDO

    Glob_FuncNum(fbn:fe2) = FuncNumTemp(1:k)


  END SUBROUTINE PermuteFunctions2


  SUBROUTINE PermuteMatrixElements(fb, fe, TempR)
    !==================================================================
    ! Subroutine PermuteMatrixElements
    !==================================================================
    ! Applies to Glob_H, Glob_S, Glob_diagH and Glob_diagS the permutation
    ! PermuteFunctions applies to the functions (block fb..fe moved to the
    ! end). Call the two together. Workspace: TempR (>= fe-fb+1).
    ! Only the LOWER triangles (excluding the diagonals) are meaningful;
    ! parts of the UPPER triangles serve as scratch, which makes the
    ! permutation possible in place. The H diagonal is in Glob_diagH for
    ! 'G' and in Glob_H(i,i) for 'I'. With nfco = fe-fb+1 and
    ! k = Glob_CurrBasisSize-fe the work runs in six stages: shift within
    ! columns 1..fb-1; stage the tail x tail block and the transposed
    ! tail x block block into the upper triangle; move the within-block
    ! triangle IN PLACE; drop both staged blocks into their new places. For
    ! 'I' the upper triangle of Glob_S is rebuilt at the end.
    ! The in-place stage
    !     DO i=fe-1,fb,-1
    !       Glob_H(i+k+1:N,i+k) = Glob_H(i+1:fe,i)
    !     ENDDO
    ! MUST run backwards: it writes column i+k (k >= 0) and reads column i,
    ! so descending i never reads a column an earlier iteration replaced
    ! (Fortran protects the right-hand side only within one statement).
    ! Ascending, it failed whenever k <= nfco-2 with silently wrong H and S;
    ! both directions were checked by simulation for N=4..12.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: fb, fe    ! block to move to the end
    REAL(wp) :: TempR(*)  ! scratch, at least fe-fb+1 long

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j        ! loop counters
    INTEGER :: fbn         ! where the moved block starts
    INTEGER :: nfco        ! size of the moved block
    INTEGER :: fep         ! fe+1, first function of the tail
    INTEGER :: k           ! size of the tail that slides down
    INTEGER :: q, p, r, s  ! index offsets for the staged blocks


    nfco = fe-fb+1
    fbn = fb+Glob_CurrBasisSize-fe


    !==================================================================
    ! Diagonals of H and S
    !==================================================================
    ! 'G' keeps the H diagonal in Glob_diagH; 'I' keeps it inside
    ! Glob_H itself. Glob_diagS is permuted in both cases.
    !------------------------------------------------------------------
    ! First we do the diagonals of H and S
    IF (Glob_GSEPSolutionMethod == 'G') THEN

      TempR(1:nfco) = Glob_diagH(fb:fe)
      DO i = fb, fbn-1
        Glob_diagH(i) = Glob_diagH(nfco+i)
      ENDDO
      Glob_diagH(fbn:Glob_CurrBasisSize) = TempR(1:nfco)

    ELSE

      DO i = fb, fe
        TempR(i-fb+1) = Glob_H(i, i)
      ENDDO
      DO i = fb, fbn-1
        Glob_H(i, i) = Glob_H(nfco+i, nfco+i)
      ENDDO
      DO i = fbn, Glob_CurrBasisSize
        Glob_H(i, i) = TempR(i-fbn+1)
      ENDDO

    ENDIF

    TempR(1:nfco) = Glob_diagS(fb:fe)
    DO i = fb, fbn-1
      Glob_diagS(i) = Glob_diagS(nfco+i)
    ENDDO
    Glob_diagS(fbn:Glob_CurrBasisSize) = TempR(1:nfco)


    !==================================================================
    ! Off-diagonal elements
    !==================================================================
    ! The six stages described in the header, applied first to H and
    ! then, identically, to S.
    !------------------------------------------------------------------
    ! Now we permute off-diagonal elements
    fep = fe+1
    k = Glob_CurrBasisSize-fe
    q = Glob_CurrBasisSize+fep
    p = fb+Glob_CurrBasisSize-fe-1
    s = Glob_CurrBasisSize+fb
    r = s-1

    ! -- H, stage 1: columns left of the block, plain shift ----------
    DO i = 1, fb-1
      TempR(1:nfco) = Glob_H(fb:fe, i)
      DO j = fb, fbn-1
        Glob_H(j, i) = Glob_H(j+nfco, i)
      ENDDO
      Glob_H(fbn:Glob_CurrBasisSize, i) = TempR(1:nfco)
    ENDDO

    ! -- H, stage 2: stage the tail x tail block into the upper ------
    DO i = fep, Glob_CurrBasisSize-1
      Glob_H(fep:q-i-1, q-i) = Glob_H(i+1:Glob_CurrBasisSize, i)
    ENDDO

    ! -- H, stage 3: stage the tail x block block, transposed --------
    DO i = fe+1, Glob_CurrBasisSize
      DO j = fb, fe
        Glob_H(j, i) = Glob_H(i, j)
      ENDDO
    ENDDO

    ! -- H, stage 4: the within-block triangle, moved in place -------
    !    DOWNWARDS, which is what keeps it from reading what it has
    !    already written - see the header.
    DO i = fe-1, fb, -1
      Glob_H(i+k+1:Glob_CurrBasisSize, i+k) = Glob_H(i+1:fe, i)
    ENDDO

    ! -- H, stage 5: drop the stage-3 data into its new home ---------
    Glob_H(fbn:Glob_CurrBasisSize, fb:fb+k-1) = Glob_H(fb:fe, fe+1:Glob_CurrBasisSize)

    ! -- H, stage 6: recover the stage-2 data ------------------------
    DO i = fb, p-1
      Glob_H(i+1:p, i) = Glob_H(fep:r-i, s-i)
    ENDDO

    ! -- S, stage 1 --------------------------------------------------
    DO i = 1, fb-1
      TempR(1:nfco) = Glob_S(fb:fe, i)
      DO j = fb, fbn-1
        Glob_S(j, i) = Glob_S(j+nfco, i)
      ENDDO
      Glob_S(fbn:Glob_CurrBasisSize, i) = TempR(1:nfco)
    ENDDO

    ! -- S, stage 2 --------------------------------------------------
    DO i = fep, Glob_CurrBasisSize-1
      Glob_S(fep:q-i-1, q-i) = Glob_S(i+1:Glob_CurrBasisSize, i)
    ENDDO

    ! -- S, stage 3 --------------------------------------------------
    DO i = fe+1, Glob_CurrBasisSize
      DO j = fb, fe
        Glob_S(j, i) = Glob_S(i, j)
      ENDDO
    ENDDO

    ! -- S, stage 4 (downwards, for the same reason as H) ------------
    DO i = fe-1, fb, -1
      Glob_S(i+k+1:Glob_CurrBasisSize, i+k) = Glob_S(i+1:fe, i)
    ENDDO

    ! -- S, stage 5 --------------------------------------------------
    Glob_S(fbn:Glob_CurrBasisSize, fb:fb+k-1) = Glob_S(fb:fe, fe+1:Glob_CurrBasisSize)

    ! -- S, stage 6 --------------------------------------------------
    DO i = fb, p-1
      Glob_S(i+1:p, i) = Glob_S(fep:r-i, s-i)
    ENDDO


    !==================================================================
    ! Rebuild the upper triangle of S for the 'I' method
    !==================================================================
    ! Inverse iteration reads Glob_S symmetrically, so the upper
    ! triangle - used as scratch above - has to be restored.
    !------------------------------------------------------------------
    IF (Glob_GSEPSolutionMethod == 'I') THEN
      ! In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      ! upper triangle of Glob_S
      DO i = fb, fe
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE PermuteMatrixElements


  SUBROUTINE PermuteMatrixElements2(fb1, fe1, fe2, TempR)
    !==================================================================
    ! Subroutine PermuteMatrixElements2
    !==================================================================
    ! Applies to Glob_H, Glob_S, Glob_diagH and Glob_diagS the block
    ! exchange PermuteFunctions2 applies to the functions (fb1..fe1 with
    ! fe1+1..fe2). Call the two together. Workspace: TempR (>= fe1-fb1+1).
    ! Same storage convention and staged scheme as PermuteMatrixElements,
    ! with t the size of the first block and k = fe2-fe1 that of the second,
    ! plus one stage for the rows BELOW fe2, since the affected region does
    ! not reach the end of the basis. The in-place stage
    !     DO i=fe1-1,fb1,-1
    !       Glob_H(i+k+1:fe2,i+k) = Glob_H(i+1:fe1,i)
    !     ENDDO
    ! must run backwards for the reason given there; ascending it failed
    ! whenever k <= t-2 (checked by simulation for N=4..11).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: fb1, fe1  ! first block
    INTEGER :: fe2       ! end of the second block, which
    ! begins at fe1+1
    REAL(wp) :: TempR(*)  ! scratch, at least fe1-fb1+1 long

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j        ! loop counters
    INTEGER :: fn          ! where the first block lands
    INTEGER :: t           ! size of the first block
    INTEGER :: fe1p        ! fe1+1, same as fb2
    INTEGER :: k           ! size of the second block
    INTEGER :: q, p, r, s  ! index offsets for the staged blocks
    INTEGER :: fb2         ! start of the second block


    t = fe1-fb1+1
    fn = fb1+fe2-fe1
    fb2 = fe1+1


    !==================================================================
    ! Diagonals of H and S
    !==================================================================
    ! First we do the diagonals of H and S
    IF (Glob_GSEPSolutionMethod == 'G') THEN

      TempR(1:t) = Glob_diagH(fb1:fe1)
      DO i = fb1, fn-1
        Glob_diagH(i) = Glob_diagH(t+i)
      ENDDO
      Glob_diagH(fn:fe2) = TempR(1:t)

    ELSE

      DO i = fb1, fe1
        TempR(i-fb1+1) = Glob_H(i, i)
      ENDDO
      DO i = fb1, fn-1
        Glob_H(i, i) = Glob_H(t+i, t+i)
      ENDDO
      DO i = fn, fe2
        Glob_H(i, i) = TempR(i-fn+1)
      ENDDO

    ENDIF

    TempR(1:t) = Glob_diagS(fb1:fe1)
    DO i = fb1, fn-1
      Glob_diagS(i) = Glob_diagS(t+i)
    ENDDO
    Glob_diagS(fn:fe2) = TempR(1:t)


    !==================================================================
    ! Off-diagonal elements
    !==================================================================
    ! Applied first to H and then, identically, to S.
    !------------------------------------------------------------------
    ! Now we permute off-diagonal elements
    fe1p = fe1+1
    k = fe2-fe1
    q = fe2+fe1p
    p = fb1+fe2-fe1-1
    s = fe2+fb1
    r = s-1

    ! -- H, columns left of the blocks, plain shift ------------------
    DO i = 1, fb1-1
      TempR(1:t) = Glob_H(fb1:fe1, i)
      DO j = fb1, fn-1
        Glob_H(j, i) = Glob_H(t+j, i)
      ENDDO
      Glob_H(fn:fe2, i) = TempR(1:t)
    ENDDO

    ! -- H, stage the second block into the upper triangle -----------
    DO i = fe1p, fe2-1
      Glob_H(fe1p:q-i-1, q-i) = Glob_H(i+1:fe2, i)
    ENDDO

    ! -- H, stage the cross block, transposed ------------------------
    DO i = fe1+1, fe2
      DO j = fb1, fe1
        Glob_H(j, i) = Glob_H(i, j)
      ENDDO
    ENDDO

    ! -- H, the within-block triangle, moved in place ----------------
    !    DOWNWARDS, which is what keeps it from reading what it has
    !    already written - see the header.
    DO i = fe1-1, fb1, -1
      Glob_H(i+k+1:fe2, i+k) = Glob_H(i+1:fe1, i)
    ENDDO

    ! -- H, drop the staged cross block into its new home ------------
    Glob_H(fn:fe2, fb1:fb1+k-1) = Glob_H(fb1:fe1, fe1+1:fe2)

    ! -- H, recover the staged second block --------------------------
    DO i = fb1, p-1
      Glob_H(i+1:p, i) = Glob_H(fe1p:r-i, s-i)
    ENDDO

    ! -- H, rows below fe2: they see the swap as a column shift ------
    DO i = fe2+1, Glob_CurrBasisSize
      TempR(1:t) = Glob_H(i, fb1:fe1)
      DO j = fb1, fn-1
        Glob_H(i, j) = Glob_H(i, t+j)
      ENDDO
      Glob_H(i, fn:fe2) = TempR(1:t)
    ENDDO

    ! -- S, columns left of the blocks -------------------------------
    DO i = 1, fb1-1
      TempR(1:t) = Glob_S(fb1:fe1, i)
      DO j = fb1, fn-1
        Glob_S(j, i) = Glob_S(t+j, i)
      ENDDO
      Glob_S(fn:fe2, i) = TempR(1:t)
    ENDDO

    ! -- S, stage the second block -----------------------------------
    DO i = fe1p, fe2-1
      Glob_S(fe1p:q-i-1, q-i) = Glob_S(i+1:fe2, i)
    ENDDO

    ! -- S, stage the cross block, transposed ------------------------
    DO i = fe1+1, fe2
      DO j = fb1, fe1
        Glob_S(j, i) = Glob_S(i, j)
      ENDDO
    ENDDO

    ! -- S, the within-block triangle (downwards, as for H) ----------
    DO i = fe1-1, fb1, -1
      Glob_S(i+k+1:fe2, i+k) = Glob_S(i+1:fe1, i)
    ENDDO

    ! -- S, drop the staged cross block ------------------------------
    Glob_S(fn:fe2, fb1:fb1+k-1) = Glob_S(fb1:fe1, fe1+1:fe2)

    ! -- S, recover the staged second block --------------------------
    DO i = fb1, p-1
      Glob_S(i+1:p, i) = Glob_S(fe1p:r-i, s-i)
    ENDDO

    ! -- S, rows below fe2 -------------------------------------------
    DO i = fe2+1, Glob_CurrBasisSize
      TempR(1:t) = Glob_S(i, fb1:fe1)
      DO j = fb1, fn-1
        Glob_S(i, j) = Glob_S(i, t+j)
      ENDDO
      Glob_S(i, fn:fe2) = TempR(1:t)
    ENDDO


    !==================================================================
    ! Rebuild the upper triangle of S for the 'I' method
    !==================================================================
    ! Inverse iteration reads Glob_S symmetrically, so the upper
    ! triangle - used as scratch above - has to be restored.
    !------------------------------------------------------------------
    IF (Glob_GSEPSolutionMethod == 'I') THEN
      ! In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      ! upper triangle of Glob_S
      DO i = fb1, fe2
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE PermuteMatrixElements2


  SUBROUTINE ReverseFuncOrder(fb, fe)
    !==================================================================
    ! Subroutine ReverseFuncOrder
    !==================================================================
    ! Reverses the order of basis functions fb..fe (Glob_NonlinParam,
    ! Glob_PWR and Glob_FuncNum together); ReverseMatElemOrder does the
    ! matrix elements. f = (fe-fb+1)/2 (integer division) is the number of
    ! pairs to swap - an odd block leaves its middle function in place -
    ! and fbm = fb-1, fep = fe+1 address the i-th pair from both ends. temp
    ! is a local array of n*(n+1) elements, of which Glob_npt are used.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: fb, fe  ! range of functions to reverse

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    REAL(wp) :: temp(Glob_n*(Glob_n+1))

    INTEGER :: i         ! loop counter over the pairs
    INTEGER :: j         ! (declared, currently unused)
    INTEGER :: f         ! number of pairs to swap
    INTEGER :: t         ! integer swap temporary
    INTEGER :: fbm, fep  ! fb-1 and fe+1, the pair addresses


    f = (fe-fb+1)/2  ! integer division!
    fbm = fb-1
    fep = fe+1


    !==================================================================
    ! Nonlinear parameters
    !==================================================================
    DO i = 1, f
      temp(1:Glob_npt) = Glob_NonlinParam(1:Glob_npt, fbm+i)
      Glob_NonlinParam(1:Glob_npt, fbm+i) = Glob_NonlinParam(1:Glob_npt, fep-i)
      Glob_NonlinParam(1:Glob_npt, fep-i) = temp(1:Glob_npt)
    ENDDO


    !==================================================================
    ! Premultiplier powers
    !==================================================================
    DO i = 1, f
      t = Glob_PWR(fbm+i)
      Glob_PWR(fbm+i) = Glob_PWR(fep-i)
      Glob_PWR(fep-i) = t
    ENDDO


    !==================================================================
    ! Function numbers
    !==================================================================
    DO i = 1, f
      t = Glob_FuncNum(fbm+i)
      Glob_FuncNum(fbm+i) = Glob_FuncNum(fep-i)
      Glob_FuncNum(fep-i) = t
    ENDDO


  END SUBROUTINE ReverseFuncOrder


  SUBROUTINE ReverseMatElemOrder(fb, fe)
    !==================================================================
    ! Subroutine ReverseMatElemOrder
    !==================================================================
    ! Applies to Glob_H, Glob_S, Glob_diagH and Glob_diagS the reversal
    ! ReverseFuncOrder applies to functions fb..fe; call the two together.
    ! Only the LOWER triangles are used and no staging is needed: a
    ! reversal is a set of independent pair swaps done in place. The H
    ! diagonal is in Glob_diagH for 'G' and in Glob_H for 'I'.
    ! f = (fe-fb+1)/2 pairs; fbm = fb-1, fep = fe+1; ff = fe+fb is the
    ! reflection constant (index x maps to ff-x). The off-diagonal work
    ! splits into columns left of the block (swap rows), the block itself
    ! (reflect both indices) and rows below the block (swap columns). For
    ! 'I' the upper triangle of Glob_S is filled in at the end. Verified by
    ! simulation for N = 3..10.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: fb, fe  ! range of functions to reverse

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j      ! loop counters
    INTEGER :: f         ! number of pairs to swap
    INTEGER :: fbm, fep  ! fb-1 and fe+1, the pair addresses
    INTEGER :: ff        ! fe+fb, the reflection constant

    REAL(wp)    :: r  ! swap temporary for the diagonals
    COMPLEX(wp) :: c  ! swap temporary for the off-diagonals


    f = (fe-fb+1)/2  ! integer division!
    fbm = fb-1
    fep = fe+1
    ff = fe+fb


    !==================================================================
    ! Diagonal elements
    !==================================================================
    ! 'G' keeps the H diagonal in Glob_diagH; 'I' keeps it inside
    ! Glob_H. Glob_diagS is swapped in both cases.
    !------------------------------------------------------------------
    ! Diagonal elements
    IF (Glob_GSEPSolutionMethod == 'G') THEN

      DO i = 1, f
        r = Glob_diagH(fbm+i)
        Glob_diagH(fbm+i) = Glob_diagH(fep-i)
        Glob_diagH(fep-i) = r
      ENDDO

    ELSE

      DO i = 1, f
        r = Glob_H(fbm+i, fbm+i)
        Glob_H(fbm+i, fbm+i) = Glob_H(fep-i, fep-i)
        Glob_H(fep-i, fep-i) = r
      ENDDO

    ENDIF

    DO i = 1, f
      r = Glob_diagS(fbm+i)
      Glob_diagS(fbm+i) = Glob_diagS(fep-i)
      Glob_diagS(fep-i) = r
    ENDDO


    !==================================================================
    ! Off-diagonal elements
    !==================================================================
    ! The three regions described in the header, applied first to H and
    ! then, identically, to S.
    !------------------------------------------------------------------
    ! Off-diagonal elements

    ! -- H, columns left of the block: swap rows within the block ----
    DO i = 1, fbm
      DO j = 1, f
        c = Glob_H(fbm+j, i)
        Glob_H(fbm+j, i) = Glob_H(fep-j, i)
        Glob_H(fep-j, i) = c
      ENDDO
    ENDDO

    ! -- H, inside the block: reflect both indices -------------------
    DO i = 1, f
      DO j = fb+i, fe-i
        c = Glob_H(j, fbm+i)
        Glob_H(j, fbm+i) = Glob_H(fep-i, ff-j)
        Glob_H(fep-i, ff-j) = c
      ENDDO
      Glob_H(j, fbm+i) = Glob_H(j, fbm+i)
    ENDDO

    ! -- H, rows below the block: swap columns within the block ------
    DO i = fb, fb+f-1
      DO j = fep, Glob_CurrBasisSize
        c = Glob_H(j, i)
        Glob_H(j, i) = Glob_H(j, ff-i)
        Glob_H(j, ff-i) = c
      ENDDO
    ENDDO

    ! -- S, columns left of the block --------------------------------
    DO i = 1, fbm
      DO j = 1, f
        c = Glob_S(fbm+j, i)
        Glob_S(fbm+j, i) = Glob_S(fep-j, i)
        Glob_S(fep-j, i) = c
      ENDDO
    ENDDO

    ! -- S, inside the block -----------------------------------------
    DO i = 1, f
      DO j = fb+i, fe-i
        c = Glob_S(j, fbm+i)
        Glob_S(j, fbm+i) = Glob_S(fep-i, ff-j)
        Glob_S(fep-i, ff-j) = c
      ENDDO
      Glob_S(j, fbm+i) = Glob_S(j, fbm+i)
    ENDDO

    ! -- S, rows below the block -------------------------------------
    DO i = fb, fb+f-1
      DO j = fep, Glob_CurrBasisSize
        c = Glob_S(j, i)
        Glob_S(j, i) = Glob_S(j, ff-i)
        Glob_S(j, ff-i) = c
      ENDDO
    ENDDO


    !==================================================================
    ! Rebuild the upper triangle of S for the 'I' method
    !==================================================================
    ! Inverse iteration reads Glob_S symmetrically, so the upper
    ! triangle has to be filled in to match the new lower triangle.
    !------------------------------------------------------------------
    IF (Glob_GSEPSolutionMethod == 'I') THEN
      ! In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      ! upper triangle of Glob_S
      DO i = fb, fe
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE ReverseMatElemOrder


  SUBROUTINE SortBasisFuncAndMatElem(fb, fe, FuncNumTemp, NonlinParamTemp, TempR)
    !==================================================================
    ! Subroutine SortBasisFuncAndMatElem
    !==================================================================
    ! Sorts basis functions fb..fe into DECREASING order of Glob_FuncNum
    ! and permutes Glob_NonlinParam, Glob_PWR, Glob_FuncNum AND the matrix
    ! elements (Glob_H, Glob_S, Glob_diagH, Glob_diagS) in one call.
    !*** PRECONDITION (unchecked): Glob_FuncNum(fb:fe) has NO GAPS. A
    ! function numbered F goes to position k-F with k = fb-1+nfco+mf
    ! (mf = smallest number, nfco = fe-fb+1), so a gap writes out of range.
    ! Workspace (caller): NonlinParamTemp (>= (fe-fb+1)*Glob_npt), TempR
    ! (>= fe-fb+1); FuncNumTemp is in the argument list but NOT REFERENCED
    ! (Glob_PWR is round-tripped through the REAL array TempR, exact for
    ! the allowed powers). Only the LOWER triangles are meaningful; every
    ! off-diagonal region is first written into the UPPER triangle and then
    ! mirrored down, so nothing is read after being overwritten. For 'I'
    ! the upper triangle of Glob_S is rebuilt at the end.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: fb, fe                        ! range to sort
    INTEGER  :: FuncNumTemp(*)                ! not referenced - see header
    REAL(wp) :: NonlinParamTemp(Glob_npt, *)  ! parameter scratch
    REAL(wp) :: TempR(*)                      ! scratch, also used for Glob_PWR

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j    ! loop counters
    INTEGER :: k       ! destination offset: pos = k-FuncNum
    INTEGER :: mf      ! smallest function number in the range
    INTEGER :: fep     ! fe+1, first function past the range
    INTEGER :: nfco    ! number of functions in the range
    INTEGER :: fbm     ! fb-1, last function before the range
    INTEGER :: cbs     ! Glob_CurrBasisSize
    INTEGER :: fi, fj  ! function numbers of i and j


    cbs = Glob_CurrBasisSize
    nfco = fe-fb+1
    fep = fe+1
    fbm = fb-1
    mf = MINVAL(Glob_FuncNum(fb:fe))
    k = fbm+nfco+mf


    !==================================================================
    ! Nonlinear parameters and premultiplier powers
    !==================================================================
    ! Gather into the workspace at the destination slot, then copy the
    ! whole block back.
    !------------------------------------------------------------------
    ! First we sort out nonlinear parameters and Z-indices
    DO i = fb, fe
      NonlinParamTemp(1:Glob_npt, nfco+mf-Glob_FuncNum(i)) = Glob_NonlinParam(1:Glob_npt, i)
    ENDDO
    Glob_NonlinParam(1:Glob_npt, fb:fe) = NonlinParamTemp(1:Glob_npt, 1:nfco)

    DO i = fb, fe
      TempR(nfco+mf-Glob_FuncNum(i)) = Glob_PWR(i)
    ENDDO
    Glob_PWR(fb:fe) = TempR(1:nfco)


    !==================================================================
    ! Diagonal matrix elements
    !==================================================================
    ! 'G' keeps the H diagonal in Glob_diagH; 'I' keeps it inside
    ! Glob_H. Glob_diagS is sorted either way.
    !------------------------------------------------------------------
    ! Then we sort out the diagonal matrix elements
    IF (Glob_GSEPSolutionMethod == 'G') THEN

      DO i = fb, fe
        TempR(nfco+mf-Glob_FuncNum(i)) = Glob_diagH(i)
      ENDDO
      Glob_diagH(fb:fe) = TempR(1:nfco)

    ELSE

      DO i = fb, fe
        TempR(nfco+mf-Glob_FuncNum(i)) = Glob_H(i, i)
      ENDDO
      DO i = fb, fe
        Glob_H(i, i) = TempR(i-fb+1)
      ENDDO

    ENDIF

    DO i = fb, fe
      TempR(nfco+mf-Glob_FuncNum(i)) = Glob_diagS(i)
    ENDDO
    Glob_diagS(fb:fe) = TempR(1:nfco)


    !==================================================================
    ! Off-diagonal matrix elements
    !==================================================================
    ! Three regions, each staged through the upper triangle and then
    ! mirrored back down. Applied first to H, then identically to S.
    !------------------------------------------------------------------
    ! Sorting out off-diagonal matrix elements

    ! -- H, columns left of the range --------------------------------
    DO i = fb, fe
      Glob_H(1:fbm, k-Glob_FuncNum(i)) = Glob_H(i, 1:fbm)
    ENDDO
    DO i = fb, fe
      Glob_H(i, 1:fbm) = Glob_H(1:fbm, i)
    ENDDO

    ! -- H, rows below the range -------------------------------------
    DO i = fb, fe
      Glob_H(k-Glob_FuncNum(i), fep:cbs) = Glob_H(fep:cbs, i)
    ENDDO
    DO i = fb, fe
      Glob_H(fep:cbs, i) = Glob_H(i, fep:cbs)
    ENDDO

    ! -- H, within the range -----------------------------------------
    !    Both branches land in the UPPER triangle; the test just orders
    !    the destination pair so the smaller index becomes the row.
    DO i = fb, fe-1
      DO j = i+1, fe
        fi = Glob_FuncNum(i)
        fj = Glob_FuncNum(j)
        IF (fi < fj) THEN
          Glob_H(k-fj, k-fi) = Glob_H(j, i)
        ELSE
          Glob_H(k-fi, k-fj) = Glob_H(j, i)
        ENDIF
      ENDDO
    ENDDO
    DO i = fb, fe-1
      Glob_H(i+1:fe, i) = Glob_H(i, i+1:fe)
    ENDDO

    ! -- S, columns left of the range --------------------------------
    DO i = fb, fe
      Glob_S(1:fbm, k-Glob_FuncNum(i)) = Glob_S(i, 1:fbm)
    ENDDO
    DO i = fb, fe
      Glob_S(i, 1:fbm) = Glob_S(1:fbm, i)
    ENDDO

    ! -- S, rows below the range -------------------------------------
    DO i = fb, fe
      Glob_S(k-Glob_FuncNum(i), fep:cbs) = Glob_S(fep:cbs, i)
    ENDDO
    DO i = fb, fe
      Glob_S(fep:cbs, i) = Glob_S(i, fep:cbs)
    ENDDO

    ! -- S, within the range -----------------------------------------
    DO i = fb, fe-1
      DO j = i+1, fe
        fi = Glob_FuncNum(i)
        fj = Glob_FuncNum(j)
        IF (fi < fj) THEN
          Glob_S(k-fj, k-fi) = Glob_S(j, i)
        ELSE
          Glob_S(k-fi, k-fj) = Glob_S(j, i)
        ENDIF
      ENDDO
    ENDDO
    DO i = fb, fe-1
      Glob_S(i+1:fe, i) = Glob_S(i, i+1:fe)
    ENDDO


    !==================================================================
    ! Rebuild the upper triangle of S for the 'I' method
    !==================================================================
    ! Inverse iteration reads Glob_S symmetrically, so the upper
    ! triangle - used as scratch above - has to be restored.
    !------------------------------------------------------------------
    IF (Glob_GSEPSolutionMethod == 'I') THEN
      ! In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      ! upper triangle of Glob_S
      DO i = fb, fe
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
      ENDDO
    ENDIF


    !==================================================================
    ! Renumber the sorted functions
    !==================================================================
    ! Done LAST, because every mapping above reads the OLD numbering.
    ! After this the range carries mf+fe-i, i.e. decreasing from fb.
    !------------------------------------------------------------------
    ! At last we change the order of basis functions
    DO i = fb, fe
      Glob_FuncNum(i) = mf+fe-i
    ENDDO


  END SUBROUTINE SortBasisFuncAndMatElem


  SUBROUTINE GetOverlapStatistics(Nmin, Nmax, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
    !==================================================================
    ! Subroutine GetOverlapStatistics
    !==================================================================
    ! Reports how close to linearly dependent the basis is: the largest,
    ! smallest and mean magnitude of the pair overlaps that involve
    ! functions Nmin..Nmax (pairs entirely below Nmin cannot have changed).
    ! MaxAbsOverlap and MinAbsOverlap keep their SIGN (the magnitudes drive
    ! the comparisons through absMaxAbsOverlap/absMinAbsOverlap); only
    ! AverageAbsOverlap is a magnitude, with C(Nmax,2)-C(Nmin-1,2) pairs in
    ! the denominator. Only the LOWER triangle of Glob_S is read. No MPI:
    ! every rank holds the same Glob_S.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER  :: Nmin, Nmax         ! range of functions to look at
    REAL(wp) :: MaxAbsOverlap      ! out: signed, largest in magnitude
    REAL(wp) :: MinAbsOverlap      ! out: signed, smallest in magnitude
    REAL(wp) :: AverageAbsOverlap  ! out: mean magnitude

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j              ! loop counters
    INTEGER :: k                 ! (declared, currently unused)
    INTEGER :: nbands, leftover  ! (declared, currently unused)
    LOGICAL :: oddband           ! (declared, currently unused)

    REAL(wp) :: absMaxAbsOverlap  ! running largest magnitude
    REAL(wp) :: absMinAbsOverlap  ! running smallest magnitude
    REAL(wp) :: absSji            ! magnitude of the current element
    REAL(wp) :: t                 ! running sum of magnitudes


    MaxAbsOverlap = ZERO
    absMaxAbsOverlap = ZERO
    MinAbsOverlap = HUGE(MinAbsOverlap)
    absMinAbsOverlap = MinAbsOverlap
    t = ZERO


    !==================================================================
    ! Pairs with one function below Nmin and one inside the range
    !==================================================================
    DO i = 1, Nmin-1
      DO j = Nmin, Nmax

        absSji = ABS(Glob_S(j, i))

        IF (absSji > absMaxAbsOverlap) THEN
          absMaxAbsOverlap = absSji
          MaxAbsOverlap = Glob_S(j, i)
        ENDIF

        IF (absSji < absMinAbsOverlap) THEN
          absMinAbsOverlap = absSji
          MinAbsOverlap = Glob_S(j, i)
        ENDIF

        t = t+absSji

      ENDDO
    ENDDO


    !==================================================================
    ! Pairs with both functions inside the range
    !==================================================================
    DO i = Nmin, Nmax
      DO j = i+1, Nmax

        absSji = ABS(Glob_S(j, i))

        IF (absSji > absMaxAbsOverlap) THEN
          absMaxAbsOverlap = absSji
          MaxAbsOverlap = Glob_S(j, i)
        ENDIF

        IF (absSji < absMinAbsOverlap) THEN
          absMinAbsOverlap = absSji
          MinAbsOverlap = Glob_S(j, i)
        ENDIF

        t = t+absSji

      ENDDO
    ENDDO


    !==================================================================
    ! The mean magnitude
    !==================================================================
    ! The denominator counts the pairs visited, so it is zero exactly
    ! when neither loop ran - a range holding a single function, such as
    ! Nmin=Nmax=1. Without this guard that is a 0/0, and MinAbsOverlap
    ! would additionally come back as the HUGE() it was primed with. All
    ! three statistics are undefined in that case, so report zeros.
    !
    ! Returning early is safe: this routine performs no collectives.
    !------------------------------------------------------------------
    IF (Nmax*(Nmax-1)-(Nmin-1)*(Nmin-2) == 0) THEN
      MaxAbsOverlap = ZERO
      MinAbsOverlap = ZERO
      AverageAbsOverlap = ZERO
      RETURN
    ENDIF

    AverageAbsOverlap = 2*t/(Nmax*(Nmax-1)-(Nmin-1)*(Nmin-2))


  END SUBROUTINE GetOverlapStatistics

  FUNCTION NumOfRowsToPermForUnitShift(j)
    !==================================================================
    ! Function NumOfRowsToPermForUnitShift
    !==================================================================
    ! Returns the largest power of two dividing j (the ruler sequence
    ! 1,2,1,4,1,2,1,8,...). The cyclic optimization uses it to decide how
    ! many functions to permute when the window shifts by one step, so the
    ! cost of keeping the basis ordered is spread over the steps. k is
    ! doubled while it divides j, then halved once.
    !*** PRECONDITION: j >= 1; with j=0 the loop never ends. Both call
    ! sites pass OptIterCounter-1 inside IF (OptIterCounter/=1), so the
    ! argument is always at least 1.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments and result
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: NumOfRowsToPermForUnitShift  ! function result
    INTEGER :: j                            ! step number, >= 1

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: k  ! running power of two


    k = 1

    DO WHILE (MOD(j, k) == 0)
      k = k*2
    ENDDO

    NumOfRowsToPermForUnitShift = k/2


  END FUNCTION NumOfRowsToPermForUnitShift


  SUBROUTINE GenerateRndIntSeq(n, s)
    !==================================================================
    ! Subroutine GenerateRndIntSeq
    !==================================================================
    ! Returns a random permutation of 1..n in s(1:n), used to visit the
    ! premultiplier powers of a window in random order. s starts as the
    ! identity and is shuffled with DO i=1,n: j = random in 1..i, swap
    ! s(i), s(j) - a valid Fisher-Yates shuffle (every permutation equally
    ! likely; the i=1 pass is a no-op). REAL(8) is used for the draw like
    ! the other random numbers in this module.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: n     ! length of the sequence
    INTEGER :: s(n)  ! out: a permutation of 1..n

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j  ! loop counter and draw
    INTEGER :: k     ! swap temporary

    REAL(8) :: r8  ! the random draw


    IF (n == 1) THEN

      s(1) = 1

    ELSE

      ! Start from the identity
      DO i = 1, n
        s(i) = i
      ENDDO

      ! Shuffle it
      DO i = 1, n
        CALL RANDOM_NUMBER(r8)
        j = INT(r8*i)+1
        k = s(i); s(i) = s(j); s(j) = k
      ENDDO

    ENDIF


  END SUBROUTINE GenerateRndIntSeq


  SUBROUTINE ReallocateBasisFuncData(FinalSize, NumOfFuncToKeep)
    !==================================================================
    ! Subroutine ReallocateBasisFuncData
    !==================================================================
    ! Resizes Glob_History, Glob_FuncNum, Glob_PWR and Glob_NonlinParam to
    ! FinalSize, keeping the first NumOfFuncToKeep functions and clearing
    ! the rest; called at the top of BasisEnlG/BasisEnlI to grow the basis
    ! to Kstop. NumOfFuncToKeep <= FinalSize is checked (EC0120).
    ! Fortran 90 cannot grow an allocatable array in place, so the data is
    ! parked, the arrays are reallocated and the data is copied back;
    ! Glob_History is copied whole so NumOfEnergyEvalDuringFullOpt survives.
    ! Glob_UseReallocFile=.FALSE. (the branch that runs) parks it in memory;
    ! .TRUE. parks it in Glob_ReallocFileName via rank 0. *** THAT BRANCH IS
    ! UNREACHABLE (nothing sets the flag) AND DEFECTIVE: its MPI_BCAST of
    ! OpenFileErr sits inside IF (Glob_ProcID==0). Left as is on purpose.
    ! The new slots get Glob_FuncNum=0 until the enlargement loop numbers
    ! them.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: FinalSize        ! new size of the arrays
    INTEGER :: NumOfFuncToKeep  ! how many functions to carry over

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i            ! loop counter
    INTEGER :: OpenFileErr  ! IOSTAT of the OPEN and of the last READ

    ! Buffers used to broadcast Glob_History field by field. An MPI
    ! broadcast needs a contiguous block of ONE datatype, and
    ! Glob_HistoryStep mixes a real with three integers, so each
    ! field is packed into a plain array, broadcast, and unpacked.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: WorkBuffReal
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: WorkBuffInt

    ! Temporaries that hold the kept data across the reallocation
    TYPE(Glob_HistoryStep), ALLOCATABLE, DIMENSION(:) :: TempHistory
    REAL(wp), ALLOCATABLE, DIMENSION(:, :)            :: TempParam
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempFunc
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempZInd


    !==================================================================
    ! Guard: the kept data has to fit in the new arrays
    !==================================================================
    ! Every rank evaluates the same condition on the same values, so
    ! either all of them abort or none does.
    !------------------------------------------------------------------
    IF (NumOfFuncToKeep > FinalSize) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0120 in ReallocateBasisFuncData:'
        WRITE(*, *) 'NumOfFuncToKeep must be smaller or equal than FinalSize'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    ENDIF


    IF (Glob_UseReallocFile) THEN

      !==================================================================
      ! Route 1: park the data in an external file  *** DEAD BRANCH ***
      !==================================================================
      ! See the header: Glob_UseReallocFile is never .TRUE., and the
      ! OpenFileErr broadcast below is misplaced. Left untouched.
      !------------------------------------------------------------------

      !--------------------------------------------------------------
      ! Rank 0 writes the kept data out
      !--------------------------------------------------------------
      ! Only rank 0 writes: every rank holds the same values, and they
      ! would otherwise all write the same file at the same time.
      ! Temporarily store the information in a file
      IF (Glob_ProcID == 0) THEN
        IF (NumOfFuncToKeep > 0) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a47)', ADVANCE='no') 'Reallocating some arrays using external file...'
          OPEN(1, FILE=Glob_ReallocFileName, FORM='unformatted', STATUS='replace')
          WRITE(1) Glob_History(1:NumOfFuncToKeep)
          WRITE(1) Glob_FuncNum(1:NumOfFuncToKeep)
          WRITE(1) Glob_PWR(1:NumOfFuncToKeep)
          WRITE(1) Glob_NonlinParam(1:Glob_npt, 1:NumOfFuncToKeep)
          CLOSE(1)
        ENDIF
      ENDIF

      !--------------------------------------------------------------
      ! Every rank resizes its own copy of the arrays
      !--------------------------------------------------------------
      ! Released in the reverse of the order they are allocated in.
      ! The contents are gone after this point on EVERY rank,
      ! including rank 0 - the file is the only copy left.
      DEALLOCATE(Glob_NonlinParam)
      DEALLOCATE(Glob_PWR)
      DEALLOCATE(Glob_FuncNum)
      DEALLOCATE(Glob_History)

      ALLOCATE(Glob_History(FinalSize))
      ALLOCATE(Glob_FuncNum(FinalSize))
      ALLOCATE(Glob_PWR(FinalSize))
      ALLOCATE(Glob_NonlinParam(Glob_npt, FinalSize))

      !--------------------------------------------------------------
      ! Rank 0 reads the data back
      !--------------------------------------------------------------
      ! Only the last READ carries IOSTAT, so a failure in one of the
      ! three earlier records is not caught here.
      IF (Glob_ProcID == 0) THEN
        IF (NumOfFuncToKeep > 0) THEN
          OPEN(1, FILE=Glob_ReallocFileName, FORM='unformatted', STATUS='old', IOSTAT=OpenFileErr)
          IF (OpenFileErr == 0) THEN
            READ(1) Glob_History(1:NumOfFuncToKeep)
            READ(1) Glob_FuncNum(1:NumOfFuncToKeep)
            READ(1) Glob_PWR(1:NumOfFuncToKeep)
            READ(1, IOSTAT=OpenFileErr) Glob_NonlinParam(1:Glob_npt, 1:NumOfFuncToKeep)
          ENDIF
          CLOSE(1)
          CALL MPI_BCAST(OpenFileErr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        ENDIF
      ENDIF

      !--------------------------------------------------------------
      ! Abort if the data could not be read back
      !--------------------------------------------------------------
      ! Unrecoverable: the arrays were already deallocated above, so
      ! there is no copy of the basis left in memory to fall back on.
      IF ((OpenFileErr /= 0) .AND. (NumOfFuncToKeep > 0)) THEN
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 2) WRITE(*, *)
          WRITE(*, *) 'Error EC0121 in ReallocateBasisFuncData:'
          WRITE(*, *) 'cannot read data from file', Glob_ReallocFileName
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
      ENDIF

      !--------------------------------------------------------------
      ! Blank the temporary file
      !--------------------------------------------------------------
      ! status='replace' truncates it, and a short string is written
      ! so that what is left on disk is small and obviously not live
      ! data.
      IF (Glob_ProcID == 0) THEN
        OPEN(1, FILE=Glob_ReallocFileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)
        WRITE(1) 'This temporary file is empty'
        CLOSE(1)
      ENDIF

      !--------------------------------------------------------------
      ! Distribute Glob_History, one field at a time
      !--------------------------------------------------------------
      ! Pack into a plain array, broadcast, unpack. Four fields, four
      ! broadcasts; the two buffers are allocated once and reused.
      ALLOCATE(WorkBuffReal(NumOfFuncToKeep))
      ALLOCATE(WorkBuffInt(NumOfFuncToKeep))

      ! -- Energy ----------------------------------------------------
      DO i = 1, NumOfFuncToKeep
        WorkBuffReal(i) = Glob_History(i)%Energy
      ENDDO
      CALL MPI_BCAST(WorkBuffReal, NumOfFuncToKeep, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      DO i = 1, NumOfFuncToKeep
        Glob_History(i)%Energy = WorkBuffReal(i)
      ENDDO

      ! -- CyclesDone ------------------------------------------------
      DO i = 1, NumOfFuncToKeep
        WorkBuffInt(i) = Glob_History(i)%CyclesDone
      ENDDO
      CALL MPI_BCAST(WorkBuffInt, NumOfFuncToKeep, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      DO i = 1, NumOfFuncToKeep
        Glob_History(i)%CyclesDone = WorkBuffInt(i)
      ENDDO

      ! -- InitFuncAtLastStep ----------------------------------------
      DO i = 1, NumOfFuncToKeep
        WorkBuffInt(i) = Glob_History(i)%InitFuncAtLastStep
      ENDDO
      CALL MPI_BCAST(WorkBuffInt, NumOfFuncToKeep, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      DO i = 1, NumOfFuncToKeep
        Glob_History(i)%InitFuncAtLastStep = WorkBuffInt(i)
      ENDDO

      ! -- NumOfEnergyEvalDuringFullOpt ------------------------------
      DO i = 1, NumOfFuncToKeep
        WorkBuffInt(i) = Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      ENDDO
      CALL MPI_BCAST(WorkBuffInt, NumOfFuncToKeep, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      DO i = 1, NumOfFuncToKeep
        Glob_History(i)%NumOfEnergyEvalDuringFullOpt = WorkBuffInt(i)
      ENDDO

      DEALLOCATE(WorkBuffReal)
      DEALLOCATE(WorkBuffInt)

      !--------------------------------------------------------------
      ! Distribute the remaining three arrays
      !--------------------------------------------------------------
      ! These are plain arrays of a single type, so each goes out in
      ! one broadcast. Glob_NonlinParam is contiguous and column
      ! major, so its leading Glob_npt*NumOfFuncToKeep elements are
      ! exactly the kept columns.
      CALL MPI_BCAST(Glob_FuncNum, NumOfFuncToKeep, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_PWR, NumOfFuncToKeep, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_NonlinParam, Glob_npt*NumOfFuncToKeep, &
                     MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      ! Closes the line opened with ADVANCE='no' above
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) 'done'
      ENDIF

    ELSE  ! if (Glob_UseReallocFile)

      !==================================================================
      ! Route 2: park the data in memory
      !==================================================================
      ! The branch that actually runs. Purely local: every rank does the
      ! same thing to its own copy, so there is nothing to communicate.
      !
      ! Peak memory is the old arrays plus the temporaries plus the new
      ! arrays. At the sizes this code is used at, the per-function data
      ! is tiny next to Glob_H and Glob_S, which are not touched here.
      !------------------------------------------------------------------

      !--------------------------------------------------------------
      ! Copy the kept data into temporaries
      !--------------------------------------------------------------
      ALLOCATE(TempHistory(NumOfFuncToKeep))
      ALLOCATE(TempFunc(NumOfFuncToKeep))
      ALLOCATE(TempZInd(NumOfFuncToKeep))
      ALLOCATE(TempParam(Glob_npt, NumOfFuncToKeep))

      TempHistory(1:NumOfFuncToKeep) = Glob_History(1:NumOfFuncToKeep)
      TempFunc(1:NumOfFuncToKeep) = Glob_FuncNum(1:NumOfFuncToKeep)
      TempZInd(1:NumOfFuncToKeep) = Glob_PWR(1:NumOfFuncToKeep)
      TempParam(1:Glob_npt, 1:NumOfFuncToKeep) = Glob_NonlinParam(1:Glob_npt, 1:NumOfFuncToKeep)

      !--------------------------------------------------------------
      ! Resize
      !--------------------------------------------------------------
      DEALLOCATE(Glob_NonlinParam)
      DEALLOCATE(Glob_PWR)
      DEALLOCATE(Glob_FuncNum)
      DEALLOCATE(Glob_History)

      ALLOCATE(Glob_History(FinalSize))
      ALLOCATE(Glob_FuncNum(FinalSize))
      ALLOCATE(Glob_PWR(FinalSize))
      ALLOCATE(Glob_NonlinParam(Glob_npt, FinalSize))

      !--------------------------------------------------------------
      ! Copy the kept data back and release the temporaries
      !--------------------------------------------------------------
      Glob_History(1:NumOfFuncToKeep) = TempHistory(1:NumOfFuncToKeep)
      Glob_FuncNum(1:NumOfFuncToKeep) = TempFunc(1:NumOfFuncToKeep)
      Glob_PWR(1:NumOfFuncToKeep) = TempZInd(1:NumOfFuncToKeep)
      Glob_NonlinParam(1:Glob_npt, 1:NumOfFuncToKeep) = TempParam(1:Glob_npt, 1:NumOfFuncToKeep)

      DEALLOCATE(TempParam)
      DEALLOCATE(TempZInd)
      DEALLOCATE(TempFunc)
      DEALLOCATE(TempHistory)

    ENDIF


    !==================================================================
    ! Clear the slots beyond the kept functions
    !==================================================================
    ! Common to both routes. ALLOCATE does not initialize, so without
    ! this the new slots would hold whatever happened to be on the
    ! heap. The loop is empty when NumOfFuncToKeep==FinalSize, i.e.
    ! when the arrays were resized to exactly the data they held.
    !
    ! Glob_FuncNum=0 here is a placeholder - see the header note.
    !------------------------------------------------------------------
    DO i = NumOfFuncToKeep+1, FinalSize
      Glob_History(i)%Energy = ZERO
      Glob_History(i)%CyclesDone = 0
      Glob_History(i)%InitFuncAtLastStep = 0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = 0
      Glob_FuncNum(i) = 0
      Glob_PWR(i) = 0
      Glob_NonlinParam(1:Glob_npt, i) = ZERO
    ENDDO


  END SUBROUTINE ReallocateBasisFuncData


  !==================================================================
  ! The QR solution method ('Q'): workspace and transaction helpers,
  ! the solve, the overlap penalty and statistics, the energy routines
  ! EnergyQA/EnergyQAM/EnergyQB, the common solver of the elimination
  ! routines and the Q cleanup, then the drivers BasisEnlQ, OptCycleQ
  ! and FullOpt1Q, the Q twins of the G and I drivers below.
  !==================================================================

  SUBROUTINE PrepareQWorkspace(MatrixOrder, Capacity, MaxActive, ErrorCode)
    ! Subroutine PrepareQWorkspace allocates the workproc-owned storage shared by
    ! the Q energy routines and one Q BBOP driver. It does not allocate Glob_H or
    ! Glob_S and does not initialize qrlinalg. The BBOP driver owns those lifetimes
    ! and must call this routine only after Glob_npt has been initialized.
    !
    ! The physical matrices can have capacity greater than MatrixOrder during basis
    ! enlargement. Only indices 1 through MatrixOrder belong to the represented
    ! problem. MaxActive is normally Kstep for BASIS_ENL, NumOfFuncToOpt for
    ! OPT_CYCLE, and FinalFunc-InitFunc+1 for FULL_OPT1.
    !
    ! Arguments:
    INTEGER, INTENT(IN)  :: MatrixOrder, Capacity, MaxActive
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: AllocationStatus

    CALL ClearQWorkspace()
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (MatrixOrder < 0) RETURN
    IF (Capacity < MAX(1, MatrixOrder)) RETURN
    IF (MaxActive < 1) RETURN
    IF (Glob_npt < 1) RETURN

    ALLOCATE(Q_Workspace%ActiveFunction(MaxActive), &
             Q_Workspace%ActivePosition(Capacity), &
             Q_Workspace%PreviousH(Capacity, MaxActive), &
             Q_Workspace%PreviousS(Capacity, MaxActive), &
             Q_Workspace%TrialH(Capacity, MaxActive), &
             Q_Workspace%TrialS(Capacity, MaxActive), &
             Q_Workspace%PreviousDiagS(MaxActive), &
             Q_Workspace%TrialDiagS(MaxActive), &
             Q_Workspace%MatrixParam(Glob_npt, MaxActive), &
             Q_Workspace%PreviousParam(Glob_npt, MaxActive), &
             Q_Workspace%AcceptedParam(Glob_npt, MaxActive), &
             Q_Workspace%InitialVector(Capacity), &
             Q_Workspace%SolvedVector(Capacity), &
             Q_Workspace%DeltaH(Capacity), &
             Q_Workspace%DeltaS(Capacity), &
             STAT = AllocationStatus)
    IF (AllocationStatus /= 0) THEN
      CALL ClearQWorkspace()
      ErrorCode = Q_METHOD_ALLOCATION_ERROR
      RETURN
    ENDIF

    Q_Workspace%MatrixOrder = MatrixOrder
    Q_Workspace%Capacity = Capacity
    Q_Workspace%MaxActive = MaxActive
    Q_Workspace%NumActive = 0
    Q_Workspace%ActiveFunction = 0
    Q_Workspace%ActivePosition = 0
    Q_Workspace%PreviousH = ZERO
    Q_Workspace%PreviousS = ZERO
    Q_Workspace%TrialH = ZERO
    Q_Workspace%TrialS = ZERO
    Q_Workspace%PreviousDiagS = ZERO
    Q_Workspace%TrialDiagS = ZERO
    Q_Workspace%MatrixParam = ZERO
    Q_Workspace%PreviousParam = ZERO
    Q_Workspace%AcceptedParam = ZERO
    Q_Workspace%InitialVector = ONE
    Q_Workspace%SolvedVector = ZERO
    Q_Workspace%DeltaH = ZERO
    Q_Workspace%DeltaS = ZERO
    Q_Workspace%AcceptedEnergy = ZERO
    Q_Workspace%AcceptedPointIsStored = .FALSE.
    Q_Workspace%LastEigenpairResidual = HUGE(ONE)
    Q_Workspace%LastFactorResidual = HUGE(ONE)
    Q_Workspace%FreshFactorizations = 0
    Q_Workspace%MatricesAreCanonical = .FALSE.
    Q_Workspace%FactorsMatchMatrices = .FALSE.
    Q_Workspace%MatrixParametersAreStored = .FALSE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE PrepareQWorkspace

  SUBROUTINE EnsureQActiveCapacity(RequiredActive, ErrorCode)
    ! Subroutine EnsureQActiveCapacity enlarges only the transaction part of an
    ! existing Q workspace. The QR factors and the full-order vector workspace are
    ! deliberately preserved. Cleanup routines initially reserve one active column
    ! because that is sufficient for elimination, but a separation test discovers
    ! the number of columns to replace only after the first eigenvector or overlap
    ! matrix has been inspected. Rebuilding the complete Q workspace at that point
    ! would discard the factorization whose reuse is the purpose of the Q method.
    !
    ! This operation is valid only between transactions. No previous, trial, or
    ! accepted point is copied: cleanup routines call it before selecting their
    ! first active set. Allocating every replacement array before move_alloc makes
    ! allocation failure transactional as well; the original workspace remains
    ! usable when any allocation fails.
    !
    ! Arguments:
    INTEGER, INTENT(IN)  :: RequiredActive
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER               :: AllocationStatus
    INTEGER, ALLOCATABLE  :: NewActiveFunction(:)
    REAL(wp), ALLOCATABLE :: NewPreviousH(:, :), NewPreviousS(:, :)
    REAL(wp), ALLOCATABLE :: NewTrialH(:, :), NewTrialS(:, :)
    REAL(wp), ALLOCATABLE :: NewPreviousDiagS(:), NewTrialDiagS(:)
    REAL(wp), ALLOCATABLE :: NewMatrixParam(:, :), NewPreviousParam(:, :)
    REAL(wp), ALLOCATABLE :: NewAcceptedParam(:, :)

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (RequiredActive < 1) RETURN
    IF (RequiredActive > Q_Workspace%Capacity) RETURN
    IF (Q_Workspace%MatrixOrder < 1) RETURN
    IF (Q_Workspace%NumActive /= 0) RETURN
    IF (Q_Workspace%MatrixParametersAreStored) RETURN
    IF (Q_Workspace%TrialIsReady) RETURN
    IF (Q_Workspace%AcceptedPointIsStored) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (RequiredActive <= Q_Workspace%MaxActive) THEN
      ErrorCode = Q_METHOD_SUCCESS
      RETURN
    ENDIF

    ALLOCATE(NewActiveFunction(RequiredActive), &
             NewPreviousH(Q_Workspace%Capacity, RequiredActive), &
             NewPreviousS(Q_Workspace%Capacity, RequiredActive), &
             NewTrialH(Q_Workspace%Capacity, RequiredActive), &
             NewTrialS(Q_Workspace%Capacity, RequiredActive), &
             NewPreviousDiagS(RequiredActive), &
             NewTrialDiagS(RequiredActive), &
             NewMatrixParam(Glob_npt, RequiredActive), &
             NewPreviousParam(Glob_npt, RequiredActive), &
             NewAcceptedParam(Glob_npt, RequiredActive), &
             STAT = AllocationStatus)
    IF (AllocationStatus /= 0) THEN
      ErrorCode = Q_METHOD_ALLOCATION_ERROR
      RETURN
    ENDIF

    NewActiveFunction = 0
    NewPreviousH = ZERO
    NewPreviousS = ZERO
    NewTrialH = ZERO
    NewTrialS = ZERO
    NewPreviousDiagS = ZERO
    NewTrialDiagS = ZERO
    NewMatrixParam = ZERO
    NewPreviousParam = ZERO
    NewAcceptedParam = ZERO

    CALL move_alloc(NewActiveFunction, Q_Workspace%ActiveFunction)
    CALL move_alloc(NewPreviousH, Q_Workspace%PreviousH)
    CALL move_alloc(NewPreviousS, Q_Workspace%PreviousS)
    CALL move_alloc(NewTrialH, Q_Workspace%TrialH)
    CALL move_alloc(NewTrialS, Q_Workspace%TrialS)
    CALL move_alloc(NewPreviousDiagS, Q_Workspace%PreviousDiagS)
    CALL move_alloc(NewTrialDiagS, Q_Workspace%TrialDiagS)
    CALL move_alloc(NewMatrixParam, Q_Workspace%MatrixParam)
    CALL move_alloc(NewPreviousParam, Q_Workspace%PreviousParam)
    CALL move_alloc(NewAcceptedParam, Q_Workspace%AcceptedParam)
    Q_Workspace%MaxActive = RequiredActive
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE EnsureQActiveCapacity

  SUBROUTINE ClearQWorkspace()
    ! Subroutine ClearQWorkspace releases workproc-owned Q storage. It intentionally
    ! does not deallocate the existing Glob_ arrays, because their lifetime belongs
    ! to the BBOP driver. The qrlinalg state is cleared here before the remaining
    ! metadata are reset.

    IF (ALLOCATED(Q_Workspace%DeltaS)) DEALLOCATE(Q_Workspace%DeltaS)
    IF (ALLOCATED(Q_Workspace%DeltaH)) DEALLOCATE(Q_Workspace%DeltaH)
    IF (ALLOCATED(Q_Workspace%SolvedVector)) DEALLOCATE(Q_Workspace%SolvedVector)
    IF (ALLOCATED(Q_Workspace%InitialVector)) DEALLOCATE(Q_Workspace%InitialVector)
    IF (ALLOCATED(Q_Workspace%AcceptedParam)) DEALLOCATE(Q_Workspace%AcceptedParam)
    IF (ALLOCATED(Q_Workspace%PreviousParam)) DEALLOCATE(Q_Workspace%PreviousParam)
    IF (ALLOCATED(Q_Workspace%MatrixParam)) DEALLOCATE(Q_Workspace%MatrixParam)
    IF (ALLOCATED(Q_Workspace%TrialDiagS)) DEALLOCATE(Q_Workspace%TrialDiagS)
    IF (ALLOCATED(Q_Workspace%PreviousDiagS)) DEALLOCATE(Q_Workspace%PreviousDiagS)
    IF (ALLOCATED(Q_Workspace%TrialS)) DEALLOCATE(Q_Workspace%TrialS)
    IF (ALLOCATED(Q_Workspace%TrialH)) DEALLOCATE(Q_Workspace%TrialH)
    IF (ALLOCATED(Q_Workspace%PreviousS)) DEALLOCATE(Q_Workspace%PreviousS)
    IF (ALLOCATED(Q_Workspace%PreviousH)) DEALLOCATE(Q_Workspace%PreviousH)
    IF (ALLOCATED(Q_Workspace%ActivePosition)) DEALLOCATE(Q_Workspace%ActivePosition)
    ! clear is valid for an initialized or empty state. Calling it on every MPI
    ! rank is therefore safe even though only rank zero will own allocated QR
    ! factors once FactorizeQFresh is implemented.
    CALL Q_Workspace%Factors%clear()

    IF (ALLOCATED(Q_Workspace%ActiveFunction)) DEALLOCATE(Q_Workspace%ActiveFunction)

    Q_Workspace%MatrixOrder = 0
    Q_Workspace%Capacity = 0
    Q_Workspace%MaxActive = 0
    Q_Workspace%NumActive = 0
    Q_Workspace%MatricesAreCanonical = .FALSE.
    Q_Workspace%FactorsMatchMatrices = .FALSE.
    Q_Workspace%MatrixParametersAreStored = .FALSE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    Q_Workspace%AcceptedEnergy = ZERO
    Q_Workspace%AcceptedPointIsStored = .FALSE.
    Q_Workspace%LastEigenpairResidual = HUGE(ONE)
    Q_Workspace%LastFactorResidual = HUGE(ONE)
    Q_Workspace%FreshFactorizations = 0

  END SUBROUTINE ClearQWorkspace

  SUBROUTINE SetQActiveFunctions(ActiveFunction, ErrorCode)
    ! Subroutine SetQActiveFunctions defines optimizer block order without changing
    ! the physical order of basis functions. ActiveFunction must contain distinct
    ! canonical indices in the range 1:Q_Workspace%MatrixOrder. The order supplied
    ! here is also the order of nonlinear parameter blocks, gradient blocks, and
    ! saved Hessian rows and columns.
    !
    ! The routine updates both maps only after the complete input has been checked,
    ! so an invalid selection leaves the previous active set unchanged.
    !
    ! Arguments:
    INTEGER, INTENT(IN)  :: ActiveFunction(:)
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: i, j, NumActive

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActivePosition)) RETURN
    NumActive = SIZE(ActiveFunction)
    IF ((NumActive < 1) .OR. (NumActive > Q_Workspace%MaxActive)) RETURN
    DO i = 1, NumActive
      IF ((ActiveFunction(i) < 1) .OR. &
          (ActiveFunction(i) > Q_Workspace%MatrixOrder)) RETURN
      DO j = 1, i-1
        IF (ActiveFunction(i) == ActiveFunction(j)) RETURN
      ENDDO
    ENDDO

    Q_Workspace%ActiveFunction = 0
    Q_Workspace%ActivePosition = 0
    Q_Workspace%ActiveFunction(1:NumActive) = ActiveFunction
    DO i = 1, NumActive
      Q_Workspace%ActivePosition(ActiveFunction(i)) = i
    ENDDO
    Q_Workspace%NumActive = NumActive
    Q_Workspace%AcceptedPointIsStored = .FALSE.
    Q_Workspace%MatrixParametersAreStored = .FALSE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE SetQActiveFunctions

  SUBROUTINE CaptureQMatrixParameters(ErrorCode)
    ! Subroutine CaptureQMatrixParameters records the nonlinear parameters that
    ! belong to the currently represented canonical H/S matrices. A Q driver calls
    ! this immediately after selecting a new active set, before it copies a DRMNG
    ! trial point into Glob_NonlinParam. The separate copy is essential: after that
    ! copy Glob_NonlinParam describes the requested trial, while Glob_H and Glob_S
    ! still describe the preceding point until ApplyQTrial commits successfully.
    !
    ! This routine cannot prove that matrix elements were calculated from the
    ! current parameters. MatricesAreCanonical is therefore an explicit caller
    ! precondition set only after a full Q assembly or a successful swap restore.
    !
    ! Arguments:
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: a, FunctionIndex

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (Q_Workspace%NumActive < 1) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%MatrixParam)) RETURN
    IF (.NOT. ALLOCATED(Glob_NonlinParam)) RETURN
    IF (SIZE(Glob_NonlinParam, 1) < Glob_npt) RETURN
    IF (SIZE(Glob_NonlinParam, 2) < Q_Workspace%MatrixOrder) RETURN

    DO a = 1, Q_Workspace%NumActive
      FunctionIndex = Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt, a) = &
        Glob_NonlinParam(1:Glob_npt, FunctionIndex)
    ENDDO
    IF (Q_Workspace%NumActive < Q_Workspace%MaxActive) THEN
      Q_Workspace%MatrixParam(1:Glob_npt, Q_Workspace%NumActive+1:Q_Workspace%MaxActive) = ZERO
    ENDIF
    Q_Workspace%MatrixParametersAreStored = .TRUE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE CaptureQMatrixParameters

  SUBROUTINE TrimQFactors(TargetOrder, ErrorCode)
    ! Subroutine TrimQFactors removes a rejected canonical suffix from the QR state.
    ! BASIS_ENL repeatedly tests new functions in positions TargetOrder+1 onward;
    ! deleting those last rows and columns restores the accepted prefix without any
    ! basis permutation or cubic refactorization. The physical suffix may retain
    ! obsolete trial values because it lies outside MatrixOrder and is overwritten
    ! before it can become active again.
    !
    ! qrlinalg intentionally has no valid order-zero factorization. When the first
    ! basis block is being selected, trimming to zero therefore reinitializes an
    ! empty capacity reservation. The next trial is constructed by a fresh
    ! factorization after its complete physical block has been staged.
    !
    ! Arguments:
    INTEGER, INTENT(IN)  :: TargetOrder
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: CurrentOrder, RootError, RecoveryError

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    CurrentOrder = Q_Workspace%MatrixOrder
    IF ((TargetOrder < 0) .OR. (TargetOrder > CurrentOrder)) RETURN
    IF (Q_Workspace%Capacity < MAX(1, CurrentOrder)) RETURN

    RootError = Q_METHOD_SUCCESS
    IF (Glob_ProcID == 0) THEN
      IF (TargetOrder == 0) THEN
        CALL Q_Workspace%Factors%clear()
        CALL Q_Workspace%Factors%initialize(Q_Workspace%Capacity, RootError)
      ELSE
        IF (.NOT. Q_Workspace%Factors%is_valid()) THEN
          RootError = Q_METHOD_INVALID_STATE
        ELSE IF (Q_Workspace%Factors%order() /= CurrentOrder) THEN
          RootError = Q_METHOD_INVALID_STATE
        ENDIF
        IF (RootError == Q_METHOD_SUCCESS) THEN
          DO WHILE (Q_Workspace%Factors%order() > TargetOrder)
            CALL Q_Workspace%Factors%delete_symmetric(&
              Q_Workspace%Factors%order(), RootError)
            IF (RootError /= Q_METHOD_SUCCESS) EXIT
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    ! Publish the smaller physical order before recovery. Its leading principal
    ! block was never modified by appending or deleting a suffix, so it is a
    ! sound source for a fresh factorization if a structural operation failed.
    Q_Workspace%MatrixOrder = TargetOrder
    Q_Workspace%NumActive = 0
    Q_Workspace%ActiveFunction = 0
    Q_Workspace%ActivePosition = 0
    Q_Workspace%MatrixParametersAreStored = .FALSE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    Q_Workspace%AcceptedPointIsStored = .FALSE.
    Q_Workspace%FactorsMatchMatrices = &
      (RootError == Q_METHOD_SUCCESS) .AND. (TargetOrder > 0)

    IF ((RootError /= Q_METHOD_SUCCESS) .AND. (TargetOrder > 0)) THEN
      CALL FactorizeQFresh(RecoveryError)
      IF (RecoveryError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RecoveryError
      ELSE
        ErrorCode = RootError
      ENDIF
      RETURN
    ENDIF

    ErrorCode = RootError

  END SUBROUTINE TrimQFactors

  FUNCTION QCanonicalMatrixElement(Matrix, i, j) RESULT(MatrixElement)
    ! Function QCanonicalMatrixElement reads a conceptual symmetric matrix element
    ! from a matrix whose lower triangle, including the diagonal, is authoritative.
    ! No Q routine may read the upper triangle directly because it may contain old
    ! G-solver workspace or undefined values.
    !
    ! Arguments:
    REAL(wp), INTENT(IN) :: Matrix(:, :)
    INTEGER, INTENT(IN)  :: i, j
    REAL(wp)             :: MatrixElement

    IF (i >= j) THEN
      MatrixElement = Matrix(i, j)
    ELSE
      MatrixElement = Matrix(j, i)
    ENDIF

  END FUNCTION QCanonicalMatrixElement

  SUBROUTINE GatherQCanonicalColumn(Matrix, MatrixOrder, ColumnIndex, Column, ErrorCode)
    ! Subroutine GatherQCanonicalColumn constructs the complete conceptual column
    ! required by qrlinalg replace_symmetric and append_symmetric. The source matrix
    ! retains only its canonical lower triangle. This O(n) gather is negligible
    ! beside the O(n*n) QR update and avoids maintaining two writable triangles.
    !
    ! Arguments:
    REAL(wp), INTENT(IN)  :: Matrix(:, :)
    INTEGER, INTENT(IN)   :: MatrixOrder, ColumnIndex
    REAL(wp), INTENT(OUT) :: Column(:)
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER :: i

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (MatrixOrder < 1) RETURN
    IF ((SIZE(Matrix, 1) < MatrixOrder) .OR. (SIZE(Matrix, 2) < MatrixOrder)) RETURN
    IF ((ColumnIndex < 1) .OR. (ColumnIndex > MatrixOrder)) RETURN
    IF (SIZE(Column) /= MatrixOrder) RETURN

    DO i = 1, MatrixOrder
      Column(i) = QCanonicalMatrixElement(Matrix, i, ColumnIndex)
    ENDDO
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE GatherQCanonicalColumn

  SUBROUTINE StoreQCanonicalColumn(Matrix, MatrixOrder, ColumnIndex, Column, ErrorCode)
    ! Subroutine StoreQCanonicalColumn commits one complete symmetric column to the
    ! canonical lower triangle. The upper triangle is intentionally untouched.
    ! For a replacement this routine is called only after qrlinalg has accepted the
    ! matching delta; during a multi-column transaction it is called after every
    ! successful sequential update so intersections use the progressively updated
    ! matrix and are not counted twice.
    !
    ! Arguments:
    REAL(wp), INTENT(INOUT) :: Matrix(:, :)
    INTEGER, INTENT(IN)     :: MatrixOrder, ColumnIndex
    REAL(wp), INTENT(IN)    :: Column(:)
    INTEGER, INTENT(OUT)    :: ErrorCode
    ! Local variables:
    INTEGER :: i

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (MatrixOrder < 1) RETURN
    IF ((SIZE(Matrix, 1) < MatrixOrder) .OR. (SIZE(Matrix, 2) < MatrixOrder)) RETURN
    IF ((ColumnIndex < 1) .OR. (ColumnIndex > MatrixOrder)) RETURN
    IF (SIZE(Column) /= MatrixOrder) RETURN

    DO i = 1, ColumnIndex
      Matrix(ColumnIndex, i) = Column(i)
    ENDDO
    DO i = ColumnIndex+1, MatrixOrder
      Matrix(i, ColumnIndex) = Column(i)
    ENDDO
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE StoreQCanonicalColumn

  SUBROUTINE AssembleQTrial(AreDerivativesNeeded, ErrorCode)
    ! Subroutine AssembleQTrial calculates every unordered matrix-element pair
    ! that touches an active function. Active diagonals must be calculated first so
    ! all raw norms are available before normalized off-diagonal elements are
    ! formed. The final target columns are staged in Q_Workspace%TrialH and TrialS;
    ! Glob_H and Glob_S remain unchanged until ApplyQTrial succeeds.
    !
    ! When derivatives are requested, Glob_D(:,a,j) will mean the G-compatible
    ! scaled derivatives with respect to canonical function
    ! Q_Workspace%ActiveFunction(a), paired with canonical function j. This removes
    ! the old assumption that the differentiated functions occupy a trailing block.
    !
    ! Arguments:
    LOGICAL, INTENT(IN)  :: AreDerivativesNeeded
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: a, b, i, j, PairNumber, MatrixOrder, NumActive
    INTEGER :: ActiveIndex, NumMatrixEntries, NumDerivativeEntries, npt, npt2
    INTEGER :: mActive, mOther
    INTEGER :: gradflag  ! 0 no derivatives, 1 active side only, 2 both sides
    ! Bra side (the active function): Lk, Ak = Lk*Lk' and M*Ak, built once
    ! per active function
    REAL(wp) :: PreLk(Glob_n, Glob_n, 1), PreAk(Glob_n, Glob_n, 1), PreMAk(Glob_n, Glob_n, 1)
    REAL(wp) :: ParamOther(Glob_n*(Glob_n+1)/2)
    REAL(wp) :: Hsum, Ssum, ActiveNorm, OtherNorm, Normalization
    REAL(wp) :: DActiveSum(2*Glob_npt_MaxAllowed), DOtherSum(2*Glob_npt_MaxAllowed)
    LOGICAL  :: OtherDerivativeNeeded
    ! Work arrays of MatrixElementsOpt, one entry per symmetry term (as in
    ! ComputeMatElem and ComputeMatElemAndDeriv of matform)
    REAL(wp) :: SymMatrixBuf(Glob_NumYHYTerms, Glob_n, Glob_n)
    REAL(wp) :: SklBuf(Glob_NumYHYTerms), HklBuf(Glob_NumYHYTerms)
    REAL(wp) :: dSkBuf(Glob_NumYHYTerms, Glob_npt), dSlBuf(Glob_NumYHYTerms, Glob_npt)
    REAL(wp) :: dHkBuf(Glob_NumYHYTerms, Glob_npt), dHlBuf(Glob_NumYHYTerms, Glob_npt)

    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    NumActive = Q_Workspace%NumActive
    npt = Glob_npt
    npt2 = 2*Glob_npt
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF ((MatrixOrder < 1) .OR. (NumActive < 1)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActivePosition)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialS)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousS)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialDiagS)) RETURN
    IF (.NOT. ALLOCATED(Glob_NonlinParam)) RETURN
    IF (.NOT. ALLOCATED(Glob_diagS)) RETURN
    IF (SIZE(Glob_NonlinParam, 1) < Glob_npt) RETURN
    IF (SIZE(Glob_NonlinParam, 2) < MatrixOrder) RETURN
    IF (SIZE(Glob_diagS) < MatrixOrder) RETURN
    IF (Glob_NumYHYTerms < 1) RETURN
    IF (AreDerivativesNeeded) THEN
      IF (.NOT. ALLOCATED(Glob_D)) RETURN
      IF (SIZE(Glob_D, 1) < npt2) RETURN
      IF (SIZE(Glob_D, 2) < NumActive) RETURN
      IF (SIZE(Glob_D, 3) < MatrixOrder) RETURN
      Glob_D = ZERO
    ENDIF

    ! Each physical pair that touches the active set is evaluated exactly once,
    ! by ONE rank (pair p goes to rank MOD(p-1,Glob_NumOfProcs)); the other
    ! ranks leave zeros and the MPI_ALLREDUCE below assembles the columns.
    ! MatrixElementsOpt returns all Glob_NumYHYTerms symmetry terms of a pair
    ! in one call, so the split is per pair, as in ComputeMatElem. For an
    ! active-active pair the value is copied into both conceptual trial
    ! columns. ActivePosition supplies an ordering independent of canonical
    ! basis indices, so this remains correct for descending and noncontiguous
    ! active lists.
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,l1)') 'AssembleQTrial: order ', MatrixOrder, ', active functions ', &
                                       NumActive, ', derivatives ', AreDerivativesNeeded
    ENDIF
    DO j = 1, Glob_NumYHYTerms
      SymMatrixBuf(j, :, :) = Glob_YHYMatr(1:Glob_n, 1:Glob_n, j)
    ENDDO
    Q_Workspace%TrialH = ZERO
    Q_Workspace%TrialS = ZERO
    PairNumber = 0
    DO a = 1, NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      CALL PrecomputeMatrices_MA_L_A(Glob_npt, 1, Glob_NonlinParam(1:Glob_npt, ActiveIndex:ActiveIndex), &
                                     Glob_MassMatrix(1:Glob_n, 1:Glob_n), PreLk, PreAk, PreMAk)
      mActive = Glob_PWR(ActiveIndex)
      DO i = 1, MatrixOrder
        b = Q_Workspace%ActivePosition(i)
        IF ((b > 0) .AND. (b < a)) CYCLE

        PairNumber = PairNumber+1
        ParamOther(1:Glob_npt) = Glob_NonlinParam(1:Glob_npt, i)
        mOther = Glob_PWR(i)
        Hsum = ZERO
        Ssum = ZERO
        IF (AreDerivativesNeeded) THEN
          DActiveSum(1:npt2) = ZERO
          DOtherSum(1:npt2) = ZERO
        ENDIF
        OtherDerivativeNeeded = AreDerivativesNeeded .AND. (b > 0) .AND. (b /= a)
        IF (MOD(PairNumber-1, Glob_NumOfProcs) == Glob_ProcID) THEN
          gradflag = 0
          IF (AreDerivativesNeeded) gradflag = 1
          IF (OtherDerivativeNeeded) gradflag = 2
          ! The active function is the k side of the kernel, the other one
          ! the l side; dHk/dSk and dHl/dSl are their derivatives. The layout
          ! of the sums is that of StoreHSD: dH first, then dS.
          CALL MatrixElementsOpt(mActive, PreLk(:, :, 1), PreAk(:, :, 1), PreMAk(:, :, 1), &
                                 mOther, ParamOther, SymMatrixBuf, &
                                 SklBuf, HklBuf, dSkBuf, dSlBuf, dHkBuf, dHlBuf, gradflag)
          DO j = 1, Glob_NumYHYTerms
            Hsum = Hsum+Glob_YHYCoeff(j)*HklBuf(j)
            Ssum = Ssum+Glob_YHYCoeff(j)*SklBuf(j)
          ENDDO
          IF (AreDerivativesNeeded) THEN
            DO j = 1, Glob_NumYHYTerms
              DActiveSum(1:npt) = DActiveSum(1:npt)+Glob_YHYCoeff(j)*dHkBuf(j, 1:npt)
              DActiveSum(npt+1:npt2) = DActiveSum(npt+1:npt2)+Glob_YHYCoeff(j)*dSkBuf(j, 1:npt)
            ENDDO
            IF (OtherDerivativeNeeded) THEN
              DO j = 1, Glob_NumYHYTerms
                DOtherSum(1:npt) = DOtherSum(1:npt)+Glob_YHYCoeff(j)*dHlBuf(j, 1:npt)
                DOtherSum(npt+1:npt2) = DOtherSum(npt+1:npt2)+Glob_YHYCoeff(j)*dSlBuf(j, 1:npt)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        Q_Workspace%TrialH(i, a) = Hsum
        Q_Workspace%TrialS(i, a) = Ssum
        IF (AreDerivativesNeeded) &
          Glob_D(1:npt2, a, i) = DActiveSum(1:npt2)
        IF ((b > 0) .AND. (b /= a)) THEN
          Q_Workspace%TrialH(ActiveIndex, b) = Hsum
          Q_Workspace%TrialS(ActiveIndex, b) = Ssum
          IF (AreDerivativesNeeded) &
            Glob_D(1:npt2, b, ActiveIndex) = DOtherSum(1:npt2)
        ENDIF
      ENDDO
    ENDDO

    ! The first n rows and m columns are not contiguous when Capacity is larger
    ! than MatrixOrder. Reducing the complete allocated arrays preserves their
    ! physical leading dimensions and avoids an incorrectly packed MPI count.
    ! PreviousH/S are only reduction receive buffers here; ApplyQTrial replaces
    ! them with the actual pre-transaction physical columns before any update.
    NumMatrixEntries = SIZE(Q_Workspace%TrialH)
    CALL MPI_ALLREDUCE(Q_Workspace%TrialH, Q_Workspace%PreviousH, &
      NumMatrixEntries, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_ALLREDUCE(Q_Workspace%TrialS, Q_Workspace%PreviousS, &
      NumMatrixEntries, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    Q_Workspace%TrialH = Q_Workspace%PreviousH
    Q_Workspace%TrialS = Q_Workspace%PreviousS
    IF (AreDerivativesNeeded) THEN
      ! Every derivative element is owned by the process that evaluated its
      ! symmetry term. An in-place collective avoids a second replicated
      ! 2*npt-by-active-by-order tensor solely for the reduction receive side.
      NumDerivativeEntries = SIZE(Glob_D)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Glob_D, NumDerivativeEntries, &
        MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF

    ! All active raw self-overlaps must be known before any column is
    ! normalized, because an active-active element depends on both new norms.
    ! The form .not.(x>tiny) also rejects a NaN, for which the comparison is
    ! false, without requiring an additional IEEE module dependency.
    DO a = 1, NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      Q_Workspace%TrialDiagS(a) = Q_Workspace%TrialS(ActiveIndex, a)
      IF (.NOT. (Q_Workspace%TrialDiagS(a) > TINY(ONE))) RETURN
    ENDDO

    DO a = 1, NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      ActiveNorm = Q_Workspace%TrialDiagS(a)
      DO i = 1, MatrixOrder
        b = Q_Workspace%ActivePosition(i)
        IF (b > 0) THEN
          OtherNorm = Q_Workspace%TrialDiagS(b)
        ELSE
          OtherNorm = Glob_diagS(i)
        ENDIF
        IF (.NOT. (OtherNorm > TINY(ONE))) RETURN
        Normalization = ONE/SQRT(ActiveNorm*OtherNorm)
        Q_Workspace%TrialH(i, a) = Q_Workspace%TrialH(i, a)*Normalization
        Q_Workspace%TrialS(i, a) = Q_Workspace%TrialS(i, a)*Normalization
        IF (AreDerivativesNeeded) THEN
          IF (i == ActiveIndex) THEN
            ! The kernel differentiates one side of an identical bra/ket pair.
            ! Both sides contribute equally to the physical diagonal.
            Glob_D(1:npt2, a, i) = TWO*Glob_D(1:npt2, a, i)/ActiveNorm
          ELSE
            ! As in StoreHSD, keep derivatives of the raw matrix element scaled
            ! by the two basis-function norms. The normalization derivative is
            ! subtracted once, explicitly, in the energy-gradient contraction.
            Glob_D(1:npt2, a, i) = Glob_D(1:npt2, a, i)*Normalization
          ENDIF
        ENDIF
      ENDDO
      ! Set the analytically normalized diagonal exactly. This avoids allowing
      ! roundoff in Sii/Sii to enter overlap tests or qrlinalg's solve path.
      Q_Workspace%TrialS(ActiveIndex, a) = ONE
    ENDDO

    Q_Workspace%TrialIsReady = .TRUE.
    Q_Workspace%TrialHasDerivatives = AreDerivativesNeeded
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE AssembleQTrial

  SUBROUTINE FactorizeQFresh(ErrorCode)
    ! Subroutine FactorizeQFresh initializes or refreshes the root-owned qrlinalg
    ! state from canonical, normalized, unshifted Glob_H and Glob_S. The represented
    ! shift is fixed to Glob_ApproxEnergy for one BBOP step. Rank zero broadcasts the
    ! status before any process continues to a collective gradient calculation.
    !
    ! Arguments:
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: RootError

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    Q_Workspace%FactorsMatchMatrices = .FALSE.
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (Q_Workspace%MatrixOrder < 1) RETURN
    IF (Q_Workspace%Capacity < Q_Workspace%MatrixOrder) RETURN
    IF (.NOT. ALLOCATED(Glob_H)) RETURN
    IF (.NOT. ALLOCATED(Glob_S)) RETURN
    IF ((SIZE(Glob_H, 1) < Q_Workspace%MatrixOrder) .OR. &
        (SIZE(Glob_H, 2) < Q_Workspace%MatrixOrder)) RETURN
    IF ((SIZE(Glob_S, 1) < Q_Workspace%MatrixOrder) .OR. &
        (SIZE(Glob_S, 2) < Q_Workspace%MatrixOrder)) RETURN

    RootError = Q_METHOD_SUCCESS
    IF (Glob_ProcID == 0) THEN
      ! Initialization is separated from fresh factorization in qrlinalg. Reuse
      ! the allocated state whenever its capacity is already correct; this
      ! preserves the lifetime update counter across ordinary refreshes.
      IF (Q_Workspace%Factors%get_capacity() /= Q_Workspace%Capacity) THEN
        CALL Q_Workspace%Factors%initialize(Q_Workspace%Capacity, RootError)
      ENDIF
      IF (RootError == Q_METHOD_SUCCESS) THEN
        CALL Q_Workspace%Factors%factorize_fresh(Glob_H, Glob_S, &
          Glob_ApproxEnergy, RootError, active_order = Q_Workspace%MatrixOrder)
      ENDIF
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    ErrorCode = RootError
    Q_Workspace%FactorsMatchMatrices = (ErrorCode == Q_METHOD_SUCCESS)
    IF (ErrorCode == Q_METHOD_SUCCESS) &
      Q_Workspace%FreshFactorizations = Q_Workspace%FreshFactorizations+1

  END SUBROUTINE FactorizeQFresh

  SUBROUTINE AppendQTrial(ErrorCode)
    ! Subroutine AppendQTrial commits a staged suffix during BASIS_ENL. The active
    ! map must be the consecutive canonical suffix BaseOrder+1:MatrixOrder. When an
    ! accepted prefix exists, qrlinalg append_symmetric grows its factors one row
    ! and column at a time. Each physical column is committed only after the
    ! matching structural update succeeds.
    !
    ! For BaseOrder=0 there is no valid factorization that can be appended to. The
    ! routine first commits the complete staged block and then constructs the first
    ! fresh factorization. This special case is confined to the first enlargement
    ! step and therefore does not affect the asymptotic candidate-selection cost.
    !
    ! On any append failure, the represented prefix is reconstructed from its
    ! untouched leading physical block. The caller receives the original append
    ! status when recovery succeeds, and a recovery status otherwise.
    !
    ! Arguments:
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: a, BaseOrder, FunctionIndex, MatrixOrder, NumActive
    INTEGER :: RootError, StoreError, RecoveryError

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    NumActive = Q_Workspace%NumActive
    BaseOrder = MatrixOrder-NumActive
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (.NOT. Q_Workspace%TrialIsReady) RETURN
    IF ((MatrixOrder < 1) .OR. (NumActive < 1) .OR. (BaseOrder < 0)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialS)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%TrialDiagS)) RETURN
    DO a = 1, NumActive
      IF (Q_Workspace%ActiveFunction(a) /= BaseOrder+a) RETURN
    ENDDO

    RootError = Q_METHOD_SUCCESS
    IF (BaseOrder == 0) THEN
      ! Every pair in the new leading block touches an active function, so the
      ! staged columns together contain a complete symmetric problem.
      DO a = 1, NumActive
        FunctionIndex = a
        CALL StoreQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
          Q_Workspace%TrialH(1:MatrixOrder, a), StoreError)
        IF (StoreError /= Q_METHOD_SUCCESS) THEN
          ErrorCode = StoreError
          RETURN
        ENDIF
        CALL StoreQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
          Q_Workspace%TrialS(1:MatrixOrder, a), StoreError)
        IF (StoreError /= Q_METHOD_SUCCESS) THEN
          ErrorCode = StoreError
          RETURN
        ENDIF
        Glob_diagS(FunctionIndex) = Q_Workspace%TrialDiagS(a)
      ENDDO
      ! An empty input basis carries a large placeholder CURRENT_ENERGY, not a
      ! meaningful inverse-iteration target. Replace such a remote shift by the
      ! lowest normalized diagonal estimate before constructing the first QR
      ! factorization. For the usual Kstep=1 start this is the exact energy.
      IF (ABS(Glob_ApproxEnergy) > &
          1000000*MAX(ONE, MAXVAL(ABS([(Glob_H(a, a), a=1, MatrixOrder)])))) THEN
        Glob_ApproxEnergy = MINVAL([(Glob_H(a, a), a=1, MatrixOrder)])* &
          Glob_InvItParameter
      ENDIF
      CALL FactorizeQFresh(RootError)
    ELSE
      IF (Glob_ProcID == 0) THEN
        IF (.NOT. Q_Workspace%Factors%is_valid()) THEN
          RootError = Q_METHOD_INVALID_STATE
        ELSE IF (Q_Workspace%Factors%order() /= BaseOrder) THEN
          RootError = Q_METHOD_INVALID_STATE
        ENDIF
      ENDIF
      CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      IF (RootError == Q_METHOD_SUCCESS) THEN
        DO a = 1, NumActive
          FunctionIndex = BaseOrder+a
          IF (Glob_ProcID == 0) THEN
            CALL Q_Workspace%Factors%append_symmetric(&
              Q_Workspace%TrialH(1:FunctionIndex, a), &
              Q_Workspace%TrialS(1:FunctionIndex, a), RootError)
          ENDIF
          CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          IF (RootError /= Q_METHOD_SUCCESS) EXIT

          CALL StoreQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
            Q_Workspace%TrialH(1:MatrixOrder, a), StoreError)
          CALL StoreQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
            Q_Workspace%TrialS(1:MatrixOrder, a), StoreError)
          Glob_diagS(FunctionIndex) = Q_Workspace%TrialDiagS(a)
        ENDDO
      ENDIF
    ENDIF

    IF (RootError /= Q_METHOD_SUCCESS) THEN
      ! Only suffix rows and columns may have been committed. Re-expose the
      ! accepted leading block and rebuild its factors instead of retaining a
      ! partially appended generation.
      Q_Workspace%MatrixOrder = BaseOrder
      Q_Workspace%FactorsMatchMatrices = .FALSE.
      IF (BaseOrder > 0) THEN
        CALL FactorizeQFresh(RecoveryError)
      ELSE
        RecoveryError = Q_METHOD_SUCCESS
        IF (Glob_ProcID == 0) THEN
          CALL Q_Workspace%Factors%clear()
          CALL Q_Workspace%Factors%initialize(Q_Workspace%Capacity, RecoveryError)
        ENDIF
        CALL MPI_BCAST(RecoveryError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
      IF (RecoveryError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RecoveryError
      ELSE
        ErrorCode = RootError
      ENDIF
      Q_Workspace%TrialIsReady = .FALSE.
      Q_Workspace%TrialHasDerivatives = .FALSE.
      RETURN
    ENDIF

    DO a = 1, NumActive
      FunctionIndex = Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt, a) = &
        Glob_NonlinParam(1:Glob_npt, FunctionIndex)
    ENDDO
    Q_Workspace%FactorsMatchMatrices = .TRUE.
    Q_Workspace%MatrixParametersAreStored = .TRUE.
    Q_Workspace%TrialIsReady = .FALSE.
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE AppendQTrial

  SUBROUTINE EvaluateQAppendedTrial(BaseOrder, TargetOrder, Evalue, ErrorCode)
    ! Subroutine EvaluateQAppendedTrial performs one complete BASIS_ENL candidate
    ! transaction. Any preceding candidate suffix is deleted, the new consecutive
    ! suffix is assembled in canonical order, appended to the accepted QR prefix,
    ! and solved. On success the active map remains installed so the selected
    ! candidate can immediately enter ordinary replacement-based optimization.
    !
    ! Arguments:
    INTEGER, INTENT(IN)   :: BaseOrder, TargetOrder
    REAL(wp), INTENT(OUT) :: Evalue
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER :: a, NumActive
    INTEGER :: ActiveFunction(Q_Workspace%MaxActive)

    Evalue = HUGE(Evalue)
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    NumActive = TargetOrder-BaseOrder
    IF ((BaseOrder < 0) .OR. (TargetOrder > Q_Workspace%Capacity)) RETURN
    IF ((NumActive < 1) .OR. (NumActive > Q_Workspace%MaxActive)) RETURN

    CALL TrimQFactors(BaseOrder, ErrorCode)
    IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN

    ! The factors still represent BaseOrder while the assembly routines need to
    ! see the complete target problem. FactorsMatchMatrices remains false until
    ! AppendQTrial has installed every staged suffix column.
    Q_Workspace%MatrixOrder = TargetOrder
    Q_Workspace%FactorsMatchMatrices = .FALSE.
    DO a = 1, NumActive
      ActiveFunction(a) = BaseOrder+a
    ENDDO
    CALL SetQActiveFunctions(ActiveFunction(1:NumActive), ErrorCode)
    IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
    CALL AssembleQTrial(.FALSE., ErrorCode)
    IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
    CALL AppendQTrial(ErrorCode)
    IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
    CALL SolveQ(Evalue, ErrorCode)

  END SUBROUTINE EvaluateQAppendedTrial

  SUBROUTINE ApplyQTrial(ErrorCode)
    ! Subroutine ApplyQTrial updates active columns in the exact order stored in
    ! Q_Workspace%ActiveFunction. For each column it gathers the currently
    ! represented physical column, subtracts it from the staged target, calls
    ! replace_symmetric on rank zero, broadcasts the status, and only then commits
    ! the physical column on every rank. A preflight check must make all ordinary
    ! argument failures impossible before the first factor is changed.
    !
    ! Arguments:
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER :: a, FunctionIndex, MatrixOrder, NumActive, RootError, StoreError
    INTEGER :: RecoveryError

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    NumActive = Q_Workspace%NumActive
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (.NOT. Q_Workspace%FactorsMatchMatrices) RETURN
    IF (.NOT. Q_Workspace%MatrixParametersAreStored) RETURN
    IF (.NOT. Q_Workspace%TrialIsReady) RETURN
    IF ((MatrixOrder < 1) .OR. (NumActive < 1)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousS)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousDiagS)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%PreviousParam)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%DeltaH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%DeltaS)) RETURN
    IF (SIZE(Q_Workspace%DeltaH) < MatrixOrder) RETURN
    IF (SIZE(Q_Workspace%DeltaS) < MatrixOrder) RETURN

    ! Preflight the private factor metadata on rank zero. Every subsequent
    ! replace_symmetric call then has a valid index and vectors of exact order;
    ! the library documents that these validated updates have no numerical
    ! failure return.
    RootError = Q_METHOD_SUCCESS
    IF (Glob_ProcID == 0) THEN
      IF (.NOT. Q_Workspace%Factors%is_valid()) RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%order() /= MatrixOrder) RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%get_capacity() /= Q_Workspace%Capacity) &
        RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%get_shift() /= Glob_ApproxEnergy) &
        RootError = Q_METHOD_INVALID_STATE
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (RootError /= Q_METHOD_SUCCESS) THEN
      ErrorCode = RootError
      Q_Workspace%FactorsMatchMatrices = .FALSE.
      RETURN
    ENDIF

    ! Snapshot every original column before the first commit. In particular,
    ! an active-active intersection must be captured before either endpoint is
    ! changed. These copies also define the nonlinear-parameter generation to
    ! which a rare transaction recovery must return.
    Q_Workspace%PreviousParam = Q_Workspace%MatrixParam
    DO a = 1, NumActive
      FunctionIndex = Q_Workspace%ActiveFunction(a)
      CALL GatherQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
        Q_Workspace%PreviousH(1:MatrixOrder, a), StoreError)
      IF (StoreError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = StoreError
        RETURN
      ENDIF
      CALL GatherQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
        Q_Workspace%PreviousS(1:MatrixOrder, a), StoreError)
      IF (StoreError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = StoreError
        RETURN
      ENDIF
      Q_Workspace%PreviousDiagS(a) = Glob_diagS(FunctionIndex)
    ENDDO

    DO a = 1, NumActive
      FunctionIndex = Q_Workspace%ActiveFunction(a)
      ! Gather the progressively updated column, not the original snapshot.
      ! Earlier active columns have already installed their shared intersection,
      ! so the later replacement sees a zero delta at that location.
      CALL GatherQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
        Q_Workspace%DeltaH(1:MatrixOrder), StoreError)
      CALL GatherQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
        Q_Workspace%DeltaS(1:MatrixOrder), StoreError)
      Q_Workspace%DeltaH(1:MatrixOrder) = &
        Q_Workspace%TrialH(1:MatrixOrder, a)-Q_Workspace%DeltaH(1:MatrixOrder)
      Q_Workspace%DeltaS(1:MatrixOrder) = &
        Q_Workspace%TrialS(1:MatrixOrder, a)-Q_Workspace%DeltaS(1:MatrixOrder)

      RootError = Q_METHOD_SUCCESS
      IF (Glob_ProcID == 0) THEN
        CALL Q_Workspace%Factors%replace_symmetric(FunctionIndex, &
          Q_Workspace%DeltaH(1:MatrixOrder), &
          Q_Workspace%DeltaS(1:MatrixOrder), RootError)
      ENDIF
      CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      IF (RootError /= Q_METHOD_SUCCESS) EXIT

      CALL StoreQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
        Q_Workspace%TrialH(1:MatrixOrder, a), StoreError)
      CALL StoreQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
        Q_Workspace%TrialS(1:MatrixOrder, a), StoreError)
      Glob_diagS(FunctionIndex) = Q_Workspace%TrialDiagS(a)
    ENDDO

    IF (RootError /= Q_METHOD_SUCCESS) THEN
      ! A library failure is not expected after the preflight, but fail safely:
      ! restore the complete physical generation and construct fresh factors
      ! instead of trying to reason about a possibly partial QR update.
      DO a = 1, NumActive
        FunctionIndex = Q_Workspace%ActiveFunction(a)
        CALL StoreQCanonicalColumn(Glob_H, MatrixOrder, FunctionIndex, &
          Q_Workspace%PreviousH(1:MatrixOrder, a), StoreError)
        CALL StoreQCanonicalColumn(Glob_S, MatrixOrder, FunctionIndex, &
          Q_Workspace%PreviousS(1:MatrixOrder, a), StoreError)
        Glob_diagS(FunctionIndex) = Q_Workspace%PreviousDiagS(a)
        Glob_NonlinParam(1:Glob_npt, FunctionIndex) = &
          Q_Workspace%PreviousParam(1:Glob_npt, a)
      ENDDO
      Q_Workspace%MatrixParam = Q_Workspace%PreviousParam
      CALL FactorizeQFresh(RecoveryError)
      IF (RecoveryError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RecoveryError
      ELSE
        ErrorCode = RootError
      ENDIF
      Q_Workspace%TrialIsReady = .FALSE.
      Q_Workspace%TrialHasDerivatives = .FALSE.
      RETURN
    ENDIF

    DO a = 1, NumActive
      FunctionIndex = Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt, a) = &
        Glob_NonlinParam(1:Glob_npt, FunctionIndex)
    ENDDO
    Q_Workspace%FactorsMatchMatrices = .TRUE.
    Q_Workspace%MatrixParametersAreStored = .TRUE.
    Q_Workspace%TrialIsReady = .FALSE.
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE ApplyQTrial

  SUBROUTINE ComputeQEigenpairResidual(Evalue, Eigenvector, AbsoluteResidual, &
                                       RelativeResidual, ErrorCode)
    ! Subroutine ComputeQEigenpairResidual measures the residual of the physical
    ! generalized eigenproblem, independently of qrlinalg's direction-change test:
    !
    !                  r = H*c-Evalue*S*c .
    !
    ! The returned relative value is ||r||/(||H*c||+|E|*||S*c||+tiny). Only the
    ! canonical lower triangles are read. This routine is serial and is called on
    ! rank zero; it performs no allocation and does not modify H, S, or c.
    !
    ! Arguments:
    REAL(wp), INTENT(IN)  :: Evalue, Eigenvector(:)
    REAL(wp), INTENT(OUT) :: AbsoluteResidual, RelativeResidual
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER  :: MatrixOrder
    REAL(wp) :: HNorm2, SNorm2, ResidualNorm2

    AbsoluteResidual = ZERO
    RelativeResidual = ZERO
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    IF (MatrixOrder < 1) RETURN
    IF (SIZE(Eigenvector) /= MatrixOrder) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%DeltaH)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%DeltaS)) RETURN

    ! The earlier scalar double loop performed the same two symmetric products
    ! but left optimized BLAS performance unused. DSYMV reads exactly the
    ! authoritative lower triangles, handles the physical leading dimension,
    ! and is substantially faster for the residual that follows every Q solve.
    ! This helper is called only on rank zero, so calling BLAS directly is also
    ! important: the MPI-routing MTMVL wrapper is collective for some calibrated
    ! matrix sizes and therefore cannot be entered by the root alone.
    CALL DSYMV('L', MatrixOrder, ONE, Glob_H, SIZE(Glob_H, 1), Eigenvector, 1, &
      ZERO, Q_Workspace%DeltaH, 1)
    CALL DSYMV('L', MatrixOrder, ONE, Glob_S, SIZE(Glob_S, 1), Eigenvector, 1, &
      ZERO, Q_Workspace%DeltaS, 1)
    HNorm2 = DOT_PRODUCT(Q_Workspace%DeltaH(1:MatrixOrder), &
                      Q_Workspace%DeltaH(1:MatrixOrder))
    SNorm2 = DOT_PRODUCT(Q_Workspace%DeltaS(1:MatrixOrder), &
                      Q_Workspace%DeltaS(1:MatrixOrder))
    Q_Workspace%DeltaH(1:MatrixOrder) = Q_Workspace%DeltaH(1:MatrixOrder)- &
      Evalue*Q_Workspace%DeltaS(1:MatrixOrder)
    ResidualNorm2 = DOT_PRODUCT(Q_Workspace%DeltaH(1:MatrixOrder), &
                             Q_Workspace%DeltaH(1:MatrixOrder))
    AbsoluteResidual = SQRT(ResidualNorm2)
    RelativeResidual = AbsoluteResidual/ &
      (SQRT(HNorm2)+ABS(Evalue)*SQRT(SNorm2)+TINY(ONE))
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE ComputeQEigenpairResidual

  SUBROUTINE SolveQ(Evalue, ErrorCode)
    ! Subroutine SolveQ calls qrlinalg inverse iteration on rank zero with
    ! distinct input and output vectors and normalization mode zero, then broadcast
    ! the physical energy, S-normalized Glob_c, convergence diagnostics, and status.
    ! A QR_ERR_NO_CONVERGENCE result contains a usable approximation. It is accepted
    ! only when the independent physical generalized-eigenpair residual satisfies
    ! the requested tolerance; otherwise one fresh-factorization retry is made.
    !
    ! Arguments:
    REAL(wp), INTENT(OUT) :: Evalue
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER        :: i, MatrixOrder, NumOfIterations, RootError, RefreshError, SolveStatus
    INTEGER        :: RetryMaxIterations
    INTEGER(int64) :: UpdatesSinceFresh, RefreshLimit
    REAL(wp)       :: AbsoluteResidual, RelativeResidual, FactorResidual
    REAL(wp)       :: RefreshTolerance, SolveTolerance
    LOGICAL        :: RefreshNeeded

    Evalue = HUGE(Evalue)
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (.NOT. Q_Workspace%FactorsMatchMatrices) RETURN
    IF (MatrixOrder < 1) RETURN
    IF (.NOT. ALLOCATED(Glob_c)) RETURN
    IF (SIZE(Glob_c) < MatrixOrder) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%InitialVector)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%SolvedVector)) RETURN

    ! A one-function generalized problem has an exact closed-form solution.
    ! Besides avoiding unnecessary inverse iteration, this is essential when a
    ! calculation starts from the conventional enormous CURRENT_ENERGY sentinel:
    ! forming lambda as shift+(lambda-shift) would otherwise lose the physical
    ! diagonal energy by catastrophic cancellation. BASIS_ENL retargets the QR
    ! shift to this first accepted energy before attempting order two.
    IF (MatrixOrder == 1) THEN
      IF (.NOT. (Glob_S(1, 1) > TINY(ONE))) THEN
        ErrorCode = QR_ERR_NONPOSITIVE_OVERLAP
        RETURN
      ENDIF
      Evalue = Glob_H(1, 1)/Glob_S(1, 1)
      Glob_c(1) = ONE/SQRT(Glob_S(1, 1))
      Q_Workspace%LastEigenpairResidual = ZERO
      Q_Workspace%LastFactorResidual = ZERO
      Glob_LastEigvalTol = ZERO
      Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2 = Glob_InvItTempCounter2+1
      ErrorCode = Q_METHOD_SUCCESS
      RETURN
    ENDIF

    ! A full refresh after O(n) replacements keeps the amortized cost of fresh
    ! O(n**3) factorizations at O(n**2) per replacement. The residual check
    ! below can request an earlier refresh when accumulated rotations drift.
    RefreshNeeded = .FALSE.
    IF (Glob_ProcID == 0) THEN
      IF (.NOT. Q_Workspace%Factors%is_valid()) THEN
        RootError = Q_METHOD_INVALID_STATE
      ELSE
        RootError = Q_METHOD_SUCCESS
        UpdatesSinceFresh = Q_Workspace%Factors%get_updates_since_fresh()
        RefreshLimit = INT(MAX(64, 8*MatrixOrder), int64)
        RefreshNeeded = (UpdatesSinceFresh >= RefreshLimit)
      ENDIF
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(RefreshNeeded, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (RootError /= Q_METHOD_SUCCESS) THEN
      ErrorCode = RootError
      Q_Workspace%FactorsMatchMatrices = .FALSE.
      RETURN
    ENDIF
    IF (RefreshNeeded) THEN
      CALL FactorizeQFresh(RefreshError)
      IF (RefreshError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RefreshError
        RETURN
      ENDIF
    ENDIF

    RootError = Q_METHOD_SUCCESS
    AbsoluteResidual = ZERO
    RelativeResidual = ZERO
    NumOfIterations = 0
    IF (Glob_ProcID == 0) THEN
      Q_Workspace%InitialVector(1:MatrixOrder) = Glob_c(1:MatrixOrder)
      IF (.NOT. (MAXVAL(ABS(Q_Workspace%InitialVector(1:MatrixOrder))) > TINY(ONE))) THEN
        DO i = 1, MatrixOrder
          Q_Workspace%InitialVector(i) = ONE
        ENDDO
      ENDIF
      CALL Q_Workspace%Factors%solve(Glob_S, &
        Q_Workspace%InitialVector(1:MatrixOrder), &
        Q_Workspace%SolvedVector(1:MatrixOrder), Evalue, Glob_EigvalTol, &
        Glob_MaxIterForGSEPIIS, 0, RelativeResidual, NumOfIterations, RootError)
      Glob_c(1:MatrixOrder) = Q_Workspace%SolvedVector(1:MatrixOrder)
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(RelativeResidual, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(NumOfIterations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_c, MatrixOrder, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    SolveStatus = RootError
    Glob_LastEigvalTol = RelativeResidual
    Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
    IF ((SolveStatus /= Q_METHOD_SUCCESS) .AND. &
        (SolveStatus /= QR_ERR_NO_CONVERGENCE)) THEN
      ErrorCode = RootError
      RETURN
    ENDIF

    ! Measure both independent numerical relationships. The physical residual
    ! checks H*c=E*S*c. The factor residual checks Q*R=(H-shift*S) with a fixed
    ! dense probe. The eigenvector must not be used for the latter because it
    ! is nearly a null vector of the deliberately near-eigenvalue shifted
    ! matrix, which amplifies harmless factorization roundoff in the relative
    ! action ratio. QR_ERR_NO_CONVERGENCE still contains an approximation; it is
    ! accepted only when this independent physical residual satisfies the
    ! caller tolerance. Either diagnostic may request one fresh-factorization
    ! retry.
    IF (Glob_ProcID == 0) THEN
      RootError = Q_METHOD_SUCCESS
      CALL ComputeQEigenpairResidual(Evalue, Glob_c(1:MatrixOrder), &
        AbsoluteResidual, Q_Workspace%LastEigenpairResidual, RootError)
      IF (RootError == Q_METHOD_SUCCESS) THEN
        DO i = 1, MatrixOrder
          Q_Workspace%InitialVector(i) = ONE+ &
            REAL(MOD(17*i, 23), wp)/REAL(23, wp)
          IF (MOD(i, 2) == 0) Q_Workspace%InitialVector(i) = &
            -Q_Workspace%InitialVector(i)
        ENDDO
        CALL Q_Workspace%Factors%factorization_residual(Glob_H, Glob_S, &
          Q_Workspace%InitialVector(1:MatrixOrder), AbsoluteResidual, &
          FactorResidual, RootError)
        Q_Workspace%LastFactorResidual = FactorResidual
      ENDIF
      RefreshTolerance = MAX(1000*EPSILON(ONE)*MatrixOrder, &
                             10*ABS(Glob_EigvalTol))
      SolveTolerance = MAX(1000*EPSILON(ONE)*MatrixOrder, &
                           100*ABS(Glob_EigvalTol))
      RefreshNeeded = (RootError == Q_METHOD_SUCCESS) .AND. &
        ((.NOT. (Q_Workspace%LastFactorResidual <= RefreshTolerance)) .OR. &
         (.NOT. (Q_Workspace%LastEigenpairResidual <= SolveTolerance)))
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(RefreshNeeded, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Q_Workspace%LastEigenpairResidual, 1, MPI_WP, 0, &
      MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Q_Workspace%LastFactorResidual, 1, MPI_WP, 0, &
      MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (RootError /= Q_METHOD_SUCCESS) THEN
      ErrorCode = RootError
      RETURN
    ENDIF

    IF (RefreshNeeded) THEN
      CALL FactorizeQFresh(RefreshError)
      IF (RefreshError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RefreshError
        RETURN
      ENDIF
      ! The normal iteration cap is tuned for small optimizer displacements.
      ! A newly appended basis function can rotate the eigenvector much farther.
      ! Only after the independent residual has rejected the first approximation
      ! do we allow this larger cap; every additional iteration remains O(n**2).
      RetryMaxIterations = MAX(120, MAX(Glob_MaxIterForGSEPIIS, 4*MatrixOrder))
      IF (Glob_ProcID == 0) THEN
        Q_Workspace%InitialVector(1:MatrixOrder) = Glob_c(1:MatrixOrder)
        CALL Q_Workspace%Factors%solve(Glob_S, &
          Q_Workspace%InitialVector(1:MatrixOrder), &
          Q_Workspace%SolvedVector(1:MatrixOrder), Evalue, Glob_EigvalTol, &
          RetryMaxIterations, 0, RelativeResidual, NumOfIterations, RootError)
        Glob_c(1:MatrixOrder) = Q_Workspace%SolvedVector(1:MatrixOrder)
      ENDIF
      CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(RelativeResidual, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(NumOfIterations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_c, MatrixOrder, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      Glob_LastEigvalTol = RelativeResidual
      Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations

      IF ((RootError == Q_METHOD_SUCCESS) .OR. &
          (RootError == QR_ERR_NO_CONVERGENCE)) THEN
        IF (Glob_ProcID == 0) THEN
          RootError = Q_METHOD_SUCCESS
          CALL ComputeQEigenpairResidual(Evalue, Glob_c(1:MatrixOrder), &
            AbsoluteResidual, Q_Workspace%LastEigenpairResidual, RootError)
          IF (RootError == Q_METHOD_SUCCESS) THEN
            DO i = 1, MatrixOrder
              Q_Workspace%InitialVector(i) = ONE+ &
                REAL(MOD(17*i, 23), wp)/REAL(23, wp)
              IF (MOD(i, 2) == 0) Q_Workspace%InitialVector(i) = &
                -Q_Workspace%InitialVector(i)
            ENDDO
            CALL Q_Workspace%Factors%factorization_residual(Glob_H, Glob_S, &
              Q_Workspace%InitialVector(1:MatrixOrder), AbsoluteResidual, &
              FactorResidual, RootError)
            Q_Workspace%LastFactorResidual = FactorResidual
          ENDIF
          IF ((RootError == Q_METHOD_SUCCESS) .AND. &
              (.NOT. (Q_Workspace%LastEigenpairResidual <= SolveTolerance))) &
            RootError = QR_ERR_NO_CONVERGENCE
          IF ((RootError == Q_METHOD_SUCCESS) .AND. &
              (.NOT. (Q_Workspace%LastFactorResidual <= RefreshTolerance))) &
            RootError = Q_METHOD_INVALID_STATE
        ENDIF
        CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Q_Workspace%LastEigenpairResidual, 1, MPI_WP, 0, &
          MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Q_Workspace%LastFactorResidual, 1, MPI_WP, 0, &
          MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
    ENDIF
    ErrorCode = RootError

  END SUBROUTINE SolveQ

  SUBROUTINE ComputeQOverlapPenalty(TotalPenalty, ErrorCode)
    ! Subroutine ComputeQOverlapPenalty evaluates the FULL_OPT1 smooth pair-overlap
    ! penalty without assuming that optimized functions form a trailing block. Each
    ! unordered pair touching the explicit active map is visited exactly once and
    ! assigned to one MPI rank. Glob_S is read only through its canonical lower
    ! triangle.
    !
    ! Arguments:
    REAL(wp), INTENT(OUT) :: TotalPenalty
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER  :: a, b, i, PairNumber, ActiveIndex
    REAL(wp) :: PairOverlap, PairPenalty, LocalPenalty, PenaltyCoefficient

    TotalPenalty = ZERO
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (Q_Workspace%NumActive < 1) RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActiveFunction)) RETURN
    IF (.NOT. ALLOCATED(Q_Workspace%ActivePosition)) RETURN
    IF (.NOT. (Glob_OverlapPenaltyThreshold2 < ONE)) RETURN

    PenaltyCoefficient = Glob_MaxOverlapPenalty/ &
      (ONE-Glob_OverlapPenaltyThreshold2)
    LocalPenalty = ZERO
    PairNumber = 0
    DO a = 1, Q_Workspace%NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      DO i = 1, Q_Workspace%MatrixOrder
        IF (i == ActiveIndex) CYCLE
        b = Q_Workspace%ActivePosition(i)
        IF ((b > 0) .AND. (b < a)) CYCLE
        PairNumber = PairNumber+1
        IF (MOD(PairNumber-1, Glob_NumOfProcs) /= Glob_ProcID) CYCLE
        PairOverlap = QCanonicalMatrixElement(Glob_S, ActiveIndex, i)
        IF (PairOverlap*PairOverlap > Glob_OverlapPenaltyThreshold2) THEN
          PairPenalty = PenaltyCoefficient* &
            (PairOverlap*PairOverlap-Glob_OverlapPenaltyThreshold2)
          LocalPenalty = LocalPenalty+PairPenalty
        ENDIF
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(LocalPenalty, TotalPenalty, 1, MPI_WP, MPI_SUM, &
      MPI_COMM_WORLD, Glob_MPIErrCode)
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE ComputeQOverlapPenalty

  SUBROUTINE ComputeQOverlapPenaltyAndAddGradient(TotalPenalty, WkGR, ErrorCode)
    ! Subroutine ComputeQOverlapPenaltyAndAddGradient evaluates the same mapped
    ! penalty as ComputeQOverlapPenalty and adds its analytic derivative to WkGR.
    ! The derivative tensor produced by AssembleQTrial stores
    !
    !  d<Sraw_ij>/sqrt(Sraw_ii*Sraw_jj)
    !
    ! for an active endpoint, while its active diagonal stores
    ! dSraw_ii/Sraw_ii. Consequently the derivative of a normalized overlap is the
    ! cross derivative minus one half of the overlap times the diagonal derivative.
    ! For an active-active pair this expression is applied independently to both
    ! endpoints. Pair ownership follows MPI rank, so the caller's subsequent
    ! gradient MPI_ALLREDUCE combines both the energy-gradient and penalty pieces.
    !
    ! Arguments:
    REAL(wp), INTENT(OUT)   :: TotalPenalty
    REAL(wp), INTENT(INOUT) :: WkGR(:)
    INTEGER, INTENT(OUT)    :: ErrorCode
    ! Local variables:
    INTEGER  :: a, b, i, m, PairNumber, ActiveIndex
    REAL(wp) :: PairOverlap, LocalPenalty, PenaltyCoefficient
    REAL(wp) :: GradientCoefficient, OverlapDerivative

    TotalPenalty = ZERO
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (Q_Workspace%NumActive < 1) RETURN
    IF (.NOT. Q_Workspace%TrialHasDerivatives) RETURN
    IF (SIZE(WkGR) < Q_Workspace%NumActive*Glob_npt) RETURN
    IF (.NOT. ALLOCATED(Glob_D)) RETURN
    IF (.NOT. (Glob_OverlapPenaltyThreshold2 < ONE)) RETURN

    PenaltyCoefficient = Glob_MaxOverlapPenalty/ &
      (ONE-Glob_OverlapPenaltyThreshold2)
    LocalPenalty = ZERO
    PairNumber = 0
    DO a = 1, Q_Workspace%NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      DO i = 1, Q_Workspace%MatrixOrder
        IF (i == ActiveIndex) CYCLE
        b = Q_Workspace%ActivePosition(i)
        IF ((b > 0) .AND. (b < a)) CYCLE
        PairNumber = PairNumber+1
        IF (MOD(PairNumber-1, Glob_NumOfProcs) /= Glob_ProcID) CYCLE
        PairOverlap = QCanonicalMatrixElement(Glob_S, ActiveIndex, i)
        IF (PairOverlap*PairOverlap > Glob_OverlapPenaltyThreshold2) THEN
          LocalPenalty = LocalPenalty+PenaltyCoefficient* &
            (PairOverlap*PairOverlap-Glob_OverlapPenaltyThreshold2)
          GradientCoefficient = TWO*PenaltyCoefficient*PairOverlap
          DO m = 1, Glob_npt
            OverlapDerivative = Glob_D(Glob_npt+m, a, i)-ONEHALF* &
              PairOverlap*Glob_D(Glob_npt+m, a, ActiveIndex)
            WkGR((a-1)*Glob_npt+m) = WkGR((a-1)*Glob_npt+m)+ &
              GradientCoefficient*OverlapDerivative
          ENDDO
          IF (b > 0) THEN
            DO m = 1, Glob_npt
              OverlapDerivative = Glob_D(Glob_npt+m, b, ActiveIndex)-ONEHALF* &
                PairOverlap*Glob_D(Glob_npt+m, b, i)
              WkGR((b-1)*Glob_npt+m) = WkGR((b-1)*Glob_npt+m)+ &
                GradientCoefficient*OverlapDerivative
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(LocalPenalty, TotalPenalty, 1, MPI_WP, MPI_SUM, &
      MPI_COMM_WORLD, Glob_MPIErrCode)
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE ComputeQOverlapPenaltyAndAddGradient

  SUBROUTINE GetQOverlapStatistics(MaxAbsOverlap, MinAbsOverlap, &
                                   AverageAbsOverlap, ErrorCode)
    ! Subroutine GetQOverlapStatistics reports statistics for all unordered pairs
    ! touching the active map. It is the canonical-index counterpart of
    ! GetOverlapStatistics, whose Nmin:Nmax interface assumes a trailing block.
    !
    ! Arguments:
    REAL(wp), INTENT(OUT) :: MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER  :: a, b, i, NumPairs, ActiveIndex
    REAL(wp) :: PairOverlap, AbsPairOverlap

    MaxAbsOverlap = ZERO
    MinAbsOverlap = HUGE(ONE)
    AverageAbsOverlap = ZERO
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (Q_Workspace%NumActive < 1) RETURN

    NumPairs = 0
    DO a = 1, Q_Workspace%NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      DO i = 1, Q_Workspace%MatrixOrder
        IF (i == ActiveIndex) CYCLE
        b = Q_Workspace%ActivePosition(i)
        IF ((b > 0) .AND. (b < a)) CYCLE
        PairOverlap = QCanonicalMatrixElement(Glob_S, ActiveIndex, i)
        AbsPairOverlap = ABS(PairOverlap)
        IF (AbsPairOverlap > ABS(MaxAbsOverlap)) MaxAbsOverlap = PairOverlap
        IF (AbsPairOverlap < ABS(MinAbsOverlap)) MinAbsOverlap = PairOverlap
        AverageAbsOverlap = AverageAbsOverlap+AbsPairOverlap
        NumPairs = NumPairs+1
      ENDDO
    ENDDO
    IF (NumPairs > 0) THEN
      AverageAbsOverlap = AverageAbsOverlap/NumPairs
    ELSE
      MinAbsOverlap = ZERO
    ENDIF
    ErrorCode = Q_METHOD_SUCCESS

  END SUBROUTINE GetQOverlapStatistics

  FUNCTION EnergyQA(AreMatElemNeeded, ErrorCode)
    ! Function EnergyQA provides the Q counterpart of EnergyGA. The active set
    ! is taken from Q_Workspace rather than encoded as a trailing Nmin:Nmax range.
    ! It assembles and applies a trial when AreMatElemNeeded is true, solves the
    ! represented problem, and returns the requested physical eigenvalue.
    !
    ! Arguments:
    LOGICAL, INTENT(IN)  :: AreMatElemNeeded
    INTEGER, INTENT(OUT) :: ErrorCode
    REAL(wp)             :: EnergyQA

    EnergyQA = HUGE(EnergyQA)
    ErrorCode = Q_METHOD_SUCCESS
    IF (AreMatElemNeeded) THEN
      CALL AssembleQTrial(.FALSE., ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
      CALL ApplyQTrial(ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
    ENDIF
    CALL SolveQ(EnergyQA, ErrorCode)

    IF ((ErrorCode == Q_METHOD_SUCCESS) .AND. Glob_OverlapPenaltyAllowed) THEN
      CALL ComputeQOverlapPenalty(Glob_TotalOverlapPenalty, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) &
        EnergyQA = EnergyQA+Glob_TotalOverlapPenalty
    ENDIF
    Glob_EnergyGACounter = Glob_EnergyGACounter+1

  END FUNCTION EnergyQA

  FUNCTION EnergyQAM(AreMatElemNeeded, ErrorCode)
    ! Function EnergyQAM provides the Q counterpart of EnergyGAM. qrlinalg
    ! always computes an eigenvector, so EnergyQA and EnergyQAM may share one solve;
    ! the separate entry point keeps G's candidate-acceptance call structure clear.
    !
    ! Arguments:
    LOGICAL, INTENT(IN)  :: AreMatElemNeeded
    INTEGER, INTENT(OUT) :: ErrorCode
    REAL(wp)             :: EnergyQAM

    EnergyQAM = EnergyQA(AreMatElemNeeded, ErrorCode)

  END FUNCTION EnergyQAM

  SUBROUTINE EnergyQB(Evalue, Gradient, AreMatElemNeeded, ErrorCode)
    ! Subroutine EnergyQB provides the Q counterpart of EnergyGB. Gradient block
    ! a corresponds to Q_Workspace%ActiveFunction(a); contractions must use the
    ! coefficient of that canonical function instead of Glob_c(a+Glob_nfru).
    !
    ! Arguments:
    REAL(wp), INTENT(OUT) :: Evalue
    REAL(wp), INTENT(OUT) :: Gradient(:)
    LOGICAL, INTENT(IN)   :: AreMatElemNeeded
    INTEGER, INTENT(OUT)  :: ErrorCode
    ! Local variables:
    INTEGER  :: a, l, m, npt, MatrixOrder, NumActive, ActiveIndex
    REAL(wp) :: W(Glob_npt_MaxAllowed), t, t2

    Evalue = HUGE(Evalue)
    Gradient = HUGE(Evalue)
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    MatrixOrder = Q_Workspace%MatrixOrder
    NumActive = Q_Workspace%NumActive
    npt = Glob_npt
    IF (SIZE(Gradient) < NumActive*npt) RETURN
    IF (.NOT. ALLOCATED(Glob_D)) RETURN
    IF (.NOT. ALLOCATED(Glob_WkGR)) RETURN
    IF (SIZE(Glob_WkGR) < NumActive*npt) RETURN

    IF (AreMatElemNeeded) THEN
      CALL AssembleQTrial(.TRUE., ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
      CALL ApplyQTrial(ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN
    ELSE
      IF (.NOT. Q_Workspace%TrialHasDerivatives) RETURN
    ENDIF

    CALL SolveQ(Evalue, ErrorCode)
    IF (ErrorCode /= Q_METHOD_SUCCESS) RETURN

    ! The derivative tensor stores the derivative of one raw matrix element,
    ! scaled by both raw basis norms. Contracting the full conceptual row gives
    ! both symmetric H/S contributions. The final diagonal subtraction removes
    ! the duplicate diagonal and, through the eigenvalue equation, accounts for
    ! the derivative of normalization of every element touching this function.
    Glob_WkGR(1:NumActive*npt) = ZERO
    DO a = 1, NumActive
      ActiveIndex = Q_Workspace%ActiveFunction(a)
      W(1:npt) = ZERO
      DO l = 1+Glob_ProcID, MatrixOrder, Glob_NumOfProcs
        t = Glob_c(l)
        DO m = 1, npt
          W(m) = W(m)+t*(Glob_D(m, a, l)-Evalue*Glob_D(m+npt, a, l))
        ENDDO
      ENDDO
      t = Glob_c(ActiveIndex)
      t2 = t*t
      DO m = 1, npt
        Glob_WkGR((a-1)*npt+m) = TWO*t*W(m)
      ENDDO
      DO m = 1+Glob_ProcID, npt, Glob_NumOfProcs
        Glob_WkGR((a-1)*npt+m) = Glob_WkGR((a-1)*npt+m)- &
          t2*(Glob_D(m, a, ActiveIndex)- &
              Evalue*Glob_D(m+npt, a, ActiveIndex))
      ENDDO
    ENDDO

    IF (Glob_OverlapPenaltyAllowed) THEN
      CALL ComputeQOverlapPenaltyAndAddGradient(Glob_TotalOverlapPenalty, &
        Glob_WkGR, ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) THEN
        Evalue = HUGE(Evalue)
        Gradient = HUGE(Evalue)
        RETURN
      ENDIF
      Evalue = Evalue+Glob_TotalOverlapPenalty
    ENDIF
    CALL MPI_ALLREDUCE(Glob_WkGR, Gradient, NumActive*npt, MPI_WP, &
      MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    Glob_EnergyGBCounter = Glob_EnergyGBCounter+1

  END SUBROUTINE EnergyQB

  SUBROUTINE SolveEliminationGSEP(GSEPSolMethod, MatrixOrder, Evalue, ErrorCode)
    ! Subroutine SolveEliminationGSEP provides the common eigensolver boundary used
    ! by the elimination and separation BBOP routines. It constructs the initial
    ! factorization. After the basis change, Q cleanup drivers preserve this state
    ! through delete_symmetric or replace_symmetric; only G rebuilds and calls this
    ! routine a second time.
    !
    ! The G path preserves the historical DSYGVX layout: it materializes the upper
    ! triangle and restores Hamiltonian and overlap diagonals before calling LAPACK.
    ! The Q path must not do that. StoreHS has already placed normalized, unshifted H
    ! and S in their canonical lower triangles, including their physical diagonals.
    ! qrlinalg receives that representation directly and owns only its factors.
    !
    ! Arguments:
    CHARACTER(1), INTENT(IN) :: GSEPSolMethod
    INTEGER, INTENT(IN)      :: MatrixOrder
    REAL(wp), INTENT(OUT)    :: Evalue
    INTEGER, INTENT(OUT)     :: ErrorCode
    ! Local variables:
    ! LAPACK specifies IFAIL as dimension (N): on failure DSYGVX writes the
    ! indices of the eigenvectors that did not converge into it.
    INTEGER  :: i, j, IFAIL(Glob_HSLeadDim), NumOfEigvalsFound
    REAL(wp) :: EVs(1)

    Evalue = HUGE(Evalue)
    ErrorCode = Q_METHOD_INVALID_ARGUMENT

    SELECT CASE (GSEPSolMethod)
    CASE ('G')
      DO i = 1, MatrixOrder
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO
      DO i = 1, MatrixOrder
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO

      IF (Glob_ProcID == 0) THEN
        CALL DSYGVX(1, 'V', 'I', 'U', MatrixOrder, Glob_H, Glob_HSLeadDim, &
          Glob_S, Glob_HSLeadDim, ZERO, ZERO, Glob_WhichEigenvalue, &
          Glob_WhichEigenvalue, Glob_AbsTolForDSYGVX, NumOfEigvalsFound, &
          EVs, Glob_c, MatrixOrder, Glob_WorkForDSYGVX, &
          Glob_LWorkForDSYGVX, Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
        Evalue = EVs(1)
      ENDIF
      CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_c, MatrixOrder, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    CASE ('Q')
      ! There is no useful previous vector after a structural elimination or a
      ! random separation. ONE gives inverse iteration a deterministic nonzero
      ! starting vector and SolveQ returns an S-normalized coefficient vector.
      Glob_c(1:MatrixOrder) = ONE
      CALL PrepareQWorkspace(MatrixOrder, MatrixOrder, 1, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) THEN
        Q_Workspace%MatricesAreCanonical = .TRUE.
        CALL FactorizeQFresh(ErrorCode)
      ENDIF
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL SolveQ(Evalue, ErrorCode)
    ENDSELECT

  END SUBROUTINE SolveEliminationGSEP

  SUBROUTINE DeleteQMaskedFunctions(RemoveMask, ErrorCode)
    ! Subroutine DeleteQMaskedFunctions removes every marked canonical basis
    ! function from both the qrlinalg state and the physical matrix representation.
    ! Deletions are submitted in descending canonical order, so an original index
    ! continues to identify the same row and column after every preceding deletion.
    ! For r removed functions this costs O(r*n**2), while the historical cleanup
    ! path recalculated O(n**2) matrix elements and constructed another O(n**3)
    ! factorization even though every survivor-survivor element was unchanged.
    !
    ! The physical lower triangles are compacted only after all root-owned QR
    ! operations succeed. Their in-place ascending survivor copy is safe: each
    ! source row and column has an original index not smaller than its destination,
    ! and no write can destroy a matrix element needed by a later survivor. The
    ! upper triangles remain deliberately unspecified under the Q canonical-layout
    ! contract. Raw overlap norms and the current eigenvector are compacted by the
    ! same survivor map.
    !
    ! A public qrlinalg deletion cannot fail after the complete metadata preflight,
    ! but the recovery path is still explicit. If a future library implementation
    ! introduces a failure, the untouched physical matrices reconstruct the old
    ! factorization before this routine reports the original error.
    !
    ! Arguments:
    INTEGER, INTENT(IN)  :: RemoveMask(:)
    INTEGER, INTENT(OUT) :: ErrorCode
    ! Local variables:
    INTEGER              :: i, j, OldI, OldJ, OldOrder, NewOrder, RootError, RecoveryError
    INTEGER, ALLOCATABLE :: Survivor(:)

    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    OldOrder = Q_Workspace%MatrixOrder
    IF (Glob_GSEPSolutionMethod /= 'Q') RETURN
    IF (.NOT. Q_Workspace%MatricesAreCanonical) RETURN
    IF (.NOT. Q_Workspace%FactorsMatchMatrices) RETURN
    IF (OldOrder < 2) RETURN
    IF (SIZE(RemoveMask) /= OldOrder) RETURN
    NewOrder = COUNT(RemoveMask == 0)
    IF ((NewOrder < 1) .OR. (NewOrder >= OldOrder)) RETURN
    IF (.NOT. ALLOCATED(Glob_H)) RETURN
    IF (.NOT. ALLOCATED(Glob_S)) RETURN
    IF (.NOT. ALLOCATED(Glob_diagS)) RETURN
    IF (.NOT. ALLOCATED(Glob_c)) RETURN
    IF ((SIZE(Glob_H, 1) < OldOrder) .OR. (SIZE(Glob_H, 2) < OldOrder)) RETURN
    IF ((SIZE(Glob_S, 1) < OldOrder) .OR. (SIZE(Glob_S, 2) < OldOrder)) RETURN
    IF (SIZE(Glob_diagS) < OldOrder) RETURN
    IF (SIZE(Glob_c) < OldOrder) RETURN

    ALLOCATE(Survivor(NewOrder), STAT=RootError)
    IF (RootError /= 0) THEN
      ErrorCode = Q_METHOD_ALLOCATION_ERROR
      RETURN
    ENDIF
    j = 0
    DO i = 1, OldOrder
      IF (RemoveMask(i) == 0) THEN
        j = j+1
        Survivor(j) = i
      ENDIF
    ENDDO

    RootError = Q_METHOD_SUCCESS
    IF (Glob_ProcID == 0) THEN
      IF (.NOT. Q_Workspace%Factors%is_valid()) RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%order() /= OldOrder) RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%get_capacity() /= Q_Workspace%Capacity) &
        RootError = Q_METHOD_INVALID_STATE
      IF (Q_Workspace%Factors%get_shift() /= Glob_ApproxEnergy) &
        RootError = Q_METHOD_INVALID_STATE
      IF (RootError == Q_METHOD_SUCCESS) THEN
        DO i = OldOrder, 1, -1
          IF (RemoveMask(i) /= 0) THEN
            CALL Q_Workspace%Factors%delete_symmetric(i, RootError)
            IF (RootError /= Q_METHOD_SUCCESS) EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    CALL MPI_BCAST(RootError, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    IF (RootError /= Q_METHOD_SUCCESS) THEN
      Q_Workspace%FactorsMatchMatrices = .FALSE.
      CALL FactorizeQFresh(RecoveryError)
      IF (RecoveryError /= Q_METHOD_SUCCESS) THEN
        ErrorCode = RecoveryError
      ELSE
        ErrorCode = RootError
      ENDIF
      DEALLOCATE(Survivor)
      RETURN
    ENDIF

    DO i = 1, NewOrder
      OldI = Survivor(i)
      Glob_c(i) = Glob_c(OldI)
      Glob_diagS(i) = Glob_diagS(OldI)
      DO j = 1, i
        OldJ = Survivor(j)
        Glob_H(i, j) = Glob_H(OldI, OldJ)
        Glob_S(i, j) = Glob_S(OldI, OldJ)
      ENDDO
    ENDDO

    Q_Workspace%MatrixOrder = NewOrder
    Q_Workspace%NumActive = 0
    Q_Workspace%ActiveFunction = 0
    Q_Workspace%ActivePosition = 0
    Q_Workspace%MatrixParametersAreStored = .FALSE.
    Q_Workspace%TrialIsReady = .FALSE.
    Q_Workspace%TrialHasDerivatives = .FALSE.
    Q_Workspace%AcceptedPointIsStored = .FALSE.
    Q_Workspace%FactorsMatchMatrices = .TRUE.
    Q_Workspace%LastEigenpairResidual = HUGE(ONE)
    Q_Workspace%LastFactorResidual = HUGE(ONE)
    ErrorCode = Q_METHOD_SUCCESS
    DEALLOCATE(Survivor)

  END SUBROUTINE DeleteQMaskedFunctions


  SUBROUTINE BasisEnlQ(Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval, &
                       OverlapThreshold, LinCoeffThreshold, ErrorCode)
    !==================================================================
    ! Subroutine BasisEnlQ
    !==================================================================
    ! Enlarges the basis from Kstart-1 to Kstop functions, Kstep at a time,
    ! with the GSEP solved by the QR method ('Q'); the Q twin of BasisEnlG
    ! and BasisEnlI, with the same candidate generation, premultiplier-power
    ! scan, DRMNG optimization, acceptance tests (energy, pair overlap,
    ! linear coefficient, shape of the new functions), history and output.
    ! The basis keeps its canonical order: a candidate block is APPENDED to
    ! the QR factors of the accepted prefix (EvaluateQAppendedTrial), a
    ! rejected one is deleted again (TrimQFactors), and the optimizer
    ! replaces the block's columns (EnergyQA/EnergyQB through ApplyQTrial)
    ! instead of recomputing and refactorizing the whole problem.
    ! Arguments as BasisEnlG, plus ErrorCode: Q_METHOD_SUCCESS on return,
    ! Q_METHOD_INVALID_ARGUMENT when an argument is out of range (nothing
    ! is done then). A fatal Q status inside the step aborts the run.
    !==================================================================
    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval
    REAL(wp), INTENT(IN) :: OverlapThreshold, LinCoeffThreshold
    INTEGER, INTENT(OUT) :: ErrorCode
    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER  :: i, j, K, AttemptToGetGoodFunc, ii, jj, jbest
    INTEGER  :: npt, nfo, nfru, nfrup1, nvmax, nv
    INTEGER  :: ErrCode, NumOfFailures, NumOfEnergyEval, NumOfGradEval
    LOGICAL  :: IsSwapFileOK, IsEnergyImproved, ExitNeeded
    LOGICAL  :: IsOverlapBad, IsAnyLinCoeffBad, IsEnergyBad
    LOGICAL  :: IsShapeBad
    INTEGER  :: NumOfShapeRedraws
    REAL(wp) :: Cfac, ShapeSsum, ShapeSabs
    INTEGER  :: wbfu_t, wmu_t, wbfu, wmu, rgm1_counter, rgm2_counter
    REAL(wp) :: ms1, ms2
    REAL(wp) :: Evalue, E_init, E_best
    REAL(wp) :: t
    ! Largest legal premultiplier power: the greatest EVEN value not
    ! exceeding Glob_MaxPowerAllowed (GenerateTrialParam only ever produces
    ! even powers), so the power scan steps through 2,4,...,PWRMax.
    INTEGER, PARAMETER :: PWRMax = 2*(Glob_MaxPowerAllowed/2)
    ! Candidate block and the best candidate block found so far
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: ParSet, ParSetBest
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: ZIndSet, ZIndSetBest
    ! The optimization variables: the nonlinear parameters of the block
    ! laid out as one flat vector of nfo*npt elements
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_best, grad
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: ZIndOptSequence
    ! Arrays used by DRMNG
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D, V, V_init
    INTEGER, PARAMETER                  :: LIV = 60
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV
    INTEGER                             :: ALG
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    IF (Kstart /= Glob_CurrBasisSize+1) RETURN
    IF ((Kstart < 1) .OR. (Kstop < Kstart) .OR. (Kstep < 1) .OR. (NTrials < 1)) RETURN
    IF (MaxEnergyEval < 0) RETURN
    !==================================================================
    ! Announce the step
    !==================================================================
    wbfu_t = 0
    wmu_t = 0
    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlQ started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstart =', Kstart
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstop =', Kstop
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstep =', Kstep
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'OptimizationType =', OptimizationType
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval =', MaxEnergyEval
    ENDIF
    !==================================================================
    ! Global state this routine works under
    !==================================================================
    ! Glob_nfru/Glob_nfo/Glob_nfa define the block being added; the Q
    ! energy routines take the active functions from Q_Workspace instead,
    ! but the swap-file and matrix-element routines still read them.
    ! Overlap penalties are off: an overlap violation REJECTS the block.
    !------------------------------------------------------------------
    Glob_GSEPSolutionMethod = 'Q'
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_nfa = Kstart+Kstep
    Glob_nfo = Kstep
    Glob_HSLeadDim = Kstop
    Glob_HSBuffLen = Kstop*Kstep
    npt = Glob_npt
    nfo = Glob_nfo
    nvmax = Kstep*npt
    rgm1_counter = 0
    rgm2_counter = 0
    ms1 = ZERO
    ms2 = ZERO
    ! Reallocate arrays that contain the information about basis
    ! functions and optimization process to the final capacity; the
    ! accepted prefix keeps its canonical order.
    CALL ReallocateBasisFuncData(Kstop, Glob_CurrBasisSize)
    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! Q stores both diagonals inside Glob_H and Glob_S (lower triangles
    ! authoritative, upper ones unused), so there is no Glob_diagH and no
    ! DSYGVX workspace. Glob_D holds the derivatives of the matrix elements
    ! with respect to the nonlinear parameters of the active functions:
    ! 2*npt of them per (active function, basis function) pair.
    !------------------------------------------------------------------
    ALLOCATE(Glob_H(Kstop, Kstop))
    ALLOCATE(Glob_S(Kstop, Kstop))
    ALLOCATE(Glob_diagS(Kstop))
    ALLOCATE(Glob_c(Kstop))
    ALLOCATE(Glob_D(2*npt, Kstep, Kstop))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ! Allocate workspace for EnergyQB
    ALLOCATE(Glob_WkGR(nvmax))
    ! Local workspace
    ALLOCATE(ParSet(npt, Kstep))
    ALLOCATE(ParSetBest(npt, Kstep))
    ALLOCATE(ZIndSet(Kstep))
    ALLOCATE(ZIndSetBest(Kstep))
    ALLOCATE(x(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(ZIndOptSequence(Kstep))
    ! Arrays used by DRMNG; LV is the documented size of V plus one
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))
    ! The Q workspace: matrix order = the accepted prefix, capacity Kstop,
    ! at most Kstep functions active at a time
    CALL PrepareQWorkspace(Glob_CurrBasisSize, Kstop, Kstep, ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0125 in BasisEnlQ: Q workspace cannot be allocated'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next (1 energy,
    ! 2 gradient) and is called again, which suits an energy that is a
    ! collective operation over all processes. The settings are those of
    ! BasisEnlG, so a G/Q comparison changes the eigensolver only.
    !------------------------------------------------------------------
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): evaluation and iteration limits, set out of the way
    ! because the budget is enforced by MaxEnergyEval below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed
    !==================================================================
    ! Initial state and energy
    !==================================================================
    ! Only the accepted prefix can be in the swap file. The unused
    ! capacity is zeroed before the full-capacity broadcast of the shared
    ! swap reader, then ONE fresh factorization of the prefix is built.
    ! With Kstart==1 there is no basis yet: Glob_CurrEnergy is primed with
    ! HUGE() and the first candidate to produce a finite energy wins.
    !------------------------------------------------------------------
    Glob_H = ZERO
    Glob_S = ZERO
    Glob_diagS = ZERO
    Glob_c = ONE
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)
    Q_Workspace%MatricesAreCanonical = .TRUE.
    ErrCode = Q_METHOD_SUCCESS
    IF (Glob_CurrBasisSize > 0) THEN
      IF (.NOT. IsSwapFileOK) THEN
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem...'
        CALL ComputeMatElem(1, Glob_CurrBasisSize)
      ELSE
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Solving eigenvalue problem...'
      ENDIF
      CALL FactorizeQFresh(ErrCode)
      IF (ErrCode == Q_METHOD_SUCCESS) CALL SolveQ(Glob_CurrEnergy, ErrCode)
    ELSE
      Glob_CurrEnergy = HUGE(Glob_CurrEnergy)
      CALL TrimQFactors(0, ErrCode)
    ENDIF
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, '(1x,a,1x,i0)') &
        'Error EC0126 in BasisEnlQ: initial Q state cannot be constructed, status', ErrCode
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy
    K = Kstart-1
    !==================================================================
    ! MAIN LOOP - one block of up to Kstep functions per iteration
    !==================================================================
    DO WHILE (K < Kstop)
      !--------------------------------------------------------------
      ! Size and place the block: nfru accepted functions stay, nfo are
      ! added, K is the new basis size (the last block may be short)
      !--------------------------------------------------------------
      nfru = K
      nfo = MIN(Kstep, Kstop-K)
      K = K+nfo
      CALL linalg_setparam(K)  ! reset linalg flags to account for changes in the basis size
      Glob_nfa = K
      Glob_nfru = nfru
      Glob_nfo = nfo
      nfrup1 = nfru+1
      nv = nfo*npt
      E_init = Glob_CurrEnergy
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Current basis size is', Glob_CurrBasisSize
        IF (nfo > 1) THEN
          WRITE(*, '(1x,a,1x,i0,a,i0)') 'Selecting functions', nfrup1, '-', K
        ELSE
          WRITE(*, '(1x,a,1x,i0)') 'Selecting function', K
        ENDIF
      ENDIF
      !==================================================================
      ! ACCEPTANCE LOOP - regenerate the block until it is acceptable
      !==================================================================
      ! The flags are primed so that the loop always runs at least once.
      ! It ends when the block passes the energy, overlap and linear
      ! coefficient tests, or when the attempt budget runs out - in which
      ! case the run stops (EC0130), as in the reference. A block in which the
      ! Young operator nearly annihilates a new function (IsShapeBad) is
      ! redrawn without limit and without spending the attempt budget.
      ! Every rejected candidate suffix is removed from the QR factors by
      ! the next EvaluateQAppendedTrial, so trials never pile up.
      !------------------------------------------------------------------
      IsOverlapBad = .TRUE.
      IsAnyLinCoeffBad = .TRUE.
      IsEnergyBad = .FALSE.
      IsShapeBad = .FALSE.
      NumOfShapeRedraws = 0
      AttemptToGetGoodFunc = 1
      DO WHILE (((IsOverlapBad .OR. IsAnyLinCoeffBad .OR. IsEnergyBad) .AND. &
                 (AttemptToGetGoodFunc <= Glob_BadOverlapOrLinCoeffLim)) .OR. IsShapeBad)
        !------------------------------------------------------------------
        ! Step 1: stochastic selection of the block
        !------------------------------------------------------------------
        ! GenerateTrialParam runs on rank 0 (it consumes the random stream)
        ! and the result is broadcast; every rank evaluates the SAME
        ! candidate, appended to the accepted prefix. A candidate the Q solve
        ! cannot handle is counted, not fatal; only an excessive FRACTION of
        ! failures is.
        !------------------------------------------------------------------
        NumOfFailures = 0
        IsEnergyImproved = .FALSE.
        wbfu = 0
        wmu = 0
        DO i = 1, NTrials
          IF (Glob_ProcID == 0) CALL GenerateTrialParam(nfo, ParSet, ZIndSet, wbfu_t, wmu_t)
          CALL MPI_BCAST(ParSet, npt*nfo, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_BCAST(ZIndSet, nfo, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          Glob_NonlinParam(1:npt, nfrup1:K) = ParSet(1:npt, 1:nfo)
          Glob_PWR(nfrup1:K) = ZIndSet(1:nfo)
          CALL EvaluateQAppendedTrial(nfru, K, Evalue, ErrCode)
          IF (ErrCode == Q_METHOD_SUCCESS) THEN
            IF (Evalue < Glob_CurrEnergy) THEN
              Glob_CurrEnergy = Evalue
              ParSetBest(1:npt, 1:nfo) = ParSet(1:npt, 1:nfo)
              ZIndSetBest(1:nfo) = ZIndSet(1:nfo)
              IsEnergyImproved = .TRUE.
              wbfu = wbfu_t
              wmu = wmu_t
            ENDIF
          ELSE
            NumOfFailures = NumOfFailures+1
            IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a,1x,i0,1x,a,1x,i0)') &
              'Warning WC0114 in BasisEnlQ: candidate', i, 'failed with Q status', ErrCode
          ENDIF
        ENDDO
        ! Too many candidates the solver could not handle at all: the
        ! basis is in a state the eigensolver cannot work with, and more
        ! trials will not fix it.
        IF (NumOfFailures*ONE/NTrials > Glob_MaxFracOfTrialFailsAllowed) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0127 in BasisEnlQ: the number of eigenvalue problem solution failures'
            WRITE(*, *) 'in random selection process exceeded limit'
            WRITE(*, '(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials, '%'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        ! The candidates were all solvable but none lowered the energy.
        IF (.NOT. (IsEnergyImproved)) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0128 in BasisEnlQ: random selection did not result'
            WRITE(*, *) 'in any energy improvement'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        ! Put the winning candidate back in place: the last candidate
        ! evaluated is not in general the best one. A successful append
        ! also records the matrix parameters the replacement transactions
        ! of EnergyQA/EnergyQB start from.
        Glob_NonlinParam(1:npt, nfrup1:K) = ParSetBest(1:npt, 1:nfo)
        Glob_PWR(nfrup1:K) = ZIndSetBest(1:nfo)
        CALL EvaluateQAppendedTrial(nfru, K, Glob_CurrEnergy, ErrCode)
        IF (ErrCode /= Q_METHOD_SUCCESS) THEN
          IF (Glob_ProcID == 0) WRITE(*, '(1x,a,1x,i0)') &
            'Error EC0129 in BasisEnlQ: selected Q candidate cannot be reconstructed, status', ErrCode
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        IF (Glob_ProcID == 0) THEN
          WRITE (*, '(1x,a)', ADVANCE='no') 'E='
          CALL writereal(6, Glob_CurrEnergy)
          IF (Verbose >= 2) WRITE (*, '(5x,a,1x,i0)') 'prototype function is', wbfu
          DO i = 1, nfo
            WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', ZIndSetBest(i)
            CALL writerealarradv(6, ParSetBest(1:npt, i), npt)
          ENDDO
          IF (Verbose >= 1) WRITE (*, *) 'Optimizing nonlinear parameters'
        ENDIF
        !------------------------------------------------------------------
        ! Prime the best point found
        !------------------------------------------------------------------
        ! Every branch of the SELECT below must leave x_best defined, because
        ! the energy and the linear coefficients are recomputed at x_best
        ! right after it; the randomly selected point is also what
        ! OptimizationType 0 keeps.
        !------------------------------------------------------------------
        DO i = 1, nfo
          x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
        ENDDO
        E_best = Glob_CurrEnergy
        x_best(1:nv) = x(1:nv)
        NumOfEnergyEval = 0
        NumOfGradEval = 0
        SELECT CASE (OptimizationType)
        !------------------------------------------------------------------
        ! OptimizationType 0: keep the randomly selected block
        !------------------------------------------------------------------
        CASE (0)
        ! x and x_best already hold the selected point - nothing to do.
        !------------------------------------------------------------------
        ! OptimizationType 1: optimize powers, then nonlinear parameters
        !------------------------------------------------------------------
        CASE (1)
          !------------------------------------------------------------------
          ! Step 2: premultiplier powers, one function at a time
          !------------------------------------------------------------------
          ! Each function of the block is tried with every legal EVEN power
          ! 2,4,...,PWRMax and the best is kept; the functions are visited
          ! in random order. The power is part of the basis function but not
          ! a DRMNG variable, so every alternative is tested with the same
          ! append/delete transaction as a random candidate, and the winning
          ! powers are rebuilt explicitly at the end (the nonlinear
          ! parameters alone cannot tell the Q transaction that this metadata
          ! changed). Skipped when the power is pinned by Glob_IsIndexFixed.
          !------------------------------------------------------------------
          NumOfFailures = 0
          IF (.NOT. Glob_IsIndexFixed) THEN
            ! Generate a random sequence which will define the order in which
            ! Z-indices should be optimized (one index at a time)
            CALL GenerateRndIntSeq(nfo, ZIndOptSequence)
            DO i = 1, nfo
              ii = ZIndOptSequence(i)
              j = Glob_PWR(nfru+ii)
              jbest = j
              DO jj = 2, PWRMax, 2
                IF (jj /= j) THEN
                  Glob_PWR(nfru+ii) = jj
                  CALL EvaluateQAppendedTrial(nfru, K, Evalue, ErrCode)
                  IF (ErrCode /= Q_METHOD_SUCCESS) THEN
                    ! Restore the best power known so far before carrying on,
                    ! so a failed trial never leaves a bad power behind.
                    NumOfFailures = NumOfFailures+1
                    Glob_PWR(nfru+ii) = jbest
                    IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
                      IF (Glob_ProcID == 0) THEN
                        WRITE(*, *) 'Error EC0128 in BasisEnlQ: number of failures in energy calculations'
                        WRITE(*, *) 'during the optimization of Z-indicies exceeded limit'
                      ENDIF
                      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
                    ENDIF
                  ELSE
                    IF (Evalue < Glob_CurrEnergy) THEN
                      Glob_CurrEnergy = Evalue
                      jbest = jj
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              Glob_PWR(nfru+ii) = jbest
            ENDDO
            CALL EvaluateQAppendedTrial(nfru, K, Glob_CurrEnergy, ErrCode)
            IF (ErrCode /= Q_METHOD_SUCCESS) THEN
              IF (Glob_ProcID == 0) WRITE(*, '(1x,a,1x,i0)') &
                'Error EC0129 in BasisEnlQ: optimized Z-indices cannot be reconstructed, status', ErrCode
              CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
            ENDIF
          ENDIF
          !------------------------------------------------------------------
          ! Step 3: nonlinear parameters, with DRMNG
          !------------------------------------------------------------------
          ! IV and V are restored from the copies made before the main loop,
          ! so each block starts the minimizer from its documented default
          ! state. Every variable gets the same scale t, taken from the
          ! RELATIVE energy gain of the block so far, floored at
          ! 10000*epsilon; at the very start of a basis t=1.
          !------------------------------------------------------------------
          IV(1:LIV) = IV_init(1:LIV)
          V(1:LV) = V_init(1:LV)
          DO i = 1, nfo
            x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
          ENDDO
          IF (nfru >= nfo) THEN
            t = MAX(ABS((E_init-Glob_CurrEnergy))/(ABS(E_init)+ABS(Glob_CurrEnergy)), &
                    10000*EPSILON(Glob_CurrEnergy))
          ELSE
            t = ONE
          ENDIF
          D(1:nv) = t
          ExitNeeded = .FALSE.
          NumOfFailures = 0
          IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.
          E_best = Glob_CurrEnergy
          x_best(1:nv) = x(1:nv)
          !------------------------------------------------------------------
          ! The reverse-communication loop
          !------------------------------------------------------------------
          ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it
          ! wants: 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its
          ! evaluation limit. A failed evaluation is reported with IV(2)=1
          ! (TOOBIG), which makes DRMNG shrink the step. The best point is
          ! tracked here because the last point DRMNG visits is not
          ! necessarily the lowest.
          !------------------------------------------------------------------
          DO WHILE (.NOT. (ExitNeeded))
            IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
            CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            SELECT CASE (IV(1))
            CASE (1)  ! Only energy is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              Evalue = EnergyQA(.TRUE., ErrCode)
              NumOfEnergyEval = NumOfEnergyEval+1
              IF (ErrCode /= Q_METHOD_SUCCESS) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
              ELSE
                Glob_CurrEnergy = Evalue
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nv) = x(1:nv)
                ENDIF
              ENDIF
            CASE (2)  ! Only gradient is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              CALL EnergyQB(Evalue, grad, .TRUE., ErrCode)
              NumOfGradEval = NumOfGradEval+1
              IF (ErrCode /= Q_METHOD_SUCCESS) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
              ELSE
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nv) = x(1:nv)
                ENDIF
              ENDIF
            CASE (3:8)  ! Some kind of convergence has been reached
              ExitNeeded = .TRUE.
            CASE (9:10)  ! Function evaluation limit has been reached.
              ! This is never supposed to happen because we
              ! count the number of function evaluations ourselves.
              ExitNeeded = .TRUE.
            CASE DEFAULT
              ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
              ! and >=14 for a bad input. None of those match a case above,
              ! so without this the loop would call DRMNG again for ever.
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                  'Warning WC0137 in BasisEnlQ: DRMNG returned IV(1) =', IV(1)
                IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
              ENDIF
              ExitNeeded = .TRUE.
            ENDSELECT
            ! A warning, not an abort: the best point found so far is still
            ! usable, and the acceptance test below decides what to do with
            ! the block.
            IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                  'Warning WC0112 in BasisEnlQ: number of failures in energy or gradient', &
                  'calculations during the optimization of nonlinear parameters', &
                  'reached the limit of', Glob_MaxEnergyFailsAllowed
              ENDIF
            ENDIF
            IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.
          ENDDO
        !------------------------------------------------------------------
        ! Anything else is a programming error, not an input choice
        !------------------------------------------------------------------
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0131 in BasisEnlQ: unsupported value of OptimizationType', OptimizationType
            IF (Verbose >= 1) WRITE(*, *) 'Allowed values are 0 (no optimization) and 1 (powers and nonlinear parameters)'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
        ENDSELECT  ! (OptimizationType)
        !------------------------------------------------------------------
        ! Step 4a: re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! The Q solve always returns the eigenvector, so EnergyQAM is one
        ! replacement transaction at x_best; it makes the matrices, the
        ! factors, the energy and Glob_c describe the same point. A failure
        ! is NOT fatal: the block is rejected and regenerated, the same
        ! response as a bad overlap.
        !------------------------------------------------------------------
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, nfru+i) = x_best((i-1)*npt+1:i*npt)
        ENDDO
        IsEnergyBad = .FALSE.
        Evalue = EnergyQAM(.TRUE., ErrCode)
        IF (ErrCode == Q_METHOD_SUCCESS) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          IsEnergyBad = .TRUE.
          Glob_CurrEnergy = E_init
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0113 in BasisEnlQ: failed to evaluate energy after the optimization', &
              'of nonlinear parameters. Generated basis function(s) are rejected', &
              'and a new attempt to generate them will be made'
          ENDIF
        ENDIF
        !------------------------------------------------------------------
        ! Step 4b: pair overlaps
        !------------------------------------------------------------------
        ! Every pair involving a NEW function is checked against the whole
        ! basis; only the canonical lower triangle of Glob_S is read.
        ! Glob_CurrEnergy is rolled back to E_init so the next attempt
        ! measures its improvement against the basis before this block.
        ! Disabled by OverlapThreshold <= 0.
        !------------------------------------------------------------------
        IsOverlapBad = .FALSE.
        IF ((ErrCode == Q_METHOD_SUCCESS) .AND. (OverlapThreshold > ZERO)) THEN
          ii = 0
          DO i = nfrup1, K
            DO j = 1, i-1
              IF (ABS(QCanonicalMatrixElement(Glob_S, i, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                Glob_CurrEnergy = E_init
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0110: overlap of the following functions exceeds threshold. ', &
                      'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') ii, ':', i, j, '    S='
                  CALL writerealadv(6, QCanonicalMatrixElement(Glob_S, i, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        !------------------------------------------------------------------
        ! Step 4c: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..K, not just the new functions:
        ! adding a function can blow up the coefficient of an old one, and
        ! that is exactly the near-linear-dependence this test is meant to
        ! catch. Disabled by LinCoeffThreshold <= 0.
        !------------------------------------------------------------------
        IsAnyLinCoeffBad = .FALSE.
        IF ((ErrCode == Q_METHOD_SUCCESS) .AND. (LinCoeffThreshold > ZERO)) THEN
          ii = 0
          DO i = 1, K
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              Glob_CurrEnergy = E_init
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0111: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', i, '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        !------------------------------------------------------------------
        ! Step 4d: shape of the new functions after the Young operator
        !------------------------------------------------------------------
        ! C = sum|c_k S_k| / |<phi|Y+Y|phi>| of each new function. Above
        ! Glob_MaxSelfOverlapCancel the operator has almost annihilated the
        ! function: what survives is round-off, which the energy test cannot
        ! tell from a genuine improvement. Such a block is redrawn as often
        ! as necessary; these redraws do not count against the attempt budget.
        !------------------------------------------------------------------
        IsShapeBad = .FALSE.
        ii = 0
        DO i = nfrup1, K
          Cfac = SelfOverlapCancellation(Glob_PWR(i), Glob_NonlinParam(1:npt, i), ShapeSsum, ShapeSabs)
          IF (Cfac > Glob_MaxSelfOverlapCancel) THEN
            ii = ii+1
            IsShapeBad = .TRUE.
            Glob_CurrEnergy = E_init
            IF (Glob_ProcID == 0) THEN
              IF ((ii == 1) .AND. (Verbose >= 1)) THEN
                WRITE(*, '(1x,a,a,es9.2,a)') 'Warning WC0116 in BasisEnlQ: the Young operator nearly annihilates ', &
                  'the following function(s), C > ', Glob_MaxSelfOverlapCancel, &
                  '. Generated basis function(s) are rejected and a new attempt to generate them will be made'
              ENDIF
              IF (Verbose >= 1) WRITE(*, '(1x,i6,a1,i6,a8,i5,a5,es10.3,a10,es10.3)') ii, ':', i, '   power', &
                Glob_PWR(i), '   C=', Cfac, '   S_raw=', ShapeSsum
            ENDIF
          ENDIF
        ENDDO
        IF (IsShapeBad) THEN
          NumOfShapeRedraws = NumOfShapeRedraws+1
          IF ((Glob_ProcID == 0) .AND. (Verbose >= 1) .AND. (MOD(NumOfShapeRedraws, 50) == 0)) &
            WRITE(*, '(1x,a,i0,a)') 'BasisEnlQ: ', NumOfShapeRedraws, ' blocks redrawn so far because of the shape test'
        ENDIF
        IF (.NOT. IsShapeBad) AttemptToGetGoodFunc = AttemptToGetGoodFunc+1
      ENDDO  ! acceptance loop
      ! The attempt budget ran out: stop, as the reference does. The last
      ! accepted block is already in the data file (SaveResults per block).
      IF (IsOverlapBad .OR. IsAnyLinCoeffBad .OR. IsEnergyBad) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'Error EC0130 in BasisEnlQ: unable to construct an acceptable candidate block'
          WRITE(*, '(1x,a,1x,i0,1x,a)') 'No block passed the energy, overlap and linear coefficient tests in', &
            Glob_BadOverlapOrLinCoeffLim, 'attempts'
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF
      !==================================================================
      ! Report the accepted block
      !==================================================================
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) THEN
        WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
        ENDIF
        WRITE (*, *) 'E=', Glob_CurrEnergy
        DO i = 1, nfo
          WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', Glob_PWR(nfru+i)
          CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
        ENDDO
      ENDIF
      !==================================================================
      ! Generator statistics
      !==================================================================
      ! Distance of the accepted parameters from their prototype,
      ! accumulated per generator method and averaged at the end of the
      ! run; rank 0 only. Skipped for the first blocks (nfru<=nfo), where
      ! prototype and new window overlap.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN
        IF (wmu == 1) THEN
          rgm1_counter = rgm1_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms1 = ms1+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF (wmu == 2) THEN
          rgm2_counter = rgm2_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms2 = ms2+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      !==================================================================
      ! Commit the block
      !==================================================================
      ! The history of a newly added function starts empty. SaveResults
      ! runs on rank 0 only and after every block, so an interrupted run
      ! can be resumed from the last completed block. The canonical order
      ! is the user-visible order, so no sorting workspace is needed.
      !------------------------------------------------------------------
      Glob_CurrBasisSize = K
      DO i = 1, nfo
        Glob_History(nfru+i)%Energy = Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone = 0
        Glob_History(nfru+i)%InitFuncAtLastStep = 0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt = 0
        Glob_FuncNum(nfru+i) = nfru+i
      ENDDO
      IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')
    ENDDO  ! main loop
    !==================================================================
    ! Hand H and S to the next step and release everything
    !==================================================================
    CALL StoreMatricesInSwapFile()
    CALL ClearQWorkspace()
    ! Deallocate arrays used by DRMNG
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)
    ! Deallocate workspace
    DEALLOCATE(ZIndOptSequence)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x)
    DEALLOCATE(ZIndSetBest)
    DEALLOCATE(ZIndSet)
    DEALLOCATE(ParSetBest)
    DEALLOCATE(ParSet)
    ! Deallocate workspace for EnergyQB
    DEALLOCATE(Glob_WkGR)
    ! Deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)
    !==================================================================
    ! Closing summary
    !==================================================================
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *) 'Random selection statistics:'
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 1 of generating basis functions was used', rgm1_counter, 'times'
      IF ((rgm1_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms1/(npt*rgm1_counter)
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 2 of generating basis functions was used', rgm2_counter, 'times'
      IF ((rgm2_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms2/(npt*rgm2_counter)
      IF (Verbose >= 2) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlQ finished'
    ENDIF
    ErrorCode = Q_METHOD_SUCCESS
  END SUBROUTINE BasisEnlQ


  SUBROUTINE OptCycleQ(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                       NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold, &
                       SavingFreq, ErrorCode)
    !==================================================================
    ! Subroutine OptCycleQ
    !==================================================================
    ! Improves an EXISTING basis (no functions added): a window of
    ! NumOfFuncToOpt functions sweeps FuncBegin..FuncEnd, advancing by
    ! NumOfFuncToShift per step, NumCycles times, optimizing the nonlinear
    ! parameters in the window with the QR method ('Q') and DRMNG; the Q
    ! twin of OptCycleG and OptCycleI, with their scheduling, acceptance
    ! tests, failure limits, saving and history. The basis is NEVER
    ! permuted: the window is an explicit active map (SetQActiveFunctions)
    ! and the optimizer replaces its columns in the QR factors
    ! (EnergyQA/EnergyQB through ApplyQTrial). G reverses the window before
    ! moving it to the end; the descending active map gives DRMNG the same
    ! block order. A step that ends with an unusable energy, an overlap
    ! above OverlapThreshold or a linear coefficient above
    ! LinCoeffThreshold is UNDONE (parameters back to x_init).
    ! Arguments as OptCycleG (K is NOT REFERENCED; MaxEnergyEval <= 0 ->
    ! Glob_MaxFuncEvalForCyclOpt; SavingFreq < 1 -> 1), plus ErrorCode:
    ! Q_METHOD_SUCCESS on return, Q_METHOD_INVALID_ARGUMENT when an
    ! argument is out of range (nothing is done then).
    !==================================================================
    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: K                   ! NOT REFERENCED - see the header
    INTEGER, INTENT(IN)  :: FuncBegin, FuncEnd  ! range of functions to sweep
    INTEGER, INTENT(IN)  :: NumOfFuncToOpt      ! window size
    INTEGER, INTENT(IN)  :: NumOfFuncToShift    ! how far the window advances per step
    INTEGER, INTENT(IN)  :: NumCycles           ! sweeps over the range
    INTEGER, INTENT(IN)  :: MaxEnergyEval       ! evaluations per step, <=0 = use default
    REAL(wp), INTENT(IN) :: OverlapThreshold    ! pair-overlap rejection, <=0 = off
    REAL(wp), INTENT(IN) :: LinCoeffThreshold   ! linear-coefficient rejection, <=0 = off
    INTEGER, INTENT(IN)  :: SavingFreq          ! save every SavingFreq steps, <1 = 1
    INTEGER, INTENT(OUT) :: ErrorCode           ! Q status of the step
    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i, j, a          ! loop counters
    INTEGER :: ii               ! counts the violations an acceptance test finds
    INTEGER :: CurrCycle        ! sweep number, 1..NumCycles
    INTEGER :: CurrFunc         ! first function of the current window
    INTEGER :: CurrFuncBegin    ! where this cycle's sweep starts
    INTEGER :: totsteps         ! steps done in this call; drives the saving
    INTEGER :: cbs              ! Glob_CurrBasisSize, the basis worked on
    INTEGER :: npt              ! Glob_npt, nonlinear parameters per function
    INTEGER :: nfo              ! functions in the current window
    INTEGER :: nv               ! nfo*npt, optimization variables this step
    INTEGER :: nvmax            ! NumOfFuncToOpt*npt, the largest nv can be
    INTEGER :: ActiveI          ! canonical index of an active function
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ActiveFunction  ! the window as an active map
    REAL(wp) :: Evalue  ! value returned by the latest solve
    REAL(wp) :: E_best  ! lowest energy this step has reached
    REAL(wp) :: E_prev  ! energy before the step, to fall back on
    REAL(wp) :: t       ! scale factor for the DRMNG vector D
    LOGICAL :: IsSwapFileOK      ! .true. when H and S came from the swap file
    LOGICAL :: ExitNeeded        ! ends the reverse-communication loop
    LOGICAL :: LastIter          ! this is the final step of the cycle
    LOGICAL :: IsOverlapBad      ! an overlap exceeded OverlapThreshold
    LOGICAL :: IsAnyLinCoeffBad  ! a coefficient exceeded LinCoeffThreshold
    INTEGER :: ErrCode          ! Q status of the latest transaction
    INTEGER :: NumOfFailures    ! failed solves in this step
    INTEGER :: NumOfEnergyEval  ! energy evaluations in this step
    INTEGER :: NumOfGradEval    ! gradient evaluations in this step
    INTEGER :: MaxEvalToUse  ! energy evaluations allowed per step
    INTEGER :: SaveEvery     ! save every SaveEvery steps
    ! The optimization variables: the nonlinear parameters of the window,
    ! flattened. x_init is the point the step started from.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_init, x_best, grad
    ! Arrays and settings used by DRMNG
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D                      ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V, V_init              ! work array and its copy
    INTEGER, PARAMETER                  :: LIV = 60               ! length of IV
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV                     ! length of V
    INTEGER                             :: ALG                    ! 2 = unconstrained minimization
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    cbs = Glob_CurrBasisSize
    IF (K /= cbs) RETURN
    IF ((FuncBegin < 1) .OR. (FuncEnd < FuncBegin) .OR. (FuncEnd > cbs)) RETURN
    IF ((NumOfFuncToOpt < 1) .OR. (NumOfFuncToShift < 1)) RETURN
    IF (NumCycles < 0) RETURN
    !==================================================================
    ! Global state this routine works under
    !==================================================================
    Glob_GSEPSolutionMethod = 'Q'
    Glob_OverlapPenaltyAllowed = .FALSE.
    npt = Glob_npt
    nvmax = NumOfFuncToOpt*npt
    Glob_HSLeadDim = cbs
    Glob_HSBuffLen = cbs*NumOfFuncToOpt
    Glob_nfa = cbs
    !==================================================================
    ! Nothing to do?
    !==================================================================
    ! Both conditions have to hold: the required number of cycles is
    ! done AND the last of them reached the end of the range.
    !------------------------------------------------------------------
    IF ((Glob_History(cbs)%CyclesDone >= NumCycles) .AND. &
        (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleQ started'
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
        IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
          FuncBegin, '-', FuncEnd, 'is already completed'
        IF (Verbose >= 1) WRITE(*, *) 'Exiting OptCycleQ...'
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleQ finished'
      ENDIF
      ErrorCode = Q_METHOD_SUCCESS
      RETURN
    ENDIF
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleQ started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
      WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin, '-', FuncEnd, 'will be performed'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval', MaxEnergyEval
    ENDIF
    !==================================================================
    ! Turn the two input limits into usable values
    !==================================================================
    ! MaxEnergyEval <= 0 would make the DRMNG loop exit before its first
    ! iteration and leave every function untouched, so it falls back on
    ! Glob_MaxFuncEvalForCyclOpt; SavingFreq <= 0 would reach a MOD by zero.
    !------------------------------------------------------------------
    IF (MaxEnergyEval > 0) THEN
      MaxEvalToUse = MaxEnergyEval
    ELSE
      MaxEvalToUse = Glob_MaxFuncEvalForCyclOpt
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a,1x,i0)') &
          'Warning WC0131 in OptCycleQ: MaxEnergyEval is', MaxEnergyEval, ',', &
          'using the default limit of', Glob_MaxFuncEvalForCyclOpt
      ENDIF
    ENDIF
    SaveEvery = MAX(SavingFreq, 1)
    IF ((SavingFreq < 1) .AND. (Glob_ProcID == 0)) THEN
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a)') &
        'Warning WC0132 in OptCycleQ: SavingFreq is', SavingFreq, ',', &
        'results will be saved after every step'
    ENDIF
    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! Q stores both diagonals inside Glob_H and Glob_S; no Glob_diagH, no
    ! DSYGVX workspace, no permutation or sorting workspace.
    !------------------------------------------------------------------
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_D(2*npt, NumOfFuncToOpt, cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ! Allocate workspace for EnergyQB
    ALLOCATE(Glob_WkGR(nvmax))
    ! Allocate arrays used by DRMNG and the active-index schedule
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))
    ALLOCATE(x(nvmax))
    ALLOCATE(x_init(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(ActiveFunction(NumOfFuncToOpt))
    CALL PrepareQWorkspace(cbs, cbs, NumOfFuncToOpt, ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0150 in OptCycleQ: Q workspace cannot be allocated'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! The settings are those of OptCycleG, so a G/Q comparison changes the
    ! eigensolver only; set once, outside the main loop.
    !------------------------------------------------------------------
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): evaluation and iteration limits, set out of the way
    ! because the budget is enforced by MaxEvalToUse below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed
    !==================================================================
    ! Initial state and energy
    !==================================================================
    ! The canonical matrices come from the swap file when it is usable;
    ! otherwise one full assembly. After that every optimizer evaluation
    ! changes only the active rows and columns.
    !------------------------------------------------------------------
    Glob_H = ZERO
    Glob_S = ZERO
    Glob_diagS = ZERO
    Glob_c = ONE
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)
    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) &
        'Computing matrix elements and constructing fresh QR factors...'
      CALL ComputeMatElem(1, cbs)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Constructing fresh QR factors...'
    ENDIF
    Q_Workspace%MatricesAreCanonical = .TRUE.
    CALL FactorizeQFresh(ErrCode)
    IF (ErrCode == Q_METHOD_SUCCESS) CALL SolveQ(Glob_CurrEnergy, ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, '(1x,a,1x,i0)') &
        'Error EC0151 in OptCycleQ: initial Q energy cannot be computed, status', ErrCode
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy
    !------------------------------------------------------------------
    ! Normalize the restart position, exactly as in OptCycleG. Function
    ! numbers stay in their input positions, so no preliminary or restart
    ! permutation is required.
    !------------------------------------------------------------------
    IF (Glob_History(cbs)%InitFuncAtLastStep < FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
    IF (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd) THEN
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
    ENDIF
    !==================================================================
    ! OUTER LOOP - one pass over FuncBegin..FuncEnd per cycle
    !==================================================================
    totsteps = 0
    DO CurrCycle = Glob_History(cbs)%CyclesDone+1, NumCycles
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Cycle', CurrCycle, 'began'
      ENDIF
      CurrFuncBegin = Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      !==================================================================
      ! INNER LOOP - one window per step
      !==================================================================
      DO CurrFunc = CurrFuncBegin, FuncEnd, NumOfFuncToShift
        totsteps = totsteps+1
        ! The window is short only on the final step of a cycle, when
        ! fewer than NumOfFuncToOpt functions are left in the range.
        nfo = MIN(FuncEnd-CurrFunc+1, NumOfFuncToOpt)
        Glob_nfo = nfo
        Glob_nfru = cbs-nfo
        nv = nfo*npt
        !------------------------------------------------------------------
        ! Install the window as the active map
        !------------------------------------------------------------------
        ! G reverses the requested range before moving it to the trailing
        ! block; this descending map gives DRMNG the same block order while
        ! every physical basis function remains at its canonical index.
        !------------------------------------------------------------------
        DO a = 1, nfo
          ActiveFunction(a) = CurrFunc+nfo-a
        ENDDO
        CALL SetQActiveFunctions(ActiveFunction(1:nfo), ErrCode)
        IF (ErrCode == Q_METHOD_SUCCESS) CALL CaptureQMatrixParameters(ErrCode)
        IF (ErrCode /= Q_METHOD_SUCCESS) THEN
          IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0152 in OptCycleQ: active Q transaction cannot be initialized'
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *)
          IF (nfo > 1) THEN
            IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Optimizing functions', CurrFunc, '-', CurrFunc+nfo-1
          ELSE
            WRITE(*, '(1x,a,1x,i0)') 'Optimizing function', CurrFunc
          ENDIF
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 2) WRITE (*, *) 'Nonlinear parameters before optimization:'
            DO a = 1, nfo
              ActiveI = ActiveFunction(a)
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(ActiveI), ':', Glob_PWR(ActiveI)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, ActiveI), npt)
            ENDDO
          ENDIF
        ENDIF
        !------------------------------------------------------------------
        ! Set up the step
        !------------------------------------------------------------------
        ! x holds the window's parameters flattened, x_init a copy to
        ! restore from if the step is rejected. D is the DRMNG scale vector:
        ! the 1/(cbs^2*sqrt(cbs)) form shrinks the steps as the basis grows,
        ! floored at 10000*epsilon.
        !------------------------------------------------------------------
        IV(1:LIV) = IV_init(1:LIV)
        V(1:LV) = V_init(1:LV)
        DO a = 1, nfo
          ActiveI = ActiveFunction(a)
          x((a-1)*npt+1:a*npt) = Glob_NonlinParam(1:npt, ActiveI)
          x_init((a-1)*npt+1:a*npt) = Glob_NonlinParam(1:npt, ActiveI)
        ENDDO
        t = MAX(ONE/(cbs*cbs*SQRT(ONE*cbs)), 10000*EPSILON(Glob_CurrEnergy))
        D(1:nv) = t
        ExitNeeded = .FALSE.
        NumOfFailures = 0
        NumOfEnergyEval = 0
        NumOfGradEval = 0
        IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.
        E_best = Glob_CurrEnergy
        x_best(1:nv) = x(1:nv)
        ! Remember the energy of the unchanged parameters, needed in case
        ! the optimization of this window has to be abandoned
        E_prev = Glob_CurrEnergy
        !------------------------------------------------------------------
        ! The reverse-communication loop
        !------------------------------------------------------------------
        ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it
        ! wants: 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its
        ! evaluation limit. A failed evaluation is reported with IV(2)=1
        ! (TOOBIG), which makes DRMNG shrink the step. The best point is
        ! tracked here because the last point DRMNG visits is not
        ! necessarily the lowest.
        !------------------------------------------------------------------
        DO WHILE (.NOT. (ExitNeeded))
          IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
          CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          SELECT CASE (IV(1))
          CASE (1)  ! Only energy is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO a = 1, nfo
              ActiveI = ActiveFunction(a)
              Glob_NonlinParam(1:npt, ActiveI) = x((a-1)*npt+1:a*npt)
            ENDDO
            Evalue = EnergyQA(.TRUE., ErrCode)
            NumOfEnergyEval = NumOfEnergyEval+1
            IF (ErrCode /= Q_METHOD_SUCCESS) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
            ELSE
              Glob_CurrEnergy = Evalue
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nv) = x(1:nv)
              ENDIF
            ENDIF
          CASE (2)  ! Only gradient is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO a = 1, nfo
              ActiveI = ActiveFunction(a)
              Glob_NonlinParam(1:npt, ActiveI) = x((a-1)*npt+1:a*npt)
            ENDDO
            CALL EnergyQB(Evalue, grad, .TRUE., ErrCode)
            NumOfGradEval = NumOfGradEval+1
            IF (ErrCode /= Q_METHOD_SUCCESS) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
            ELSE
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nv) = x(1:nv)
              ENDIF
            ENDIF
          CASE (3:8)  ! Some kind of convergence has been reached
            ExitNeeded = .TRUE.
          CASE (9:10)  ! Function evaluation limit has been reached.
            ! This is never supposed to happen because we
            ! count the number of function evaluations ourselves.
            ExitNeeded = .TRUE.
          CASE DEFAULT
            ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
            ! and >=14 for a bad input. None of those match a case above,
            ! so without this the loop would call DRMNG again for ever.
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                'Warning WC0139 in OptCycleQ: DRMNG returned IV(1) =', IV(1)
              IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
            ENDIF
            ExitNeeded = .TRUE.
          ENDSELECT
          ! A warning, not an abort: the best point found so far is
          ! still usable, and the acceptance tests below decide what
          ! becomes of this step.
          IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0123 in OptCycleQ: number of failures in energy or gradient', &
                'calculations during the optimization of nonlinear parameters', &
                'reached the limit of', Glob_MaxEnergyFailsAllowed
            ENDIF
          ENDIF
          IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.
        ENDDO  ! while
        !------------------------------------------------------------------
        ! Re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! DRMNG's last requested point need not be its lowest-energy point.
        ! The replacement transaction at x_best makes the parameters, the
        ! matrices, the factors, the energy and Glob_c describe one state.
        !------------------------------------------------------------------
        DO a = 1, nfo
          ActiveI = ActiveFunction(a)
          Glob_NonlinParam(1:npt, ActiveI) = x_best((a-1)*npt+1:a*npt)
        ENDDO
        Evalue = EnergyQAM(.TRUE., ErrCode)
        IF (ErrCode == Q_METHOD_SUCCESS) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          ! Glob_CurrEnergy is intentionally left unchanged here: the value
          ! returned by a failed evaluation is meaningless
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0120 in OptCycleQ: failed to evaluate energy after optimization', &
              'of nonlinear parameters. The values of the nonlinear parameters', &
              'will be left unchanged.'
          ENDIF
        ENDIF
        !------------------------------------------------------------------
        ! Acceptance test: pair overlaps
        !------------------------------------------------------------------
        ! Every unordered pair touching the window is checked once, against
        ! the whole basis - including functions with canonical indices
        ! greater than an active one, which a lower-triangle prefix loop
        ! would miss. Only the canonical lower triangle of Glob_S is read.
        ! Function NUMBERS are printed, as in OptCycleG.
        !------------------------------------------------------------------
        IsOverlapBad = .FALSE.
        IF ((ErrCode == Q_METHOD_SUCCESS) .AND. (OverlapThreshold > ZERO)) THEN
          ii = 0
          DO a = 1, nfo
            ActiveI = ActiveFunction(a)
            DO j = 1, cbs
              IF (j == ActiveI) CYCLE
              IF ((Q_Workspace%ActivePosition(j) > 0) .AND. (Q_Workspace%ActivePosition(j) < a)) CYCLE
              IF (ABS(QCanonicalMatrixElement(Glob_S, ActiveI, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0121: overlap of the following functions exceeds threshold. ', &
                      'Nonlinear parameters will be left unchanged'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') &
                    ii, ':', Glob_FuncNum(ActiveI), Glob_FuncNum(j), '    S='
                  CALL writerealadv(6, QCanonicalMatrixElement(Glob_S, ActiveI, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        !------------------------------------------------------------------
        ! Acceptance test: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..cbs, not just the window:
        ! moving one function can blow up the coefficient of another.
        !------------------------------------------------------------------
        IsAnyLinCoeffBad = .FALSE.
        IF ((ErrCode == Q_METHOD_SUCCESS) .AND. (LinCoeffThreshold > ZERO)) THEN
          ii = 0
          DO i = 1, cbs
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0122: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', Glob_FuncNum(i), '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        ! Reported only for a step that was actually kept
        IF ((Glob_ProcID == 0) .AND. (ErrCode == Q_METHOD_SUCCESS) .AND. (.NOT. IsOverlapBad) .AND. (.NOT. IsAnyLinCoeffBad)) THEN
          IF (Verbose >= 1) THEN
          WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
          ENDIF
          WRITE (*, *) 'E=', Glob_CurrEnergy
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 1) WRITE (*, *) 'Nonlinear parameters after optimization:'
            DO a = 1, nfo
              ActiveI = ActiveFunction(a)
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(ActiveI), ':', Glob_PWR(ActiveI)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, ActiveI), npt)
            ENDDO
          ENDIF
        ENDIF
        !------------------------------------------------------------------
        ! Undo the step if it was not acceptable
        !------------------------------------------------------------------
        ! The parameters go back to x_init and their columns are replaced
        ! again, so the matrices and the factors match the basis before the
        ! window moves on. ApplyQTrial restores a provable matrix/factor
        ! generation on an update failure; a solve failure leaves the
        ! restored physical point represented, so keeping E_prev is safe.
        !------------------------------------------------------------------
        IF ((ErrCode /= Q_METHOD_SUCCESS) .OR. IsOverlapBad .OR. IsAnyLinCoeffBad) THEN
          DO a = 1, nfo
            ActiveI = ActiveFunction(a)
            Glob_NonlinParam(1:npt, ActiveI) = x_init((a-1)*npt+1:a*npt)
          ENDDO
          Evalue = EnergyQA(.TRUE., ErrCode)
          IF (ErrCode == Q_METHOD_SUCCESS) THEN
            Glob_CurrEnergy = Evalue
          ELSE
            IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE(*, '(1x,a,1x,a)') &
              'Warning WC0124 in OptCycleQ: energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_CurrEnergy = E_prev
          ENDIF
        ENDIF
        !------------------------------------------------------------------
        ! Record the step and save
        !------------------------------------------------------------------
        ! On the last step of a cycle the position is reset to 0 and the
        ! cycle counter advances, which tells a resumed run that this cycle
        ! is finished. Saving happens on the first few steps whatever
        ! SaveEvery says, then every SaveEvery steps, and always on the
        ! final step of a cycle. The canonical order is the user-visible
        ! order: no sorting.
        !------------------------------------------------------------------
        IF (CurrFunc > FuncEnd-NumOfFuncToShift) THEN
          LastIter = .TRUE.
        ELSE
          LastIter = .FALSE.
        ENDIF
        Glob_History(cbs)%Energy = Glob_CurrEnergy
        IF (LastIter) THEN
          Glob_History(cbs)%InitFuncAtLastStep = 0
          Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
        ELSE
          Glob_History(cbs)%InitFuncAtLastStep = CurrFunc
        ENDIF
        IF (Glob_ProcID == 0) THEN
          IF ((totsteps <= Glob_MinMandSavSteps) .OR. (MOD(totsteps, SaveEvery) == 0) .OR. &
              (CurrFunc+NumOfFuncToShift >= FuncEnd)) THEN
            CALL SaveResults(Sort='no')
          ENDIF
        ENDIF
      ENDDO  ! end cycle CurrCycle
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Cycle', CurrCycle, ' finished'
      ENDIF
      IF (CurrCycle /= NumCycles) Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
    ENDDO  ! End of main optimization cycle
    !==================================================================
    ! Hand H and S to the next BBOP step and release everything
    !==================================================================
    CALL StoreMatricesInSwapFile()
    CALL ClearQWorkspace()
    DEALLOCATE(ActiveFunction)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x_init)
    DEALLOCATE(x)
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)
    DEALLOCATE(Glob_WkGR)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)
    ErrorCode = Q_METHOD_SUCCESS
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine OptCycleQ finished'
  END SUBROUTINE OptCycleQ


  SUBROUTINE FullOpt1Q(InitFunc, FinalFunc, MaxEnergyEval, OverlapThreshold, MaxOverlapPenalty, &
                       DataSaveMinTimeInterv, HessianSaveMinTimeInterv, HessFileName, ErrorCode)
    !==================================================================
    ! Subroutine FullOpt1Q
    !==================================================================
    ! Optimizes the nonlinear parameters of functions InitFunc..FinalFunc
    ! SIMULTANEOUSLY with the QR method ('Q') and DRMNG; the Q twin of
    ! FullOpt1G and FullOpt1I, with their DRMNG controls, smooth overlap
    ! PENALTY (OverlapThreshold >= 1.0 turns it off; the energy printed
    ! during the optimization includes it), Hessian restart and save policy,
    ! timed data saves and history accounting. The range is an explicit
    ! active map and is never moved to the end of the basis, so QR factors,
    ! H, S, coefficients, function identities and Hessian blocks keep one
    ! stable ordering. The Hessian is read from HessFileName on entry;
    ! after an improving evaluation the data file is written when
    ! DataSaveMinTimeInterv seconds have passed since its last save and the
    ! Hessian when HessianSaveMinTimeInterv seconds have; both are written
    ! once more at the end. A HessFileName of ' ', 'none', 'NONE' or 'None'
    ! disables the Hessian file. MaxEnergyEval is the total budget of
    ! energy evaluations of this basis size, the history count included;
    ! DRMNG enforces the remainder (IV(17)). The step ends at the best
    ! accepted point, never on a rejected trial point. ErrorCode:
    ! Q_METHOD_SUCCESS on return, Q_METHOD_INVALID_ARGUMENT when an
    ! argument is out of range (nothing is done then).
    !==================================================================
    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN)                        :: InitFunc, FinalFunc       ! range to optimize
    INTEGER, INTENT(IN)                        :: MaxEnergyEval             ! total budget, history included
    REAL(wp), INTENT(IN)                       :: OverlapThreshold          ! penalty threshold, >=1 = off
    REAL(wp), INTENT(IN)                       :: MaxOverlapPenalty         ! penalty magnitude
    REAL(4), INTENT(IN)                        :: DataSaveMinTimeInterv     ! seconds between data saves
    REAL(4), INTENT(IN)                        :: HessianSaveMinTimeInterv  ! seconds between Hessian saves
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: HessFileName
    INTEGER, INTENT(OUT)                       :: ErrorCode                 ! Q status of the step
    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER :: i                   ! loop counter
    INTEGER :: npt                 ! Glob_npt, parameters per function
    INTEGER :: nfa                 ! Glob_CurrBasisSize
    INTEGER :: nfo                 ! functions being optimized
    INTEGER :: nv                  ! nfo*npt, optimization variables
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ActiveFunction  ! the range as an active map
    REAL(wp) :: Evalue             ! value returned by the latest solve
    REAL(wp) :: CurrentEnergy      ! the value DRMNG is driven with
    REAL(wp) :: t                  ! scale factor for the DRMNG vector D
    REAL(wp) :: MaxAbsOverlap      ! reported by GetQOverlapStatistics
    REAL(wp) :: MinAbsOverlap      !   "
    REAL(wp) :: AverageAbsOverlap  !   "
    LOGICAL :: IsSwapFileOK       ! H and S came from the swap file
    LOGICAL :: ExitNeeded         ! ends the reverse-communication loop
    LOGICAL :: SaveHessian        ! HessFileName names a real file
    LOGICAL :: IsHessFileOK       ! a usable Hessian was read back
    LOGICAL :: IsHessSaveSuccess  ! set by SaveHessianFile, not read
    INTEGER :: ErrCode                            ! Q status of the latest transaction
    INTEGER :: NumOfFailures                      ! failed solves so far
    INTEGER :: NumOfEnergyEval                    ! energy evaluations so far
    INTEGER :: NumOfGradEval                      ! gradient evaluations so far
    INTEGER :: NumOfEnergyEvalDuringFullOpt_Init  ! count carried in from history
    REAL(4) :: TimeOfLastSave      ! CPU time of the last data save
    REAL(4) :: TimeOfLastHessSave  ! CPU time of the last Hessian save
    ! The optimization variables
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, grad
    ! Arrays and settings used by DRMNG (rank 0 only)
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D         ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V         ! work array
    INTEGER, PARAMETER                  :: LIV = 60  ! length of IV
    INTEGER                             :: IV(LIV)
    INTEGER                             :: LV        ! length of V
    INTEGER                             :: ALG       ! 2 = unconstrained minimization
    INTEGER                             :: IVLMAT    ! IV(42), where V holds the Hessian
    ErrorCode = Q_METHOD_INVALID_ARGUMENT
    nfa = Glob_CurrBasisSize
    IF ((InitFunc < 1) .OR. (FinalFunc < InitFunc) .OR. (FinalFunc > nfa)) RETURN
    IF (MaxEnergyEval < 0) RETURN
    !==================================================================
    ! Overlap penalty
    !==================================================================
    IsHessFileOK = .FALSE.
    IF (OverlapThreshold >= ONE) THEN
      Glob_OverlapPenaltyAllowed = .FALSE.
    ELSE
      Glob_OverlapPenaltyAllowed = .TRUE.
      Glob_OverlapPenaltyThreshold2 = OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty = MaxOverlapPenalty
    ENDIF
    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine FullOpt1Q started'
      IF (Verbose >= 1) WRITE(*, *) 'Simultaneous optimization of nonlinear parameters of basis functions'
      IF (Verbose >= 1) WRITE(*, *) InitFunc, '  through', FinalFunc, '  will be attempted'
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Overlap threshold is ', ABS(OverlapThreshold)
        IF (Verbose >= 1) WRITE(*, *) 'Max value of a pair overlap penalty is ', Glob_MaxOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Warning! The energy value that will be shown during the optimization'
        WRITE(*, *) 'may differ from the actual energy'
      ELSE
        IF (Verbose >= 1) WRITE(*, *) 'No constraints on overlaps will be imposed'
      ENDIF
    ENDIF
    !==================================================================
    ! Global state and array allocation
    !==================================================================
    Glob_GSEPSolutionMethod = 'Q'
    Glob_nfa = nfa
    Glob_nfo = FinalFunc-InitFunc+1
    Glob_nfru = nfa-Glob_nfo
    Glob_HSLeadDim = nfa
    npt = Glob_npt
    nfo = Glob_nfo
    nv = nfo*npt
    Glob_HSBuffLen = MAX(MIN(nfa*(nfa+1)/2, 1000), 30*nfa)
    ALLOCATE(Glob_H(nfa, nfa))
    ALLOCATE(Glob_S(nfa, nfa))
    ALLOCATE(Glob_diagS(nfa))
    ALLOCATE(Glob_D(2*npt, nfo, nfa))
    ALLOCATE(Glob_c(nfa))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ! Allocate workspace for EnergyQB
    ALLOCATE(Glob_WkGR(nv))
    ! D and V are only ever touched on rank 0, which is the only rank
    ! that runs DRMNG, so they are allocated there alone.
    LV = 71 + nv*(nv+13)/2 + 1
    IF (Glob_ProcID == 0) THEN
      ALLOCATE(D(nv))
      ALLOCATE(V(LV))
    ENDIF
    ALLOCATE(x(nv))
    ALLOCATE(grad(nv))
    ALLOCATE(ActiveFunction(nfo))
    CALL PrepareQWorkspace(nfa, nfa, nfo, ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0160 in FullOpt1Q: Q workspace cannot be allocated'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    DO i = 1, nfo
      ActiveFunction(i) = InitFunc+i-1
    ENDDO
    ! Setting up a logical variable that determines
    ! whether the hessian should be saved from time to time
    IF ((HessFileName == ' ') .OR. (HessFileName == 'none') .OR. &
        (HessFileName == 'NONE') .OR. (HessFileName == 'None')) THEN
      SaveHessian = .FALSE.
    ELSE
      SaveHessian = .TRUE.
    ENDIF
    !==================================================================
    ! Initial energy
    !==================================================================
    ! The full canonical problem is restored from the swap file or
    ! assembled once; then the active map is installed and the initial
    ! factors are built. No permutation and no sorting workspace.
    !------------------------------------------------------------------
    Glob_H = ZERO
    Glob_S = ZERO
    Glob_diagS = ZERO
    Glob_c = ONE
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)
    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a59)', ADVANCE='no') &
        'Computing matrix elements and solving eigenvalue problem...'
      CALL ComputeMatElem(1, nfa)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    ENDIF
    Q_Workspace%MatricesAreCanonical = .TRUE.
    CALL SetQActiveFunctions(ActiveFunction, ErrCode)
    IF (ErrCode == Q_METHOD_SUCCESS) CALL CaptureQMatrixParameters(ErrCode)
    IF (ErrCode == Q_METHOD_SUCCESS) CALL FactorizeQFresh(ErrCode)
    IF (ErrCode == Q_METHOD_SUCCESS) Glob_CurrEnergy = EnergyQA(.FALSE., ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) WRITE(*, '(1x,a,1x,i0)') &
        'Error EC0161 in FullOpt1Q: initial Q energy cannot be computed, status', ErrCode
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'
    CALL GetQOverlapStatistics(MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap, ErrCode)
    IF (Glob_ProcID == 0) THEN
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                           ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (including overlap penalty ', Glob_CurrEnergy
      ELSE
        WRITE(*, *) 'Initial energy                            ', Glob_CurrEnergy
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                           ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                           ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap              ', AverageAbsOverlap
    ENDIF
    CALL CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave = Glob_TimeSinceStart
    TimeOfLastHessSave = Glob_TimeSinceStart
    !==================================================================
    ! Set up DRMNG and load the Hessian
    !==================================================================
    ! The packed Hessian ordering is [parameters of InitFunc, parameters
    ! of InitFunc+1, ...], which is exactly the active-map order and
    ! therefore stable across restarts. IV(25)=0 tells DRMNG not to
    ! overwrite the Hessian it was handed. The scale vector D is only built
    ! when there is no usable stored one to reuse.
    !------------------------------------------------------------------
    ALG = 2
    IF (Glob_ProcID == 0) THEN
      CALL DIVSET(ALG, IV, LIV, LV, V)
      IV(18) = 1000000  ! iteration limit; the evaluation limit IV(17) is set from the budget below
      IV(19) = -1  ! set summary print format
      IV(20) = 0; IV(22) = 0; IV(23) = -1; IV(24) = 0
      V(31) = 0.0_wp
      V(32) = 2*EPSILON(V(32))
      V(37) = 2*EPSILON(V(37))
      ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
      ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
      ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
      V(35) = Glob_MaxScStepAllowedInOpt
      IV(1) = 12  ! DIVSET has been called and some default values were changed
    ENDIF
    DO i = 1, nfo
      x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, ActiveFunction(i))
    ENDDO
    IF (Glob_ProcID == 0) THEN
      ! If SaveHessian=.true. then try to read the Hessian from the file
      IF (SaveHessian) THEN
        IVLMAT = IV(42)
        CALL ReadHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessFileOK)
        IF (IsHessFileOK) IV(25) = 0
      ENDIF
      IF ((.NOT. Glob_FullOptSaveD) .OR. (.NOT. IsHessFileOK) .OR. (.NOT. SaveHessian)) THEN
        ! Set the scaling vector
        t = MAX(ONE/(nfa*nfa*SQRT(ONE*nfa)), 10000*EPSILON(Glob_CurrEnergy))
        D(1:nv) = t
      ENDIF
    ENDIF
    !==================================================================
    ! The reverse-communication loop
    !==================================================================
    ! IV(1) on return: 1 = wants an energy, 2 = wants a gradient,
    ! 3..8 = converged, 9,10 = its own limit reached. A failed energy
    ! is reported with IV(2)=1, which makes DRMNG shrink the step.
    !------------------------------------------------------------------
    ExitNeeded = .FALSE.
    NumOfFailures = 0
    NumOfEnergyEval = 0
    NumOfGradEval = 0
    NumOfEnergyEvalDuringFullOpt_Init = Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt
    ! The remaining budget of energy evaluations of this basis size is
    ! DRMNG's own function-evaluation limit (checked after every
    ! evaluation, so the count is exact)
    IF (Glob_ProcID == 0) IV(17) = MAX(1, MaxEnergyEval-NumOfEnergyEvalDuringFullOpt_Init)
    ! DRMNG does not read FX on its first (IV(1)=12) entry, but giving
    ! it a defined value keeps -finit-real=nan builds quiet.
    CurrentEnergy = Glob_CurrEnergy
    DO WHILE (.NOT. (ExitNeeded))
      IF (Glob_ProcID == 0) CALL DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
      CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      SELECT CASE (IV(1))
      CASE (1)  ! Only energy is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, ActiveFunction(i)) = x((i-1)*npt+1:i*npt)
        ENDDO
        Evalue = EnergyQA(.TRUE., ErrCode)
        NumOfEnergyEval = NumOfEnergyEval+1
        IF (ErrCode /= Q_METHOD_SUCCESS) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
        ELSE
          CurrentEnergy = Evalue
        ENDIF
        !--------------------------------------------------------------
        ! Periodic saving, on rank 0, only when the energy improved
        !--------------------------------------------------------------
        IF ((Glob_ProcID == 0) .AND. (ErrCode == Q_METHOD_SUCCESS)) THEN
          IF (Evalue < Glob_CurrEnergy) THEN
            CALL CPU_TIME(Glob_TimeSinceStart)
            IF (Glob_TimeSinceStart-TimeOfLastSave > DataSaveMinTimeInterv) THEN
              IF (Glob_OverlapPenaltyAllowed) THEN
                Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
              ELSE
                Glob_CurrEnergy = Evalue
              ENDIF
              Glob_History(nfa)%Energy = Glob_CurrEnergy
              Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt = &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              CALL SaveResults(Sort='no')
              WRITE(*, *) 'Data file has been updated'
              CALL GetQOverlapStatistics(MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap, ErrCode)
              IF (Verbose >= 1) WRITE(*, *) 'Some current statistics:'
              IF (Glob_OverlapPenaltyAllowed) THEN
                IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
              ELSE
                WRITE(*, *) 'Energy                            ', Evalue
              ENDIF
              IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
              TimeOfLastSave = Glob_TimeSinceStart
            ENDIF
            IF ((Glob_TimeSinceStart-TimeOfLastHessSave > HessianSaveMinTimeInterv) &
                .AND. (SaveHessian)) THEN
              CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
              TimeOfLastHessSave = Glob_TimeSinceStart
            ENDIF
          ENDIF
        ENDIF
      CASE (2)  ! Only gradient is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, ActiveFunction(i)) = x((i-1)*npt+1:i*npt)
        ENDDO
        CALL EnergyQB(Evalue, grad, .TRUE., ErrCode)
        NumOfGradEval = NumOfGradEval+1
        IF (ErrCode /= Q_METHOD_SUCCESS) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
        ENDIF
      CASE (3:8)  ! Some kind of convergence has been reached
        ExitNeeded = .TRUE.
      CASE (9:10)  ! DRMNG reached its evaluation (or iteration) limit
        ! IV(17) holds the remaining energy-evaluation budget of this basis
        ! size, so this is the normal end of a step that exhausts it. x may
        ! still sit on a rejected trial point; the code after the loop goes
        ! back to the accepted point when that is the case.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *) 'Warning WC0130 in FullOpt1Q: number of energy evaluations reached limit'
          IF (Verbose >= 1) WRITE(*, '(1x,a,i0,a,i0,a)') '(', NumOfEnergyEval, ' in this step, ', &
            NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval, ' in total)'
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.
      CASE DEFAULT
        ! DRMNG returns 63 or 65 when it gives up on an uncomputable
        ! value, and >=14 for a bad input. None of those match a case
        ! above, so without this the loop would call DRMNG again for
        ! ever. Leave with the best point found so far.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
            'Warning WC0136 in FullOpt1Q: DRMNG returned IV(1) =', IV(1)
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.
      ENDSELECT
      IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,a,1x,a)') &
            'Error EC0162 in FullOpt1Q: number of failures in energy or gradient', &
            'calculations during the optimization of nonlinear parameters', &
            'exceeded limit'
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF
    ENDDO
    !==================================================================
    ! Final energy at the best point found
    !==================================================================
    ! DRMNG can leave x on a rejected trial point of the last line search
    ! (evaluation limit, false convergence). Its accepted iterate x0 is in
    ! V(IV(43):IV(43)+nv-1) with f0 = V(13); V(10) is the last f evaluated.
    ! Whenever the point in x is not better than x0, go back to x0, so
    ! that the final energy and the data file describe the best point.
    ! The replacement transaction at that point makes the parameters, the
    ! matrices, the factors, the energy and Glob_c one final generation.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      IF ((IV(31) >= 1) .AND. (V(10) >= V(13))) THEN
        IF ((V(10) > V(13)) .AND. (Verbose >= 2)) &
          WRITE(*, *) 'The last trial point is discarded: back to the last accepted point'
        x(1:nv) = V(IV(43):IV(43)+nv-1)
      ENDIF
    ENDIF
    CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    DO i = 1, nfo
      Glob_NonlinParam(1:npt, ActiveFunction(i)) = x((i-1)*npt+1:i*npt)
    ENDDO
    Evalue = EnergyQA(.TRUE., ErrCode)
    IF (ErrCode /= Q_METHOD_SUCCESS) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, '(1x,a,1x,i0)') 'Error EC0163 in FullOpt1Q: failed to evaluate energy after the optimization, status', ErrCode
        IF (Verbose >= 1) WRITE(*, *) 'of nonlinear parameters'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    ! Glob_CurrEnergy is the physical energy, without the penalty
    IF (Glob_OverlapPenaltyAllowed) THEN
      Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
    ELSE
      Glob_CurrEnergy = Evalue
    ENDIF
    CALL GetQOverlapStatistics(MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap, ErrCode)
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
      WRITE(*, *) 'Final energy and overlap statistics:'
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
      ELSE
        WRITE(*, *) 'Energy                            ', Evalue
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
    ENDIF
    ! Adding data to history
    Glob_History(nfa)%Energy = Glob_CurrEnergy
    Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt = &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
    !==================================================================
    ! Save and release everything
    !==================================================================
    ! The Hessian goes to disk with the data file, so a restart finds
    ! the file even when no evaluation improved the energy
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      CALL SaveResults(Sort='no')
      IF (SaveHessian) CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
    ENDIF
    CALL StoreMatricesInSwapFile()
    IF (Glob_OverlapPenaltyAllowed) Glob_OverlapPenaltyAllowed = .FALSE.
    CALL ClearQWorkspace()
    DEALLOCATE(ActiveFunction)
    DEALLOCATE(grad)
    DEALLOCATE(x)
    IF (Glob_ProcID == 0) THEN
      DEALLOCATE(V)
      DEALLOCATE(D)
    ENDIF
    DEALLOCATE(Glob_WkGR)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)
    ErrorCode = Q_METHOD_SUCCESS
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine FullOpt1Q finished'
  END SUBROUTINE FullOpt1Q


  SUBROUTINE BasisEnlG(Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval, &
                       OverlapThreshold, LinCoeffThreshold)
    !==================================================================
    ! Subroutine BasisEnlG
    !==================================================================
    ! Enlarges the basis from Kstart-1 to Kstop functions, Kstep at a time,
    ! with the GSEP solved by DSYGVX ('G'); BasisEnlI is the inverse-
    ! iteration twin. For each block: (1) stochastic selection - NTrials
    ! candidate blocks drawn from the distribution of the existing
    ! parameters, the best kept; (2) premultiplier powers, scanned one
    ! function at a time in random order (skipped with Glob_IsIndexFixed);
    ! (3) nonlinear parameters optimized with DRMNG through the
    ! reverse-communication loop (OptimizationType 1); (4) acceptance - the
    ! block is regenerated when the energy cannot be evaluated, a new pair
    ! overlap exceeds OverlapThreshold or a linear coefficient exceeds
    ! LinCoeffThreshold, up to Glob_BadOverlapOrLinCoeffLim attempts (then
    ! the last block is kept). A threshold <= 0 disables its test. The
    ! names ZIndSet/ZIndSetBest/ZIndOptSequence and the "Z-index" messages
    ! refer to the premultiplier power in Glob_PWR.
    ! Arguments: Kstart, Kstop, Kstep (the last block may be smaller),
    ! NTrials, OptimizationType (0 keep the selected parameters, 1 optimize;
    ! else EC0129), MaxEnergyEval (per block), OverlapThreshold,
    ! LinCoeffThreshold. Precondition: the first Kstart-1 functions are in
    ! Glob_NonlinParam and Glob_PWR.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval
    REAL(wp), INTENT(IN) :: OverlapThreshold, LinCoeffThreshold

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    ! Local variables:
    INTEGER  :: i, j, K, AttemptToGetGoodFunc, ii, jj, jbest
    INTEGER  :: np, npt, nfo, nfa, nfru, nfrup1, nvmax, nv
    INTEGER  :: OpenFileErr, ErrCode, NumOfFailures, NumOfEnergyEval, NumOfGradEval
    LOGICAL  :: IsSwapFileOK, IsEnergyImproved, ExitNeeded
    LOGICAL  :: IsOverlapBad, IsAnyLinCoeffBad, IsEnergyBad
    LOGICAL  :: IsShapeBad
    INTEGER  :: NumOfShapeRedraws
    REAL(wp) :: Cfac, ShapeSsum, ShapeSabs
    INTEGER  :: wbfu_t, wmu_t, wbfu, wmu, rgm1_counter, rgm2_counter, BlockSizeForDSYGVX
    REAL(wp) :: ms1, ms2
    REAL(wp) :: Evalue, E_init, E_best
    REAL(wp) :: t

    ! Largest legal premultiplier power: the greatest EVEN value not
    ! exceeding Glob_MaxPowerAllowed. Powers are always even in this
    ! frame - GenerateTrialParam only ever produces even ones - so the
    ! power scan below steps through 2,4,...,PWRMax.
    INTEGER, PARAMETER :: PWRMax = 2*(Glob_MaxPowerAllowed/2)

    ! Candidate block and the best candidate block found so far
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: ParSet, ParSetBest
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: ZIndSet, ZIndSetBest

    ! The optimization variables: the nonlinear parameters of the block
    ! laid out as one flat vector of nfo*npt elements
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_best, grad
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: ZIndOptSequence

    ! Arrays used by DRMNG
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D, V, V_init
    INTEGER, PARAMETER                  :: LIV = 60
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV
    INTEGER                             :: ALG

    ! Allocatable work space
    ! *** These six are NOT REFERENCED anywhere in this routine. They are
    ! left over from when the reallocation of the per-function arrays was
    ! written out inline here; it now lives in ReallocateBasisFuncData,
    ! which declares its own copies. Same for OpenFileErr above, which was
    ! used by the inline swap-file handling now in
    ! ReadSwapFileAndDistributeData.
    REAL(wp), ALLOCATABLE, DIMENSION(:)               :: WorkBuffReal
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: WorkBuffInt
    TYPE(Glob_HistoryStep), ALLOCATABLE, DIMENSION(:) :: TempHistory
    REAL(wp), ALLOCATABLE, DIMENSION(:, :)            :: TempParam
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempZInd
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempFunc
    !====================================================
    ! These variables are used when a finite difference gradient is computed
    ! real(wp),allocatable,dimension(:)     ::    fx,fgrad
    ! real(wp)                                    deltax,Evalue1
    !====================================================


    !==================================================================
    ! Announce the step
    !==================================================================
    ! wbfu_t/wmu_t are filled by GenerateTrialParam on rank 0 and carry
    ! which existing function was used as prototype and which generator
    ! method produced the candidate. They feed the statistics printed at
    ! the very end, so they are cleared once here.
    !------------------------------------------------------------------
    wbfu_t = 0
    wmu_t = 0

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlG started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstart =', Kstart
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstop =', Kstop
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstep =', Kstep
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'OptimizationType =', OptimizationType
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval =', MaxEnergyEval
    ENDIF


    !==================================================================
    ! Global state this routine works under
    !==================================================================
    ! Glob_nfru/Glob_nfo/Glob_nfa define the optimization window (nfru
    ! frozen functions, nfo being optimized, nfa the last function) that
    ! the energy routines and Glob_D are indexed against; they are reset
    ! per iteration in the main loop. Overlap penalties are off: an overlap
    ! violation REJECTS the block instead of penalizing the energy.
    !------------------------------------------------------------------
    Glob_GSEPSolutionMethod = 'G'
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_nfa = Kstart+Kstep
    Glob_nfo = Kstep
    Glob_HSLeadDim = Kstop
    Glob_HSBuffLen = Kstop*Kstep
    np = Glob_np
    npt = Glob_npt
    nfo = Glob_nfo
    nfa = Glob_nfa
    nvmax = Kstep*Glob_npt

    ! Reallocate arrays that contain the information
    ! about basis functions and optimization process.
    CALL ReallocateBasisFuncData(Kstop, Glob_CurrBasisSize)


    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! Only the LOWER triangles of Glob_H and Glob_S are meaningful; the
    ! upper ones are used as scratch by the permutation routines.
    !
    ! Glob_D holds the derivatives of the matrix elements with respect
    ! to the nonlinear parameters of the nfo functions in the window -
    ! 2*npt of them per (window function, basis function) pair.
    !------------------------------------------------------------------
    ! Allocate some global arrays
    ALLOCATE(Glob_H(Kstop, Kstop))
    ALLOCATE(Glob_S(Kstop, Kstop))
    ALLOCATE(Glob_diagH(Kstop))
    ALLOCATE(Glob_diagS(Kstop))
    ALLOCATE(Glob_c(Kstop))
    ALLOCATE(Glob_D(2*npt, Kstep, Kstop))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    !==================================================================
    ! LAPACK workspace for DSYGVX
    !==================================================================
    ! (BlockSize+3)*Kstop is the OPTIMAL size ILAENV suggests; 8*Kstop is
    ! the minimum DSYGVX requires. Taking the larger of the two is what
    ! keeps the call valid when ILAENV returns a small block size, which
    ! it does for different values in Netlib and MKL.
    !------------------------------------------------------------------
    ! Allocate workspace for DSYGVX
    BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', Kstop, Kstop, Kstop, Kstop)
    Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*Kstop, 8*Kstop)
    ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    ALLOCATE(Glob_IWorkForDSYGVX(5*Kstop))

    ! Allocate workspace for EnergyGB
    ALLOCATE(Glob_WkGR(Kstep*npt))

    !------------------------------------------------------------------
    ! Local workspace
    !------------------------------------------------------------------
    ! Allocate workspace
    ALLOCATE(ParSet(npt, Kstep))
    ALLOCATE(ParSetBest(npt, Kstep))
    ALLOCATE(ZIndSet(Kstep))
    ALLOCATE(ZIndSetBest(Kstep))
    ALLOCATE(x(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(ZIndOptSequence(Kstep))


    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next (1 energy,
    ! 2 gradient) and is called again, which suits an energy that is a
    ! collective operation over all processes.
    ! LV is the documented size of the V work array plus one.
    !------------------------------------------------------------------
    nvmax = npt*Kstep
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))

    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next (1 energy,
    ! 2 gradient) and is called again, which suits an energy that is a
    ! collective operation over all processes.
    ! The parameters are set once, outside the main loop.
    !------------------------------------------------------------------

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): iteration and function-evaluation limits, set out of
    ! the way because the budget is enforced by MaxEnergyEval below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    ! V(35)=0.1*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed


    !==================================================================
    ! Initial energy
    !==================================================================
    ! The swap file, when it is valid, carries H and S for the basis we
    ! start from, so the matrix elements do not have to be recomputed -
    ! hence .false. for AreMatElemNeeded in that branch.
    !
    ! With Kstart==1 there is no basis yet and no energy to compute, so
    ! Glob_CurrEnergy is primed with HUGE() and the first candidate to
    ! produce a finite energy wins.
    !------------------------------------------------------------------
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    ! Calculating the initial energy
    ErrCode = 0
    IF (Kstart > 1) THEN
      IF (IsSwapFileOK) THEN
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Solving eigenvalue problem...'
        Glob_CurrEnergy = EnergyGA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
      ELSE
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem...'
        Glob_CurrEnergy = EnergyGA(1, Glob_CurrBasisSize, .TRUE., ErrCode)
      ENDIF
    ELSE
      Glob_CurrEnergy = HUGE(Glob_CurrEnergy)
    ENDIF
    IF (ErrCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0125 in BasisEnlG: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy

    rgm1_counter = 0
    rgm2_counter = 0
    ms1 = ZERO
    ms2 = ZERO
    K = Kstart-1


    !==================================================================
    ! MAIN LOOP - one block of up to Kstep functions per iteration
    !==================================================================
    ! Main loop begins here
    DO WHILE (K < Kstop)

      !--------------------------------------------------------------
      ! Size and place the window for this block
      !--------------------------------------------------------------
      ! nfru = functions that stay frozen, nfo = functions being added,
      ! K = the new basis size. The last block is short when Kstop-K is
      ! less than Kstep.
      IF (K+Kstep <= Kstop) THEN
        nfo = Kstep
        nfru = K
        K = K+Kstep
      ELSE
        nfo = Kstop-K
        nfru = K
        K = Kstop
      ENDIF

      CALL linalg_setparam(K)  ! reset linalg flags to account for changes in the basis size

      Glob_nfa = K
      Glob_nfru = nfru
      Glob_nfo = nfo
      nfrup1 = nfru+1
      nv = nfo*npt
      E_init = Glob_CurrEnergy

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Current basis size is', Glob_CurrBasisSize
        IF (nfo > 1) THEN
          WRITE(*, '(1x,a,1x,i0,a,i0)') 'Selecting functions', nfrup1, '-', K
        ELSE
          WRITE(*, '(1x,a,1x,i0)') 'Selecting function', K
        ENDIF
      ENDIF


      !==================================================================
      ! ACCEPTANCE LOOP - regenerate the block until it is acceptable
      !==================================================================
      ! The three flags are primed so that the loop always runs at least
      ! once. It ends when the block passes the energy, overlap and linear
      ! coefficient tests, or when the attempt budget runs out - in which
      ! case the last block is kept regardless, which is deliberate: a
      ! basis that is slightly too linearly dependent is better than no
      ! progress at all. A block in which the Young operator nearly
      ! annihilates a new function (IsShapeBad, C > Glob_MaxSelfOverlapCancel)
      ! is redrawn without limit and without spending the attempt budget.
      !------------------------------------------------------------------
      IsOverlapBad = .TRUE.
      IsAnyLinCoeffBad = .TRUE.
      IsEnergyBad = .FALSE.
      IsShapeBad = .FALSE.
      NumOfShapeRedraws = 0
      AttemptToGetGoodFunc = 1

      DO WHILE (((IsOverlapBad .OR. IsAnyLinCoeffBad .OR. IsEnergyBad) .AND. &
                 (AttemptToGetGoodFunc <= Glob_BadOverlapOrLinCoeffLim)) .OR. IsShapeBad)

        !------------------------------------------------------------------
        ! Step 1: stochastic selection of the block
        !------------------------------------------------------------------
        ! GenerateTrialParam runs on rank 0 (it consumes the random stream) and
        ! the result is broadcast; every rank evaluates the SAME candidate. A
        ! candidate that makes DSYGVX fail is counted, not fatal; only an
        ! excessive FRACTION of failures is.
        !------------------------------------------------------------------
        NumOfFailures = 0
        IsEnergyImproved = .FALSE.
        wbfu = 0
        wmu = 0

        DO i = 1, NTrials
          IF (Glob_ProcID == 0) CALL GenerateTrialParam(nfo, ParSet, ZIndSet, wbfu_t, wmu_t)
          CALL MPI_BCAST(ParSet, npt*nfo, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_BCAST(ZIndSet, nfo, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          Glob_NonlinParam(1:npt, nfrup1:K) = ParSet(1:npt, 1:nfo)
          Glob_PWR(nfrup1:K) = ZIndSet(1:nfo)
          Evalue = EnergyGA(nfrup1, K, .TRUE., ErrCode)
          IF (ErrCode == 0) THEN
            IF (Evalue < Glob_CurrEnergy) THEN
              Glob_CurrEnergy = Evalue
              ParSetBest(1:npt, 1:nfo) = ParSet(1:npt, 1:nfo)
              ZIndSetBest(1:nfo) = ZIndSet(1:nfo)
              IsEnergyImproved = .TRUE.
              wbfu = wbfu_t
              wmu = wmu_t
            ENDIF
          ELSE
            NumOfFailures = NumOfFailures+1
          ENDIF
        ENDDO

        ! Too many candidates the solver could not handle at all: the
        ! basis is in a state the eigensolver cannot work with, and more
        ! trials will not fix it.
        IF (NumOfFailures*ONE/NTrials > Glob_MaxFracOfTrialFailsAllowed) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0126 in BasisEnlG: the number of eigenvalue problem solution failures'
            WRITE(*, *) 'in random selection process exceeded limit'
            WRITE(*, '(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials, '%'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! The candidates were all solvable but none lowered the energy.
        IF (.NOT. (IsEnergyImproved)) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0127 in BasisEnlG: random selection did not result'
            WRITE(*, *) 'in any energy improvement'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! Put the winning candidate back in place. The re-solve is needed
        ! because the last candidate evaluated is not in general the best
        ! one, so Glob_H/Glob_S hold the wrong block.
        Glob_NonlinParam(1:npt, nfrup1:K) = ParSetBest(1:npt, 1:nfo)
        Glob_PWR(nfrup1:K) = ZIndSetBest(1:nfo)
        Glob_CurrEnergy = EnergyGA(nfrup1, K, .TRUE., ErrCode)

        IF (Glob_ProcID == 0) THEN
          WRITE (*, '(1x,a)', ADVANCE='no') 'E='
          CALL writereal(6, Glob_CurrEnergy)
          IF (Verbose >= 2) WRITE (*, '(5x,a,1x,i0)') 'prototype function is', wbfu
          DO i = 1, nfo
            WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', ZIndSetBest(i)
            CALL writerealarradv(6, ParSetBest(1:npt, i), npt)
          ENDDO
          IF (Verbose >= 1) WRITE (*, *) 'Optimizing nonlinear parameters'
        ENDIF


        !------------------------------------------------------------------
        ! Prime the best point found
        !------------------------------------------------------------------
        ! Every branch of the SELECT below must leave x_best defined, because
        ! the energy and the linear coefficients are recomputed at x_best right
        ! after it; the randomly selected point is also what OptimizationType 0
        ! keeps.
        !------------------------------------------------------------------
        DO i = 1, nfo
          x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
        ENDDO
        E_best = Glob_CurrEnergy
        x_best(1:nfo*npt) = x(1:nfo*npt)


        ! Optimization of the nonlinear parameters
        SELECT CASE (OptimizationType)

        !------------------------------------------------------------------
        ! OptimizationType 0: keep the randomly selected block
        !------------------------------------------------------------------
        CASE (0)  ! No optimization
        ! x and x_best already hold the selected point - nothing to do.
        ! Note that this is NOT what OptimizationType==0 meant in
        ! workproc.f90, where anything other than 1 selected the DMNG
        ! optimization. A driver that wants the parameters optimized
        ! must now pass 1.

        !------------------------------------------------------------------
        ! OptimizationType 1: optimize powers, then nonlinear parameters
        !------------------------------------------------------------------
        CASE (1)

          !------------------------------------------------------------------
          ! Step 2: premultiplier powers, one function at a time
          !------------------------------------------------------------------
          ! Each function in the window is tried with every legal EVEN power
          ! 2,4,...,PWRMax and the best is kept; the functions are visited in
          ! random order. This exhaustive scan costs PWRMax/2 - 1 energy
          ! evaluations per function (a 1-D minimizer or a window around the
          ! current power would be cheaper). Skipped when the power is pinned
          ! by Glob_IsIndexFixed.
          !------------------------------------------------------------------
          NumOfFailures = 0
          IF (.NOT. Glob_IsIndexFixed) THEN
            ! Generate a random sequence which will define the order in which
            ! Z-indices should be optimized (one index at a time)
            CALL GenerateRndIntSeq(nfo, ZIndOptSequence)
            ! Loop where Z-indices are optimized. Note: there is some room for improvement
            ! here as I programmed it in a simple way when all matrix element of functions
            ! nfrup1 through K are computed each time while it is not always necessary.
            DO i = 1, nfo
              ii = ZIndOptSequence(i)
              j = Glob_PWR(nfru+ii)
              jbest = j
              DO jj = 2, PWRMax, 2
                IF (jj /= j) THEN
                  Glob_PWR(nfru+ii) = jj
                  Evalue = EnergyGA(nfrup1, K, .TRUE., ErrCode)
                  IF (ErrCode /= 0) THEN
                    ! Restore the best power known so far before carrying on,
                    ! so a failed trial never leaves a bad power behind.
                    NumOfFailures = NumOfFailures+1
                    Glob_PWR(nfru+ii) = jbest
                    IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
                      IF (Glob_ProcID == 0) THEN
                        WRITE(*, *) 'Error EC0128 in BasisEnlG: number of failures in energy calculations'
                        WRITE(*, *) 'during the optimization of Z-indicies exceeded limit'
                      ENDIF
                      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
                    ENDIF
                  ELSE
                    IF (Evalue < Glob_CurrEnergy) THEN
                      Glob_CurrEnergy = Evalue
                      jbest = jj
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              Glob_PWR(nfru+ii) = jbest
            ENDDO
          ENDIF

          !------------------------------------------------------------------
          ! Step 3: nonlinear parameters, with DRMNG
          !------------------------------------------------------------------
          ! IV and V are restored from the copies made before the main loop, so
          ! each block starts the minimizer from its documented default state.
          !------------------------------------------------------------------
          ! Now we optimize nonlinear parameters

          ! Setting IV and V values as was in their initial copies
          IV(1:LIV) = IV_init(1:LIV)
          V(1:LV) = V_init(1:LV)

          nv = nfo*npt
          DO i = 1, nfo
            x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
          ENDDO

          !--------------------------------------------------------------
          ! Scale vector D for DRMNG
          !--------------------------------------------------------------
          ! Every variable gets the same scale t, taken from the RELATIVE energy
          ! gain of the block so far (a block that already moved the energy a lot
          ! may take larger steps), floored at 10000*epsilon. At the very start
          ! of a basis (nfru>=nfo fails, E_init is the HUGE() sentinel) t=1.
          IF (nfru >= nfo) THEN
            t = MAX(ABS((E_init-Glob_CurrEnergy))/(ABS(E_init)+ABS(Glob_CurrEnergy)), &
                    10000*EPSILON(Glob_CurrEnergy))
          ELSE
            t = ONE
          ENDIF
          ! if (Glob_ProcID==0) write(*,*) 'scaling coeff=',t !remove later
          DO i = 1, nfo
            ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
            DO j = 1, npt
              ! Make sure none of the D(i) will be zero or smaller than the threshold
              ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
              D(npt*(i-1)+j) = t
              ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
            ENDDO
          ENDDO

          ExitNeeded = .FALSE.
          NumOfFailures = 0
          NumOfEnergyEval = 0
          NumOfGradEval = 0
          IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.
          E_best = Glob_CurrEnergy
          x_best(1:nfo*npt) = x(1:nfo*npt)

          !------------------------------------------------------------------
          ! The reverse-communication loop
          !------------------------------------------------------------------
          ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it wants:
          ! 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its evaluation
          ! limit. On a failed evaluation the energy handed to DRMNG keeps its
          ! previous value, so the step counts as no reduction and the trust
          ! radius shrinks (IV(2), the TOOBIG flag, is deliberately left 0). The
          ! best point is tracked here because the last point DRMNG visits is not
          ! necessarily the lowest.
          !------------------------------------------------------------------
          DO WHILE (.NOT. (ExitNeeded))

            IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
            CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

            SELECT CASE (IV(1))

            CASE (1)  ! Only energy is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              Evalue = EnergyGA(nfrup1, K, .TRUE., ErrCode)
              NumOfEnergyEval = NumOfEnergyEval+1
              IF (ErrCode /= 0) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
              ELSE
                Glob_CurrEnergy = Evalue
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nfo*npt) = x(1:nfo*npt)
                ENDIF
              ENDIF
              ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
              IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

            CASE (2)  ! Only gradient is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              CALL EnergyGB(Evalue, grad, .TRUE., ErrCode)
              NumOfGradEval = NumOfGradEval+1
              IF (ErrCode /= 0) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
              ELSE
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nfo*npt) = x(1:nfo*npt)
                ENDIF
              ENDIF
              ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
              IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1
            !===================================
            ! Finite-difference check of the gradient (debug code) was removed
            ! here; see NEW_workproc.f90.bak4_20260924 if it is needed again.
            !===================================

            CASE (3:8)  ! Some kind of convergence has been reached
              ExitNeeded = .TRUE.

            CASE (9:10)  ! Function evaluation limit has been reached.
              ! This never suppose to happen because we
              ! count the number of function evaluations ourselves.
              ExitNeeded = .TRUE.


            CASE DEFAULT
              ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
              ! and >=14 for a bad input. None of those match a case above,
              ! so without this the loop would call DRMNG again for ever.
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                  'Warning WC0137 in BasisEnlG: DRMNG returned IV(1) =', IV(1)
                IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
              ENDIF
              ExitNeeded = .TRUE.
            ENDSELECT

            ! A warning, not an abort: the best point found so far is still
            ! usable, and the acceptance test below decides what to do with
            ! the block.
            IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                  'Warning WC0112 in BasisEnlG: number of failures in energy or gradient', &
                  'calculations during the optimization of nonlinear parameters', &
                  'reached the limit of', Glob_MaxEnergyFailsAllowed
              ENDIF
              ! call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
            ENDIF

            IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.

          ENDDO

        !------------------------------------------------------------------
        ! Anything else is a programming error, not an input choice
        !------------------------------------------------------------------
        ! Without this the SELECT fell through silently and the block was
        ! left at whatever point x_best happened to hold.
        !------------------------------------------------------------------
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0129 in BasisEnlG: unsupported value of OptimizationType', OptimizationType
            IF (Verbose >= 1) WRITE(*, *) 'Allowed values are 0 (no optimization) and 1 (powers and nonlinear parameters)'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)

        ENDSELECT  ! (OptimizationType)


        !------------------------------------------------------------------
        ! Step 4a: re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! EnergyGAM is used rather than EnergyGA because it also produces
        ! Glob_c, which the linear-coefficient test below needs.
        !
        ! A failure here is NOT fatal: the block is simply rejected and
        ! regenerated, the same response as a bad overlap.
        !------------------------------------------------------------------
        ! Calculate the energy and the linear coefficients at the best point found
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, nfru+i) = x_best((i-1)*npt+1:i*npt)
        ENDDO

        IsEnergyBad = .FALSE.
        Evalue = EnergyGAM(nfrup1, K, .TRUE., ErrCode)
        !!We run EnergyGA again because EnergyGAM might give slightly different
        !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
        ! Glob_CurrEnergy=EnergyGA(nfrup1,K,.false.,ErrCode)
        IF (ErrCode == 0) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          ! Reject the generated basis function(s) and make another attempt to
          ! generate them, just like it is done when the overlap or the linear
          ! parameters are found to be bad
          IsEnergyBad = .TRUE.
          Glob_CurrEnergy = E_init
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0113 in BasisEnlG: failed to evaluate energy after the optimization', &
              'of nonlinear parameters. Generated basis function(s) are rejected', &
              'and a new attempt to generate them will be made'
          ENDIF
        ENDIF


        !------------------------------------------------------------------
        ! Step 4b: pair overlaps
        !------------------------------------------------------------------
        ! Every pair involving a NEW function is checked against the whole
        ! basis (lower triangle of Glob_S). Glob_CurrEnergy is rolled back to
        ! E_init so the next attempt measures its improvement against the basis
        ! before this block. Disabled by OverlapThreshold <= 0.
        !------------------------------------------------------------------
        IsOverlapBad = .FALSE.
        IF (OverlapThreshold > ZERO) THEN
          ii = 0
          DO i = nfrup1, K
            DO j = 1, i-1
              IF (ABS(Glob_S(i, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                Glob_CurrEnergy = E_init
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0110: overlap of the following functions exceeds threshold. ', &
                      'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') ii, ':', i, j, '    S='
                  CALL writerealadv(6, Glob_S(i, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Step 4c: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..K, not just the new functions:
        ! adding a function can blow up the coefficient of an old one, and
        ! that is exactly the near-linear-dependence this test is meant to
        ! catch.
        !
        ! Disabled by passing LinCoeffThreshold <= 0.
        !------------------------------------------------------------------
        ! Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad = .FALSE.
        IF (LinCoeffThreshold > ZERO) THEN
          ii = 0
          DO i = 1, K
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              Glob_CurrEnergy = E_init
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0111: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', i, '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Step 4d: shape of the new functions after the Young operator
        !------------------------------------------------------------------
        ! C = sum|c_k S_k| / |<phi|Y+Y|phi>| of each new function. Above
        ! Glob_MaxSelfOverlapCancel the operator has almost annihilated the
        ! function: what survives is round-off, which the energy test cannot
        ! tell from a genuine improvement. Such a block is redrawn as often as
        ! necessary; these redraws do not count against the attempt budget.
        !------------------------------------------------------------------
        IsShapeBad = .FALSE.
        ii = 0
        DO i = nfrup1, K
          Cfac = SelfOverlapCancellation(Glob_PWR(i), Glob_NonlinParam(1:npt, i), ShapeSsum, ShapeSabs)
          IF (Cfac > Glob_MaxSelfOverlapCancel) THEN
            ii = ii+1
            IsShapeBad = .TRUE.
            Glob_CurrEnergy = E_init
            IF (Glob_ProcID == 0) THEN
              IF ((ii == 1) .AND. (Verbose >= 1)) THEN
                WRITE(*, '(1x,a,a,es9.2,a)') 'Warning WC0114 in BasisEnlG: the Young operator nearly annihilates ', &
                  'the following function(s), C > ', Glob_MaxSelfOverlapCancel, &
                  '. Generated basis function(s) are rejected and a new attempt to generate them will be made'
              ENDIF
              IF (Verbose >= 1) WRITE(*, '(1x,i6,a1,i6,a8,i5,a5,es10.3,a10,es10.3)') ii, ':', i, '   power', &
                Glob_PWR(i), '   C=', Cfac, '   S_raw=', ShapeSsum
            ENDIF
          ENDIF
        ENDDO
        IF (IsShapeBad) THEN
          NumOfShapeRedraws = NumOfShapeRedraws+1
          IF ((Glob_ProcID == 0) .AND. (Verbose >= 1) .AND. (MOD(NumOfShapeRedraws, 50) == 0)) &
            WRITE(*, '(1x,a,i0,a)') 'BasisEnlG: ', NumOfShapeRedraws, ' blocks redrawn so far because of the shape test'
        ENDIF

        IF (.NOT. IsShapeBad) AttemptToGetGoodFunc = AttemptToGetGoodFunc+1

      ENDDO  ! (IsOverlapBad.or.IsAnyLinCoeffBad).and. &
      ! (AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)


      !==================================================================
      ! Report the accepted block
      !==================================================================
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) THEN
        WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
        ENDIF
        WRITE (*, *) 'E=', Glob_CurrEnergy
        DO i = 1, nfo
          WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', Glob_PWR(nfru+i)
          CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
        ENDDO
      ENDIF

      !==================================================================
      ! Generator statistics
      !==================================================================
      ! Distance of the accepted parameters from their prototype, accumulated
      ! per generator method and averaged at the end of the run. wbfu (first
      ! prototype) and wmu (method) come from GenerateTrialParam for the
      ! winning candidate; rank 0 only. Skipped for the first blocks
      ! (nfru<=nfo), where prototype and new window overlap.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN

        IF (wmu == 1) THEN
          rgm1_counter = rgm1_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms1 = ms1+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF (wmu == 2) THEN
          rgm2_counter = rgm2_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms2 = ms2+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

      !==================================================================
      ! Commit the block
      !==================================================================
      ! The history of a newly added function starts empty: no cycles done
      ! and no energy evaluations spent on full optimization yet. The
      ! numbering is the identity here, which is the condition
      ! SortBasisFuncAndMatElem requires.
      !
      ! SaveResults runs on rank 0 only and every iteration, so an
      ! interrupted run can be resumed from the last completed block.
      !------------------------------------------------------------------
      Glob_CurrBasisSize = K

      DO i = 1, nfo
        Glob_History(nfru+i)%Energy = Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone = 0
        Glob_History(nfru+i)%InitFuncAtLastStep = 0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt = 0
      ENDDO

      DO i = 1, nfo
        Glob_FuncNum(nfru+i) = nfru+i
      ENDDO

      IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')

    ENDDO
    ! Main loop ends here


    !==================================================================
    ! Hand H and S to the next step and release everything
    !==================================================================
    ! The swap file lets the next BBOP step skip recomputing the matrix
    ! elements of the basis just built.
    !
    ! Deallocation is in the reverse of the allocation order throughout.
    !------------------------------------------------------------------
    CALL StoreMatricesInSwapFile()

    ! Deallocate arrays used by DRMNG
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)

    ! Deallocate workspace
    DEALLOCATE(ZIndOptSequence)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x)
    DEALLOCATE(ZIndSetBest)
    DEALLOCATE(ZIndSet)
    DEALLOCATE(ParSetBest)
    DEALLOCATE(ParSet)

    ! Deallocate workspace for EnergyGB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for DSYGVX
    DEALLOCATE(Glob_IWorkForDSYGVX)
    DEALLOCATE(Glob_WorkForDSYGVX)

    ! Deallocate global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_diagH)
    DEALLOCATE(GLob_S)
    DEALLOCATE(Glob_H)


    !==================================================================
    ! Closing summary
    !==================================================================
    ! The averages are per nonlinear parameter, hence the division by
    ! npt as well as by the number of times the method was used. Guarded
    ! against the counter being zero, which is the case whenever one of
    ! the two methods never ran.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *) 'Random selection statistics:'
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 1 of generating basis functions was used', rgm1_counter, 'times'
      IF ((rgm1_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms1/(npt*rgm1_counter)
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 2 of generating basis functions was used', rgm2_counter, 'times'
      IF ((rgm2_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms2/(npt*rgm2_counter)
      IF (Verbose >= 2) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlG finished'
    ENDIF


  END SUBROUTINE BasisEnlG


  SUBROUTINE BasisEnlI(Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval, &
                       OverlapThreshold, LinCoeffThreshold)
    !==================================================================
    ! Subroutine BasisEnlI
    !==================================================================
    ! Enlarges the basis from Kstart-1 to Kstop functions, Kstep at a time,
    ! with the GSEP solved by INVERSE ITERATION ('I'); a twin of BasisEnlG
    ! (same stages and acceptance tests, documented there). Differences:
    ! Glob_H holds the SHIFTED matrix H - Glob_ApproxEnergy*S and Glob_invD
    ! the LDL' diagonal (no Glob_diagH); the workspace is Glob_WorkForGSEPIIS
    ! and Glob_LastEigvector; EnergyIA/EnergyIAM/EnergyIB replace
    ! EnergyGA/EnergyGAM/EnergyGB; v_good keeps the last converged
    ! eigenvector (see its allocation); the average number of inverse
    ! iterations per solve is reported per block. The shift is re-anchored
    ! per block (RefreshINVITShift) and placed on eigenvalue
    ! Glob_WhichEigenvalue at the start of the step
    ! (RetargetShiftToEigenvalue) when Glob_EigIdxTargeting is 1.
    ! OptimizationType: 0 keep the selected parameters, 1 optimize; else
    ! EC0139.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: Kstart, Kstop, Kstep, NTrials, OptimizationType, MaxEnergyEval
    REAL(wp), INTENT(IN) :: OverlapThreshold, LinCoeffThreshold

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    ! Local variables:
    INTEGER  :: i, j, K, AttemptToGetGoodFunc, ii, jj, jbest
    INTEGER  :: np, npt, nfo, nfa, nfru, nfrup1, nvmax, nv
    INTEGER  :: OpenFileErr, ErrCode, NumOfFailures, NumOfEnergyEval, NumOfGradEval
    LOGICAL  :: IsSwapFileOK, IsEnergyImproved, ExitNeeded
    LOGICAL  :: IsOverlapBad, IsAnyLinCoeffBad, IsEnergyBad
    LOGICAL  :: IsShapeBad
    INTEGER  :: NumOfShapeRedraws
    REAL(wp) :: Cfac, ShapeSsum, ShapeSabs
    INTEGER  :: wbfu_t, wmu_t, wbfu, wmu, rgm1_counter, rgm2_counter
    REAL(wp) :: ms1, ms2
    REAL(wp) :: Evalue, E_init, E_best
    REAL(wp) :: t

    ! Largest legal premultiplier power: the greatest EVEN value not
    ! exceeding Glob_MaxPowerAllowed. Powers are always even in this
    ! frame - GenerateTrialParam only ever produces even ones - so the
    ! power scan below steps through 2,4,...,PWRMax.
    INTEGER, PARAMETER :: PWRMax = 2*(Glob_MaxPowerAllowed/2)

    ! Candidate block and the best candidate block found so far
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: ParSet, ParSetBest
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: ZIndSet, ZIndSetBest

    ! The optimization variables: the nonlinear parameters of the block
    ! laid out as one flat vector of nfo*npt elements
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_best, grad

    ! Last eigenvector that came out of a SUCCESSFUL solve - see the
    ! note where it is allocated
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: v_good

    INTEGER, ALLOCATABLE, DIMENSION(:) :: ZIndOptSequence

    ! Arrays used by DRMNG
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D, V, V_init
    INTEGER, PARAMETER                  :: LIV = 60
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV
    INTEGER                             :: ALG

    ! Allocatable work space
    ! *** These six are NOT REFERENCED anywhere in this routine, and
    ! neither is OpenFileErr above. They are left over from when the
    ! reallocation and the swap-file handling were written out inline
    ! here; both now live in their own routines, which declare their
    ! own copies.
    REAL(wp), ALLOCATABLE, DIMENSION(:)               :: WorkBuffReal
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: WorkBuffInt
    TYPE(Glob_HistoryStep), ALLOCATABLE, DIMENSION(:) :: TempHistory
    REAL(wp), ALLOCATABLE, DIMENSION(:, :)            :: TempParam
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempZInd
    INTEGER, ALLOCATABLE, DIMENSION(:)                :: TempFunc
    !====================================================
    ! These variables are used when a finite difference gradient is computed
    ! real(wp),allocatable,dimension(:)     ::    fx,fgrad
    ! real(wp)                                    deltax,Evalue1
    !====================================================


    !==================================================================
    ! Announce the step
    !==================================================================
    ! wbfu_t/wmu_t are filled by GenerateTrialParam on rank 0 and carry
    ! which existing function was used as prototype and which generator
    ! method produced the candidate. They feed the statistics printed at
    ! the very end, so they are cleared once here.
    !------------------------------------------------------------------
    wbfu_t = 0
    wmu_t = 0

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlI started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstart =', Kstart
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstop =', Kstop
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Kstep =', Kstep
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'OptimizationType =', OptimizationType
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval =', MaxEnergyEval
    ENDIF


    !==================================================================
    ! Global state this routine works under
    !==================================================================
    ! Glob_nfru/Glob_nfo/Glob_nfa define the optimization window (nfru
    ! frozen functions, nfo being optimized, nfa the last function) that
    ! the energy routines and Glob_D are indexed against; they are reset
    ! per iteration in the main loop. Overlap penalties are off: an overlap
    ! violation REJECTS the block instead of penalizing the energy.
    !------------------------------------------------------------------
    Glob_GSEPSolutionMethod = 'I'
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_nfa = Kstart+Kstep
    Glob_nfo = Kstep
    Glob_HSLeadDim = Kstop
    Glob_HSBuffLen = Kstop*Kstep
    np = Glob_np
    npt = Glob_npt
    nfo = Glob_nfo
    nfa = Glob_nfa
    nvmax = Kstep*Glob_npt

    ! Reallocate arrays that contain the information
    ! about basis functions and optimization process.
    CALL ReallocateBasisFuncData(Kstop, Glob_CurrBasisSize)


    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! Only the LOWER triangles of Glob_H and Glob_S are meaningful; the
    ! upper ones are used as scratch by the permutation routines.
    !
    ! No Glob_diagH here: on the 'I' path the diagonal of H stays inside
    ! Glob_H. Glob_invD takes its place, holding the diagonal of the
    ! LDL' factorization that inverse iteration works with.
    !------------------------------------------------------------------
    ! Allocate some global arrays
    ALLOCATE(Glob_H(Kstop, Kstop))
    ALLOCATE(Glob_S(Kstop, Kstop))
    ALLOCATE(Glob_diagS(Kstop))
    ALLOCATE(Glob_invD(Kstop))
    ALLOCATE(Glob_c(Kstop))
    ALLOCATE(Glob_D(2*npt, Kstep, Kstop))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    !==================================================================
    ! Inverse-iteration workspace
    !==================================================================
    ! Glob_LastEigvector is the STARTING VECTOR of the next inverse
    ! iteration. Carrying the previous solution over is what makes the
    ! method cheap: consecutive evaluations differ only slightly, so a
    ! few iterations suffice. It is primed with all ones because there
    ! is nothing better to start from.
    !------------------------------------------------------------------
    ! Allocate workspace for subroutine GSEPIIS, which is called
    ! inside EnergyIA, EnergyIAM, and EnergyIB
    ALLOCATE(Glob_WorkForGSEPIIS(Kstop))
    ALLOCATE(Glob_LastEigvector(Kstop))
    Glob_LastEigvector(1:Kstop) = ONE

    !------------------------------------------------------------------
    ! ... and the fallback copy of it
    !------------------------------------------------------------------
    ! A failed (diverged or half-converged) inverse iteration leaves
    ! Glob_LastEigvector in a state that poisons every later solve, so the
    ! vector of the last SUCCESSFUL solve is kept in v_good and restored
    ! whenever ErrCode comes back nonzero. No counterpart on the 'G' path,
    ! where DSYGVX starts from scratch every time.
    !------------------------------------------------------------------
    ALLOCATE(v_good(Kstop))
    v_good(1:Kstop) = ONE

    ! Allocate workspace for EnergyIB
    ALLOCATE(Glob_WkGR(Kstep*npt))

    !------------------------------------------------------------------
    ! Local workspace
    !------------------------------------------------------------------
    ! Allocate workspace
    ALLOCATE(ParSet(npt, Kstep))
    ALLOCATE(ParSetBest(npt, Kstep))
    ALLOCATE(ZIndSet(Kstep))
    ALLOCATE(ZIndSetBest(Kstep))
    ALLOCATE(x(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(ZIndOptSequence(Kstep))


    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next (1 energy,
    ! 2 gradient) and is called again, which suits an energy that is a
    ! collective operation over all processes.
    ! LV is the documented size of the V work array plus one.
    !------------------------------------------------------------------
    nvmax = npt*Kstep
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))

    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next (1 energy,
    ! 2 gradient) and is called again, which suits an energy that is a
    ! collective operation over all processes.
    ! The parameters are set once, outside the main loop.
    !------------------------------------------------------------------

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): iteration and function-evaluation limits, set out of
    ! the way because the budget is enforced by MaxEnergyEval below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    ! V(35)=0.1*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed


    !==================================================================
    ! Initial energy
    !==================================================================
    ! The swap file, when it is valid, carries H and S for the basis we
    ! start from, so the matrix elements do not have to be recomputed -
    ! hence .false. for AreMatElemNeeded in that branch.
    !
    ! With Kstart==1 there is no basis yet and no energy to compute, so
    ! Glob_CurrEnergy is primed with HUGE() and the first candidate to
    ! produce a finite energy wins.
    !------------------------------------------------------------------
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    WrkP_WrongStateCount = 0

    ! Calculating the initial energy
    ErrCode = 0
    IF (Kstart > 1) THEN
      IF (IsSwapFileOK) THEN
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Solving eigenvalue problem...'
        Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
      ELSE
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem...'
        Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .TRUE., ErrCode)
      ENDIF
    ELSE
      Glob_CurrEnergy = HUGE(Glob_CurrEnergy)
    ENDIF
    !==================================================================
    ! Put the inverse-iteration shift on the requested eigenvalue
    !==================================================================
    ! Once per BBOP step; inert unless Glob_EigIdxTargeting==1. It runs
    ! BEFORE the fatal check below: a shift sitting between two
    ! eigenvalues at nearly equal distance is the usual reason the first
    ! solve does not converge, and moving the shift is precisely the
    ! cure.
    !------------------------------------------------------------------
    IF ((Kstart > 1) .AND. (Glob_EigIdxTargeting == 1)) THEN
      CALL RetargetShiftToEigenvalue(Glob_CurrBasisSize, 'BasisEnlI')
      Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
    ENDIF

    ! The second test catches a solve that converged on a level other
    ! than WHICH_EIGENVALUE: EnergyIA returns it with ErrCode=0 and the
    ! rejection sentinel as the energy - see IsRequestedEigenstate.
    IF ((ErrCode /= 0) .OR. ((Kstart > 1) .AND. (ABS(Glob_CurrEnergy) > 1.0E10_wp))) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0135 in BasisEnlI: initial energy cannot be computed'
        IF (ErrCode == 0) WRITE(*, *) '(inverse iteration converged on a level other than WHICH_EIGENVALUE)'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy

    rgm1_counter = 0
    rgm2_counter = 0
    ms1 = ZERO
    ms2 = ZERO
    K = Kstart-1


    !==================================================================
    ! MAIN LOOP - one block of up to Kstep functions per iteration
    !==================================================================
    ! Main loop begins here
    DO WHILE (K < Kstop)

      !--------------------------------------------------------------
      ! Size and place the window for this block
      !--------------------------------------------------------------
      ! nfru = functions that stay frozen, nfo = functions being added,
      ! K = the new basis size. The last block is short when Kstop-K is
      ! less than Kstep.
      IF (K+Kstep <= Kstop) THEN
        nfo = Kstep
        nfru = K
        K = K+Kstep
      ELSE
        nfo = Kstop-K
        nfru = K
        K = Kstop
      ENDIF

      CALL linalg_setparam(K)  ! reset linalg flags to account for changes in the basis size

      Glob_nfa = K
      Glob_nfru = nfru
      Glob_nfo = nfo
      nfrup1 = nfru+1
      nv = nfo*npt
      E_init = Glob_CurrEnergy

      ! Per-block counters behind the "average iterations in GSEPIIS"
      ! line printed further down. Counter1 counts the solves and
      ! Counter2 accumulates their iterations; both are bumped inside
      ! EnergyIA, EnergyIAM and EnergyIB.
      Glob_InvItTempCounter1 = 0
      Glob_InvItTempCounter2 = 0

      !--------------------------------------------------------------
      ! Re-anchor the inverse-iteration shift on the energy reached so
      ! far, once per block. The shift must
      ! stay fixed WHILE a block is being added - the trial solves reuse
      ! the LDL^T of the leading nfru block, valid for one shift only -
      ! but between blocks there is nothing to preserve, and a shift
      ! frozen for the whole step is what made every solve run to the
      ! iteration cap by the tenth added function. Skipped on the first
      ! block of the step, where RetargetShiftToEigenvalue has just
      ! placed the shift. The refresh invalidates the stored
      ! factorization, so the solve that follows starts from row 1.
      !--------------------------------------------------------------
      IF (nfru > Kstart-1) THEN
        CALL RefreshINVITShift(nfru)
        Evalue = EnergyIA(1, nfru, .FALSE., ErrCode)
        IF ((ErrCode == 0) .AND. (ABS(Evalue) < 1.0E10_wp)) Glob_CurrEnergy = Evalue
        E_init = Glob_CurrEnergy
      ENDIF

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Current basis size is', Glob_CurrBasisSize
        IF (nfo > 1) THEN
          WRITE(*, '(1x,a,1x,i0,a,i0)') 'Selecting functions', nfrup1, '-', K
        ELSE
          WRITE(*, '(1x,a,1x,i0)') 'Selecting function', K
        ENDIF
      ENDIF


      !==================================================================
      ! ACCEPTANCE LOOP - regenerate the block until it is acceptable
      !==================================================================
      ! The three flags are primed so that the loop always runs at least
      ! once. It ends when the block passes the energy, overlap and linear
      ! coefficient tests, or when the attempt budget runs out - in which
      ! case the last block is kept regardless, which is deliberate: a
      ! basis that is slightly too linearly dependent is better than no
      ! progress at all. A block in which the Young operator nearly
      ! annihilates a new function (IsShapeBad, C > Glob_MaxSelfOverlapCancel)
      ! is redrawn without limit and without spending the attempt budget.
      !------------------------------------------------------------------
      IsOverlapBad = .TRUE.
      IsAnyLinCoeffBad = .TRUE.
      IsEnergyBad = .FALSE.
      IsShapeBad = .FALSE.
      NumOfShapeRedraws = 0
      AttemptToGetGoodFunc = 1

      DO WHILE (((IsOverlapBad .OR. IsAnyLinCoeffBad .OR. IsEnergyBad) .AND. &
                 (AttemptToGetGoodFunc <= Glob_BadOverlapOrLinCoeffLim)) .OR. IsShapeBad)

        !------------------------------------------------------------------
        ! Step 1: stochastic selection of the block
        !------------------------------------------------------------------
        ! GenerateTrialParam runs on rank 0 (it consumes the random stream) and
        ! the result is broadcast; every rank evaluates the SAME candidate.
        ! Every successful candidate refreshes v_good and every failed one
        ! restores from it, so the next inverse iteration always starts from a
        ! converged vector.
        !------------------------------------------------------------------
        NumOfFailures = 0
        IsEnergyImproved = .FALSE.
        wbfu = 0
        wmu = 0

        DO i = 1, NTrials
          IF (Glob_ProcID == 0) CALL GenerateTrialParam(nfo, ParSet, ZIndSet, wbfu_t, wmu_t)
          CALL MPI_BCAST(ParSet, npt*nfo, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_BCAST(ZIndSet, nfo, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          Glob_NonlinParam(1:npt, nfrup1:K) = ParSet(1:npt, 1:nfo)
          Glob_PWR(nfrup1:K) = ZIndSet(1:nfo)
          Evalue = EnergyIA(nfrup1, K, .TRUE., ErrCode)
          IF (ErrCode == 0) THEN
            v_good(1:K) = Glob_LastEigvector(1:K)
            IF (Evalue < Glob_CurrEnergy) THEN
              Glob_CurrEnergy = Evalue
              ParSetBest(1:npt, 1:nfo) = ParSet(1:npt, 1:nfo)
              ZIndSetBest(1:nfo) = ZIndSet(1:nfo)
              IsEnergyImproved = .TRUE.
              wbfu = wbfu_t
              wmu = wmu_t
            ENDIF
          ELSE
            NumOfFailures = NumOfFailures+1
            ! Restore the last good eigenvector as the vector left by the failed
            ! inverse iteration process may be unusable as a starting vector
            Glob_LastEigvector(1:K) = v_good(1:K)
          ENDIF
        ENDDO

        ! Too many candidates the solver could not handle at all: the
        ! basis is in a state inverse iteration cannot work with, and
        ! more trials will not fix it.
        IF (NumOfFailures*ONE/NTrials > Glob_MaxFracOfTrialFailsAllowed) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0136 in BasisEnlI: the number of eigenvalue problem solution failures'
            WRITE(*, *) 'in random selection process exceeded limit'
            WRITE(*, '(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials, '%'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! The candidates were all solvable but none lowered the energy.
        IF (.NOT. (IsEnergyImproved)) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0137 in BasisEnlI: random selection did not result'
            WRITE(*, *) 'in any energy improvement'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! Put the winning candidate back in place. The re-solve is needed
        ! because the last candidate evaluated is not in general the best
        ! one, so Glob_H/Glob_S hold the wrong block.
        Glob_NonlinParam(1:npt, nfrup1:K) = ParSetBest(1:npt, 1:nfo)
        Glob_PWR(nfrup1:K) = ZIndSetBest(1:nfo)
        Evalue = EnergyIA(nfrup1, K, .TRUE., ErrCode)
        IF (ErrCode == 0) THEN
          Glob_CurrEnergy = Evalue
          v_good(1:K) = Glob_LastEigvector(1:K)
        ELSE
          Glob_LastEigvector(1:K) = v_good(1:K)
        ENDIF

        IF (Glob_ProcID == 0) THEN
          WRITE (*, '(1x,a)', ADVANCE='no') 'E='
          CALL writereal(6, Glob_CurrEnergy)
          IF (Verbose >= 2) WRITE (*, '(5x,a,1x,i0)') 'prototype function is', wbfu
          DO i = 1, nfo
            WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', ZIndSetBest(i)
            CALL writerealarradv(6, ParSetBest(1:npt, i), npt)
          ENDDO
          IF (Verbose >= 1) WRITE (*, *) 'Optimizing nonlinear parameters'
        ENDIF


        !------------------------------------------------------------------
        ! Prime the best point found
        !------------------------------------------------------------------
        ! Every branch of the SELECT below must leave x_best defined, because
        ! the energy and the linear coefficients are recomputed at x_best right
        ! after it; the randomly selected point is also what OptimizationType 0
        ! keeps.
        !------------------------------------------------------------------
        DO i = 1, nfo
          x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
        ENDDO
        E_best = Glob_CurrEnergy
        x_best(1:nfo*npt) = x(1:nfo*npt)


        ! Optimization of the nonlinear parameters
        SELECT CASE (OptimizationType)

        !------------------------------------------------------------------
        ! OptimizationType 0: keep the randomly selected block
        !------------------------------------------------------------------
        CASE (0)  ! No optimization
        ! x and x_best already hold the selected point - nothing to do.
        ! Note that this is NOT what OptimizationType==0 meant in
        ! workproc.f90, where anything other than 1 selected the DMNG
        ! optimization. A driver that wants the parameters optimized
        ! must now pass 1.

        !------------------------------------------------------------------
        ! OptimizationType 1: optimize powers, then nonlinear parameters
        !------------------------------------------------------------------
        CASE (1)

          !------------------------------------------------------------------
          ! Step 2: premultiplier powers, one function at a time
          !------------------------------------------------------------------
          ! Each function in the window is tried with every legal EVEN power
          ! 2,4,...,PWRMax and the best is kept; the functions are visited in
          ! random order. This exhaustive scan costs PWRMax/2 - 1 energy
          ! evaluations per function (a 1-D minimizer or a window around the
          ! current power would be cheaper). Skipped when the power is pinned
          ! by Glob_IsIndexFixed.
          !------------------------------------------------------------------
          NumOfFailures = 0
          IF (.NOT. Glob_IsIndexFixed) THEN
            ! Generate a random sequence which will define the order in which
            ! Z-indices should be optimized (one index at a time)
            CALL GenerateRndIntSeq(nfo, ZIndOptSequence)
            ! Loop where Z-indices are optimized. Note: there is some room for improvement
            ! here as I programmed it in a simple way when all matrix element of functions
            ! nfrup1 through K are computed each time while it is not always necessary.
            DO i = 1, nfo
              ii = ZIndOptSequence(i)
              j = Glob_PWR(nfru+ii)
              jbest = j
              DO jj = 2, PWRMax, 2
                IF (jj /= j) THEN
                  Glob_PWR(nfru+ii) = jj
                  Evalue = EnergyIA(nfrup1, K, .TRUE., ErrCode)
                  IF (ErrCode /= 0) THEN
                    ! Restore the best power known so far before carrying on,
                    ! so a failed trial never leaves a bad power behind.
                    NumOfFailures = NumOfFailures+1
                    Glob_PWR(nfru+ii) = jbest
                    IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
                      IF (Glob_ProcID == 0) THEN
                        WRITE(*, *) 'Error EC0138 in BasisEnlI: number of failures in energy calculations'
                        WRITE(*, *) 'during the optimization of Z-indicies exceeded limit'
                      ENDIF
                      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
                    ENDIF
                  ELSE
                    IF (Evalue < Glob_CurrEnergy) THEN
                      Glob_CurrEnergy = Evalue
                      jbest = jj
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
              Glob_PWR(nfru+ii) = jbest
            ENDDO
          ENDIF

          !------------------------------------------------------------------
          ! Step 3: nonlinear parameters, with DRMNG
          !------------------------------------------------------------------
          ! IV and V are restored from the copies made before the main loop, so
          ! each block starts the minimizer from its documented default state.
          !------------------------------------------------------------------
          ! Now we optimize nonlinear parameters

          ! Setting IV and V values as was in their initial copies
          IV(1:LIV) = IV_init(1:LIV)
          V(1:LV) = V_init(1:LV)

          nv = nfo*npt
          DO i = 1, nfo
            x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
          ENDDO

          !--------------------------------------------------------------
          ! Scale vector D for DRMNG
          !--------------------------------------------------------------
          ! Every variable gets the same scale t, taken from the RELATIVE energy
          ! gain of the block so far (a block that already moved the energy a lot
          ! may take larger steps), floored at 10000*epsilon. At the very start
          ! of a basis (nfru>=nfo fails, E_init is the HUGE() sentinel) t=1.
          IF (nfru >= nfo) THEN
            t = MAX(ABS((E_init-Glob_CurrEnergy))/(ABS(E_init)+ABS(Glob_CurrEnergy)), &
                    10000*EPSILON(Glob_CurrEnergy))
          ELSE
            t = ONE
          ENDIF
          ! if (Glob_ProcID==0) write(*,*) 'scaling coeff=',t !remove later
          DO i = 1, nfo
            ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
            DO j = 1, npt
              ! Make sure none of the D(i) will be zero or smaller than the threshold
              ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
              D(npt*(i-1)+j) = t
              ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
            ENDDO
          ENDDO

          ExitNeeded = .FALSE.
          NumOfFailures = 0
          NumOfEnergyEval = 0
          NumOfGradEval = 0
          IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.
          E_best = Glob_CurrEnergy
          x_best(1:nfo*npt) = x(1:nfo*npt)

          !------------------------------------------------------------------
          ! The reverse-communication loop
          !------------------------------------------------------------------
          ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it wants:
          ! 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its evaluation
          ! limit. On a failed evaluation the energy handed to DRMNG keeps its
          ! previous value, so the step counts as no reduction and the trust
          ! radius shrinks (IV(2), the TOOBIG flag, is deliberately left 0). The
          ! best point is tracked here because the last point DRMNG visits is not
          ! necessarily the lowest.
          ! On this path a failure also restores v_good.
          !------------------------------------------------------------------
          DO WHILE (.NOT. (ExitNeeded))

            IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
            CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

            SELECT CASE (IV(1))

            CASE (1)  ! Only energy is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              Evalue = EnergyIA(nfrup1, K, .TRUE., ErrCode)
              NumOfEnergyEval = NumOfEnergyEval+1
              IF (ErrCode /= 0) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
                ! Restore the last good eigenvector as the vector left by the failed
                ! inverse iteration process may be unusable as a starting vector
                Glob_LastEigvector(1:K) = v_good(1:K)
              ELSE
                Glob_CurrEnergy = Evalue
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nfo*npt) = x(1:nfo*npt)
                ENDIF
              ENDIF
              ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
              IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

            CASE (2)  ! Only gradient is needed
              CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
              DO i = 1, nfo
                Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
              ENDDO
              CALL EnergyIB(Evalue, grad, .TRUE., ErrCode)
              NumOfGradEval = NumOfGradEval+1
              IF (ErrCode /= 0) THEN
                NumOfFailures = NumOfFailures+1
                IV(2) = 1
                ! Restore the last good eigenvector as the vector left by the failed
                ! inverse iteration process may be unusable as a starting vector
                Glob_LastEigvector(1:K) = v_good(1:K)
              ELSE
                IF (Evalue < E_best) THEN
                  E_best = Evalue
                  x_best(1:nfo*npt) = x(1:nfo*npt)
                ENDIF
              ENDIF
              ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
              IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1
            !===================================
            ! Finite-difference check of the gradient (debug code) was removed
            ! here; see NEW_workproc.f90.bak4_20260924 if it is needed again.
            !===================================

            CASE (3:8)  ! Some kind of convergence has been reached
              ExitNeeded = .TRUE.

            CASE (9:10)  ! Function evaluation limit has been reached.
              ! This never suppose to happen because we
              ! count the number of function evaluations ourselves.
              ExitNeeded = .TRUE.


            CASE DEFAULT
              ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
              ! and >=14 for a bad input. None of those match a case above,
              ! so without this the loop would call DRMNG again for ever.
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                  'Warning WC0138 in BasisEnlI: DRMNG returned IV(1) =', IV(1)
                IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
              ENDIF
              ExitNeeded = .TRUE.
            ENDSELECT

            ! A warning, not an abort: the best point found so far is still
            ! usable, and the acceptance test below decides what to do with
            ! the block.
            IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
              IF (Glob_ProcID == 0) THEN
                IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                  'Warning WC0117 in BasisEnlI: number of failures in energy or gradient', &
                  'calculations during the optimization of nonlinear parameters', &
                  'reached the limit of', Glob_MaxEnergyFailsAllowed
              ENDIF
              ! call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
            ENDIF

            IF (NumOfEnergyEval >= MaxEnergyEval) ExitNeeded = .TRUE.

          ENDDO

        !------------------------------------------------------------------
        ! Anything else is a programming error, not an input choice
        !------------------------------------------------------------------
        ! Without this the SELECT fell through silently and the block was
        ! left at whatever point x_best happened to hold.
        !------------------------------------------------------------------
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Error EC0139 in BasisEnlI: unsupported value of OptimizationType', OptimizationType
            IF (Verbose >= 1) WRITE(*, *) 'Allowed values are 0 (no optimization) and 1 (powers and nonlinear parameters)'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)

        ENDSELECT  ! (OptimizationType)


        !------------------------------------------------------------------
        ! Step 4a: re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! EnergyIAM is used rather than EnergyIA because it also produces
        ! Glob_c, which the linear-coefficient test below needs.
        !
        ! A failure here is NOT fatal: the block is simply rejected and
        ! regenerated, the same response as a bad overlap - and the starting
        ! vector is restored from v_good so the next attempt is not poisoned.
        !------------------------------------------------------------------
        ! Calculate the energy and the linear coefficients at the best point found
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, nfru+i) = x_best((i-1)*npt+1:i*npt)
        ENDDO

        IsEnergyBad = .FALSE.
        Evalue = EnergyIAM(nfrup1, K, .TRUE., ErrCode)
        IF (ErrCode == 0) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          ! Reject the generated basis function(s) and make another attempt to
          ! generate them, just like it is done when the overlap or the linear
          ! parameters are found to be bad
          IsEnergyBad = .TRUE.
          Glob_CurrEnergy = E_init
          Glob_LastEigvector(1:K) = v_good(1:K)
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0118 in BasisEnlI: failed to evaluate energy after the optimization', &
              'of nonlinear parameters. Generated basis function(s) are rejected', &
              'and a new attempt to generate them will be made'
          ENDIF
        ENDIF


        !------------------------------------------------------------------
        ! Step 4b: pair overlaps
        !------------------------------------------------------------------
        ! Every pair involving a NEW function is checked against the whole
        ! basis (lower triangle of Glob_S). Glob_CurrEnergy is rolled back to
        ! E_init so the next attempt measures its improvement against the basis
        ! before this block. Disabled by OverlapThreshold <= 0.
        !------------------------------------------------------------------
        IsOverlapBad = .FALSE.
        IF (OverlapThreshold > ZERO) THEN
          ii = 0
          DO i = nfrup1, K
            DO j = 1, i-1
              IF (ABS(Glob_S(i, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                Glob_CurrEnergy = E_init
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0115: overlap of the following functions exceeds threshold. ', &
                      'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') ii, ':', i, j, '    S='
                  CALL writerealadv(6, Glob_S(i, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Step 4c: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..K, not just the new functions:
        ! adding a function can blow up the coefficient of an old one, and
        ! that is exactly the near-linear-dependence this test is meant to
        ! catch.
        !
        ! Disabled by passing LinCoeffThreshold <= 0.
        !------------------------------------------------------------------
        ! Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad = .FALSE.
        IF (LinCoeffThreshold > ZERO) THEN
          ii = 0
          DO i = 1, K
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              Glob_CurrEnergy = E_init
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0116: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', i, '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Step 4d: shape of the new functions after the Young operator
        !------------------------------------------------------------------
        ! C = sum|c_k S_k| / |<phi|Y+Y|phi>| of each new function. Above
        ! Glob_MaxSelfOverlapCancel the operator has almost annihilated the
        ! function: what survives is round-off, which the energy test cannot
        ! tell from a genuine improvement. Such a block is redrawn as often as
        ! necessary; these redraws do not count against the attempt budget.
        !------------------------------------------------------------------
        IsShapeBad = .FALSE.
        ii = 0
        DO i = nfrup1, K
          Cfac = SelfOverlapCancellation(Glob_PWR(i), Glob_NonlinParam(1:npt, i), ShapeSsum, ShapeSabs)
          IF (Cfac > Glob_MaxSelfOverlapCancel) THEN
            ii = ii+1
            IsShapeBad = .TRUE.
            Glob_CurrEnergy = E_init
            IF (Glob_ProcID == 0) THEN
              IF ((ii == 1) .AND. (Verbose >= 1)) THEN
                WRITE(*, '(1x,a,a,es9.2,a)') 'Warning WC0119 in BasisEnlI: the Young operator nearly annihilates ', &
                  'the following function(s), C > ', Glob_MaxSelfOverlapCancel, &
                  '. Generated basis function(s) are rejected and a new attempt to generate them will be made'
              ENDIF
              IF (Verbose >= 1) WRITE(*, '(1x,i6,a1,i6,a8,i5,a5,es10.3,a10,es10.3)') ii, ':', i, '   power', &
                Glob_PWR(i), '   C=', Cfac, '   S_raw=', ShapeSsum
            ENDIF
          ENDIF
        ENDDO
        IF (IsShapeBad) THEN
          NumOfShapeRedraws = NumOfShapeRedraws+1
          IF ((Glob_ProcID == 0) .AND. (Verbose >= 1) .AND. (MOD(NumOfShapeRedraws, 50) == 0)) &
            WRITE(*, '(1x,a,i0,a)') 'BasisEnlI: ', NumOfShapeRedraws, ' blocks redrawn so far because of the shape test'
        ENDIF

        IF (.NOT. IsShapeBad) AttemptToGetGoodFunc = AttemptToGetGoodFunc+1

      ENDDO  ! (IsOverlapBad.or.IsAnyLinCoeffBad).and. &
      ! (AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)


      !==================================================================
      ! Report the accepted block
      !==================================================================
      ! The last line is the average number of inverse iterations per
      ! solve over this block - the figure to watch if the fixed shift
      ! starts costing convergence. Counter1 is at least NTrials here,
      ! since every trial ran one solve, so the division is safe for any
      ! sensible NTrials.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) THEN
        WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
        ENDIF
        WRITE (*, *) 'E=', Glob_CurrEnergy
        DO i = 1, nfo
          WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') nfru+i, ':', Glob_PWR(nfru+i)
          CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
        ENDDO
        IF (Verbose >= 1) WRITE (*, *) 'Average number of iterations in GSEPIIS: ', &
          (Glob_InvItTempCounter2*ONE)/Glob_InvItTempCounter1
      ENDIF

      !==================================================================
      ! Generator statistics
      !==================================================================
      ! Distance of the accepted parameters from their prototype, accumulated
      ! per generator method and averaged at the end of the run. wbfu (first
      ! prototype) and wmu (method) come from GenerateTrialParam for the
      ! winning candidate; rank 0 only. Skipped for the first blocks
      ! (nfru<=nfo), where prototype and new window overlap.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN

        IF (wmu == 1) THEN
          rgm1_counter = rgm1_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms1 = ms1+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF (wmu == 2) THEN
          rgm2_counter = rgm2_counter+1
          IF (nfru > nfo) THEN
            DO i = 1, nfo
              DO j = 1, npt
                t = (Glob_NonlinParam(j, wbfu+i-1)-Glob_NonlinParam(j, nfru+i)) &
                   /Glob_NonlinParam(j, wbfu+i-1)
                ms2 = ms2+ABS(t)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

      !==================================================================
      ! Commit the block
      !==================================================================
      ! The history of a newly added function starts empty: no cycles done
      ! and no energy evaluations spent on full optimization yet. The
      ! numbering is the identity here, which is the condition
      ! SortBasisFuncAndMatElem requires.
      !
      ! SaveResults runs on rank 0 only and every iteration, so an
      ! interrupted run can be resumed from the last completed block.
      !------------------------------------------------------------------
      Glob_CurrBasisSize = K

      DO i = 1, nfo
        Glob_History(nfru+i)%Energy = Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone = 0
        Glob_History(nfru+i)%InitFuncAtLastStep = 0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt = 0
      ENDDO

      DO i = 1, nfo
        Glob_FuncNum(nfru+i) = nfru+i
      ENDDO

      IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')

    ENDDO
    ! Main loop ends here


    !==================================================================
    ! Hand H and S to the next step and release everything
    !==================================================================
    ! The swap file lets the next BBOP step skip recomputing the matrix
    ! elements of the basis just built. On this path Glob_H holds the
    ! SHIFTED matrix, and StoreMatricesInSwapFile adds the shift back
    ! before writing - see the note there.
    !
    ! Deallocation is in the reverse of the allocation order throughout.
    !------------------------------------------------------------------
    CALL StoreMatricesInSwapFile()

    ! Deallocate arrays used by DRMNG
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)

    ! Deallocate workspace
    DEALLOCATE(ZIndOptSequence)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x)
    DEALLOCATE(ZIndSetBest)
    DEALLOCATE(ZIndSet)
    DEALLOCATE(ParSetBest)
    DEALLOCATE(ParSet)

    ! Deallocate workspace for EnergyIB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for subroutine GSEPIIS, which is called
    ! inside EnergyIA, EnergyIAM, and EnergyIB
    DEALLOCATE(v_good)
    DEALLOCATE(Glob_LastEigvector)
    DEALLOCATE(Glob_WorkForGSEPIIS)

    ! Deallocate global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_invD)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(GLob_S)
    DEALLOCATE(Glob_H)


    !==================================================================
    ! Closing summary
    !==================================================================
    ! The averages are per nonlinear parameter, hence the division by
    ! npt as well as by the number of times the method was used. Guarded
    ! against the counter being zero, which is the case whenever one of
    ! the two methods never ran.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *) 'Random selection statistics:'
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 1 of generating basis functions was used', rgm1_counter, 'times'
      IF ((rgm1_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms1/(npt*rgm1_counter)
      IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Method 2 of generating basis functions was used', rgm2_counter, 'times'
      IF ((rgm2_counter /= 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ', ms2/(npt*rgm2_counter)
      IF (Verbose >= 2) WRITE(*, *)
      IF (WrkP_WrongStateCount > 0) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Trial points refused because inverse iteration landed on'
        WRITE(*, *) 'a level other than WHICH_EIGENVALUE =', Glob_WhichEigenvalue, &
                   ' :', WrkP_WrongStateCount
        WRITE(*, *)
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Routine BasisEnlI finished'
    ENDIF


  END SUBROUTINE BasisEnlI


  SUBROUTINE OptCycleG(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                       NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold, SavingFreq)
    !==================================================================
    ! Subroutine OptCycleG
    !==================================================================
    ! Improves an EXISTING basis (no functions added): a window of
    ! NumOfFuncToOpt functions sweeps FuncBegin..FuncEnd, advancing by
    ! NumOfFuncToShift per step (normally equal, so the windows tile the
    ! range), NumCycles times, optimizing the nonlinear parameters in the
    ! window with DSYGVX ('G') and DRMNG; OptCycleI is the inverse-iteration
    ! twin. The window is physically MOVED to the end of the basis
    ! (PermuteFunctions/PermuteMatrixElements and their 2-variants), so the
    ! energy routines only recompute the trailing block, and
    ! SortBasisFuncAndMatElem restores the order after each cycle.
    ! Glob_History(cbs)%CyclesDone and %InitFuncAtLastStep are written after
    ! every step and read on entry, so an interrupted run resumes. A step
    ! that ends with an unusable energy, an overlap above OverlapThreshold
    ! or a linear coefficient above LinCoeffThreshold is UNDONE (parameters
    ! back to x_init) and the sweep moves on; a threshold <= 0 disables
    ! its test. Premultiplier powers are not optimized here.
    ! Arguments: K (NOT REFERENCED; Glob_CurrBasisSize is used), FuncBegin,
    ! FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, NumCycles, MaxEnergyEval
    ! (per step; <= 0 -> Glob_MaxFuncEvalForCyclOpt), OverlapThreshold,
    ! LinCoeffThreshold, SavingFreq (save every SavingFreq steps; < 1 -> 1).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: K                   ! NOT REFERENCED - see the header
    INTEGER, INTENT(IN)  :: FuncBegin, FuncEnd  ! range of functions to sweep
    INTEGER, INTENT(IN)  :: NumOfFuncToOpt      ! window size
    INTEGER, INTENT(IN)  :: NumOfFuncToShift    ! how far the window advances per step
    INTEGER, INTENT(IN)  :: NumCycles           ! sweeps over the range
    INTEGER, INTENT(IN)  :: MaxEnergyEval       ! evaluations per step, <=0 = use default
    REAL(wp), INTENT(IN) :: OverlapThreshold    ! pair-overlap rejection, <=0 = off
    REAL(wp), INTENT(IN) :: LinCoeffThreshold   ! linear-coefficient rejection, <=0 = off
    INTEGER, INTENT(IN)  :: SavingFreq          ! save every SavingFreq steps, <1 = 1

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and loop indices -----------------------------------
    INTEGER :: i, j  ! loop counters, and the range bounds
    ! handed to the permutation routines
    INTEGER :: ii              ! counts the violations an acceptance test finds
    INTEGER :: m               ! upper bound passed to PermuteFunctions2
    INTEGER :: CurrCycle       ! sweep number, 1..NumCycles
    INTEGER :: CurrFunc        ! first function of the current window
    INTEGER :: CurrFuncBegin   ! where this cycle's sweep starts
    INTEGER :: OptIterCounter  ! step number within the current cycle
    INTEGER :: totsteps        ! steps done in this call; drives the saving

    ! -- basis and window geometry -----------------------------------
    INTEGER :: cbs     ! Glob_CurrBasisSize, the basis worked on
    INTEGER :: npt     ! Glob_npt, nonlinear parameters per function
    INTEGER :: nfo     ! functions in the current window
    INTEGER :: nfru    ! frozen functions ahead of the window, cbs-nfo
    INTEGER :: nfrup1  ! nfru+1, first function of the window
    INTEGER :: nfs     ! min(NumOfFuncToShift,nfo)
    INTEGER :: nv      ! nfo*npt, optimization variables this step
    INTEGER :: nvmax   ! NumOfFuncToOpt*npt, the largest nv can be
    INTEGER :: fbn     ! start of the swept range once moved to the end
    INTEGER :: ip      ! functions of the range a previous run finished
    INTEGER :: q       ! functions left in the range at the cycle start
    INTEGER :: nr      ! rows to permute this step (ruler sequence)

    ! -- energies ----------------------------------------------------
    REAL(wp) :: Evalue  ! value returned by the latest solve
    REAL(wp) :: E_best  ! lowest energy this step has reached
    REAL(wp) :: E_prev  ! energy before the step, to fall back on
    REAL(wp) :: t       ! scale factor for the DRMNG vector D

    ! -- status flags ------------------------------------------------
    LOGICAL :: IsSwapFileOK      ! .true. when H and S came from the swap file
    LOGICAL :: ExitNeeded        ! ends the reverse-communication loop
    LOGICAL :: LastIter          ! this is the final step of the cycle
    LOGICAL :: IsOverlapBad      ! an overlap exceeded OverlapThreshold
    LOGICAL :: IsAnyLinCoeffBad  ! a coefficient exceeded LinCoeffThreshold

    ! -- solver bookkeeping ------------------------------------------
    INTEGER :: ErrCode          ! non-zero when a solve failed
    INTEGER :: NumOfFailures    ! failed solves in this step
    INTEGER :: NumOfEnergyEval  ! energy evaluations in this step
    INTEGER :: NumOfGradEval    ! gradient evaluations in this step
    ! ILAENV block size; sizes the LAPACK workspace for DSYGVX
    INTEGER :: BlockSizeForDSYGVX

    ! -- effective limits --------------------------------------------
    !    Derived from the arguments above, so that an absent or nonsensical
    !    input value cannot disable the optimization or divide by zero.
    INTEGER :: MaxEvalToUse  ! energy evaluations allowed per step
    INTEGER :: SaveEvery     ! save every SaveEvery steps

    ! -- declared but NOT REFERENCED ---------------------------------
    !    np and nfco are assigned once and never read; nfa and
    !    AttemptToGetGoodOverlap are never touched at all. The last one is
    !    a leftover from BasisEnlG's retry loop, which this routine has not.
    INTEGER :: np
    INTEGER :: nfa
    INTEGER :: nfco
    INTEGER :: AttemptToGetGoodOverlap

    ! -- workspace for the permutation and sorting routines ----------
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: NonlinParamTemp
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: FuncNumTemp
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TempR

    ! -- the optimization variables ----------------------------------
    !    The nonlinear parameters of the window, flattened. x_init is the
    !    point the step started from, to restore if the step is rejected.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_init, x_best, grad

    ! -- arrays and settings used by DRMNG ---------------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D                      ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V, V_init              ! work array and its copy
    INTEGER, PARAMETER                  :: LIV = 60               ! length of IV
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV                     ! length of V
    INTEGER                             :: ALG                    ! 2 = unconstrained minimization


    !==================================================================
    ! Global state this routine works under
    !==================================================================
    ! The basis size is Glob_CurrBasisSize (not K). fbn is where the swept
    ! range starts once moved to the end: FuncBegin shifted down by the
    ! length of the tail FuncEnd..cbs. Overlap penalties are off: a
    ! violation REJECTS the step.
    !------------------------------------------------------------------
    cbs = Glob_CurrBasisSize
    Glob_GSEPSolutionMethod = 'G'
    Glob_OverlapPenaltyAllowed = .FALSE.
    np = Glob_np
    npt = Glob_npt
    nvmax = NumOfFuncToOpt*npt
    Glob_HSLeadDim = Glob_CurrBasisSize
    Glob_HSBuffLen = Glob_CurrBasisSize*NumOfFuncToOpt
    Glob_nfa = Glob_CurrBasisSize
    nfco = FuncEnd-FuncBegin+1
    fbn = FuncBegin+Glob_CurrBasisSize-FuncEnd


    !==================================================================
    ! Nothing to do?
    !==================================================================
    ! Both conditions have to hold: the required number of cycles is
    ! done AND the last of them reached the end of the range.
    !
    ! Printed on rank 0 only. Without the guard every process wrote
    ! these six lines and the output came out interleaved.
    !------------------------------------------------------------------
    ! Checking if cyclic optimization is already completed for this basis size
    IF ((Glob_History(cbs)%CyclesDone >= NumCycles) .AND. &
        (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleG started'
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
        IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
          FuncBegin, '-', FuncEnd, 'is already completed'
        IF (Verbose >= 1) WRITE(*, *) 'Exiting OptCycleG...'
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleG finished'
      ENDIF
      RETURN
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleG started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
      WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin, '-', FuncEnd, 'will be performed'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval', MaxEnergyEval
    ENDIF


    !==================================================================
    ! Turn the two input limits into usable values
    !==================================================================
    ! MaxEnergyEval (field G of OPT_CYCLE) and SavingFreq (field H) are not
    ! validated by ReadIOFile. MaxEnergyEval <= 0 would make the DRMNG loop
    ! exit before its first iteration and leave every function untouched,
    ! so it falls back on Glob_MaxFuncEvalForCyclOpt; SavingFreq <= 0
    ! would reach a MOD by zero.
    !------------------------------------------------------------------
    IF (MaxEnergyEval > 0) THEN
      MaxEvalToUse = MaxEnergyEval
    ELSE
      MaxEvalToUse = Glob_MaxFuncEvalForCyclOpt
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a,1x,i0)') &
          'Warning WC0131 in OptCycleG: MaxEnergyEval is', MaxEnergyEval, ',', &
          'using the default limit of', Glob_MaxFuncEvalForCyclOpt
      ENDIF
    ENDIF

    SaveEvery = MAX(SavingFreq, 1)
    IF ((SavingFreq < 1) .AND. (Glob_ProcID == 0)) THEN
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a)') &
        'Warning WC0132 in OptCycleG: SavingFreq is', SavingFreq, ',', &
        'results will be saved after every step'
    ENDIF


    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! Everything is sized by cbs, the basis this routine works on - it
    ! does not change here.
    !------------------------------------------------------------------
    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_D(2*npt, NumOfFuncToOpt, cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    !==================================================================
    ! LAPACK workspace for DSYGVX
    !==================================================================
    ! (BlockSize+3)*cbs is the OPTIMAL size ILAENV suggests; 8*cbs is
    ! the minimum DSYGVX requires. The larger of the two keeps the call
    ! valid whichever library ILAENV comes from.
    !------------------------------------------------------------------
    ! Allocate workspace for DSYGVX
    BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
    Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
    ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))

    ! Allocate workspace for EnergyGB
    ALLOCATE(Glob_WkGR(NumOfFuncToOpt*npt))

    ! Allocate workspace for SaveResults (we will use Sort='yes'
    ! option, which requires workspace)
    ALLOCATE(Glob_IntWorkArrForSaveResults(cbs))

    ! Allocate arrays used by DRMNG
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))

    !------------------------------------------------------------------
    ! Local workspace
    !------------------------------------------------------------------
    ! The three Temp arrays are the scratch the permutation and sorting
    ! routines need; cbs-FuncBegin+1 is the longest range any of them is
    ! ever asked to handle.
    !------------------------------------------------------------------
    ! Allocate workspace
    ALLOCATE(x(nvmax))
    ALLOCATE(x_init(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(NonlinParamTemp(npt, cbs-FuncBegin+1))
    ALLOCATE(FuncNumTemp(cbs-FuncBegin+1))
    ALLOCATE(TempR(cbs-FuncBegin+1))


    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next - 1 for
    ! an energy, 2 for a gradient - and we supply it and call again.
    ! That is what makes it usable here, where evaluating the energy is
    ! itself a collective operation over all the processes.
    !------------------------------------------------------------------
    ! Setting some parameters for DRMNG
    ! We do it outside of the main loop so that no time
    ! is wasted for doing exactly the same operation over
    ! and over again

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): iteration and function-evaluation limits, set out
    ! of the way because the budget is enforced by MaxEvalToUse below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    ! V(35)=0.1*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed


    !==================================================================
    ! Bring the basis into the layout the sweep expects
    !==================================================================
    ! Three moves: REVERSE FuncBegin..FuncEnd (so the function the sweep
    ! starts on ends up last), move the range to the END of the basis, and
    ! rotate out any part a previous run already finished. The matrix
    ! elements are permuted along only when the swap file supplied usable
    ! ones; otherwise they are recomputed below anyway.
    !------------------------------------------------------------------
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    ! Changing the order of the functions to be optimized and, if necessary,
    ! the corresponding matrix elements to reverse
    CALL ReverseFuncOrder(FuncBegin, FuncEnd)
    IF (IsSwapFileOK) CALL ReverseMatElemOrder(FuncBegin, FuncEnd)

    ! Shifting the set of basis functions to be optimized to
    ! the very end and, if necessary, doing the permutation of matrix elements
    ! to reflect this change.
    CALL PermuteFunctions(FuncBegin, FuncEnd, FuncNumTemp, NonlinParamTemp)
    IF (IsSwapFileOK) CALL PermuteMatrixElements(FuncBegin, FuncEnd, TempR)

    !------------------------------------------------------------------
    ! Normalize the restart position
    !------------------------------------------------------------------
    ! A recorded position outside FuncBegin..FuncEnd-1 means this range
    ! is not where the previous run stopped, so the sweep starts from
    ! FuncBegin. Landing at or past FuncEnd means that cycle finished,
    ! so it is counted and the next one starts over.
    !------------------------------------------------------------------
    ! If the last optimized function number is greater than FuncEnd-1 or smaller than FuncBegin
    ! we need to change it so that the optimization begins from FuncBegin
    IF (Glob_History(cbs)%InitFuncAtLastStep < FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
    IF (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd) THEN
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
    ENDIF

    ! We need to make an initial permutation to shift the functions that were
    ! already optimized (if any) in the current optimization cycle
    ip = Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
    IF (ip > 0) THEN
      CALL PermuteFunctions(fbn, cbs-ip, FuncNumTemp, NonlinParamTemp)
      IF (IsSwapFileOK) CALL PermuteMatrixElements(fbn, cbs-ip, TempR)
    ENDIF


    !==================================================================
    ! Initial energy
    !==================================================================
    ! With usable matrix elements only the eigenvalue problem has to be
    ! solved; without them everything is recomputed, which is what the
    ! .true. asks for.
    !------------------------------------------------------------------
    ! Calculating the initial energy
    IF (IsSwapFileOK) THEN
      ! Getting initial energy
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyGA(1, cbs, .FALSE., ErrCode)
    ELSE
      ! Getting initial energy
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyGA(1, cbs, .TRUE., ErrCode)
    ENDIF
    IF (ErrCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0145 in OptCycleG: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy


    !==================================================================
    ! OUTER LOOP - one pass over FuncBegin..FuncEnd per cycle
    !==================================================================
    ! Starts at CyclesDone+1, so a resumed run does not repeat cycles it
    ! already finished.
    !------------------------------------------------------------------
    ! Here comes main optimization cycle
    totsteps = 0

    DO CurrCycle = Glob_History(cbs)%CyclesDone+1, NumCycles

      ! Doing cycle number CurrCycle
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Cycle', CurrCycle, 'began'
      ENDIF

      CurrFuncBegin = Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      q = FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
      OptIterCounter = 0


      !==================================================================
      ! INNER LOOP - one window per step
      !==================================================================
      DO CurrFunc = CurrFuncBegin, FuncEnd, NumOfFuncToShift

        ! Note that CurrFunc counts functions as they were not shifted back
        ! by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
        ! initial input file.
        totsteps = totsteps+1

        ! The window is the LAST nfo functions of the basis; it is
        ! short only on the final step of a cycle, when fewer than
        ! NumOfFuncToOpt functions are left in the range.
        nfo = MIN(FuncEnd-CurrFunc+1, NumOfFuncToOpt)
        Glob_nfo = nfo
        nv = nfo*npt
        nfru = cbs-nfo
        Glob_nfru = nfru
        nfrup1 = nfru+1
        OptIterCounter = OptIterCounter+1

        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *)
          IF (nfo > 1) THEN
            IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Optimizing functions', CurrFunc, '-', CurrFunc+nfo-1
          ELSE
            WRITE(*, '(1x,a,1x,i0)') 'Optimizing function', CurrFunc
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Advance the window
        !------------------------------------------------------------------
        ! Skipped on the first step of a cycle. PermuteFunctions/
        ! PermuteMatrixElements push the block just optimized back in front of
        ! the untouched ones, then PermuteFunctions2/PermuteMatrixElements2
        ! swap two adjacent blocks of nr functions to bring the next window
        ! forward. nr follows the ruler sequence NumOfRowsToPermForUnitShift
        ! (1,2,1,4,...), spreading the cost of keeping the basis ordered; the
        ! "abnormal" branch is the tail of a cycle, where fewer than 2*nr
        ! functions remain.
        ! The permutations are applied whenever the first block is non-empty;
        ! the matrix elements are real by now.
        !------------------------------------------------------------------
        nfs = MIN(NumOfFuncToShift, nfo)

        IF (OptIterCounter /= 1) THEN

          i = cbs-nfo-nfs+1
          j = cbs-nfs
          CALL PermuteFunctions(i, j, FuncNumTemp, NonlinParamTemp)
          CALL PermuteMatrixElements(i, j, TempR)

          nr = NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)

          IF (q-nfo >= nr*2) THEN
            ! normal permutation, there is sufficient number of functions left
            m = cbs-nfo
            j = m-nr
            i = j-nr+1
          ELSE
            ! abnormal permutation, not sufficient number of functions left
            m = cbs-nfo
            j = m-nr
            i = cbs-q+1
          ENDIF

          ! At the tail of a cycle the first block can be empty (i = j+1: the ruler
          ! width equals the number of functions not yet reordered). Nothing is
          ! left to swap then; the reference code falls through as a no-op.
          IF (i <= j) THEN
            CALL PermuteFunctions2(i, j, m, FuncNumTemp, NonlinParamTemp)
            CALL PermuteMatrixElements2(i, j, m, TempR)
          ENDIF

        ENDIF

        IF (Glob_ProcID == 0) THEN
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 2) WRITE (*, *) 'Nonlinear parameters before optimization:'
            DO i = 1, nfo
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(nfru+i), ':', Glob_PWR(nfru+i)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
            ENDDO
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Set up the step
        !------------------------------------------------------------------
        ! x holds the window's parameters flattened, x_init a copy to restore
        ! from if the step is rejected. D is the DRMNG scale vector: the
        ! 1/(cbs^2*sqrt(cbs)) form shrinks the steps as the basis grows (a
        ! stiffer problem), floored at 10000*epsilon.
        !------------------------------------------------------------------
        IV(1:LIV) = IV_init(1:LIV)
        V(1:LV) = V_init(1:LV)

        DO i = 1, nfo
          x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
          x_init((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
        ENDDO

        t = MAX(ONE/(cbs*cbs*SQRT(ONE*cbs)), 10000*EPSILON(Glob_CurrEnergy))
        DO i = 1, nfo
          ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          DO j = 1, npt
            ! Make sure none of the D(i) will be zero or smaller than the threshold
            ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j) = t
            ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          ENDDO
        ENDDO

        ExitNeeded = .FALSE.
        NumOfFailures = 0
        NumOfEnergyEval = 0
        NumOfGradEval = 0
        IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.
        E_best = Glob_CurrEnergy
        x_best(1:nfo*npt) = x(1:nfo*npt)
        ! Remember the energy that corresponds to the current (i.e. unchanged) values
        ! of the nonlinear parameters. It is needed in case the optimization of this
        ! function/set of functions has to be abandoned
        E_prev = Glob_CurrEnergy

        !------------------------------------------------------------------
        ! The reverse-communication loop
        !------------------------------------------------------------------
        ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it wants:
        ! 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its evaluation
        ! limit. On a failed evaluation the energy handed to DRMNG keeps its
        ! previous value, so the step counts as no reduction and the trust
        ! radius shrinks (IV(2), the TOOBIG flag, is deliberately left 0). The
        ! best point is tracked here because the last point DRMNG visits is not
        ! necessarily the lowest.
        !------------------------------------------------------------------
        DO WHILE (.NOT. (ExitNeeded))

          IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
          CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

          SELECT CASE (IV(1))

          CASE (1)  ! Only energy is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO i = 1, nfo
              Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
            ENDDO
            Evalue = EnergyGA(nfrup1, cbs, .TRUE., ErrCode)
            NumOfEnergyEval = NumOfEnergyEval+1
            IF (ErrCode /= 0) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
            ELSE
              Glob_CurrEnergy = Evalue
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nfo*npt) = x(1:nfo*npt)
              ENDIF
            ENDIF
            ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
            IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

          CASE (2)  ! Only gradient is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO i = 1, nfo
              Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
            ENDDO
            CALL EnergyGB(Evalue, grad, .TRUE., ErrCode)
            NumOfGradEval = NumOfGradEval+1
            IF (ErrCode /= 0) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
            ELSE
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nfo*npt) = x(1:nfo*npt)
              ENDIF
            ENDIF
            ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
            IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

          CASE (3:8)  ! Some kind of convergence has been reached
            ExitNeeded = .TRUE.

          CASE (9:10)  ! Function evaluation limit has been reached.
            ! This is never supposed to happen because we
            ! count the number of function evaluations ourselves.
            ExitNeeded = .TRUE.


          CASE DEFAULT
            ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
            ! and >=14 for a bad input. None of those match a case above,
            ! so without this the loop would call DRMNG again for ever.
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                'Warning WC0139 in OptCycleG: DRMNG returned IV(1) =', IV(1)
              IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
            ENDIF
            ExitNeeded = .TRUE.
          ENDSELECT

          ! A warning, not an abort: the best point found so far is
          ! still usable, and the acceptance tests below decide what
          ! becomes of this step.
          IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0123 in OptCycleG: number of failures in energy or gradient', &
                'calculations during the optimization of nonlinear parameters', &
                'reached the limit of', Glob_MaxEnergyFailsAllowed
            ENDIF
            ! call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
          ENDIF

          IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.

        ENDDO  ! while


        !------------------------------------------------------------------
        ! Re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! EnergyGAM rather than EnergyGA because it also produces Glob_c,
        ! which the linear-coefficient test needs.
        !------------------------------------------------------------------
        ! Compute the energy and the linear coefficients at the best point found
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, nfru+i) = x_best((i-1)*npt+1:i*npt)
        ENDDO

        Evalue = EnergyGAM(nfrup1, cbs, .TRUE., ErrCode)
        !!We run EnergyGA again because EnergyGAM might give slightly different
        !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
        ! Glob_CurrEnergy=EnergyGA(nfrup1,cbs,.false.,ErrCode)
        IF (ErrCode == 0) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          ! Note that Glob_CurrEnergy is intentionally left unchanged here as
          ! the value returned by a failed energy evaluation is meaningless
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0120 in OptCycleG: failed to evaluate energy after optimization', &
              'of nonlinear parameters. The values of the nonlinear parameters', &
              'will be left unchanged.'
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Acceptance test: pair overlaps
        !------------------------------------------------------------------
        ! Every pair involving a function in the window, against the whole
        ! basis. Only the LOWER triangle of Glob_S is read. Function NUMBERS
        ! are printed, not positions, because the basis has been permuted
        ! and the positions would mean nothing to the reader.
        !------------------------------------------------------------------
        ! Checking if overlap is OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad = .FALSE.
        IF (OverlapThreshold > ZERO) THEN
          ii = 0
          DO i = nfrup1, cbs
            DO j = 1, i-1
              IF (ABS(Glob_S(i, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0121: overlap of the following functions exceeds threshold. ', &
                      'Nonlinear parameters will be left unchanged'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') &
                    ii, ':', Glob_FuncNum(i), Glob_FuncNum(j), '    S='
                  CALL writerealadv(6, Glob_S(i, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Acceptance test: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..cbs, not just the window:
        ! moving one function can blow up the coefficient of another, and
        ! that is exactly the near-linear-dependence this test is for.
        !------------------------------------------------------------------
        ! Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad = .FALSE.
        IF (LinCoeffThreshold > ZERO) THEN
          ii = 0
          DO i = 1, cbs
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0122: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', i, '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        ! Reported only for a step that was actually kept
        IF ((Glob_ProcID == 0) .AND. (ErrCode == 0) .AND. (.NOT. IsOverlapBad) .AND. (.NOT. IsAnyLinCoeffBad)) THEN
          IF (Verbose >= 1) THEN
          WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
          ENDIF
          WRITE (*, *) 'E=', Glob_CurrEnergy
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 1) WRITE (*, *) 'Nonlinear parameters after optimization:'
            DO i = 1, nfo
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(nfru+i), ':', Glob_PWR(nfru+i)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
            ENDDO
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Undo the step if it was not acceptable
        !------------------------------------------------------------------
        ! The parameters go back to x_init and the matrix elements are
        ! recomputed for them, so H and S match the basis again before the
        ! window moves on. If even THAT re-solve fails the parameters are
        ! still the original ones and their matrix elements have just been
        ! computed, so it is safe to carry on with the recorded energy.
        !------------------------------------------------------------------
        IF ((ErrCode /= 0) .OR. IsOverlapBad .OR. IsAnyLinCoeffBad) THEN
          ! restore initial values of nonlinear parameters
          DO i = 1, nfo
            Glob_NonlinParam(1:npt, nfru+i) = x_init((i-1)*npt+1:i*npt)
          ENDDO
          Evalue = EnergyGA(nfrup1, cbs, .TRUE., ErrCode)
          IF (ErrCode == 0) THEN
            Glob_CurrEnergy = Evalue
          ELSE
            ! Leave this function (or set of functions) as it was and proceed to
            ! the next one. The nonlinear parameters are the original ones and
            ! the matrix elements that correspond to them have been computed,
            ! so it is safe to continue
            IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE(*, '(1x,a,1x,a)') &
              'Warning WC0124 in OptCycleG: energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_CurrEnergy = E_prev
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Record the step and save
        !------------------------------------------------------------------
        ! On the last step of a cycle the position is reset to 0 and the
        ! cycle counter advances, which is what tells a resumed run that
        ! this cycle is finished.
        !
        ! Saving happens on the first few steps whatever SaveEvery says, so
        ! that a run interrupted early still leaves something behind, then
        ! every SaveEvery steps, and always on the final step of a cycle.
        !------------------------------------------------------------------
        IF (CurrFunc > FuncEnd-NumOfFuncToShift) THEN
          LastIter = .TRUE.
        ELSE
          LastIter = .FALSE.
        ENDIF

        Glob_History(cbs)%Energy = Glob_CurrEnergy
        IF (LastIter) THEN
          Glob_History(cbs)%InitFuncAtLastStep = 0
          Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
        ELSE
          Glob_History(cbs)%InitFuncAtLastStep = CurrFunc
        ENDIF

        IF (Glob_ProcID == 0) THEN
          IF ((totsteps <= Glob_MinMandSavSteps) .OR. (MOD(totsteps, SaveEvery) == 0) .OR. &
              (CurrFunc+NumOfFuncToShift >= FuncEnd)) THEN
            CALL SaveResults(Sort='yes')
          ENDIF
        ENDIF

      ENDDO  ! end cycle CurrCycle


      !------------------------------------------------------------------
      ! End of a cycle
      !------------------------------------------------------------------
      ! Between cycles the basis is sorted back into function-number order
      ! so the next sweep starts from a known layout. Not done after the
      ! LAST cycle, because the final ordering below covers it.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Cycle', CurrCycle, ' finished'
      ENDIF

      IF (CurrCycle /= NumCycles) THEN
        Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
        IF (Glob_ProcID == 0) WRITE(*, '(1x,a47)', ADVANCE='no') &
          'Ordering basis functions and matrix elements...'
        CALL SortBasisFuncAndMatElem(fbn, cbs, FuncNumTemp, NonlinParamTemp, TempR)
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
      ENDIF

    ENDDO  ! End of main optimization cycle


    !==================================================================
    ! Put the basis back the way it came in
    !==================================================================
    ! Sort restores function-number order and the reverse undoes the
    ! reversal applied at the start, so the basis leaves this routine in
    ! the same layout it arrived in - only with better parameters. Both
    ! run on every rank, since every rank holds its own copy.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) WRITE(*, '(1x,a53)', ADVANCE='no') &
      'Final ordering basis functions and matrix elements...'
    CALL SortBasisFuncAndMatElem(FuncBegin, cbs, FuncNumTemp, NonlinParamTemp, TempR)
    CALL ReverseFuncOrder(FuncBegin, cbs)
    CALL ReverseMatElemOrder(FuncBegin, cbs)
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'

    ! Hand H and S to the next BBOP step
    CALL StoreMatricesInSwapFile()


    !==================================================================
    ! Release everything, in the reverse of the allocation order
    !==================================================================
    ! deallocate workspace
    DEALLOCATE(TempR)
    DEALLOCATE(FuncNumTemp)
    DEALLOCATE(NonlinParamTemp)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x_init)
    DEALLOCATE(x)

    ! deallocate arrays used by DRMNG
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)

    ! deallocate workspace for SaveResults
    DEALLOCATE(Glob_IntWorkArrForSaveResults)

    ! deallocate workspace for EnergyGB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for DSYGVX
    DEALLOCATE(Glob_IWorkForDSYGVX)
    DEALLOCATE(Glob_WorkForDSYGVX)

    ! Deallocate some global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine OptCycleG finished'


  END SUBROUTINE OptCycleG


  SUBROUTINE OptCycleI(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                       NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold, SavingFreq)
    !==================================================================
    ! Subroutine OptCycleI
    !==================================================================
    ! Improves an EXISTING basis by sweeping a window of NumOfFuncToOpt
    ! functions over FuncBegin..FuncEnd with INVERSE ITERATION ('I'); a
    ! twin of OptCycleG (sweep, window schedule, acceptance tests and
    ! restart bookkeeping are documented there). Differences: Glob_H holds
    ! the SHIFTED matrix H - Glob_ApproxEnergy*S and Glob_invD the LDL'
    ! diagonal; the workspace is Glob_WorkForGSEPIIS and Glob_LastEigvector;
    ! EnergyIA/EnergyIAM/EnergyIB replace EnergyGA/EnergyGAM/EnergyGB;
    ! v_good keeps the last converged eigenvector; two REFACTORIZATION
    ! solves (EnergyIA with AreMatElemNeeded=.false., fatal on failure:
    ! EC0151, EC0154) rebuild the factorization after the window advances
    ! and after the between-cycle sort, since any permutation invalidates
    ! it; one extra retry on the rejection path; the average number of
    ! inverse iterations per solve is reported. The shift is refreshed per
    ! window (RefreshINVITShift) and, with Glob_EigIdxTargeting=1, placed
    ! on eigenvalue Glob_WhichEigenvalue at the start of the step
    ! (RetargetShiftToEigenvalue). Arguments as OptCycleG; K is NOT
    ! REFERENCED.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: K                   ! NOT REFERENCED - see the header
    INTEGER, INTENT(IN)  :: FuncBegin, FuncEnd  ! range of functions to sweep
    INTEGER, INTENT(IN)  :: NumOfFuncToOpt      ! window size
    INTEGER, INTENT(IN)  :: NumOfFuncToShift    ! how far the window advances per step
    INTEGER, INTENT(IN)  :: NumCycles           ! sweeps over the range
    INTEGER, INTENT(IN)  :: MaxEnergyEval       ! evaluations per step, <=0 = use default
    REAL(wp), INTENT(IN) :: OverlapThreshold    ! pair-overlap rejection, <=0 = off
    REAL(wp), INTENT(IN) :: LinCoeffThreshold   ! linear-coefficient rejection, <=0 = off
    INTEGER, INTENT(IN)  :: SavingFreq          ! save every SavingFreq steps, <1 = 1

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and loop indices -----------------------------------
    INTEGER :: i, j  ! loop counters, and the range bounds
    ! handed to the permutation routines
    INTEGER :: ii              ! counts the violations an acceptance test finds
    INTEGER :: m               ! upper bound passed to PermuteFunctions2
    INTEGER :: CurrCycle       ! sweep number, 1..NumCycles
    INTEGER :: CurrFunc        ! first function of the current window
    INTEGER :: CurrFuncBegin   ! where this cycle's sweep starts
    INTEGER :: OptIterCounter  ! step number within the current cycle
    INTEGER :: totsteps        ! steps done in this call; drives the saving

    ! -- basis and window geometry -----------------------------------
    INTEGER :: cbs     ! Glob_CurrBasisSize, the basis worked on
    INTEGER :: npt     ! Glob_npt, nonlinear parameters per function
    INTEGER :: nfo     ! functions in the current window
    INTEGER :: nfru    ! frozen functions ahead of the window, cbs-nfo
    INTEGER :: nfrup1  ! nfru+1, first function of the window
    INTEGER :: nfs     ! min(NumOfFuncToShift,nfo)
    INTEGER :: nv      ! nfo*npt, optimization variables this step
    INTEGER :: nvmax   ! NumOfFuncToOpt*npt, the largest nv can be
    INTEGER :: fbn     ! start of the swept range once moved to the end
    INTEGER :: ip      ! functions of the range a previous run finished
    INTEGER :: q       ! functions left in the range at the cycle start
    INTEGER :: nr      ! rows to permute this step (ruler sequence)

    ! -- energies ----------------------------------------------------
    REAL(wp) :: Evalue  ! value returned by the latest solve
    REAL(wp) :: E_best  ! lowest energy this step has reached
    REAL(wp) :: E_prev  ! energy before the step, to fall back on
    REAL(wp) :: t       ! scale factor for the DRMNG vector D

    ! -- status flags ------------------------------------------------
    LOGICAL :: IsSwapFileOK      ! .true. when H and S came from the swap file
    LOGICAL :: ExitNeeded        ! ends the reverse-communication loop
    LOGICAL :: LastIter          ! this is the final step of the cycle
    LOGICAL :: IsOverlapBad      ! an overlap exceeded OverlapThreshold
    LOGICAL :: IsAnyLinCoeffBad  ! a coefficient exceeded LinCoeffThreshold

    ! -- solver bookkeeping ------------------------------------------
    INTEGER :: ErrCode          ! non-zero when a solve failed
    INTEGER :: NumOfFailures    ! failed solves in this step
    INTEGER :: NumOfEnergyEval  ! energy evaluations in this step
    INTEGER :: NumOfGradEval    ! gradient evaluations in this step

    ! -- effective limits --------------------------------------------
    !    Derived from the arguments above, so that an absent or nonsensical
    !    input value cannot disable the optimization or divide by zero.
    INTEGER :: MaxEvalToUse  ! energy evaluations allowed per step
    INTEGER :: SaveEvery     ! save every SaveEvery steps

    ! -- declared but NOT REFERENCED ---------------------------------
    !    np and nfco are assigned once and never read; nfa and
    !    AttemptToGetGoodOverlap are never touched at all. The last one is
    !    a leftover from BasisEnlI's retry loop, which this routine has not.
    INTEGER :: np
    INTEGER :: nfa
    INTEGER :: nfco
    INTEGER :: AttemptToGetGoodOverlap

    ! -- workspace for the permutation and sorting routines ----------
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: NonlinParamTemp
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: FuncNumTemp
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TempR

    ! -- the optimization variables ----------------------------------
    !    The nonlinear parameters of the window, flattened. x_init is the
    !    point the step started from, to restore if the step is rejected.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, x_init, x_best, grad

    ! -- inverse-iteration fallback ----------------------------------
    !    Last eigenvector that came out of a SUCCESSFUL solve - see the
    !    note where it is allocated.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: v_good

    ! -- arrays and settings used by DRMNG ---------------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D                      ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V, V_init              ! work array and its copy
    INTEGER, PARAMETER                  :: LIV = 60               ! length of IV
    INTEGER                             :: IV(LIV), IV_init(LIV)
    INTEGER                             :: LV                     ! length of V
    INTEGER                             :: ALG                    ! 2 = unconstrained minimization


    !==================================================================
    ! Global state this routine works under
    !==================================================================
    ! The basis size is taken from Glob_CurrBasisSize, not from the K
    ! argument. fbn is where the swept range starts once it has been
    ! moved to the end.
    !
    ! Overlap penalties are off: an overlap violation REJECTS the step
    ! here rather than being folded into the energy.
    !------------------------------------------------------------------
    ! Setting the values of some global variables
    cbs = Glob_CurrBasisSize
    Glob_GSEPSolutionMethod = 'I'
    Glob_OverlapPenaltyAllowed = .FALSE.
    np = Glob_np
    npt = Glob_npt
    nvmax = NumOfFuncToOpt*npt
    Glob_HSLeadDim = Glob_CurrBasisSize
    Glob_HSBuffLen = Glob_CurrBasisSize*NumOfFuncToOpt
    Glob_nfa = Glob_CurrBasisSize
    nfco = FuncEnd-FuncBegin+1
    fbn = FuncBegin+Glob_CurrBasisSize-FuncEnd


    !==================================================================
    ! Nothing to do?
    !==================================================================
    ! Both conditions have to hold: the required number of cycles is
    ! done AND the last of them reached the end of the range.
    !
    ! Printed on rank 0 only. Without the guard every process wrote
    ! these six lines and the output came out interleaved.
    !------------------------------------------------------------------
    ! Checking if cyclic optimization is already completed for this basis size
    IF ((Glob_History(cbs)%CyclesDone >= NumCycles) .AND. &
        (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleI started'
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
        IF (Verbose >= 2) WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
          FuncBegin, '-', FuncEnd, 'is already completed'
        IF (Verbose >= 1) WRITE(*, *) 'Exiting OptCycleI...'
        IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleI finished'
      ENDIF
      RETURN
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine OptCycleI started'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'Basis size is', cbs
      WRITE(*, '(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin, '-', FuncEnd, 'will be performed'
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') 'MaxEnergyEval', MaxEnergyEval
    ENDIF


    !==================================================================
    ! Turn the two input limits into usable values
    !==================================================================
    ! MaxEnergyEval (field G of OPT_CYCLE) and SavingFreq (field H) are not
    ! validated by ReadIOFile. MaxEnergyEval <= 0 would make the DRMNG loop
    ! exit before its first iteration and leave every function untouched,
    ! so it falls back on Glob_MaxFuncEvalForCyclOpt; SavingFreq <= 0
    ! would reach a MOD by zero.
    !------------------------------------------------------------------
    IF (MaxEnergyEval > 0) THEN
      MaxEvalToUse = MaxEnergyEval
    ELSE
      MaxEvalToUse = Glob_MaxFuncEvalForCyclOpt
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a,1x,i0)') &
          'Warning WC0133 in OptCycleI: MaxEnergyEval is', MaxEnergyEval, ',', &
          'using the default limit of', Glob_MaxFuncEvalForCyclOpt
      ENDIF
    ENDIF

    SaveEvery = MAX(SavingFreq, 1)
    IF ((SavingFreq < 1) .AND. (Glob_ProcID == 0)) THEN
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,1x,a)') &
        'Warning WC0134 in OptCycleI: SavingFreq is', SavingFreq, ',', &
        'results will be saved after every step'
    ENDIF


    !==================================================================
    ! Allocate the matrices, the derivative store and the MPI buffers
    !==================================================================
    ! No Glob_diagH here: on the 'I' path the diagonal of H stays inside
    ! Glob_H. Glob_invD takes its place, holding the diagonal of the
    ! LDL' factorization that inverse iteration works with.
    !------------------------------------------------------------------
    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_invD(cbs))
    ALLOCATE(Glob_D(2*npt, NumOfFuncToOpt, cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    !==================================================================
    ! Inverse-iteration workspace
    !==================================================================
    ! Glob_LastEigvector is the STARTING VECTOR of the next inverse
    ! iteration. Carrying the previous solution over is what makes the
    ! method cheap: consecutive evaluations differ only slightly, so a
    ! few iterations suffice. It is primed with all ones because there
    ! is nothing better to start from.
    !------------------------------------------------------------------
    ! Allocate workspace for subroutine GSEPIIS, which is called
    ! inside EnergyIA, EnergyIAM, and EnergyIB
    ALLOCATE(Glob_WorkForGSEPIIS(cbs))
    ALLOCATE(Glob_LastEigvector(cbs))
    Glob_LastEigvector(1:cbs) = ONE

    !------------------------------------------------------------------
    ! ... and the fallback copy of it
    !------------------------------------------------------------------
    ! A failed (diverged or half-converged) inverse iteration leaves
    ! Glob_LastEigvector in a state that poisons every later solve, so the
    ! vector of the last SUCCESSFUL solve is kept in v_good and restored
    ! whenever ErrCode comes back nonzero. No counterpart on the 'G' path,
    ! where DSYGVX starts from scratch every time.
    ! Primed with ONE to match Glob_LastEigvector.
    !------------------------------------------------------------------
    ALLOCATE(v_good(cbs))
    v_good(1:cbs) = ONE

    ! Allocate workspace for EnergyIB
    ALLOCATE(Glob_WkGR(NumOfFuncToOpt*npt))

    ! Allocate workspace for SaveResults (we will use Sort='yes'
    ! option, which requires workspace)
    ALLOCATE(Glob_IntWorkArrForSaveResults(cbs))

    ! Allocate arrays used by DRMNG
    ALLOCATE(D(nvmax))
    LV = 71+nvmax*(nvmax+13)/2 + 1
    ALLOCATE(V(LV))
    ALLOCATE(V_init(LV))

    !------------------------------------------------------------------
    ! Local workspace
    !------------------------------------------------------------------
    ! The three Temp arrays are the scratch the permutation and sorting
    ! routines need; cbs-FuncBegin+1 is the longest range any of them is
    ! ever asked to handle.
    !------------------------------------------------------------------
    ! Allocate workspace
    ALLOCATE(x(nvmax))
    ALLOCATE(x_init(nvmax))
    ALLOCATE(x_best(nvmax))
    ALLOCATE(grad(nvmax))
    ALLOCATE(NonlinParamTemp(npt, cbs-FuncBegin+1))
    ALLOCATE(FuncNumTemp(cbs-FuncBegin+1))
    ALLOCATE(TempR(cbs-FuncBegin+1))


    !==================================================================
    ! Set up DRMNG
    !==================================================================
    ! DRMNG is the REVERSE-COMMUNICATION form of the SUMSL quasi-Newton
    ! minimizer: it returns with IV(1) saying what it wants next - 1 for
    ! an energy, 2 for a gradient - and we supply it and call again.
    ! That is what makes it usable here, where evaluating the energy is
    ! itself a collective operation over all the processes.
    !------------------------------------------------------------------
    ! Setting some parameters for DRMNG
    ! We do it outside of the main loop so that no time
    ! is wasted for doing exactly the same operation over
    ! and over again

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    CALL DIVSET(ALG, IV_init, LIV, LV, V_init)
    ! IV(17)/IV(18): iteration and function-evaluation limits, set out
    ! of the way because the budget is enforced by MaxEvalToUse below.
    IV_init(17) = 1000000
    IV_init(18) = 1000000
    IV_init(19) = 0  ! set summary print format
    ! Silence every report SUMSL would print by itself
    IV_init(20) = 0; IV_init(22) = 0; IV_init(23) = -1; IV_init(24) = 0
    V_init(31) = 0.0_wp
    V_init(32) = 2*EPSILON(V_init(32))
    V_init(37) = 2*EPSILON(V_init(37))
    ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
    ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
    ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35) = Glob_MaxScStepAllowedInOpt*ONE
    ! V(35)=0.1*ONE
    IV_init(1) = 12  ! DIVSET has been called and some default values were changed


    !==================================================================
    ! Bring the basis into the layout the sweep expects
    !==================================================================
    ! Three moves: reverse FuncBegin..FuncEnd so the function the sweep
    ! starts on ends up last, move that whole range to the END of the
    ! basis, and rotate out any part of it that a previous run already
    ! finished.
    !
    ! The matrix elements are only permuted along when the swap file
    ! gave us usable ones. Without it H and S hold nothing worth
    ! preserving and are recomputed below.
    !------------------------------------------------------------------
    ! Read swap file (if necessary) and distribute the data
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    ! Changing the order of the functions to be optimized and, if necessary,
    ! the corresponding matrix elements to reverse
    CALL ReverseFuncOrder(FuncBegin, FuncEnd)
    IF (IsSwapFileOK) CALL ReverseMatElemOrder(FuncBegin, FuncEnd)

    ! Shifting the set of basis functions to be optimized to
    ! the very end and, if necessary, doing the permutation of matrix elements
    ! to reflect this change.
    CALL PermuteFunctions(FuncBegin, FuncEnd, FuncNumTemp, NonlinParamTemp)
    IF (IsSwapFileOK) CALL PermuteMatrixElements(FuncBegin, FuncEnd, TempR)

    !------------------------------------------------------------------
    ! Normalize the restart position
    !------------------------------------------------------------------
    ! A recorded position outside FuncBegin..FuncEnd-1 means this range
    ! is not where the previous run stopped, so the sweep starts from
    ! FuncBegin. Landing at or past FuncEnd means that cycle finished,
    ! so it is counted and the next one starts over.
    !------------------------------------------------------------------
    ! If the last optimized function number is greater than FuncEnd-1 or smaller than FuncBegin
    ! we need to change it so that the optimization begins from FuncBegin
    IF (Glob_History(cbs)%InitFuncAtLastStep < FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
    IF (Glob_History(cbs)%InitFuncAtLastStep >= FuncEnd) THEN
      Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
    ENDIF

    ! We need to make an initial permutation to shift the functions that were
    ! already optimized (if any) in the current optimization cycle
    ip = Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
    IF (ip > 0) THEN
      CALL PermuteFunctions(fbn, cbs-ip, FuncNumTemp, NonlinParamTemp)
      IF (IsSwapFileOK) CALL PermuteMatrixElements(fbn, cbs-ip, TempR)
    ENDIF


    !==================================================================
    ! Initial energy
    !==================================================================
    ! With usable matrix elements only the eigenvalue problem has to be
    ! solved; without them everything is recomputed, which is what the
    ! .true. asks for.
    !------------------------------------------------------------------
    WrkP_WrongStateCount = 0

    ! Calculating the initial energy
    IF (IsSwapFileOK) THEN
      ! Getting initial energy
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyIA(1, cbs, .FALSE., ErrCode)
    ELSE
      ! Getting initial energy
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyIA(1, cbs, .TRUE., ErrCode)
    ENDIF
    !==================================================================
    ! Put the inverse-iteration shift on the requested eigenvalue
    !==================================================================
    ! Once per BBOP step; inert unless Glob_EigIdxTargeting==1. It runs
    ! BEFORE the fatal check below: a shift sitting between two
    ! eigenvalues at nearly equal distance is the usual reason the first
    ! solve does not converge, and moving the shift is precisely the
    ! cure.
    !------------------------------------------------------------------
    IF (Glob_EigIdxTargeting == 1) THEN
      CALL RetargetShiftToEigenvalue(cbs, 'OptCycleI')
      Glob_CurrEnergy = EnergyIA(1, cbs, .FALSE., ErrCode)
    ENDIF

    ! The second test catches a solve that converged on a level other
    ! than WHICH_EIGENVALUE: EnergyIA returns it with ErrCode=0 and the
    ! rejection sentinel as the energy - see IsRequestedEigenstate.
    IF ((ErrCode /= 0) .OR. (ABS(Glob_CurrEnergy) > 1.0E10_wp)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0150 in OptCycleI: initial energy cannot be computed'
        IF (ErrCode == 0) WRITE(*, *) '(inverse iteration converged on a level other than WHICH_EIGENVALUE)'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) WRITE(*, *) 'Initial energy ', Glob_CurrEnergy


    !==================================================================
    ! OUTER LOOP - one pass over FuncBegin..FuncEnd per cycle
    !==================================================================
    ! Starts at CyclesDone+1, so a resumed run does not repeat cycles it
    ! already finished.
    !------------------------------------------------------------------
    ! Here comes main optimization cycle
    totsteps = 0

    DO CurrCycle = Glob_History(cbs)%CyclesDone+1, NumCycles

      ! Doing cycle number CurrCycle
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,1x,a)') 'Cycle', CurrCycle, 'began'
      ENDIF

      CurrFuncBegin = Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      q = FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
      OptIterCounter = 0


      !==================================================================
      ! INNER LOOP - one window per step
      !==================================================================
      DO CurrFunc = CurrFuncBegin, FuncEnd, NumOfFuncToShift

        ! Note that CurrFunc counts functions as they were not shifted back
        ! by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
        ! initial input file.
        totsteps = totsteps+1

        ! The window is the LAST nfo functions of the basis; it is
        ! short only on the final step of a cycle, when fewer than
        ! NumOfFuncToOpt functions are left in the range.
        nfo = MIN(FuncEnd-CurrFunc+1, NumOfFuncToOpt)
        Glob_nfo = nfo
        nv = nfo*npt
        nfru = cbs-nfo
        Glob_nfru = nfru
        nfrup1 = nfru+1
        OptIterCounter = OptIterCounter+1

        ! Counter1 counts the solves done in this step and Counter2
        ! accumulates their inverse iterations; their ratio is the
        ! average reported below. Reset per step, so the figure tracks
        ! how the shift is holding up as the sweep proceeds.
        Glob_InvItTempCounter1 = 0
        Glob_InvItTempCounter2 = 0

        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *)
          IF (nfo > 1) THEN
            IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Optimizing functions', CurrFunc, '-', CurrFunc+nfo-1
          ELSE
            WRITE(*, '(1x,a,1x,i0)') 'Optimizing function', CurrFunc
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Advance the window
        !------------------------------------------------------------------
        ! Skipped on the first step of a cycle. PermuteFunctions/
        ! PermuteMatrixElements push the block just optimized back in front of
        ! the untouched ones, then PermuteFunctions2/PermuteMatrixElements2
        ! swap two adjacent blocks of nr functions to bring the next window
        ! forward. nr follows the ruler sequence NumOfRowsToPermForUnitShift
        ! (1,2,1,4,...), spreading the cost of keeping the basis ordered; the
        ! "abnormal" branch is the tail of a cycle, where fewer than 2*nr
        ! functions remain.
        ! The refactorization that follows (AreMatElemNeeded=.false.) rebuilds
        ! the factorization of H - sigma*S after the permutation; the matrix
        ! elements themselves are already correct.
        !------------------------------------------------------------------
        nfs = MIN(NumOfFuncToShift, nfo)

        IF (OptIterCounter /= 1) THEN

          i = cbs-nfo-nfs+1
          j = cbs-nfs
          CALL PermuteFunctions(i, j, FuncNumTemp, NonlinParamTemp)
          CALL PermuteMatrixElements(i, j, TempR)

          nr = NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)

          IF (q-nfo >= nr*2) THEN
            ! normal permutation, there is sufficient number of functions left
            m = cbs-nfo
            j = m-nr
            i = j-nr+1
          ELSE
            ! abnormal permutation, not sufficient number of functions left
            m = cbs-nfo
            j = m-nr
            i = cbs-q+1
          ENDIF

          ! At the tail of a cycle the first block can be empty (i = j+1: the ruler
          ! width equals the number of functions not yet reordered). Nothing is
          ! left to swap then; the reference code falls through as a no-op.
          IF (i <= j) THEN
            CALL PermuteFunctions2(i, j, m, FuncNumTemp, NonlinParamTemp)
            CALL PermuteMatrixElements2(i, j, m, TempR)
          ENDIF

          ! call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
          Glob_CurrEnergy = EnergyIA(i, cbs, .FALSE., ErrCode)
          IF ((ErrCode /= 0) .OR. (ABS(Glob_CurrEnergy) > 1.0E10_wp)) THEN
            IF (Glob_ProcID == 0) WRITE(*, *) &
              'Error EC0151 in OptCycleI: energy cannot be computed after permuting basis functions'
            CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
          ENDIF

        ENDIF

        !------------------------------------------------------------------
        ! Re-anchor the inverse-iteration shift on the energy reached so
        ! far, ONCE PER WINDOW. The shift must
        ! stay fixed within a window - EnergyIA and EnergyIB reuse the LDL^T
        ! of the leading nfru block, valid for one shift only - but between
        ! windows there is nothing to preserve, and a shift frozen for a
        ! whole OPT_CYCLE step is what makes the solver diverge late in a
        ! cycle. See RefreshINVITShift. The refresh invalidates the stored
        ! factorization, so the solve that follows starts from row 1.
        !------------------------------------------------------------------
        CALL RefreshINVITShift(cbs)
        Evalue = EnergyIA(1, cbs, .FALSE., ErrCode)
        IF ((ErrCode == 0) .AND. (ABS(Evalue) < 1.0E10_wp)) Glob_CurrEnergy = Evalue

        IF (Glob_ProcID == 0) THEN
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 2) WRITE (*, *) 'Nonlinear parameters before optimization:'
            DO i = 1, nfo
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(nfru+i), ':', Glob_PWR(nfru+i)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
            ENDDO
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Set up the step
        !------------------------------------------------------------------
        ! x holds the window's parameters flattened, x_init a copy to restore
        ! from if the step is rejected, and v_good the matching eigenvector. D
        ! is the DRMNG scale vector: the 1/(cbs^2*sqrt(cbs)) form shrinks the
        ! steps as the basis grows (a stiffer problem), floored at 10000*epsilon.
        !------------------------------------------------------------------
        IV(1:LIV) = IV_init(1:LIV)
        V(1:LV) = V_init(1:LV)

        DO i = 1, nfo
          x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
          x_init((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, nfru+i)
        ENDDO

        t = MAX(ONE/(cbs*cbs*SQRT(ONE*cbs)), 10000*EPSILON(Glob_CurrEnergy))
        DO i = 1, nfo
          ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          DO j = 1, npt
            ! Make sure none of the D(i) will be zero or smaller than the threshold
            ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j) = t
            ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          ENDDO
        ENDDO

        ExitNeeded = .FALSE.
        NumOfFailures = 0
        NumOfEnergyEval = 0
        NumOfGradEval = 0
        IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.
        E_best = Glob_CurrEnergy
        x_best(1:nfo*npt) = x(1:nfo*npt)
        ! Remember the energy and the eigenvector that correspond to the current
        ! (i.e. unchanged) values of the nonlinear parameters. They are needed in case
        ! the optimization of this function/set of functions has to be abandoned
        E_prev = Glob_CurrEnergy
        v_good(1:cbs) = Glob_LastEigvector(1:cbs)

        !------------------------------------------------------------------
        ! The reverse-communication loop
        !------------------------------------------------------------------
        ! DRMNG runs on rank 0 and IV is broadcast. IV(1) says what it wants:
        ! 1 an energy at x, 2 a gradient, 3..8 converged, 9,10 its evaluation
        ! limit. On a failed evaluation the energy handed to DRMNG keeps its
        ! previous value, so the step counts as no reduction and the trust
        ! radius shrinks (IV(2), the TOOBIG flag, is deliberately left 0). The
        ! best point is tracked here because the last point DRMNG visits is not
        ! necessarily the lowest.
        ! On this path a failure also restores v_good.
        !------------------------------------------------------------------
        DO WHILE (.NOT. (ExitNeeded))

          IF (Glob_ProcID == 0) CALL DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
          CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

          SELECT CASE (IV(1))

          CASE (1)  ! Only energy is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO i = 1, nfo
              Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
            ENDDO
            Evalue = EnergyIA(nfrup1, cbs, .TRUE., ErrCode)
            NumOfEnergyEval = NumOfEnergyEval+1
            IF (ErrCode /= 0) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
              ! Restore the last good eigenvector as the vector left by the failed
              ! inverse iteration process may be unusable as a starting vector
              Glob_LastEigvector(1:cbs) = v_good(1:cbs)
            ELSE
              Glob_CurrEnergy = Evalue
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nfo*npt) = x(1:nfo*npt)
              ENDIF
            ENDIF
            ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
            IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

          CASE (2)  ! Only gradient is needed
            CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
            DO i = 1, nfo
              Glob_NonlinParam(1:npt, nfru+i) = x((i-1)*npt+1:i*npt)
            ENDDO
            CALL EnergyIB(Evalue, grad, .TRUE., ErrCode)
            NumOfGradEval = NumOfGradEval+1
            IF (ErrCode /= 0) THEN
              NumOfFailures = NumOfFailures+1
              IV(2) = 1
              ! Restore the last good eigenvector as the vector left by the failed
              ! inverse iteration process may be unusable as a starting vector
              Glob_LastEigvector(1:cbs) = v_good(1:cbs)
            ELSE
              IF (Evalue < E_best) THEN
                E_best = Evalue
                x_best(1:nfo*npt) = x(1:nfo*npt)
              ENDIF
            ENDIF
            ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
            IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

          CASE (3:8)  ! Some kind of convergence has been reached
            ExitNeeded = .TRUE.

          CASE (9:10)  ! Function evaluation limit has been reached.
            ! This is never supposed to happen because we
            ! count the number of function evaluations ourselves.
            ExitNeeded = .TRUE.


          CASE DEFAULT
            ! DRMNG answers an IV(2) failure report with IV(1)=63 or 65,
            ! and >=14 for a bad input. None of those match a case above,
            ! so without this the loop would call DRMNG again for ever.
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
                'Warning WC0140 in OptCycleI: DRMNG returned IV(1) =', IV(1)
              IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
            ENDIF
            ExitNeeded = .TRUE.
          ENDSELECT

          ! A warning, not an abort: the best point found so far is
          ! still usable, and the acceptance tests below decide what
          ! becomes of this step.
          IF (NumOfFailures == Glob_MaxEnergyFailsAllowed) THEN
            IF (Glob_ProcID == 0) THEN
              IF (Verbose >= 1) WRITE(*, '(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0128 in OptCycleI: number of failures in energy or gradient', &
                'calculations during the optimization of nonlinear parameters', &
                'reached the limit of', Glob_MaxEnergyFailsAllowed
            ENDIF
            ! call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
          ENDIF

          IF (NumOfEnergyEval >= MaxEvalToUse) ExitNeeded = .TRUE.

        ENDDO  ! while


        !------------------------------------------------------------------
        ! Re-solve at the best point, for the linear coefficients
        !------------------------------------------------------------------
        ! EnergyIAM rather than EnergyIA because it also produces Glob_c,
        ! which the linear-coefficient test needs.
        !------------------------------------------------------------------
        ! Compute the energy and the linear coefficients at the best point found
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, nfru+i) = x_best((i-1)*npt+1:i*npt)
        ENDDO

        Evalue = EnergyIAM(nfrup1, cbs, .TRUE., ErrCode)
        IF (ErrCode == 0) THEN
          Glob_CurrEnergy = Evalue
        ELSE
          ! Note that Glob_CurrEnergy is intentionally left unchanged here as
          ! the value returned by a failed energy evaluation is meaningless
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,a,1x,a)') &
              'Warning WC0125 in OptCycleI: failed to evaluate energy after optimization', &
              'of nonlinear parameters. The values of the nonlinear parameters', &
              'will be left unchanged.'
          ENDIF
          Glob_LastEigvector(1:cbs) = v_good(1:cbs)
        ENDIF

        !------------------------------------------------------------------
        ! Acceptance test: pair overlaps
        !------------------------------------------------------------------
        ! Every pair involving a function in the window, against the whole
        ! basis. Only the LOWER triangle of Glob_S is read. Function NUMBERS
        ! are printed, not positions, because the basis has been permuted
        ! and the positions would mean nothing to the reader.
        !------------------------------------------------------------------
        ! Checking if overlap is OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad = .FALSE.
        IF (OverlapThreshold > ZERO) THEN
          ii = 0
          DO i = nfrup1, cbs
            DO j = 1, i-1
              IF (ABS(Glob_S(i, j)) > OverlapThreshold) THEN
                ii = ii+1
                IsOverlapBad = .TRUE.
                IF (Glob_ProcID == 0) THEN
                  IF (ii == 1) THEN
                    IF (Verbose >= 1) WRITE(*, *) 'Warning WC0126: overlap of the following functions exceeds threshold. ', &
                      'Nonlinear parameters will be left unchanged'
                  ENDIF
                  WRITE(*, '(1x,i6,a1,i6,i6,a6)', ADVANCE='no') &
                    ii, ':', Glob_FuncNum(i), Glob_FuncNum(j), '    S='
                  CALL writerealadv(6, Glob_S(i, j))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        !------------------------------------------------------------------
        ! Acceptance test: linear coefficients
        !------------------------------------------------------------------
        ! The scan covers the WHOLE basis, 1..cbs, not just the window:
        ! moving one function can blow up the coefficient of another, and
        ! that is exactly the near-linear-dependence this test is for.
        !------------------------------------------------------------------
        ! Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad = .FALSE.
        IF (LinCoeffThreshold > ZERO) THEN
          ii = 0
          DO i = 1, cbs
            IF (ABS(Glob_c(i)) > LinCoeffThreshold) THEN
              ii = ii+1
              IsAnyLinCoeffBad = .TRUE.
              IF (Glob_ProcID == 0) THEN
                IF (ii == 1) THEN
                  IF (Verbose >= 1) THEN
                  WRITE(*,*) 'Warning WC0127: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                  ENDIF
                ENDIF
                WRITE(*, '(1x,i6,a1,i6,a6)', ADVANCE='no') ii, ':', i, '    c='
                CALL writerealadv(6, Glob_c(i))
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        ! Reported only for a step that was actually kept. The GSEPIIS
        ! average is the figure to watch if the fixed shift starts
        ! costing convergence - see the header. Counter1 is at least 1
        ! by now, because the EnergyIAM solve above bumped it.
        IF ((Glob_ProcID == 0) .AND. (ErrCode == 0) .AND. (.NOT. IsOverlapBad) .AND. (.NOT. IsAnyLinCoeffBad)) THEN
          IF (Verbose >= 1) THEN
          WRITE (*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
          ENDIF
          WRITE (*, *) 'E=', Glob_CurrEnergy
          IF (Glob_AreParamPrintedInCycleOptX) THEN
            IF (Verbose >= 1) WRITE (*, *) 'Nonlinear parameters after optimization:'
            DO i = 1, nfo
              WRITE(*, '(1x,i6,a1,i6)', ADVANCE='no') Glob_FuncNum(nfru+i), ':', Glob_PWR(nfru+i)
              CALL writerealarradv(6, Glob_NonlinParam(1:npt, nfru+i), npt)
            ENDDO
            IF (Verbose >= 1) WRITE (*, '(1x,a41,f8.4)') 'Average number of iterations in GSEPIIS: ', &
              (Glob_InvItTempCounter2*ONE)/Glob_InvItTempCounter1
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Undo the step if it was not acceptable
        !------------------------------------------------------------------
        ! The parameters go back to x_init and the matrix elements are
        ! recomputed. A failure there can only come from the inverse iteration,
        ! so the starting vector is restored and the solve retried with
        ! AreMatElemNeeded=.false.; if that fails too the sweep carries on with
        ! the recorded energy.
        !------------------------------------------------------------------
        IF ((ErrCode /= 0) .OR. IsOverlapBad .OR. IsAnyLinCoeffBad) THEN
          ! restore initial values of nonlinear parameters
          DO i = 1, nfo
            Glob_NonlinParam(1:npt, nfru+i) = x_init((i-1)*npt+1:i*npt)
          ENDDO

          Evalue = EnergyIA(nfrup1, cbs, .TRUE., ErrCode)
          IF (ErrCode /= 0) THEN
            ! The matrix elements are correct, so the failure can only come from
            ! the inverse iteration process itself. Restore the last good
            ! eigenvector and try once more. No matrix elements need be recomputed
            Glob_LastEigvector(1:cbs) = v_good(1:cbs)
            Evalue = EnergyIA(nfrup1, cbs, .FALSE., ErrCode)
          ENDIF

          IF (ErrCode == 0) THEN
            Glob_CurrEnergy = Evalue
          ELSE
            ! Leave this function (or set of functions) as it was and proceed to
            ! the next one. The nonlinear parameters are the original ones and
            ! the matrix elements that correspond to them have been computed. The
            ! factorization of H-Glob_ApproxEnergy*S is redone at the next step,
            ! so it is safe to continue
            IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE(*, '(1x,a,1x,a)') &
              'Warning WC0129 in OptCycleI: energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_LastEigvector(1:cbs) = v_good(1:cbs)
            Glob_CurrEnergy = E_prev
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! Record the step and save
        !------------------------------------------------------------------
        ! On the last step of a cycle the position is reset to 0 and the
        ! cycle counter advances, which is what tells a resumed run that
        ! this cycle is finished.
        !
        ! Saving happens on the first few steps whatever SaveEvery says, so
        ! that a run interrupted early still leaves something behind, then
        ! every SaveEvery steps, and always on the final step of a cycle.
        !------------------------------------------------------------------
        IF (CurrFunc > FuncEnd-NumOfFuncToShift) THEN
          LastIter = .TRUE.
        ELSE
          LastIter = .FALSE.
        ENDIF

        Glob_History(cbs)%Energy = Glob_CurrEnergy
        IF (LastIter) THEN
          Glob_History(cbs)%InitFuncAtLastStep = 0
          Glob_History(cbs)%CyclesDone = Glob_History(cbs)%CyclesDone+1
        ELSE
          Glob_History(cbs)%InitFuncAtLastStep = CurrFunc
        ENDIF

        IF (Glob_ProcID == 0) THEN
          IF ((totsteps <= Glob_MinMandSavSteps) .OR. (MOD(totsteps, SaveEvery) == 0) .OR. &
              (CurrFunc+NumOfFuncToShift >= FuncEnd)) THEN
            CALL SaveResults(Sort='yes')
          ENDIF
        ENDIF

      ENDDO  ! end cycle CurrCycle


      !------------------------------------------------------------------
      ! End of a cycle
      !------------------------------------------------------------------
      ! Between cycles the basis is sorted back into function-number order
      ! so the next sweep starts from a known layout, and the
      ! factorization is rebuilt for that layout. Not done after the LAST
      ! cycle, because the final ordering below covers it.
      !------------------------------------------------------------------
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 1) WRITE(*, *)
        IF (Verbose >= 1) WRITE(*, *) 'Cycle', CurrCycle, ' finished'
      ENDIF

      IF (CurrCycle /= NumCycles) THEN
        Glob_History(cbs)%InitFuncAtLastStep = FuncBegin-NumOfFuncToShift
        IF (Glob_ProcID == 0) WRITE(*, '(1x,a47)', ADVANCE='no') &
          'Ordering basis functions and matrix elements...'
        CALL SortBasisFuncAndMatElem(fbn, cbs, FuncNumTemp, NonlinParamTemp, TempR)
        ! call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
        Glob_CurrEnergy = EnergyIA(fbn, cbs, .FALSE., ErrCode)
        IF ((ErrCode /= 0) .OR. (ABS(Glob_CurrEnergy) > 1.0E10_wp)) THEN
          IF (Glob_ProcID == 0) WRITE(*, *) &
            'Error EC0154 in OptCycleI: energy cannot be computed after sorting basis functions'
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
      ENDIF

    ENDDO  ! End of main optimization cycle


    !==================================================================
    ! Put the basis back the way it came in
    !==================================================================
    ! Sort restores function-number order and the reverse undoes the
    ! reversal applied at the start, so the basis leaves this routine in
    ! the same layout it arrived in - only with better parameters. Both
    ! run on every rank, since every rank holds its own copy.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) WRITE(*, '(1x,a53)', ADVANCE='no') &
      'Final ordering basis functions and matrix elements...'
    CALL SortBasisFuncAndMatElem(FuncBegin, cbs, FuncNumTemp, NonlinParamTemp, TempR)
    CALL ReverseFuncOrder(FuncBegin, cbs)
    CALL ReverseMatElemOrder(FuncBegin, cbs)
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'

    ! Hand H and S to the next BBOP step. On this path Glob_H holds the
    ! SHIFTED matrix, and StoreMatricesInSwapFile adds the shift back
    ! before writing - see the note there.
    CALL StoreMatricesInSwapFile()


    !==================================================================
    ! Release everything, in the reverse of the allocation order
    !==================================================================
    ! deallocate workspace
    DEALLOCATE(TempR)
    DEALLOCATE(FuncNumTemp)
    DEALLOCATE(NonlinParamTemp)
    DEALLOCATE(grad)
    DEALLOCATE(x_best)
    DEALLOCATE(x_init)
    DEALLOCATE(x)

    ! deallocate arrays used by DRMNG
    DEALLOCATE(V_init)
    DEALLOCATE(V)
    DEALLOCATE(D)

    ! deallocate workspace for SaveResults
    DEALLOCATE(Glob_IntWorkArrForSaveResults)

    ! deallocate workspace for EnergyIB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for GSEPIIS
    DEALLOCATE(Glob_WorkForGSEPIIS)
    DEALLOCATE(v_good)
    DEALLOCATE(Glob_LastEigvector)

    ! Deallocate some global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_invD)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (WrkP_WrongStateCount > 0)) THEN
      IF (Verbose >= 1) WRITE(*, *) 'Trial points refused because inverse iteration landed on'
      WRITE(*, *) 'a level other than WHICH_EIGENVALUE =', Glob_WhichEigenvalue, &
                 ' :', WrkP_WrongStateCount
    ENDIF
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine OptCycleI finished'


  END SUBROUTINE OptCycleI


  SUBROUTINE FullOpt1G(InitFunc, FinalFunc, MaxEnergyEval, OverlapThreshold, MaxOverlapPenalty, &
                       DataSaveMinTimeInterv, HessianSaveMinTimeInterv, HessFileName)
    !==================================================================
    ! Subroutine FullOpt1G
    !==================================================================
    ! Optimizes the nonlinear parameters of functions InitFunc..FinalFunc
    ! SIMULTANEOUSLY with DSYGVX ('G') and DRMNG; FullOpt1I is the
    ! inverse-iteration twin. The block is moved to the end of the basis on
    ! entry and back on exit. Overlaps are handled by a smooth quadratic
    ! PENALTY added to the energy (OverlapThreshold >= 1.0 turns it off);
    ! the energy printed during the optimization includes it. The Hessian
    ! is read from HessFileName on entry; after an improving evaluation
    ! the data file is written when DataSaveMinTimeInterv seconds have
    ! passed since its last save and the Hessian when
    ! HessianSaveMinTimeInterv seconds have; both are written once more
    ! at the end. A HessFileName of ' ', 'none', 'NONE' or 'None'
    ! disables the Hessian file. MaxEnergyEval is the total budget of
    ! energy evaluations of this basis size, the history count included;
    ! DRMNG enforces the remainder (IV(17)). The step ends at the best
    ! accepted point, never on a rejected trial point.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)                        :: InitFunc, FinalFunc       ! block to optimize
    INTEGER, INTENT(IN)                        :: MaxEnergyEval             ! total budget, history included
    REAL(wp), INTENT(IN)                       :: OverlapThreshold          ! penalty threshold, >=1 = off
    REAL(wp), INTENT(IN)                       :: MaxOverlapPenalty         ! penalty magnitude
    REAL(4), INTENT(IN)                        :: DataSaveMinTimeInterv     ! seconds between data saves
    REAL(4), INTENT(IN)                        :: HessianSaveMinTimeInterv  ! seconds between Hessian saves
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: HessFileName

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and geometry ---------------------------------------
    INTEGER :: i, j                ! loop counters
    INTEGER :: npt                 ! Glob_npt, parameters per function
    INTEGER :: nfa                 ! Glob_CurrBasisSize
    INTEGER :: nfo                 ! functions being optimized
    INTEGER :: nv                  ! nfo*npt, optimization variables
    INTEGER :: tas                 ! size of the block handed to the permutations
    INTEGER :: InitFuncNew         ! where InitFunc sits after the shift to the end
    INTEGER :: BlockSizeForDSYGVX  ! ILAENV block size for the LAPACK workspace

    ! -- energies ----------------------------------------------------
    REAL(wp) :: Evalue             ! value returned by the latest solve
    REAL(wp) :: CurrentEnergy      ! the value DRMNG is driven with
    REAL(wp) :: t                  ! scale factor for the DRMNG vector D
    REAL(wp) :: MaxAbsOverlap      ! reported by GetOverlapStatistics
    REAL(wp) :: MinAbsOverlap      !   "
    REAL(wp) :: AverageAbsOverlap  !   "

    ! -- status flags ------------------------------------------------
    LOGICAL :: IsSwapFileOK       ! H and S came from the swap file
    LOGICAL :: ExitNeeded         ! ends the reverse-communication loop
    LOGICAL :: SaveHessian        ! HessFileName names a real file
    LOGICAL :: IsHessFileOK       ! a usable Hessian was read back
    LOGICAL :: IsHessSaveSuccess  ! set by SaveHessianFile, not read

    ! -- solver bookkeeping ------------------------------------------
    INTEGER :: ErrCode                            ! non-zero when a solve failed
    INTEGER :: NumOfFailures                      ! failed solves so far
    INTEGER :: NumOfEnergyEval                    ! energy evaluations so far
    INTEGER :: NumOfGradEval                      ! gradient evaluations so far
    INTEGER :: NumOfEnergyEvalDuringFullOpt_Init  ! count carried in from history

    ! -- timing ------------------------------------------------------
    REAL(4) :: TimeOfLastSave      ! CPU time of the last data save
    REAL(4) :: TimeOfLastHessSave  ! CPU time of the last Hessian save

    ! -- declared but NOT REFERENCED ---------------------------------
    INTEGER :: OpenFileErr
    INTEGER :: np
    INTEGER :: nfru

    ! -- workspace for the permutation routines ----------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: NonlinParamTemp
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: FuncNumTemp
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TempR

    ! -- the optimization variables ----------------------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, grad

    ! -- arrays and settings used by DRMNG (rank 0 only) -------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D         ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V         ! work array
    INTEGER, PARAMETER                  :: LIV = 60  ! length of IV
    INTEGER                             :: IV(LIV)
    INTEGER                             :: LV        ! length of V
    INTEGER                             :: ALG       ! 2 = unconstrained minimization
    INTEGER                             :: IVLMAT    ! IV(42), where V holds the Hessian
    !!====================================================
    !!These variables are used when a finite difference gradient is computed
    ! real(wp)                                    deltax,Evalue1
    !!====================================================


    !==================================================================
    ! Overlap penalty
    !==================================================================
    IsHessFileOK = .FALSE.

    IF (OverlapThreshold >= ONE) THEN
      Glob_OverlapPenaltyAllowed = .FALSE.
    ELSE
      Glob_OverlapPenaltyAllowed = .TRUE.
      Glob_OverlapPenaltyThreshold2 = OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty = MaxOverlapPenalty
    ENDIF

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine FullOpt1G started'
      IF (Verbose >= 1) WRITE(*, *) 'Simultaneous optimization of nonlinear parameters of basis functions'
      IF (Verbose >= 1) WRITE(*, *) InitFunc, '  through', FinalFunc, '  will be attempted'
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Overlap threshold is ', ABS(OverlapThreshold)
        IF (Verbose >= 1) WRITE(*, *) 'Max value of a pair overlap penalty is ', Glob_MaxOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Warning! The energy value that will be shown during the optimization'
        WRITE(*, *) 'may differ from the actual energy'
      ELSE
        IF (Verbose >= 1) WRITE(*, *) 'No constraints on overlaps will be imposed'
      ENDIF
    ENDIF


    !==================================================================
    ! Global state and array allocation
    !==================================================================
    ! The whole basis is the problem here; nfo of it is being varied.
    !------------------------------------------------------------------
    ! Setting the values of some global variables
    Glob_GSEPSolutionMethod = 'G'
    Glob_nfa = Glob_CurrBasisSize
    Glob_nfo = FinalFunc-InitFunc+1
    Glob_nfru = Glob_CurrBasisSize-Glob_nfo
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    nfa = Glob_nfa
    nfo = Glob_nfo
    nfru = Glob_nfru
    nv = nfo*npt
    InitFuncNew = nfa+InitFunc-FinalFunc
    Glob_HSBuffLen = MAX(MIN(nfa*(nfa+1)/2, 1000), 30*nfa)

    ! Allocate some global arrays
    ALLOCATE(Glob_H(nfa, nfa))
    ALLOCATE(Glob_S(nfa, nfa))
    ALLOCATE(Glob_diagH(nfa))
    ALLOCATE(Glob_diagS(nfa))
    ALLOCATE(Glob_D(2*npt, nfo, nfa))
    ALLOCATE(Glob_c(nfa))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', nfa, nfa, nfa, nfa)
    Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*nfa, 8*nfa)
    ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    ALLOCATE(Glob_IWorkForDSYGVX(5*nfa))

    ! Allocate workspace for EnergyGB
    ALLOCATE(Glob_WkGR(nfo*npt))

    ! D and V are only ever touched on rank 0, which is the only rank
    ! that runs DRMNG, so they are allocated there alone.
    ! Allocate arrays used by DRMNG
    LV = 71 + nv*(nv+13)/2 + 1
    IF (Glob_ProcID == 0) THEN
      ALLOCATE(D(nv))
      ALLOCATE(V(LV))
    ENDIF

    ! Allocate workspace
    ALLOCATE(x(nv))
    ALLOCATE(grad(nv))

    ! Setting up a logical variable that determines
    ! whether the hessian should be saved from time to time
    IF ((HessFileName == ' ') .OR. (HessFileName == 'none') .OR. &
        (HessFileName == 'NONE') .OR. (HessFileName == 'None')) THEN
      SaveHessian = .FALSE.
    ELSE
      SaveHessian = .TRUE.
    ENDIF


    !==================================================================
    ! Move the block to the end of the basis
    !==================================================================
    ! The matrix elements are permuted along only when the swap file
    ! supplied usable ones; otherwise EnergyGA recomputes them below,
    ! already in the permuted order.
    !------------------------------------------------------------------
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    ! Shifting basis functions InitFunc through FinalFunc to
    ! the very end and, if necessary, doing the permutation of
    ! matrix elements to reflect this change in function order.
    IF (FinalFunc /= nfa) THEN
      ! allocate space
      tas = FinalFunc-InitFunc+1
      ALLOCATE(NonlinParamTemp(1:npt, 1:tas))
      ALLOCATE(FuncNumTemp(1:tas))
      ALLOCATE(TempR(1:tas))
      CALL PermuteFunctions(InitFunc, FinalFunc, FuncNumTemp, NonlinParamTemp)
      IF (IsSwapFileOK) CALL PermuteMatrixElements(InitFunc, FinalFunc, TempR)
      DEALLOCATE(TempR)
      DEALLOCATE(FuncNumTemp)
      DEALLOCATE(NonlinParamTemp)
      ! Allocate workspace for SaveResults (it uses Sort='yes'
      ! option, which requires workspace)
      ALLOCATE(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
    ENDIF


    !==================================================================
    ! Initial energy
    !==================================================================
    ! Calculating the initial energy
    IF (IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyGA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a59)', ADVANCE='no') &
        'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyGA(1, Glob_CurrBasisSize, .TRUE., ErrCode)
    ENDIF
    IF (ErrCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0160 in FullOpt1G: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'

    IF (Glob_ProcID == 0) THEN
      CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                           ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (including overlap penalty ', Glob_CurrEnergy
      ELSE
        WRITE(*, *) 'Initial energy                            ', Glob_CurrEnergy
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                           ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                           ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap              ', AverageAbsOverlap
    ENDIF

    CALL CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave = Glob_TimeSinceStart
    TimeOfLastHessSave = Glob_TimeSinceStart


    !==================================================================
    ! Set up DRMNG and load the Hessian
    !==================================================================
    ! IV(25)=0 tells DRMNG not to overwrite the Hessian it was handed.
    ! The scale vector D is only built when there is no usable stored
    ! one to reuse.
    !------------------------------------------------------------------
    ! Setting parameters for DRMNG

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    IF (Glob_ProcID == 0) THEN
      CALL DIVSET(ALG, IV, LIV, LV, V)
      IV(18) = 1000000  ! iteration limit; the evaluation limit IV(17) is set from the budget below
      IV(19) = -1  ! set summary print format
      IV(20) = 0; IV(22) = 0; IV(23) = -1; IV(24) = 0
      V(31) = 0.0_wp
      V(32) = 2*EPSILON(V(32))
      V(37) = 2*EPSILON(V(37))
      ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
      ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
      ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
      V(35) = Glob_MaxScStepAllowedInOpt
      IV(1) = 12  ! DIVSET has been called and some default values were changed
      nv = nfo*npt
    ENDIF

    DO i = 1, nfo
      x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, InitFuncNew+i-1)
    ENDDO

    IF (Glob_ProcID == 0) THEN
      ! If SaveHessian=.true. then try to read the Hessian
      ! from the file
      IF (SaveHessian) THEN
        IVLMAT = IV(42)
        CALL ReadHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessFileOK)
        IF (IsHessFileOK) IV(25) = 0
      ENDIF
      IF ((.NOT. Glob_FullOptSaveD) .OR. (.NOT. IsHessFileOK) .OR. (.NOT. SaveHessian)) THEN
        !    !Set the scaling vector
        !    do i=1,nfo
        !      t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
        !      do j=1,np
        !        !Make sure none of the D(i) will be zero or smaller than the threshold
        !        D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        !      enddo
        !    enddo
        ! Set the scaling vector
        t = MAX(ONE/(nfa*nfa*SQRT(ONE*nfa)), 10000*EPSILON(Glob_CurrEnergy))
        DO i = 1, nfo
          ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          DO j = 1, npt
            ! Make sure none of the D(i) will be zero or smaller than the threshold
            ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j) = t
            ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF


    !==================================================================
    ! The reverse-communication loop
    !==================================================================
    ! IV(1) on return: 1 = wants an energy, 2 = wants a gradient,
    ! 3..8 = converged, 9,10 = its own limit reached. A failed energy
    ! is reported with IV(2)=1, which makes DRMNG shrink the step.
    !------------------------------------------------------------------
    ExitNeeded = .FALSE.
    NumOfFailures = 0
    NumOfEnergyEval = 0
    NumOfGradEval = 0
    NumOfEnergyEvalDuringFullOpt_Init = Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt
    ! The remaining budget of energy evaluations of this basis size is
    ! DRMNG's own function-evaluation limit (checked after every
    ! evaluation, so the count is exact)
    IF (Glob_ProcID == 0) IV(17) = MAX(1, MaxEnergyEval-NumOfEnergyEvalDuringFullOpt_Init)
    ! DRMNG does not read FX on its first (IV(1)=12) entry, but giving
    ! it a defined value keeps -finit-real=nan builds quiet.
    CurrentEnergy = Glob_CurrEnergy

    DO WHILE (.NOT. (ExitNeeded))

      IF (Glob_ProcID == 0) CALL DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
      CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      SELECT CASE (IV(1))

      CASE (1)  ! Only energy is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
        ENDDO
        Evalue = EnergyGA(InitFuncNew, nfa, .TRUE., ErrCode)
        NumOfEnergyEval = NumOfEnergyEval+1
        IF (ErrCode /= 0) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
        ELSE
          CurrentEnergy = Evalue
        ENDIF
        ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
        IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

        !--------------------------------------------------------------
        ! Periodic saving, on rank 0, only when the energy improved
        !--------------------------------------------------------------
        IF (Glob_ProcID == 0) THEN
          IF (Evalue < Glob_CurrEnergy) THEN
            ! Save data and Hessian if necessary
            CALL CPU_TIME(Glob_TimeSinceStart)

            IF (Glob_TimeSinceStart-TimeOfLastSave > DataSaveMinTimeInterv) THEN
              ! Save the results if more than DataSaveMinTimeInterv seconds
              ! have passed since the last save
              IF (Glob_OverlapPenaltyAllowed) THEN
                Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
              ELSE
                Glob_CurrEnergy = Evalue
              ENDIF
              ! Changing history
              Glob_History(Glob_CurrBasisSize)%Energy = Glob_CurrEnergy
              Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt = &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              ! Sort='yes' needs the workspace that is allocated only when
              ! the block was actually moved, i.e. when FinalFunc/=nfa.
              IF (FinalFunc == nfa) THEN
                CALL SaveResults(Sort='no')
              ELSE
                CALL SaveResults(Sort='yes')
              ENDIF
              WRITE(*, *) 'Data file has been updated'
              CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
              IF (Verbose >= 1) WRITE(*, *) 'Some current statistics:'
              IF (Glob_OverlapPenaltyAllowed) THEN
                IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
              ELSE
                WRITE(*, *) 'Energy                            ', Evalue
              ENDIF
              IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
              TimeOfLastSave = Glob_TimeSinceStart
            ENDIF

            IF ((Glob_TimeSinceStart-TimeOfLastHessSave > HessianSaveMinTimeInterv) &
                .AND. (SaveHessian)) THEN
              ! Save the Hessian if more than HessianSaveMinTimeInterv seconds
              ! have passed since the last save
              CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
              TimeOfLastHessSave = Glob_TimeSinceStart
            ENDIF

          ENDIF
        ENDIF

      CASE (2)  ! Only gradient is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
        ENDDO
        CALL EnergyGB(Evalue, grad, .TRUE., ErrCode)
        NumOfGradEval = NumOfGradEval+1
        IF (ErrCode /= 0) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
        ENDIF
        ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
        IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1
      !===================================
      ! Finite-difference check of the gradient (debug code) was removed
      ! here; see NEW_workproc.f90.bak4_20260924 if it is needed again.
      !===================================

      CASE (3:8)  ! Some kind of convergence has been reached
        ExitNeeded = .TRUE.

      CASE (9:10)  ! DRMNG reached its evaluation (or iteration) limit
        ! IV(17) holds the remaining energy-evaluation budget of this basis
        ! size, so this is the normal end of a step that exhausts it. x may
        ! still sit on a rejected trial point; the code after the loop goes
        ! back to the accepted point when that is the case.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *) 'Warning WC0130 in FullOpt1G: number of energy evaluations reached limit'
          IF (Verbose >= 1) WRITE(*, '(1x,a,i0,a,i0,a)') '(', NumOfEnergyEval, ' in this step, ', &
            NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval, ' in total)'
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.

      CASE DEFAULT
        ! DRMNG returns 63 or 65 when it gives up on an uncomputable
        ! value, and >=14 for a bad input. None of those match a case
        ! above, so without this the loop would call DRMNG again for
        ! ever. Leave with the best point found so far.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
            'Warning WC0136 in FullOpt1G: DRMNG returned IV(1) =', IV(1)
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.

      ENDSELECT

      IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,a,1x,a)') &
            'Error EC0161 in FullOpt1G: number of failures in energy or gradient', &
            'calculations during the optimization of nonlinear parameters', &
            'exceeded limit'
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF

    ENDDO


    !==================================================================
    ! Final energy at the best point found
    !==================================================================
    ! DRMNG can leave x on a rejected trial point of the last line search
    ! (evaluation limit, false convergence). Its accepted iterate x0 is in
    ! V(IV(43):IV(43)+nv-1) with f0 = V(13); V(10) is the last f evaluated.
    ! Whenever the point in x is not better than x0, go back to x0, so
    ! that the final energy and the data file describe the best point.
    IF (Glob_ProcID == 0) THEN
      IF ((IV(31) >= 1) .AND. (V(10) >= V(13))) THEN
        IF ((V(10) > V(13)) .AND. (Verbose >= 2)) &
          WRITE(*, *) 'The last trial point is discarded: back to the last accepted point'
        x(1:nv) = V(IV(43):IV(43)+nv-1)
      ENDIF
    ENDIF
    CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    DO i = 1, nfo
      Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
    ENDDO
    Evalue = EnergyGA(InitFuncNew, nfa, .TRUE., ErrCode)
    IF (ErrCode /= 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0162 in FullOpt1G: failed to evaluate energy after the optimization'
        IF (Verbose >= 1) WRITE(*, *) 'of nonlinear parameters'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! Glob_CurrEnergy is the physical energy, without the penalty
    IF (Glob_OverlapPenaltyAllowed) THEN
      Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
    ELSE
      Glob_CurrEnergy = Evalue
    ENDIF

    ! Printing the number of energy/gradient evaluations
    ! and the energy after the optimization
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
      WRITE(*, *) 'Final energy and overlap statistics:'
      CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
      ELSE
        WRITE(*, *) 'Energy                            ', Evalue
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
    ENDIF

    ! Adding data to history
    Glob_History(Glob_CurrBasisSize)%Energy = Glob_CurrEnergy
    Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt = &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval


    !==================================================================
    ! Move the block back and release everything
    !==================================================================
    ! Shifting basis functions InitFunc through FinalFunc back from
    ! the very end to the middle, where they initially were, and, if
    ! necessary, doing the permutation of matrix elements to reflect
    ! this change in function order.
    IF (FinalFunc /= nfa) THEN
      ! allocate space
      tas = nfa-FinalFunc
      ALLOCATE(NonlinParamTemp(1:npt, 1:tas))
      ALLOCATE(FuncNumTemp(1:tas))
      ALLOCATE(TempR(1:tas))
      CALL PermuteFunctions(InitFunc, InitFunc+tas-1, FuncNumTemp, NonlinParamTemp)
      IF (Glob_UseSwapFile) CALL PermuteMatrixElements(InitFunc, InitFunc+tas-1, TempR)
      DEALLOCATE(TempR)
      DEALLOCATE(FuncNumTemp)
      DEALLOCATE(NonlinParamTemp)
      ! deallocate workspace for SaveResults
      DEALLOCATE(Glob_IntWorkArrForSaveResults)
    ENDIF

    ! saving results; the Hessian goes to disk with them, so a restart
    ! finds the file even when no evaluation improved the energy
    IF (Glob_ProcID == 0) THEN
      CALL SaveResults(Sort='no')
      IF (SaveHessian) CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
    ENDIF

    CALL StoreMatricesInSwapFile()

    IF (Glob_OverlapPenaltyAllowed) Glob_OverlapPenaltyAllowed = .FALSE.

    ! deallocate workspace
    DEALLOCATE(grad)
    DEALLOCATE(x)

    ! deallocate arrays used by DRMNG
    IF (Glob_ProcID == 0) THEN
      DEALLOCATE(V)
      DEALLOCATE(D)
    ENDIF

    ! deallocate workspace for EnergyGB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for DSYGVX
    DEALLOCATE(Glob_IWorkForDSYGVX)
    DEALLOCATE(Glob_WorkForDSYGVX)

    ! Deallocate some global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine FullOpt1G finished'


  END SUBROUTINE FullOpt1G


  SUBROUTINE FullOpt1I(InitFunc, FinalFunc, MaxEnergyEval, OverlapThreshold, MaxOverlapPenalty, &
                       DataSaveMinTimeInterv, HessianSaveMinTimeInterv, HessFileName)
    !==================================================================
    ! Subroutine FullOpt1I
    !==================================================================
    ! Optimizes the nonlinear parameters of functions InitFunc..FinalFunc
    ! SIMULTANEOUSLY with INVERSE ITERATION ('I') and DRMNG; a twin of
    ! FullOpt1G (block moved to the end of the basis, quadratic overlap
    ! PENALTY, Hessian and data file saved periodically - see there). v_good
    ! holds Glob_LastEigvector from the last successful solve and is
    ! restored when one fails (read only on a failure path). While the
    ! WHOLE basis is optimized (nfru==0) EnergyIB re-anchors the shift on
    ! the best energy so far (WrkP_RefreshShiftInIB, WrkP_LastINVITEnergy).
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN)                        :: InitFunc, FinalFunc       ! block to optimize
    INTEGER, INTENT(IN)                        :: MaxEnergyEval             ! total budget, history included
    REAL(wp), INTENT(IN)                       :: OverlapThreshold          ! penalty threshold, >=1 = off
    REAL(wp), INTENT(IN)                       :: MaxOverlapPenalty         ! penalty magnitude
    REAL(4), INTENT(IN)                        :: DataSaveMinTimeInterv     ! seconds between data saves
    REAL(4), INTENT(IN)                        :: HessianSaveMinTimeInterv  ! seconds between Hessian saves
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: HessFileName

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and geometry ---------------------------------------
    INTEGER :: i, j         ! loop counters
    INTEGER :: npt          ! Glob_npt, parameters per function
    INTEGER :: nfa          ! Glob_CurrBasisSize
    INTEGER :: nfo          ! functions being optimized
    INTEGER :: nv           ! nfo*npt, optimization variables
    INTEGER :: tas          ! size of the block handed to the permutations
    INTEGER :: InitFuncNew  ! where InitFunc sits after the shift to the end

    ! -- energies ----------------------------------------------------
    REAL(wp) :: Evalue             ! value returned by the latest solve
    REAL(wp) :: CurrentEnergy      ! the value DRMNG is driven with
    REAL(wp) :: t                  ! scale factor for the DRMNG vector D
    REAL(wp) :: MaxAbsOverlap      ! reported by GetOverlapStatistics
    REAL(wp) :: MinAbsOverlap      !   "
    REAL(wp) :: AverageAbsOverlap  !   "

    ! -- status flags ------------------------------------------------
    LOGICAL :: IsSwapFileOK       ! H and S came from the swap file
    LOGICAL :: ExitNeeded         ! ends the reverse-communication loop
    LOGICAL :: SaveHessian        ! HessFileName names a real file
    LOGICAL :: IsHessFileOK       ! a usable Hessian was read back
    LOGICAL :: IsHessSaveSuccess  ! set by SaveHessianFile, not read

    ! -- solver bookkeeping ------------------------------------------
    INTEGER :: ErrCode                            ! non-zero when a solve failed
    INTEGER :: NumOfFailures                      ! failed solves so far
    INTEGER :: NumOfEnergyEval                    ! energy evaluations so far
    INTEGER :: NumOfGradEval                      ! gradient evaluations so far
    INTEGER :: NumOfEnergyEvalDuringFullOpt_Init  ! count carried in from history

    ! -- timing ------------------------------------------------------
    REAL(4) :: TimeOfLastSave      ! CPU time of the last data save
    REAL(4) :: TimeOfLastHessSave  ! CPU time of the last Hessian save

    ! -- declared but NOT REFERENCED ---------------------------------
    INTEGER :: OpenFileErr
    INTEGER :: np
    INTEGER :: nfru

    ! -- workspace for the permutation routines ----------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: NonlinParamTemp
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: FuncNumTemp
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TempR

    ! -- the optimization variables ----------------------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: x, grad

    ! -- inverse-iteration fallback ----------------------------------
    !    Copy of Glob_LastEigvector from the last SUCCESSFUL solve.
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: v_good

    ! -- arrays and settings used by DRMNG (rank 0 only) -------------
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: D         ! scale vector
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: V         ! work array
    INTEGER, PARAMETER                  :: LIV = 60  ! length of IV
    INTEGER                             :: IV(LIV)
    INTEGER                             :: LV        ! length of V
    INTEGER                             :: ALG       ! 2 = unconstrained minimization
    INTEGER                             :: IVLMAT    ! IV(42), where V holds the Hessian


    !==================================================================
    ! Overlap penalty
    !==================================================================
    IsHessFileOK = .FALSE.

    IF (OverlapThreshold >= ONE) THEN
      Glob_OverlapPenaltyAllowed = .FALSE.
    ELSE
      Glob_OverlapPenaltyAllowed = .TRUE.
      Glob_OverlapPenaltyThreshold2 = OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty = MaxOverlapPenalty
    ENDIF

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine FullOpt1I started'
      IF (Verbose >= 1) WRITE(*, *) 'Simultaneous optimization of nonlinear parameters of basis functions'
      IF (Verbose >= 1) WRITE(*, *) InitFunc, '  through', FinalFunc, '  will be attempted'
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Overlap threshold is ', ABS(OverlapThreshold)
        IF (Verbose >= 1) WRITE(*, *) 'Max value of a pair overlap penalty is ', Glob_MaxOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Warning! The energy value that will be shown during the optimization'
        WRITE(*, *) 'may differ from the actual energy'
      ELSE
        IF (Verbose >= 1) WRITE(*, *) 'No constraints on overlaps will be imposed'
      ENDIF
    ENDIF


    !==================================================================
    ! Global state and array allocation
    !==================================================================
    ! The whole basis is the problem here; nfo of it is being varied.
    !------------------------------------------------------------------
    ! Setting the values of some global variables
    Glob_GSEPSolutionMethod = 'I'
    Glob_nfa = Glob_CurrBasisSize
    Glob_nfo = FinalFunc-InitFunc+1
    Glob_nfru = Glob_CurrBasisSize-Glob_nfo
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    nfa = Glob_nfa
    nfo = Glob_nfo
    nfru = Glob_nfru
    nv = nfo*npt
    InitFuncNew = nfa+InitFunc-FinalFunc
    Glob_HSBuffLen = MAX(MIN(nfa*(nfa+1)/2, 1000), 30*nfa)

    ! Allocate some global arrays
    ALLOCATE(Glob_H(nfa, nfa))
    ALLOCATE(Glob_S(nfa, nfa))
    ALLOCATE(Glob_diagS(nfa))
    ALLOCATE(Glob_invD(nfa))
    ALLOCATE(Glob_D(2*npt, nfo, nfa))
    ALLOCATE(Glob_c(nfa))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DkBuff2(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff1(2*npt, Glob_HSBuffLen))
    ALLOCATE(Glob_DlBuff2(2*npt, Glob_HSBuffLen))

    ! Allocate workspace for subroutine GSEPIIS, which is called
    ! inside EnergyIA and EnergyIB
    ALLOCATE(Glob_WorkForGSEPIIS(nfa))
    ALLOCATE(Glob_LastEigvector(nfa))
    Glob_LastEigvector(1:nfa) = ONE

    ! A failed inverse iteration leaves Glob_LastEigvector unusable as
    ! a starting vector, so the last good one is kept here.
    ALLOCATE(v_good(nfa))
    v_good(1:nfa) = ONE

    ! Allocate workspace for EnergyIB
    ALLOCATE(Glob_WkGR(nfo*npt))

    ! D and V are only ever touched on rank 0, which is the only rank
    ! that runs DRMNG, so they are allocated there alone.
    ! Allocate arrays used by DRMNG
    LV = 71 + nv*(nv+13)/2 + 1
    IF (Glob_ProcID == 0) THEN
      ALLOCATE(D(nv))
      ALLOCATE(V(LV))
    ENDIF

    ! Allocate workspace
    ALLOCATE(x(nv))
    ALLOCATE(grad(nv))

    ! Setting up a logical variable that determines
    ! whether the hessian should be saved from time to time
    IF ((HessFileName == ' ') .OR. (HessFileName == 'none') .OR. &
        (HessFileName == 'NONE') .OR. (HessFileName == 'None')) THEN
      SaveHessian = .FALSE.
    ELSE
      SaveHessian = .TRUE.
    ENDIF


    !==================================================================
    ! Move the block to the end of the basis
    !==================================================================
    ! The matrix elements are permuted along only when the swap file
    ! supplied usable ones; otherwise EnergyIA recomputes them below,
    ! already in the permuted order.
    !------------------------------------------------------------------
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    ! Shifting basis functions InitFunc through FinalFunc to
    ! the very end and, if necessary, doing the permutation of
    ! matrix elements to reflect this change in function order.
    IF (FinalFunc /= nfa) THEN
      ! allocate space
      tas = FinalFunc-InitFunc+1
      ALLOCATE(NonlinParamTemp(1:npt, 1:tas))
      ALLOCATE(FuncNumTemp(1:tas))
      ALLOCATE(TempR(1:tas))
      CALL PermuteFunctions(InitFunc, FinalFunc, FuncNumTemp, NonlinParamTemp)
      IF (IsSwapFileOK) CALL PermuteMatrixElements(InitFunc, FinalFunc, TempR)
      DEALLOCATE(TempR)
      DEALLOCATE(FuncNumTemp)
      DEALLOCATE(NonlinParamTemp)
      ! Allocate workspace for SaveResults (it uses Sort='yes'
      ! option, which requires workspace)
      ALLOCATE(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
    ENDIF


    !==================================================================
    ! Initial energy
    !==================================================================
    ! Calculating the initial energy
    IF (IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a59)', ADVANCE='no') &
        'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .TRUE., ErrCode)
    ENDIF
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'

    !==================================================================
    ! Put the inverse-iteration shift on the requested eigenvalue
    !==================================================================
    ! Once per BBOP step; inert unless Glob_EigIdxTargeting==1. It runs
    ! BEFORE the fatal check below: a shift sitting between two
    ! eigenvalues at nearly equal distance is the usual reason the first
    ! solve does not converge, and moving the shift is precisely the
    ! cure.
    !------------------------------------------------------------------
    IF (Glob_EigIdxTargeting == 1) THEN
      CALL RetargetShiftToEigenvalue(Glob_CurrBasisSize, 'FullOpt1I')
      Glob_CurrEnergy = EnergyIA(1, Glob_CurrBasisSize, .FALSE., ErrCode)
    ENDIF

    ! The second test catches a solve that converged on a level other
    ! than WHICH_EIGENVALUE: EnergyIA returns it with ErrCode=0 and the
    ! rejection sentinel as the energy - see IsRequestedEigenstate.
    IF ((ErrCode /= 0) .OR. (ABS(Glob_CurrEnergy) > 1.0E10_wp)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0165 in FullOpt1I: initial energy cannot be computed'
        IF (ErrCode == 0) WRITE(*, *) '(inverse iteration converged on a level other than WHICH_EIGENVALUE)'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                           ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Initial energy (including overlap penalty)', Glob_CurrEnergy
      ELSE
        WRITE(*, *) 'Initial energy                            ', Glob_CurrEnergy
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                           ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                           ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap              ', AverageAbsOverlap
    ENDIF

    CALL CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave = Glob_TimeSinceStart
    TimeOfLastHessSave = Glob_TimeSinceStart


    !==================================================================
    ! Set up DRMNG and load the Hessian
    !==================================================================
    ! IV(25)=0 tells DRMNG not to overwrite the Hessian it was handed.
    ! The scale vector D is only built when there is no usable stored
    ! one to reuse.
    !------------------------------------------------------------------
    ! Setting parameters for DRMNG

    ! Call DIVSET to get default values in IV and V arrays
    ! ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG = 2
    IF (Glob_ProcID == 0) THEN
      CALL DIVSET(ALG, IV, LIV, LV, V)
      IV(18) = 1000000  ! iteration limit; the evaluation limit IV(17) is set from the budget below
      IV(19) = -1  ! set summary print format
      IV(20) = 0; IV(22) = 0; IV(23) = -1; IV(24) = 0
      V(31) = 0.0_wp
      V(32) = 2*EPSILON(V(32))
      V(37) = 2*EPSILON(V(37))
      ! V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
      ! VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
      ! MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
      V(35) = Glob_MaxScStepAllowedInOpt
      IV(1) = 12  ! DIVSET has been called and some default values were changed
      nv = nfo*npt
    ENDIF

    DO i = 1, nfo
      x((i-1)*npt+1:i*npt) = Glob_NonlinParam(1:npt, InitFuncNew+i-1)
    ENDDO

    IF (Glob_ProcID == 0) THEN
      ! If SaveHessian=.true. then try to read the Hessian
      ! from the file
      IF (SaveHessian) THEN
        IVLMAT = IV(42)
        CALL ReadHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessFileOK)
        IF (IsHessFileOK) IV(25) = 0
      ENDIF
      IF ((.NOT. Glob_FullOptSaveD) .OR. (.NOT. IsHessFileOK) .OR. (.NOT. SaveHessian)) THEN
        !    !Set the scaling vector
        !    do i=1,nfo
        !      t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
        !      do j=1,np
        !        !Make sure none of the D(i) will be zero or smaller than the threshold
        !        D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
        !      enddo
        !    enddo
        ! Set the scaling vector
        t = MAX(ONE/(nfa*nfa*SQRT(ONE*nfa)), 10000*EPSILON(Glob_CurrEnergy))
        DO i = 1, nfo
          ! t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          DO j = 1, npt
            ! Make sure none of the D(i) will be zero or smaller than the threshold
            ! D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j) = t
            ! write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF


    !==================================================================
    ! The reverse-communication loop
    !==================================================================
    ! IV(1) on return: 1 = wants an energy, 2 = wants a gradient,
    ! 3..8 = converged, 9,10 = its own limit reached. A failed energy
    ! is reported with IV(2)=1, which makes DRMNG shrink the step.
    !------------------------------------------------------------------
    ExitNeeded = .FALSE.
    NumOfFailures = 0
    NumOfEnergyEval = 0
    NumOfGradEval = 0
    NumOfEnergyEvalDuringFullOpt_Init = Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt
    ! The remaining budget of energy evaluations of this basis size is
    ! DRMNG's own function-evaluation limit (checked after every
    ! evaluation, so the count is exact)
    IF (Glob_ProcID == 0) IV(17) = MAX(1, MaxEnergyEval-NumOfEnergyEvalDuringFullOpt_Init)
    ! DRMNG does not read FX on its first (IV(1)=12) entry, but giving
    ! it a defined value keeps -finit-real=nan builds quiet.
    !------------------------------------------------------------------
    ! Shift refresh inside EnergyIB: the best energy seen so far seeds it,
    ! and the switch is on only while this optimization runs and only pays
    ! off when the whole basis is optimized (nfru==0) - see EnergyIB. The
    ! wrong-level tally restarts with the step.
    !------------------------------------------------------------------
    WrkP_LastINVITEnergy = Glob_CurrEnergy
    IF (Glob_OverlapPenaltyAllowed) WrkP_LastINVITEnergy = Glob_CurrEnergy-Glob_TotalOverlapPenalty
    WrkP_RefreshShiftInIB = (Glob_nfru == 0)
    WrkP_WrongStateCount = 0

    CurrentEnergy = Glob_CurrEnergy
    v_good(1:nfa) = Glob_LastEigvector(1:nfa)

    DO WHILE (.NOT. (ExitNeeded))

      IF (Glob_ProcID == 0) CALL DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
      CALL MPI_BCAST(IV, LIV, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      SELECT CASE (IV(1))

      CASE (1)  ! Only energy is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
        ENDDO
        Evalue = EnergyIA(InitFuncNew, nfa, .TRUE., ErrCode)
        NumOfEnergyEval = NumOfEnergyEval+1
        IF (ErrCode /= 0) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
          ! Restore the last good eigenvector as the vector left by the failed
          ! inverse iteration process may be unusable as a starting vector
          Glob_LastEigvector(1:nfa) = v_good(1:nfa)
        ELSE
          CurrentEnergy = Evalue
          v_good(1:nfa) = Glob_LastEigvector(1:nfa)
        ENDIF
        ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
        IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

        !--------------------------------------------------------------
        ! Periodic saving, on rank 0, only when the energy improved
        !--------------------------------------------------------------
        IF (Glob_ProcID == 0) THEN
          IF (Evalue < Glob_CurrEnergy) THEN
            ! Save data and Hessian if necessary
            CALL CPU_TIME(Glob_TimeSinceStart)

            IF (Glob_TimeSinceStart-TimeOfLastSave > DataSaveMinTimeInterv) THEN
              ! Save the results if more than DataSaveMinTimeInterv seconds
              ! have passed since the last save
              IF (Glob_OverlapPenaltyAllowed) THEN
                Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
              ELSE
                Glob_CurrEnergy = Evalue
              ENDIF
              ! Changing history
              Glob_History(Glob_CurrBasisSize)%Energy = Glob_CurrEnergy
              Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt = &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              ! Sort='yes' needs the workspace that is allocated only when
              ! the block was actually moved, i.e. when FinalFunc/=nfa.
              IF (FinalFunc == nfa) THEN
                CALL SaveResults(Sort='no')
              ELSE
                CALL SaveResults(Sort='yes')
              ENDIF
              WRITE(*, *) 'Data file has been updated'
              CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
              IF (Verbose >= 1) WRITE(*, *) 'Some current statistics:'
              IF (Glob_OverlapPenaltyAllowed) THEN
                IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
                IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
              ELSE
                WRITE(*, *) 'Energy                            ', Evalue
              ENDIF
              IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
              IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
              TimeOfLastSave = Glob_TimeSinceStart
            ENDIF

            IF ((Glob_TimeSinceStart-TimeOfLastHessSave > HessianSaveMinTimeInterv) &
                .AND. (SaveHessian)) THEN
              ! Save the Hessian if more than HessianSaveMinTimeInterv seconds
              ! have passed since the last save
              CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
              TimeOfLastHessSave = Glob_TimeSinceStart
            ENDIF

          ENDIF
        ENDIF

      CASE (2)  ! Only gradient is needed
        CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        DO i = 1, nfo
          Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
        ENDDO
        CALL EnergyIB(Evalue, grad, .TRUE., ErrCode)
        NumOfGradEval = NumOfGradEval+1
        IF (ErrCode /= 0) THEN
          NumOfFailures = NumOfFailures+1
          IV(2) = 1
          ! Restore the last good eigenvector as the vector left by the failed
          ! inverse iteration process may be unusable as a starting vector
          Glob_LastEigvector(1:nfa) = v_good(1:nfa)
        ENDIF
        ! The rejection sentinel (wrong state, ErrCode = 0) is reported to DRMNG like a failure
        IF ((ErrCode == 0) .AND. (ABS(Evalue) > 1.0E30_wp)) IV(2) = 1

      CASE (3:8)  ! Some kind of convergence has been reached
        ExitNeeded = .TRUE.

      CASE (9:10)  ! DRMNG reached its evaluation (or iteration) limit
        ! IV(17) holds the remaining energy-evaluation budget of this basis
        ! size, so this is the normal end of a step that exhausts it. x may
        ! still sit on a rejected trial point; the code after the loop goes
        ! back to the accepted point when that is the case.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, *) 'Warning WC0135 in FullOpt1I: number of energy evaluations reached limit'
          IF (Verbose >= 1) WRITE(*, '(1x,a,i0,a,i0,a)') '(', NumOfEnergyEval, ' in this step, ', &
            NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval, ' in total)'
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.

      CASE DEFAULT
        ! DRMNG returns 63 or 65 when it gives up on an uncomputable
        ! value, and >=14 for a bad input. None of those match a case
        ! above, so without this the loop would call DRMNG again for
        ! ever. Leave with the best point found so far.
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0)') &
            'Warning WC0141 in FullOpt1I: DRMNG returned IV(1) =', IV(1)
          IF (Verbose >= 1) WRITE(*, *) 'Optimization is terminated'
        ENDIF
        ExitNeeded = .TRUE.

      ENDSELECT

      IF (NumOfFailures > Glob_MaxEnergyFailsAllowed) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,a,1x,a)') &
            'Error EC0166 in FullOpt1I: number of failures in energy or gradient', &
            'calculations during the optimization of nonlinear parameters', &
            'exceeded limit'
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF

    ENDDO


    !==================================================================
    ! Final energy at the best point found
    !==================================================================
    ! DRMNG can leave x on a rejected trial point of the last line search
    ! (evaluation limit, false convergence). Its accepted iterate x0 is in
    ! V(IV(43):IV(43)+nv-1) with f0 = V(13); V(10) is the last f evaluated.
    ! Whenever the point in x is not better than x0, go back to x0, so
    ! that the final energy and the data file describe the best point.
    IF (Glob_ProcID == 0) THEN
      IF ((IV(31) >= 1) .AND. (V(10) >= V(13))) THEN
        IF ((V(10) > V(13)) .AND. (Verbose >= 2)) &
          WRITE(*, *) 'The last trial point is discarded: back to the last accepted point'
        x(1:nv) = V(IV(43):IV(43)+nv-1)
      ENDIF
    ENDIF
    CALL MPI_BCAST(x, nv, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    DO i = 1, nfo
      Glob_NonlinParam(1:npt, InitFuncNew+i-1) = x((i-1)*npt+1:i*npt)
    ENDDO
    Evalue = EnergyIA(InitFuncNew, nfa, .TRUE., ErrCode)
    IF (ErrCode /= 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0167 in FullOpt1I: failed to evaluate energy after the optimization'
        IF (Verbose >= 1) WRITE(*, *) 'of nonlinear parameters'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! Glob_CurrEnergy is the physical energy, without the penalty
    IF (Glob_OverlapPenaltyAllowed) THEN
      Glob_CurrEnergy = Evalue-Glob_TotalOverlapPenalty
    ELSE
      Glob_CurrEnergy = Evalue
    ENDIF

    ! Printing the number of energy/gradient evaluations
    ! and the energy after the optimization
    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 1) WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, '(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations', NumOfEnergyEval, '/', NumOfGradEval
      WRITE(*, *) 'Final energy and overlap statistics:'
      CALL GetOverlapStatistics(InitFuncNew, nfa, MaxAbsOverlap, MinAbsOverlap, AverageAbsOverlap)
      IF (Glob_OverlapPenaltyAllowed) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Energy (without overlap penalty)  ', Evalue-Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Overlap penalty                   ', Glob_TotalOverlapPenalty
        IF (Verbose >= 1) WRITE(*, *) 'Energy (including overlap penalty)', Evalue
      ELSE
        WRITE(*, *) 'Energy                            ', Evalue
      ENDIF
      IF (Verbose >= 1) WRITE(*, *) 'Maximal overlap                   ', MaxAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Minimal overlap                   ', MinAbsOverlap
      IF (Verbose >= 1) WRITE(*, *) 'Average abs value of overlap      ', AverageAbsOverlap
    ENDIF

    ! Adding data to history
    Glob_History(Glob_CurrBasisSize)%Energy = Glob_CurrEnergy
    Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt = &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval


    !==================================================================
    ! Move the block back and release everything
    !==================================================================
    ! Shifting basis functions InitFunc through FinalFunc back from
    ! the very end to the middle, where they initially were, and, if
    ! necessary, doing the permutation of matrix elements to reflect
    ! this change in function order.
    IF (FinalFunc /= nfa) THEN
      ! allocate space
      tas = nfa-FinalFunc
      ALLOCATE(NonlinParamTemp(1:npt, 1:tas))
      ALLOCATE(FuncNumTemp(1:tas))
      ALLOCATE(TempR(1:tas))
      CALL PermuteFunctions(InitFunc, InitFunc+tas-1, FuncNumTemp, NonlinParamTemp)
      IF (Glob_UseSwapFile) CALL PermuteMatrixElements(InitFunc, InitFunc+tas-1, TempR)
      DEALLOCATE(TempR)
      DEALLOCATE(FuncNumTemp)
      DEALLOCATE(NonlinParamTemp)
      ! deallocate workspace for SaveResults
      DEALLOCATE(Glob_IntWorkArrForSaveResults)
    ENDIF

    ! The shift refresh in EnergyIB is for this optimization only
    WrkP_RefreshShiftInIB = .FALSE.

    ! saving results; the Hessian goes to disk with them, so a restart
    ! finds the file even when no evaluation improved the energy
    IF (Glob_ProcID == 0) THEN
      CALL SaveResults(Sort='no')
      IF (SaveHessian) CALL SaveHessianFile(V, IVLMAT, D, nv, HessFileName, IsHessSaveSuccess)
    ENDIF

    CALL StoreMatricesInSwapFile()

    IF (Glob_OverlapPenaltyAllowed) Glob_OverlapPenaltyAllowed = .FALSE.

    ! deallocate workspace
    DEALLOCATE(grad)
    DEALLOCATE(x)

    ! deallocate arrays used by DRMNG
    IF (Glob_ProcID == 0) THEN
      DEALLOCATE(V)
      DEALLOCATE(D)
    ENDIF

    ! deallocate workspace for EnergyIB
    DEALLOCATE(Glob_WkGR)

    ! Deallocate workspace for GSEPIIS
    DEALLOCATE(Glob_LastEigvector)
    DEALLOCATE(v_good)
    DEALLOCATE(Glob_WorkForGSEPIIS)

    ! Deallocate some global arrays
    DEALLOCATE(Glob_DlBuff2)
    DEALLOCATE(Glob_DlBuff1)
    DEALLOCATE(Glob_DkBuff2)
    DEALLOCATE(Glob_DkBuff1)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_D)
    DEALLOCATE(Glob_invD)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (WrkP_WrongStateCount > 0)) THEN
      IF (Verbose >= 1) WRITE(*, *) 'Trial points refused because inverse iteration landed on'
      WRITE(*, *) 'a level other than WHICH_EIGENVALUE =', Glob_WhichEigenvalue, &
                 ' :', WrkP_WrongStateCount
    ENDIF
    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine FullOpt1I finished'


  END SUBROUTINE FullOpt1I


  SUBROUTINE EliminateLittleContribFunc(LinCoeffThreshold, FileName, PrintInfoSpec, GSEPSolMethod)
    ! Subroutine EliminateLittleContribFunc eliminates basis
    ! functions whose contribution to the energy is small. More
    ! precisely, it eliminates functions that have linear coefficients
    ! whose absolute values are smaller than LinCoeffThreshold. It is
    ! important to note that this subroutine uses coefficients in
    ! front of normalized functions (so that Glob_S is the overlap
    ! matrix of normalized functions, with Glob_S(i,i)=1).
    ! The result is stored in file whose name is defined by parameter
    ! FileName. After saving the results the subroutine terminates
    ! the program.
    ! Parameter PrintInfoSpec specifies what information should be
    ! shown regarding the linear coefficients:
    ! PrintInfoSpec=0,1 : the subroutine does not print any specific info
    !                   regarding linear coefficients.
    ! PrintInfoSpec=2   : the subroutine prints linear coefficients
    !                     of all functions
    ! GSEPSolMethod     : optional solver selection. It defaults to G so existing
    !                     callers retain their behavior; main passes Q explicitly.

    ! Arguments:
    REAL(wp), INTENT(IN)                       :: LinCoeffThreshold
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName
    INTEGER, INTENT(IN)                        :: PrintInfoSpec
    CHARACTER(1), INTENT(IN), OPTIONAL         :: GSEPSolMethod

    ! Local variables:
    INTEGER                                :: i, j
    INTEGER                                :: np, npt, cbs
    INTEGER                                :: OpenFileErr, ErrorCode
    LOGICAL                                :: IsSwapFileOK
    INTEGER                                :: BlockSizeForDSYGVX
    REAL(wp)                               :: Evalue
    REAL(wp)                               :: Min_c, Max_c
    REAL(wp)                               :: Aver_c
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: NonlinParamTemp
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: MaskArray
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: ZIndTemp
    CHARACTER(Glob_FileNameLength)         :: ch_temp
    CHARACTER(1)                           :: Method

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine EliminateLittleContribFunc started'
      WRITE(*, *) 'Linear coefficient threshold value is', LinCoeffThreshold
    ENDIF

    ! Setting the values of some global variables
    Method = 'G'
    IF (PRESENT(GSEPSolMethod)) Method = GSEPSolMethod
    Glob_GSEPSolutionMethod = Method
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize

    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (Method == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (Method == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Reading data from swap file
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    ENDIF

    CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0170 in EliminateLittleContribFunc: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *) 'Basis size before elimination', cbs
      WRITE(*, *) 'Energy before elimination    ', Evalue
    ENDIF

    ALLOCATE(NonlinParamTemp(1:npt, cbs))
    ALLOCATE(MaskArray(1:cbs))
    ALLOCATE(ZIndTemp(cbs))
    MaskArray = 0

    IF ((PrintInfoSpec > 1) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, *) 'List of all linear coefficients:'
      DO i = 1, cbs
        WRITE(*, '(i6,a4,f19.12)') i, '  c=', Glob_c(i)
      ENDDO
      WRITE(*, *)
    ENDIF

    Min_c = Glob_c(1)
    Max_c = Glob_c(1)
    Aver_c = ZERO
    j = 0
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) >= LinCoeffThreshold) THEN
        NonlinParamTemp(1:npt, i-j) = Glob_NonlinParam(1:npt, i)
        ZIndTemp(i-j) = Glob_PWR(i)
      ELSE
        IF ((j == 0) .AND. (Glob_ProcID == 0)) WRITE(*, *) 'Little contributing function list:'
        j = j+1
        MaskArray(i) = 1
        IF (Glob_ProcID == 0) &
          WRITE(*, '(i6,a1,i6,a4,f19.12)') j, ':', i, '  c=', Glob_c(i)
      ENDIF
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Aver_c = Aver_c+ABS(Glob_c(i))
    ENDDO
    Aver_c = Aver_c/cbs

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      WRITE(*, '(1x,a41,e16.9)') 'Smallest by magnitude linear coeff.     =', Min_c
      WRITE(*, '(1x,a41,e16.9)') 'Largest by magnitude linear coeff.      =', Max_c
      WRITE(*, '(1x,a41,e16.9)') 'Average absolute value of linear coeff. =', Aver_c
      WRITE(*, *)
    ENDIF

    IF (j == 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'There are no functions with the contribution'
        WRITE(*, *) 'smaller than ', LinCoeffThreshold
        WRITE(*, *) 'No output file have been written. Program will now stop'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! An empty basis has no generalized eigenproblem and qrlinalg intentionally
    ! has no valid order-zero state. Fail before changing the published basis
    ! size or overwriting the input/output data with an unusable result.
    IF (j == cbs) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'All basis functions are below the coefficient threshold'
        WRITE(*, *) 'No output file has been written. Program will now stop'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    ENDIF

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *) 'Basis size before elimination', cbs
      WRITE(*, *) 'Energy before elimination    ', Evalue
    ENDIF

    Glob_NonlinParam(1:npt, 1:cbs-j) = NonlinParamTemp(1:npt, 1:cbs-j)
    Glob_PWR(1:cbs-j) = ZIndTemp(1:cbs-j)

    IF (Glob_ProcID == 0) THEN
      IF (Method == 'Q') THEN
        WRITE(*, *) 'Deleting selected rows and columns from the QR factorization...'
      ELSE
        WRITE(*, *) 'Computing matrix elements and solving eigenvalue problem with the'
        WRITE(*, *) 'basis where little contributing functions are eliminated...'
      ENDIF
    ENDIF
    IF (Method == 'Q') THEN
      ! Every survivor-survivor matrix element is unchanged. Preserve its
      ! canonical value and delete the matching factor rows and columns instead
      ! of evaluating the complete smaller matrix for a second time.
      CALL DeleteQMaskedFunctions(MaskArray, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) THEN
        Glob_CurrBasisSize = cbs-j
        cbs = Glob_CurrBasisSize
        CALL SolveQ(Evalue, ErrorCode)
      ENDIF
    ELSE
      Glob_CurrBasisSize = cbs-j
      cbs = Glob_CurrBasisSize
      CALL ComputeMatElem(1, cbs)
      CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    ENDIF
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0171 in EliminateLittleContribFunc: energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *) 'Basis size after elimination ', cbs
      WRITE(*, *) 'Energy after elimination     ', Evalue
    ENDIF

    DO i = 1, cbs
      Glob_History(i)%Energy = ZERO
      Glob_History(i)%CyclesDone = 0
      Glob_History(i)%InitFuncAtLastStep = 0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = 0
    ENDDO
    Glob_History(cbs)%Energy = Evalue

    Glob_CurrEnergy = Evalue
    Glob_LastEigvalTol = 1.0E+35_wp
    Glob_BestEigvalTol = 1.0E+35_wp
    Glob_WorstEigvalTol = 1.0E-35_wp

    ch_temp = Glob_DataFileName
    Glob_DataFileName = FileName
    IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')
    Glob_DataFileName = ch_temp

    DEALLOCATE(MaskArray)
    DEALLOCATE(ZIndTemp)
    DEALLOCATE(NonlinParamTemp)

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    IF (Method == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    ! Deallocate workspace for DSYGVX
    IF (Method == 'G') THEN
      DEALLOCATE(Glob_IWorkForDSYGVX)
      DEALLOCATE(Glob_WorkForDSYGVX)
    ELSE
      CALL ClearQWorkspace()
    ENDIF

    IF (Glob_ProcID == 0) THEN
      i = LEN_TRIM(FileName)
      WRITE(*, *) 'New basis has been saved in file', FileName(1:i)
      WRITE(*, *) 'Program will now stop'
    ENDIF

    CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

  END SUBROUTINE EliminateLittleContribFunc

  SUBROUTINE EliminateLinDepFunc(LinDepThreshold, FileName, PrintInfoSpec, GSEPSolMethod)
    ! Subroutine EliminateLinDepFunc eliminates linearly dependent
    ! functions. It checks for pair linear dependency only. It
    ! removes those functions from the basis whose overlap (absolute value)
    ! with any of other basis function with a smaller number is greater
    ! than LinDepThreshold. It is important to note that this subroutine
    ! uses normalized functions (so that Glob_S is the overlap
    ! matrix of normalized functions, with Glob_S(i,i)=1).
    ! The result is stored in a file whose name is defined by parameter
    ! FileName. After saving the results the program is terminated.
    ! Parameter PrintInfoSpec specifies what information should be
    ! printed during linear dependency check:
    ! PrintInfoSpec=0 : the subroutine does not print any info
    !                   regarding basis functions that are linearly dependent.
    ! PrintInfoSpec=1 : the subroutine prints the overlap values of
    !                   linearly dependent functions.
    ! PrintInfoSpec=2 : same as the previous case, but in addition it also prints the
    !                   nonlinear parameters of linearly dependent functions.
    ! GSEPSolMethod   : optional solver selection. It defaults to G so existing
    !                   callers retain their behavior; main passes Q explicitly.

    ! Arguments:
    REAL(wp), INTENT(IN)                       :: LinDepThreshold
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName
    INTEGER, INTENT(IN)                        :: PrintInfoSpec
    CHARACTER(1), INTENT(IN), OPTIONAL         :: GSEPSolMethod

    ! Local variables:
    INTEGER                            :: i, j, k
    INTEGER                            :: np, npt, cbs
    INTEGER                            :: OpenFileErr, ErrorCode
    LOGICAL                            :: IsSwapFileOK
    INTEGER                            :: BlockSizeForDSYGVX
    REAL(wp)                           :: Evalue
    REAL(wp)                           :: MaxOverlap, MinOverlap
    REAL(wp)                           :: AverOverlap
    REAL(wp)                           :: Min_c, Max_c
    REAL(wp)                           :: Average_c
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MaskArray
    CHARACTER(Glob_FileNameLength)     :: ch_temp
    CHARACTER(1)                       :: Method

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine EliminateLinDepFunc started'
      WRITE(*, *) 'Overlap threshold value is', LinDepThreshold
    ENDIF

    ! Setting the values of some global variables
    Method = 'G'
    IF (PRESENT(GSEPSolMethod)) Method = GSEPSolMethod
    Glob_GSEPSolutionMethod = Method
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize

    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (Method == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (Method == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Allocate local workspace
    ALLOCATE(MaskArray(1:cbs))

    ! Reading data from swap file

    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'
    ENDIF

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0175 in EliminateLinDepFunc: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *) 'Pair linear dependency check:'
      WRITE(*, *)
    ENDIF
    ! Check overlap
    MaskArray(1:cbs) = 0
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > LinDepThreshold) THEN
          MaskArray(j) = MaskArray(j)+1
          k = k+1
          IF (Glob_ProcID == 0) THEN
            IF (PrintInfoSpec > 0) WRITE(*, '(i6,a1,i6,i6,a5,f17.14)') k, ':', i, j, '   S=', Glob_S(j, i)
            IF (PrintInfoSpec == 2) THEN
              WRITE(*, '(1x,i6,1x,i6)', ADVANCE='no') i, Glob_PWR(i)
              CALL writerealarradv(6, Glob_NonlinParam(1:Glob_npt, i), Glob_npt)
              WRITE(*, *) '      c=', Glob_c(i)
              WRITE(*, '(1x,i6,1x,i6)', ADVANCE='no') j, Glob_PWR(j)
              CALL writerealarradv(6, Glob_NonlinParam(1:Glob_npt, j), Glob_npt)
              WRITE(*, *) '      c=', Glob_c(j)
            ENDIF
            WRITE(*, *)
          ENDIF
        ENDIF
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      ! An order-one basis has no off-diagonal pair statistics.
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients:
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      WRITE(*, *) '========== Summary before elimination: =========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size before elimination =', cbs
      WRITE(*, *) 'Energy before elimination     =', Evalue
      WRITE(*, *) '================================================'
      WRITE(*, *)
    ENDIF

    IF (k == 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'No linearly dependent functions have been found'
        WRITE(*, *) 'No file have been written. Program will now stop'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    j = 0
    i = 1
    DO WHILE (i+j <= cbs)
      IF (MaskArray(i+j) > 0) THEN
        j = j+1
      ELSE
        Glob_NonlinParam(1:npt, i) = Glob_NonlinParam(1:npt, i+j)
        Glob_PWR(i) = Glob_PWR(i+j)
        i = i+1
      ENDIF
    ENDDO

    ! k counts offending pairs for the diagnostic list. One later function can
    ! overlap several earlier functions, but MaskArray removes that function
    ! only once. j is the number of unique masked functions actually skipped by
    ! the compaction loop and therefore defines the new basis order.
    IF (Method == 'Q') THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a)', ADVANCE='no') &
        'Deleting selected rows and columns from QR factors...'
      CALL DeleteQMaskedFunctions(MaskArray, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) THEN
        Glob_CurrBasisSize = cbs-j
        cbs = Glob_CurrBasisSize
        CALL SolveQ(Evalue, ErrorCode)
      ENDIF
    ELSE
      Glob_CurrBasisSize = cbs-j
      cbs = Glob_CurrBasisSize
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) ' done'
        IF (Verbose >= 2) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      ENDIF
      CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    ENDIF
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0176 in EliminateLinDepFunc: energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! Check overlap
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      ! Elimination can legitimately leave one surviving basis function.
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients:
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *)
      WRITE(*, *) '========== Summary after elimination: =========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size after elimination  =', cbs
      WRITE(*, *) 'Energy after elimination      =', Evalue
      WRITE(*, *) '================================================'
      WRITE(*, *)
    ENDIF

    DO i = 1, cbs
      Glob_History(i)%Energy = ZERO
      Glob_History(i)%CyclesDone = 0
      Glob_History(i)%InitFuncAtLastStep = 0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = 0
    ENDDO
    Glob_History(cbs)%Energy = Evalue

    Glob_CurrEnergy = Evalue
    Glob_LastEigvalTol = 1.0E+35_wp
    Glob_BestEigvalTol = 1.0E+35_wp
    Glob_WorstEigvalTol = 1.0E-35_wp

    ch_temp = Glob_DataFileName
    Glob_DataFileName = FileName
    IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')
    Glob_DataFileName = ch_temp

    ! deallocate local workspace
    DEALLOCATE(MaskArray)

    ! Deallocate workspace for DSYGVX
    IF (Method == 'G') THEN
      DEALLOCATE(Glob_IWorkForDSYGVX)
      DEALLOCATE(Glob_WorkForDSYGVX)
    ELSE
      CALL ClearQWorkspace()
    ENDIF

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    IF (Method == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF (Glob_ProcID == 0) THEN
      i = LEN_TRIM(FileName)
      WRITE(*, *) 'New basis has been saved in file', FileName(1:i)
      WRITE(*, *) 'Program will now stop'
    ENDIF

    CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

  END SUBROUTINE EliminateLinDepFunc

  SUBROUTINE SeparateLinDepFunc(LinDepThreshold, SeparationParam, FileName, PrintInfoSpec, GSEPSolMethod)
    ! Subroutine SeparateLinDepFunc does exactly the same thing as
    ! subroutine EliminateLinDepFunc does, but without throwing away
    ! linearly dependent functions. Instead, it changes the parameters
    ! of such functions randomly (the random shift is controlled by
    ! argument SeparationParam, so that a_new lies within
    ! interval [a_old*(1-SeparationParam),a_old*(1+SeparationParam)]).
    ! Optional GSEPSolMethod defaults to G for source compatibility. The Q path
    ! keeps the basis canonical and consumes only the lower H/S triangles.

    ! Arguments:
    REAL(wp), INTENT(IN)                       :: LinDepThreshold
    REAL(wp), INTENT(IN)                       :: SeparationParam
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName
    INTEGER, INTENT(IN)                        :: PrintInfoSpec
    CHARACTER(1), INTENT(IN), OPTIONAL         :: GSEPSolMethod

    ! Local variables:
    INTEGER                            :: i, j, k, NumActive
    INTEGER                            :: np, npt, cbs
    INTEGER                            :: OpenFileErr, ErrorCode
    LOGICAL                            :: IsSwapFileOK
    INTEGER                            :: BlockSizeForDSYGVX
    REAL(wp)                           :: Evalue, r
    REAL(8)                            :: r8
    REAL(wp)                           :: MaxOverlap, MinOverlap
    REAL(wp)                           :: AverOverlap
    REAL(wp)                           :: Min_c, Max_c
    REAL(wp)                           :: Average_c
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MaskArray, ActiveFunction
    CHARACTER(Glob_FileNameLength)     :: ch_temp
    CHARACTER(1)                       :: Method

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine SeparateLinDepFunc started'
      WRITE(*, *) 'Overlap threshold value is', LinDepThreshold
      WRITE(*, *) 'Separation paramameter is ', SeparationParam
    ENDIF

    ! Setting the values of some global variables
    Method = 'G'
    IF (PRESENT(GSEPSolMethod)) Method = GSEPSolMethod
    Glob_GSEPSolutionMethod = Method
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize

    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (Method == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (Method == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Allocate local workspace
    ALLOCATE(MaskArray(1:cbs))

    ! Reading data from swap file

    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'
    ENDIF

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0180 in SeparateLinDepFunc: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *) 'Pair linear dependency check:'
      WRITE(*, *)
    ENDIF
    ! Check overlap
    MaskArray(1:cbs) = 0
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > LinDepThreshold) THEN
          MaskArray(j) = MaskArray(j)+1
          k = k+1
          IF (Glob_ProcID == 0) THEN
            IF (PrintInfoSpec > 0) WRITE(*, '(i6,a1,i6,i6,a5,f17.14)') k, ':', i, j, '   S=', Glob_S(j, i)
            IF (PrintInfoSpec == 2) THEN
              WRITE(*, '(1x,i6,1x,i6)', ADVANCE='no') i, Glob_PWR(i)
              CALL writerealarradv(6, Glob_NonlinParam(1:Glob_npt, i), Glob_npt)
              WRITE(*, *) '      c=', Glob_c(i)
              WRITE(*, '(1x,i6,1x,i6)', ADVANCE='no') j, Glob_PWR(j)
              CALL writerealarradv(6, Glob_NonlinParam(1:Glob_npt, j), Glob_npt)
              WRITE(*, *) '      c=', Glob_c(j)
            ENDIF
            WRITE(*, *)
          ENDIF
        ENDIF
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients:
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      WRITE(*, *) '========== Summary before separation: =========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size before separation =', cbs
      WRITE(*, *) 'Energy before separation     =', Evalue
      WRITE(*, *) '================================================'
      WRITE(*, *)
    ENDIF

    IF (k == 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'No linearly dependent functions have been found'
        WRITE(*, *) 'No file have been written. Program will now stop'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Method == 'Q') THEN
      ! The initial cleanup workspace holds one transaction column. Expand that
      ! storage only to the number of unique functions selected by the overlap
      ! mask, while retaining the already computed full QR factorization.
      NumActive = COUNT(MaskArray > 0)
      CALL EnsureQActiveCapacity(NumActive, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) THEN
        ALLOCATE(ActiveFunction(NumActive))
        j = 0
        DO i = 1, cbs
          IF (MaskArray(i) > 0) THEN
            j = j+1
            ActiveFunction(j) = i
          ENDIF
        ENDDO
        CALL SetQActiveFunctions(ActiveFunction, ErrorCode)
      ENDIF
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL CaptureQMatrixParameters(ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) THEN
        IF (Glob_ProcID == 0) WRITE(*, *) &
          'Error EC0181 in SeparateLinDepFunc: Q transaction cannot be prepared'
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
      ENDIF
    ENDIF

    j = 0
    i = 1
    DO i = 1, cbs
      IF (MaskArray(i) > 0) THEN
        DO j = 1, npt
          CALL RANDOM_NUMBER(r8)
          r = TWO*(r8-ONEHALF)*SeparationParam
          Glob_NonlinParam(j, i) = Glob_NonlinParam(j, i)*(1+r)
        ENDDO
      ENDIF
    ENDDO
    CALL MPI_BCAST(Glob_NonlinParam, cbs*npt, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (Method == 'Q') THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a)', ADVANCE='no') &
        'Computing selected matrix columns and updating QR...'
      CALL AssembleQTrial(.FALSE., ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL ApplyQTrial(ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL SolveQ(Evalue, ErrorCode)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) ' done'
        IF (Verbose >= 2) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      ENDIF
      CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    ENDIF
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0181 in EliminateLinDepFunc: energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! Check overlap
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients:
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *)
      WRITE(*, *) '========== Summary after separation: ==========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size after separation  =', cbs
      WRITE(*, *) 'Energy after separation      =', Evalue
      WRITE(*, *) '================================================'
      WRITE(*, *)
    ENDIF

    DO i = 1, cbs
      Glob_History(i)%Energy = ZERO
      Glob_History(i)%CyclesDone = 0
      Glob_History(i)%InitFuncAtLastStep = 0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = 0
    ENDDO
    Glob_History(cbs)%Energy = Evalue

    Glob_CurrEnergy = Evalue
    Glob_LastEigvalTol = 1.0E+35_wp
    Glob_BestEigvalTol = 1.0E+35_wp
    Glob_WorstEigvalTol = 1.0E-35_wp

    ch_temp = Glob_DataFileName
    Glob_DataFileName = FileName
    IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')
    Glob_DataFileName = ch_temp

    ! deallocate local workspace
    IF (ALLOCATED(ActiveFunction)) DEALLOCATE(ActiveFunction)
    DEALLOCATE(MaskArray)

    ! Deallocate workspace for DSYGVX
    IF (Method == 'G') THEN
      DEALLOCATE(Glob_IWorkForDSYGVX)
      DEALLOCATE(Glob_WorkForDSYGVX)
    ELSE
      CALL ClearQWorkspace()
    ENDIF

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    IF (Method == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF (Glob_ProcID == 0) THEN
      i = LEN_TRIM(FileName)
      WRITE(*, *) 'New basis has been saved in file', FileName(1:i)
      WRITE(*, *) 'Program will now stop'
    ENDIF

    CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

  END SUBROUTINE SeparateLinDepFunc

  SUBROUTINE SeparateFuncLargeCoeff(LCThreshold, SeparationParam, FileName, PrintInfoSpec, GSEPSolMethod)
    ! Subroutine SeparateFuncLargeCoeff changes linear parameters of
    ! those basis functions whose linear coefficients (more precisely their
    ! absolute values) exceed LCThreshold. It is important to note that this
    ! subroutine uses normalized functions (so that Glob_S is the overlap
    ! matrix of normalized functions, with Glob_S(i,i)=1). The change
    ! in nonlinear parameters is controlled by argument SeparationParam.
    ! Basically the nonlinear parameters of bad functions are chosen
    ! randomly from the interval [a_old*(1-SeparationParam), a_old*(1+SeparationParam)].
    ! After saving the results the program is terminated.
    ! Parameter PrintInfoSpec specifies what information should be
    ! shown:
    ! PrintInfoSpec=0 : the subroutine does not print any info
    !                   about basis functions.
    ! PrintInfoSpec=1 : the subroutine prints the values of the linear and
    !                   nonlinear parameters of bad functions.
    ! PrintInfoSpec=2 : same as the previous case, but in addition it also prints
    !                   the linear parameters of all basis functions.
    ! GSEPSolMethod   : optional solver selection. It defaults to G so existing
    !                   callers retain their behavior; main passes Q explicitly.

    ! Arguments:
    REAL(wp), INTENT(IN)                       :: LCThreshold
    REAL(wp), INTENT(IN)                       :: SeparationParam
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName
    INTEGER, INTENT(IN)                        :: PrintInfoSpec
    CHARACTER(1), INTENT(IN), OPTIONAL         :: GSEPSolMethod

    ! Local variables:
    INTEGER                            :: i, j, k, NumActive
    INTEGER                            :: np, npt, cbs
    INTEGER                            :: OpenFileErr, ErrorCode
    LOGICAL                            :: IsSwapFileOK
    INTEGER                            :: BlockSizeForDSYGVX
    REAL(wp)                           :: Evalue, r
    REAL(8)                            :: r8
    REAL(wp)                           :: MaxOverlap, MinOverlap
    REAL(wp)                           :: AverOverlap
    REAL(wp)                           :: Min_c, Max_c
    REAL(wp)                           :: Average_c
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ActiveFunction
    CHARACTER(Glob_FileNameLength)     :: ch_temp
    CHARACTER(1)                       :: Method

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine SeparateFuncLargeCoeff started'
      WRITE(*, *) 'Linear coefficient threshold value is', LCThreshold
      WRITE(*, *) 'Separation parameter is              ', SeparationParam
    ENDIF

    ! Setting the values of some global variables
    Method = 'G'
    IF (PRESENT(GSEPSolMethod)) Method = GSEPSolMethod
    Glob_GSEPSolutionMethod = Method
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    np = Glob_np
    npt = Glob_npt
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize

    ! Allocate some global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (Method == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (Method == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Reading data from swap file
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) ' done'
    ENDIF

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
    CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) &
        'Error EC0185 in SeparateFuncLargeCoeff: initial energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *)
      IF (PrintInfoSpec >= 2) THEN
        WRITE(*, *) 'Linear coefficients of all basis functions before separation:'
        DO i = 1, cbs
          WRITE(*, '(1x,i6,a4,f19.12)') i, '  c=', Glob_c(i)
        ENDDO
      ENDIF
      WRITE(*, *)
    ENDIF

    ! Check overlap
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients
    k = 0
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > LCThreshold) THEN
        k = k+1
        IF ((Glob_ProcID == 0) .AND. (PrintInfoSpec >= 1)) THEN
          IF (k == 1) WRITE(*, *) 'Functions whose linear coefficients exceed threshold'
          WRITE(*, '(1x,i5,a3,i6,a4,f19.12)') k, ':  ', i, '  c=', Glob_c(i)
          WRITE(*, *) Glob_NonlinParam(1:Glob_npt, i)
        ENDIF
      ENDIF
    ENDDO

    IF (k == 0) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'There are no functions whose linear coefficients exceed threshold'
        WRITE(*, *) 'No file have been written. Program will now stop'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (Method == 'Q') THEN
      ! The coefficient scan already provides the exact active count. Reserve
      ! only those transaction columns and preserve the initial factorization.
      NumActive = k
      CALL EnsureQActiveCapacity(NumActive, ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) THEN
        ALLOCATE(ActiveFunction(NumActive))
        k = 0
        DO i = 1, cbs
          IF (ABS(Glob_c(i)) > LCThreshold) THEN
            k = k+1
            ActiveFunction(k) = i
          ENDIF
        ENDDO
        CALL SetQActiveFunctions(ActiveFunction, ErrorCode)
      ENDIF
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL CaptureQMatrixParameters(ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) THEN
        IF (Glob_ProcID == 0) WRITE(*, *) &
          'Error EC0186 in SeparateFuncLargeCoeff: Q transaction cannot be prepared'
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
      ENDIF
    ENDIF
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *) '========== Summary before separation: =========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size before separation =', cbs
      WRITE(*, *) 'Energy before separation     =', Evalue
      WRITE(*, *) '================================================'
    ENDIF

    ! now change the nonlinear parameters of bad functions
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > LCThreshold) THEN
        DO j = 1, npt
          CALL RANDOM_NUMBER(r8)
          r = (r8-ONEHALF)*2*SeparationParam
          Glob_NonlinParam(j, i) = Glob_NonlinParam(j, i)*(1+r)
        ENDDO
      ENDIF
    ENDDO
    CALL MPI_BCAST(Glob_NonlinParam, cbs*Glob_npt, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (Method == 'Q') THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a)', ADVANCE='no') &
        'Computing selected matrix columns and updating QR...'
      CALL AssembleQTrial(.FALSE., ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL ApplyQTrial(ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL SolveQ(Evalue, ErrorCode)
    ELSE
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a28)', ADVANCE='no') 'Computing matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) ' done'
        IF (Verbose >= 2) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      ENDIF
      CALL SolveEliminationGSEP(Method, cbs, Evalue, ErrorCode)
    ENDIF
    IF (ErrorCode /= 0) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'Error EC0186 in SeparateFuncLargeCoeff: energy cannot be computed'
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    ! Check overlap
    MaxOverlap = ZERO
    MinOverlap = HUGE(MinOverlap)/2
    AverOverlap = ZERO
    k = 0
    DO i = 1, cbs
      DO j = i+1, cbs
        IF (ABS(Glob_S(j, i)) > ABS(MaxOverlap)) MaxOverlap = Glob_S(j, i)
        IF (ABS(Glob_S(j, i)) < ABS(MinOverlap)) MinOverlap = Glob_S(j, i)
        AverOverlap = AverOverlap+ABS(Glob_S(j, i))
      ENDDO
    ENDDO
    IF (cbs > 1) THEN
      AverOverlap = AverOverlap/(cbs*(cbs-1)/TWO)
    ELSE
      MaxOverlap = ZERO
      MinOverlap = ZERO
      AverOverlap = ZERO
    ENDIF

    ! Check linear coefficients
    Min_c = HUGE(Min_c)/2
    Max_c = ZERO
    Average_c = ZERO
    DO i = 1, cbs
      IF (ABS(Glob_c(i)) > ABS(Max_c)) Max_c = Glob_c(i)
      IF (ABS(Glob_c(i)) < ABS(Min_c)) Min_c = Glob_c(i)
      Average_c = Average_c+ABS(Glob_c(i))
    ENDDO
    Average_c = Average_c/cbs

    IF (Glob_ProcID == 0) THEN
      IF (Verbose >= 2) WRITE(*, *) ' done'
      WRITE(*, *) '========== Summary after separation: ==========='
      WRITE(*, *) 'Maximal by magnitude overlap         =', MaxOverlap
      WRITE(*, *) 'Minimal by magnitude overlap         =', MinOverlap
      WRITE(*, *) 'Average absolute value of overlap    =', AverOverlap
      WRITE(*, *) 'Maximal by magnitude lin. coeff.     =', Max_c
      WRITE(*, *) 'Minimal by magnitude lin. coeff.     =', Min_c
      WRITE(*, *) 'Average absolute value of lin. coeff.=', Average_c
      WRITE(*, *) 'Basis size after separation', cbs
      WRITE(*, *) 'Energy after separation    ', Evalue
      WRITE(*, *) '================================================'
    ENDIF

    DO i = 1, cbs
      Glob_History(i)%Energy = ZERO
      Glob_History(i)%CyclesDone = 0
      Glob_History(i)%InitFuncAtLastStep = 0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt = 0
    ENDDO
    Glob_History(cbs)%Energy = Evalue

    Glob_CurrEnergy = Evalue
    Glob_LastEigvalTol = 1.0E+35_wp
    Glob_BestEigvalTol = 1.0E+35_wp
    Glob_WorstEigvalTol = 1.0E-35_wp

    ch_temp = Glob_DataFileName
    Glob_DataFileName = FileName
    IF (Glob_ProcID == 0) CALL SaveResults(Sort='no')
    Glob_DataFileName = ch_temp

    IF (ALLOCATED(ActiveFunction)) DEALLOCATE(ActiveFunction)

    ! Deallocate workspace for DSYGVX
    IF (Method == 'G') THEN
      DEALLOCATE(Glob_IWorkForDSYGVX)
      DEALLOCATE(Glob_WorkForDSYGVX)
    ELSE
      CALL ClearQWorkspace()
    ENDIF

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    DEALLOCATE(Glob_diagS)
    IF (Method == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF (Glob_ProcID == 0) THEN
      i = LEN_TRIM(FileName)
      WRITE(*, *) 'New basis has been saved in file', FileName(1:i)
      WRITE(*, *) 'Program will now stop'
    ENDIF

    CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop

  END SUBROUTINE SeparateFuncLargeCoeff


  SUBROUTINE SaveHSWF(FileName1, FileName2, FileName3, FileName4, GSEPSolMethod)
    !==================================================================
    ! Subroutine SaveHSWF
    !==================================================================
    ! Writes the Hamiltonian, the overlap matrix, the eigenvector and the
    ! wave function to four files (a name of ' ', 'none', 'NONE' or 'None'
    ! skips that file). GSEPSolMethod picks the solver ('G' DSYGVX, 'I'
    ! inverse iteration, 'Q' QR factorization); the eigenproblem is solved only when file 3 or 4
    ! is wanted. The routine RETURNS normally.
    ! File formats: H and S - one element per line preceded by its two
    ! indices, FULL matrix (both triangles); eigenvector - one entry per
    ! line with its index; wave function - a header (particles, masses,
    ! charges, Young operator, basis size, energy), then one line per
    ! function: index, linear coefficient, a colon, the power, the
    ! parameters. In 'I' mode Glob_H holds H - Glob_ApproxEnergy*S, so the
    ! shift is added back per element and the unit diagonal of S is written
    ! as ONE.
    !==================================================================

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName1      ! Hamiltonian
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName2      ! overlap matrix
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName3      ! eigenvector
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName4      ! whole wave function
    CHARACTER(1)                               :: GSEPSolMethod  ! 'G' = DSYGVX, 'I' = inverse iteration

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and sizes ------------------------------------------
    INTEGER :: i, j                ! loop counters
    INTEGER :: cbs                 ! Glob_CurrBasisSize
    INTEGER :: n                   ! Glob_n, pseudoparticles
    INTEGER :: npt                 ! Glob_npt, parameters per function
    INTEGER :: BlockSizeForDSYGVX  ! ILAENV block size ('G' only)

    ! -- which files were asked for ----------------------------------
    LOGICAL :: IsHNeeded     ! FileName1 names a real file
    LOGICAL :: IsSNeeded     ! FileName2   "
    LOGICAL :: IsEVNeeded    ! FileName3   "
    LOGICAL :: IsWFNeeded    ! FileName4   "
    LOGICAL :: IsSwapFileOK  ! H and S came from the swap file

    ! -- eigensolver -------------------------------------------------
    INTEGER  :: ErrorCode          ! DSYGVX INFO, or GSEPIIS status
    INTEGER  :: NumOfEigvecs       ! eigenvectors asked of DSYGVX
    INTEGER  :: NumOfEigvalsFound  ! DSYGVX M
    INTEGER  :: NumOfIterations    ! inverse iterations used
    REAL(wp) :: Evalue             ! the eigenvalue, on every rank

    ! -- declared but NOT REFERENCED ---------------------------------
    INTEGER :: np

    ! -- DSYGVX output, 'G' path only --------------------------------
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: Eigvals
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: Eigvecs
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: IFAIL


    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine SaveHSWF started'
      IF (Verbose >= 1) WRITE(*, *) 'Number of basis functions', Glob_CurrBasisSize
      IF (Verbose >= 1) WRITE(*, *) 'GSEP solution method ', GSEPsolMethod
    ENDIF

    IF ((GSEPsolMethod /= 'G') .AND. (GSEPsolMethod /= 'I') .AND. (GSEPsolMethod /= 'Q')) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0190 in SaveHSWF: wrong GSEP solution method'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF


    !==================================================================
    ! Global state and what has been asked for
    !==================================================================
    ! DSYGVX is asked for a few eigenvectors beyond the wanted one so
    ! that Glob_WhichEigenvalue can be picked out of the set; inverse
    ! iteration returns the one vector it converges to.
    !------------------------------------------------------------------
    ! Setting the values of some global and local variables
    Glob_GSEPSolutionMethod = GSEPsolMethod
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    n = Glob_n
    np = Glob_np
    npt = Glob_npt
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize
    IF (GSEPsolMethod == 'G') NumOfEigvecs = MIN(cbs, Glob_WhichEigenvalue+10)
    IF (GSEPsolMethod == 'I') NumOfEigvecs = 1
    IF (GSEPsolMethod == 'Q') NumOfEigvecs = 1

    ! Setting logical variables that determine if everything (H, S, eigenvector, wave function)
    ! needs to be saved
    IF ((FileName1 == ' ') .OR. (FileName1 == 'none') .OR. &
        (FileName1 == 'NONE') .OR. (FileName1 == 'None')) THEN
      IsHNeeded = .FALSE.
    ELSE
      IsHNeeded = .TRUE.
    ENDIF
    IF ((FileName2 == ' ') .OR. (FileName2 == 'none') .OR. &
        (FileName2 == 'NONE') .OR. (FileName2 == 'None')) THEN
      IsSNeeded = .FALSE.
    ELSE
      IsSNeeded = .TRUE.
    ENDIF
    IF ((FileName3 == ' ') .OR. (FileName3 == 'none') .OR. &
        (FileName3 == 'NONE') .OR. (FileName3 == 'None')) THEN
      IsEVNeeded = .FALSE.
    ELSE
      IsEVNeeded = .TRUE.
    ENDIF
    IF ((FileName4 == ' ') .OR. (FileName4 == 'none') .OR. &
        (FileName4 == 'NONE') .OR. (FileName4 == 'None')) THEN
      IsWFNeeded = .FALSE.
    ELSE
      IsWFNeeded = .TRUE.
    ENDIF

    ! Nothing asked for at all. Leaving now skips allocating the
    ! matrices and computing every matrix element for output that
    ! would never be written. The same check is in workproc.f90.
    IF (.NOT. (IsHNeeded .OR. IsSNeeded .OR. IsEVNeeded .OR. IsWFNeeded)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'SaveHSWF: all four file names are none - nothing to save.'
        IF (Verbose >= 1) WRITE(*, *) 'Routine SaveHSWF finished'
      ENDIF
      RETURN
    ENDIF


    !==================================================================
    ! Allocate, by solution method
    !==================================================================
    ! 'G' keeps the H diagonal in Glob_diagH; 'I' keeps it inside
    ! Glob_H and needs Glob_invD for the LDL' factorization instead.
    !------------------------------------------------------------------
    ! Allocate global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (GSEPsolMethod == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    IF (GSEPsolMethod == 'I') ALLOCATE(Glob_invD(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (GSEPsolMethod == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Allocate workspace for subroutine GSEPIIS
    IF (GSEPsolMethod == 'I') THEN
      ALLOCATE(Glob_WorkForGSEPIIS(cbs))
      ALLOCATE(Glob_LastEigvector(cbs))
      Glob_LastEigvector(1:cbs) = ONE
    ENDIF

    ! Allocate local arrays
    IF (GSEPsolMethod == 'G') THEN
      ALLOCATE(Eigvals(NumOfEigvecs))
      ALLOCATE(Eigvecs(cbs, NumOfEigvecs))
      ALLOCATE(IFAIL(cbs))
    ENDIF


    !==================================================================
    ! Get the matrix elements
    !==================================================================
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a52)', ADVANCE='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
    ENDIF


    !==================================================================
    ! METHOD 'G' - DSYGVX
    !==================================================================
    IF (GSEPSolMethod == 'G') THEN

      ! Both triangles are filled: DSYGVX reads them, and the files
      ! below carry the full matrix.
      DO i = 1, cbs
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO
      DO i = 1, cbs
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO

      ! Saving matrices H and S:
      IF (Glob_ProcID == 0) THEN

        IF (IsHNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving the Hamiltonian matrix...'
          OPEN(2, FILE=FileName1)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              CALL writerealadv(2, Glob_H(i, j))
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF

        IF (IsSNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving the overlap matrix...'
          OPEN(2, FILE=FileName2)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              CALL writerealadv(2, Glob_S(i, j))
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF

      ENDIF

      ! Solved on rank 0 only, then broadcast. Skipped entirely when
      ! neither the eigenvector nor the wave function was asked for.
      IF ((IsEVNeeded) .OR. (IsWFNeeded)) THEN

        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
          CALL DSYGVX(1, 'V', 'I', 'U', cbs, Glob_H, Glob_HSLeadDim, Glob_S, Glob_HSLeadDim, &
                      ZERO, ZERO, 1, NumOfEigvecs, Glob_AbsTolForDSYGVX, &
                      NumOfEigvalsFound, Eigvals, Eigvecs, cbs, Glob_WorkForDSYGVX, Glob_LWorkForDSYGVX, &
                      Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
          ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
                      !$        VL, VU, IL, IU, ABSTOL,
                      !$        M, W, Z, LDZ, WORK, LWORK,
                      !$        IWORK, IFAIL, INFO )
        ENDIF

        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        IF (ErrorCode /= 0) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'failed'
            WRITE(*, *) &
              'Error EC0191 in SaveHSWF: routine DSYGVX failed with error code', ErrorCode
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! Glob_WhichEigenvalue picks the wanted level out of the set
        ! sending the eigenvalue and the eigenvector to all processes
        IF (Glob_ProcID == 0) THEN
          Evalue = Eigvals(Glob_WhichEigenvalue)
          Glob_c(1:cbs) = Eigvecs(1:cbs, Glob_WhichEigenvalue)
        ENDIF
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Glob_c, cbs, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        Glob_CurrEnergy = Evalue

        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 2) WRITE(*, *) 'done'
          WRITE(*, *) 'Energy: ', Evalue
        ENDIF

      ENDIF

    ENDIF  ! if (GSEPSolMethod=='G')


    !==================================================================
    ! METHOD 'I' - inverse iteration
    !==================================================================
    ! Glob_H holds H - Glob_ApproxEnergy*S here, so the shift is added
    ! back element by element as the file is written. Only the LOWER
    ! triangle is stored, hence the three-way index test.
    !------------------------------------------------------------------
    IF (GSEPSolMethod == 'I') THEN

      ! Saving matrices H and S:
      IF (Glob_ProcID == 0) THEN

        IF (IsHNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving Hamiltonian matrix...'
          OPEN(2, FILE=FileName1)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              IF (i == j) CALL writerealadv(2, Glob_H(i, j)+Glob_ApproxEnergy)
              IF (i > j) CALL writerealadv(2, Glob_H(i, j)+Glob_ApproxEnergy*Glob_S(i, j))
              IF (i < j) CALL writerealadv(2, Glob_H(j, i)+Glob_ApproxEnergy*Glob_S(j, i))
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF

        IF (IsSNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving overlap matrix...'
          OPEN(2, FILE=FileName2)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              IF (i == j) THEN
                CALL writerealadv(2, ONE)
              ELSE
                CALL writerealadv(2, Glob_S(i, j))
              ENDIF
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF

      ENDIF

      IF ((IsEVNeeded) .OR. (IsWFNeeded)) THEN

        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'

        ! A one-function basis has nothing to iterate on: the energy
        ! is the single diagonal element.
        IF (cbs == 1) THEN
          Glob_CurrEnergy = Glob_diagH(1)
          NumOfIterations = 1
          ErrorCode = 0
        ELSE
          CALL GSEPIIS(1, cbs, Glob_H, Glob_HSLeadDim, Glob_invD, Glob_S, Glob_HSLeadDim, &
                       Glob_ApproxEnergy, Glob_LastEigvector, Glob_WorkForGSEPIIS, Glob_EigvalTol, &
                       Evalue, Glob_c, Glob_LastEigvalTol, Glob_MaxIterForGSEPIIS, &
                       0, NumOfIterations, ErrorCode)
          ! GSEPIIS(k,n,M,nM,invD,B,nB, &
          !        apprlambda,v,w,Tol, &
          !        lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
          IF (Glob_LastEigvalTol > Glob_WorstEigvalTol) Glob_WorstEigvalTol = Glob_LastEigvalTol
          IF (Glob_LastEigvalTol > Glob_BestEigvalTol) Glob_BestEigvalTol = Glob_LastEigvalTol
          CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_BCAST(Glob_c, cbs, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        ENDIF

        Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
        Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
        Glob_CurrEnergy = Evalue

        IF (ErrorCode /= 0) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'failed'
            WRITE(*, *) 'Error EC0192 in SaveHSWF: the energy cannot be computed'
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF

        ! print the energy
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 2) WRITE(*, *) 'done'
          WRITE(*, *) 'Energy: ', Evalue
        ENDIF

      ENDIF

    ENDIF

    !==================================================================
    ! METHOD 'Q' - QR factorization
    !==================================================================
    ! The physical matrices already use the canonical normalized lower
    ! triangles with both diagonals in place; the files carry the full
    ! matrix, formed element by element while writing, and the unused
    ! upper triangles are left untouched.
    !------------------------------------------------------------------
    IF (GSEPSolMethod == 'Q') THEN
      IF (Glob_ProcID == 0) THEN
        IF (IsHNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving Hamiltonian matrix...'
          OPEN(2, FILE=FileName1)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              CALL writerealadv(2, QCanonicalMatrixElement(Glob_H, i, j))
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF
        IF (IsSNeeded) THEN
          IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving overlap matrix...'
          OPEN(2, FILE=FileName2)
          DO i = 1, cbs
            DO j = 1, cbs
              WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
              CALL writerealadv(2, QCanonicalMatrixElement(Glob_S, i, j))
            ENDDO
          ENDDO
          CLOSE(2)
          IF (Verbose >= 2) WRITE(*, *) 'done'
        ENDIF
      ENDIF
      IF ((IsEVNeeded) .OR. (IsWFNeeded)) THEN
        IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
        ! ONE gives inverse iteration a deterministic nonzero starting
        ! vector; SolveQ returns an S-normalized coefficient vector
        Glob_c = ONE
        CALL PrepareQWorkspace(cbs, cbs, 1, ErrorCode)
        Q_Workspace%MatricesAreCanonical = .TRUE.
        IF (ErrorCode == Q_METHOD_SUCCESS) CALL FactorizeQFresh(ErrorCode)
        IF (ErrorCode == Q_METHOD_SUCCESS) CALL SolveQ(Evalue, ErrorCode)
        IF (ErrorCode /= Q_METHOD_SUCCESS) THEN
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'failed'
            WRITE(*, '(1x,a,1x,i0)') 'Error EC0193 in SaveHSWF: Q energy cannot be computed, status', ErrorCode
          ENDIF
          CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
        ENDIF
        Glob_CurrEnergy = Evalue
        IF (Glob_ProcID == 0) THEN
          IF (Verbose >= 2) WRITE(*, *) 'done'
          WRITE(*, *) 'Energy: ', Evalue
        ENDIF
      ENDIF
    ENDIF


    !==================================================================
    ! Write the eigenvector and the wave function
    !==================================================================
    ! Both on rank 0 only. The wave-function header carries enough to
    ! rebuild the problem: particle count, masses, charges, the Young
    ! operator, the basis size and the energy.
    !------------------------------------------------------------------
    ! Saving the eigenvector
    IF ((IsEVNeeded) .AND. (Glob_ProcID == 0)) THEN
      IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving eigenvector...'
      OPEN(2, FILE=FileName3)
      DO i = 1, cbs
        WRITE(2, '(1x,i6,1x)', ADVANCE='no') i
        CALL writerealadv(2, Glob_c(i))
      ENDDO
      CLOSE(2)
      IF (Verbose >= 2) WRITE(*, *) 'done'
    ENDIF

    ! Saving the wave function
    IF ((IsWFNeeded) .AND. (Glob_ProcID == 0)) THEN
      IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving wave function...'
      OPEN(2, FILE=FileName4)
      WRITE(2, '(1x,a)') 'PG_0S WAVE FUNCTION FILE'
      WRITE(2, '(1x,a9,1x,i6)') 'PARTICLES', Glob_n+1
      WRITE(2, '(1x,a6)', ADVANCE='no') 'MASSES'
      CALL writerealarradv(2, Glob_Mass, Glob_n+1)
      WRITE(2, '(1x,a7)', ADVANCE='no') 'CHARGES'
      CALL writereal(2, Glob_PseudoCharge0)
      CALL writerealarradv(2, Glob_PseudoCharge, Glob_n)
      j = LEN_TRIM(Glob_YOperatorString)
      WRITE(2, '(1x,a8)', ADVANCE='no') 'SYMMETRY'
      CALL writestringadv(2, Glob_YOperatorString, j)
      WRITE(2, '(1x,a10,1x,i6)') 'BASIS_SIZE', Glob_CurrBasisSize
      WRITE(2, '(1x,a14)', ADVANCE='no') 'CURRENT_ENERGY'
      CALL writerealadv(2, Glob_CurrEnergy)
      WRITE(2, *) '=============================='
      ! One line per function: index, coefficient, colon, power, parameters
      DO i = 1, cbs
        WRITE(2, '(1x,i6,1x)', ADVANCE='no') i
        CALL writereal(2, Glob_c(i))
        WRITE(2, '(1x,a1,1x)', ADVANCE='no') ':'
        WRITE(2, '(i6,1x)', ADVANCE='no') Glob_PWR(i)
        CALL writerealarradv(2, Glob_NonlinParam(1:Glob_npt, i), Glob_npt)
      ENDDO
      CLOSE(2)
      IF (Verbose >= 2) WRITE(*, *) 'done'
    ENDIF


    !==================================================================
    ! Release everything, by solution method
    !==================================================================
    IF (GSEPSolMethod == 'G') THEN
      DEALLOCATE(IFAIL)
      DEALLOCATE(Eigvecs)
      DEALLOCATE(Eigvals)
    ENDIF

    IF (GSEPsolMethod == 'I') THEN
      DEALLOCATE(Glob_LastEigvector)
      DEALLOCATE(Glob_WorkForGSEPIIS)
    ENDIF
    IF (GSEPSolMethod == 'Q') CALL ClearQWorkspace()

    ! deallocate workspace for DSYGVX
    IF (GSEPSolMethod == 'G') THEN
      DEALLOCATE(Glob_WorkForDSYGVX)
      DEALLOCATE(Glob_IWorkForDSYGVX)
    ENDIF

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    IF (GSEPSolMethod == 'I') DEALLOCATE(Glob_invD)
    DEALLOCATE(Glob_diagS)
    IF (GSEPSolMethod == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine SaveHSWF finished'


  END SUBROUTINE SaveHSWF


  SUBROUTINE SaveHSRaw(FileName1, FileName2, GSEPSolMethod)
    !==================================================================
    ! Subroutine SaveHSRaw
    !==================================================================
    ! Writes the UNNORMALIZED Hamiltonian and overlap matrices:
    !     SAVE_HS_R  <Method>  <BasisSize>  <File1> <File2>
    ! File1 receives H, File2 S, in the 'i j value' layout of SaveHSWF
    ! (' ', 'none', 'NONE', 'None' skips a file). RETURNS normally.
    ! StoreHS keeps only normalized elements plus the raw diagonal in
    ! Glob_diagS, so the raw values are reconstructed at write time:
    !     S_raw(i,i) = Glob_diagS(i)
    !     S_raw(i,j) = S_norm(i,j) * SQRT(diagS(i))*SQRT(diagS(j))
    !     H_raw(i,j) = H_norm(i,j) * SQRT(diagS(i))*SQRT(diagS(j))
    ! (multiplications only, so full working precision; SQRT(a)*SQRT(b)
    ! rather than SQRT(a*b) so a near-null diagS cannot underflow). In 'I'
    ! mode the shift is added back before writing. Matrix elements come from
    ! the swap file when valid, else are recomputed; no eigenproblem is
    ! solved. Afterwards the BASIS HEALTH TABLE (cancellation inside the
    ! self-overlaps) is printed on the screen and into Glob_HealthFileName.
    !==================================================================

    IMPLICIT NONE

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName1      ! unnormalized H
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName2      ! unnormalized S
    CHARACTER(1), INTENT(IN)                   :: GSEPSolMethod  ! 'G', 'I' or 'Q': how H is stored
    CHARACTER(1)                               :: Method         ! storage used: 'G' also for 'Q'

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER  :: i, j                                         ! loop counters
    INTEGER  :: cbs                                          ! Glob_CurrBasisSize
    INTEGER  :: nflag                                        ! functions with C above Glob_MaxSelfOverlapCancel
    INTEGER  :: OpenFileErr                                  ! IOSTAT of the health-file OPEN
    INTEGER  :: u, iu                                        ! output unit and its selector
    LOGICAL  :: IsHNeeded                                    ! FileName1 names a real file
    LOGICAL  :: IsSNeeded                                    ! FileName2   "
    LOGICAL  :: IsSwapFileOK                                 ! H and S came from the swap file
    REAL(wp) :: Scale, Sraw, Hraw                            ! normalization factor and raw elements
    REAL(wp) :: Ssum, Sabs, Cfac, Dlost, Dhave, Dleft, Rerr  ! health table

    ! Buffers for the diagonal-only recomputation used by the health table.
    ! Only the cbs diagonal pairs are redone here, not all cbs(cbs+1)/2, so
    ! this costs a fraction of the matrix build.


    cbs = Glob_CurrBasisSize

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine SaveHSRaw started'
      IF (Verbose >= 1) WRITE(*, *) 'Number of basis functions', cbs
      IF (Verbose >= 1) WRITE(*, *) 'GSEP solution method ', GSEPSolMethod
    ENDIF

    IF ((GSEPSolMethod /= 'G') .AND. (GSEPSolMethod /= 'I') .AND. (GSEPSolMethod /= 'Q')) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0203 in SaveHSRaw: wrong GSEP solution method'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
    ENDIF

    IF (cbs < 1) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'SaveHSRaw: empty basis, nothing to do.'
      RETURN
    ENDIF

    ! Setting logical variables that determine what needs to be saved
    IF ((FileName1 == ' ') .OR. (FileName1 == 'none') .OR. &
        (FileName1 == 'NONE') .OR. (FileName1 == 'None')) THEN
      IsHNeeded = .FALSE.
    ELSE
      IsHNeeded = .TRUE.
    ENDIF
    IF ((FileName2 == ' ') .OR. (FileName2 == 'none') .OR. &
        (FileName2 == 'NONE') .OR. (FileName2 == 'None')) THEN
      IsSNeeded = .FALSE.
    ELSE
      IsSNeeded = .TRUE.
    ENDIF

    IF ((.NOT. IsHNeeded) .AND. (.NOT. IsSNeeded)) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'SaveHSRaw: both file names are none - nothing to save.'
        IF (Verbose >= 1) WRITE(*, *) 'Routine SaveHSRaw finished'
      ENDIF
      RETURN
    ENDIF


    !==================================================================
    ! Global state and allocation, as in SaveHSWF minus the solver arrays
    !==================================================================
    ! The raw elements do not depend on the solver: 'Q' uses the 'G' storage
    Method = GSEPSolMethod
    IF (Method == 'Q') Method = 'G'
    Glob_GSEPSolutionMethod = Method
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = cbs
    Glob_HSBuffLen = MAX(MIN(cbs*(cbs+1)/2, 1000), 30*cbs)

    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (Method == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))


    !==================================================================
    ! Get the matrix elements
    !==================================================================
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a52)', ADVANCE='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
    ENDIF

    ! Mirror the stored lower triangles into the upper ones. 'G' keeps
    ! the H diagonal apart in Glob_diagH and leaves S's implicit; 'I' has
    ! already written both diagonals into the matrices themselves.
    DO i = 1, cbs
      DO j = 1, i-1
        Glob_H(j, i) = Glob_H(i, j)
        Glob_S(j, i) = Glob_S(i, j)
      ENDDO
    ENDDO
    IF (Method == 'G') THEN
      DO i = 1, cbs
        Glob_H(i, i) = Glob_diagH(i)
        Glob_S(i, i) = ONE
      ENDDO
    ENDIF


    !==================================================================
    ! Undo the normalization and write
    !==================================================================
    IF (Glob_ProcID == 0) THEN

      IF (IsHNeeded) THEN
        IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving unnormalized Hamiltonian matrix...'
        OPEN(2, FILE=FileName1, STATUS='replace')
        DO i = 1, cbs
          DO j = 1, cbs
            Scale = SQRT(Glob_diagS(i))*SQRT(Glob_diagS(j))
            IF (i == j) THEN
              Sraw = Glob_diagS(i)
            ELSE
              Sraw = Glob_S(i, j)*Scale
            ENDIF
            IF (Method == 'G') THEN
              IF (i == j) THEN
                Hraw = Glob_diagH(i)*Glob_diagS(i)
              ELSE
                Hraw = Glob_H(i, j)*Scale
              ENDIF
            ELSE
              IF (i == j) THEN
                Hraw = (Glob_H(i, i)+Glob_ApproxEnergy)*Glob_diagS(i)
              ELSE
                Hraw = Glob_H(i, j)*Scale+Glob_ApproxEnergy*Sraw
              ENDIF
            ENDIF
            WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
            CALL writerealadv(2, Hraw)
          ENDDO
        ENDDO
        CLOSE(2)
        IF (Verbose >= 2) WRITE(*, *) 'done'
      ENDIF

      IF (IsSNeeded) THEN
        IF (Verbose >= 2) WRITE(*, '(1x,a)', ADVANCE='no') 'Saving unnormalized overlap matrix...'
        OPEN(2, FILE=FileName2, STATUS='replace')
        DO i = 1, cbs
          DO j = 1, cbs
            IF (i == j) THEN
              Sraw = Glob_diagS(i)
            ELSE
              Sraw = Glob_S(i, j)*SQRT(Glob_diagS(i))*SQRT(Glob_diagS(j))
            ENDIF
            WRITE(2, '(1x,i6,1x,i6,1x)', ADVANCE='no') i, j
            CALL writerealadv(2, Sraw)
          ENDDO
        ENDDO
        CLOSE(2)
        IF (Verbose >= 2) WRITE(*, *) 'done'
      ENDIF

    ENDIF


    !==================================================================
    ! Basis health table
    !==================================================================
    ! <phi_i|Y+Y|phi_i> is a sum of Glob_NumYHYTerms signed terms; when they
    ! cancel, the result carries fewer correct digits than the working
    ! precision. C = sum|c_k S_k| / |sum c_k S_k| is the factor by which
    ! machine epsilon is amplified: log10(C) digits are lost and
    !-log10(epsilon) - log10(C) survive. Only the diagonal is recomputed
    ! (cbs pairs), on rank 0 alone.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN

      Dhave = -LOG10(EPSILON(Sraw))

      ! STATUS='replace', not 'unknown': a shorter run must not leave the
      ! tail of a longer earlier one underneath its own output.
      OPEN(3, FILE=Glob_HealthFileName, STATUS='replace', IOSTAT=OpenFileErr)
      IF ((OpenFileErr /= 0) .AND. (Verbose >= 1)) WRITE(*, *) 'Warning in SaveHSRaw: could not open ', &
                                     TRIM(Glob_HealthFileName)

      IF (Verbose >= 1) WRITE(*, *)
      WRITE(*, *) 'Basis health from the diagonal of Y+Y'
      WRITE(*, '(1x,a,f5.2,a)') 'Working precision: -log10(EPSILON) = ', Dhave, &
                               ' decimal digits'
      WRITE(*, *)
      WRITE(*, '(1x,a)') '  func  power        S_raw(i,i)        sum|terms|'// &
                        '     cancel C   lost   left      rel.err'
      IF (OpenFileErr == 0) THEN
        WRITE(3, *) 'Basis health from the diagonal of Y+Y'
        WRITE(3, '(1x,a,f5.2,a)') 'Working precision: -log10(EPSILON) = ', Dhave, &
                                 ' decimal digits'
        WRITE(3, *)
        WRITE(3, '(1x,a)') '  func  power        S_raw(i,i)        sum|terms|'// &
                          '     cancel C   lost   left      rel.err'
      ENDIF

      nflag = 0
      DO i = 1, cbs
        Cfac = SelfOverlapCancellation(Glob_PWR(i), Glob_NonlinParam(1:Glob_npt, i), Ssum, Sabs)
        Dlost = LOG10(Cfac)
        Dleft = Dhave-Dlost
        IF (Dleft < ZERO) Dleft = ZERO
        IF (Cfac > Glob_MaxSelfOverlapCancel) nflag = nflag+1

        ! The relative error carried by S_raw(i,i): machine epsilon
        ! amplified by the cancellation. Same statement as 'left', but
        ! expressed as a plain accuracy instead of a count of digits.
        Rerr = Cfac*EPSILON(Sraw)
        IF (Rerr > ONE) Rerr = ONE

        WRITE(*, '(1x,i6,1x,i6,3(1x,e17.9),2(1x,f6.1),1x,e13.5)') &
          i, Glob_PWR(i), Ssum, Sabs, Cfac, Dlost, Dleft, Rerr
        IF (OpenFileErr == 0) &
          WRITE(3, '(1x,i6,1x,i6,3(1x,e17.9),2(1x,f6.1),1x,e13.5)') &
            i, Glob_PWR(i), Ssum, Sabs, Cfac, Dlost, Dleft, Rerr
      ENDDO


      !----------------------------------------------------------------
      ! Summary and legend, sent to the screen and to the mirror file
      !----------------------------------------------------------------
      DO iu = 1, 2
        IF (iu == 1) THEN
          u = 6
        ELSE
          IF (OpenFileErr /= 0) CYCLE
          u = 3
        ENDIF

        WRITE(u, *)
        WRITE(u, '(1x,a,es8.1,a,i6)') 'Functions with cancellation above ', Glob_MaxSelfOverlapCancel, ': ', nflag
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' HOW TO READ THE TABLE'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' S_raw(i,i)  The self-overlap <phi_i|Y+Y|phi_i>, not normalized. It is'
        WRITE(u, '(1x,a)') '             what is LEFT of the basis function after the Young operator'
        WRITE(u, '(1x,a)') '             has acted on it. Same number as the diagonal of the raw S'
        WRITE(u, '(1x,a)') '             file written by this step.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' sum|terms|  Y+Y is a sum of many signed terms. This is the sum of their'
        WRITE(u, '(1x,a)') '             absolute values: how big the numbers were BEFORE they were'
        WRITE(u, '(1x,a)') '             added together.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' cancel C    = sum|terms| / |S_raw(i,i)|   -  how much cancelled.'
        WRITE(u, '(1x,a)') '               C = 1        nothing cancelled, the terms just added up'
        WRITE(u, '(1x,a)') '               C = 1000     the terms were 1000 times bigger than the'
        WRITE(u, '(1x,a)') '                            number that survived'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' lost        = log10(C).  The number of significant digits destroyed.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' left        = working precision - lost.  Significant digits still correct.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' rel.err     = C * EPSILON.  The same statement as "left", written as a'
        WRITE(u, '(1x,a)') '             plain relative accuracy instead of a count of digits.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' WHY lost = log10(C)'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' The one idea to hold on to:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     NUMBER OF CORRECT DIGITS IS ITSELF A LOGARITHM.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' "n correct digits" and "relative error 10^-n" are the same statement:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     relative error 1E-16   <->   16 correct digits'
        WRITE(u, '(1x,a)') '     relative error 1E-06   <->    6 correct digits'
        WRITE(u, '(1x,a)') '     digits = -log10(relative error)'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' That is also why the working precision above is -log10(EPSILON).'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Now, cancellation multiplies the relative error by C:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     error after  =  C  x  error before'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Take -log10 of both sides. A logarithm turns multiplication into'
        WRITE(u, '(1x,a)') ' subtraction, so:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     digits after  =  digits before  -  log10(C)'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' and the amount subtracted, log10(C), is what the "lost" column shows.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' The plain-language version: log10(C) is simply HOW MANY FACTORS OF TEN'
        WRITE(u, '(1x,a)') ' are in C, and each factor of ten eats exactly one digit.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     C = 10          1 digit lost'
        WRITE(u, '(1x,a)') '     C = 100         2 digits lost'
        WRITE(u, '(1x,a)') '     C = 1 000       3 digits lost'
        WRITE(u, '(1x,a)') '     C = 10^10      10 digits lost'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' A WORKED EXAMPLE'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Add four numbers, carrying 16 digits:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '       1000.0  -  999.0  +  500.0  -  500.9   =   0.1'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '   sum|terms| = 1000 + 999 + 500 + 500.9      =   2999.9'
        WRITE(u, '(1x,a)') '   C          = 2999.9 / 0.1                  =   30000'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Each term is stored to 16 digits, so 1000 carries an ABSOLUTE error of'
        WRITE(u, '(1x,a)') ' 1000 x 1E-16 = 1E-13. The errors add up just like the numbers do:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '   absolute error of the total  ~  2999.9 x 1E-16  =  3E-13'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' The errors did NOT cancel - only the numbers did. The answer shrank to'
        WRITE(u, '(1x,a)') ' 0.1 while its error stayed at 3E-13, so'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '   relative error  =  3E-13 / 0.1  =  3E-12   ->  about 11.5 digits left'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Started with 16 digits, ended with 11.5: 4.5 digits lost. And indeed'
        WRITE(u, '(1x,a)') ' log10(30000) = 4.5.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' WHY THE VALUES ARE NOT WHOLE NUMBERS'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Because they are logarithms, not counts. A fraction is a factor:'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     C = 1E+10   ->  lost = 10.0'
        WRITE(u, '(1x,a)') '     C = 3E+10   ->  lost = 10.5      (log10(3) = 0.5)'
        WRITE(u, '(1x,a)') '     C = 1E+11   ->  lost = 11.0'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' So "10.5 digits lost" means the cancellation is three times worse than'
        WRITE(u, '(1x,a)') ' a clean 1E+10 - it is one number, ten-and-a-half, not "10 and 5".'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Likewise "left = 5.2 digits" means a bit better than 5 correct digits'
        WRITE(u, '(1x,a)') ' and not quite 6. If whole numbers are easier, round DOWN and read the'
        WRITE(u, '(1x,a)') ' rel.err column, which says the same thing without any logarithms.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' WHY IT MATTERS'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Normalization divides the elements of row i and column i by'
        WRITE(u, '(1x,a)') ' SQRT(S_raw(i,i)*S_raw(j,j)). An error in S_raw(i,i) therefore does not'
        WRITE(u, '(1x,a)') ' stay on the diagonal - it spreads over that whole row and column, at'
        WRITE(u, '(1x,a)') ' half the relative size, because the square root halves it.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' So one bad function limits the accuracy of the entire calculation,'
        WRITE(u, '(1x,a)') ' not only of itself.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' WHAT MAKES C LARGE'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' The Young operator antisymmetrizes some exchanges and symmetrizes'
        WRITE(u, '(1x,a)') ' others. A function that is nearly SYMMETRIC under an exchange the'
        WRITE(u, '(1x,a)') ' operator ANTISYMMETRIZES is nearly annihilated: the terms almost'
        WRITE(u, '(1x,a)') ' cancel and only a tiny residual survives. Normalization then scales'
        WRITE(u, '(1x,a)') ' that residual back up to unit norm, magnifying every rounding error'
        WRITE(u, '(1x,a)') ' inside it by the same factor C.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' To find WHICH exchange, rerun with SYMMETRY (1+Pij), then with'
        WRITE(u, '(1x,a)') ' (1-Pij), for one pair ij at a time, and compare S_raw(i,i):'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '     <phi|Pij|phi>/<phi|phi> = (d_plus - d_minus)/(d_plus + d_minus)'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' A value near +1 means the function is symmetric under that exchange.'
        WRITE(u, '(1x,a)') ' If that exchange sits in an antisymmetrizer, it is the culprit.'
        WRITE(u, *)
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, '(1x,a)') ' WHAT TO DO'
        WRITE(u, '(1x,a)') '=================================================================='
        WRITE(u, *)
        WRITE(u, '(1x,a)') '   C below 1.0E+03    normal, nothing to do'
        WRITE(u, '(1x,a)') '   C up to  1.0E+06    still usable in double precision'
        WRITE(u, '(1x,a)') '   C above  1.0E+06    re-optimize that function, or run PREC=10'
        WRITE(u, '(1x,a)') '                       (extended precision buys about 3 more digits)'
        WRITE(u, *)
        WRITE(u, '(1x,a)') ' Re-optimizing is usually better than deleting: such a function often'
        WRITE(u, '(1x,a)') ' still lowers the energy, it is only its norm that is badly determined.'
        WRITE(u, *)
      ENDDO

      IF (OpenFileErr == 0) THEN
        CLOSE(3)
        WRITE(*, *) 'Table also written to ', TRIM(Glob_HealthFileName)
      ENDIF

    ENDIF


    !==================================================================
    ! Release everything
    !==================================================================
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_diagS)
    IF (Method == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE(*, *) 'Routine SaveHSRaw finished'

  END SUBROUTINE SaveHSRaw


  SUBROUTINE OverlapDiag(GSEPSolMethod, FileName)
    !==================================================================
    ! Subroutine OverlapDiag
    !==================================================================
    ! Diagonalizes the OVERLAP MATRIX ALONE (S*x = lambda*x, the normalized
    ! S with unit diagonal that the solvers use) and reports its spectrum:
    !     OVERLAP_D  <Method>  <BasisSize>  [<FileName>]
    ! (fires when BasisSize equals the current basis size; RETURNS
    ! normally). FileName receives the full spectrum; ReadIOFile fills it
    ! with Glob_OverlapFileName when the line names no file.
    ! Only 'G' (DSYEVX on rank 0) is
    ! accepted; the former 'I' path (unpivoted LDL^T inverse iteration) was
    ! removed because it gave less (no lambda_max, no condition number), was
    ! least reliable exactly when S is nearly singular, and was not faster.
    ! lambda_min -> 0 means a combination of basis functions is nearly the
    ! zero function; the eigenvalues sum to N; the condition number
    ! lambda_max/lambda_min says how many decimal digits the basis has eaten
    ! and when PREC=10 becomes necessary. Glob_WhichEigenvalue is not
    ! consulted (it names a state of the HAMILTONIAN). Matrix elements from
    ! the swap file when valid, else recomputed.
    !==================================================================

    IMPLICIT NONE

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    CHARACTER(1), INTENT(IN) :: GSEPSolMethod  ! only 'G' is accepted
    CHARACTER(*), INTENT(IN) :: FileName       ! mirror file for the spectrum

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER                             :: N                  ! Glob_CurrBasisSize
    INTEGER                             :: i, j               ! loop counters
    INTEGER                             :: ErrorCode          ! DSYEVX INFO
    INTEGER                             :: OpenFileErr        ! IOSTAT of the mirror-file OPEN
    INTEGER                             :: BlockSizeForDSYEVX ! ILAENV block size
    INTEGER                             :: LWork              ! DSYEVX workspace length
    INTEGER                             :: NumOfEigvalsFound  ! DSYEVX M (N for RANGE='A')
    INTEGER                             :: NLow, NHigh        ! eigenvalues listed in the low and high blocks
    INTEGER                             :: IHigh              ! first index of the high block
    REAL(wp)                            :: CondNum            ! lambda_max / lambda_min
    REAL(wp)                            :: Zdummy(1, 1)       ! DSYEVX Z, not referenced with JOBZ='N'
    LOGICAL                             :: IsFileOK           ! the mirror file could be opened
    LOGICAL                             :: IsSwapFileOK       ! H and S came from the swap file
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: Work               ! DSYEVX workspace
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: IWork              ! DSYEVX integer workspace (5*N)
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: IFAIL              ! DSYEVX IFAIL (used only with JOBZ='V')
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: Eigvals            ! the spectrum of S, ascending


    N = Glob_CurrBasisSize

    ! Only rank 0 opens the file, but IsFileOK is read on every rank in
    ! the CLOSE test below. Fortran does not promise to short-circuit
    ! .AND., so this has to be defined everywhere, not just where it is set.
    IsFileOK = .FALSE.

    IF (GSEPSolMethod /= 'G') THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error in OverlapDiag: eigensolver ', GSEPSolMethod, &
                   ' is not available for OVERLAP_D'
        IF (Verbose >= 2) WRITE(*, *) 'Only G is recognized. Skipping this step...'
      ENDIF
      RETURN
    ENDIF

    IF (N < 1) THEN
      IF (Glob_ProcID == 0) WRITE(*, *) 'OverlapDiag: empty basis, nothing to do.'
      RETURN
    ENDIF

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine OverlapDiag started'
      IF (Verbose >= 1) WRITE(*, *) 'Number of basis functions', N
      WRITE(*, *) 'Eigensolver G: DSYEVX, full spectrum of S'
    ENDIF

    ! Workspace size, same rule as everywhere else in this file: DSYEVX
    ! needs at least 8*N, which MAX((NB+3)*N, 8*N) always satisfies.
    BlockSizeForDSYEVX = ILAENV(1, 'DSYTRD', 'VIU', N, N, N, N)
    LWork = MAX((BlockSizeForDSYEVX+3)*N, 8*N)


    !==================================================================
    ! Global state and allocation
    !==================================================================
    ! 'G' storage puts the normalized off-diagonals of S in the lower
    ! triangle and leaves the unit diagonal implicit, which is the
    ! cheapest way to get S without the shift being folded into anything.
    ! Glob_H and Glob_diagH are allocated because StoreHS writes the
    ! Hamiltonian elements as a side effect of computing the overlap -
    ! they come out of the same matelem call and are simply not used here.
    !------------------------------------------------------------------
    Glob_GSEPSolutionMethod = 'G'
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = N
    Glob_HSBuffLen = MAX(MIN(N*(N+1)/2, 1000), 30*N)

    ALLOCATE(Glob_H(N, N))
    ALLOCATE(Glob_S(N, N))
    ALLOCATE(Glob_diagH(N))
    ALLOCATE(Glob_diagS(N))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))
    ALLOCATE(Work(LWork))
    ALLOCATE(IWork(5*N))
    ALLOCATE(IFAIL(N))
    ALLOCATE(Eigvals(N))


    !==================================================================
    ! Get the matrix elements
    !==================================================================
    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a52)', ADVANCE='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      CALL ComputeMatElem(1, N)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
    ENDIF

    ! Complete S: mirror the lower triangle into the upper and set the
    ! unit diagonal. After this Glob_S holds the full symmetric matrix.
    DO i = 1, N
      DO j = 1, i-1
        Glob_S(j, i) = Glob_S(i, j)
      ENDDO
      Glob_S(i, i) = ONE
    ENDDO


    !==================================================================
    ! Open the mirror file
    !==================================================================
    ! Everything below goes to the screen (summary plus extreme
    ! eigenvalues) and to this file (whole spectrum), as ExpectationValues
    ! mirrors into Glob_ExpValFileName. STATUS='replace', not 'unknown',
    ! so a shorter run does not leave the tail of an earlier one in place.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      OPEN(2, FILE=TRIM(FileName), STATUS='replace', IOSTAT=OpenFileErr)
      IsFileOK = (OpenFileErr == 0)
      IF (.NOT. IsFileOK) THEN
        IF (Verbose >= 1) WRITE(*, *) 'Warning in OverlapDiag: could not open ', TRIM(FileName)
        IF (Verbose >= 2) WRITE(*, *) 'Screen output only.'
      ELSE
        WRITE(2, *)
        WRITE(2, *) 'Routine OverlapDiag started'
        WRITE(2, *) 'Number of basis functions', N
        WRITE(2, *) 'Eigensolver G: DSYEVX, full spectrum of S'
      ENDIF
    ENDIF

    CondNum = ZERO


    !==================================================================
    ! Whole spectrum
    !==================================================================
    ! DSYEVX destroys the triangle it reads, which is fine - S is not
    ! needed afterwards. Eigenvalues come back in Eigvals in ASCENDING
    ! order. With JOBZ='N', RANGE='A' and ABSTOL=0 it takes the same path
    ! as DSYEV: DSYTRD, then DSTERF; VL, VU, IL and IU are not referenced.
    ! Solved on rank 0 and broadcast. A failure is not fatal: the
    ! step is a diagnostic, so it is reported and the step is skipped.
    !------------------------------------------------------------------
    IF (Glob_ProcID == 0) THEN
      CALL DSYEVX('N', 'A', 'U', N, Glob_S, Glob_HSLeadDim, ZERO, ZERO, 1, N, ZERO, &
                  NumOfEigvalsFound, Eigvals, Zdummy, 1, Work, LWork, IWork, IFAIL, ErrorCode)
    ENDIF
    CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    IF (ErrorCode /= 0) THEN

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Warning in OverlapDiag: routine DSYEVX failed with error code', ErrorCode
        WRITE(*, *) 'No eigenvalues are reported for this step.'
      ENDIF

    ELSE

      CALL MPI_BCAST(Eigvals, N, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

      IF (Eigvals(1) > ZERO) CondNum = Eigvals(N)/Eigvals(1)

      ! How many eigenvalues go in each block, and where the "highest"
      ! block starts. NLow and NHigh never overlap: on a basis smaller
      ! than 20 the two would otherwise list the same eigenvalues twice.
      NLow = MIN(10, N)
      NHigh = MIN(10, N-NLow)
      IHigh = N-NHigh+1

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        WRITE(*, *) 'Overlap matrix eigenvalues:'
        IF (Verbose >= 1) WRITE(*, '(a)') '-------------------------------------'
        ! Every '=' in this block sits in column 21: the label strings are
        ! padded to 21 characters, and the lambda lines reach the same
        ! column as 4 + len('lambda(') + 5 + len(')   ').
        WRITE(*, '(a)', ADVANCE='no') '    smallest        ='
        CALL writerealadv(6, Eigvals(1))
        WRITE(*, '(a)', ADVANCE='no') '    largest         ='
        CALL writerealadv(6, Eigvals(N))
        WRITE(*, '(a)', ADVANCE='no') '    sum (must be N) ='
        CALL writerealadv(6, SUM(Eigvals(1:N)))
        WRITE(*, '(a)', ADVANCE='no') '    condition number='
        CALL writerealadv(6, CondNum)

        WRITE(*, *)
        WRITE(*, '(a,i4,a)') '  lowest ', NLow, ' :'
        DO i = 1, NLow
          WRITE(*, '(4x,a,i5,a)', ADVANCE='no') 'lambda(', i, ')   ='
          CALL writerealadv(6, Eigvals(i))
        ENDDO

        IF (NHigh > 0) THEN
          WRITE(*, *)
          WRITE(*, '(a,i4,a)') '  highest ', NHigh, ' :'
          DO i = IHigh, N
            WRITE(*, '(4x,a,i5,a)', ADVANCE='no') 'lambda(', i, ')   ='
            CALL writerealadv(6, Eigvals(i))
          ENDDO
        ENDIF

        ! The file gets the WHOLE spectrum, in the same labelled form.
        IF (IsFileOK) THEN
          WRITE(2, *)
          WRITE(2, *) 'Overlap matrix eigenvalues:'
          WRITE(2, '(a)') '-------------------------------------'
          WRITE(2, '(a)', ADVANCE='no') '    smallest        ='
          CALL writerealadv(2, Eigvals(1))
          WRITE(2, '(a)', ADVANCE='no') '    largest         ='
          CALL writerealadv(2, Eigvals(N))
          WRITE(2, '(a)', ADVANCE='no') '    sum (must be N) ='
          CALL writerealadv(2, SUM(Eigvals(1:N)))
          WRITE(2, '(a)', ADVANCE='no') '    condition number='
          CALL writerealadv(2, CondNum)
          WRITE(2, *)
          DO i = 1, N
            WRITE(2, '(4x,a,i5,a)', ADVANCE='no') 'lambda(', i, ')   ='
            CALL writerealadv(2, Eigvals(i))
          ENDDO
          WRITE(*, *)
          WRITE(*, *) '   full spectrum written to ', TRIM(FileName)
        ENDIF
      ENDIF

    ENDIF

    IF ((Glob_ProcID == 0) .AND. IsFileOK) CLOSE(2)


    !==================================================================
    ! Release everything
    !==================================================================
    DEALLOCATE(Eigvals)
    DEALLOCATE(IFAIL)
    DEALLOCATE(IWork)
    DEALLOCATE(Work)
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_diagS)
    DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine OverlapDiag finished'
    ENDIF

  END SUBROUTINE OverlapDiag


  FUNCTION NumFieldWidth(vals) RESULT(w)
    !==================================================================
    ! Function NumFieldWidth
    !==================================================================
    ! Width of an F<w>.16 field that shows every value of vals with 16
    ! decimals, aligned on the decimal point, with two blanks in front
    ! of the widest value when the label before the field ends with one
    ! blank. Used by the screen output of ExpectationValues.
    !==================================================================
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: vals(:)
    INTEGER              :: w
    INTEGER        :: i, l
    CHARACTER(128) :: buf

    w = 18  ! 0.xxxxxxxxxxxxxxxx
    DO i = 1, SIZE(vals)
      WRITE(buf, '(f0.16)') vals(i)
      buf = ADJUSTL(buf)
      l = LEN_TRIM(buf)
      ! F0.16 may drop the zero before the point; F<w>.16 prints it
      IF (buf(1:1) == '.') l = l+1
      IF (buf(1:2) == '-.') l = l+1
      w = MAX(w, l)
    ENDDO
    w = w+1
  END FUNCTION NumFieldWidth


  SUBROUTINE PairSetDescription(iset, descr, single, haszero)
    !==================================================================
    ! Subroutine PairSetDescription
    !==================================================================
    ! Index-set notation of the equivalent-pair set iset (Glob_EqvPairList)
    ! for the screen output of ExpectationValues. A stored pair (a,a) is
    ! the distance of particle a from particle 0, i.e. the pair (0,a).
    !   'i=2,3,4,5,6,7'        all pairs (0,i): single = .TRUE., label r_i
    !   'i=1; j=2,3,4,5,6,7'   all pairs (i,j) with i in one list and j
    !                          in the other
    !   'i<j in 2,3,4,5,6,7'   all pairs within one list
    !   'ij=23,24,...'         any other set (does not occur for sets of
    !                          identical particles)
    ! Every particle number is listed. haszero: particle 0 appears in a
    ! list, i.e. r_0j stands for r_j on that line.
    !==================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN)       :: iset
    CHARACTER(*), INTENT(OUT) :: descr
    LOGICAL, INTENT(OUT)      :: single, haszero
    INTEGER      :: j, k, a, b, p, nu
    LOGICAL      :: ok
    LOGICAL      :: inU(0:Glob_n), inP(0:Glob_n), inB(0:Glob_n)
    LOGICAL      :: present(0:Glob_n, 0:Glob_n)
    CHARACTER(8) :: num

    k = Glob_NumOfPairsInEqvPairSet(iset)
    present = .FALSE.
    inU = .FALSE.
    single = .TRUE.
    DO j = 1, k
      a = Glob_EqvPairList(1, j, iset)
      b = Glob_EqvPairList(2, j, iset)
      IF (a == b) THEN
        a = 0
      ELSE
        single = .FALSE.
      ENDIF
      present(a, b) = .TRUE.
      present(b, a) = .TRUE.
      inU(a) = .TRUE.
      inU(b) = .TRUE.
    ENDDO
    nu = COUNT(inU)
    haszero = .FALSE.

    ! all pairs (0,i)
    IF (single) THEN
      inB = inU
      inB(0) = .FALSE.
      descr = 'i='//IndexList(inB)
      RETURN
    ENDIF

    ! one pair (a,b), a < b
    IF (k == 1) THEN
      a = Glob_EqvPairList(1, 1, iset)
      b = Glob_EqvPairList(2, 1, iset)
      WRITE(descr, '(a,i0,a,i0)') 'i=', a, '; j=', b
      RETURN
    ENDIF

    ! all pairs within one list
    IF (k == nu*(nu-1)/2) THEN
      ok = .TRUE.
      DO a = 0, Glob_n
        DO b = a+1, Glob_n
          IF (inU(a) .AND. inU(b) .AND. (.NOT. present(a, b))) ok = .FALSE.
        ENDDO
      ENDDO
      IF (ok) THEN
        descr = 'i<j in '//IndexList(inU)
        haszero = inU(0)
        RETURN
      ENDIF
    ENDIF

    ! all pairs (i,j), i in P and j in B: B = partners of the smallest index
    p = 0
    DO WHILE (.NOT. inU(p))
      p = p+1
    ENDDO
    inB = present(p, :)
    inP = inU .AND. (.NOT. inB)
    ok = (k == COUNT(inP)*COUNT(inB))
    DO a = 0, Glob_n
      DO b = 0, Glob_n
        IF (inP(a) .AND. inB(b) .AND. (.NOT. present(a, b))) ok = .FALSE.
      ENDDO
    ENDDO
    IF (ok) THEN
      descr = 'i='//TRIM(IndexList(inP))//'; j='//IndexList(inB)
      haszero = inU(0)
      RETURN
    ENDIF

    ! anything else: the pairs themselves
    descr = 'ij='
    DO j = 1, k
      a = Glob_EqvPairList(1, j, iset)
      b = Glob_EqvPairList(2, j, iset)
      IF (a == b) a = 0
      WRITE(num, '(i0,i0)') a, b
      IF (j == 1) THEN
        descr = TRIM(descr)//TRIM(num)
      ELSE
        descr = TRIM(descr)//','//TRIM(num)
      ENDIF
    ENDDO
    haszero = inU(0)

  CONTAINS

    FUNCTION IndexList(mask) RESULT(s)
      ! 'a,b,c' for the particle numbers with mask = .TRUE.
      LOGICAL, INTENT(IN) :: mask(0:Glob_n)
      CHARACTER(64)       :: s
      INTEGER      :: m
      CHARACTER(8) :: item
      s = ''
      DO m = 0, Glob_n
        IF (mask(m)) THEN
          WRITE(item, '(i0)') m
          IF (LEN_TRIM(s) == 0) THEN
            s = item
          ELSE
            s = TRIM(s)//','//TRIM(item)
          ENDIF
        ENDIF
      ENDDO
    END FUNCTION IndexList

  END SUBROUTINE PairSetDescription


  SUBROUTINE ExpectationValues(Action, SymmAdaptMethod, FileName1, FileName2, FileName3, FileName4, GSEPSolMethod)
    !==================================================================
    ! Subroutine ExpectationValues
    !==================================================================
    ! Computes expectation values over the current basis and prints them to
    ! the screen and to Glob_ExpValFileName; GSEPSolMethod picks the solver
    ! ('G' DSYGVX, 'I' inverse iteration, 'Q' QR factorization). RETURNS normally. Through
    ! ExpValuesMatElem of module matelem: S, T, V and H = T + V; for every
    ! pseudoparticle pair 1/r, r, r^2 and delta(r); the mass-velocity,
    ! Darwin and orbit-orbit corrections (also times alpha^2). Two-particle
    ! quantities are also averaged over the equivalent-pair sets of
    ! ProgramDataInit. SymmAdaptMethod must be 1 (Y'Y on the ket;
    ! ExpValuesMatElem takes a single symmetry matrix), else EC0194.
    !
    ! Action is 'EXPC_VALS' or 'DENSITIES'. As in the reference codes,
    ! DENSITIES additionally evaluates the correlation function on a
    ! user-supplied grid; only the nucleus-nucleus one (pair 0-1, r_1, the
    ! g_1 of the reference codes) is available in PG_0S, and the particle
    ! densities are not:
    !   FileName1 - grid for the correlation function (one value r >= 0 per
    !               line); 'none', 'NONE' or 'None' skips the correlation
    !               function. A bad or missing grid file stops the run
    !               (EC0196, EC0197).
    !   FileName2 - file where the correlation function is stored
    !   FileName3 - grid for the particle densities: not available; a name
    !               other than 'none' is ignored with warning WC0195
    !   FileName4 - file for the particle densities: NOT REFERENCED
    ! For EXPC_VALS the four names are not referenced.
    !==================================================================

    IMPLICIT NONE

    !------------------------------------------------------------------
    ! Arguments
    !------------------------------------------------------------------
    CHARACTER(9), INTENT(IN)                   :: Action           ! 'EXPC_VALS' or 'DENSITIES'
    INTEGER, INTENT(IN)                        :: SymmAdaptMethod  ! must be 1 - see the header
    CHARACTER(1)                               :: GSEPSolMethod    ! 'G' = DSYGVX, 'I' = inverse iteration
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName1        ! correlation-function grid
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName2        ! correlation-function output
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName3        ! density grid - not available
    CHARACTER(Glob_FileNameLength), INTENT(IN) :: FileName4        ! NOT REFERENCED

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------

    ! -- counters and sizes ------------------------------------------
    INTEGER :: i, j, k             ! loop counters
    INTEGER :: a, b, c             ! pair indices and position in MEkl
    INTEGER :: counter             ! basis pair number, for the MPI split
    INTEGER :: cbs                 ! Glob_CurrBasisSize
    INTEGER :: n                   ! Glob_n, pseudoparticles
    INTEGER :: npt                 ! Glob_npt, parameters per function
    INTEGER :: np2                 ! n*(n+1)/2, pairs including (i,i)
    INTEGER :: BlockSizeForDSYGVX  ! ILAENV block size ('G' only)

    ! -- eigensolver -------------------------------------------------
    INTEGER                                :: ErrorCode          ! DSYGVX INFO, or GSEPIIS status
    INTEGER                                :: NumOfEigvecs       ! eigenvectors asked of DSYGVX
    INTEGER                                :: NumOfEigvalsFound  ! DSYGVX M
    INTEGER                                :: NumOfIterations    ! inverse iterations used
    REAL(wp)                               :: Evalue             ! the eigenvalue, on every rank
    LOGICAL                                :: IsSwapFileOK       ! H and S came from the swap file
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: Eigvals
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: Eigvecs
    INTEGER, ALLOCATABLE, DIMENSION(:)     :: IFAIL

    ! -- correlation function ----------------------------------------
    LOGICAL  :: AreCorrFuncNeeded  ! DENSITIES with a grid file other than 'none'
    LOGICAL  :: IsFile1OK          ! the grid file was read without problems
    INTEGER  :: OpenFileErr        ! IOSTAT of the grid-file OPEN and READ
    INTEGER  :: NumCFGridPoints    ! number of grid points
    REAL(wp) :: temp1              ! one grid value while counting
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: CFGrid  ! grid points
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: CFkl    ! one basis pair
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: CF      ! the correlation function

    ! -- the accumulation --------------------------------------------
    INTEGER                             :: NumOfExpcVals  ! length of the packed value vector
    REAL(wp)                            :: factor         ! 2*c_i*c_j/sqrt(S_ii S_jj), or half of it on the diagonal
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: MEkl           ! one basis pair, packed
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: MEkl_s         ! this rank's partial sums
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: MEkl_r         ! the reduced sums

    ! -- one basis pair (kl) and the total (no suffix) ---------------
    REAL(wp)                               :: Skl, Tkl, Vkl, MVkl, Darwinkl, Darwin1kl, OOkl
    REAL(wp)                               :: H, S, T, V, MV, Darwin, Darwin1, OO
    ! Bra side of the pairs (i,j): Li, Ai = Li*Li' and M*Ai, built once per i
    REAL(wp) :: PreLk(Glob_n, Glob_n, 1), PreAk(Glob_n, Glob_n, 1), PreMAk(Glob_n, Glob_n, 1)
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: rmkl, rkl, r2kl, deltarkl
    REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: rm, r, r2, deltar
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: PrintedVals  ! every value of the screen block, for the field width
    INTEGER       :: wfld                                 ! field width of the value lines (NumFieldWidth)
    CHARACTER(40) :: fmte, fmtv, fmtp1, fmtp2, fmtd1, fmtd2  ! run-time formats built from wfld
    INTEGER                                   :: nsets, iq, wd  ! equivalent-pair sets, quantity, description width
    REAL(wp), ALLOCATABLE, DIMENSION(:, :)    :: SetAvg         ! (quantity, set): average over the equivalent pairs
    CHARACTER(128), ALLOCATABLE, DIMENSION(:) :: SetDescr       ! index-set notation of each set
    LOGICAL, ALLOCATABLE, DIMENSION(:)        :: SetSingle      ! pairs (0,i): label r_i instead of r_ij
    LOGICAL, ALLOCATABLE, DIMENSION(:)        :: SetHasZero     ! particle 0 listed: r_0j stands for r_j
    CHARACTER(11)                             :: lab            ! quantity label, e.g. delta(r_ij)
    CHARACTER(40)                             :: fmts           ! run-time format of the symmetrized lines

    ! -- symmetrized averages ----------------------------------------
    REAL(wp) :: beta  ! running sum over one equivalent-pair set


    IF (Glob_ProcID == 0) THEN
      WRITE(*, *)
      IF (Verbose >= 1) WRITE(*, *) 'Routine ExpectationValues started'
      IF (Verbose >= 1) WRITE(*, *) 'Number of basis functions', Glob_CurrBasisSize
      IF (Verbose >= 1) WRITE(*, *) 'GSEP solution method ', GSEPsolMethod
    ENDIF
    IF ((GSEPsolMethod /= 'G') .AND. (GSEPsolMethod /= 'I') .AND. (GSEPsolMethod /= 'Q')) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0195 in ExpectationValues: wrong GSEP solution method'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! STOP
    ENDIF
    IF (SymmAdaptMethod /= 1) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error EC0194 in ExpectationValues: SymmAdaptMethod must be 1'
        WRITE(*, *) 'ExpValuesMatElem takes a single symmetry matrix (Y''Y on the ket)'
      ENDIF
      CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! STOP
    ENDIF

    ! Setting the values of some global and local variables
    Glob_GSEPSolutionMethod = GSEPsolMethod
    Glob_OverlapPenaltyAllowed = .FALSE.
    Glob_HSLeadDim = Glob_CurrBasisSize
    n = Glob_n
    npt = Glob_npt
    np2 = n*(n+1)/2
    Glob_HSBuffLen = MAX(MIN(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2, 1000), 30*Glob_CurrBasisSize)
    cbs = Glob_CurrBasisSize
    IF (GSEPsolMethod == 'G') NumOfEigvecs = MIN(cbs, Glob_WhichEigenvalue+10)
    IF (GSEPsolMethod == 'I') NumOfEigvecs = 1
    IF (GSEPsolMethod == 'Q') NumOfEigvecs = 1


    !==================================================================
    ! Grid for the correlation function (DENSITIES only)
    !==================================================================
    ! As in the reference codes: the correlation function is needed when
    ! the action is DENSITIES and FileName1 is not 'none'. Only the
    ! nucleus-nucleus function (pair 0-1) is available; the particle
    ! densities are not, so a density grid (FileName3) other than 'none'
    ! is ignored with a warning.
    !------------------------------------------------------------------
    AreCorrFuncNeeded = .FALSE.

    IF (Action == 'DENSITIES') THEN
      IF (.NOT. ((FileName1 == ' ') .OR. (FileName1 == 'none') .OR. (FileName1 == 'NONE') .OR. &
                 (FileName1 == 'None'))) AreCorrFuncNeeded = .TRUE.
      IF (.NOT. ((FileName3 == ' ') .OR. (FileName3 == 'none') .OR. (FileName3 == 'NONE') .OR. &
                 (FileName3 == 'None'))) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'Warning WC0195 in ExpectationValues: particle densities are not available in PG_0S;'
          WRITE(*, *) 'the density grid file ', TRIM(FileName3), ' is ignored'
        ENDIF
      ENDIF
    ENDIF

    ! Here we determine the number of grid points for correlation
    ! function calculation (one value r >= 0 per line)
    NumCFGridPoints = 0
    IF (AreCorrFuncNeeded) THEN
      IF (Glob_ProcID == 0) THEN
        IsFile1OK = .TRUE.
        OPEN(1, FILE=FileName1, STATUS='old', IOSTAT=OpenFileErr)
        IF (OpenFileErr == 0) THEN
          NumCFGridPoints = 0
          temp1 = ZERO
          DO WHILE (OpenFileErr == 0)
            READ(1, *, IOSTAT=OpenFileErr) temp1
            IF (temp1 < ZERO) THEN
              OpenFileErr = 111
            ELSE
              NumCFGridPoints = NumCFGridPoints+1
            ENDIF
          ENDDO
          NumCFGridPoints = NumCFGridPoints-1
          IF ((OpenFileErr > 0) .OR. (NumCFGridPoints == 0)) THEN
            ! Improper data in the input file
            WRITE(*, *) 'Error EC0196 in ExpectationValues: improper data in file ', TRIM(FileName1)
            IsFile1OK = .FALSE.
          ENDIF
          CLOSE(1)
        ELSE
          IsFile1OK = .FALSE.
        ENDIF
      ENDIF
      CALL MPI_BCAST(IsFile1OK, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(NumCFGridPoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ! stop if there were problems with correlation function grid file
      IF (.NOT. IsFile1OK) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'Error EC0197 in ExpectationValues: cannot open CF grid file ', TRIM(FileName1)
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! stop
      ENDIF
    ENDIF

    ! Allocate the arrays for the grid points and the function values;
    ! one element when the function is not needed, to have valid arguments
    IF (AreCorrFuncNeeded) THEN
      ALLOCATE(CFGrid(NumCFGridPoints))
      ALLOCATE(CFkl(NumCFGridPoints))
      ALLOCATE(CF(NumCFGridPoints))
    ELSE
      ALLOCATE(CFGrid(1))
      ALLOCATE(CFkl(1))
    ENDIF

    ! Now we open file FileName1 again, but this time we read the data
    ! from it into array CFGrid
    IF (AreCorrFuncNeeded) THEN
      IF (Glob_ProcID == 0) THEN
        OPEN(1, FILE=FileName1, STATUS='old')
        DO i = 1, NumCFGridPoints
          READ(1, *) CFGrid(i)
        ENDDO
        CLOSE(1)
      ENDIF
      CALL MPI_BCAST(CFGrid, NumCFGridPoints, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF


    ! Allocate global arrays
    ALLOCATE(Glob_H(cbs, cbs))
    ALLOCATE(Glob_S(cbs, cbs))
    IF (GSEPsolMethod == 'G') ALLOCATE(Glob_diagH(cbs))
    ALLOCATE(Glob_diagS(cbs))
    IF (GSEPsolMethod == 'I') ALLOCATE(Glob_invD(cbs))
    ALLOCATE(Glob_c(cbs))
    ALLOCATE(Glob_HklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_HklBuff2(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff1(Glob_HSBuffLen))
    ALLOCATE(Glob_SklBuff2(Glob_HSBuffLen))

    ! Allocate workspace for DSYGVX
    IF (GSEPsolMethod == 'G') THEN
      BlockSizeForDSYGVX = ILAENV(1, 'DSYTRD', 'VIU', cbs, cbs, cbs, cbs)
      Glob_LWorkForDSYGVX = MAX((BlockSizeForDSYGVX+3)*cbs, 8*cbs)
      ALLOCATE(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      ALLOCATE(Glob_IWorkForDSYGVX(5*cbs))
    ENDIF

    ! Allocate workspace for subroutine GSEPIIS
    IF (GSEPsolMethod == 'I') THEN
      ALLOCATE(Glob_WorkForGSEPIIS(cbs))
      ALLOCATE(Glob_LastEigvector(cbs))
      Glob_LastEigvector(1:cbs) = ONE
    ENDIF

    ! Allocate local arrays
    IF (GSEPsolMethod == 'G') THEN
      ALLOCATE(Eigvals(NumOfEigvecs))
      ALLOCATE(Eigvecs(cbs, NumOfEigvecs))
      ALLOCATE(IFAIL(cbs))
    ENDIF

    ! Packed layout of MEkl (np2 = n*(n+1)/2 pairs, upper triangle a<=b):
    !  1/r       MEkl(        1 :   np2)
    !  r         MEkl(  np2 + 1 : 2*np2)
    !  r^2       MEkl(2*np2 + 1 : 3*np2)
    !  delta(r)  MEkl(3*np2 + 1 : 4*np2)
    !  S, T, V, MV, Darwin, Darwin1, OO
    !            MEkl(4*np2 + 1 : 4*np2 + 7)
    !  CF        MEkl(4*np2 + 8 : 4*np2 + 7 + NumCFGridPoints)  (only when needed)
    NumOfExpcVals = 4*np2+7
    IF (AreCorrFuncNeeded) NumOfExpcVals = NumOfExpcVals+NumCFGridPoints

    ALLOCATE(MEkl(NumOfExpcVals))
    ALLOCATE(MEkl_s(NumOfExpcVals))
    ALLOCATE(MEkl_r(NumOfExpcVals))

    ALLOCATE(rmkl(n, n))
    ALLOCATE(rm(n, n))

    ALLOCATE(rkl(n, n))
    ALLOCATE(r(n, n))

    ALLOCATE(r2kl(n, n))
    ALLOCATE(r2(n, n))

    ALLOCATE(deltarkl(n, n))
    ALLOCATE(deltar(n, n))

    CALL ReadSwapFileAndDistributeData(IsSwapFileOK)

    IF (.NOT. IsSwapFileOK) THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a52)', ADVANCE='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      CALL ComputeMatElem(1, cbs)
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, *) 'done'
    ENDIF

    IF (GSEPSolMethod == 'G') THEN
      DO i = 1, cbs
        DO j = 1, i-1
          Glob_H(j, i) = Glob_H(i, j)
        ENDDO
        Glob_H(i, i) = Glob_diagH(i)
      ENDDO
      DO i = 1, cbs
        DO j = 1, i-1
          Glob_S(j, i) = Glob_S(i, j)
        ENDDO
        Glob_S(i, i) = ONE
      ENDDO

      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
        CALL DSYGVX(1, 'V', 'I', 'U', cbs, Glob_H, Glob_HSLeadDim, Glob_S, Glob_HSLeadDim, &
                    ZERO, ZERO, 1, NumOfEigvecs, Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound, Eigvals, Eigvecs, cbs, Glob_WorkForDSYGVX, Glob_LWorkForDSYGVX, &
                    Glob_IWorkForDSYGVX, IFAIL, ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
                    !$      VL, VU, IL, IU, ABSTOL,
                    !$      M, W, Z, LDZ, WORK, LWORK,
                    !$      IWORK, IFAIL, INFO )
      ENDIF
      CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      IF (ErrorCode /= 0) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'failed'
          WRITE(*, *) &
            'Error EC0200 in ExpectationValues: routine DSYGVX failed with error code', ErrorCode
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! STOP
      ENDIF

      ! sending the eigenvalue and the eigenvector to all processes
      IF (Glob_ProcID == 0) THEN
        Evalue = Eigvals(Glob_WhichEigenvalue)
        Glob_c(1:cbs) = Eigvecs(1:cbs, Glob_WhichEigenvalue)
      ENDIF
      CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_BCAST(Glob_c, cbs, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      Glob_CurrEnergy = Evalue

      ! print the lower part of the spectrum
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) 'done'
        WRITE(*, *) 'Energy: ', Evalue
        WRITE(*, *)
        WRITE(*, *) 'Lowest eigenvalues:'
        WRITE(*, *) '------------------------------------------'
        WRITE(fmte, '(a,i0,a)') '(1x,i6,1x,f', NumFieldWidth(Eigvals(1:NumOfEigvalsFound)), '.16)'
        DO i = 1, NumOfEigvalsFound
          WRITE(*, fmte) i, Eigvals(i)
        ENDDO
        WRITE(*, *)
      ENDIF
    ENDIF  ! IF (GSEPSolMethod=='G')

    IF (GSEPSolMethod == 'I') THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      IF (cbs == 1) THEN
        ! Glob_H holds the SHIFTED matrix on this path (see StoreHS), so
        ! the shift is added back. Glob_diagH is not allocated for 'I'.
        Evalue = Glob_H(1, 1)+Glob_ApproxEnergy
        Glob_c(1) = ONE
        NumOfIterations = 1
        ErrorCode = 0
      ELSE
        CALL GSEPIIS(1, cbs, Glob_H, Glob_HSLeadDim, Glob_invD, Glob_S, Glob_HSLeadDim, &
                     Glob_ApproxEnergy, Glob_LastEigvector, Glob_WorkForGSEPIIS, Glob_EigvalTol, &
                     Evalue, Glob_c, Glob_LastEigvalTol, Glob_MaxIterForGSEPIIS, &
                     0, NumOfIterations, ErrorCode)
        ! GSEPIIS(k,n,M,nM,invD,B,nB, &
        !        apprlambda,v,w,Tol, &
        !        lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
        IF (Glob_LastEigvalTol > Glob_WorstEigvalTol) Glob_WorstEigvalTol = Glob_LastEigvalTol
        IF (Glob_LastEigvalTol > Glob_BestEigvalTol) Glob_BestEigvalTol = Glob_LastEigvalTol
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Evalue, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(Glob_c, cbs, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
      Glob_InvItTempCounter1 = Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2 = Glob_InvItTempCounter2+NumOfIterations
      Glob_CurrEnergy = Evalue
      IF (ErrorCode /= 0) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'failed'
          WRITE(*, *) 'Error EC0201 in ExpectationValues: the energy cannot be computed'
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! STOP
      ENDIF
      ! print the energy
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) 'done'
        WRITE(*, *) 'Energy: ', Evalue
      ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! 'Q': one fresh factorization of the canonical matrices and one solve
    !------------------------------------------------------------------
    IF (GSEPSolMethod == 'Q') THEN
      IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a29)', ADVANCE='no') 'Solving eigenvalue problem...'
      Glob_c = ONE
      CALL PrepareQWorkspace(cbs, cbs, 1, ErrorCode)
      Q_Workspace%MatricesAreCanonical = .TRUE.
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL FactorizeQFresh(ErrorCode)
      IF (ErrorCode == Q_METHOD_SUCCESS) CALL SolveQ(Evalue, ErrorCode)
      IF (ErrorCode /= Q_METHOD_SUCCESS) THEN
        IF (Glob_ProcID == 0) THEN
          WRITE(*, *) 'failed'
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0202 in ExpectationValues: Q energy cannot be computed, status', ErrorCode
        ENDIF
        CALL MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode)  ! STOP
      ENDIF
      Glob_CurrEnergy = Evalue
      IF (Glob_ProcID == 0) THEN
        IF (Verbose >= 2) WRITE(*, *) 'done'
        WRITE(*, *) 'Energy: ', Evalue
      ENDIF
    ENDIF

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 2)) WRITE(*, '(1x,a31)', ADVANCE='no') 'Computing expectation values...'


    !==================================================================
    ! Main loop over basis-function pairs
    !==================================================================
    ! The pairs are dealt out to the ranks round-robin. For each pair the
    ! Y'Y terms are summed with their coefficients, the result is
    ! normalized and weighted with the linear coefficients, and added to
    ! this rank's partial sums; one MPI_ALLREDUCE then combines them.
    !------------------------------------------------------------------
    MEkl_s(1:NumOfExpcVals) = ZERO
    counter = 0
    DO i = 1, cbs
      CALL PrecomputeMatrices_MA_L_A(npt, 1, Glob_NonlinParam(1:npt, i:i), &
                                     Glob_MassMatrix(1:n, 1:n), PreLk, PreAk, PreMAk)
      DO j = 1, i
        counter = counter+1
        IF (MOD(counter, Glob_NumOfProcs) == Glob_ProcID) THEN

          ! The matrix elements come back unnormalized; 1/sqrt(S_ii S_jj)
          ! normalizes them the way StoreHS does. Same range hazard and
          ! the same fallback as in matform.f90.
          factor = ONE/SQRT(Glob_diagS(i)*Glob_diagS(j))
          IF (.NOT. ((factor > ZERO) .AND. (factor <= HUGE(factor)))) &
            factor = ONE/(SQRT(Glob_diagS(i))*SQRT(Glob_diagS(j)))
          IF (i == j) THEN
            factor = factor*Glob_c(i)*Glob_c(j)
          ELSE
            factor = TWO*factor*Glob_c(i)*Glob_c(j)
          ENDIF

          DO k = 1, Glob_NumYHYTerms

            CALL ExpValuesMatElem(Glob_PWR(i), PreLk(:, :, 1), PreAk(:, :, 1), PreMAk(:, :, 1), &
                                  Glob_PWR(j), Glob_NonlinParam(1:npt, j), &
                                  Glob_YHYMatr(1:n, 1:n, k), &
                                  Skl, Tkl, Vkl, rmkl, rkl, r2kl, deltarkl, &
                                  MVkl, Darwinkl, Darwin1kl, OOkl, &
                                  NumCFGridPoints, CFGrid, CFkl, AreCorrFuncNeeded)

            c = 0
            DO a = 1, n
              DO b = a, n
                c = c+1; MEkl(c) = rmkl(a, b)
              ENDDO
            ENDDO
            DO a = 1, n
              DO b = a, n
                c = c+1; MEkl(c) = rkl(a, b)
              ENDDO
            ENDDO
            DO a = 1, n
              DO b = a, n
                c = c+1; MEkl(c) = r2kl(a, b)
              ENDDO
            ENDDO
            DO a = 1, n
              DO b = a, n
                c = c+1; MEkl(c) = deltarkl(a, b)
              ENDDO
            ENDDO
            c = c+1; MEkl(c) = Skl
            c = c+1; MEkl(c) = Tkl
            c = c+1; MEkl(c) = Vkl
            c = c+1; MEkl(c) = MVkl
            c = c+1; MEkl(c) = Darwinkl
            c = c+1; MEkl(c) = Darwin1kl
            c = c+1; MEkl(c) = OOkl
            IF (AreCorrFuncNeeded) THEN
              MEkl(c+1:c+NumCFGridPoints) = CFkl(1:NumCFGridPoints)
              c = c+NumCFGridPoints
            ENDIF

            MEkl_s(1:NumOfExpcVals) = MEkl_s(1:NumOfExpcVals) &
                                   +factor*Glob_YHYCoeff(k)*MEkl(1:NumOfExpcVals)

          ENDDO  ! k=1,Glob_NumYHYTerms

        ENDIF
      ENDDO
    ENDDO

    ! Combining the results of all processes
    CALL MPI_ALLREDUCE(MEkl_s, MEkl_r, NumOfExpcVals, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)


    ! Extracting expectation values from array MEkl_r
    c = 0
    DO a = 1, n
      DO b = a, n
        c = c+1
        rm(a, b) = MEkl_r(c); rm(b, a) = MEkl_r(c)
      ENDDO
    ENDDO
    DO a = 1, n
      DO b = a, n
        c = c+1
        r(a, b) = MEkl_r(c); r(b, a) = MEkl_r(c)
      ENDDO
    ENDDO
    DO a = 1, n
      DO b = a, n
        c = c+1
        r2(a, b) = MEkl_r(c); r2(b, a) = MEkl_r(c)
      ENDDO
    ENDDO
    DO a = 1, n
      DO b = a, n
        c = c+1
        deltar(a, b) = MEkl_r(c); deltar(b, a) = MEkl_r(c)
      ENDDO
    ENDDO
    c = c+1; S = MEkl_r(c)
    c = c+1; T = MEkl_r(c)
    c = c+1; V = MEkl_r(c)
    c = c+1; MV = MEkl_r(c)
    c = c+1; Darwin = MEkl_r(c)
    c = c+1; Darwin1 = MEkl_r(c)
    c = c+1; OO = MEkl_r(c)
    IF (AreCorrFuncNeeded) THEN
      CF(1:NumCFGridPoints) = MEkl_r(c+1:c+NumCFGridPoints)
      c = c+NumCFGridPoints
    ENDIF

    ! ExpValuesMatElem returns T and V separately; H is their sum.
    H = T+V


    !==================================================================
    ! Printing results
    !==================================================================
    IF (Glob_ProcID == 0) THEN

      ! Opening an additional file where selected expectation values will be saved
      OPEN(2, FILE=Glob_ExpValFileName, STATUS='replace')
      IF (Verbose >= 2) WRITE(*, *) 'done'
      IF (Verbose >= 2) WRITE(*, *)
      WRITE(*, *) 'Expectation values:'
      WRITE(*, *) '------------------------------------------'
      ! One field width for the whole block: 16 decimals, aligned on the
      ! decimal point, two blanks between '=' and the widest value
      ALLOCATE(PrintedVals(10+4*np2))
      PrintedVals(1:10) = (/ H, S, T, V, MV, Darwin, OO, MV*(Glob_FineStructConst**2), &
                             Darwin*(Glob_FineStructConst**2), OO*(Glob_FineStructConst**2) /)
      k = 10
      DO i = 1, n
        DO j = i, n
          PrintedVals(k+1:k+4) = (/ rm(i, j), r(i, j), r2(i, j), deltar(i, j) /)
          k = k+4
        ENDDO
      ENDDO
      wfld = NumFieldWidth(PrintedVals)
      DEALLOCATE(PrintedVals)
      WRITE(fmtv, '(a,i0,a)') '(1x,a,f', wfld, '.16)'
      WRITE(fmtp1, '(a,i0,a)') '(1x,a22,i1,a3,f', wfld, '.16)'
      WRITE(fmtp2, '(a,i0,a)') '(1x,a21,i1,i1,a3,f', wfld, '.16)'
      WRITE(fmtd1, '(a,i0,a)') '(1x,a21,i1,a1,a3,f', wfld, '.16)'
      WRITE(fmtd2, '(a,i0,a)') '(1x,a20,i1,i1,a1,a3,f', wfld, '.16)'
      WRITE(*, fmtv) '                      H = ', H
      WRITE(*, fmtv) '                      S = ', S
      WRITE(*, fmtv) '                      T = ', T
      WRITE(*, fmtv) '                      V = ', V
      WRITE(*, fmtv) '                     MV = ', MV
      WRITE(*, fmtv) '                 Darwin = ', Darwin
      ! Darwin1 is computed but, as in ExpcVals, not printed
      ! WRITE(*,*) '                Darwin1=',Darwin1
      WRITE(*, fmtv) '                     OO = ', OO
      WRITE(*, fmtv) '           (alpha^2)*MV = ', MV*(Glob_FineStructConst**2)
      WRITE(*, fmtv) '       (alpha^2)*Darwin = ', Darwin*(Glob_FineStructConst**2)
      WRITE(*, fmtv) '           (alpha^2)*OO = ', OO*(Glob_FineStructConst**2)
      WRITE(*, *)

      WRITE(2, '(a)', ADVANCE='no') '                  basis '
      WRITE(2, *) cbs
      WRITE(2, '(a)', ADVANCE='no') '                 Energy '
      CALL writerealadv(2, Evalue)
      WRITE(2, '(a)', ADVANCE='no') '                      H '
      CALL writerealadv(2, H)
      WRITE(2, '(a)', ADVANCE='no') '                      S '
      CALL writerealadv(2, S)
      WRITE(2, '(a)', ADVANCE='no') '                      T '
      CALL writerealadv(2, T)
      WRITE(2, '(a)', ADVANCE='no') '                      V '
      CALL writerealadv(2, V)
      WRITE(2, '(a)', ADVANCE='no') '                     MV '
      CALL writerealadv(2, MV)
      WRITE(2, '(a)', ADVANCE='no') '                 Darwin '
      CALL writerealadv(2, Darwin)
      WRITE(2, '(a)', ADVANCE='no') '                     OO '
      CALL writerealadv(2, OO)
      WRITE(2, '(a)', ADVANCE='no') '           (alpha^2)*MV '
      CALL writerealadv(2, MV*(Glob_FineStructConst**2))
      WRITE(2, '(a)', ADVANCE='no') '       (alpha^2)*Darwin '
      CALL writerealadv(2, Darwin*(Glob_FineStructConst**2))
      WRITE(2, '(a)', ADVANCE='no') '           (alpha^2)*OO '
      CALL writerealadv(2, OO*(Glob_FineStructConst**2))

      IF ((Glob_NumOfIdentPartSets /= Glob_n+1) .AND. (SymmAdaptMethod == 1)) THEN
        IF (Verbose >= 1) WRITE(*, *) '(Warning! These values do not account for indistinguishability of'
        WRITE(*, *) 'identical particles and other possible symmetries of the system)'
        WRITE(*, *)
      ENDIF
      DO i = 1, n
        WRITE(*, fmtp1) '                  1/r_', i, ' = ', rm(i, i)
        DO j = i+1, n
          WRITE(*, fmtp2) '                 1/r_', i, j, ' = ', rm(i, j)
        ENDDO
      ENDDO
      IF (Verbose >= 1) WRITE(*, *)
      DO i = 1, n
        WRITE(*, fmtp1) '                    r_', i, ' = ', r(i, i)
        DO j = i+1, n
          WRITE(*, fmtp2) '                   r_', i, j, ' = ', r(i, j)
        ENDDO
      ENDDO
      IF (Verbose >= 1) WRITE(*, *)
      DO i = 1, n
        WRITE(*, fmtp1) '                  r^2_', i, ' = ', r2(i, i)
        DO j = i+1, n
          WRITE(*, fmtp2) '                 r^2_', i, j, ' = ', r2(i, j)
        ENDDO
      ENDDO
      IF (Verbose >= 1) WRITE(*, *)
      DO i = 1, n
        WRITE(*, fmtd1) '            delta(r_', i, ')', ' = ', deltar(i, i)
        DO j = i+1, n
          WRITE(*, fmtd2) '            delta(r_', i, j, ')', ' = ', deltar(i, j)
        ENDDO
      ENDDO
      IF (Verbose >= 1) WRITE(*, *)

      IF (Glob_NumOfIdentPartSets /= Glob_n+1) THEN
        WRITE(*, *) 'Based on the particle mass and charge values it was determined'
        WRITE(*, *) 'that the system has the following sets of identical particles:'
        DO i = 1, Glob_NumOfIdentPartSets
          j = Glob_NumOfPartInIdentPartSet(i)
          WRITE(*, '(1x,a3,i2,a13)', ADVANCE='no') 'set', i, ' :  particles'
          ! four blanks before the first number, two between the others
          WRITE(fmts, '(a,i0,a)') '(2x,', j, '(2x,i0))'
          WRITE(*, fmts) Glob_IdentPartList(1:j, i)
        ENDDO
        WRITE(*, *)
        WRITE(*, *)
        WRITE(*, *) 'Properly symmetrized expectation values :'
        WRITE(*, *) '------------------------------------------'
        WRITE(*, *) 'Properly symmetrized expectation values of two-particle quantities'
        WRITE(*, *) 'that account for permutational symmetry of the above mentioned sets'
        WRITE(*, *) 'of identical particles are:'
        IF (Verbose >= 1) WRITE(*, *) '(Warning! An additional symmetrization might be necessary if the'
        WRITE(*, *) 'Young operator contains other types of symmetries)'
        WRITE(*, *)
        !--------------------------------------------------------------
        ! Averages over the equivalent pairs of every set, then the
        ! screen block in index-set notation (PairSetDescription) with
        ! the values aligned on the decimal point, then the file lines
        ! (label of the first pair of each set, as before).
        !--------------------------------------------------------------
        nsets = Glob_NumOfNoneqvPairSets
        ALLOCATE(SetAvg(4, nsets))
        ALLOCATE(SetDescr(nsets))
        ALLOCATE(SetSingle(nsets))
        ALLOCATE(SetHasZero(nsets))
        DO i = 1, nsets
          k = Glob_NumOfPairsInEqvPairSet(i)
          SetAvg(:, i) = ZERO
          DO j = 1, k
            a = Glob_EqvPairList(1, j, i)
            b = Glob_EqvPairList(2, j, i)
            SetAvg(1, i) = SetAvg(1, i)+rm(a, b)
            SetAvg(2, i) = SetAvg(2, i)+r(a, b)
            SetAvg(3, i) = SetAvg(3, i)+r2(a, b)
            SetAvg(4, i) = SetAvg(4, i)+deltar(a, b)
          ENDDO
          SetAvg(:, i) = SetAvg(:, i)/k
          CALL PairSetDescription(i, SetDescr(i), SetSingle(i), SetHasZero(i))
        ENDDO
        wd = MAX(15, MAXVAL(LEN_TRIM(SetDescr)))
        wfld = NumFieldWidth(RESHAPE(SetAvg, (/ 4*nsets /)))
        WRITE(fmts, '(a,i0,a)') '(1x,a,4x,a11,a3,f', wfld, '.16)'
        WRITE(*, '(1x,a)') 'Each value is the average over the equivalent pairs listed in front of it;'
        WRITE(*, '(1x,a)') 'r_i is the distance of particle i from particle 0, r_ij that between i and j'
        IF (ANY(SetHasZero)) WRITE(*, '(1x,a)') '(r_0j stands for r_j)'
        WRITE(*, *)
        DO iq = 1, 4
          DO i = 1, nsets
            SELECT CASE (iq)
            CASE (1)
              lab = '1/r_'
            CASE (2)
              lab = 'r_'
            CASE (3)
              lab = 'r^2_'
            CASE (4)
              lab = 'delta(r_'
            END SELECT
            IF (SetSingle(i)) THEN
              lab = TRIM(lab)//'i'
            ELSE
              lab = TRIM(lab)//'ij'
            ENDIF
            IF (iq == 4) lab = TRIM(lab)//')'
            WRITE(*, fmts) SetDescr(i)(1:wd), ADJUSTR(lab), ' = ', SetAvg(iq, i)
          ENDDO
          WRITE(*, *)
        ENDDO
        ! write to file
        DO i = 1, nsets
          a = Glob_EqvPairList(1, 1, i)
          b = Glob_EqvPairList(2, 1, i)
          IF (a /= b) WRITE(2, '(a,i1,i1,1x)', ADVANCE='no') '                 1/r_', a, b
          IF (a == b) WRITE(2, '(a,i1,1x)', ADVANCE='no') '                  1/r_', a
          CALL writerealadv(2, SetAvg(1, i))
        ENDDO
        DO i = 1, nsets
          a = Glob_EqvPairList(1, 1, i)
          b = Glob_EqvPairList(2, 1, i)
          IF (a /= b) WRITE(2, '(a,i1,i1,1x)', ADVANCE='no') '                   r_', a, b
          IF (a == b) WRITE(2, '(a,i1,1x)', ADVANCE='no') '                    r_', a
          CALL writerealadv(2, SetAvg(2, i))
        ENDDO
        DO i = 1, nsets
          a = Glob_EqvPairList(1, 1, i)
          b = Glob_EqvPairList(2, 1, i)
          IF (a /= b) WRITE(2, '(a,i1,i1,1x)', ADVANCE='no') '                 r^2_', a, b
          IF (a == b) WRITE(2, '(a,i1,1x)', ADVANCE='no') '                  r^2_', a
          CALL writerealadv(2, SetAvg(3, i))
        ENDDO
        DO i = 1, nsets
          a = Glob_EqvPairList(1, 1, i)
          b = Glob_EqvPairList(2, 1, i)
          IF (a /= b) WRITE(2, '(a,i1,i1,a1,1x)', ADVANCE='no') '            delta(r_', a, b, ')'
          IF (a == b) WRITE(2, '(a,i1,a1,1x)', ADVANCE='no') '             delta(r_', a, ')'
          CALL writerealadv(2, SetAvg(4, i))
        ENDDO
        DEALLOCATE(SetHasZero)
        DEALLOCATE(SetSingle)
        DEALLOCATE(SetDescr)
        DEALLOCATE(SetAvg)
      ENDIF

      CLOSE(2)

      ! Saving the correlation function in the layout of the reference
      ! codes: the grid column r and g1, the nucleus-nucleus function
      ! (pair 0-1), which is the only one available here
      IF (AreCorrFuncNeeded) THEN
        OPEN(1, FILE=FileName2, STATUS='replace')
        ! first we print titles of the data columns
        WRITE(1, '(10x,a2,1x)', ADVANCE='no') '#r'
        WRITE(1, '(21x,a1,i1,1x)', ADVANCE='no') 'g', 1
        WRITE(1, *)
        ! then we print the data columns themselves
        DO k = 1, NumCFGridPoints
          WRITE(1, '(1x,e23.16)', ADVANCE='no') CFGrid(k)
          WRITE(1, '(1x,e23.16)', ADVANCE='no') CF(k)
          WRITE(1, *)
        ENDDO
        CLOSE(1)
        WRITE(*, *) 'Correlation function has been stored in file ', TRIM(FileName2)
        WRITE(*, *)
      ENDIF

    ENDIF


    ! deallocate local arrays

    DEALLOCATE(rmkl)
    DEALLOCATE(rm)

    DEALLOCATE(rkl)
    DEALLOCATE(r)

    DEALLOCATE(r2kl)
    DEALLOCATE(r2)

    DEALLOCATE(deltarkl)
    DEALLOCATE(deltar)

    DEALLOCATE(CFGrid)
    DEALLOCATE(CFkl)
    IF (AreCorrFuncNeeded) DEALLOCATE(CF)

    DEALLOCATE(MEkl)
    DEALLOCATE(MEkl_s)
    DEALLOCATE(MEkl_r)

    IF (GSEPSolMethod == 'G') THEN
      DEALLOCATE(IFAIL)
      DEALLOCATE(Eigvecs)
      DEALLOCATE(Eigvals)
    ENDIF

    IF (GSEPsolMethod == 'I') THEN
      DEALLOCATE(Glob_LastEigvector)
      DEALLOCATE(Glob_WorkForGSEPIIS)
    ENDIF
    IF (GSEPSolMethod == 'Q') CALL ClearQWorkspace()

    ! deallocate workspace for DSYGVX
    IF (GSEPSolMethod == 'G') THEN
      DEALLOCATE(Glob_WorkForDSYGVX)
      DEALLOCATE(Glob_IWorkForDSYGVX)
    ENDIF

    ! deallocate global arrays
    DEALLOCATE(Glob_SklBuff2)
    DEALLOCATE(Glob_SklBuff1)
    DEALLOCATE(Glob_HklBuff2)
    DEALLOCATE(Glob_HklBuff1)
    DEALLOCATE(Glob_c)
    IF (GSEPSolMethod == 'I') DEALLOCATE(Glob_invD)
    DEALLOCATE(Glob_diagS)
    IF (GSEPSolMethod == 'G') DEALLOCATE(Glob_diagH)
    DEALLOCATE(Glob_S)
    DEALLOCATE(Glob_H)

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 1)) WRITE (*, *) 'Routine ExpectationValues finished'

  END SUBROUTINE ExpectationValues


END MODULE workproc
