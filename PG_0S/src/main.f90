PROGRAM main
  !==================================================================
  ! Program main
  !==================================================================
  ! Driver of the Basis Building and Optimization Program (BBOP): reads
  ! the data file, initializes the program data, seeds the random
  ! generators, executes the BBOP steps in order and empties the swap
  ! file. Besides the PG_0S frame it provides ECG_RND_SEED seeding, the
  ! ECG_EIG_IDX_TARGETING override, recovery of the last physical energy
  ! before the inverse-iteration shift, the SAVE_HS_R and OVERLAP_D steps
  ! and the FULL_OPT1 history budget test. DENSITIES and MOMT_DENS are
  ! not available and stay commented out.
  !==================================================================

  USE workproc

  IMPLICIT NONE

  !------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------
  INTEGER :: i                     ! BBOP step counter
  INTEGER :: iw                    ! trimmed file-name length
  INTEGER :: k                     ! basis size, in the energy recovery loop
  INTEGER :: Kstart, Kstop, Kstep  ! BASIS_ENL range and block size
  INTEGER :: OpenFileErr           ! IOSTAT of the swap-file OPEN
  INTEGER :: OptimizationType      ! 1 = optimize the new functions, see BasisEnlG/I/Q
  INTEGER :: QErrorCode            ! status returned by the drivers of the Q method
  REAL(8) :: r8                    ! scratch for the random generators

  ! Random-generator seeding
  INTEGER                            :: RNSeedSize
  INTEGER, ALLOCATABLE, DIMENSION(:) :: Seed

  ! Deterministic seeding (environment variable ECG_RND_SEED)
  CHARACTER(LEN=32) :: SeedEnv
  INTEGER           :: SeedVal, SeedStat, SeedErr

  ! Eigenvalue-index targeting override (environment variable ECG_EIG_IDX_TARGETING)
  CHARACTER(LEN=32) :: EigIdxEnv
  INTEGER           :: EigIdxVal, EigIdxStat, EigIdxErr


  !==================================================================
  ! MPI
  !==================================================================
  CALL MPI_INIT(Glob_MPIErrCode)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Glob_ProcID, Glob_MPIErrCode)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Glob_NumOfProcs, Glob_MPIErrCode)

  IF (Glob_ProcID == 0) THEN
    WRITE(*, '(1x,a,1x,a,1x,a)') 'Program', Glob_BasisType, 'started'
    WRITE(*, '(1x,a,1x,i0)') 'Number of parallel MPI processes running:', Glob_NumOfProcs
    WRITE(*, *)
  ENDIF


  !==================================================================
  ! Input and program data
  !==================================================================
  CALL ReadIOFile()
  IF (Glob_IsOptCycleScripted) CALL ReadBlackList()
  CALL ProgramDataInit()


  !==================================================================
  ! Random number generators
  !==================================================================
  ! Seeded from the clock unless the environment variable ECG_RND_SEED
  ! holds an integer, which makes runs reproducible for a given process
  ! count. The seed is offset by the rank so that processes drawing their
  ! own candidates do not test identical functions. DRNOR_START needs a
  ! nonzero seed, hence the two substitutions.
  !------------------------------------------------------------------
  CALL GET_ENVIRONMENT_VARIABLE('ECG_RND_SEED', SeedEnv, STATUS=SeedStat)
  SeedErr = 1
  IF (SeedStat == 0) READ(SeedEnv, *, IOSTAT=SeedErr) SeedVal

  CALL RANDOM_SEED(SIZE=RNSeedSize)
  ALLOCATE(Seed(RNSeedSize))

  IF (SeedErr == 0) THEN

    IF (SeedVal == 0) SeedVal = 519
    DO i = 1, RNSeedSize
      Seed(i) = SeedVal+37*(i-1)+7919*Glob_ProcID
    ENDDO
    CALL RANDOM_SEED(PUT=Seed(1:RNSeedSize))
    r8 = drnor_start(SeedVal+7919*Glob_ProcID)
    IF (Glob_ProcID == 0) THEN
      WRITE(*, *) 'Random generators seeded from ECG_RND_SEED = ', SeedVal
      WRITE(*, *)
    ENDIF

  ELSE

    CALL RANDOM_SEED()
    CALL RANDOM_SEED(GET=Seed(1:RNSeedSize))
    CALL SYSTEM_CLOCK(COUNT=Seed(1))
    Seed(1:RNSeedSize) = Seed(1:RNSeedSize)+Glob_ProcID
    CALL RANDOM_SEED(PUT=Seed(1:RNSeedSize))
    CALL RANDOM_NUMBER(r8)
    CALL RANDOM_NUMBER(r8)
    CALL RANDOM_NUMBER(r8)
    CALL RANDOM_NUMBER(r8)
    SeedVal = NINT(r8*25000)+Glob_ProcID
    IF (SeedVal == 0) SeedVal = 519
    r8 = drnor_start(SeedVal)

  ENDIF

  DEALLOCATE(Seed)


  !==================================================================
  ! Eigenvalue-index targeting for the inverse-iteration path
  !==================================================================
  ! On by default (Glob_EigIdxTargeting=1 in linalg): the data file
  ! states WHICH_EIGENVALUE and the DSYGVX ('G') path has always
  ! honoured it, so an 'I' step that silently tracked whatever level
  ! happened to sit nearest the shift would make 'G' and 'I' steps of
  ! the same script optimize different states. Set
  ! ECG_EIG_IDX_TARGETING=0 to get the index-blind behaviour back.
  ! See RetargetShiftToEigenvalue and IsRequestedEigenstate.
  !------------------------------------------------------------------
  CALL GET_ENVIRONMENT_VARIABLE('ECG_EIG_IDX_TARGETING', EigIdxEnv, STATUS=EigIdxStat)
  EigIdxErr = 1
  IF (EigIdxStat == 0) READ(EigIdxEnv, *, IOSTAT=EigIdxErr) EigIdxVal
  IF (EigIdxErr == 0) THEN
    IF (EigIdxVal == 0) THEN
      Glob_EigIdxTargeting = 0
    ELSE
      Glob_EigIdxTargeting = 1
    ENDIF
    IF (Glob_ProcID == 0) THEN
      WRITE(*, *) 'Eigenvalue-index targeting set from ECG_EIG_IDX_TARGETING = ', &
                 Glob_EigIdxTargeting
      WRITE(*, *)
    ENDIF
  ENDIF


  !==================================================================
  ! Swap file
  !==================================================================
  ! Emptied here because it may still hold garbage left by the last
  ! run, if that run failed.
  !------------------------------------------------------------------
  IF ((Glob_ProcID == 0) .AND. Glob_UseSwapFile) THEN
    OPEN(1, FILE=Glob_SwapFileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)
    IF (OpenFileErr == 0) THEN
      WRITE(1) 'Swap file is empty'
      CLOSE(1)
    ENDIF
  ENDIF

  ! 1 = optimize the newly selected functions (the premultiplier powers,
  ! then the nonlinear parameters with DRMNG); 0 would keep the randomly
  ! selected block as it is - see the SELECT CASE in BasisEnlG/BasisEnlI.
  ! NOTE this is not what OptimizationType meant in the old workproc.f90,
  ! where anything other than 1 selected the DMNG optimization.
  OptimizationType = 1


  !==================================================================
  ! BBOP steps
  !==================================================================
  DO i = 1, Glob_NumOfBBOPSteps

    ! The parallel/serial mode of each linear algebra routine is selected
    ! through linalg_setparam, called whenever the problem size changes.
    CALL linalg_setparam(Glob_CurrBasisSize)

    Glob_CurrBBOPStep = i

    !----------------------------------------------------------------
    ! Recover the last physical energy before forming the shift
    !----------------------------------------------------------------
    ! Glob_CurrEnergy may still hold the rejection sentinel (1e31) after a
    ! step in which every trial was rejected; used as the inverse-iteration
    ! shift it makes the solver return garbage (multiples of 2**39).
    !----------------------------------------------------------------
    IF (ABS(Glob_CurrEnergy) > 1.0E10_wp) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        WRITE(*, *) '*** WARNING: the previous step left a non-physical energy ***'
        WRITE(*, *) 'Glob_CurrEnergy = ', Glob_CurrEnergy
      ENDIF
      DO k = Glob_CurrBasisSize, 1, -1
        IF ((ABS(Glob_History(k)%Energy) < 1.0E10_wp) .AND. (Glob_History(k)%Energy /= ZERO)) THEN
          Glob_CurrEnergy = Glob_History(k)%Energy
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *) 'Recovered the energy recorded at basis size ', k
            WRITE(*, *) 'Glob_CurrEnergy = ', Glob_CurrEnergy
            WRITE(*, *)
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    Glob_ApproxEnergy = Glob_CurrEnergy*Glob_InvItParameter


    SELECT CASE (Glob_BBOP(i)%Action)

    !----------------------------------------------------------------
    ! BASIS_ENL  Method  Kstart Kstop Kstep NTrials MaxEnergyEval Q R
    !----------------------------------------------------------------
    CASE ('BASIS_ENL')
      Kstart = Glob_BBOP(i)%A
      Kstop = Glob_BBOP(i)%B
      Kstep = Glob_BBOP(i)%C
      IF (Kstop > Glob_CurrBasisSize) THEN
        IF (Kstart <= Glob_CurrBasisSize+1) THEN
          Kstart = Glob_CurrBasisSize+1
          SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
          CASE ('G')
            CALL BasisEnlG(Kstart, Kstop, Kstep, Glob_BBOP(i)%D, OptimizationType, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%Q, Glob_BBOP(i)%R)
          CASE ('I')
            CALL BasisEnlI(Kstart, Kstop, Kstep, Glob_BBOP(i)%D, OptimizationType, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%Q, Glob_BBOP(i)%R)
          CASE ('Q')
            CALL BasisEnlQ(Kstart, Kstop, Kstep, Glob_BBOP(i)%D, OptimizationType, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, QErrorCode)
            IF ((QErrorCode /= Q_METHOD_SUCCESS) .AND. (Glob_ProcID == 0)) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,i0)') 'Error EC0014 in main: BASIS_ENL Q failed at BBOP step', i, &
                'with status', QErrorCode
            ENDIF
          CASE DEFAULT
            IF (Glob_ProcID == 0) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: BASIS_ENL step', i, &
                'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
              WRITE(*, *) 'Only G, I and Q are recognized. Skipping this step...'
            ENDIF
          END SELECT
        ELSE
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0)') 'Error EC0001 in main: incorrect BBOP step', i
            WRITE(*, *) 'One or more parameters in BASIS_ENL are incorrect'
          ENDIF
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! OPT_CYCLE  Method  BasisSize FuncBegin FuncEnd NumOfFuncToOpt
    !            NumOfFuncToShift NumCycles MaxEnergyEval Q R SavingFreq
    !----------------------------------------------------------------
    CASE ('OPT_CYCLE')
      IF (Glob_ProcID == 0) THEN
        IF ((Glob_BBOP(i)%C > Glob_BBOP(i)%A) .OR. (Glob_BBOP(i)%B > Glob_BBOP(i)%C)) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0002 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in OPT_CYCLE are incorrect'
        ENDIF
      ENDIF
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%B > 0) .AND. &
          (Glob_BBOP(i)%C <= Glob_CurrBasisSize)) THEN
        IF (Glob_History(Glob_CurrBasisSize)%CyclesDone < Glob_BBOP(i)%F) THEN
          SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
          CASE ('G')
            CALL OptCycleG(Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%F, Glob_BBOP(i)%G, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          CASE ('I')
            CALL OptCycleI(Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%F, Glob_BBOP(i)%G, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          CASE ('Q')
            CALL OptCycleQ(Glob_BBOP(i)%A, Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, &
                           Glob_BBOP(i)%E, Glob_BBOP(i)%F, Glob_BBOP(i)%G, Glob_BBOP(i)%Q, Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H, QErrorCode)
            IF ((QErrorCode /= Q_METHOD_SUCCESS) .AND. (Glob_ProcID == 0)) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,i0)') 'Error EC0015 in main: OPT_CYCLE Q failed at BBOP step', i, &
                'with status', QErrorCode
            ENDIF
          CASE DEFAULT
            IF (Glob_ProcID == 0) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: OPT_CYCLE step', i, &
                'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
              WRITE(*, *) 'Only G, I and Q are recognized. Skipping this step...'
            ENDIF
          END SELECT
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! FULL_OPT1  Method  BasisSize InitFunc FinalFunc MaxEnergyEval Q R
    !            DataSaveInterval HessianSaveInterval HessianFile
    !----------------------------------------------------------------
    CASE ('FULL_OPT1')
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%B > 0) .AND. &
          (Glob_BBOP(i)%C <= Glob_CurrBasisSize)) THEN
        ! A restarted run does not repeat a full optimization whose
        ! energy-evaluation budget (field D) the history already shows
        ! as spent.
        IF (Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt < Glob_BBOP(i)%D) THEN
          SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
          CASE ('G')
            CALL FullOpt1G(Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, Glob_BBOP(i)%Q, &
                           Glob_BBOP(i)%R, REAL(Glob_BBOP(i)%E, 4), REAL(Glob_BBOP(i)%F, 4), &
                           Glob_BBOP(i)%FileName1)
          CASE ('I')
            CALL FullOpt1I(Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, Glob_BBOP(i)%Q, &
                           Glob_BBOP(i)%R, REAL(Glob_BBOP(i)%E, 4), REAL(Glob_BBOP(i)%F, 4), &
                           Glob_BBOP(i)%FileName1)
          CASE ('Q')
            CALL FullOpt1Q(Glob_BBOP(i)%B, Glob_BBOP(i)%C, Glob_BBOP(i)%D, Glob_BBOP(i)%Q, &
                           Glob_BBOP(i)%R, REAL(Glob_BBOP(i)%E, 4), REAL(Glob_BBOP(i)%F, 4), &
                           Glob_BBOP(i)%FileName1, QErrorCode)
            IF ((QErrorCode /= Q_METHOD_SUCCESS) .AND. (Glob_ProcID == 0)) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,i0)') 'Error EC0016 in main: FULL_OPT1 Q failed at BBOP step', i, &
                'with status', QErrorCode
            ENDIF
          CASE DEFAULT
            IF (Glob_ProcID == 0) THEN
              WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: FULL_OPT1 step', i, &
                'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
              WRITE(*, *) 'Only G, I and Q are recognized. Skipping this step...'
            ENDIF
          END SELECT
        ELSE
          IF (Glob_ProcID == 0) THEN
            WRITE(*, *)
            WRITE(*, '(1x,a,1x,i0,1x,a)') 'FULL_OPT1 step', i, &
              'skipped: its energy-evaluation budget was already spent (see the history)'
          ENDIF
        ENDIF
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0003 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in FULL_OPT1 are incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! ELIM_LCFN  Method  BasisSize Threshold OutFile
    !----------------------------------------------------------------
    ! Drops every basis function whose linear coefficient is smaller in
    ! magnitude than Q and writes the reduced basis to FileName1. The
    ! routine solves with DSYGVX ('G') or with the QR method ('Q');
    ! inverse iteration is not available in the four cleanup steps.
    !----------------------------------------------------------------
    CASE ('ELIM_LCFN')
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%Q > ZERO)) THEN
        SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
        CASE ('G')
          CALL EliminateLittleContribFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1, &
                                          Glob_ElimRoutPrintSpec, 'G')
        CASE ('Q')
          CALL EliminateLittleContribFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1, &
                                          Glob_ElimRoutPrintSpec, 'Q')
        CASE ('I')
          IF (Glob_ProcID == 0) WRITE(*, *) 'Sorry, GSEP solution method I does not work in ELIM_LCFN'
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: ELIM_LCFN step', i, &
              'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
            WRITE(*, *) 'Only G and Q are recognized. Skipping this step...'
          ENDIF
        END SELECT
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0004 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in ELIM_LCFN are incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! ELIM_LND1  Method  BasisSize Threshold OutFile
    !----------------------------------------------------------------
    CASE ('ELIM_LND1')
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%Q > ZERO)) THEN
        SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
        CASE ('G')
          CALL EliminateLinDepFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1, Glob_ElimRoutPrintSpec, 'G')
        CASE ('Q')
          CALL EliminateLinDepFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%FileName1, Glob_ElimRoutPrintSpec, 'Q')
        CASE ('I')
          IF (Glob_ProcID == 0) WRITE(*, *) 'Sorry, GSEP solution method I does not work in ELIM_LND1'
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: ELIM_LND1 step', i, &
              'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
            WRITE(*, *) 'Only G and Q are recognized. Skipping this step...'
          ENDIF
        END SELECT
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0005 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in ELIM_LND1 are incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! SEPR_LND1  Method  BasisSize Threshold SeparationParam OutFile
    !----------------------------------------------------------------
    CASE ('SEPR_LND1')
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%Q > ZERO)) THEN
        SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
        CASE ('G')
          CALL SeparateLinDepFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%R, Glob_BBOP(i)%FileName1, &
                                  Glob_ElimRoutPrintSpec, 'G')
        CASE ('Q')
          CALL SeparateLinDepFunc(Glob_BBOP(i)%Q, Glob_BBOP(i)%R, Glob_BBOP(i)%FileName1, &
                                  Glob_ElimRoutPrintSpec, 'Q')
        CASE ('I')
          IF (Glob_ProcID == 0) WRITE(*, *) 'Sorry, GSEP solution method I does not work in SEPR_LND1'
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: SEPR_LND1 step', i, &
              'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
            WRITE(*, *) 'Only G and Q are recognized. Skipping this step...'
          ENDIF
        END SELECT
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0006 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in SEPR_LND1 are incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! SEPR_FLCF  Method  BasisSize Threshold SeparationParam OutFile
    !----------------------------------------------------------------
    CASE ('SEPR_FLCF')
      IF ((Glob_BBOP(i)%A == Glob_CurrBasisSize) .AND. (Glob_BBOP(i)%Q > ZERO)) THEN
        SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
        CASE ('G')
          CALL SeparateFuncLargeCoeff(Glob_BBOP(i)%Q, Glob_BBOP(i)%R, Glob_BBOP(i)%FileName1, &
                                      Glob_ElimRoutPrintSpec, 'G')
        CASE ('Q')
          CALL SeparateFuncLargeCoeff(Glob_BBOP(i)%Q, Glob_BBOP(i)%R, Glob_BBOP(i)%FileName1, &
                                      Glob_ElimRoutPrintSpec, 'Q')
        CASE ('I')
          IF (Glob_ProcID == 0) WRITE(*, *) 'Sorry, GSEP solution method I does not work in SEPR_FLCF'
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: SEPR_FLCF step', i, &
              'has GSEP solution method', Glob_BBOP(i)%GSEPSolutionMethod
            WRITE(*, *) 'Only G and Q are recognized. Skipping this step...'
          ENDIF
        END SELECT
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0007 in main: incorrect BBOP step', i
          WRITE(*, *) 'One or more parameters in SEPR_FLCF are incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! EXPC_VALS  Method  BasisSize
    !----------------------------------------------------------------
    ! The G, I and Q solution methods are supported; ExpectationValues
    ! keeps the paths completely separate and stops on any other method.
    !----------------------------------------------------------------
    CASE ('EXPC_VALS')
      IF (Glob_BBOP(i)%A == Glob_CurrBasisSize) THEN
        CALL ExpectationValues(Glob_BBOP(i)%Action, 1, Glob_FileNameNone, Glob_FileNameNone, Glob_FileNameNone, &
                               Glob_FileNameNone, Glob_BBOP(i)%GSEPSolutionMethod)
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0008 in main: incorrect BBOP step', i
          WRITE(*, *) 'Second parameter in EXPC_VALS is incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! DENSITIES and MOMT_DENS - not available in this frame: the density
    ! outputs were removed from ExpectationValues and ReadIOFile rejects
    ! the two actions. The PG_0S dispatch is kept below for reference.
    !----------------------------------------------------------------
    ! CASE('DENSITIES')
    !  IF (Glob_BBOP(i)%A==Glob_CurrBasisSize) THEN
    !    CALL ExpectationValues(Glob_BBOP(i)%Action,1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
    !                           Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4,Glob_BBOP(i)%GSEPSolutionMethod)
    !  ELSE
    !    IF (Glob_ProcID==0) THEN
    !      WRITE(*,'(1x,a,1x,i0)') 'Error EC0009 in main: incorrect BBOP step',i
    !      WRITE(*,*) 'Second parameter in DENSITIES is incorrect'
    !    ENDIF
    !  ENDIF
    ! CASE('MOMT_DENS')
    !  IF (Glob_BBOP(i)%A==Glob_CurrBasisSize) THEN
    !    CALL ExpectationValues(Glob_BBOP(i)%Action,1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
    !                           Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4,Glob_BBOP(i)%GSEPSolutionMethod)
    !  ELSE
    !    IF (Glob_ProcID==0) THEN
    !      WRITE(*,'(1x,a,1x,i0)') 'Error EC0010 in main: incorrect BBOP step',i
    !      WRITE(*,*) 'Second parameter in MOMT_DENS is incorrect'
    !    ENDIF
    !  ENDIF

    !----------------------------------------------------------------
    ! SAVE_FILE  BasisSize FileName   (no solution-method field)
    !----------------------------------------------------------------
    ! Writes the current basis to FileName1 in the data-file format,
    ! leaving Glob_DataFileName untouched so the run keeps
    ! checkpointing to the main data file afterwards.
    !----------------------------------------------------------------
    CASE ('SAVE_FILE')
      IF (Glob_BBOP(i)%A == Glob_CurrBasisSize) THEN
        IF (Glob_ProcID == 0) THEN
          iw = LEN_TRIM(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          WRITE(*, *)
          WRITE(*, *) 'Saving basis in file ', Glob_BBOP(i)%FileName1(1:iw), '...'
          CALL SaveResults(Filename=Glob_BBOP(i)%FileName1, Sort='no')
          WRITE(*, *) ' done'
          WRITE(*, *)
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! SAVE_HSWF  Method  BasisSize File1 File2 File3 File4
    !----------------------------------------------------------------
    ! Writes the (normalized) Hamiltonian, the overlap, the eigenvector
    ! and the whole wave function; 'none' in a slot skips that file.
    !----------------------------------------------------------------
    CASE ('SAVE_HSWF')
      IF (Glob_BBOP(i)%A == Glob_CurrBasisSize) THEN
        CALL SaveHSWF(Glob_BBOP(i)%FileName1, Glob_BBOP(i)%FileName2, &
                      Glob_BBOP(i)%FileName3, Glob_BBOP(i)%FileName4, &
                      Glob_BBOP(i)%GSEPSolutionMethod)
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0011 in main: incorrect BBOP step', i
          WRITE(*, *) 'Second parameter in SAVE_HSWF is incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! SAVE_HS_R  Method  BasisSize File1 File2
    !----------------------------------------------------------------
    ! Writes the UNNORMALIZED H and S. SAVE_HSWF writes the normalized
    ! pair, which is what the eigenproblem uses; the raw elements are
    ! not retained anywhere, so SaveHSRaw reconstructs them from the
    ! stored normalized matrices and Glob_diagS, and adds the basis
    ! health table.
    !----------------------------------------------------------------
    CASE ('SAVE_HS_R')
      IF (Glob_BBOP(i)%A == Glob_CurrBasisSize) THEN
        CALL SaveHSRaw(Glob_BBOP(i)%FileName1, Glob_BBOP(i)%FileName2, &
                       Glob_BBOP(i)%GSEPSolutionMethod)
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0012 in main: incorrect BBOP step', i
          WRITE(*, *) 'Second parameter in SAVE_HS_R is incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! OVERLAP_D  Eigensolver  BasisSize  [FileName]
    !----------------------------------------------------------------
    ! Diagonalizes the overlap matrix on its own and reports its
    ! spectrum and condition number. Same trigger rule as EXPC_VALS;
    ! the optional 4th field names the file that receives the full
    ! spectrum (ReadIOFile substitutes overlap.txt when it is absent).
    ! Only 'G' is available: the eigensolver field is kept
    ! so that the line looks like every other BBOP action and
    ! round-trips through SaveResults - see OverlapDiag for why the
    ! inverse-iteration path was dropped.
    !----------------------------------------------------------------
    CASE ('OVERLAP_D')
      IF (Glob_BBOP(i)%A == Glob_CurrBasisSize) THEN
        SELECT CASE (Glob_BBOP(i)%GSEPSolutionMethod)
        CASE ('G')
          CALL OverlapDiag('G', Glob_BBOP(i)%FileName1)
        CASE DEFAULT
          IF (Glob_ProcID == 0) THEN
            WRITE(*, '(1x,a,1x,i0,1x,a,1x,a)') 'Warning in main: OVERLAP_D step', i, &
              'has eigensolver', Glob_BBOP(i)%GSEPSolutionMethod
            WRITE(*, *) 'Only G is available for OVERLAP_D. Skipping this step...'
          ENDIF
        END SELECT
      ELSE
        IF (Glob_ProcID == 0) THEN
          WRITE(*, '(1x,a,1x,i0)') 'Error EC0013 in main: incorrect BBOP step', i
          WRITE(*, *) 'Second parameter in OVERLAP_D is incorrect'
        ENDIF
      ENDIF

    !----------------------------------------------------------------
    ! Anything else the parser let through
    !----------------------------------------------------------------
    CASE DEFAULT
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        WRITE(*, '(1x,a,1x,a,1x,a)') 'Warning in main: BBOP action', TRIM(Glob_BBOP(i)%Action), &
          'is recognized by the input parser but not implemented here.'
        WRITE(*, *) 'Skipping this step...'
      ENDIF

    END SELECT

  ENDDO


  !==================================================================
  ! Finish
  !==================================================================
  IF (Glob_ProcID == 0) THEN
    IF (Glob_UseSwapFile) THEN
      ! Empty swap file as it may contain the data saved after
      ! the last basis building and optimization program step
      OPEN(1, FILE=Glob_SwapFileName, FORM='unformatted', STATUS='replace', IOSTAT=OpenFileErr)
      IF (OpenFileErr == 0) THEN
        WRITE(1) 'Swap file is empty'
        CLOSE(1)
        WRITE(*, *)
        WRITE(*, *) 'Swap file has been cleaned up'
      ENDIF
    ENDIF
    WRITE(*, *) ' '
    WRITE(*, *) 'Basis Building and Optimization Program is completed'
    WRITE(*, *) 'Program has stopped'
  ENDIF

  CALL MPI_FINALIZE(Glob_MPIErrCode)

END PROGRAM main
