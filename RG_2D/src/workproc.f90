module workproc
!This module contains basic work subroutines
  use misc
  use matform
  use matelem
  use linalg
  use iso_fortran_env, only: int64
  use qrlinalg, only: qr_real_state,QR_SUCCESS,QR_ERR_INVALID_ARGUMENT, &
    QR_ERR_ALLOCATION,QR_ERR_INVALID_STATE, &
    QR_ERR_DIMENSION_MISMATCH,QR_ERR_NO_CONVERGENCE, &
    QR_ERR_NONPOSITIVE_OVERLAP
  implicit none

!The Q method keeps the physical basis in its canonical order. Unlike the G
!method, it must not move the functions currently being optimized to the end
!of the basis because the factors owned by qrlinalg describe one particular
!row and column order. Q_ActiveFunction maps an optimizer block to that fixed
!basis order; Q_ActivePosition is the inverse map.
!
!This storage is deliberately kept in workproc because the Q BBOP drivers and
!their trial/acceptance policy belong here. Physical H, S, nonlinear
!parameters, raw diagonal overlaps, derivatives, and linear coefficients
!remain in their existing Glob_ arrays. Factors is initialized and used only
!on rank zero because qrlinalg is serial; the physical matrices and the small
!relationship flags remain under workproc's MPI-aware control. Each numerical
!Q entry point constructs or updates that state through the transaction
!helpers below before requesting an eigenpair.
  integer,parameter :: Q_METHOD_SUCCESS=QR_SUCCESS
  integer,parameter :: Q_METHOD_INVALID_ARGUMENT=QR_ERR_INVALID_ARGUMENT
  integer,parameter :: Q_METHOD_ALLOCATION_ERROR=QR_ERR_ALLOCATION
  integer,parameter :: Q_METHOD_INVALID_STATE=QR_ERR_INVALID_STATE
  integer,parameter :: Q_METHOD_DIMENSION_MISMATCH=QR_ERR_DIMENSION_MISMATCH

  type :: QMethodWorkspace
    !The QR state is meaningful only on rank zero. Keeping it in the same
    !workspace makes ownership and lifetime explicit without exposing private
    !Q and R storage to the remainder of ECGPACK.
    type(qr_real_state) :: Factors

    !MatrixOrder is the active order represented by the physical matrices.
    !Capacity is the leading dimension allocated for Glob_H and Glob_S.
    !MaxActive is the largest simultaneous optimization block in this BBOP.
    integer :: MatrixOrder=0
    integer :: Capacity=0
    integer :: MaxActive=0
    integer :: NumActive=0

    !These flags describe relationships owned by workproc. They are not
    !duplicates of qrlinalg's internal validity flag: qrlinalg cannot know
    !whether a caller has changed Glob_H or Glob_S since the last update.
    logical :: MatricesAreCanonical=.false.
    logical :: FactorsMatchMatrices=.false.
    logical :: MatrixParametersAreStored=.false.
    logical :: TrialIsReady=.false.
    logical :: TrialHasDerivatives=.false.

    !Q_ActiveFunction(a) is the canonical basis index represented by optimizer
    !block a. Q_ActivePosition(i) is a, or zero when function i is inactive.
    integer,allocatable :: ActiveFunction(:)
    integer,allocatable :: ActivePosition(:)

    !A trial is assembled completely before any physical matrix column or QR
    !factor is changed. Columns here are full conceptual symmetric columns,
    !even though only the lower triangles of Glob_H and Glob_S are canonical.
    !Previous columns permit a failed multi-column transaction to be restored.
    real(wp),allocatable :: PreviousH(:,:),PreviousS(:,:)
    real(wp),allocatable :: TrialH(:,:),TrialS(:,:)
    real(wp),allocatable :: PreviousDiagS(:),TrialDiagS(:)

    !MatrixParam records the nonlinear parameters represented by the current
    !physical H/S columns. This must be independent of Glob_NonlinParam because
    !DRMNG writes its next requested point there before matrix assembly starts.
    !PreviousParam belongs to PreviousH/S and permits complete transaction
    !recovery. AcceptedParam is independent of both: a successful energy
    !evaluation does not by itself accept a trial as the optimizer's best.
    real(wp),allocatable :: MatrixParam(:,:),PreviousParam(:,:),AcceptedParam(:,:)
    real(wp) :: AcceptedEnergy=ZERO
    logical :: AcceptedPointIsStored=.false.

    !These diagnostics distinguish inverse-iteration accuracy from agreement
    !between the updated factors and the physical matrices. Both are relative
    !action residuals evaluated for the most recently returned eigenvector.
    real(wp) :: LastEigenpairResidual=huge(ONE)
    real(wp) :: LastFactorResidual=huge(ONE)
    integer :: FreshFactorizations=0

    !The qrlinalg solve requires distinct input and output vectors. DeltaH and
    !DeltaS hold one replacement column and are also useful for residual and
    !fresh-factorization checks. They are replicated because the surrounding
    !transaction is collective, but only rank zero passes them to qrlinalg.
    real(wp),allocatable :: InitialVector(:),SolvedVector(:)
    real(wp),allocatable :: DeltaH(:),DeltaS(:)
  end type QMethodWorkspace

  type(QMethodWorkspace),save :: Q_Workspace

contains

  subroutine ReadIOFile()
!Subroutine ReadIOFile reads data (nonlinear variational
!parameters and other information) from the input/output
!file whose name is specified by global variable
!Glob_DataFileName. If there is no such a file in the
!current directory then the program stops.  On entry to DSYGVX

!Local variables:
    integer        OpenFileErr
    real(wp)    ReadRealA,ReadRealB
    real(wp),allocatable,dimension(:) :: ReadRealArr
    integer        ReadInt,ReadErr
    integer        WorkInt(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
    real(wp),allocatable,dimension(:) :: WorkBuffReal
    integer,allocatable,dimension(:)     :: WorkBuffInt
    integer        i,j,Line,j1,j2,j3,j4
    character(70)  ReadChar
    character(5)   ReadBasisType  !Basis type specifier from the optional BASIS_TYPE line
    character(256) ReadLine
    integer        ReadLineErr
    logical        ErrorInDataFile,IsBBOPStep

    ErrorInDataFile=.false.

    if (Glob_ProcID==0) then
      open(1,file=Glob_DataFileName,status='old',iostat=OpenFileErr)
      if (OpenFileErr/=0) then
        write (*,*) 'Error EC0100 in DataFileInit: data file not found - ',trim(adjustl(Glob_DataFileName))
        ErrorInDataFile=.true.
      endif
    endif

    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

!Reading information
    if (Glob_ProcID==0) Line=0
    if (Glob_ProcID==0) then
      write(*,*) 'Reading initial conditions from data file ',trim(adjustl(Glob_DataFileName))
!Read the optional BASIS_TYPE line that may precede the PARTICLES line. The whole
!record is read into a buffer first and parsed internally so that backspace 1
!reliably pushes the line back if it turns out not to be a BASIS_TYPE line.
      read(1,'(A)',iostat=ReadLineErr) ReadLine
      if (ReadLineErr==0) read(ReadLine,*,iostat=ReadErr) ReadChar(1:10),ReadBasisType
      if ((ReadLineErr==0).and.(ReadErr==0).and.(ReadChar(1:10)=='BASIS_TYPE')) then
        if (ReadBasisType/=Glob_BasisType) then
          write(*,*) 'Error EC0101 in data file: basis type specifier "',ReadBasisType, &
            '" does not match the basis type of this code "',Glob_BasisType,'"'
          ErrorInDataFile=.true.
        else
          write(*,'(1x,a10,1x,a5)') ReadChar(1:10),ReadBasisType
          Glob_BasisTypeSupplied=.true.
          Line=Line+1
        endif
      else
        if (ReadLineErr==0) backspace 1
      endif
      read(1,*) ReadChar(1:9),ReadInt
      write(*,'(1x,a9,1x,i6)') ReadChar(1:9),ReadInt
      Line=Line+1
      Glob_n=ReadInt-1 !Glob_n is the number of pseudoparticles
      if ((Glob_n<1).or.(ReadChar(1:9)/='PARTICLES')) then
        write(*,*) 'Error EC0102 in data file, line ',Line
        ErrorInDataFile=.true.
      endif
    endif
    call MPI_BCAST(Glob_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_BasisTypeSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (Glob_n/=Glob_AllowedNumOfPseudoParticles) then
      if (Glob_ProcID==0) then
        write (*,*) 'The version of the code you are running was compiled for the case'
        write (*,*) 'when the number of particles in the system is equal to', &
          Glob_AllowedNumOfParticles
        write (*,*) 'while the number of particles specified in the input file is',Glob_n+1
        write (*,*) 'Please make appropriate changes. Program will now stop.'
      endif
      ErrorInDataFile=.true.
    endif
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    Glob_np=Glob_n*(Glob_n+1)/2
    Glob_npt=Glob_np
    Glob_2Raised3n2=TWO**((3*Glob_n)/TWO)
    Glob_PiRaised3n2=Glob_Pi**((3*Glob_n)/TWO)

!Read the optional VECTOR_COUPLING_SCHEME line that may appear right after the
!PARTICLES line. As with the other optional descriptors, the record is read into
!a buffer first and parsed internally so that backspace 1 reliably pushes the
!line back if it turns out not to be a VECTOR_COUPLING_SCHEME line.
    if (Glob_ProcID==0) then
      read(1,'(A)',iostat=ReadLineErr) ReadLine
      if (ReadLineErr==0) then
        read(ReadLine,*,iostat=ReadErr) ReadChar(1:22),ReadInt
        if ((ReadErr==0).and.(ReadChar(1:22)=='VECTOR_COUPLING_SCHEME')) then
          Glob_IsVectorCouplingSchemeProvided=.true.
          Glob_VectorCouplingScheme=ReadInt
          write(*,'(1x,a22,1x,i6)') ReadChar(1:22),Glob_VectorCouplingScheme
          Line=Line+1
          if ((Glob_VectorCouplingScheme<0).or.(Glob_VectorCouplingScheme>2)) then
            write(*,*) 'Error EC0105 in data file, line ',Line
            write(*,*) 'VECTOR_COUPLING_SCHEME value must be 0, 1, or 2'
            ErrorInDataFile=.true.
          endif
          !Scheme 1 places two pseudoparticles in p-states, which requires the two
          !indices to be different. This is impossible when there is only one pseudoparticle.
          if ((Glob_VectorCouplingScheme==1).and.(Glob_n<2)) then
            write(*,*) 'Error EC0107 in data file, line ',Line
            write(*,*) 'VECTOR_COUPLING_SCHEME 1 requires at least two pseudoparticles'
            ErrorInDataFile=.true.
          endif
        else
          backspace 1
        endif
      endif
    endif
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    call MPI_BCAST(Glob_IsVectorCouplingSchemeProvided,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_VectorCouplingScheme,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!Read the optional FIXED_INDEX_1 and FIXED_INDEX_2 lines that may appear between
!the PARTICLES and MASSES lines (in any order). As with the other optional
!descriptors, each record is read into a buffer first and parsed internally so
!that backspace 1 reliably pushes the line back if it turns out not to be a
!FIXED_INDEX_1/FIXED_INDEX_2 line.
    if (Glob_ProcID==0) then
      do i=1,2
        read(1,'(A)',iostat=ReadLineErr) ReadLine
        if (ReadLineErr/=0) exit
        read(ReadLine,*,iostat=ReadErr) ReadChar(1:13),ReadInt
        if ((ReadErr==0).and.(ReadChar(1:13)=='FIXED_INDEX_1')) then
          Glob_IsIndexFixed(1)=.true.
          Glob_IndexFixedValue(1)=ReadInt
          write(*,'(1x,a13,1x,i6)') ReadChar(1:13),Glob_IndexFixedValue(1)
          Line=Line+1
          if ((Glob_IndexFixedValue(1)<1).or.(Glob_IndexFixedValue(1)>Glob_n)) then
            write(*,*) 'Error EC0103 in data file, line ',Line
            write(*,*) 'FIXED_INDEX_1 value must be in the range from 1 to',Glob_n
            ErrorInDataFile=.true.
          endif
        elseif ((ReadErr==0).and.(ReadChar(1:13)=='FIXED_INDEX_2')) then
          Glob_IsIndexFixed(2)=.true.
          Glob_IndexFixedValue(2)=ReadInt
          write(*,'(1x,a13,1x,i6)') ReadChar(1:13),Glob_IndexFixedValue(2)
          Line=Line+1
          if ((Glob_IndexFixedValue(2)<1).or.(Glob_IndexFixedValue(2)>Glob_n)) then
            write(*,*) 'Error EC0104 in data file, line ',Line
            write(*,*) 'FIXED_INDEX_2 value must be in the range from 1 to',Glob_n
            ErrorInDataFile=.true.
          endif
        else
          backspace 1
          exit
        endif
      enddo
!When the vector coupling scheme is explicitly specified, only one of the two
!indices of the basis functions may be fixed.
      if (Glob_IsVectorCouplingSchemeProvided.and. &
          Glob_IsIndexFixed(1).and.Glob_IsIndexFixed(2)) then
        write(*,*) 'Error EC0106 in data file, line ',Line
        write(*,*) 'FIXED_INDEX_1 and FIXED_INDEX_2 cannot both be present'
        write(*,*) 'when VECTOR_COUPLING_SCHEME is specified'
        ErrorInDataFile=.true.
      endif
    endif
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    call MPI_BCAST(Glob_IsIndexFixed,2,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_IndexFixedValue,2,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    allocate(Glob_Mass(Glob_n+1))
    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:6),Glob_Mass(1:Glob_n+1)
      write(*,'(1x,a6)',advance='no') ReadChar(1:6)
      call writerealarradv(6,Glob_Mass,Glob_n+1)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_Mass,Glob_n+1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    allocate(Glob_PseudoCharge(Glob_n))
    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:7),Glob_PseudoCharge0,Glob_PseudoCharge(1:Glob_n)
      write(*,'(1x,a7)',advance='no') ReadChar(1:7)
      call writereal(6,Glob_PseudoCharge0)
      call writerealarradv(6,Glob_PseudoCharge,Glob_n)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_PseudoCharge0,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_PseudoCharge,Glob_n,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      Glob_RepulsionScalingParam=1.0_wp
      Glob_RepScalParamSupplied=.false.
      Glob_RepulsionScalingParamPlus=1.0_wp
      Glob_RepScalParamPlusSupplied=.false.
      Glob_RepulsionScalingParamMinus=1.0_wp
      Glob_RepScalParamMinusSupplied=.false.
      do i=1,3
        !Read exactly one record into a character buffer first, then parse it
        !with an internal read. This keeps unit 1 advancing by a single record
        !so that backspace 1 reliably pushes the line back on all compilers
        !(list-directed reads from unit 1 may span several records, which makes
        !backspace behave differently in gfortran vs ifort).
        read(1,'(A)',iostat=ReadLineErr) ReadLine
        if (ReadLineErr/=0) exit
        read(ReadLine,*,iostat=ReadErr) ReadChar(1:29),ReadRealA
        if ((ReadErr/=0).or.(ReadChar(1:23)/='REPULSION_SCALING_PARAM')) then
          backspace 1
        else
          if (ReadChar(1:28)=='REPULSION_SCALING_PARAM_PLUS') then
            Glob_RepulsionScalingParamPlus=ReadRealA
            Glob_RepScalParamPlusSupplied=.true.
            write(*,'(1x,a28)',advance='no') ReadChar(1:28)
            call writerealadv(6,Glob_RepulsionScalingParamPlus)
          elseif (ReadChar(1:29)=='REPULSION_SCALING_PARAM_MINUS') then
            Glob_RepulsionScalingParamMinus=ReadRealA
            Glob_RepScalParamMinusSupplied=.true.
            write(*,'(1x,a28)',advance='no') ReadChar(1:29)
            call writerealadv(6,Glob_RepulsionScalingParamMinus)
          else
            Glob_RepulsionScalingParam=ReadRealA
            Glob_RepScalParamSupplied=.true.
            write(*,'(1x,a23)',advance='no') ReadChar(1:23)
            call writerealadv(6,Glob_RepulsionScalingParam)
          endif
          Line=Line+1
        endif
      enddo
    endif
    call MPI_BCAST(Glob_RepScalParamSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RepulsionScalingParam,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RepScalParamPlusSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RepulsionScalingParamPlus,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RepScalParamMinusSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RepulsionScalingParamMinus,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      !Read one record into a buffer and parse it internally (see the comment in
      !the repulsion-scaling block above) so that backspace 1 is portable.
      read(1,'(A)',iostat=ReadLineErr) ReadLine
      if (ReadLineErr==0) read(ReadLine,*,iostat=ReadErr) ReadChar(1:24),Glob_AttractionScalingParam
      if ((ReadLineErr/=0).or.(ReadErr/=0).or.(ReadChar(1:24)/='ATTRACTION_SCALING_PARAM')) then
        Glob_AttractionScalingParam=1.0_wp
        Glob_AttrScalParamSupplied=.false.
        if (ReadLineErr==0) backspace 1
      else
        Glob_AttrScalParamSupplied=.true.
        write(*,'(1x,a24)',advance='no') ReadChar(1:24)
        call writerealadv(6,Glob_AttractionScalingParam)
        Line=Line+1
      endif
    endif
    call MPI_BCAST(Glob_AttrScalParamSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_AttractionScalingParam,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:8),Glob_YOperatorString
      j=len_trim(Glob_YOperatorString)
      write(*,'(1x,a8)',advance='no') ReadChar(1:8)
      call writestringadv(6,Glob_YOperatorString,j)
      Line=Line+1
    endif
    do i=1,Glob_YOperatorStringLength
      WorkInt(i)=ichar(Glob_YOperatorString(i:i))
    enddo
    call MPI_BCAST(WorkInt,Glob_YOperatorStringLength,MPI_INTEGER,0, &
                   MPI_COMM_WORLD,Glob_MPIErrCode)
    do i=1,Glob_YOperatorStringLength
      Glob_YOperatorString(i:i)=char(WorkInt(i))
    enddo

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:10),Glob_CurrBasisSize
      write(*,'(1x,a10,1x,i6)')  ReadChar(1:10),Glob_CurrBasisSize
      Line=Line+1
    endif
    call MPI_BCAST(Glob_CurrBasisSize,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:14),Glob_CurrEnergy
      write(*,'(1x,a14)',advance='no')  ReadChar(1:14)
      call writerealadv(6,Glob_CurrEnergy)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_CurrEnergy,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:16),Glob_WhichEigenvalue
      write(*,'(1x,a16,1x,i6)')  ReadChar(1:16),Glob_WhichEigenvalue
      Line=Line+1
    endif
    call MPI_BCAST(Glob_WhichEigenvalue,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:16),Glob_EigvalTol
      write(*,'(1x,a16)',advance='no')  ReadChar(1:16)
      call writerealadv(6,Glob_EigvalTol)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_EigvalTol,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:14),Glob_InvItParameter
      write(*,'(1x,a14)',advance='no')  ReadChar(1:14)
      call writerealadv(6,Glob_InvItParameter)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_InvItParameter,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:15),Glob_LastEigvalTol
      write(*,'(1x,a15)',advance='no')  ReadChar(1:15)
      call writerealadv(6,Glob_LastEigvalTol)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_LastEigvalTol,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:15),Glob_BestEigvalTol
      write(*,'(1x,a15)',advance='no')  ReadChar(1:15)
      call writerealadv(6,Glob_BestEigvalTol)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_BestEigvalTol,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:16),Glob_WorstEigvalTol
      write(*,'(1x,a16)',advance='no')  ReadChar(1:16)
      call writerealadv(6,Glob_WorstEigvalTol)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_WorstEigvalTol,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,*) ReadChar(1:15),Glob_RG_p1,Glob_RG_s1,Glob_RG_s2
      write(*,'(1x,a15)',advance='no') ReadChar(1:15)
      call writereal(6,Glob_RG_p1)
      call writereal(6,Glob_RG_s1)
      call writerealadv(6,Glob_RG_s2)
      Line=Line+1
    endif
    call MPI_BCAST(Glob_RG_p1,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RG_s1,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_RG_s2,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Glob_ProcID==0) then
      read(1,'(a70)')   ReadChar(1:70)
      write(*,'(a70)')  ReadChar(1:70)
      read(1,'(a70)')   ReadChar(1:70)
      write(*,'(a70)')  ReadChar(1:70)
      read(1,'(a70)')   ReadChar(1:70)
      write(*,'(a70)')  ReadChar(1:70)
      Line=Line+3
    endif

!Reading Basis Building and Optimization Program
    ReadChar(1:70)=' '
    if (Glob_ProcID==0) then
      Glob_NumOfBBOPSteps=0
      IsBBOPStep=.true.
      do while (IsBBOPStep)
        read(1,*) ReadChar(1:9)
        if ((ReadChar(1:9)=='BASIS_ENL').or.(ReadChar(1:9)=='OPT_CYCLE').or.  &
            (ReadChar(1:9)=='FULL_OPT1').or.(ReadChar(1:9)=='EXPC_VALS').or.  &
            (ReadChar(1:9)=='ELIM_LCFN').or.(ReadChar(1:9)=='ELIM_LND1').or.  &
            (ReadChar(1:9)=='SEPR_LND1').or.(ReadChar(1:9)=='SEPR_FLCF').or.  &
            (ReadChar(1:9)=='DENSITIES').or.(ReadChar(1:9)=='SAVE_FILE').or.  &
            (ReadChar(1:9)=='SAVE_HSWF')) then
          Glob_NumOfBBOPSteps=Glob_NumOfBBOPSteps+1
        else
          IsBBOPStep=.false.
        endif
      enddo
      do i=1,Glob_NumOfBBOPSteps+1
        backspace 1
      enddo
    endif
    call MPI_BCAST(Glob_NumOfBBOPSteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    allocate(Glob_BBOP(Glob_NumOfBBOPSteps))

    if (Glob_ProcID==0) then
      do i=1,Glob_NumOfBBOPSteps
        read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod
      enddo
      do i=1,Glob_NumOfBBOPSteps
        backspace 1
      enddo
      do i=1,Glob_NumOfBBOPSteps
        select case (Glob_BBOP(i)%Action(1:9))
        case('BASIS_ENL')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod , &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C, &
            Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R
          !write(*,'(1x,a9,1x,a1,5(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
          !                   Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writerealadv(6,Glob_BBOP(i)%R)
        case('OPT_CYCLE')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C, &
            Glob_BBOP(i)%D,Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G, &
            Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%H
          !write(*,'(1x,a9,1x,a1,7(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
          !                 Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E, &
          !                 Glob_BBOP(i)%F,Glob_BBOP(i)%G
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writereal(6,Glob_BBOP(i)%R)
          !write(*,'(1x,i6)') Glob_BBOP(i)%H
          Glob_IsOptCycleScripted=.true.
        case('FULL_OPT1')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D, &
            Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
            Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,4(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
          !                 Glob_BBOP(i)%C,Glob_BBOP(i)%D
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writereal(6,Glob_BBOP(i)%R)
          !write(*,'(2(1x,i6),1x)',advance='no') Glob_BBOP(i)%E,Glob_BBOP(i)%F
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('EXPC_VALS')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !write(*,'(1x,a9,1x,a1,1x,i6)') Glob_BBOP(i)%Action(1:9),  &
          !        Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        case('DENSITIES')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName3(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)
          j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,1x,i6)',advance='no')              &
          !      Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !call writestring(6,Glob_BBOP(i)%FileName1,j1)
          !call writestring(6,Glob_BBOP(i)%FileName2,j2)
          !call writestring(6,Glob_BBOP(i)%FileName3,j3)
          !call writestringadv(6,Glob_BBOP(i)%FileName4,j4)
        case('ELIM_LCFN')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('ELIM_LND1')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('SEPR_LND1')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writereal(6,Glob_BBOP(i)%R)
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('SEPR_FLCF')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
          !      Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          !call writereal(6,Glob_BBOP(i)%Q)
          !call writereal(6,Glob_BBOP(i)%R)
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('SAVE_FILE')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A, &
            Glob_BBOP(i)%FileName1(1:Glob_FileNameLength)
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          !write(*,'(1x,a9,1x,i6,1x)',advance='no') Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A
          !call writestringadv(6,Glob_BBOP(i)%FileName1,j)
        case('SAVE_HSWF')
          read(1,*) Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod, &
            Glob_BBOP(i)%A,Glob_BBOP(i)%FileName1(1:Glob_FileNameLength), &
            Glob_BBOP(i)%FileName2(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName3(1:Glob_FileNameLength),                &
            Glob_BBOP(i)%FileName4(1:Glob_FileNameLength)
          j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
        endselect
      enddo
      read(1,'(a70)')   ReadChar(1:70)
      !write(*,'(a70)')  ReadChar(1:70)
    endif
    do i=1,Glob_NumOfBBOPSteps
      do j=1,9
        WorkInt(j)=ichar(Glob_BBOP(i)%Action(j:j))
      enddo
      call MPI_BCAST(WorkInt,9,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,9
        Glob_BBOP(i)%Action(j:j)=char(WorkInt(j))
      enddo
      j=ichar(Glob_BBOP(i)%GSEPSolutionMethod)
      call MPI_BCAST(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      Glob_BBOP(i)%GSEPSolutionMethod=char(j)
      call MPI_BCAST(Glob_BBOP(i)%A,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%B,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%C,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%D,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%E,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%F,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%G,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%H,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%Q,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_BBOP(i)%R,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,Glob_FileNameLength
        WorkInt(j)=ichar(Glob_BBOP(i)%FileName1(j:j))
      enddo
      call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,Glob_FileNameLength
        Glob_BBOP(i)%FileName1(j:j)=char(WorkInt(j))
      enddo
      do j=1,Glob_FileNameLength
        WorkInt(j)=ichar(Glob_BBOP(i)%FileName2(j:j))
      enddo
      call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,Glob_FileNameLength
        Glob_BBOP(i)%FileName2(j:j)=char(WorkInt(j))
      enddo
      do j=1,Glob_FileNameLength
        WorkInt(j)=ichar(Glob_BBOP(i)%FileName3(j:j))
      enddo
      call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,Glob_FileNameLength
        Glob_BBOP(i)%FileName3(j:j)=char(WorkInt(j))
      enddo
      do j=1,Glob_FileNameLength
        WorkInt(j)=ichar(Glob_BBOP(i)%FileName4(j:j))
      enddo
      call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do j=1,Glob_FileNameLength
        Glob_BBOP(i)%FileName4(j:j)=char(WorkInt(j))
      enddo
    enddo
    call MPI_BCAST(Glob_IsOptCycleScripted,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    allocate(Glob_History(Glob_CurrBasisSize))
    allocate(Glob_FuncNum(Glob_CurrBasisSize))
    allocate(Glob_Index(Glob_CurrBasisSize,2))
    allocate(Glob_NonlinParam(Glob_npt,Glob_CurrBasisSize))

    if (Glob_CurrBasisSize==0) then
      !Stop reading if basis size is zero
      close(1)
      return
    endif

    allocate(WorkBuffReal(Glob_CurrBasisSize))
    allocate(WorkBuffInt(Glob_CurrBasisSize))

!Reading history
    if (Glob_ProcID==0) then
      do i=1,Glob_CurrBasisSize
        read(1,*)  j,Glob_History(i)%Energy,Glob_History(i)%CyclesDone, &
          Glob_History(i)%InitFuncAtLastStep, &
          Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      enddo
      read(1,*)   ReadChar(1:70)
      write(*,*)
    endif
    do i=1,Glob_CurrBasisSize
      WorkBuffReal(i)=Glob_History(i)%Energy
    enddo
    call MPI_BCAST(WorkBuffReal,Glob_CurrBasisSize,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    do i=1,Glob_CurrBasisSize
      Glob_History(i)%Energy=WorkBuffReal(i)
    enddo
    do i=1,Glob_CurrBasisSize
      WorkBuffInt(i)=Glob_History(i)%CyclesDone
    enddo
    call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    do i=1,Glob_CurrBasisSize
      Glob_History(i)%CyclesDone=WorkBuffInt(i)
    enddo
    do i=1,Glob_CurrBasisSize
      WorkBuffInt(i)=Glob_History(i)%InitFuncAtLastStep
    enddo
    call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    do i=1,Glob_CurrBasisSize
      Glob_History(i)%InitFuncAtLastStep=WorkBuffInt(i)
    enddo
    do i=1,Glob_CurrBasisSize
      WorkBuffInt(i)=Glob_History(i)%NumOfEnergyEvalDuringFullOpt
    enddo
    call MPI_BCAST(WorkBuffInt,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    do i=1,Glob_CurrBasisSize
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=WorkBuffInt(i)
    enddo

    deallocate(WorkBuffReal)
    deallocate(WorkBuffInt)

!Reading basis functions
    if (Glob_ProcID==0) then
      do i=1,Glob_CurrBasisSize
        read(1,*) j,Glob_Index(i,1),Glob_Index(i,2),Glob_NonlinParam(1:Glob_npt,i)
      enddo
    endif
    call MPI_BCAST(Glob_NonlinParam,Glob_npt*Glob_CurrBasisSize, &
                   MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_Index,2*Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!Setting function numbers as they were read
    do i=1,Glob_CurrBasisSize
      Glob_FuncNum(i)=i
    enddo

    close(1)

  end subroutine ReadIOFile

  subroutine SaveResults(FileName,Sort)
!Subroutine SaveResults saves current results of the
!calculations.
!  Input parameters :
!  Filename - Optional file name, where the results need to be saved. If it is
!         not passed then the result are saved in the file whose name is defined
!         by global variable Glob_DataFileName
!  Sort - Optional string variable. If its value is equal to 'yes'
!         or 'YES' then the basis functions are sorted out according to
!         their numbers. In this case a work array is needed. The
!         subroutine uses array Glob_IntWorkArrForSaveResults. Thus,
!         in order to do calls with Sort='yes' or without passing one first
!         needs to allocate that global array. The length of this integer
!         array should be at least Glob_CurrBasisSize. If the value of
!         Sort is not equal to 'yes' or 'YES', or it is simply not passed, then
!         the basis functions will not be sorted out. In this case the allocation
!         of global array Glob_IntWorkArrForSaveResults is not necessary.
!Arguments:
    character(*)   FileName,Sort
    optional   ::  FileName,Sort
!Local variables:
    integer i,j,j1,j2,j3,j4
    logical SortNeeded

    if (Glob_ProcID==0) then
      if (present(FileName)) then
        open(1,file=FileName,status='replace')
      else
        open(1,file=Glob_DataFileName,status='replace')
      endif
      if (Glob_BasisTypeSupplied) write(1,'(1x,a10,1x,a5)') 'BASIS_TYPE',Glob_BasisType
      write(1,'(1x,a9,1x,i6)') 'PARTICLES',Glob_n+1
      if (Glob_IsVectorCouplingSchemeProvided) &
        write(1,'(1x,a22,1x,i6)') 'VECTOR_COUPLING_SCHEME',Glob_VectorCouplingScheme
      if (Glob_IsIndexFixed(1)) write(1,'(1x,a13,1x,i6)') 'FIXED_INDEX_1',Glob_IndexFixedValue(1)
      if (Glob_IsIndexFixed(2)) write(1,'(1x,a13,1x,i6)') 'FIXED_INDEX_2',Glob_IndexFixedValue(2)
      write(1,'(1x,a6)',advance='no') 'MASSES'
      call writerealarradv(1,Glob_Mass,Glob_n+1)
      write(1,'(1x,a7)',advance='no') 'CHARGES'
      call writereal(1,Glob_PseudoCharge0)
      call writerealarradv(1,Glob_PseudoCharge,Glob_n)
      if (Glob_RepScalParamSupplied) then
        write(1,'(1x,a23)',advance='no') 'REPULSION_SCALING_PARAM'
        call writerealadv(1,Glob_RepulsionScalingParam)
      endif
      if (Glob_RepScalParamPlusSupplied) then
        write(1,'(1x,a28)',advance='no') 'REPULSION_SCALING_PARAM_PLUS'
        call writerealadv(1,Glob_RepulsionScalingParamPlus)
      endif
      if (Glob_RepScalParamMinusSupplied) then
        write(1,'(1x,a29)',advance='no') 'REPULSION_SCALING_PARAM_MINUS'
        call writerealadv(1,Glob_RepulsionScalingParamMinus)
      endif
      if (Glob_AttrScalParamSupplied) then
        write(1,'(1x,a24)',advance='no') 'ATTRACTION_SCALING_PARAM'
        call writerealadv(1,Glob_AttractionScalingParam)
      endif
      i=len_trim(Glob_YOperatorString)
      write(1,'(1x,a8)',advance='no') 'SYMMETRY'
      call writestringadv(1,Glob_YOperatorString,i)
      write(1,'(1x,a10,1x,i6)') 'BASIS_SIZE',Glob_CurrBasisSize
      write(1,'(1x,a14)',advance='no') 'CURRENT_ENERGY'
      call writerealadv(1,Glob_CurrEnergy)
      write(1,'(1x,a16,1x,i6)') 'WHICH_EIGENVALUE',Glob_WhichEigenvalue
      write(1,'(1x,a16)',advance='no') 'EIGVAL_TOLERANCE'
      call writerealadv(1,Glob_EigvalTol)
      write(1,'(1x,a14)',advance='no') 'INVITPARAMETER'
      call writerealadv(1,Glob_InvItParameter)
      write(1,'(1x,a15)',advance='no') 'LAST_EIGVAL_TOL'
      call writerealadv(1,Glob_LastEigvalTol)
      write(1,'(1x,a15)',advance='no') 'BEST_EIGVAL_TOL'
      call writerealadv(1,Glob_BestEigvalTol)
      write(1,'(1x,a16)',advance='no') 'WORST_EIGVAL_TOL'
      call writerealadv(1,Glob_WorstEigvalTol)
      write(1,'(1x,a15)',advance='no') 'GENERATOR_PARAM'
      call writereal(1,Glob_RG_p1)
      call writereal(1,Glob_RG_s1)
      call writerealadv(1,Glob_RG_s2)
      write(1,*) '=============================='
      i=Glob_CurrBasisSize
      write(1,'(1x,i6)',advance='no')  i
      call writereal(1,Glob_History(i)%Energy)
      write(1,'(3(1x,i6))') Glob_History(i)%CyclesDone,  &
        Glob_History(i)%InitFuncAtLastStep,Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      write(1,*) '=============================='
      do i=1,Glob_NumOfBBOPSteps
        select case (Glob_BBOP(i)%Action(1:9))
        case('BASIS_ENL')
          write(1,'(1x,a9,1x,a1,5(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E
          call writereal(1,Glob_BBOP(i)%Q)
          call writerealadv(1,Glob_BBOP(i)%R)
        case('OPT_CYCLE')
          write(1,'(1x,a9,1x,a1,7(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%E, &
            Glob_BBOP(i)%F,Glob_BBOP(i)%G
          call writereal(1,Glob_BBOP(i)%Q)
          call writereal(1,Glob_BBOP(i)%R)
          write(1,'(1x,i6)') Glob_BBOP(i)%H
        case('FULL_OPT1')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,4(1x,i6))',advance='no') Glob_BBOP(i)%Action(1:9), &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A,Glob_BBOP(i)%B, &
            Glob_BBOP(i)%C,Glob_BBOP(i)%D
          call writereal(1,Glob_BBOP(i)%Q)
          call writereal(1,Glob_BBOP(i)%R)
          write(1,'(2(1x,i6),1x)',advance='no') Glob_BBOP(i)%E,Glob_BBOP(i)%F
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('EXPC_VALS')
          write(1,'(1x,a9,1x,a1,1x,i6)') Glob_BBOP(i)%Action(1:9),  &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
        case('DENSITIES')
          j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,1x,i6)',advance='no')              &
            Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writestring(1,Glob_BBOP(i)%FileName1,j1)
          call writestring(1,Glob_BBOP(i)%FileName2,j2)
          call writestring(1,Glob_BBOP(i)%FileName3,j3)
          call writestringadv(1,Glob_BBOP(i)%FileName4,j4)
        case('ELIM_LCFN')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writereal(1,Glob_BBOP(i)%Q)
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('ELIM_LND1')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writereal(1,Glob_BBOP(i)%Q)
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('SEPR_LND1')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writereal(1,Glob_BBOP(i)%Q)
          call writereal(1,Glob_BBOP(i)%R)
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('SEPR_FLCF')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,1x,i6)',advance='no') Glob_BBOP(i)%Action(1:9),  &
            Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writereal(1,Glob_BBOP(i)%Q)
          call writereal(1,Glob_BBOP(i)%R)
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('SAVE_FILE')
          j=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,i6,1x)',advance='no') Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%A
          call writestringadv(1,Glob_BBOP(i)%FileName1,j)
        case('SAVE_HSWF')
          j1=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          j2=len_trim(Glob_BBOP(i)%FileName2(1:Glob_FileNameLength))
          j3=len_trim(Glob_BBOP(i)%FileName3(1:Glob_FileNameLength))
          j4=len_trim(Glob_BBOP(i)%FileName4(1:Glob_FileNameLength))
          write(1,'(1x,a9,1x,a1,1x,i6)',advance='no')              &
            Glob_BBOP(i)%Action(1:9),Glob_BBOP(i)%GSEPSolutionMethod,Glob_BBOP(i)%A
          call writestring(1,Glob_BBOP(i)%FileName1,j1)
          call writestring(1,Glob_BBOP(i)%FileName2,j2)
          call writestring(1,Glob_BBOP(i)%FileName3,j3)
          call writestringadv(1,Glob_BBOP(i)%FileName4,j4)
        endselect
      enddo
      write(1,*) '=============================='
      do i=1,Glob_CurrBasisSize
        write(1,'(1x,i6)',advance='no')  i
        call writereal(1,Glob_History(i)%Energy)
        write(1,'(3(1x,i6))') Glob_History(i)%CyclesDone,  &
          Glob_History(i)%InitFuncAtLastStep,Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      enddo
      write(1,*) '=============================='
      SortNeeded=.false.
      if (present(Sort)) then
        if ((Sort=='yes').or.(Sort=='YES')) SortNeeded=.true.
      endif
      if (SortNeeded) then
        do i=1,Glob_CurrBasisSize
          Glob_IntWorkArrForSaveResults(Glob_FuncNum(i))=i
        enddo
        do i=1,Glob_CurrBasisSize
          write(1,'(1x,i6,1x,i6,1x,i6)',advance='no') i,Glob_Index(Glob_IntWorkArrForSaveResults(i),1)&
            ,Glob_Index(Glob_IntWorkArrForSaveResults(i),2)
          call writerealarradv(1,Glob_NonlinParam(1:Glob_npt,Glob_IntWorkArrForSaveResults(i)),Glob_npt)
        enddo
      else
        do i=1,Glob_CurrBasisSize
          write(1,'(1x,i6,1x,i6,1x,i6)',advance='no') i,Glob_Index(i,1),Glob_Index(i,2)
          call writerealarradv(1,Glob_NonlinParam(1:Glob_npt,i),Glob_npt)
        enddo
      endif
      close(1)
    endif

  end subroutine SaveResults

  subroutine ReadBlackList()
!Subroutine ReadBlackList reads the file containing the list of
!basis functions that must not be optimized in cyclic optimization
!routines (the name of the file is specified by Glob_BlackListFileName)
!The file should contain function numbers (not necessarily ordered), one per
!line, no blank lines, repetitions are ok.
!If the file does not exist or does not contain anything then it is assumed
!that no functions are blacklisted.
!Local variables:
    integer        OpenFileErr,i,ReadErr
    logical        ErrorInFile

    ErrorInFile=.false.
    if (Glob_ProcID==0) then
      open(1,file=Glob_BlackListFileName,status='old',iostat=OpenFileErr)
      if (OpenFileErr/=0) then
        write(*,*)
        write(*,*) 'Black list file not found: ',trim(Glob_BlackListFileName)
        write(*,*) 'No basis functions will be blacklisted during cyclic optimization'
        write(*,*)
        ErrorInFile=.true.
      endif
    endif

    call MPI_BCAST(ErrorInFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInFile) then
      Glob_lbf=0
      return
    endif

    if (Glob_ProcID==0) then
      Glob_lbf=0
      ReadErr=0
      do while ((ReadErr==0).and.(.not.ErrorInFile))
        read(1,*,iostat=ReadErr) i
        if (ReadErr==0) then
          if ((i<=Glob_CurrBasisSize).and.(i>0)) then
            Glob_lbf=max(i,Glob_lbf)
          else
            ErrorInFile=.true.
          endif
        endif
      enddo
    endif
    call MPI_BCAST(ErrorInFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (ErrorInFile) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0110 in ReadBlackList: incorrect values in file ',trim(Glob_BlackListFileName)
        close(1)
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    call MPI_BCAST(Glob_lbf,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (Glob_lbf==0) then
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,*) 'Black list file is empty - ',trim(Glob_BlackListFileName)
        write(*,*) 'no basis functions will be blacklisted during cyclic optimization'
        write(*,*)
      endif
      return
    endif
    allocate(Glob_Blacklisted(Glob_lbf))
    if (Glob_ProcID==0) then
      Glob_Blacklisted(1:Glob_lbf)=.false.
      rewind(1)
      ReadErr=0
      do while (ReadErr==0)
        read(1,*,iostat=ReadErr) i
        if (ReadErr==0) Glob_Blacklisted(i)=.true.
      enddo
      close(1)
      write(*,*)
      write(*,*) 'Black list file has been read - ',trim(Glob_BlackListFileName)
      write(*,*) 'The following basis functions will be blacklisted during cyclic optimization:'
      do i=1,Glob_lbf
        if (Glob_Blacklisted(i)) write(*,'(1x,i6)',advance='no') i
      enddo
      write(*,*)
      write(*,*)
    endif
    call MPI_BCAST(Glob_Blacklisted,Glob_lbf,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)

  end subroutine ReadBlackList

  subroutine ProgramDataInit()
!Subroutine ProgramDataInit initializes some data needed for
!calculations. It should be called at the start of the program,
!right after reading input/output file.

    integer       n,npart,indexh
    integer       i,j,k,p,q,t,s,w,ii,jj,kk
    character(1)  c1,cc1
    integer       StrLen,NumFactY
    integer       TotNumOfYTerms,TotNumOfYHYTerms,CurrNumOfTerms
    integer       L,R,FirstLPos,LastRPos,MaxNumTermsInFact
    integer,allocatable,dimension(:)      :: TempSymCoeff,TempSymCoeff1,NumTermsInYOpFact
    integer,allocatable,dimension(:,:,:)  :: TempSymMatr,TempSymMatr1
    integer,allocatable,dimension(:,:)    :: Matr1,Matr2,Matr3,Matr4
    integer                               :: Coeff,Cf3
    character(Glob_YOperatorStringLength),allocatable,dimension(:) :: YOpStr, YHOpStr
    integer       pi,pj,pt,ps
    logical       AreTermsIdentical
    integer,allocatable,dimension(:)      :: IdentParticleSet
    integer,allocatable,dimension(:,:)    :: IdentPseudoPartPairSet
    real(wp)   qi,qj,qq,sqq
    real(wp)   mk,mi,m0,mh,ml

    if (Glob_ProcID==0) write(*,*) 'Initializing program data'

!Initialization of some machine-dependent numbers
    Glob_AbsTolForDSYGVX=2*DLAMCH('S')

    n=Glob_n
    npart=n+1

!Constructing Glob_PseudoChargeMatrix and Glob_ScaledPseudoChargeMatrix
    allocate(Glob_PseudoChargeMatrix(0:n,0:n))
    allocate(Glob_ScaledPseudoChargeMatrix(0:n,0:n))
    Glob_PseudoChargeMatrix(0:n,0:n)=ZERO
    Glob_ScaledPseudoChargeMatrix(0:n,0:n)=ZERO
    do j=0,n
      if (j==0) then
        qj=Glob_PseudoCharge0
      else
        qj=Glob_PseudoCharge(j)
      endif
      do i=0,n
        if (i==0) then
          qi=Glob_PseudoCharge0
        else
          qi=Glob_PseudoCharge(i)
        endif
        qq=qi*qj
        if (qq<0.0_wp) then
          sqq=qq*Glob_AttractionScalingParam
        else
          if ((qi>0.0_wp).and.(qj>0.0_wp)) then
            sqq=qq*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamPlus
          else
            sqq=qq*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamMinus
          endif
        endif
        Glob_PseudoChargeMatrix(i,j)=qq
        Glob_ScaledPseudoChargeMatrix(i,j)=sqq        
      enddo
    enddo    

!Constructing Glob_MassMatrix
    allocate(Glob_MassMatrix(n,n))
    Glob_MassMatrix(1:n,1:n)=ONEHALF/Glob_Mass(1)
    do i=1,n
      Glob_MassMatrix(i,i)=Glob_MassMatrix(i,i)+ONEHALF/Glob_Mass(i+1)
    enddo

!Determine the components of vector Glob_bvc
!that is used in evaluation of particle densities
    allocate(Glob_bvc(n,npart))
    Glob_MassTotal=sum(Glob_Mass(1:npart))
    do i=1,npart
      Glob_bvc(1:n,i)=-Glob_Mass(2:n+1)/Glob_MassTotal
    enddo
    do i=2,npart
      Glob_bvc(i-1,i)=Glob_bvc(i-1,i)+ONE
    enddo

!Determine the mass and the index of the lightest and heaviest particle
!(reference particle excluded).
!and the index if heaviest
    indexh=0
    mh=0
    ml=2*Glob_MassTotal
    do i=1,n
      if (Glob_Mass(i+1)<ml) then
        ml=Glob_Mass(i+1)
      endif
      if (Glob_Mass(i+1)>mh) then
        mh=Glob_Mass(i+1)
        indexh=i
      endif
    enddo

    if (abs(ml - mh) > 1.d-14) then
      Glob_ArePseudoParticleMassesTheSame = .false.
    endif

    m0=Glob_Mass(1)
    Glob_dmva21 = (m0**3 + ml**3)/(TWO*m0*ml*(m0+ml)**2)
    Glob_dmva22 = (m0**3 + mh**3)/(TWO*m0*mh*(m0+mh)**2)
!Glob_dmva2 = (m0**2)/(TWO*mh*(m0+mh)**2)
    Glob_dmvB(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO
!if (.not. Glob_ArePseudoParticleMassesTheSame) then
!  do i=1,n
!    if (i == indexh) cycle
!    mi=Glob_Mass(i+1)
!    mui = m0*mi/(m0+mi)
!    Glob_dmvB(i,i) = -ONE/mui
!  enddo
!  do i=1,n
!    do j=i+1,n
!      Glob_dmvB(i,j) = -ONE/m0
!      Glob_dmvB(j,i) = Glob_dmvB(i,j)
!    enddo
!  enddo
!endif
    Glob_dmvM(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO
    Glob_dmvM(1:n,1:n)=Glob_MassMatrix(1:n,1:n)
!Glob_dmvMB=Glob_dmvM+Glob_dmvB

!Constructing all Pij transposition matrices
!
!  P1i  (i/=1) has the following view:
!
!                    i-1
!
!      |  1   0   0  -1   0   0   0  |
!      |                             |
!      |  0   1   0  -1   0   0   0  |
!      |                             |
!      |  0   0   1  -1   0   0   0  |
!      |                             |          This particular
!P1i = |  0   0   0  -1   0   0   0  |  i-1     matrix is P15 for the
!      |                             |          case of 8 particles
!      |  0   0   0  -1   1   0   0  |
!      |                             |
!      |  0   0   0  -1   0   1   0  |
!      |                             |
!      |  0   0   0  -1   0   0   1  |
!
!
!  Pij has the following form:
!
!                i-1         j-1
!
!      |  1   0   0   0   0   0   0  |
!      |                             |
!      |  0   1   0   0   0   0   0  |
!      |                             |
!      |  0   0   0   0   0   1   0  |  i-1
!      |                             |
!Pij = |  0   0   0   1   0   0   0  |          This particular
!      |                             |          matrix is P47 for the
!      |  0   0   0   0   1   0   0  |          case of 8 particles
!      |                             |
!      |  0   0   1   0   0   0   0  |  j-1
!      |                             |
!      |  0   0   0   0   0   0   1  |
!
    allocate(Glob_Transposit(n,n,npart,npart))
!First set all of them to be unit matrices
    Glob_Transposit(1:n,1:n,1:npart,1:npart)=0
    do i=1,npart
      do j=1,npart
        do k=1,n
          Glob_Transposit(k,k,i,j)=1
        enddo
      enddo
    enddo
!Now continue depending on type of transposition (P1i or Pij)
    do i=2,npart
      Glob_Transposit(1:n,i-1,1,i)=-1
    enddo
    do i=2,npart
      do j=i+1,npart
        Glob_Transposit(i-1,i-1,i,j)=0
        Glob_Transposit(j-1,j-1,i,j)=0
        Glob_Transposit(j-1,i-1,i,j)=1
        Glob_Transposit(i-1,j-1,i,j)=1
        Glob_Transposit(i-1,i-1,j,i)=0
        Glob_Transposit(j-1,j-1,j,i)=0
        Glob_Transposit(j-1,i-1,j,i)=1
        Glob_Transposit(i-1,j-1,j,i)=1
      enddo
    enddo

!Constructing the Young operator based on the content of
!a string variable Glob_YOperatorString

!First we throw away spaces and multiplication signs
!from Glob_YOperatorString
    StrLen=len_trim(Glob_YOperatorString)
    do i=1,StrLen
      c1=Glob_YOperatorString(i:i)
      if ((c1==' ').or.(c1=='*')) then
        do j=i,StrLen-1
          Glob_YOperatorString(j:j)=Glob_YOperatorString(j+1:j+1)
        enddo
        Glob_YOperatorString(j:j)=' '
      endif
    enddo
    StrLen=len_trim(Glob_YOperatorString)

!Checking for wrong symbols in Glob_YOperatorString
    do i=1,StrLen
      c1=Glob_YOperatorString(i:i)
      if ((c1/='1').and.(c1/='2').and.(c1/='3').and.(c1/='4').and.(c1/='5').and. &
          (c1/='6').and.(c1/='7').and.(c1/='8').and.(c1/='9').and.(c1/='0').and. &
          (c1/='P').and.(c1/='+').and.(c1/='-').and.(c1/='*').and.(c1/=')').and. &
          (c1/='(')) then
        write(*,*) 'Error EC0115 in ProgramDataInit: the Young operator expression'
        write(*,*) 'contains wrong symbols'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
    enddo

!Checking if the number of left and right brackets is the same
!and counting how many brackets there are
    L=0
    R=0
    do i=1,StrLen
      if (Glob_YOperatorString(i:i)==')') R=R+1
      if (Glob_YOperatorString(i:i)=='(') L=L+1
    enddo
    if (R/=L) then
      write(*,*) 'Error EC0116 in ProgramDataInit: the numer of left and right brackets in the'
      write(*,*) 'Young operator is different'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!NumFactY is the number of factors in the Young operator,
!FirstLPos is the position of the first left bracket,
!LastRPos is the position of the last right bracket,
    if (R/=0) then
      FirstLPos=scan(Glob_YOperatorString(1:StrLen),'(')
      LastRPos=scan(Glob_YOperatorString(1:StrLen),')',back=.true.)
      NumFactY=R
      i=0
      do j=1,R
        k=0
        i=i+1
        c1=Glob_YOperatorString(i:i)
        do while (c1/='(')
          if (c1/='*') k=1
          i=i+1
          c1=Glob_YOperatorString(i:i)
        enddo
        if (k==1) NumFactY=NumFactY+1
        do while (c1/=')')
          i=i+1
          c1=Glob_YOperatorString(i:i)
        enddo
      enddo
      if (Glob_YOperatorString(StrLen:StrLen)/=')') NumFactY=NumFactY+1
    else
      NumFactY=1
    endif

!Splitting Glob_YOperatorString into an array of smaller
!strings, YOpStr. Each column of this array will contain just
!one factor, with no brackets. A '+' or a '-' sign is added in
!front of the first term in a factor if needed. Multiplication
!signs are dropped .
    allocate(YOpStr(NumFactY))
    do i=1,NumFactY
      YOpStr(i)=' '
    enddo
    if (R==0) then
      c1=Glob_YOperatorString(1:1)
      if ((c1/='+').or.(c1/='-')) then
        YOpStr(1)(1:1)='+'
        YOpStr(1)(2:StrLen+1)=Glob_YOperatorString(1:StrLen)
      else
        YOpStr(1)(1:StrLen)=Glob_YOperatorString(1:StrLen)
      endif
    else
      i=1
      k=1
      p=i
      q=0
      c1=Glob_YOperatorString(i:i)
      if ((c1/='(').and.(c1/='+').and.(c1/='-')) then
        q=1
        YOpStr(k)(1:1)='+'
      endif
      do while (Glob_YOperatorString(i:i)/='(')
        i=i+1
      enddo
      if (i>1) then
        YOpStr(k)(p+q:i-1+q)=Glob_YOperatorString(p:i-1)
        !if (YOpStr(k)(i-1+q:i-1+q)=='*') YOpStr(k)(i-1+q:i-1+q)=' '
        k=k+1
      endif
      do j=1,R
        i=i+1
        p=i
        q=0
        c1=Glob_YOperatorString(i:i)
        if ((c1/=')').and.(c1/='+').and.(c1/='-')) then
          q=1
          YOpStr(k)(1:1)='+'
        endif
        do while (Glob_YOperatorString(i:i)/=')')
          i=i+1
        enddo
        YOpStr(k)(1+q:i+q-p)=Glob_YOperatorString(p:i-1)
        k=k+1
        i=i+1
        c1=Glob_YOperatorString(i:i)
        if (c1=='*') then
          i=i+1
          c1=Glob_YOperatorString(i:i)
        endif
        if ((c1/='(').and.(c1/=' '))  then
          p=i
          YOpStr(k)(1:1)='+'
          do while ((c1/='(').and.(c1/=' '))
            i=i+1
            c1=Glob_YOperatorString(i:i)
          enddo
          YOpStr(k)(2:i+1-p)=Glob_YOperatorString(p:i-1)
          k=k+1
        endif
      enddo
    endif

!Print all factors in the Young operator
!j=StrLen+1
!do i=1,NumFactY
!  write (*,'(1x,i3,1x,a3,a<j>)') i,':  ',YOpStr(i)(1:j)
!enddo

!Creating an array that contains all the factors of the
!Y^{\dagger} operator. Basically, Y^{\dagger} is the reversed Y (i.e.
!the order of all factors is reversed as well as
!permutation products (if there are any) in each factor come
!in reverse order.
    allocate(YHOpStr(NumFactY))
    do i=1,NumFactY
      YHOpStr(i)=' '
    enddo
    do i=NumFactY,1,-1
      s=NumFactY-i+1
      j=1
      c1=YOpStr(s)(j:j)
      do while (c1/=' ')
        if (c1=='P') then
          k=0 !k counts the number of Permutations in the current term
          t=0
          do while ((c1/='+').and.(c1/='-').and.(c1/=' '))
            if (c1=='P') k=k+1
            t=t+1
            c1=YOpStr(s)(j+t:j+t)
          enddo
          do t=1,k
            YHOpStr(i)(j+3*(k-t):j+3*(k-t)+2)=YOpStr(s)(j+3*(t-1):j+3*(t-1)+2)
          enddo
          j=j+3*k
        else
          YHOpStr(i)(j:j)=c1
          j=j+1
        endif
        c1=YOpStr(s)(j:j)
      enddo
    enddo

!Counting how many terms there are in each factor of the Young
!operator, as well as the total number of terms in the nonsimplified
!Young operator
    allocate(NumTermsInYOpFact(NumFactY))
    TotNumOfYTerms=1
    do k=1,NumFactY
      j=0
      do i=1,Glob_YOperatorStringLength
        if ((YOpStr(k)(i:i)=='+').or.(YOpStr(k)(i:i)=='-')) j=j+1
      enddo
      NumTermsInYOpFact(k)=j
      TotNumOfYTerms=TotNumOfYTerms*j
      TotNumOfYHYTerms=TotNumOfYTerms*TotNumOfYTerms
    enddo
    if (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y operator:     ',TotNumOfYTerms
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y^{+}Y operator:',TotNumOfYHYTerms
    endif

    allocate(Matr1(1:n,1:n))
    allocate(Matr2(1:n,1:n))
    allocate(Matr3(1:n,1:n))
    allocate(Matr4(1:n,1:n))

!Multiplying all factors in YOpStr and placing actual matrices
!and coefficients in arrays Glob_YMatr and Glob_YCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.
    CurrNumOfTerms=NumTermsInYOpFact(NumFactY)
    allocate(TempSymCoeff(CurrNumOfTerms))
    allocate(TempSymMatr(n,n,CurrNumOfTerms))
    do j=NumFactY,1,-1
      !reading the current factor
      k=0
      i=1
      c1=YOpStr(j)(i:i)
      p=i
      do while (c1/=' ')
        i=i+1
        c1=YOpStr(j)(i:i)
        do while ((c1/='P').and.(c1/='+').and.(c1/='-').and.(i<Glob_YOperatorStringLength))
          i=i+1
          c1=YOpStr(j)(i:i)
        enddo
        if (i-p>1) then
          read(YOpStr(j)(p:i-1),*) Coeff
        else
          if (YOpStr(j)(i-1:i-1)=='+') then
            Coeff=1
          else
            Coeff=-1
          endif
        endif
        Matr1=Glob_Transposit(1:n,1:n,1,1)
        do while (c1=='P')
          read(YOpStr(j)(i+1:i+1),*) p
          read(YOpStr(j)(i+2:i+2),*) q
          Matr2(1:n,1:n)=Glob_Transposit(1:n,1:n,p,q)
          Matr4(1:n,1:n)=Matr1(1:n,1:n)
          do ii=1,n
            do jj=1,n
              w=0
              do kk=1,n
                w=w+Matr2(ii,kk)*Matr4(kk,jj)
              enddo
              Matr1(ii,jj)=w
            enddo
          enddo
          i=i+3
          c1=YOpStr(j)(i:i)
        enddo
        k=k+1
        p=i
        if (j/=NumFactY) then
          if (k==1) then
            Matr3(1:n,1:n)=Matr1(1:n,1:n)
            Cf3=Coeff
          else
            do s=1,t
              Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
              q=t*(k-1)+s
              do ii=1,n
                do jj=1,n
                  w=0
                  do kk=1,n
                    w=w+Matr2(ii,kk)*Matr1(kk,jj)
                  enddo
                  TempSymMatr(ii,jj,q)=w
                enddo
              enddo
              TempSymCoeff(q)=Coeff*TempSymCoeff(s)
            enddo
          endif
        else
          TempSymMatr(1:n,1:n,k)=Matr1(1:n,1:n)
          TempSymCoeff(k)=Coeff
        endif
      enddo
      if (j/=NumFactY) then
        do s=1,t
          Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
          do ii=1,n
            do jj=1,n
              w=0
              do kk=1,n
                w=w+Matr2(ii,kk)*Matr3(kk,jj)
              enddo
              TempSymMatr(ii,jj,s)=w
            enddo
          enddo
          TempSymCoeff(s)=Cf3*TempSymCoeff(s)
        enddo
      endif
      !mark the identical terms (adding their coefficients
      !and setting all of them but one to zero)
      t=CurrNumOfTerms
      do i=1,CurrNumOfTerms
        if (TempSymCoeff(i)==0) cycle
        do s=i+1,CurrNumOfTerms
          if (TempSymCoeff(s)==0) cycle
          if (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
            TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
            if (TempSymCoeff(i)==0) t=t-1
            TempSymCoeff(s)=0
            t=t-1
          endif
        enddo
      enddo
      !reallocate arrays containing symmetry terms
      !to allow for multiplication by the next factor
      if (j/=1) then
        allocate(TempSymCoeff1(t))
        allocate(TempSymMatr1(n,n,t))
        s=0
        do i=1,CurrNumOfTerms
          if (TempSymCoeff(i)/=0) then
            s=s+1
            TempSymCoeff1(s)=TempSymCoeff(i)
            TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
          endif
        enddo
        CurrNumOfTerms=t*NumTermsInYOpFact(j-1)
        deallocate(TempSymCoeff)
        deallocate(TempSymMatr)
        allocate(TempSymCoeff(CurrNumOfTerms))
        allocate(TempSymMatr(n,n,CurrNumOfTerms))
        TempSymCoeff(1:t)=TempSymCoeff1(1:t)
        TempSymMatr(1:n,1:n,1:t)=TempSymMatr1(1:n,1:n,1:t)
        deallocate(TempSymCoeff1)
        deallocate(TempSymMatr1)
      endif
    enddo

    Glob_NumYTerms=t
    allocate(Glob_YCoeff(Glob_NumYTerms))
    allocate(Glob_YMatr(n,n,Glob_NumYTerms))
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        Glob_YCoeff(s)=TempSymCoeff(i)
        Glob_YMatr(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      endif
    enddo
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    if (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in simplified Y operator:        ',Glob_NumYTerms
    endif

!Now doing the same thing for Y^{\dagger}Y operator, that
!is expanding it and collecting identical terms

!Multiplying all factors in YHOpStr by already existing
!matrices and coefficients of Y. and placing actual matrices
!and coefficients in arrays Glob_YHYMatr and Glob_YHYCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.
    CurrNumOfTerms=NumTermsInYOpFact(1)*Glob_NumYTerms
    allocate(TempSymCoeff(CurrNumOfTerms))
    allocate(TempSymMatr(n,n,CurrNumOfTerms))
!TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
!TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)
    TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
    TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)
    t=Glob_NumYTerms
    do j=NumFactY,1,-1
      !reading the current factor
      k=0
      i=1
      c1=YHOpStr(j)(i:i)
      p=i
      do while (c1/=' ')
        i=i+1
        c1=YHOpStr(j)(i:i)
        do while ((c1/='P').and.(c1/='+').and.(c1/='-').and.(i<Glob_YOperatorStringLength))
          i=i+1
          c1=YHOpStr(j)(i:i)
        enddo
        if (i-p>1) then
          read(YHOpStr(j)(p:i-1),*) Coeff
        else
          if (YHOpStr(j)(i-1:i-1)=='+') then
            Coeff=1
          else
            Coeff=-1
          endif
        endif
        Matr1=Glob_Transposit(1:n,1:n,1,1)
        do while (c1=='P')
          read(YHOpStr(j)(i+1:i+1),*) p
          read(YHOpStr(j)(i+2:i+2),*) q
          Matr2(1:n,1:n)=Glob_Transposit(1:n,1:n,p,q)
          Matr4(1:n,1:n)=Matr1(1:n,1:n)
          do ii=1,n
            do jj=1,n
              w=0
              do kk=1,n
                w=w+Matr2(ii,kk)*Matr4(kk,jj)
              enddo
              Matr1(ii,jj)=w
            enddo
          enddo
          i=i+3
          c1=YHOpStr(j)(i:i)
        enddo
        k=k+1
        p=i
        if (k==1) then
          Matr3(1:n,1:n)=Matr1(1:n,1:n)
          Cf3=Coeff
        else
          do s=1,t
            Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
            q=t*(k-1)+s
            do ii=1,n
              do jj=1,n
                w=0
                do kk=1,n
                  w=w+Matr2(ii,kk)*Matr1(kk,jj)
                enddo
                TempSymMatr(ii,jj,q)=w
              enddo
            enddo
            TempSymCoeff(q)=Coeff*TempSymCoeff(s)
          enddo
        endif
      enddo
      do s=1,t
        Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
        do ii=1,n
          do jj=1,n
            w=0
            do kk=1,n
              w=w+Matr2(ii,kk)*Matr3(kk,jj)
            enddo
            TempSymMatr(ii,jj,s)=w
          enddo
        enddo
        TempSymCoeff(s)=Cf3*TempSymCoeff(s)
      enddo
      !mark the identical terms (adding their coefficients
      !and setting all of them but one to zero)
      t=CurrNumOfTerms
      do i=1,CurrNumOfTerms
        if (TempSymCoeff(i)==0) cycle
        do s=i+1,CurrNumOfTerms
          if (TempSymCoeff(s)==0) cycle
          if (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
            TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
            if (TempSymCoeff(i)==0) t=t-1
            TempSymCoeff(s)=0
            t=t-1
          endif
        enddo
      enddo
      !reallocate arrays containing symmetry terms
      !to allow for multiplication by the next factor
      if (j/=1) then
        allocate(TempSymCoeff1(t))
        allocate(TempSymMatr1(n,n,t))
        s=0
        do i=1,CurrNumOfTerms
          if (TempSymCoeff(i)/=0) then
            s=s+1
            TempSymCoeff1(s)=TempSymCoeff(i)
            TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
          endif
        enddo
        CurrNumOfTerms=t*NumTermsInYOpFact(NumFactY-j+2)
        deallocate(TempSymCoeff)
        deallocate(TempSymMatr)
        allocate(TempSymCoeff(CurrNumOfTerms))
        allocate(TempSymMatr(n,n,CurrNumOfTerms))
        TempSymCoeff(1:t)=TempSymCoeff1(1:t)
        TempSymMatr(1:n,1:n,1:t)=TempSymMatr1(1:n,1:n,1:t)
        deallocate(TempSymCoeff1)
        deallocate(TempSymMatr1)
      endif
    enddo

    Glob_NumYHYTerms=t
    allocate(Glob_YHYCoeff(Glob_NumYHYTerms))
    allocate(Glob_YHYMatr(n,n,Glob_NumYHYTerms))
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        Glob_YHYCoeff(s)=TempSymCoeff(i)
        Glob_YHYMatr(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      endif
    enddo
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    if (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in simplified Y^{+}Y operator:   ',Glob_NumYHYTerms
    endif

!Print all independent Y operator matrices and coefficients
!open(1,file='symterms_new.txt',status='replace')
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYTerms
!    write(1,'(10x,a2,i3,a4,i3)') 'i=',i,'  C=',Glob_YCoeff(i)
!    do j=1,n
!      write(1,'(<n>(1x,i2))') Glob_YMatr(j,1:n,i)
!    enddo
!    write(1,*)
!  enddo
!endif
!write(1,*)
!write(1,*) '-----'
!Print all independent Y^{\dagger}Y operator matrices and coefficients
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYHYTerms
!    write(1,'(10x,a2,i3,a4,i3)') 'i=',i,'  C=',Glob_YHYCoeff(i)
!    do j=1,n
!      write(1,'(<n>(1x,i2))') Glob_YHYMatr(j,1:n,i)
!    enddo
!    write(1,*)
!  enddo
!endif
!close(1)
!stop

!open(1,file='symterms_new.txt',status='replace')
!if (Glob_ProcID==0) then
!  write(1,*)
!  do i=1,Glob_NumYHYTerms
!    write(1,'(1x,a14,i2,a2,E41.34)') 'Glob_YHYCoeff(',i,')=',Glob_YHYCoeff(i)
!  enddo
!  do i=1,Glob_NumYHYTerms
!    do j=1,n
!      do k=1,n
!        write(1,'(1x,a13,i2,a1,i2,a1,i2,a2,e23.16)') &
!              'Glob_YHYMatr(', k , ',' , j , ',' , i, ')=', Glob_YHYMatr(k,j,i)
!      enddo
!    enddo
!  enddo
!endif
!close(1)
!stop

    deallocate(Matr4)
    deallocate(Matr3)
    deallocate(Matr2)
    deallocate(Matr1)
    deallocate(NumTermsInYOpFact)
    deallocate(YHOpStr)
    deallocate(YOpStr)

!Now we determine which particles are identical. This determination
!is based on the input values of masses and charges only. The information
!about the sets of identical particles may be needed for proper
!symmetrization of expectation values of operators that involve
!two-particle quantities (such as interparticle distances).
!The set of particles to which particle i belongs is labelled by IdentParticleSet(i)
!If IdentParticleSet(i)=IdentParticleSet(j) then it means that
!particles i and j are identical. The largest value in IdentParticleSet
!gives the total number of identical particle sets
    allocate(IdentParticleSet(npart))
    IdentParticleSet(1)=1
    k=1
    do i=2,npart
      s=0
      j=0
      do while ((j<i-1).and.(s==0))
        j=j+1
        if (j>1) then
          if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge(j-1)==Glob_PseudoCharge(i-1))) then
            IdentParticleSet(i)=IdentParticleSet(j)
            s=1
          endif
        else
          !j=1 case
          if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge0==Glob_PseudoCharge(i-1))) then
            IdentParticleSet(i)=IdentParticleSet(j)
            s=1
          endif
        endif
      enddo
      if (s==0) then
        k=k+1
        IdentParticleSet(i)=k
      endif
    enddo
    Glob_NumOfIdentPartSets=maxval(IdentParticleSet(1:npart))

!Below we determine which pairs of pseudoparticles are identical.
!The information about this is stored in array IdentPseudoPartPairSet(1:n,1:n)
!Diagonal elements do not actually designate pairs of pseudoparticles but
!rather a single pseudoparticle, which corresponds to a certain pair of particles.
!If IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(k,l) then
!it means these pairs ij and kl should be equivalent.
!The largest value of array IdentPseudoPartPairSet gives the number
!of nonequivalent pairs.
    allocate(IdentPseudoPartPairSet(1:n,1:n))
    IdentPseudoPartPairSet(1:n,1:n)=0
    k=0
    do i=1,n
      do j=i,n
        if (i==j) then
          pi=1; pj=j+1
        else
          pi=i+1; pj=j+1
        endif
        w=0
        do s=1,i
          if (s==i) then
            q=j-1
          else
            q=n
          endif
          do t=s,q
            if (w==1) cycle
            if (s==t) then
              ps=1; pt=t+1
            else
              ps=s+1; pt=t+1
            endif
            if ((IdentParticleSet(ps)==IdentParticleSet(pi)).and. &
                (IdentParticleSet(pt)==IdentParticleSet(pj))) then
              w=1
              IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(s,t)
            endif
          enddo
        enddo
        if (w==0) then
          k=k+1
          IdentPseudoPartPairSet(i,j)=k
        endif
      enddo
    enddo
    Glob_NumOfNoneqvPairSets=maxval(IdentPseudoPartPairSet(1:n,1:n))

!Now we create arrays Glob_NumOfPartInIdentPartSet
!and Glob_IdentPartList
    allocate(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
    allocate(Glob_IdentPartList(npart,Glob_NumOfIdentPartSets))
    Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)=0
    Glob_IdentPartList(1:npart,1:Glob_NumOfIdentPartSets)=0
    do i=1,npart
      k=IdentParticleSet(i)
      Glob_NumOfPartInIdentPartSet(k)=Glob_NumOfPartInIdentPartSet(k)+1
      Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k),k)=i
    enddo

!Create arrays Glob_NumOfPairsInEqvPairSet and
!Glob_EqvPairList
    allocate(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
    allocate(Glob_EqvPairList(2,n*(n+1)/2,Glob_NumOfNoneqvPairSets))
    Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets)=0
    Glob_EqvPairList(1:2,1:n*(n+1)/2,1:Glob_NumOfNoneqvPairSets)=0
    do i=1,n
      do j=i,n
        k=IdentPseudoPartPairSet(i,j)
        Glob_NumOfPairsInEqvPairSet(k)=Glob_NumOfPairsInEqvPairSet(k)+1
        Glob_EqvPairList(1,Glob_NumOfPairsInEqvPairSet(k),k)=i
        Glob_EqvPairList(2,Glob_NumOfPairsInEqvPairSet(k),k)=j
      enddo
    enddo

    deallocate(IdentPseudoPartPairSet)
    deallocate(IdentParticleSet)

  end subroutine ProgramDataInit

  subroutine GenerateTrialParam(nfun,x,m,k,method_used)
!Subroutine GenerateTrialParam generates nonlinear
!parameters for trial functions. It uses three global variables
!to do that: Glob_RG_p1, Glob_RG_s1, Glob_RG_s2.
!Glob_RG_p1 is the probability of using the first generation method,
!so 1-Glob_RG_p1 is the probability of using the second one. Both
!methods select a random function from the existing basis (or several
!functions if nfun>1). Then, the first method generates nonlinear
!parameters that are independently normally distributed around the
!corresponding parameters of the selected function with standard
!deviation Glob_RG_s1*x(i), where x(i) is the i-th nonlinear parameter
!of the selected function. The second method generates nonlinear
!parameters in the following way: it generates a random number, r, from
!a normal (0,Glob_RG_s2) distribution and then multiplies all x(i) by (1+r).
!It also makes sure that 0.8<|1+r|>1.2 as we do not want almost linearly
!dependent functions in the basis. The generation of Z-indices is done in the
!following manner. If the first generation method is used then ZIndexChangeProb*100%
!probability is that a Z-index will be the same as in the prototype function,
!the rest 100%-ZIndexChangeProb*100% of the probability is uniformly distributed
!among all other possible values of Z-index. If the second generation method is used
!then all Z-indices are always the same as in the prototype functions.
!  Input parameter:
!   nfun - number of trial functions whose parameters
!          are to be generated;
!  Output parameters :
!    x(1:Glob_npt,1:nfun) - a 2D-array containing generated
!          nonlinear parameters of the trial functions;
!    m(1:nfun,2) - an 2D-array containing generated
!          x,y-indices of the trial functions;
!      k - specifies which function in the existing
!          basis was selected to generate a candidate. In case when
!          nfun>1 k specifies which was the first function (out of
!          a set) selected. In case Glob_CurrBasisSize==0 k is set to 0.
!    method_used - specifies which generation method was used (1 or 2).
!          In case Glob_CurrBasisSize==0 method_used is set to 0.

!Arguments :
    integer        nfun
    real(wp)    x(Glob_npt,nfun)
    integer        m(nfun,2)
    integer        k, method_used

!Local variables :
    integer        i,j,jj,p
    real(8)    r,sumf,r1
!Constants that define the uniform distribution
!in the case of Glob_CurrBasisSize==0
    real(8)  :: Lmin=-0.5_8
    real(8)  :: Lmax= 0.5_8
    real(8)  :: IndexChangeProb=0.5_8

    if (Glob_CurrBasisSize==0) then !making uniform distribution
      do i=1,nfun
        do j=1,Glob_np
          call random_number(r)
          x(j,i)=r*(Lmax-Lmin)+Lmin
        enddo
        call random_number(r)
        call random_number(r1)
        r=r*Glob_n
        r1=r1*Glob_n
        m(i,1)=1+int(r)
        m(i,2)=1+int(r1)
        if (m(i,1)>Glob_n) m(i,1)=Glob_n
        if (m(i,2)>Glob_n) m(i,2)=Glob_n
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if (m(i,1)==m(i,2)) then
!         if (m(i,1)==Glob_n) then
!             m(i,2)=m(i,1)-1
!         elseif (m(i,1)<Glob_n) then
!             m(i,2)=m(i,1)+1
!         endif
!    end if
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      method_used=0
      k=0
    else
      if (Glob_CurrBasisSize<nfun) then
        call random_number(r)
        if (r<GLob_RG_p1) then
          !method 1
          do i=1,nfun
            call random_number(r)
            k=int(r*(Glob_CurrBasisSize))+1
            do j=1,Glob_npt
              x(j,i)=(Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j,k)
            enddo
            call random_number(r)
            if (r>IndexChangeProb) then
              m(i,1)=Glob_Index(k,1)
              m(i,2)=Glob_Index(k,2)
            else
              j=Glob_Index(k,1)
              jj=Glob_Index(k,2)
              m(i,1)=j
              m(i,2)=jj
              do while (m(i,1)==j .and. m(i,2)==jj) !.and. m(i,1)==m(i,2)
                call random_number(r)
                r=r*Glob_n
                m(i,1)=1+int(r)
                call random_number(r)
                r=r*Glob_n
                m(i,2)=1+int(r)
              enddo
            endif
          enddo
          method_used=1
        else
          !method 2
          do i=1,nfun
            call random_number(r)
            k=int(r*(Glob_CurrBasisSize))+1
            r=Glob_RG_s2*drnor()+ONE
            do while ((abs(r)>0.8E0_wp).and.(abs(r)<1.2E0_wp))
              r=Glob_RG_s2*drnor()+ONE
            enddo
            do j=1,Glob_npt
              x(j,i)=r*Glob_NonlinParam(j,k)
            enddo
            m(i,1)=Glob_Index(k,1)
            m(i,2)=Glob_Index(k,2)
          enddo
          method_used=2
        endif
      else
        call random_number(r)
        k=int(r*(Glob_CurrBasisSize-nfun+1))+1
        call random_number(r)
        if (r<GLob_RG_p1) then
          !method 1
          do i=1,nfun
            do j=1,Glob_npt
              x(j,i)=(Glob_RG_s1*drnor()+ONE)*Glob_NonlinParam(j,k+i-1)
            enddo
            call random_number(r)
            if ((r>IndexChangeProb).or.(Glob_n==1)) then
              m(i,1)=Glob_Index(k+i-1,1)
              m(i,2)=Glob_Index(k+i-1,2)
            else
              j=Glob_Index(k+i-1,1)
              jj=Glob_Index(k+i-1,2)
              m(i,1)=j
              m(i,2)=jj
              do while (m(i,1)==j .and. m(i,2)==jj)  ! .and. m(i,1)==m(i,2)
                call random_number(r)
                r=r*Glob_n
                m(i,1)=1+int(r)
                call random_number(r)
                r=r*Glob_n
                m(i,2)=1+int(r)
              enddo
            endif
          enddo
          method_used=1
        else
          !method 2  !in method two z indices will not changed
          r=Glob_RG_s2*drnor()+ONE
          do while ((abs(r)>0.8E0_wp).and.(abs(r)<1.2E0_wp))
            r=Glob_RG_s2*drnor()+ONE
          enddo
          do i=1,nfun
            do j=1,Glob_npt
              x(j,i)=r*Glob_NonlinParam(j,k+i-1)
            enddo
            m(i,1)=Glob_Index(k+i-1,1)
            m(i,2)=Glob_Index(k+i-1,2)
          enddo
          method_used=2
        endif
      endif
    endif

    !If any of the indices are fixed then we set them to the fixed values
    if (Glob_IsIndexFixed(1)) then
      m(1:nfun,1)=Glob_IndexFixedValue(1)
    endif
    if (Glob_IsIndexFixed(2)) then
      m(1:nfun,2)=Glob_IndexFixedValue(2)
    endif

    !If Glob_VectorCouplingScheme is nonzero then we need to make sure that the generated indices
    !are compatible with the coupling scheme. Note that ReadIOFile does not allow both indices to
    !be fixed when a coupling scheme is specified, so there is always at least one index that is
    !free to be adjusted here. When neither index is fixed we adjust the first one.
    !Glob_VectorCouplingScheme==1 : two particles in p-states, so the two indices must be different
    if (Glob_VectorCouplingScheme==1) then
      do i=1,nfun
        do while (m(i,1)==m(i,2))
          call random_number(r)
          j=1+int(r*Glob_n)
          if (j>Glob_n) j=Glob_n
          if (Glob_IsIndexFixed(1)) then
            m(i,2)=j
          else
            m(i,1)=j
          endif
        enddo
      enddo
    endif
    !Glob_VectorCouplingScheme==2 : a single particle in a d-state, so the two indices must be the same
    if (Glob_VectorCouplingScheme==2) then
      if (Glob_IsIndexFixed(2)) then
        m(1:nfun,1)=m(1:nfun,2)
      else
        m(1:nfun,2)=m(1:nfun,1)
      endif
    endif

  end subroutine GenerateTrialParam

  subroutine ComputeOverlapPenalty(MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty)
!Subroutine ComputeOverlapPenalty computes overlap penalty. Overlap
!penalty may be used during full optimization of all nonlinear parameters
!in order to prevent severe pair linear dependencies of basis functions.
!The penalty function is P = sum_{i=1,K; j=Glob_nfru+1,K; i<j} Pij, where
!  Pij=(abs(Sij)^2-t^2)*b/(1-t^2),  Sij^2>t^2
!  Pij=0, Sij<=t^2
!Here
!  b=MaxPairOverlapPenalty
!  t^2=OverlapThreshold2
!  Sij is the pair overlap value
!The result is returned in TotalPenalty
!Arguments:
    real(wp)  MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty
!Local variables
    integer      i,j,k,nbands,leftover
    real(wp)  pen_coeff,tp
    logical      oddband

    tp=ZERO
    pen_coeff=MaxPairOverlapPenalty/(ONE-OverlapThreshold2)
    if (Glob_NumOfProcs==1) then
      !In case of a single MPI process we do not split work
      do i=1,Glob_nfru
        do j=Glob_nfru+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) &
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
        enddo
      enddo
      do i=Glob_nfru+1,Glob_nfa
        do j=i+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) &
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
        enddo
      enddo
    else
      !In case of more than one MPI process we split the work
      !more or less evenly between processes
      do i=1+Glob_ProcID,Glob_nfru,Glob_NumOfProcs
        do j=Glob_nfru+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) &
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
        enddo
      enddo
      nbands=Glob_nfo/Glob_NumOfProcs
      leftover=mod(Glob_nfo,Glob_NumOfProcs)
      oddband=.true.
      do k=1,nbands
        if (oddband) then
          i=Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
          oddband=.false.
        else
          i=Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
          oddband=.true.
        endif
        do j=i+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) &
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
        enddo
      enddo !k
      if (leftover>1) then
        if (oddband) then
          i=Glob_nfa-Glob_ProcID+1
        else
          i=Glob_nfa-leftover+1+Glob_ProcID
        endif
        if ((i<Glob_nfa).and.(i>Glob_nfa-leftover)) then
          do j=i+1,Glob_nfa
            if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) &
              tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
          enddo
        endif
      endif
    endif
    call MPI_ALLREDUCE(tp,TotalPenalty,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

  end subroutine ComputeOverlapPenalty

  subroutine ComputeOverlapPenaltyAndAddGradient(MaxPairOverlapPenalty,OverlapThreshold2, &
                                                 TotalPenalty,WkGR)
!Subroutine ComputeOverlapPenalty computes overlap penalty as well as the
!addition to the energy gradient due to the penalty. Overlap
!penalty may be used during full optimization of all nonlinear parameters
!in order to prevent severe pair linear dependencies of basis functions.
!The penalty function is P = sum_{i=1,K; j=Glob_nfru+1,K; i<j} Pij, where
!  Pij=(abs(Sij)^2-t^2)*b/(1-t^2),  Sij^2>t^2
!  Pij=0, Sij<=t^2
!Here
!  b=MaxPairOverlapPenalty
!  t^2=OverlapThreshold2
!  Sij is the pair overlap value
!The overlap penalty is returned in TotalPenalty. The gradient is added to WkGR.
!Note that even though the subroutine divides the work between MPI processes, it does
!not combine the results for the gradient. It is assumed that such an operation will be done
!manually right after ComputeOverlapPenaltyAndAddGradient is called. This is done to avoid
!combining the gradient components twice (once after the gradient is computed and once
!after the gradient addition due to the overlap penalties is computed).
!Arguments:
    real(wp)  MaxPairOverlapPenalty,OverlapThreshold2,TotalPenalty,WkGR(:)
!Local variables
    integer      i,j,k,m,nbands,leftover
    real(wp)  pen_coeff,tp,temp1,temp2
    logical      oddband

    tp=ZERO
    pen_coeff=MaxPairOverlapPenalty/(ONE-OverlapThreshold2)
    if (Glob_NumOfProcs==1) then
      !In case of a single MPI process we do not split work
      do i=1,Glob_nfru
        do j=Glob_nfru+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) then
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
            temp1=2*pen_coeff*Glob_S(j,i)
            temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
            do m=1,Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,j-Glob_nfru,i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2)
            enddo
          endif
        enddo
      enddo
      do i=Glob_nfru+1,Glob_nfa
        do j=i+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) then
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
            temp1=2*pen_coeff*Glob_S(j,i)
            temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
            do m=1,Glob_npt
              WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,i-Glob_nfru,j) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2)
            enddo
            do m=1,Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,j-Glob_nfru,i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2)
            enddo
          endif
        enddo
      enddo
    else
      !In case of more than one MPI process we split the work
      !more or less evenly between processes
      do i=1+Glob_ProcID,Glob_nfru,Glob_NumOfProcs
        do j=Glob_nfru+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) then
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
            temp1=2*pen_coeff*Glob_S(j,i)
            temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
            do m=1,Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,j-Glob_nfru,i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2)
            enddo
          endif
        enddo
      enddo
      nbands=Glob_nfo/Glob_NumOfProcs
      leftover=mod(Glob_nfo,Glob_NumOfProcs)
      oddband=.true.
      do k=1,nbands
        if (oddband) then
          i=Glob_nfru+k*Glob_NumOfProcs-Glob_ProcID
          oddband=.false.
        else
          i=Glob_nfru+1+(k-1)*Glob_NumOfProcs+Glob_ProcID
          oddband=.true.
        endif
        do j=i+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) then
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
            temp1=2*pen_coeff*Glob_S(j,i)
            temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
            do m=1,Glob_npt
              WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,i-Glob_nfru,j) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2)
            enddo
            do m=1,Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,j-Glob_nfru,i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2)
            enddo
          endif
        enddo
      enddo !k
      if (leftover>1) then
        if (oddband) then
          i=Glob_nfa-Glob_ProcID+1
        else
          i=Glob_nfa-leftover+1+Glob_ProcID
        endif
        if ((i<Glob_nfa).and.(i>Glob_nfa-leftover)) then
          do j=i+1,Glob_nfa
          if (Glob_S(j,i)*Glob_S(j,i)>OverlapThreshold2) then
            tp=tp+pen_coeff*(Glob_S(j,i)*Glob_S(j,i)-OverlapThreshold2)
            temp1=2*pen_coeff*Glob_S(j,i)
            temp2=sqrt(Glob_diagS(i)/Glob_diagS(j))
            do m=1,Glob_npt
              WkGR((i-Glob_nfru-1)*Glob_npt+m)=WkGR((i-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,i-Glob_nfru,j) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,i-Glob_nfru,i)*Glob_S(j,i)/temp2)
            enddo
            do m=1,Glob_npt
              WkGR((j-Glob_nfru-1)*Glob_npt+m)=WkGR((j-Glob_nfru-1)*Glob_npt+m) &
                                                +temp1*(Glob_D(Glob_npt+m,j-Glob_nfru,i) &
                                                        -ONEHALF*Glob_D(Glob_npt+m,j-Glob_nfru,j)*Glob_S(j,i)*temp2)
            enddo
          endif
          enddo
        endif
      endif
    endif
    call MPI_ALLREDUCE(tp,TotalPenalty,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

  end subroutine ComputeOverlapPenaltyAndAddGradient

  function EnergyGA(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyGA computes the energy of the system under
!consideration. Routine DSYGVX from LAPACK is used to solve
!GSEP. The function finds the energy level defined by the
!global variable Glob_WhichEigenvalue. The calculations are
!done with a basis of Nmax functions. It is assumed that
!nonlinear parameters of the basis functions are stored in
!global array Glob_NonlinParam and all the matrix elements of the
!first Nmin-1 functions have already been calculated and stored in
!proper global arrays. If no matrix elements have been calculated,
!one should set Nmin=0. In the case when all matrix elements have been
!computed and it is only necessary to solve GSEP one
!needs to set AreMatElemNeeded=.false.
!Upon successful exit ErrorCode should be equal to 0. If
!ErrorCode = Nmax + i then the leading minor of S of size i is
!not positive definite.

    real(wp)    EnergyGA
!Arguments:
    integer        Nmin,Nmax
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables:
    integer        i,j
    real(wp)    Evalue, EVs(1)
    real(wp)    Z(1)
    integer        NumOfEigvalsFound
    integer        IFAIL(1)

    if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
    if (Nmax==1) then
      EnergyGA=Glob_diagH(1)
      ErrorCode=0
    else
      !Copying H and S matrix elements from the lower triangles
      !to the upper ones. The diagonals are copied from the global
      !arrays where they are stored
      do i=1,Nmax
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,Nmax
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo
      if (Glob_ProcID==0) then
        call DSYGVX(1,'N','I','U',Nmax,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,   &
                    ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound,EVs,Z,Nmax,Glob_WorkForDSYGVX,  &
                    Glob_LWorkForDSYGVX,Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
!       LWORK, IWORK, IFAIL, INFO )
        Evalue=EVs(1)
      endif
      if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_OverlapPenaltyAllowed) then
        EnergyGA=Evalue+Glob_TotalOverlapPenalty
      else
        EnergyGA=Evalue
      endif
    endif

    Glob_EnergyGACounter=Glob_EnergyGACounter+1

  end function EnergyGA

  function EnergyGAM(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyGAM is essentially the same as EnergyGA.
!The only difference is that EnergyGAM also computes the
!linear coefficients (stored in Glob_c).

    real(wp)    EnergyGAM
!Arguments:
    integer        Nmin,Nmax
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables:
    integer        i,j
    real(wp)    Evalue, EVs(1)
    integer        NumOfEigvalsFound
    integer        IFAIL(1)

    if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
    if (Nmax==1) then
      EnergyGAM=Glob_diagH(1)
      Glob_c(1)=ONE
      ErrorCode=0
    else
      !Copying H and S matrix elements from the lower triangles
      !to the upper ones. The diagonals are copied from the global
      !arrays where they are stored
      do i=1,Nmax
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,Nmax
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo
      if (Glob_ProcID==0) then
        call   DSYGVX(1,'V','I','U',Nmax,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,   &
                      ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForDSYGVX, &
                      NumOfEigvalsFound,EVs,Glob_c,Nmax,Glob_WorkForDSYGVX,  &
                      Glob_LWorkForDSYGVX,Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
!       LWORK, IWORK, IFAIL, INFO )
        Evalue=EVs(1)
      endif
      if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_c,Nmax,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_OverlapPenaltyAllowed) then
        EnergyGAM=Evalue+Glob_TotalOverlapPenalty
      else
        EnergyGAM=Evalue
      endif
    endif

    Glob_EnergyGACounter=Glob_EnergyGACounter+1

  end function EnergyGAM

  subroutine EnergyGB(Evalue,Gradient,AreMatElemNeeded,ErrorCode)
!subroutine EnergyGB computes the energy value and the
!gradient with respect to the nonlinear parameters of the last
!Glob_nfo functions of the system under consideration. Routine
!DSYGVX from LAPACK is used to solve GSEP. The subroutine finds
!the energy level defined by the  global variable
!Glob_WhichEigenvalue, as well as the corresponding gradient. It
!is assumed that nonlinear parameters of the basis functions are
!stored in global array Glob_NonlinParam and all the matrix
!elements of the first Glob_nfru functions have been computed
!and stored in proper global arrays. If no matrix elements have
!been calculated, one should first call ComputeMatElemAndDeriv
!and calculate the matrix elements of the first Glob_nfru
!functions. In case when all matrix elements have been computed
!and it is only necessary to solve GSEP and find gradient, one
!needs to set AreMatElemNeeded=.false.
!Upon successful exit ErrorCode should be equal to 0. If
!ErrorCode <= Glob_nfa then the eigenvector failed to
!converge. If ErrorCode = Glob_nfa + i then the
!leading minor of S of size i is not positive definite.
!The data in Gradient is ordered in the following manner:
!Gradient=(dEdvechLi,dEdvechL{i+1},...,dEdvechLj),
!where i=Glob_nfru+1, j=Glob_nfa

!Arguments:
    real(wp)    Evalue
    real(wp)    Gradient(Glob_npt*Glob_nfo)
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables
    integer        nfo,nfa,nfru,npt
    integer        i,j,k,l,m,nbands,leftover
    logical        oddband
    integer        N,NumOfEigvalsFound,IFAIL(1)
    real(wp)    EVs(1)
    real(wp)    W(Glob_npt_MaxAllowed),t,t2
    real(wp)    pen_coeff

    nfo=Glob_nfo
    nfa=Glob_nfa
    npt=Glob_npt
    nfru=Glob_nfru

    if (AreMatElemNeeded) call ComputeMatElemAndDeriv(nfru+1,nfa)

    if (nfa==1) then
      Evalue=Glob_diagH(1)/Glob_diagS(1)
      Glob_c(1)=ONE
      ErrorCode=0
    else
      !Copying H and S matrix elements from the lower triangles
      !to the upper ones. The diagonals are copied from the global
      !arrays where they are stored
      do i=1,nfa
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,nfa
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo

      if (Glob_ProcID==0) then
        call   DSYGVX(1,'V','I','U',nfa,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,   &
                      ZERO,ZERO,Glob_WhichEigenvalue,Glob_WhichEigenvalue,Glob_AbsTolForDSYGVX, &
                      NumOfEigvalsFound,EVs,Glob_c,nfa,Glob_WorkForDSYGVX,  &
                      Glob_LWorkForDSYGVX,Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
        ! SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
!       LWORK, IWORK, IFAIL, INFO )
        Evalue=EVs(1)
      endif
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_c,nfa,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

!Computing gradient
    do k=1,nfo
      W(1:npt)=ZERO
      do l=1+Glob_ProcID,nfa,Glob_NumOfProcs
        t=Glob_c(l)
        do m=1,npt
          W(m)=W(m)+t*(Glob_D(m,k,l)-Evalue*Glob_D(m+npt,k,l))
        enddo
      enddo
      t=Glob_c(k+nfru)
      t2=t*t
      do m=1,npt
        Glob_WkGR((k-1)*npt+m)=TWO*t*W(m)
      enddo
      do m=1+Glob_ProcID,npt,Glob_NumOfProcs
        Glob_WkGR((k-1)*npt+m)=Glob_WkGR((k-1)*npt+m)-t2*(Glob_D(m,k,k+nfru) &
                                                          -Evalue*Glob_D(m+npt,k,k+nfru))
      enddo
    enddo

    if ((Glob_OverlapPenaltyAllowed).and.(nfa/=1)) then
      call ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2, &
                                               Glob_TotalOverlapPenalty,Glob_WkGR)
      Evalue=Evalue+Glob_TotalOverlapPenalty
    endif
    call MPI_ALLREDUCE(Glob_WkGR,Gradient,nfo*npt,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

    Glob_EnergyGBCounter=Glob_EnergyGBCounter+1

  end subroutine EnergyGB

  function EnergyIA(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyIA computes the energy of the system under
!consideration. Routine GSEPIIS from module linalg is used to
!solve GSEP. The function computes the energy level which is
!the closest one to the value of the global variable
!Glob_ApproxEnergy. The calculations are done with a basis
!of Nmax functions. It is assumed that nonlinear parameters
!of the basis functions are stored in global array
!Glob_NonlinParam and all matrix elements of the
!first Nmin-1 functions have already been calculated and stored in
!proper global arrays. If no matrix elements have been calculated,
!one should set Nmin=0. In the case when all matrix elements have
!been computed and only the solution of GSEP is needed, one
!should set AreMatElemNeeded=.false.
!Upon successful exit ErrorCode should be equal to 0. In the
!case of nonzero ErrorCode please see the description of GSEPIIS
!as ErrorCode is simply passed from that routine to EnergyIA

    real(wp)    EnergyIA
!Arguments:
    integer        Nmin,Nmax
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables:
    real(wp)    Evalue
    integer        NumOfIterations

    if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
    if (Nmax==1) then
      EnergyIA=Glob_diagH(1)
      NumOfIterations=1
      ErrorCode=0
    else
      Glob_c(1:Nmax)=Glob_LastEigvector(1:Nmax)
      call GSEPIIS(Nmin,Nmax,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim, &
                   Glob_ApproxEnergy,Glob_c,Glob_WorkForGSEPIIS,Glob_EigvalTol, &
                   Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGSEPIIS, &
                   -1,NumOfIterations,ErrorCode)
      !GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
      if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol
      if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol
      if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_OverlapPenaltyAllowed) then
        EnergyIA=Evalue+Glob_TotalOverlapPenalty
      else
        EnergyIA=Evalue
      endif
    endif

    Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIACounter=Glob_EnergyIACounter+1

  end function EnergyIA

  function EnergyIAM(Nmin,Nmax,AreMatElemNeeded,ErrorCode)
!Function EnergyIAM is pretty much the same function as
!EnergyIAM. The only exception is that it also generates
!a properly normalized eigenvector, which gives correct
!linear coefficients (stored in Glob_c).

    real(wp)    EnergyIAM
!Arguments:
    integer        Nmin,Nmax
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables:
    real(wp)    Evalue
    integer        NumOfIterations

    if (AreMatElemNeeded) call ComputeMatElem(Nmin,Nmax)
    if (Nmax==1) then
      EnergyIAM=Glob_H(1,1)+Glob_ApproxEnergy
      Glob_c(1)=ONE
      NumOfIterations=1
      ErrorCode=0
    else
      Glob_c(1:Nmax)=Glob_LastEigvector(1:Nmax)
      call GSEPIIS(Nmin,Nmax,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim, &
                   Glob_ApproxEnergy,Glob_c,Glob_WorkForGSEPIIS,Glob_EigvalTol, &
                   Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGSEPIIS, &
                   0,NumOfIterations,ErrorCode)
      !GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
      Glob_c(1:Nmax)=Glob_LastEigvector(1:Nmax)
      if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol
      if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol
      if (Glob_OverlapPenaltyAllowed) call ComputeOverlapPenalty(Glob_MaxOverlapPenalty, &
                                                                 Glob_OverlapPenaltyThreshold2,Glob_TotalOverlapPenalty)
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_OverlapPenaltyAllowed) then
        EnergyIAM=Evalue+Glob_TotalOverlapPenalty
      else
        EnergyIAM=Evalue
      endif
    endif

    Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIACounter=Glob_EnergyIACounter+1

  end function EnergyIAM

  subroutine EnergyIB(Evalue,Gradient,AreMatElemNeeded,ErrorCode)
!subroutine EnergyIB computes the energy value and the
!gradient with respect to the nonlinear parameters of the last
!Glob_nfo functions of the system under consideration.
!Routine GSEPIIS from module linalg is used to
!solve GSEP. The subroutine computes the energy level which is
!the closest one to the value of the global variable
!Glob_ApproxEnergy, as well as the corresponding gradient.
!It is assumed that nonlinear parameters of the basis functions are
!stored in global array Glob_NonlinParam and all the matrix
!elements of the first Glob_nfru functions have been computed
!and stored in proper global arrays. If no matrix elements have
!been calculated, one should first call ComputeMatElemAndDeriv
!and calculate the matrix elements of the first Glob_nfru
!functions. In case when all matrix elements have been computed
!and it is only necessary to solve GSEP and find gradient, one
!needs to set AreMatElemNeeded=.false.
!Upon successful exit ErrorCode should be equal to 0. In the
!case of nonzero ErrorCode please see the description of GSEPIIS
!as ErrorCode is simply passed from that routine to EnergyIA.
!The data in Gradient is ordered in the following manner:
!Gradient=(dEdvechLi,dEdvechL{i+1},...,dEdvechLj),
!where i=Glob_nfru+1, j=Glob_nfa

!Arguments:
    real(wp)    Evalue
    real(wp)    Gradient(Glob_npt*Glob_nfo)
    logical        AreMatElemNeeded
    integer        ErrorCode

!Local variables:
    integer        nfo,nfa,nfru,npt
    integer        i,j,k,l,m,nbands,leftover
    logical        oddband
    integer        NumOfIterations
    real(wp)    W(Glob_npt_MaxAllowed),t,t2
    real(wp)    pen_coeff

    nfo=Glob_nfo
    nfa=Glob_nfa
    npt=Glob_npt
    nfru=Glob_nfru

    if (AreMatElemNeeded) call ComputeMatElemAndDeriv(nfru+1,nfa)

    if (nfa==1) then
      Evalue=Glob_H(1,1)+Glob_ApproxEnergy
      Glob_c(1)=ONE
      NumOfIterations=1
      ErrorCode=0
    else
      Glob_c(1:nfa)=Glob_LastEigvector(1:nfa)
      call GSEPIIS(nfru+1,nfa,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim, &
                   Glob_ApproxEnergy,Glob_c,Glob_WorkForGSEPIIS,Glob_EigvalTol, &
                   Evalue,Glob_LastEigvector,Glob_LastEigvalTol,Glob_MaxIterForGSEPIIS, &
                   0,NumOfIterations,ErrorCode)
      !GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
      !lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
      Glob_c(1:nfa)=Glob_LastEigvector(1:nfa)
      if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol
      if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

!Computing gradient
    do k=1,nfo
      W(1:npt)=ZERO
      do l=1+Glob_ProcID,nfa,Glob_NumOfProcs
        t=Glob_c(l)
        do m=1,npt
          W(m)=W(m)+t*(Glob_D(m,k,l)-Evalue*Glob_D(m+npt,k,l))
        enddo
      enddo
      t=Glob_c(k+nfru)
      t2=t*t
      do m=1,npt
        Glob_WkGR((k-1)*npt+m)=2*t*W(m)
      enddo
      do m=1+Glob_ProcID,npt,Glob_NumOfProcs
        Glob_WkGR((k-1)*npt+m)=Glob_WkGR((k-1)*npt+m)-t2*(Glob_D(m,k,k+nfru) &
                                                          -Evalue*Glob_D(m+npt,k,k+nfru))
      enddo
    enddo

    if ((Glob_OverlapPenaltyAllowed).and.(nfa/=1)) then
      call ComputeOverlapPenaltyAndAddGradient(Glob_MaxOverlapPenalty,Glob_OverlapPenaltyThreshold2, &
                                               Glob_TotalOverlapPenalty,Glob_WkGR)
      Evalue=Evalue+Glob_TotalOverlapPenalty
    endif
    call MPI_ALLREDUCE(Glob_WkGR,Gradient,nfo*npt,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)

    Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
    Glob_EnergyIBCounter=Glob_EnergyIBCounter+1

  end subroutine EnergyIB

  subroutine ReadSwapFileAndDistributeData(IsSwapFileOK)
!This subroutine reads data from a swap file (if it is
!available, and if global variable Glob_UseSwapFile=.true.)
!and sends the read data to all processes. If the result is
!success then the value of a logical variable IsSwapFileOK
!is .true. on exit. It is assumed that the current basis
!size is equal to Glob_CurrBasisSize.

!Arguments:
    logical IsSwapFileOK
!Local variables:
    integer i,j,OpenFileErr

    IsSwapFileOK=.false.
    if (Glob_UseSwapFile) then
      if (Glob_ProcID==0) then
        open(1,file=Glob_SwapFileName,form='unformatted',status='old',iostat=OpenFileErr)
        if (OpenFileErr==0) then
          read(1) i
          if (i==Glob_CurrBasisSize) then
            !Reading H and S from a matrix stored in file
            !Remember that the lower part (including the diagonal)
            !contains elements of H, while the upper part contains S
            write(*,'(1x,a33)',advance='no') 'Reading H and S from swap file...'
            do j=1,Glob_CurrBasisSize
              if (OpenFileErr==0) then
                read(1,iostat=OpenFileErr) Glob_H(1:Glob_CurrBasisSize,j)
              endif
            enddo
            !reading diagS
            if (OpenFileErr==0) read(1,iostat=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize)
            if (OpenFileErr==0) then
              IsSwapFileOK=.true.
              write(*,*) 'completed'
              close(1)
            else
              write(*,*) 'failed'
            endif
            !erase information in swap file to free disc space
            open(1,file=Glob_SwapFileName,form='unformatted',status='replace',iostat=OpenFileErr)
            write(1) 'Swap file is empty'
            close(1)
          endif
        endif
      endif
    endif
    call MPI_BCAST(IsSwapFileOK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!If swap file is OK then send the data to all processes
    if (IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a35)',advance='no') 'Sending H and S to all processes...'
      call MPI_BCAST(Glob_H,Glob_HSLeadDim*Glob_HSLeadDim,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_diagS,Glob_CurrBasisSize,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      !Remember that the lower part (including the diagonal)
      !contains elements of H, while the upper part contains S

      !Restoring proper storage of data for Glob_GSEPSolutionMethod='G'
      if (Glob_GSEPSolutionMethod=='G') then
        do i=1,Glob_CurrBasisSize
          do j=1,i-1
            Glob_S(i,j)=Glob_H(j,i)
          enddo
          Glob_diagH(i)=Glob_H(i,i)
        enddo
      endif
      !Restoring proper storage of data for Glob_GSEPSolutionMethod='I'
      if (Glob_GSEPSolutionMethod=='I') then
        do i=1,Glob_CurrBasisSize
          do j=1,i-1
            Glob_S(i,j)=Glob_H(j,i)
            Glob_S(j,i)=Glob_H(j,i)
          enddo
          Glob_S(i,i)=ONE
        enddo
        do i=1,Glob_CurrBasisSize
          Glob_H(i,i)=Glob_H(i,i)-Glob_ApproxEnergy
          do j=i+1,Glob_CurrBasisSize
            Glob_H(j,i)=Glob_H(j,i)-Glob_ApproxEnergy*Glob_S(j,i)
          enddo
        enddo
      endif

      !Q uses the same compact swap-file representation as G, but restores
      !unshifted physical matrices with both diagonals stored in place. The
      !lower triangles remain the only authoritative matrix storage.
      if (Glob_GSEPSolutionMethod=='Q') then
        do i=1,Glob_CurrBasisSize
          do j=1,i-1
            Glob_S(i,j)=Glob_H(j,i)
          enddo
          Glob_S(i,i)=ONE
        enddo
      endif

      if (Glob_ProcID==0) write(*,*) 'completed'
    else
      if (Glob_ProcID==0) then
        write(*,*) 'Matrices H and S were not read from swap file'
        write(*,*) 'All H and S matrix elements need be (re)computed'
      endif
    endif

  end subroutine ReadSwapFileAndDistributeData

  subroutine StoreMatricesInSwapFile()
!This subroutine stores H and S matrices, including
!their diagonals, in a swap file. This is done only
!if the value of global variable Glob_UseSwapFile is .true.

!Local variables:
    integer i,j,OpenFileErr

    if (Glob_UseSwapFile) then
!If swap file is allowed to use then write H and S matrix
!elements into it
      if (Glob_CurrBBOPStep/=Glob_NumOfBBOPSteps) then
        if (Glob_ProcID==0) then
          open(1,file=Glob_SwapFileName,form='unformatted',status='replace',iostat=OpenFileErr)
          if (OpenFileErr==0) then
            write(1) Glob_CurrBasisSize
            !We store a matrix whose lower part (including the diagonal)
            !contains elements of H, while the upper part contains S.
            !Also, we store the diagonal of S.

            if (Glob_GSEPSolutionMethod=='G') then
              do i=1,Glob_CurrBasisSize
                do j=1,i-1
                  Glob_H(j,i)=Glob_S(i,j)
                enddo
                Glob_H(i,i)=Glob_diagH(i)
              enddo
            endif
            if (Glob_GSEPSolutionMethod=='I') then
              do i=1,Glob_CurrBasisSize
                do j=1,i-1
                  Glob_H(j,i)=Glob_S(i,j)
                enddo
                Glob_H(i,i)=Glob_H(i,i)+Glob_ApproxEnergy
                do j=i+1,Glob_CurrBasisSize
                  Glob_H(j,i)=Glob_H(j,i)+Glob_ApproxEnergy*Glob_S(j,i)
                enddo
              enddo
            endif

            if (Glob_GSEPSolutionMethod=='Q') then
              !Borrow only the unused upper triangle of H for serialized S.
              !The canonical lower triangles and the physical H diagonal are
              !unchanged, so live factors still match after this routine.
              do i=1,Glob_CurrBasisSize
                do j=1,i-1
                  Glob_H(j,i)=Glob_S(i,j)
                enddo
              enddo
            endif

            write(*,*)
            write(*,'(1x,a33)',advance='no') 'Writing H and S into swap file... '
            do i=1,Glob_CurrBasisSize
              if (OpenFileErr==0) then
                write(1,iostat=OpenFileErr) Glob_H(1:Glob_CurrBasisSize,i)
              endif
            enddo
            if (OpenFileErr==0) then
              write(1,iostat=OpenFileErr) Glob_diagS(1:Glob_CurrBasisSize)
            endif
            if (OpenFileErr==0) then
              write(*,*) 'completed'
              write(*,*)
            else
              write(*,*) 'failed to complete'
              write(*,*)
            endif
          endif
          close(1)
        endif
      endif
    endif

  end subroutine StoreMatricesInSwapFile

  subroutine ReadHessianFile(V,IVLMAT,D,nvar,FileName,IsHessFileOK)
!This subroutine reads Hessian from file FileName (if it exists
!and contains appropriate amount of data). nvar stands for the
!total number of variables in the problem. The Hessian is
!placed in array V, starting from element V(IVLMAT).
!In the case when the global parameter Glob_FullOptSaveD=.true.
!the subroutine also reads scaling vector D from the same
!file. It is assumed that this scaling vector is placed after
!Hessian data in the file.
!If the reading is successful then the value of logical variable
!IsHessFileOK is .true. on exit.

!Arguments:
    real(wp) V(*),D(*)
    integer     IVLMAT,nvar
    character(Glob_FileNameLength) FileName
    logical     IsHessFileOK
!Local variables:
    integer i,j,OpenFileErr

    IsHessFileOK=.false.
    if (Glob_ProcID==0) then
      open(1,file=FileName,form='unformatted',status='old',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        read(1,iostat=OpenFileErr) j
        if ((j==nvar).and.(OpenFileErr==0)) then
          write(*,'(1x,a28)',advance='no') 'Reading Hessian from file...'
          read(1,iostat=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)
          if ((OpenFileErr==0).and.(Glob_FullOptSaveD)) read(1,iostat=OpenFileErr) D(1:nvar)
          if (OpenFileErr==0) then
            IsHessFileOK=.true.
            write(*,*) 'done'
            close(1)
          else
            write(*,*) 'failed'
            write(*,*) 'Warning WC0100: Default Hessian initialization must be used'
          endif
        else
          write(*,*) 'Warning WC0101: Hessian file ',FileName
          write(*,*) 'is inconsistent with the dimension of the current optimization problem'
          write(*,*) 'and will not be read'
          write(*,*) 'Default Hessian initialization must be used'
        endif
      else
        write(*,*) 'Warning WC0102: Cannot open Hessian file ',FileName
        write(*,*) 'Default Hessian initialization must be used'
      endif
    endif

  end subroutine ReadHessianFile

  subroutine SaveHessianFile(V,IVLMAT,D,nvar,FileName,IsSuccess)
!This subroutine saves Hessian to file FileName. nvar stands for the
!total number of variables in the problem. The Hessian is taken from
!V, starting from element V(IVLMAT).
!In case when the global parameter Glob_FullOptSaveD=.true. the
!subroutine also saves the scaling vector D (used by the optimization
!routine). Vector D is placed after V, i.e. in the end of the file, so
!that even if user happens to continue the same optimization problem with
!Glob_FullOptSaveD changed from .true. to .false., the routine
!will be able to run.
!If the file is successfully saved then the value of logical variable
!IsSuccess is .true. on exit.

!Arguments:
    real(wp) V(*),D(*)
    integer     IVLMAT,nvar
    character(Glob_FileNameLength) FileName
    logical     IsSuccess
!Local variables:
    integer i,OpenFileErr

    IsSuccess=.false.
    if (Glob_ProcID==0) then
      open(1,file=FileName,form='unformatted',status='replace',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        write(1,iostat=OpenFileErr) nvar
        if (OpenFileErr==0) then
          write(*,'(1x,a17)',advance='no') 'Saving Hessian...'
          write(1,iostat=OpenFileErr) V(IVLMAT:IVLMAT+nvar*(nvar+1)/2-1)
          if (Glob_FullOptSaveD) write(1,iostat=OpenFileErr) D(1:nvar)
          close(1)
          if (OpenFileErr==0) then
            IsSuccess=.true.
            write(*,*) 'done'
          else
            write(*,*) 'failed'
          endif
        endif
      endif
      if (.not.IsSuccess) write(*,*) 'Warning WC0105: Hessian was not saved'
    endif

  end subroutine SaveHessianFile

  subroutine PermuteFunctions(fb,fe,FuncNumTemp,NonlinParamTemp)
!Subroutine PermuteFunctions permutes basis functions in such
!a way that functions fb through fe go to the very end of the
!basis (which makes their optimization more time efficient).
!The subroutine permutes the nonlinear parameters (Glob_NonlinParam)
!x,y,z-indices (Glob_Index), and function numbers (Glob_FuncNum).
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp and NonlinParamTemp. The length of FuncNumTemp
!should be at least fe-fb+1, while the length of NonlinParamTemp must be
!at least (fe-fb+1)*Glob_npt

!Arguments:
    integer      fb,fe
    integer      FuncNumTemp(*)
    real(wp)  NonlinParamTemp(Glob_npt,*)
!Local variables:
    integer      i,fbn,nfco

    nfco=fe-fb+1
    fbn=fb+Glob_CurrBasisSize-fe

    NonlinParamTemp(1:Glob_npt,1:nfco)=Glob_NonlinParam(1:Glob_npt,fb:fe)
    do i=fb,fbn-1
      Glob_NonlinParam(1:Glob_npt,i)=Glob_NonlinParam(1:Glob_npt,nfco+i)
    enddo
    Glob_NonlinParam(1:Glob_npt,fbn:Glob_CurrBasisSize)=NonlinParamTemp(1:Glob_npt,1:nfco)

    FuncNumTemp(1:nfco)=Glob_Index(fb:fe,1)
    do i=fb,fbn-1
      Glob_Index(i,1)=Glob_Index(nfco+i,1)
    enddo
    Glob_Index(fbn:Glob_CurrBasisSize,1)=FuncNumTemp(1:nfco)

    FuncNumTemp(1:nfco)=Glob_Index(fb:fe,2)
    do i=fb,fbn-1
      Glob_Index(i,2)=Glob_Index(nfco+i,2)
    enddo
    Glob_Index(fbn:Glob_CurrBasisSize,2)=FuncNumTemp(1:nfco)

    FuncNumTemp(1:nfco)=Glob_FuncNum(fb:fe)
    do i=fb,fbn-1
      Glob_FuncNum(i)=Glob_FuncNum(nfco+i)
    enddo
    Glob_FuncNum(fbn:Glob_CurrBasisSize)=FuncNumTemp(1:nfco)

  end subroutine PermuteFunctions

  subroutine PermuteFunctions2(fb1,fe1,fe2,FuncNumTemp,NonlinParamTemp)
!Subroutine PermuteFunctions2 permutes basis functions in such
!a way that the set of functions fb1 through fe1 exchanges its positions
!with the set of function fb2 through fe2, where fb2=fe1+1.
!The subroutine permutes the nonlinear parameters (Glob_NonlinParam),
!x,y,z-indices (Glob_Index), and function numbers (Glob_FuncNum).
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp and NonlinParamTemp. The length of FuncNumTemp
!should be at least fe1-fb1+1, while the length of NonlinParamTemp must be
!at least (fe1-fb1+1)*Glob_npt

!Arguments:
    integer      fb1,fe1,fe2
    integer      FuncNumTemp(*)
    real(wp)      NonlinParamTemp(Glob_npt,*)
!Local variables:
    integer  i,fbn,k

    k=fe1-fb1+1
    fbn=fb1+fe2-fe1

    NonlinParamTemp(1:Glob_npt,1:k)=Glob_NonlinParam(1:Glob_npt,fb1:fe1)
    do i=fb1,fbn-1
      Glob_NonlinParam(1:Glob_npt,i)=Glob_NonlinParam(1:Glob_npt,k+i)
    enddo
    Glob_NonlinParam(1:Glob_npt,fbn:fe2)=NonlinParamTemp(1:Glob_npt,1:k)

    FuncNumTemp(1:k)=Glob_Index(fb1:fe1,1)
    do i=fb1,fbn-1
      Glob_Index(i,1)=Glob_Index(k+i,1)
    enddo
    Glob_Index(fbn:fe2,1)=FuncNumTemp(1:k)

    FuncNumTemp(1:k)=Glob_Index(fb1:fe1,2)
    do i=fb1,fbn-1
      Glob_Index(i,2)=Glob_Index(k+i,2)
    enddo
    Glob_Index(fbn:fe2,2)=FuncNumTemp(1:k)

    FuncNumTemp(1:k)=Glob_FuncNum(fb1:fe1)
    do i=fb1,fbn-1
      Glob_FuncNum(i)=Glob_FuncNum(k+i)
    enddo
    Glob_FuncNum(fbn:fe2)=FuncNumTemp(1:k)

  end subroutine PermuteFunctions2

  subroutine PermuteMatrixElements(fb,fe,TempR)
!This subroutine permutes matrix elements of H and S in such a
!way that functions fb and fe go to the very end of the basis
!(which makes their optimization more time efficient). The subroutine
!permutes elements of arrays Glob_H, Glob_S, Glob_diagH, Glob_diagS
!This subroutine requires workspace, which must be provided by
!by array TempR. The length of this array should be at
!least fe-fb+1.
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper
!triangles (excluding the diagonal) are used as workspace.

!Arguments:
    integer        fb,fe
    real(wp)    TempR(*)
!Local variables:
    integer        i,j,fbn,nfco,fep,k,q,p,r,s

    nfco=fe-fb+1
    fbn=fb+Glob_CurrBasisSize-fe

!First we do the diagonals of H and S
    if (Glob_GSEPSolutionMethod=='G') then
      TempR(1:nfco)=Glob_diagH(fb:fe)
      do i=fb,fbn-1
        Glob_diagH(i)=Glob_diagH(nfco+i)
      enddo
      Glob_diagH(fbn:Glob_CurrBasisSize)=TempR(1:nfco)
    else
      do i=fb,fe
        TempR(i-fb+1)=Glob_H(i,i)
      enddo
      do i=fb,fbn-1
        Glob_H(i,i)=Glob_H(nfco+i,nfco+i)
      enddo
      do i=fbn,Glob_CurrBasisSize
        Glob_H(i,i)=TempR(i-fbn+1)
      enddo
    endif
    TempR(1:nfco)=Glob_diagS(fb:fe)
    do i=fb,fbn-1
      Glob_diagS(i)=Glob_diagS(nfco+i)
    enddo
    Glob_diagS(fbn:Glob_CurrBasisSize)=TempR(1:nfco)

!Now we permute off-diagonal elements
    fep=fe+1
    k=Glob_CurrBasisSize-fe
    q=Glob_CurrBasisSize+fep
    p=fb+Glob_CurrBasisSize-fe-1
    s=Glob_CurrBasisSize+fb
    r=s-1
    do i=1,fb-1
      TempR(1:nfco)=Glob_H(fb:fe,i)
      do j=fb,fbn-1
        Glob_H(j,i)=Glob_H(j+nfco,i)
      enddo
      Glob_H(fbn:Glob_CurrBasisSize,i)=TempR(1:nfco)
    enddo
    do i=fep,Glob_CurrBasisSize-1
      Glob_H(fep:q-i-1,q-i)=Glob_H(i+1:Glob_CurrBasisSize,i)
    enddo
    do i=fe+1,Glob_CurrBasisSize
      do j=fb,fe
        Glob_H(j,i)=Glob_H(i,j)
      enddo
    enddo
    do i=fb,fe-1
      Glob_H(i+k+1:Glob_CurrBasisSize,i+k)=Glob_H(i+1:fe,i)
    enddo
    Glob_H(fbn:Glob_CurrBasisSize,fb:fb+k-1)=Glob_H(fb:fe,fe+1:Glob_CurrBasisSize)
    do i=fb,p-1
      Glob_H(i+1:p,i)=Glob_H(fep:r-i,s-i)
    enddo
    do i=1,fb-1
      TempR(1:nfco)=Glob_S(fb:fe,i)
      do j=fb,fbn-1
        Glob_S(j,i)=Glob_S(j+nfco,i)
      enddo
      Glob_S(fbn:Glob_CurrBasisSize,i)=TempR(1:nfco)
    enddo
    do i=fep,Glob_CurrBasisSize-1
      Glob_S(fep:q-i-1,q-i)=Glob_S(i+1:Glob_CurrBasisSize,i)
    enddo
    do i=fe+1,Glob_CurrBasisSize
      do j=fb,fe
        Glob_S(j,i)=Glob_S(i,j)
      enddo
    enddo
    do i=fb,fe-1
      Glob_S(i+k+1:Glob_CurrBasisSize,i+k)=Glob_S(i+1:fe,i)
    enddo
    Glob_S(fbn:Glob_CurrBasisSize,fb:fb+k-1)=Glob_S(fb:fe,fe+1:Glob_CurrBasisSize)
    do i=fb,p-1
      Glob_S(i+1:p,i)=Glob_S(fep:r-i,s-i)
    enddo
    if (Glob_GSEPSolutionMethod=='I') then
      !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      !upper triangle of Glob_S
      do i=fb,fe
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
      enddo
    endif

  end subroutine PermuteMatrixElements

  subroutine PermuteMatrixElements2(fb1,fe1,fe2,TempR)
!This subroutine permutes matrix elements of H and S in such a
!way that the set of functions fb1 through fe1 exchanges its
!positions with the set of function fb2 through fe2, where
!fb2=fe1+1.
!The subroutine permutes elements of arrays Glob_H, Glob_S,
!Glob_diagH, Glob_diagS. This subroutine requires workspace,
!which must be provided by array TempR. The length of
!this arrays should be at least fe1-fb1+1.
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper
!triangles (excluding the diagonal) are used as workspace.
!Arguments:
    integer      fb1,fe1,fe2
    real(wp)  TempR(*)
!Local variables:
    integer      i,j,fn,t,fe1p,k,q,p,r,s,fb2

    t=fe1-fb1+1
    fn=fb1+fe2-fe1
    fb2=fe1+1

!First we do the diagonals of H and S
    if (Glob_GSEPSolutionMethod=='G') then
      TempR(1:t)=Glob_diagH(fb1:fe1)
      do i=fb1,fn-1
        Glob_diagH(i)=Glob_diagH(t+i)
      enddo
      Glob_diagH(fn:fe2)=TempR(1:t)
    else
      do i=fb1,fe1
        TempR(i-fb1+1)=Glob_H(i,i)
      enddo
      do i=fb1,fn-1
        Glob_H(i,i)=Glob_H(t+i,t+i)
      enddo
      do i=fn,fe2
        Glob_H(i,i)=TempR(i-fn+1)
      enddo
    endif
    TempR(1:t)=Glob_diagS(fb1:fe1)
    do i=fb1,fn-1
      Glob_diagS(i)=Glob_diagS(t+i)
    enddo
    Glob_diagS(fn:fe2)=TempR(1:t)

!Now we permute off-diagonal elements
    fe1p=fe1+1
    k=fe2-fe1
    q=fe2+fe1p
    p=fb1+fe2-fe1-1
    s=fe2+fb1
    r=s-1
    do i=1,fb1-1
      TempR(1:t)=Glob_H(fb1:fe1,i)
      do j=fb1,fn-1
        Glob_H(j,i)=Glob_H(t+j,i)
      enddo
      Glob_H(fn:fe2,i)=TempR(1:t)
    enddo
    do i=fe1p,fe2-1
      Glob_H(fe1p:q-i-1,q-i)=Glob_H(i+1:fe2,i)
    enddo
    do i=fe1+1,fe2
      do j=fb1,fe1
        Glob_H(j,i)=Glob_H(i,j)
      enddo
    enddo
    do i=fb1,fe1-1
      Glob_H(i+k+1:fe2,i+k)=Glob_H(i+1:fe1,i)
    enddo
    Glob_H(fn:fe2,fb1:fb1+k-1)=Glob_H(fb1:fe1,fe1+1:fe2)
    do i=fb1,p-1
      Glob_H(i+1:p,i)=Glob_H(fe1p:r-i,s-i)
    enddo
    do i=fe2+1,Glob_CurrBasisSize
      TempR(1:t)=Glob_H(i,fb1:fe1)
      do j=fb1,fn-1
        Glob_H(i,j)=Glob_H(i,t+j)
      enddo
      Glob_H(i,fn:fe2)=TempR(1:t)
    enddo
    do i=1,fb1-1
      TempR(1:t)=Glob_S(fb1:fe1,i)
      do j=fb1,fn-1
        Glob_S(j,i)=Glob_S(t+j,i)
      enddo
      Glob_S(fn:fe2,i)=TempR(1:t)
    enddo
    do i=fe1p,fe2-1
      Glob_S(fe1p:q-i-1,q-i)=Glob_S(i+1:fe2,i)
    enddo
    do i=fe1+1,fe2
      do j=fb1,fe1
        Glob_S(j,i)=Glob_S(i,j)
      enddo
    enddo
    do i=fb1,fe1-1
      Glob_S(i+k+1:fe2,i+k)=Glob_S(i+1:fe1,i)
    enddo
    Glob_S(fn:fe2,fb1:fb1+k-1)=Glob_S(fb1:fe1,fe1+1:fe2)
    do i=fb1,p-1
      Glob_S(i+1:p,i)=Glob_S(fe1p:r-i,s-i)
    enddo
    do i=fe2+1,Glob_CurrBasisSize
      TempR(1:t)=Glob_S(i,fb1:fe1)
      do j=fb1,fn-1
        Glob_S(i,j)=Glob_S(i,t+j)
      enddo
      Glob_S(i,fn:fe2)=TempR(1:t)
    enddo
    if (Glob_GSEPSolutionMethod=='I') then
      !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      !upper triangle of Glob_S
      do i=fb1,fe2
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
      enddo
    endif

  end subroutine PermuteMatrixElements2

  subroutine ReverseFuncOrder(fb,fe)
!Subroutine ReverseFuncOrder changes the order of basis functions
!fb through fe to reverse.
!The subroutine permutes nonlinear parameters (Glob_NonlinParam),
!x,y,z-indices (Glob_Index), and function numbers (Glob_FuncNum).

!Arguments:
    integer  fb,fe
!Local variables:
    real(wp) temp(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1))
    integer i,j,f,t,fbm,fep

    f=(fe-fb+1)/2 !integer division!
    fbm=fb-1
    fep=fe+1
    do i=1,f
      temp(1:Glob_npt)=Glob_NonlinParam(1:Glob_npt,fbm+i)
      Glob_NonlinParam(1:Glob_npt,fbm+i)=Glob_NonlinParam(1:Glob_npt,fep-i)
      Glob_NonlinParam(1:Glob_npt,fep-i)=temp(1:Glob_npt)
    enddo

    do i=1,f
      t=Glob_Index(fbm+i,1)
      Glob_Index(fbm+i,1)=Glob_Index(fep-i,1)
      Glob_Index(fep-i,1)=t
    enddo

    do i=1,f
      t=Glob_Index(fbm+i,2)
      Glob_Index(fbm+i,2)=Glob_Index(fep-i,2)
      Glob_Index(fep-i,2)=t
    enddo

    do i=1,f
      t=Glob_FuncNum(fbm+i)
      Glob_FuncNum(fbm+i)=Glob_FuncNum(fep-i)
      Glob_FuncNum(fep-i)=t
    enddo

  end subroutine ReverseFuncOrder

  subroutine ReverseMatElemOrder(fb,fe)
!Subroutine ReverseMatElemOrder changes the order of matrix elements
!corresponding to basis functions fb through fe to reverse.
!The subroutine permutes elements of arrays Glob_H, Glob_S, Glob_diagH,
!Glob_diagS. Note that it uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonal). Upper triangles are not
!referenced.

!Arguments:
    integer  fb,fe
!Local variables:
    integer i,j,f,fbm,fep,ff
    real(wp)     r
    complex(wp)  c

    f=(fe-fb+1)/2 !integer division!
    fbm=fb-1
    fep=fe+1
    ff=fe+fb
!Diagonal elements
    if (Glob_GSEPSolutionMethod=='G') then
      do i=1,f
        r=Glob_diagH(fbm+i)
        Glob_diagH(fbm+i)=Glob_diagH(fep-i)
        Glob_diagH(fep-i)=r
      enddo
    else
      do i=1,f
        r=Glob_H(fbm+i,fbm+i)
        Glob_H(fbm+i,fbm+i)=Glob_H(fep-i,fep-i)
        Glob_H(fep-i,fep-i)=r
      enddo
    endif
    do i=1,f
      r=Glob_diagS(fbm+i)
      Glob_diagS(fbm+i)=Glob_diagS(fep-i)
      Glob_diagS(fep-i)=r
    enddo
!Off-diagonal elements
    do i=1,fbm
      do j=1,f
        c=Glob_H(fbm+j,i)
        Glob_H(fbm+j,i)=Glob_H(fep-j,i)
        Glob_H(fep-j,i)=c
      enddo
    enddo
    do i=1,f
      do j=fb+i,fe-i
        c=Glob_H(j,fbm+i)
        Glob_H(j,fbm+i)=Glob_H(fep-i,ff-j)
        Glob_H(fep-i,ff-j)=c
      enddo
      Glob_H(j,fbm+i)=Glob_H(j,fbm+i)
    enddo
    do i=fb,fb+f-1
      do j=fep,Glob_CurrBasisSize
        c=Glob_H(j,i)
        Glob_H(j,i)=Glob_H(j,ff-i)
        Glob_H(j,ff-i)=c
      enddo
    enddo
    do i=1,fbm
      do j=1,f
        c=Glob_S(fbm+j,i)
        Glob_S(fbm+j,i)=Glob_S(fep-j,i)
        Glob_S(fep-j,i)=c
      enddo
    enddo
    do i=1,f
      do j=fb+i,fe-i
        c=Glob_S(j,fbm+i)
        Glob_S(j,fbm+i)=Glob_S(fep-i,ff-j)
        Glob_S(fep-i,ff-j)=c
      enddo
      Glob_S(j,fbm+i)=Glob_S(j,fbm+i)
    enddo
    do i=fb,fb+f-1
      do j=fep,Glob_CurrBasisSize
        c=Glob_S(j,i)
        Glob_S(j,i)=Glob_S(j,ff-i)
        Glob_S(j,ff-i)=c
      enddo
    enddo
    if (Glob_GSEPSolutionMethod=='I') then
      !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      !upper triangle of Glob_S
      do i=fb,fe
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
      enddo
    endif

  end subroutine ReverseMatElemOrder

  subroutine MakeFuncAndMEPermForCyclicOpt(fb,fe,blacklisted,blsize,FuncNumTemp,NonlinParamTemp, &
                                           TempR,fbnew,Permute_ME)
!Subroutine MakeFuncAndMEPermForCyclicOpt changes the order of basis functions in
!such a way that it is suitable for start of optimization cycles
!in cyclic optimization routines. Basically it moves functions fb through fe to
!the very end of the list while also reversing their order and excluding the
!functions that are blacklisted. Functions 1 through fb-1 remain completely
!untouched. The blacklist (logical array blacklisted) has only blsize elements,
!which can be smaller than the current basis size. Workspace must be supplied
!in arrays FuncNumTemp and NonlinParamTemp. The length of FuncNumTemp
!should be at least Glob_CurrBasisSize-fb+1, while the length of NonlinParamTemp
!must be at least (Glob_CurrBasisSize-fb+1)*Glob_npt. If some functions from fb
!to fe happen to be blacklisted then the new number of functions that need to be
!optimized is smaller than fe-fb+1. Upon exit the subroutine provides a new start
!function for cyclic optimization, fbnew. The new end function is, of course,
!the last function in the basis, Glob_CurrBasisSize.
!When Permute_ME=.true. this subroutine also changes the order of matrix elements
!accordingly. Otherwise nothing is done with matrix elements. When Permute_ME=.true.
!the user needs to supply additional workspace in array TempR. The length of this array
!should be at least Glob_CurrBasisSize-fb+1. If Permute_ME=.false. this array is not
!referenced.
!Arguments:
    integer        fb,fe
    logical        blacklisted(*)
    integer        blsize
    integer        FuncNumTemp(*)
    real(wp)    NonlinParamTemp(Glob_npt,*)
    real(wp)    TempR(*)
    integer        fbnew
    logical        Permute_ME
!Local variables:
    integer        i,j,cbs,cbs1,cbs1mi,q,p,fbm1,gfni,gfnj

    cbs=Glob_CurrBasisSize
    cbs1=cbs+1
    fbm1=fb-1
    q=cbs-fe
    p=cbs-fb+1
    j=q

!First we permute function numbers:
    do i=1,q
      FuncNumTemp(i)=Glob_FuncNum(cbs1-i)
    enddo
    do i=q+1,p
      cbs1mi=cbs1-i
      if (cbs1mi<=blsize) then
        if (blacklisted(cbs1mi)) then
          j=j+1
          FuncNumTemp(j)=Glob_FuncNum(cbs1mi)
        endif
      endif
    enddo
    fbnew=fb+j
    do i=q+1,p
      cbs1mi=cbs1-i
      if (cbs1mi<=blsize) then
        if (.not.blacklisted(cbs1mi)) then
          j=j+1
          FuncNumTemp(j)=Glob_FuncNum(cbs1mi)
        endif
      else
        j=j+1
        FuncNumTemp(j)=Glob_FuncNum(cbs1mi)
      endif
    enddo
    Glob_FuncNum(fb:cbs)=FuncNumTemp(1:p)
!Then we permute nonlinear parameters
    do i=fb,cbs
      NonlinParamTemp(1:Glob_npt,i-fbm1)=Glob_NonlinParam(1:Glob_npt,Glob_FuncNum(i))
    enddo
    do i=fb,cbs
      Glob_NonlinParam(1:Glob_npt,i)=NonlinParamTemp(1:Glob_npt,i-fbm1)
    enddo

!We also permute x,y,z-indices of functions
    do i=fb,cbs
      FuncNumTemp(i-fbm1)=Glob_Index(Glob_FuncNum(i),1)
    enddo
    do i=fb,cbs
      Glob_Index(i,1)=FuncNumTemp(i-fbm1)
    enddo

    do i=fb,cbs
      FuncNumTemp(i-fbm1)=Glob_Index(Glob_FuncNum(i),2)
    enddo
    do i=fb,cbs
      Glob_Index(i,2)=FuncNumTemp(i-fbm1)
    enddo

!Exit if matrix elements are not needed to be permuted
    if (.not.Permute_ME) return

!Permute diagonal matrix elements
    if (Glob_GSEPSolutionMethod=='G') then
      do i=fb,cbs
        TempR(i-fbm1)=Glob_diagH(Glob_FuncNum(i))
      enddo
      do i=fb,cbs
        Glob_diagH(i)=TempR(i-fbm1)
      enddo
    else
      do i=fb,cbs
        TempR(i-fbm1)=Glob_H(Glob_FuncNum(i),Glob_FuncNum(i))
      enddo
      do i=fb,cbs
        Glob_H(i,i)=TempR(i-fbm1)
      enddo
    endif
    do i=fb,cbs
      TempR(i-fbm1)=Glob_diagS(Glob_FuncNum(i))
    enddo
    do i=fb,cbs
      Glob_diagS(i)=TempR(i-fbm1)
    enddo
!Permute off-diagonal matrix elements
    do i=fb,cbs
      do j=1,i-1
        gfni=Glob_FuncNum(i)
        gfnj=Glob_FuncNum(j)
        if (gfni>gfnj) then
          Glob_H(j,i)=Glob_H(gfni,gfnj)
          Glob_S(j,i)=Glob_S(gfni,gfnj)
        else
          Glob_H(j,i)=Glob_H(gfnj,gfni)
          Glob_S(j,i)=Glob_S(gfnj,gfni)
        endif
      enddo
    enddo
    do i=fb,cbs
      do j=1,i-1
        Glob_H(i,j)=Glob_H(j,i)
        Glob_S(i,j)=Glob_S(j,i)
      enddo
    enddo

  end subroutine MakeFuncAndMEPermForCyclicOpt

  subroutine SortBasisFuncAndMatElem(fb,fe,FuncNumTemp,NonlinParamTemp,TempR)
!Subroutine SortBasisFuncAndMatElem permutes basis functions
!fb through fe so that they are sorted in decreasing order
!(that is the values of Glob_FuncNum(i) decrease). It also permutes
!the corresponding matrix elements. The subroutine permutes the
!nonlinear parameters (Glob_NonlinParam), x,y,z-indices (Glob_Index),
!function numbers (Glob_FuncNum), the elements of arrays Glob_H, Glob_S,
!Glob_diagH, and glob_diagS. It is assumed that the set of
!values Glob_FuncNum(fb:fe) ranges from some minimal value
!to the maximal one with no gaps (for example:  5,8,6,7,9). If
!this condition is not satisfied the subroutine will fail without
!warning.
!This subroutine requires workspace, which must be provided by
!arrays FuncNumTemp, NonlinParamTemp, and TempR. The length of
!FuncNumTemp, TempR, TempC should be at least fe-fb+1, while the
!length of NonlinParamTemp must be at least (fe-fb+1)*Glob_npt.
!Note that the subroutine uses only the lower triangles of matrices
!Glob_H and Glob_S (excluding the diagonals). Parts of the upper
!triangles (excluding the diagonal) are used as workspace.

!Arguments:
    integer       fb,fe
    integer       FuncNumTemp(*)
    real(wp)       NonlinParamTemp(Glob_npt,*)
    real(wp)       TempR(*)
!Local variables:
    integer  i,j,k,mf,fep,nfco,fbm,cbs,fi,fj

    cbs=Glob_CurrBasisSize
    nfco=fe-fb+1
    fep=fe+1
    fbm=fb-1
    mf=minval(Glob_FuncNum(fb:fe))
    k=fbm+nfco+mf
!First we sort out nonlinear parameters and x,y-indices
    do i=fb,fe
      NonlinParamTemp(1:Glob_npt,nfco+mf-Glob_FuncNum(i))=Glob_NonlinParam(1:Glob_npt,i)
    enddo
    Glob_NonlinParam(1:Glob_npt,fb:fe)=NonlinParamTemp(1:Glob_npt,1:nfco)

    do i=fb,fe
      TempR(nfco+mf-Glob_FuncNum(i))=Glob_Index(i,1)
    enddo
    Glob_Index(fb:fe,1)=TempR(1:nfco)

    do i=fb,fe
      TempR(nfco+mf-Glob_FuncNum(i))=Glob_Index(i,2)
    enddo
    Glob_Index(fb:fe,2)=TempR(1:nfco)

!Then we sort out the diagonal matrix elements
    if (Glob_GSEPSolutionMethod=='G') then
      do i=fb,fe
        TempR(nfco+mf-Glob_FuncNum(i))=Glob_diagH(i)
      enddo
      Glob_diagH(fb:fe)=TempR(1:nfco)
    else
      do i=fb,fe
        TempR(nfco+mf-Glob_FuncNum(i))=Glob_H(i,i)
      enddo
      do i=fb,fe
        Glob_H(i,i)=TempR(i-fb+1)
      enddo
    endif
    do i=fb,fe
      TempR(nfco+mf-Glob_FuncNum(i))=Glob_diagS(i)
    enddo
    Glob_diagS(fb:fe)=TempR(1:nfco)
!Sorting out off-diagonal matrix elements
    do i=fb,fe
      Glob_H(1:fbm,k-Glob_FuncNum(i))=Glob_H(i,1:fbm)
    enddo
    do i=fb,fe
      Glob_H(i,1:fbm)=Glob_H(1:fbm,i)
    enddo
    do i=fb,fe
      Glob_H(k-Glob_FuncNum(i),fep:cbs)=Glob_H(fep:cbs,i)
    enddo
    do i=fb,fe
      Glob_H(fep:cbs,i)=Glob_H(i,fep:cbs)
    enddo
    do i=fb,fe-1
      do j=i+1,fe
        fi=Glob_FuncNum(i)
        fj=Glob_FuncNum(j)
        if (fi<fj) then
          Glob_H(k-fj,k-fi)=Glob_H(j,i)
        else
          Glob_H(k-fi,k-fj)=Glob_H(j,i)
        endif
      enddo
    enddo
    do i=fb,fe-1
      Glob_H(i+1:fe,i)=Glob_H(i,i+1:fe)
    enddo
    do i=fb,fe
      Glob_S(1:fbm,k-Glob_FuncNum(i))=Glob_S(i,1:fbm)
    enddo
    do i=fb,fe
      Glob_S(i,1:fbm)=Glob_S(1:fbm,i)
    enddo
    do i=fb,fe
      Glob_S(k-Glob_FuncNum(i),fep:cbs)=Glob_S(fep:cbs,i)
    enddo
    do i=fb,fe
      Glob_S(fep:cbs,i)=Glob_S(i,fep:cbs)
    enddo
    do i=fb,fe-1
      do j=i+1,fe
        fi=Glob_FuncNum(i)
        fj=Glob_FuncNum(j)
        if (fi<fj) then
          Glob_S(k-fj,k-fi)=Glob_S(j,i)
        else
          Glob_S(k-fi,k-fj)=Glob_S(j,i)
        endif
      enddo
    enddo
    do i=fb,fe-1
      Glob_S(i+1:fe,i)=Glob_S(i,i+1:fe)
    enddo
    if (Glob_GSEPSolutionMethod=='I') then
      !In case Glob_GSEPSolutionMethod=='I' we need to fill out the
      !upper triangle of Glob_S
      do i=fb,fe
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
      enddo
    endif
!At last we change the order of basis functions
    do i=fb,fe
      Glob_FuncNum(i)=mf+fe-i
    enddo

  end subroutine SortBasisFuncAndMatElem

  subroutine GetOverlapStatistics(Nmin,Nmax,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
!Subroutine GetOverlapStatistics determines the largest by magnitude overlap, the smallest
!by magnitude overlap, and the average overlap magnitude for basis functions
!ranging from Nmin to Nmax.
!Arguments:
    integer      Nmin,Nmax
    real(wp)  MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap
!Local variables
    integer      i,j,k,nbands,leftover
    logical      oddband
    real(wp)  absMaxAbsOverlap,absMinAbsOverlap,absSji,t

    MaxAbsOverlap=ZERO
    absMaxAbsOverlap=ZERO
    MinAbsOverlap=huge(MinAbsOverlap)
    absMinAbsOverlap=MinAbsOverlap
    t=ZERO

    do i=1,Nmin-1
      do j=Nmin,Nmax
        absSji=abs(Glob_S(j,i))
        if (absSji>absMaxAbsOverlap) then
          absMaxAbsOverlap=absSji
          MaxAbsOverlap=Glob_S(j,i)
        endif
        if (absSji<absMinAbsOverlap) then
          absMinAbsOverlap=absSji
          MinAbsOverlap=Glob_S(j,i)
        endif
        t=t+absSji
      enddo
    enddo
    do i=Nmin,Nmax
      do j=i+1,Nmax
        absSji=abs(Glob_S(j,i))
        if (absSji>absMaxAbsOverlap) then
          absMaxAbsOverlap=absSji
          MaxAbsOverlap=Glob_S(j,i)
        endif
        if (absSji<absMinAbsOverlap) then
          absMinAbsOverlap=absSji
          MinAbsOverlap=Glob_S(j,i)
        endif
        t=t+absSji
      enddo
    enddo
    AverageAbsOverlap=2*t/(Nmax*(Nmax-1)-(Nmin-1)*(Nmin-2))

  end subroutine GetOverlapStatistics

  function NumOfRowsToPermForUnitShift(j)
!Function NumOfRowsToPerm is used in cyclic optimization
!procedure to find the number of function that must be
!permuted. What it does, it returns the following values:
!
!      j  NumOfRowsToPermForUnitShift(j)
!
!      1           1
!      2           2
!      3           1
!      4           4
!      5           1
!      6           2
!      7           1
!      8           8
!      9           1
!     10           2
!     11           1
!     12           4
!     13           1
!     14           2
!     15           1
!     16          16
!     17           1
!      .           .
!      .           .
!      .           .       .

    integer NumOfRowsToPermForUnitShift,j,k
    k=1
    do while (mod(j,k)==0)
      k=k*2
    enddo
    NumOfRowsToPermForUnitShift=k/2

  end function NumOfRowsToPermForUnitShift

  subroutine GenerateRndIntSeq(n,s)
!Routine GenerateRndIntSeq generates a random sequence (permutation)
!of integers from 1 to n. The result is returned in array s(1:n)
!Arguments:
    integer      n,s(n)
!Local variables:
    integer      i,j,k
    real(8)  r8

    if (n==1) then
      s(1)=1
    else
      do i=1,n
        s(i)=i
      enddo
      do i=1,n
        call random_number(r8)
        j=int(r8*i)+1
        k=s(i); s(i)=s(j); s(j)=k
      enddo
    endif

  end subroutine GenerateRndIntSeq

  subroutine ReallocateBasisFuncData(FinalSize,NumOfFuncToKeep)
!Subroutine ReallocateBasisFuncData reallocates arrays that contain
!basis function data. The final arrays may be either larger than
!the initial ones or smaller. The affected arrays are:
!   Glob_NonlinParam
!   Glob_FuncNum
!   Glob_History
!The parameters of the subroutine are the following:
!   FinalSize - the final size of the arrays
!   NumOfFuncToKeep - the number of functions in the initial array
!       whose data will be copied to the reallocated arrays. More
!       specifically, the data of first NumOfFuncToKeep functions
!       is copied. Notice that NumOfFuncToKeep<=FinalSize.
!Parameters:
    integer       FinalSize,NumOfFuncToKeep
!Local variables:
    integer                                         :: i,OpenFileErr
    real(wp),allocatable,dimension(:)            :: WorkBuffReal
    integer,allocatable,dimension(:)                :: WorkBuffInt
    type(Glob_HistoryStep),allocatable,dimension(:) :: TempHistory
    real(wp),allocatable,dimension(:,:)          :: TempParam
    integer,allocatable,dimension(:)                :: TempFunc
    integer,allocatable,dimension(:,:)              :: TempInd

    if (NumOfFuncToKeep>FinalSize) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0120 in ReallocateBasisFuncData:'
        write(*,*) 'NumOfFuncToKeep must be smaller or equal than FinalSize'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    if (Glob_UseReallocFile) then
      !Temporarily store the information in a file
      if (Glob_ProcID==0) then
        if (NumOfFuncToKeep>0) then
          write(*,'(1x,a47)',advance='no') 'Reallocating some arrays using external file...'
          open(1,file=Glob_ReallocFileName,form='unformatted',status='replace')
          write(1) Glob_History(1:NumOfFuncToKeep)
          write(1) Glob_FuncNum(1:NumOfFuncToKeep)
          write(1) Glob_Index(1:NumOfFuncToKeep,1)
          write(1) Glob_Index(1:NumOfFuncToKeep,2)
          write(1) Glob_NonlinParam(1:Glob_npt,1:NumOfFuncToKeep)
          close(1)
        endif
      endif
      deallocate(Glob_NonlinParam)
      deallocate(Glob_Index)
      deallocate(Glob_FuncNum)
      deallocate(Glob_History)
      allocate(Glob_History(FinalSize))
      allocate(Glob_FuncNum(FinalSize))
      allocate(Glob_Index(FinalSize,2))
      allocate(Glob_NonlinParam(Glob_npt,FinalSize))
      if (Glob_ProcID==0) then
        if (NumOfFuncToKeep>0) then
          open(1,file=Glob_ReallocFileName,form='unformatted',status='old',iostat=OpenFileErr)
          if (OpenFileErr==0) then
            read(1) Glob_History(1:NumOfFuncToKeep)
            read(1) Glob_FuncNum(1:NumOfFuncToKeep)
            read(1) Glob_Index(1:NumOfFuncToKeep,1)
            read(1) Glob_Index(1:NumOfFuncToKeep,2)
            read(1,iostat=OpenFileErr) Glob_NonlinParam(1:Glob_npt,1:NumOfFuncToKeep)
          endif
          close(1)
          call MPI_BCAST(OpenFileErr,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        endif
      endif
      if ((OpenFileErr/=0).and.(NumOfFuncToKeep>0)) then
        if (Glob_ProcID==0) then
          write(*,*)
          write(*,*) 'Error EC0121 in ReallocateBasisFuncData:'
          write(*,*) 'cannot read data from file',Glob_ReallocFileName
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
      if (Glob_ProcID==0) then
        open(1,file=Glob_ReallocFileName,form='unformatted',status='replace',iostat=OpenFileErr)
        write(1)  'This temporary file is empty'
        close(1)
      endif
      allocate(WorkBuffReal(NumOfFuncToKeep))
      allocate(WorkBuffInt(NumOfFuncToKeep))
      do i=1,NumOfFuncToKeep
        WorkBuffReal(i)=Glob_History(i)%Energy
      enddo
      call MPI_BCAST(WorkBuffReal,NumOfFuncToKeep,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do i=1,NumOfFuncToKeep
        Glob_History(i)%Energy=WorkBuffReal(i)
      enddo
      do i=1,NumOfFuncToKeep
        WorkBuffInt(i)=Glob_History(i)%CyclesDone
      enddo
      call MPI_BCAST(WorkBuffInt,NumOfFuncToKeep,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do i=1,NumOfFuncToKeep
        Glob_History(i)%CyclesDone=WorkBuffInt(i)
      enddo
      do i=1,NumOfFuncToKeep
        WorkBuffInt(i)=Glob_History(i)%InitFuncAtLastStep
      enddo
      call MPI_BCAST(WorkBuffInt,NumOfFuncToKeep,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do i=1,NumOfFuncToKeep
        Glob_History(i)%InitFuncAtLastStep=WorkBuffInt(i)
      enddo
      do i=1,NumOfFuncToKeep
        WorkBuffInt(i)=Glob_History(i)%NumOfEnergyEvalDuringFullOpt
      enddo
      call MPI_BCAST(WorkBuffInt,NumOfFuncToKeep,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      do i=1,NumOfFuncToKeep
        Glob_History(i)%NumOfEnergyEvalDuringFullOpt=WorkBuffInt(i)
      enddo
      deallocate(WorkBuffReal)
      deallocate(WorkBuffInt)
      call MPI_BCAST(Glob_FuncNum,NumOfFuncToKeep,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_Index,2*NumOfFuncToKeep,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_NonlinParam,Glob_npt*NumOfFuncToKeep, &
                     MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_ProcID==0) then
        write(*,*) 'done'
      endif
    else !if (Glob_UseReallocFile)
      allocate(TempHistory(NumOfFuncToKeep))
      allocate(TempFunc(NumOfFuncToKeep))
      allocate(TempInd(NumOfFuncToKeep,2))
      allocate(TempParam(Glob_npt,NumOfFuncToKeep))
      TempHistory(1:NumOfFuncToKeep)=Glob_History(1:NumOfFuncToKeep)
      TempFunc(1:NumOfFuncToKeep)=Glob_FuncNum(1:NumOfFuncToKeep)
      TempInd(1:NumOfFuncToKeep,1)=Glob_Index(1:NumOfFuncToKeep,1)
      TempInd(1:NumOfFuncToKeep,2)=Glob_Index(1:NumOfFuncToKeep,2)
      TempParam(1:Glob_npt,1:NumOfFuncToKeep)=Glob_NonlinParam(1:Glob_npt,1:NumOfFuncToKeep)
      deallocate(Glob_NonlinParam)
      deallocate(Glob_Index)
      deallocate(Glob_FuncNum)
      deallocate(Glob_History)
      allocate(Glob_History(FinalSize))
      allocate(Glob_FuncNum(FinalSize))
      allocate(Glob_Index(FinalSize,2))
      allocate(Glob_NonlinParam(Glob_npt,FinalSize))
      Glob_History(1:NumOfFuncToKeep)=TempHistory(1:NumOfFuncToKeep)
      Glob_FuncNum(1:NumOfFuncToKeep)=TempFunc(1:NumOfFuncToKeep)
      Glob_Index(1:NumOfFuncToKeep,1)=TempInd(1:NumOfFuncToKeep,1)
      Glob_Index(1:NumOfFuncToKeep,2)=TempInd(1:NumOfFuncToKeep,2)
      Glob_NonlinParam(1:Glob_npt,1:NumOfFuncToKeep)=TempParam(1:Glob_npt,1:NumOfFuncToKeep)
      deallocate(TempParam)
      deallocate(TempInd)
      deallocate(TempFunc)
      deallocate(TempHistory)
    endif

    do i=NumOfFuncToKeep+1,FinalSize
      Glob_History(i)%Energy=ZERO
      Glob_History(i)%CyclesDone=0
      Glob_History(i)%InitFuncAtLastStep=0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
      Glob_FuncNum(i)=0
      Glob_Index(i,1)=0
      Glob_Index(i,2)=0
      Glob_NonlinParam(1:Glob_npt,i)=ZERO
    enddo

  end subroutine ReallocateBasisFuncData

  subroutine PrepareQWorkspace(MatrixOrder,Capacity,MaxActive,ErrorCode)
!Subroutine PrepareQWorkspace allocates the workproc-owned storage shared by
!the Q energy routines and one Q BBOP driver. It does not allocate Glob_H or
!Glob_S and does not initialize qrlinalg. The BBOP driver owns those lifetimes
!and must call this routine only after Glob_npt has been initialized.
!
!The physical matrices can have capacity greater than MatrixOrder during basis
!enlargement. Only indices 1 through MatrixOrder belong to the represented
!problem. MaxActive is normally Kstep for BASIS_ENL, NumOfFuncToOpt for
!OPT_CYCLE, and FinalFunc-InitFunc+1 for FULL_OPT1.
!
!Arguments:
    integer,intent(in)  :: MatrixOrder,Capacity,MaxActive
    integer,intent(out) :: ErrorCode
!Local variables:
    integer AllocationStatus

    call ClearQWorkspace()
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (MatrixOrder<0) return
    if (Capacity<max(1,MatrixOrder)) return
    if (MaxActive<1) return
    if (Glob_npt<1) return

    allocate(Q_Workspace%ActiveFunction(MaxActive), &
             Q_Workspace%ActivePosition(Capacity), &
             Q_Workspace%PreviousH(Capacity,MaxActive), &
             Q_Workspace%PreviousS(Capacity,MaxActive), &
             Q_Workspace%TrialH(Capacity,MaxActive), &
             Q_Workspace%TrialS(Capacity,MaxActive), &
             Q_Workspace%PreviousDiagS(MaxActive), &
             Q_Workspace%TrialDiagS(MaxActive), &
             Q_Workspace%MatrixParam(Glob_npt,MaxActive), &
             Q_Workspace%PreviousParam(Glob_npt,MaxActive), &
             Q_Workspace%AcceptedParam(Glob_npt,MaxActive), &
             Q_Workspace%InitialVector(Capacity), &
             Q_Workspace%SolvedVector(Capacity), &
             Q_Workspace%DeltaH(Capacity), &
             Q_Workspace%DeltaS(Capacity), &
             stat=AllocationStatus)
    if (AllocationStatus/=0) then
      call ClearQWorkspace()
      ErrorCode=Q_METHOD_ALLOCATION_ERROR
      return
    endif

    Q_Workspace%MatrixOrder=MatrixOrder
    Q_Workspace%Capacity=Capacity
    Q_Workspace%MaxActive=MaxActive
    Q_Workspace%NumActive=0
    Q_Workspace%ActiveFunction=0
    Q_Workspace%ActivePosition=0
    Q_Workspace%PreviousH=ZERO
    Q_Workspace%PreviousS=ZERO
    Q_Workspace%TrialH=ZERO
    Q_Workspace%TrialS=ZERO
    Q_Workspace%PreviousDiagS=ZERO
    Q_Workspace%TrialDiagS=ZERO
    Q_Workspace%MatrixParam=ZERO
    Q_Workspace%PreviousParam=ZERO
    Q_Workspace%AcceptedParam=ZERO
    Q_Workspace%InitialVector=ONE
    Q_Workspace%SolvedVector=ZERO
    Q_Workspace%DeltaH=ZERO
    Q_Workspace%DeltaS=ZERO
    Q_Workspace%AcceptedEnergy=ZERO
    Q_Workspace%AcceptedPointIsStored=.false.
    Q_Workspace%LastEigenpairResidual=huge(ONE)
    Q_Workspace%LastFactorResidual=huge(ONE)
    Q_Workspace%FreshFactorizations=0
    Q_Workspace%MatricesAreCanonical=.false.
    Q_Workspace%FactorsMatchMatrices=.false.
    Q_Workspace%MatrixParametersAreStored=.false.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine PrepareQWorkspace

  subroutine EnsureQActiveCapacity(RequiredActive,ErrorCode)
!Subroutine EnsureQActiveCapacity enlarges only the transaction part of an
!existing Q workspace. The QR factors and the full-order vector workspace are
!deliberately preserved. Cleanup routines initially reserve one active column
!because that is sufficient for elimination, but a separation test discovers
!the number of columns to replace only after the first eigenvector or overlap
!matrix has been inspected. Rebuilding the complete Q workspace at that point
!would discard the factorization whose reuse is the purpose of the Q method.
!
!This operation is valid only between transactions. No previous, trial, or
!accepted point is copied: cleanup routines call it before selecting their
!first active set. Allocating every replacement array before move_alloc makes
!allocation failure transactional as well; the original workspace remains
!usable when any allocation fails.
!
!Arguments:
    integer,intent(in)  :: RequiredActive
    integer,intent(out) :: ErrorCode
!Local variables:
    integer AllocationStatus
    integer,allocatable :: NewActiveFunction(:)
    real(wp),allocatable :: NewPreviousH(:,:),NewPreviousS(:,:)
    real(wp),allocatable :: NewTrialH(:,:),NewTrialS(:,:)
    real(wp),allocatable :: NewPreviousDiagS(:),NewTrialDiagS(:)
    real(wp),allocatable :: NewMatrixParam(:,:),NewPreviousParam(:,:)
    real(wp),allocatable :: NewAcceptedParam(:,:)

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (RequiredActive<1) return
    if (RequiredActive>Q_Workspace%Capacity) return
    if (Q_Workspace%MatrixOrder<1) return
    if (Q_Workspace%NumActive/=0) return
    if (Q_Workspace%MatrixParametersAreStored) return
    if (Q_Workspace%TrialIsReady) return
    if (Q_Workspace%AcceptedPointIsStored) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (RequiredActive<=Q_Workspace%MaxActive) then
      ErrorCode=Q_METHOD_SUCCESS
      return
    endif

    allocate(NewActiveFunction(RequiredActive), &
             NewPreviousH(Q_Workspace%Capacity,RequiredActive), &
             NewPreviousS(Q_Workspace%Capacity,RequiredActive), &
             NewTrialH(Q_Workspace%Capacity,RequiredActive), &
             NewTrialS(Q_Workspace%Capacity,RequiredActive), &
             NewPreviousDiagS(RequiredActive), &
             NewTrialDiagS(RequiredActive), &
             NewMatrixParam(Glob_npt,RequiredActive), &
             NewPreviousParam(Glob_npt,RequiredActive), &
             NewAcceptedParam(Glob_npt,RequiredActive), &
             stat=AllocationStatus)
    if (AllocationStatus/=0) then
      ErrorCode=Q_METHOD_ALLOCATION_ERROR
      return
    endif

    NewActiveFunction=0
    NewPreviousH=ZERO
    NewPreviousS=ZERO
    NewTrialH=ZERO
    NewTrialS=ZERO
    NewPreviousDiagS=ZERO
    NewTrialDiagS=ZERO
    NewMatrixParam=ZERO
    NewPreviousParam=ZERO
    NewAcceptedParam=ZERO

    call move_alloc(NewActiveFunction,Q_Workspace%ActiveFunction)
    call move_alloc(NewPreviousH,Q_Workspace%PreviousH)
    call move_alloc(NewPreviousS,Q_Workspace%PreviousS)
    call move_alloc(NewTrialH,Q_Workspace%TrialH)
    call move_alloc(NewTrialS,Q_Workspace%TrialS)
    call move_alloc(NewPreviousDiagS,Q_Workspace%PreviousDiagS)
    call move_alloc(NewTrialDiagS,Q_Workspace%TrialDiagS)
    call move_alloc(NewMatrixParam,Q_Workspace%MatrixParam)
    call move_alloc(NewPreviousParam,Q_Workspace%PreviousParam)
    call move_alloc(NewAcceptedParam,Q_Workspace%AcceptedParam)
    Q_Workspace%MaxActive=RequiredActive
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine EnsureQActiveCapacity

  subroutine ClearQWorkspace()
!Subroutine ClearQWorkspace releases workproc-owned Q storage. It intentionally
!does not deallocate the existing Glob_ arrays, because their lifetime belongs
!to the BBOP driver. The qrlinalg state is cleared here before the remaining
!metadata are reset.

    if (allocated(Q_Workspace%DeltaS)) deallocate(Q_Workspace%DeltaS)
    if (allocated(Q_Workspace%DeltaH)) deallocate(Q_Workspace%DeltaH)
    if (allocated(Q_Workspace%SolvedVector)) deallocate(Q_Workspace%SolvedVector)
    if (allocated(Q_Workspace%InitialVector)) deallocate(Q_Workspace%InitialVector)
    if (allocated(Q_Workspace%AcceptedParam)) deallocate(Q_Workspace%AcceptedParam)
    if (allocated(Q_Workspace%PreviousParam)) deallocate(Q_Workspace%PreviousParam)
    if (allocated(Q_Workspace%MatrixParam)) deallocate(Q_Workspace%MatrixParam)
    if (allocated(Q_Workspace%TrialDiagS)) deallocate(Q_Workspace%TrialDiagS)
    if (allocated(Q_Workspace%PreviousDiagS)) deallocate(Q_Workspace%PreviousDiagS)
    if (allocated(Q_Workspace%TrialS)) deallocate(Q_Workspace%TrialS)
    if (allocated(Q_Workspace%TrialH)) deallocate(Q_Workspace%TrialH)
    if (allocated(Q_Workspace%PreviousS)) deallocate(Q_Workspace%PreviousS)
    if (allocated(Q_Workspace%PreviousH)) deallocate(Q_Workspace%PreviousH)
    if (allocated(Q_Workspace%ActivePosition)) deallocate(Q_Workspace%ActivePosition)
    !clear is valid for an initialized or empty state. Calling it on every MPI
    !rank is therefore safe even though only rank zero will own allocated QR
    !factors once FactorizeQFresh is implemented.
    call Q_Workspace%Factors%clear()

    if (allocated(Q_Workspace%ActiveFunction)) deallocate(Q_Workspace%ActiveFunction)

    Q_Workspace%MatrixOrder=0
    Q_Workspace%Capacity=0
    Q_Workspace%MaxActive=0
    Q_Workspace%NumActive=0
    Q_Workspace%MatricesAreCanonical=.false.
    Q_Workspace%FactorsMatchMatrices=.false.
    Q_Workspace%MatrixParametersAreStored=.false.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    Q_Workspace%AcceptedEnergy=ZERO
    Q_Workspace%AcceptedPointIsStored=.false.
    Q_Workspace%LastEigenpairResidual=huge(ONE)
    Q_Workspace%LastFactorResidual=huge(ONE)
    Q_Workspace%FreshFactorizations=0

  end subroutine ClearQWorkspace

  subroutine SetQActiveFunctions(ActiveFunction,ErrorCode)
!Subroutine SetQActiveFunctions defines optimizer block order without changing
!the physical order of basis functions. ActiveFunction must contain distinct
!canonical indices in the range 1:Q_Workspace%MatrixOrder. The order supplied
!here is also the order of nonlinear parameter blocks, gradient blocks, and
!saved Hessian rows and columns.
!
!The routine updates both maps only after the complete input has been checked,
!so an invalid selection leaves the previous active set unchanged.
!
!Arguments:
    integer,intent(in)  :: ActiveFunction(:)
    integer,intent(out) :: ErrorCode
!Local variables:
    integer i,j,NumActive

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%ActivePosition)) return
    NumActive=size(ActiveFunction)
    if ((NumActive<1).or.(NumActive>Q_Workspace%MaxActive)) return
    do i=1,NumActive
      if ((ActiveFunction(i)<1).or. &
          (ActiveFunction(i)>Q_Workspace%MatrixOrder)) return
      do j=1,i-1
        if (ActiveFunction(i)==ActiveFunction(j)) return
      enddo
    enddo

    Q_Workspace%ActiveFunction=0
    Q_Workspace%ActivePosition=0
    Q_Workspace%ActiveFunction(1:NumActive)=ActiveFunction
    do i=1,NumActive
      Q_Workspace%ActivePosition(ActiveFunction(i))=i
    enddo
    Q_Workspace%NumActive=NumActive
    Q_Workspace%AcceptedPointIsStored=.false.
    Q_Workspace%MatrixParametersAreStored=.false.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine SetQActiveFunctions

  subroutine CaptureQMatrixParameters(ErrorCode)
!Subroutine CaptureQMatrixParameters records the nonlinear parameters that
!belong to the currently represented canonical H/S matrices. A Q driver calls
!this immediately after selecting a new active set, before it copies a DRMNG
!trial point into Glob_NonlinParam. The separate copy is essential: after that
!copy Glob_NonlinParam describes the requested trial, while Glob_H and Glob_S
!still describe the preceding point until ApplyQTrial commits successfully.
!
!This routine cannot prove that matrix elements were calculated from the
!current parameters. MatricesAreCanonical is therefore an explicit caller
!precondition set only after a full Q assembly or a successful swap restore.
!
!Arguments:
    integer,intent(out) :: ErrorCode
!Local variables:
    integer a,FunctionIndex

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (Q_Workspace%NumActive<1) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%MatrixParam)) return
    if (.not.allocated(Glob_NonlinParam)) return
    if (size(Glob_NonlinParam,1)<Glob_npt) return
    if (size(Glob_NonlinParam,2)<Q_Workspace%MatrixOrder) return

    do a=1,Q_Workspace%NumActive
      FunctionIndex=Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt,a)= &
        Glob_NonlinParam(1:Glob_npt,FunctionIndex)
    enddo
    if (Q_Workspace%NumActive<Q_Workspace%MaxActive) then
      Q_Workspace%MatrixParam(1:Glob_npt,Q_Workspace%NumActive+1:Q_Workspace%MaxActive)=ZERO
    endif
    Q_Workspace%MatrixParametersAreStored=.true.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine CaptureQMatrixParameters

  subroutine TrimQFactors(TargetOrder,ErrorCode)
!Subroutine TrimQFactors removes a rejected canonical suffix from the QR state.
!BASIS_ENL repeatedly tests new functions in positions TargetOrder+1 onward;
!deleting those last rows and columns restores the accepted prefix without any
!basis permutation or cubic refactorization. The physical suffix may retain
!obsolete trial values because it lies outside MatrixOrder and is overwritten
!before it can become active again.
!
!qrlinalg intentionally has no valid order-zero factorization. When the first
!basis block is being selected, trimming to zero therefore reinitializes an
!empty capacity reservation. The next trial is constructed by a fresh
!factorization after its complete physical block has been staged.
!
!Arguments:
    integer,intent(in)  :: TargetOrder
    integer,intent(out) :: ErrorCode
!Local variables:
    integer CurrentOrder,RootError,RecoveryError

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    CurrentOrder=Q_Workspace%MatrixOrder
    if ((TargetOrder<0).or.(TargetOrder>CurrentOrder)) return
    if (Q_Workspace%Capacity<max(1,CurrentOrder)) return

    RootError=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) then
      if (TargetOrder==0) then
        call Q_Workspace%Factors%clear()
        call Q_Workspace%Factors%initialize(Q_Workspace%Capacity,RootError)
      else
        if (.not.Q_Workspace%Factors%is_valid()) then
          RootError=Q_METHOD_INVALID_STATE
        else if (Q_Workspace%Factors%order()/=CurrentOrder) then
          RootError=Q_METHOD_INVALID_STATE
        endif
        if (RootError==Q_METHOD_SUCCESS) then
          do while (Q_Workspace%Factors%order()>TargetOrder)
            call Q_Workspace%Factors%delete_symmetric( &
              Q_Workspace%Factors%order(),RootError)
            if (RootError/=Q_METHOD_SUCCESS) exit
          enddo
        endif
      endif
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    !Publish the smaller physical order before recovery. Its leading principal
    !block was never modified by appending or deleting a suffix, so it is a
    !sound source for a fresh factorization if a structural operation failed.
    Q_Workspace%MatrixOrder=TargetOrder
    Q_Workspace%NumActive=0
    Q_Workspace%ActiveFunction=0
    Q_Workspace%ActivePosition=0
    Q_Workspace%MatrixParametersAreStored=.false.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    Q_Workspace%AcceptedPointIsStored=.false.
    Q_Workspace%FactorsMatchMatrices= &
      (RootError==Q_METHOD_SUCCESS).and.(TargetOrder>0)

    if ((RootError/=Q_METHOD_SUCCESS).and.(TargetOrder>0)) then
      call FactorizeQFresh(RecoveryError)
      if (RecoveryError/=Q_METHOD_SUCCESS) then
        ErrorCode=RecoveryError
      else
        ErrorCode=RootError
      endif
      return
    endif

    ErrorCode=RootError

  end subroutine TrimQFactors

  function QCanonicalMatrixElement(Matrix,i,j) result(MatrixElement)
!Function QCanonicalMatrixElement reads a conceptual symmetric matrix element
!from a matrix whose lower triangle, including the diagonal, is authoritative.
!No Q routine may read the upper triangle directly because it may contain old
!G-solver workspace or undefined values.
!
!Arguments:
    real(wp),intent(in) :: Matrix(:,:)
    integer,intent(in)  :: i,j
    real(wp) MatrixElement

    if (i>=j) then
      MatrixElement=Matrix(i,j)
    else
      MatrixElement=Matrix(j,i)
    endif

  end function QCanonicalMatrixElement

  subroutine GatherQCanonicalColumn(Matrix,MatrixOrder,ColumnIndex,Column,ErrorCode)
!Subroutine GatherQCanonicalColumn constructs the complete conceptual column
!required by qrlinalg replace_symmetric and append_symmetric. The source matrix
!retains only its canonical lower triangle. This O(n) gather is negligible
!beside the O(n*n) QR update and avoids maintaining two writable triangles.
!
!Arguments:
    real(wp),intent(in)  :: Matrix(:,:)
    integer,intent(in)   :: MatrixOrder,ColumnIndex
    real(wp),intent(out) :: Column(:)
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer i

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (MatrixOrder<1) return
    if ((size(Matrix,1)<MatrixOrder).or.(size(Matrix,2)<MatrixOrder)) return
    if ((ColumnIndex<1).or.(ColumnIndex>MatrixOrder)) return
    if (size(Column)/=MatrixOrder) return

    do i=1,MatrixOrder
      Column(i)=QCanonicalMatrixElement(Matrix,i,ColumnIndex)
    enddo
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine GatherQCanonicalColumn

  subroutine StoreQCanonicalColumn(Matrix,MatrixOrder,ColumnIndex,Column,ErrorCode)
!Subroutine StoreQCanonicalColumn commits one complete symmetric column to the
!canonical lower triangle. The upper triangle is intentionally untouched.
!For a replacement this routine is called only after qrlinalg has accepted the
!matching delta; during a multi-column transaction it is called after every
!successful sequential update so intersections use the progressively updated
!matrix and are not counted twice.
!
!Arguments:
    real(wp),intent(inout) :: Matrix(:,:)
    integer,intent(in)     :: MatrixOrder,ColumnIndex
    real(wp),intent(in)    :: Column(:)
    integer,intent(out)    :: ErrorCode
!Local variables:
    integer i

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (MatrixOrder<1) return
    if ((size(Matrix,1)<MatrixOrder).or.(size(Matrix,2)<MatrixOrder)) return
    if ((ColumnIndex<1).or.(ColumnIndex>MatrixOrder)) return
    if (size(Column)/=MatrixOrder) return

    do i=1,ColumnIndex
      Matrix(ColumnIndex,i)=Column(i)
    enddo
    do i=ColumnIndex+1,MatrixOrder
      Matrix(i,ColumnIndex)=Column(i)
    enddo
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine StoreQCanonicalColumn

  subroutine AssembleQTrial(AreDerivativesNeeded,ErrorCode)
!Subroutine AssembleQTrial calculates every unordered matrix-element pair
!that touches an active function. Active diagonals must be calculated first so
!all raw norms are available before normalized off-diagonal elements are
!formed. The final target columns are staged in Q_Workspace%TrialH and TrialS;
!Glob_H and Glob_S remain unchanged until ApplyQTrial succeeds.
!
!When derivatives are requested, Glob_D(:,a,j) will mean the G-compatible
!scaled derivatives with respect to canonical function
!Q_Workspace%ActiveFunction(a), paired with canonical function j. This removes
!the old assumption that the differentiated functions occupy a trailing block.
!
!Arguments:
    logical,intent(in)  :: AreDerivativesNeeded
    integer,intent(out) :: ErrorCode
!Local variables:
    integer a,b,i,j,q,PairNumber,MatrixOrder,NumActive
    integer ActiveIndex,NumMatrixEntries,NumDerivativeEntries,npt2
    integer mActive,mmActive,mOther,mmOther
    real(wp) ParamActive(Glob_AllowedNumOfPseudoParticles* &
                         (Glob_AllowedNumOfPseudoParticles+1)/2)
    real(wp) ParamOther(Glob_AllowedNumOfPseudoParticles* &
                        (Glob_AllowedNumOfPseudoParticles+1)/2)
    real(wp) Hkl,Skl,Tkl,Vkl,Hsum,Ssum,ActiveNorm,OtherNorm,Normalization
    real(wp) DActive(2*Glob_npt_MaxAllowed),DOther(2*Glob_npt_MaxAllowed)
    real(wp) DActiveSum(2*Glob_npt_MaxAllowed),DOtherSum(2*Glob_npt_MaxAllowed)
    logical OtherDerivativeNeeded

    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    NumActive=Q_Workspace%NumActive
    npt2=2*Glob_npt
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if ((MatrixOrder<1).or.(NumActive<1)) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%ActivePosition)) return
    if (.not.allocated(Q_Workspace%TrialH)) return
    if (.not.allocated(Q_Workspace%TrialS)) return
    if (.not.allocated(Q_Workspace%PreviousH)) return
    if (.not.allocated(Q_Workspace%PreviousS)) return
    if (.not.allocated(Q_Workspace%TrialDiagS)) return
    if (.not.allocated(Glob_NonlinParam)) return
    if (.not.allocated(Glob_diagS)) return
    if (size(Glob_NonlinParam,1)<Glob_npt) return
    if (size(Glob_NonlinParam,2)<MatrixOrder) return
    if (size(Glob_diagS)<MatrixOrder) return
    if (Glob_NumYHYTerms<1) return
    if (AreDerivativesNeeded) then
      if (.not.allocated(Glob_D)) return
      if (size(Glob_D,1)<npt2) return
      if (size(Glob_D,2)<NumActive) return
      if (size(Glob_D,3)<MatrixOrder) return
      Glob_D=ZERO
    endif

    !Each physical pair that touches the active set is evaluated exactly once.
    !For an active-active pair, the value is copied into both conceptual trial
    !columns. ActivePosition supplies an ordering independent of canonical
    !basis indices, so this remains correct for descending and noncontiguous
    !active lists.
    Q_Workspace%TrialH=ZERO
    Q_Workspace%TrialS=ZERO
    PairNumber=0
    do a=1,NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      ParamActive(1:Glob_npt)=Glob_NonlinParam(1:Glob_npt,ActiveIndex)
      mActive=Glob_Index(ActiveIndex,1)
      mmActive=Glob_Index(ActiveIndex,2)
      do i=1,MatrixOrder
        b=Q_Workspace%ActivePosition(i)
        if ((b>0).and.(b<a)) cycle

        PairNumber=PairNumber+1
        ParamOther(1:Glob_npt)=Glob_NonlinParam(1:Glob_npt,i)
        mOther=Glob_Index(i,1)
        mmOther=Glob_Index(i,2)
        Hsum=ZERO
        Ssum=ZERO
        if (AreDerivativesNeeded) then
          DActiveSum(1:npt2)=ZERO
          DOtherSum(1:npt2)=ZERO
        endif
        OtherDerivativeNeeded=AreDerivativesNeeded.and.(b>0).and.(b/=a)
        q=(PairNumber-1)*Glob_NumYHYTerms-1
        do j=1,Glob_NumYHYTerms
          if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then
            call MatrixElementsHS_RG_2D(mActive,mmActive,mOther,mmOther, &
              ParamActive,ParamOther, &
              Glob_YHYMatr(1:Glob_n,1:Glob_n,j),Hkl,Skl, &
              Tkl,Vkl,DActive,DOther,AreDerivativesNeeded,OtherDerivativeNeeded)
            Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
            Ssum=Ssum+Glob_YHYCoeff(j)*Skl
            if (AreDerivativesNeeded) then
              DActiveSum(1:npt2)=DActiveSum(1:npt2)+ &
                Glob_YHYCoeff(j)*DActive(1:npt2)
              if (OtherDerivativeNeeded) then
                DOtherSum(1:npt2)=DOtherSum(1:npt2)+ &
                  Glob_YHYCoeff(j)*DOther(1:npt2)
              endif
            endif
          endif
        enddo
        Q_Workspace%TrialH(i,a)=Hsum
        Q_Workspace%TrialS(i,a)=Ssum
        if (AreDerivativesNeeded) &
          Glob_D(1:npt2,a,i)=DActiveSum(1:npt2)
        if ((b>0).and.(b/=a)) then
          Q_Workspace%TrialH(ActiveIndex,b)=Hsum
          Q_Workspace%TrialS(ActiveIndex,b)=Ssum
          if (AreDerivativesNeeded) &
            Glob_D(1:npt2,b,ActiveIndex)=DOtherSum(1:npt2)
        endif
      enddo
    enddo

    !The first n rows and m columns are not contiguous when Capacity is larger
    !than MatrixOrder. Reducing the complete allocated arrays preserves their
    !physical leading dimensions and avoids an incorrectly packed MPI count.
    !PreviousH/S are only reduction receive buffers here; ApplyQTrial replaces
    !them with the actual pre-transaction physical columns before any update.
    NumMatrixEntries=size(Q_Workspace%TrialH)
    call MPI_ALLREDUCE(Q_Workspace%TrialH,Q_Workspace%PreviousH, &
      NumMatrixEntries,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_ALLREDUCE(Q_Workspace%TrialS,Q_Workspace%PreviousS, &
      NumMatrixEntries,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    Q_Workspace%TrialH=Q_Workspace%PreviousH
    Q_Workspace%TrialS=Q_Workspace%PreviousS
    if (AreDerivativesNeeded) then
      !Every derivative element is owned by the process that evaluated its
      !symmetry term. An in-place collective avoids a second replicated
      !2*npt-by-active-by-order tensor solely for the reduction receive side.
      NumDerivativeEntries=size(Glob_D)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Glob_D,NumDerivativeEntries, &
        MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

    !All active raw self-overlaps must be known before any column is
    !normalized, because an active-active element depends on both new norms.
    !The form .not.(x>tiny) also rejects a NaN, for which the comparison is
    !false, without requiring an additional IEEE module dependency.
    do a=1,NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      Q_Workspace%TrialDiagS(a)=Q_Workspace%TrialS(ActiveIndex,a)
      if (.not.(Q_Workspace%TrialDiagS(a)>tiny(ONE))) return
    enddo

    do a=1,NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      ActiveNorm=Q_Workspace%TrialDiagS(a)
      do i=1,MatrixOrder
        b=Q_Workspace%ActivePosition(i)
        if (b>0) then
          OtherNorm=Q_Workspace%TrialDiagS(b)
        else
          OtherNorm=Glob_diagS(i)
        endif
        if (.not.(OtherNorm>tiny(ONE))) return
        Normalization=ONE/sqrt(ActiveNorm*OtherNorm)
        Q_Workspace%TrialH(i,a)=Q_Workspace%TrialH(i,a)*Normalization
        Q_Workspace%TrialS(i,a)=Q_Workspace%TrialS(i,a)*Normalization
        if (AreDerivativesNeeded) then
          if (i==ActiveIndex) then
            !The kernel differentiates one side of an identical bra/ket pair.
            !Both sides contribute equally to the physical diagonal.
            Glob_D(1:npt2,a,i)=TWO*Glob_D(1:npt2,a,i)/ActiveNorm
          else
            !As in StoreHSD, keep derivatives of the raw matrix element scaled
            !by the two basis-function norms. The normalization derivative is
            !subtracted once, explicitly, in the energy-gradient contraction.
            Glob_D(1:npt2,a,i)=Glob_D(1:npt2,a,i)*Normalization
          endif
        endif
      enddo
      !Set the analytically normalized diagonal exactly. This avoids allowing
      !roundoff in Sii/Sii to enter overlap tests or qrlinalg's solve path.
      Q_Workspace%TrialS(ActiveIndex,a)=ONE
    enddo

    Q_Workspace%TrialIsReady=.true.
    Q_Workspace%TrialHasDerivatives=AreDerivativesNeeded
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine AssembleQTrial

  subroutine FactorizeQFresh(ErrorCode)
!Subroutine FactorizeQFresh initializes or refreshes the root-owned qrlinalg
!state from canonical, normalized, unshifted Glob_H and Glob_S. The represented
!shift is fixed to Glob_ApproxEnergy for one BBOP step. Rank zero broadcasts the
!status before any process continues to a collective gradient calculation.
!
!Arguments:
    integer,intent(out) :: ErrorCode
!Local variables:
    integer RootError

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    Q_Workspace%FactorsMatchMatrices=.false.
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (Q_Workspace%MatrixOrder<1) return
    if (Q_Workspace%Capacity<Q_Workspace%MatrixOrder) return
    if (.not.allocated(Glob_H)) return
    if (.not.allocated(Glob_S)) return
    if ((size(Glob_H,1)<Q_Workspace%MatrixOrder).or. &
        (size(Glob_H,2)<Q_Workspace%MatrixOrder)) return
    if ((size(Glob_S,1)<Q_Workspace%MatrixOrder).or. &
        (size(Glob_S,2)<Q_Workspace%MatrixOrder)) return

    RootError=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) then
      !Initialization is separated from fresh factorization in qrlinalg. Reuse
      !the allocated state whenever its capacity is already correct; this
      !preserves the lifetime update counter across ordinary refreshes.
      if (Q_Workspace%Factors%get_capacity()/=Q_Workspace%Capacity) then
        call Q_Workspace%Factors%initialize(Q_Workspace%Capacity,RootError)
      endif
      if (RootError==Q_METHOD_SUCCESS) then
        call Q_Workspace%Factors%factorize_fresh(Glob_H,Glob_S, &
          Glob_ApproxEnergy,RootError,active_order=Q_Workspace%MatrixOrder)
      endif
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    ErrorCode=RootError
    Q_Workspace%FactorsMatchMatrices=(ErrorCode==Q_METHOD_SUCCESS)
    if (ErrorCode==Q_METHOD_SUCCESS) &
      Q_Workspace%FreshFactorizations=Q_Workspace%FreshFactorizations+1

  end subroutine FactorizeQFresh

  subroutine AppendQTrial(ErrorCode)
!Subroutine AppendQTrial commits a staged suffix during BASIS_ENL. The active
!map must be the consecutive canonical suffix BaseOrder+1:MatrixOrder. When an
!accepted prefix exists, qrlinalg append_symmetric grows its factors one row
!and column at a time. Each physical column is committed only after the
!matching structural update succeeds.
!
!For BaseOrder=0 there is no valid factorization that can be appended to. The
!routine first commits the complete staged block and then constructs the first
!fresh factorization. This special case is confined to the first enlargement
!step and therefore does not affect the asymptotic candidate-selection cost.
!
!On any append failure, the represented prefix is reconstructed from its
!untouched leading physical block. The caller receives the original append
!status when recovery succeeds, and a recovery status otherwise.
!
!Arguments:
    integer,intent(out) :: ErrorCode
!Local variables:
    integer a,BaseOrder,FunctionIndex,MatrixOrder,NumActive
    integer RootError,StoreError,RecoveryError

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    NumActive=Q_Workspace%NumActive
    BaseOrder=MatrixOrder-NumActive
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (.not.Q_Workspace%TrialIsReady) return
    if ((MatrixOrder<1).or.(NumActive<1).or.(BaseOrder<0)) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%TrialH)) return
    if (.not.allocated(Q_Workspace%TrialS)) return
    if (.not.allocated(Q_Workspace%TrialDiagS)) return
    do a=1,NumActive
      if (Q_Workspace%ActiveFunction(a)/=BaseOrder+a) return
    enddo

    RootError=Q_METHOD_SUCCESS
    if (BaseOrder==0) then
      !Every pair in the new leading block touches an active function, so the
      !staged columns together contain a complete symmetric problem.
      do a=1,NumActive
        FunctionIndex=a
        call StoreQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
          Q_Workspace%TrialH(1:MatrixOrder,a),StoreError)
        if (StoreError/=Q_METHOD_SUCCESS) then
          ErrorCode=StoreError
          return
        endif
        call StoreQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
          Q_Workspace%TrialS(1:MatrixOrder,a),StoreError)
        if (StoreError/=Q_METHOD_SUCCESS) then
          ErrorCode=StoreError
          return
        endif
        Glob_diagS(FunctionIndex)=Q_Workspace%TrialDiagS(a)
      enddo
      !An empty input basis carries a large placeholder CURRENT_ENERGY, not a
      !meaningful inverse-iteration target. Replace such a remote shift by the
      !lowest normalized diagonal estimate before constructing the first QR
      !factorization. For the usual Kstep=1 start this is the exact energy.
      if (abs(Glob_ApproxEnergy)> &
          1000000*max(ONE,maxval(abs([(Glob_H(a,a),a=1,MatrixOrder)])))) then
        Glob_ApproxEnergy=minval([(Glob_H(a,a),a=1,MatrixOrder)])* &
          Glob_InvItParameter
      endif
      call FactorizeQFresh(RootError)
    else
      if (Glob_ProcID==0) then
        if (.not.Q_Workspace%Factors%is_valid()) then
          RootError=Q_METHOD_INVALID_STATE
        else if (Q_Workspace%Factors%order()/=BaseOrder) then
          RootError=Q_METHOD_INVALID_STATE
        endif
      endif
      call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (RootError==Q_METHOD_SUCCESS) then
        do a=1,NumActive
          FunctionIndex=BaseOrder+a
          if (Glob_ProcID==0) then
            call Q_Workspace%Factors%append_symmetric( &
              Q_Workspace%TrialH(1:FunctionIndex,a), &
              Q_Workspace%TrialS(1:FunctionIndex,a),RootError)
          endif
          call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (RootError/=Q_METHOD_SUCCESS) exit

          call StoreQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
            Q_Workspace%TrialH(1:MatrixOrder,a),StoreError)
          call StoreQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
            Q_Workspace%TrialS(1:MatrixOrder,a),StoreError)
          Glob_diagS(FunctionIndex)=Q_Workspace%TrialDiagS(a)
        enddo
      endif
    endif

    if (RootError/=Q_METHOD_SUCCESS) then
      !Only suffix rows and columns may have been committed. Re-expose the
      !accepted leading block and rebuild its factors instead of retaining a
      !partially appended generation.
      Q_Workspace%MatrixOrder=BaseOrder
      Q_Workspace%FactorsMatchMatrices=.false.
      if (BaseOrder>0) then
        call FactorizeQFresh(RecoveryError)
      else
        RecoveryError=Q_METHOD_SUCCESS
        if (Glob_ProcID==0) then
          call Q_Workspace%Factors%clear()
          call Q_Workspace%Factors%initialize(Q_Workspace%Capacity,RecoveryError)
        endif
        call MPI_BCAST(RecoveryError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      endif
      if (RecoveryError/=Q_METHOD_SUCCESS) then
        ErrorCode=RecoveryError
      else
        ErrorCode=RootError
      endif
      Q_Workspace%TrialIsReady=.false.
      Q_Workspace%TrialHasDerivatives=.false.
      return
    endif

    do a=1,NumActive
      FunctionIndex=Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt,a)= &
        Glob_NonlinParam(1:Glob_npt,FunctionIndex)
    enddo
    Q_Workspace%FactorsMatchMatrices=.true.
    Q_Workspace%MatrixParametersAreStored=.true.
    Q_Workspace%TrialIsReady=.false.
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine AppendQTrial

  subroutine EvaluateQAppendedTrial(BaseOrder,TargetOrder,Evalue,ErrorCode)
!Subroutine EvaluateQAppendedTrial performs one complete BASIS_ENL candidate
!transaction. Any preceding candidate suffix is deleted, the new consecutive
!suffix is assembled in canonical order, appended to the accepted QR prefix,
!and solved. On success the active map remains installed so the selected
!candidate can immediately enter ordinary replacement-based optimization.
!
!Arguments:
    integer,intent(in)  :: BaseOrder,TargetOrder
    real(wp),intent(out) :: Evalue
    integer,intent(out) :: ErrorCode
!Local variables:
    integer a,NumActive
    integer ActiveFunction(Q_Workspace%MaxActive)

    Evalue=huge(Evalue)
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    NumActive=TargetOrder-BaseOrder
    if ((BaseOrder<0).or.(TargetOrder>Q_Workspace%Capacity)) return
    if ((NumActive<1).or.(NumActive>Q_Workspace%MaxActive)) return

    call TrimQFactors(BaseOrder,ErrorCode)
    if (ErrorCode/=Q_METHOD_SUCCESS) return

    !The factors still represent BaseOrder while the assembly routines need to
    !see the complete target problem. FactorsMatchMatrices remains false until
    !AppendQTrial has installed every staged suffix column.
    Q_Workspace%MatrixOrder=TargetOrder
    Q_Workspace%FactorsMatchMatrices=.false.
    do a=1,NumActive
      ActiveFunction(a)=BaseOrder+a
    enddo
    call SetQActiveFunctions(ActiveFunction(1:NumActive),ErrorCode)
    if (ErrorCode/=Q_METHOD_SUCCESS) return
    call AssembleQTrial(.false.,ErrorCode)
    if (ErrorCode/=Q_METHOD_SUCCESS) return
    call AppendQTrial(ErrorCode)
    if (ErrorCode/=Q_METHOD_SUCCESS) return
    call SolveQ(Evalue,ErrorCode)

  end subroutine EvaluateQAppendedTrial

  subroutine ApplyQTrial(ErrorCode)
!Subroutine ApplyQTrial updates active columns in the exact order stored in
!Q_Workspace%ActiveFunction. For each column it gathers the currently
!represented physical column, subtracts it from the staged target, calls
!replace_symmetric on rank zero, broadcasts the status, and only then commits
!the physical column on every rank. A preflight check must make all ordinary
!argument failures impossible before the first factor is changed.
!
!Arguments:
    integer,intent(out) :: ErrorCode
!Local variables:
    integer a,FunctionIndex,MatrixOrder,NumActive,RootError,StoreError
    integer RecoveryError

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    NumActive=Q_Workspace%NumActive
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (.not.Q_Workspace%FactorsMatchMatrices) return
    if (.not.Q_Workspace%MatrixParametersAreStored) return
    if (.not.Q_Workspace%TrialIsReady) return
    if ((MatrixOrder<1).or.(NumActive<1)) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%PreviousH)) return
    if (.not.allocated(Q_Workspace%PreviousS)) return
    if (.not.allocated(Q_Workspace%PreviousDiagS)) return
    if (.not.allocated(Q_Workspace%PreviousParam)) return
    if (.not.allocated(Q_Workspace%DeltaH)) return
    if (.not.allocated(Q_Workspace%DeltaS)) return
    if (size(Q_Workspace%DeltaH)<MatrixOrder) return
    if (size(Q_Workspace%DeltaS)<MatrixOrder) return

    !Preflight the private factor metadata on rank zero. Every subsequent
    !replace_symmetric call then has a valid index and vectors of exact order;
    !the library documents that these validated updates have no numerical
    !failure return.
    RootError=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) then
      if (.not.Q_Workspace%Factors%is_valid()) RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%order()/=MatrixOrder) RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%get_capacity()/=Q_Workspace%Capacity) &
        RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%get_shift()/=Glob_ApproxEnergy) &
        RootError=Q_METHOD_INVALID_STATE
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (RootError/=Q_METHOD_SUCCESS) then
      ErrorCode=RootError
      Q_Workspace%FactorsMatchMatrices=.false.
      return
    endif

    !Snapshot every original column before the first commit. In particular,
    !an active-active intersection must be captured before either endpoint is
    !changed. These copies also define the nonlinear-parameter generation to
    !which a rare transaction recovery must return.
    Q_Workspace%PreviousParam=Q_Workspace%MatrixParam
    do a=1,NumActive
      FunctionIndex=Q_Workspace%ActiveFunction(a)
      call GatherQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
        Q_Workspace%PreviousH(1:MatrixOrder,a),StoreError)
      if (StoreError/=Q_METHOD_SUCCESS) then
        ErrorCode=StoreError
        return
      endif
      call GatherQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
        Q_Workspace%PreviousS(1:MatrixOrder,a),StoreError)
      if (StoreError/=Q_METHOD_SUCCESS) then
        ErrorCode=StoreError
        return
      endif
      Q_Workspace%PreviousDiagS(a)=Glob_diagS(FunctionIndex)
    enddo

    do a=1,NumActive
      FunctionIndex=Q_Workspace%ActiveFunction(a)
      !Gather the progressively updated column, not the original snapshot.
      !Earlier active columns have already installed their shared intersection,
      !so the later replacement sees a zero delta at that location.
      call GatherQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
        Q_Workspace%DeltaH(1:MatrixOrder),StoreError)
      call GatherQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
        Q_Workspace%DeltaS(1:MatrixOrder),StoreError)
      Q_Workspace%DeltaH(1:MatrixOrder)= &
        Q_Workspace%TrialH(1:MatrixOrder,a)-Q_Workspace%DeltaH(1:MatrixOrder)
      Q_Workspace%DeltaS(1:MatrixOrder)= &
        Q_Workspace%TrialS(1:MatrixOrder,a)-Q_Workspace%DeltaS(1:MatrixOrder)

      RootError=Q_METHOD_SUCCESS
      if (Glob_ProcID==0) then
        call Q_Workspace%Factors%replace_symmetric(FunctionIndex, &
          Q_Workspace%DeltaH(1:MatrixOrder), &
          Q_Workspace%DeltaS(1:MatrixOrder),RootError)
      endif
      call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (RootError/=Q_METHOD_SUCCESS) exit

      call StoreQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
        Q_Workspace%TrialH(1:MatrixOrder,a),StoreError)
      call StoreQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
        Q_Workspace%TrialS(1:MatrixOrder,a),StoreError)
      Glob_diagS(FunctionIndex)=Q_Workspace%TrialDiagS(a)
    enddo

    if (RootError/=Q_METHOD_SUCCESS) then
      !A library failure is not expected after the preflight, but fail safely:
      !restore the complete physical generation and construct fresh factors
      !instead of trying to reason about a possibly partial QR update.
      do a=1,NumActive
        FunctionIndex=Q_Workspace%ActiveFunction(a)
        call StoreQCanonicalColumn(Glob_H,MatrixOrder,FunctionIndex, &
          Q_Workspace%PreviousH(1:MatrixOrder,a),StoreError)
        call StoreQCanonicalColumn(Glob_S,MatrixOrder,FunctionIndex, &
          Q_Workspace%PreviousS(1:MatrixOrder,a),StoreError)
        Glob_diagS(FunctionIndex)=Q_Workspace%PreviousDiagS(a)
        Glob_NonlinParam(1:Glob_npt,FunctionIndex)= &
          Q_Workspace%PreviousParam(1:Glob_npt,a)
      enddo
      Q_Workspace%MatrixParam=Q_Workspace%PreviousParam
      call FactorizeQFresh(RecoveryError)
      if (RecoveryError/=Q_METHOD_SUCCESS) then
        ErrorCode=RecoveryError
      else
        ErrorCode=RootError
      endif
      Q_Workspace%TrialIsReady=.false.
      Q_Workspace%TrialHasDerivatives=.false.
      return
    endif

    do a=1,NumActive
      FunctionIndex=Q_Workspace%ActiveFunction(a)
      Q_Workspace%MatrixParam(1:Glob_npt,a)= &
        Glob_NonlinParam(1:Glob_npt,FunctionIndex)
    enddo
    Q_Workspace%FactorsMatchMatrices=.true.
    Q_Workspace%MatrixParametersAreStored=.true.
    Q_Workspace%TrialIsReady=.false.
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine ApplyQTrial

  subroutine ComputeQEigenpairResidual(Evalue,Eigenvector,AbsoluteResidual, &
                                       RelativeResidual,ErrorCode)
!Subroutine ComputeQEigenpairResidual measures the residual of the physical
!generalized eigenproblem, independently of qrlinalg's direction-change test:
!
!                  r = H*c-Evalue*S*c .
!
!The returned relative value is ||r||/(||H*c||+|E|*||S*c||+tiny). Only the
!canonical lower triangles are read. This routine is serial and is called on
!rank zero; it performs no allocation and does not modify H, S, or c.
!
!Arguments:
    real(wp),intent(in)  :: Evalue,Eigenvector(:)
    real(wp),intent(out) :: AbsoluteResidual,RelativeResidual
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer MatrixOrder
    real(wp) HNorm2,SNorm2,ResidualNorm2

    AbsoluteResidual=ZERO
    RelativeResidual=ZERO
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    if (MatrixOrder<1) return
    if (size(Eigenvector)/=MatrixOrder) return
    if (.not.allocated(Q_Workspace%DeltaH)) return
    if (.not.allocated(Q_Workspace%DeltaS)) return

    !The earlier scalar double loop performed the same two symmetric products
    !but left optimized BLAS performance unused. DSYMV reads exactly the
    !authoritative lower triangles, handles the physical leading dimension,
    !and is substantially faster for the residual that follows every Q solve.
    !This helper is called only on rank zero, so calling BLAS directly is also
    !important: the MPI-routing MTMVL wrapper is collective for some calibrated
    !matrix sizes and therefore cannot be entered by the root alone.
    call DSYMV('L',MatrixOrder,ONE,Glob_H,size(Glob_H,1),Eigenvector,1, &
      ZERO,Q_Workspace%DeltaH,1)
    call DSYMV('L',MatrixOrder,ONE,Glob_S,size(Glob_S,1),Eigenvector,1, &
      ZERO,Q_Workspace%DeltaS,1)
    HNorm2=dot_product(Q_Workspace%DeltaH(1:MatrixOrder), &
                      Q_Workspace%DeltaH(1:MatrixOrder))
    SNorm2=dot_product(Q_Workspace%DeltaS(1:MatrixOrder), &
                      Q_Workspace%DeltaS(1:MatrixOrder))
    Q_Workspace%DeltaH(1:MatrixOrder)=Q_Workspace%DeltaH(1:MatrixOrder)- &
      Evalue*Q_Workspace%DeltaS(1:MatrixOrder)
    ResidualNorm2=dot_product(Q_Workspace%DeltaH(1:MatrixOrder), &
                             Q_Workspace%DeltaH(1:MatrixOrder))
    AbsoluteResidual=sqrt(ResidualNorm2)
    RelativeResidual=AbsoluteResidual/ &
      (sqrt(HNorm2)+abs(Evalue)*sqrt(SNorm2)+tiny(ONE))
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine ComputeQEigenpairResidual

  subroutine SolveQ(Evalue,ErrorCode)
!Subroutine SolveQ calls qrlinalg inverse iteration on rank zero with
!distinct input and output vectors and normalization mode zero, then broadcast
!the physical energy, S-normalized Glob_c, convergence diagnostics, and status.
!A QR_ERR_NO_CONVERGENCE result contains a usable approximation. It is accepted
!only when the independent physical generalized-eigenpair residual satisfies
!the requested tolerance; otherwise one fresh-factorization retry is made.
!
!Arguments:
    real(wp),intent(out) :: Evalue
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer i,MatrixOrder,NumOfIterations,RootError,RefreshError,SolveStatus
    integer RetryMaxIterations
    integer(int64) UpdatesSinceFresh,RefreshLimit
    real(wp) AbsoluteResidual,RelativeResidual,FactorResidual
    real(wp) RefreshTolerance,SolveTolerance
    logical RefreshNeeded

    Evalue=huge(Evalue)
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (.not.Q_Workspace%FactorsMatchMatrices) return
    if (MatrixOrder<1) return
    if (.not.allocated(Glob_c)) return
    if (size(Glob_c)<MatrixOrder) return
    if (.not.allocated(Q_Workspace%InitialVector)) return
    if (.not.allocated(Q_Workspace%SolvedVector)) return

    !A one-function generalized problem has an exact closed-form solution.
    !Besides avoiding unnecessary inverse iteration, this is essential when a
    !calculation starts from the conventional enormous CURRENT_ENERGY sentinel:
    !forming lambda as shift+(lambda-shift) would otherwise lose the physical
    !diagonal energy by catastrophic cancellation. BASIS_ENL retargets the QR
    !shift to this first accepted energy before attempting order two.
    if (MatrixOrder==1) then
      if (.not.(Glob_S(1,1)>tiny(ONE))) then
        ErrorCode=QR_ERR_NONPOSITIVE_OVERLAP
        return
      endif
      Evalue=Glob_H(1,1)/Glob_S(1,1)
      Glob_c(1)=ONE/sqrt(Glob_S(1,1))
      Q_Workspace%LastEigenpairResidual=ZERO
      Q_Workspace%LastFactorResidual=ZERO
      Glob_LastEigvalTol=ZERO
      Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2=Glob_InvItTempCounter2+1
      ErrorCode=Q_METHOD_SUCCESS
      return
    endif

    !A full refresh after O(n) replacements keeps the amortized cost of fresh
    !O(n**3) factorizations at O(n**2) per replacement. The residual check
    !below can request an earlier refresh when accumulated rotations drift.
    RefreshNeeded=.false.
    if (Glob_ProcID==0) then
      if (.not.Q_Workspace%Factors%is_valid()) then
        RootError=Q_METHOD_INVALID_STATE
      else
        RootError=Q_METHOD_SUCCESS
        UpdatesSinceFresh=Q_Workspace%Factors%get_updates_since_fresh()
        RefreshLimit=int(max(64,8*MatrixOrder),int64)
        RefreshNeeded=(UpdatesSinceFresh>=RefreshLimit)
      endif
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(RefreshNeeded,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (RootError/=Q_METHOD_SUCCESS) then
      ErrorCode=RootError
      Q_Workspace%FactorsMatchMatrices=.false.
      return
    endif
    if (RefreshNeeded) then
      call FactorizeQFresh(RefreshError)
      if (RefreshError/=Q_METHOD_SUCCESS) then
        ErrorCode=RefreshError
        return
      endif
    endif

    RootError=Q_METHOD_SUCCESS
    AbsoluteResidual=ZERO
    RelativeResidual=ZERO
    NumOfIterations=0
    if (Glob_ProcID==0) then
      Q_Workspace%InitialVector(1:MatrixOrder)=Glob_c(1:MatrixOrder)
      if (.not.(maxval(abs(Q_Workspace%InitialVector(1:MatrixOrder)))>tiny(ONE))) then
        do i=1,MatrixOrder
          Q_Workspace%InitialVector(i)=ONE
        enddo
      endif
      call Q_Workspace%Factors%solve(Glob_S, &
        Q_Workspace%InitialVector(1:MatrixOrder), &
        Q_Workspace%SolvedVector(1:MatrixOrder),Evalue,Glob_EigvalTol, &
        Glob_MaxIterForGSEPIIS,0,RelativeResidual,NumOfIterations,RootError)
      Glob_c(1:MatrixOrder)=Q_Workspace%SolvedVector(1:MatrixOrder)
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(RelativeResidual,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(NumOfIterations,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_c,MatrixOrder,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    SolveStatus=RootError
    Glob_LastEigvalTol=RelativeResidual
    Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
    Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
    if ((SolveStatus/=Q_METHOD_SUCCESS).and. &
        (SolveStatus/=QR_ERR_NO_CONVERGENCE)) then
      ErrorCode=RootError
      return
    endif

    !Measure both independent numerical relationships. The physical residual
    !checks H*c=E*S*c. The factor residual checks Q*R=(H-shift*S) with a fixed
    !dense probe. The eigenvector must not be used for the latter because it
    !is nearly a null vector of the deliberately near-eigenvalue shifted
    !matrix, which amplifies harmless factorization roundoff in the relative
    !action ratio. QR_ERR_NO_CONVERGENCE still contains an approximation; it is
    !accepted only when this independent physical residual satisfies the
    !caller tolerance. Either diagnostic may request one fresh-factorization
    !retry.
    if (Glob_ProcID==0) then
      RootError=Q_METHOD_SUCCESS
      call ComputeQEigenpairResidual(Evalue,Glob_c(1:MatrixOrder), &
        AbsoluteResidual,Q_Workspace%LastEigenpairResidual,RootError)
      if (RootError==Q_METHOD_SUCCESS) then
        do i=1,MatrixOrder
          Q_Workspace%InitialVector(i)=ONE+ &
            real(mod(17*i,23),wp)/real(23,wp)
          if (mod(i,2)==0) Q_Workspace%InitialVector(i)= &
            -Q_Workspace%InitialVector(i)
        enddo
        call Q_Workspace%Factors%factorization_residual(Glob_H,Glob_S, &
          Q_Workspace%InitialVector(1:MatrixOrder),AbsoluteResidual, &
          FactorResidual,RootError)
        Q_Workspace%LastFactorResidual=FactorResidual
      endif
      RefreshTolerance=max(1000*epsilon(ONE)*MatrixOrder, &
                           10*abs(Glob_EigvalTol))
      SolveTolerance=max(1000*epsilon(ONE)*MatrixOrder, &
                         100*abs(Glob_EigvalTol))
      RefreshNeeded=(RootError==Q_METHOD_SUCCESS).and. &
        ((.not.(Q_Workspace%LastFactorResidual<=RefreshTolerance)).or. &
         (.not.(Q_Workspace%LastEigenpairResidual<=SolveTolerance)))
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(RefreshNeeded,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Q_Workspace%LastEigenpairResidual,1,MPI_WP,0, &
      MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Q_Workspace%LastFactorResidual,1,MPI_WP,0, &
      MPI_COMM_WORLD,Glob_MPIErrCode)
    if (RootError/=Q_METHOD_SUCCESS) then
      ErrorCode=RootError
      return
    endif

    if (RefreshNeeded) then
      call FactorizeQFresh(RefreshError)
      if (RefreshError/=Q_METHOD_SUCCESS) then
        ErrorCode=RefreshError
        return
      endif
      !The normal iteration cap is tuned for small optimizer displacements.
      !A newly appended basis function can rotate the eigenvector much farther.
      !Only after the independent residual has rejected the first approximation
      !do we allow this larger cap; every additional iteration remains O(n**2).
      RetryMaxIterations=max(120,max(Glob_MaxIterForGSEPIIS,4*MatrixOrder))
      if (Glob_ProcID==0) then
        Q_Workspace%InitialVector(1:MatrixOrder)=Glob_c(1:MatrixOrder)
        call Q_Workspace%Factors%solve(Glob_S, &
          Q_Workspace%InitialVector(1:MatrixOrder), &
          Q_Workspace%SolvedVector(1:MatrixOrder),Evalue,Glob_EigvalTol, &
          RetryMaxIterations,0,RelativeResidual,NumOfIterations,RootError)
        Glob_c(1:MatrixOrder)=Q_Workspace%SolvedVector(1:MatrixOrder)
      endif
      call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(RelativeResidual,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(NumOfIterations,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_c,MatrixOrder,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      Glob_LastEigvalTol=RelativeResidual
      Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations

      if ((RootError==Q_METHOD_SUCCESS).or. &
          (RootError==QR_ERR_NO_CONVERGENCE)) then
        if (Glob_ProcID==0) then
          RootError=Q_METHOD_SUCCESS
          call ComputeQEigenpairResidual(Evalue,Glob_c(1:MatrixOrder), &
            AbsoluteResidual,Q_Workspace%LastEigenpairResidual,RootError)
          if (RootError==Q_METHOD_SUCCESS) then
            do i=1,MatrixOrder
              Q_Workspace%InitialVector(i)=ONE+ &
                real(mod(17*i,23),wp)/real(23,wp)
              if (mod(i,2)==0) Q_Workspace%InitialVector(i)= &
                -Q_Workspace%InitialVector(i)
            enddo
            call Q_Workspace%Factors%factorization_residual(Glob_H,Glob_S, &
              Q_Workspace%InitialVector(1:MatrixOrder),AbsoluteResidual, &
              FactorResidual,RootError)
            Q_Workspace%LastFactorResidual=FactorResidual
          endif
          if ((RootError==Q_METHOD_SUCCESS).and. &
              (.not.(Q_Workspace%LastEigenpairResidual<=SolveTolerance))) &
            RootError=QR_ERR_NO_CONVERGENCE
          if ((RootError==Q_METHOD_SUCCESS).and. &
              (.not.(Q_Workspace%LastFactorResidual<=RefreshTolerance))) &
            RootError=Q_METHOD_INVALID_STATE
        endif
        call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(Q_Workspace%LastEigenpairResidual,1,MPI_WP,0, &
          MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(Q_Workspace%LastFactorResidual,1,MPI_WP,0, &
          MPI_COMM_WORLD,Glob_MPIErrCode)
      endif
    endif
    ErrorCode=RootError

  end subroutine SolveQ

  subroutine ComputeQOverlapPenalty(TotalPenalty,ErrorCode)
!Subroutine ComputeQOverlapPenalty evaluates the FULL_OPT1 smooth pair-overlap
!penalty without assuming that optimized functions form a trailing block. Each
!unordered pair touching the explicit active map is visited exactly once and
!assigned to one MPI rank. Glob_S is read only through its canonical lower
!triangle.
!
!Arguments:
    real(wp),intent(out) :: TotalPenalty
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer a,b,i,PairNumber,ActiveIndex
    real(wp) PairOverlap,PairPenalty,LocalPenalty,PenaltyCoefficient

    TotalPenalty=ZERO
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (Q_Workspace%NumActive<1) return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (.not.allocated(Q_Workspace%ActiveFunction)) return
    if (.not.allocated(Q_Workspace%ActivePosition)) return
    if (.not.(Glob_OverlapPenaltyThreshold2<ONE)) return

    PenaltyCoefficient=Glob_MaxOverlapPenalty/ &
      (ONE-Glob_OverlapPenaltyThreshold2)
    LocalPenalty=ZERO
    PairNumber=0
    do a=1,Q_Workspace%NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      do i=1,Q_Workspace%MatrixOrder
        if (i==ActiveIndex) cycle
        b=Q_Workspace%ActivePosition(i)
        if ((b>0).and.(b<a)) cycle
        PairNumber=PairNumber+1
        if (mod(PairNumber-1,Glob_NumOfProcs)/=Glob_ProcID) cycle
        PairOverlap=QCanonicalMatrixElement(Glob_S,ActiveIndex,i)
        if (PairOverlap*PairOverlap>Glob_OverlapPenaltyThreshold2) then
          PairPenalty=PenaltyCoefficient* &
            (PairOverlap*PairOverlap-Glob_OverlapPenaltyThreshold2)
          LocalPenalty=LocalPenalty+PairPenalty
        endif
      enddo
    enddo
    call MPI_ALLREDUCE(LocalPenalty,TotalPenalty,1,MPI_WP,MPI_SUM, &
      MPI_COMM_WORLD,Glob_MPIErrCode)
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine ComputeQOverlapPenalty

  subroutine ComputeQOverlapPenaltyAndAddGradient(TotalPenalty,WkGR,ErrorCode)
!Subroutine ComputeQOverlapPenaltyAndAddGradient evaluates the same mapped
!penalty as ComputeQOverlapPenalty and adds its analytic derivative to WkGR.
!The derivative tensor produced by AssembleQTrial stores
!
!  d<Sraw_ij>/sqrt(Sraw_ii*Sraw_jj)
!
!for an active endpoint, while its active diagonal stores
!dSraw_ii/Sraw_ii. Consequently the derivative of a normalized overlap is the
!cross derivative minus one half of the overlap times the diagonal derivative.
!For an active-active pair this expression is applied independently to both
!endpoints. Pair ownership follows MPI rank, so the caller's subsequent
!gradient MPI_ALLREDUCE combines both the energy-gradient and penalty pieces.
!
!Arguments:
    real(wp),intent(out)   :: TotalPenalty
    real(wp),intent(inout) :: WkGR(:)
    integer,intent(out)    :: ErrorCode
!Local variables:
    integer a,b,i,m,PairNumber,ActiveIndex
    real(wp) PairOverlap,LocalPenalty,PenaltyCoefficient
    real(wp) GradientCoefficient,OverlapDerivative

    TotalPenalty=ZERO
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (Q_Workspace%NumActive<1) return
    if (.not.Q_Workspace%TrialHasDerivatives) return
    if (size(WkGR)<Q_Workspace%NumActive*Glob_npt) return
    if (.not.allocated(Glob_D)) return
    if (.not.(Glob_OverlapPenaltyThreshold2<ONE)) return

    PenaltyCoefficient=Glob_MaxOverlapPenalty/ &
      (ONE-Glob_OverlapPenaltyThreshold2)
    LocalPenalty=ZERO
    PairNumber=0
    do a=1,Q_Workspace%NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      do i=1,Q_Workspace%MatrixOrder
        if (i==ActiveIndex) cycle
        b=Q_Workspace%ActivePosition(i)
        if ((b>0).and.(b<a)) cycle
        PairNumber=PairNumber+1
        if (mod(PairNumber-1,Glob_NumOfProcs)/=Glob_ProcID) cycle
        PairOverlap=QCanonicalMatrixElement(Glob_S,ActiveIndex,i)
        if (PairOverlap*PairOverlap>Glob_OverlapPenaltyThreshold2) then
          LocalPenalty=LocalPenalty+PenaltyCoefficient* &
            (PairOverlap*PairOverlap-Glob_OverlapPenaltyThreshold2)
          GradientCoefficient=TWO*PenaltyCoefficient*PairOverlap
          do m=1,Glob_npt
            OverlapDerivative=Glob_D(Glob_npt+m,a,i)-ONEHALF* &
              PairOverlap*Glob_D(Glob_npt+m,a,ActiveIndex)
            WkGR((a-1)*Glob_npt+m)=WkGR((a-1)*Glob_npt+m)+ &
              GradientCoefficient*OverlapDerivative
          enddo
          if (b>0) then
            do m=1,Glob_npt
              OverlapDerivative=Glob_D(Glob_npt+m,b,ActiveIndex)-ONEHALF* &
                PairOverlap*Glob_D(Glob_npt+m,b,i)
              WkGR((b-1)*Glob_npt+m)=WkGR((b-1)*Glob_npt+m)+ &
                GradientCoefficient*OverlapDerivative
            enddo
          endif
        endif
      enddo
    enddo
    call MPI_ALLREDUCE(LocalPenalty,TotalPenalty,1,MPI_WP,MPI_SUM, &
      MPI_COMM_WORLD,Glob_MPIErrCode)
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine ComputeQOverlapPenaltyAndAddGradient

  subroutine GetQOverlapStatistics(MaxAbsOverlap,MinAbsOverlap, &
                                   AverageAbsOverlap,ErrorCode)
!Subroutine GetQOverlapStatistics reports statistics for all unordered pairs
!touching the active map. It is the canonical-index counterpart of
!GetOverlapStatistics, whose Nmin:Nmax interface assumes a trailing block.
!
!Arguments:
    real(wp),intent(out) :: MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer a,b,i,NumPairs,ActiveIndex
    real(wp) PairOverlap,AbsPairOverlap

    MaxAbsOverlap=ZERO
    MinAbsOverlap=huge(ONE)
    AverageAbsOverlap=ZERO
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (Q_Workspace%NumActive<1) return

    NumPairs=0
    do a=1,Q_Workspace%NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      do i=1,Q_Workspace%MatrixOrder
        if (i==ActiveIndex) cycle
        b=Q_Workspace%ActivePosition(i)
        if ((b>0).and.(b<a)) cycle
        PairOverlap=QCanonicalMatrixElement(Glob_S,ActiveIndex,i)
        AbsPairOverlap=abs(PairOverlap)
        if (AbsPairOverlap>abs(MaxAbsOverlap)) MaxAbsOverlap=PairOverlap
        if (AbsPairOverlap<abs(MinAbsOverlap)) MinAbsOverlap=PairOverlap
        AverageAbsOverlap=AverageAbsOverlap+AbsPairOverlap
        NumPairs=NumPairs+1
      enddo
    enddo
    if (NumPairs>0) then
      AverageAbsOverlap=AverageAbsOverlap/NumPairs
    else
      MinAbsOverlap=ZERO
    endif
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine GetQOverlapStatistics

  function EnergyQA(AreMatElemNeeded,ErrorCode)
!Function EnergyQA provides the Q counterpart of EnergyGA. The active set
!is taken from Q_Workspace rather than encoded as a trailing Nmin:Nmax range.
!It assembles and applies a trial when AreMatElemNeeded is true, solves the
!represented problem, and returns the requested physical eigenvalue.
!
!Arguments:
    logical,intent(in)  :: AreMatElemNeeded
    integer,intent(out) :: ErrorCode
    real(wp) EnergyQA

    EnergyQA=huge(EnergyQA)
    ErrorCode=Q_METHOD_SUCCESS
    if (AreMatElemNeeded) then
      call AssembleQTrial(.false.,ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) return
      call ApplyQTrial(ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) return
    endif
    call SolveQ(EnergyQA,ErrorCode)

    if ((ErrorCode==Q_METHOD_SUCCESS).and.Glob_OverlapPenaltyAllowed) then
      call ComputeQOverlapPenalty(Glob_TotalOverlapPenalty,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) &
        EnergyQA=EnergyQA+Glob_TotalOverlapPenalty
    endif
    Glob_EnergyGACounter=Glob_EnergyGACounter+1

  end function EnergyQA

  function EnergyQAM(AreMatElemNeeded,ErrorCode)
!Function EnergyQAM provides the Q counterpart of EnergyGAM. qrlinalg
!always computes an eigenvector, so EnergyQA and EnergyQAM may share one solve;
!the separate entry point keeps G's candidate-acceptance call structure clear.
!
!Arguments:
    logical,intent(in)  :: AreMatElemNeeded
    integer,intent(out) :: ErrorCode
    real(wp) EnergyQAM

    EnergyQAM=EnergyQA(AreMatElemNeeded,ErrorCode)

  end function EnergyQAM

  subroutine EnergyQB(Evalue,Gradient,AreMatElemNeeded,ErrorCode)
!Subroutine EnergyQB provides the Q counterpart of EnergyGB. Gradient block
!a corresponds to Q_Workspace%ActiveFunction(a); contractions must use the
!coefficient of that canonical function instead of Glob_c(a+Glob_nfru).
!
!Arguments:
    real(wp),intent(out) :: Evalue
    real(wp),intent(out) :: Gradient(:)
    logical,intent(in)   :: AreMatElemNeeded
    integer,intent(out)  :: ErrorCode
!Local variables:
    integer a,l,m,npt,MatrixOrder,NumActive,ActiveIndex
    real(wp) W(Glob_npt_MaxAllowed),t,t2

    Evalue=huge(Evalue)
    Gradient=huge(Evalue)
    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    MatrixOrder=Q_Workspace%MatrixOrder
    NumActive=Q_Workspace%NumActive
    npt=Glob_npt
    if (size(Gradient)<NumActive*npt) return
    if (.not.allocated(Glob_D)) return
    if (.not.allocated(Glob_WkGR)) return
    if (size(Glob_WkGR)<NumActive*npt) return

    if (AreMatElemNeeded) then
      call AssembleQTrial(.true.,ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) return
      call ApplyQTrial(ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) return
    else
      if (.not.Q_Workspace%TrialHasDerivatives) return
    endif

    call SolveQ(Evalue,ErrorCode)
    if (ErrorCode/=Q_METHOD_SUCCESS) return

    !The derivative tensor stores the derivative of one raw matrix element,
    !scaled by both raw basis norms. Contracting the full conceptual row gives
    !both symmetric H/S contributions. The final diagonal subtraction removes
    !the duplicate diagonal and, through the eigenvalue equation, accounts for
    !the derivative of normalization of every element touching this function.
    Glob_WkGR(1:NumActive*npt)=ZERO
    do a=1,NumActive
      ActiveIndex=Q_Workspace%ActiveFunction(a)
      W(1:npt)=ZERO
      do l=1+Glob_ProcID,MatrixOrder,Glob_NumOfProcs
        t=Glob_c(l)
        do m=1,npt
          W(m)=W(m)+t*(Glob_D(m,a,l)-Evalue*Glob_D(m+npt,a,l))
        enddo
      enddo
      t=Glob_c(ActiveIndex)
      t2=t*t
      do m=1,npt
        Glob_WkGR((a-1)*npt+m)=TWO*t*W(m)
      enddo
      do m=1+Glob_ProcID,npt,Glob_NumOfProcs
        Glob_WkGR((a-1)*npt+m)=Glob_WkGR((a-1)*npt+m)- &
          t2*(Glob_D(m,a,ActiveIndex)- &
              Evalue*Glob_D(m+npt,a,ActiveIndex))
      enddo
    enddo

    if (Glob_OverlapPenaltyAllowed) then
      call ComputeQOverlapPenaltyAndAddGradient(Glob_TotalOverlapPenalty, &
        Glob_WkGR,ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) then
        Evalue=huge(Evalue)
        Gradient=huge(Evalue)
        return
      endif
      Evalue=Evalue+Glob_TotalOverlapPenalty
    endif
    call MPI_ALLREDUCE(Glob_WkGR,Gradient,NumActive*npt,MPI_WP, &
      MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    Glob_EnergyGBCounter=Glob_EnergyGBCounter+1

  end subroutine EnergyQB

  subroutine BasisEnlQ(Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval, &
                       OverlapThreshold,LinCoeffThreshold,ErrorCode)
!Subroutine BasisEnlQ enlarges the canonical basis using the same random
!candidate generation, optional DRMNG optimization, acceptance thresholds,
!history, and output policy as BasisEnlG. Accepted functions are appended at
!the end. Each rejected trial suffix is removed from the QR factors before the
!next candidate is appended, so repeated trials do not repeatedly increase the
!active matrix order.
!
!Arguments:
    integer,intent(in)  :: Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
    integer,intent(out) :: ErrorCode
!Local variables:
    integer i,j,K,AttemptToGetGoodFunc,ii,jj,jbest,jj2,jbest2,j2
    integer npt,nfo,nfru,nv,nvmax
    integer ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
    integer wbfu_t,wmu_t,wbfu,wmu,rgm1_counter,rgm2_counter
    real(wp) ms1,ms2,Evalue,E_init,E_best,t
    real(wp),allocatable :: ParSet(:,:),ParSetBest(:,:)
    integer,allocatable :: IndSet(:,:),IndSetBest(:,:),IndOptSequence(:)
    real(wp),allocatable :: x(:),x_best(:),grad(:)
    real(wp),allocatable :: D(:),V(:),V_init(:)
    integer,parameter :: LIV=60
    integer IV(LIV),IV_init(LIV),LV,ALG
    logical IsSwapFileOK,IsEnergyImproved,ExitNeeded
    logical IsOverlapBad,IsAnyLinCoeffBad,IsEnergyBad

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    if (Kstart/=Glob_CurrBasisSize+1) return
    if ((Kstart<1).or.(Kstop<Kstart).or.(Kstep<1).or.(NTrials<1)) return
    if ((OptimizationType<0).or.(OptimizationType>1)) return
    if (MaxEnergyEval<0) return

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine BasisEnlQ started'
      write(*,'(1x,a,1x,i0)') 'Kstart =',Kstart
      write(*,'(1x,a,1x,i0)') 'Kstop =',Kstop
      write(*,'(1x,a,1x,i0)') 'Kstep =',Kstep
      write(*,'(1x,a,1x,i0)') 'OptimizationType =',OptimizationType
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval =',MaxEnergyEval
    endif

    Glob_GSEPSolutionMethod='Q'
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Kstop
    Glob_HSBuffLen=Kstop*Kstep
    npt=Glob_npt
    nvmax=Kstep*npt
    wbfu_t=0
    wmu_t=0
    rgm1_counter=0
    rgm2_counter=0
    ms1=ZERO
    ms2=ZERO

!Reallocate basis metadata to final capacity before matrix storage is created.
!The accepted prefix remains in its original canonical order.
    call ReallocateBasisFuncData(Kstop,Glob_CurrBasisSize)

    allocate(Glob_H(Kstop,Kstop))
    allocate(Glob_S(Kstop,Kstop))
    allocate(Glob_diagS(Kstop))
    allocate(Glob_D(2*npt,Kstep,Kstop))
    allocate(Glob_c(Kstop))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_WkGR(nvmax))
    allocate(ParSet(npt,Kstep))
    allocate(ParSetBest(npt,Kstep))
    allocate(IndSet(Kstep,2))
    allocate(IndSetBest(Kstep,2))
    allocate(IndOptSequence(Kstep))
    allocate(x(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2+1
    allocate(V(LV))
    allocate(V_init(LV))

    call PrepareQWorkspace(Glob_CurrBasisSize,Kstop,Kstep,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0125 in BasisEnlQ: Q workspace cannot be allocated'
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif

!Use exactly the same optimizer controls as BasisEnlG so comparisons isolate
!the eigensolver and matrix-layout changes.
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
    IV_init(1)=12

!Only the accepted prefix can be present in the swap file. Initialize unused
!capacity before the full-capacity MPI broadcast performed by the shared swap
!reader, then construct one fresh prefix factorization.
    Glob_H=ZERO
    Glob_S=ZERO
    Glob_diagS=ZERO
    Glob_c=ONE
    call ReadSwapFileAndDistributeData(IsSwapFileOK)
    Q_Workspace%MatricesAreCanonical=.true.
    if (Glob_CurrBasisSize>0) then
      if (.not.IsSwapFileOK) call ComputeMatElem(1,Glob_CurrBasisSize)
      call FactorizeQFresh(ErrCode)
      if (ErrCode==Q_METHOD_SUCCESS) call SolveQ(Glob_CurrEnergy,ErrCode)
    else
      Glob_CurrEnergy=huge(Glob_CurrEnergy)
      call TrimQFactors(0,ErrCode)
    endif
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,'(1x,a,1x,i0)') &
        'Error EC0126 in BasisEnlQ: initial Q state cannot be constructed, status',ErrCode
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif
    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

    K=Kstart-1
    do while (K<Kstop)
      nfru=K
      nfo=min(Kstep,Kstop-K)
      K=K+nfo
      nv=nfo*npt
      Glob_nfa=K
      Glob_nfru=nfru
      Glob_nfo=nfo
      call linalg_setparam(K)
      E_init=Glob_CurrEnergy

      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0)') 'Current basis size is',Glob_CurrBasisSize
        if (nfo>1) then
          write(*,'(1x,a,1x,i0,a,i0)') 'Selecting functions',nfru+1,'-',K
        else
          write(*,'(1x,a,1x,i0)') 'Selecting function',K
        endif
      endif

      IsOverlapBad=.true.
      IsAnyLinCoeffBad=.true.
      IsEnergyBad=.false.
      AttemptToGetGoodFunc=1
      do while ((IsOverlapBad.or.IsAnyLinCoeffBad.or.IsEnergyBad).and. &
                (AttemptToGetGoodFunc<=Glob_BadOverlapOrLinCoeffLim))
        NumOfFailures=0
        IsEnergyImproved=.false.
        wbfu=0
        wmu=0

!Random candidate selection uses append/delete transactions. Glob_CurrEnergy
!tracks the best trial found so far, while ParSetBest owns its parameters.
        do i=1,NTrials
          if (Glob_ProcID==0) call GenerateTrialParam(nfo,ParSet,IndSet,wbfu_t,wmu_t)
          call MPI_BCAST(ParSet,npt*nfo,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_BCAST(IndSet,2*nfo,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          Glob_NonlinParam(1:npt,nfru+1:K)=ParSet(1:npt,1:nfo)
          Glob_Index(nfru+1:K,1)=IndSet(1:nfo,1)
          Glob_Index(nfru+1:K,2)=IndSet(1:nfo,2)
          call EvaluateQAppendedTrial(nfru,K,Evalue,ErrCode)
          if (ErrCode==Q_METHOD_SUCCESS) then
            if (Evalue<Glob_CurrEnergy) then
              Glob_CurrEnergy=Evalue
              ParSetBest(1:npt,1:nfo)=ParSet(1:npt,1:nfo)
              IndSetBest(1:nfo,1:2)=IndSet(1:nfo,1:2)
              IsEnergyImproved=.true.
              wbfu=wbfu_t
              wmu=wmu_t
            endif
          else
            NumOfFailures=NumOfFailures+1
            if (Glob_ProcID==0) write(*,'(1x,a,1x,i0,1x,a,1x,i0)') &
              'Warning WC0114 in BasisEnlQ: candidate',i,'failed with Q status',ErrCode
          endif
        enddo
        if (NumOfFailures*ONE/NTrials>Glob_MaxFracOfTrialFailsAllowed) then
          if (Glob_ProcID==0) then
            write(*,*) 'Error EC0127 in BasisEnlQ: Q candidate solve failures exceeded limit'
            write(*,'(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials,'%'
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif
        if (.not.IsEnergyImproved) then
          if (Glob_ProcID==0) write(*,*) &
            'Error EC0128 in BasisEnlQ: random selection produced no energy improvement'
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif

!Rebuild the winning candidate because the final random trial need not be the
!winner. A successful append captures the matrix parameter generation needed
!by the ordinary EnergyQA/EnergyQB replacement transactions below.
        Glob_NonlinParam(1:npt,nfru+1:K)=ParSetBest(1:npt,1:nfo)
        Glob_Index(nfru+1:K,1)=IndSetBest(1:nfo,1)
        Glob_Index(nfru+1:K,2)=IndSetBest(1:nfo,2)
        call EvaluateQAppendedTrial(nfru,K,Glob_CurrEnergy,ErrCode)
        if (ErrCode/=Q_METHOD_SUCCESS) then
          if (Glob_ProcID==0) write(*,*) &
            'Error EC0129 in BasisEnlQ: selected Q candidate cannot be reconstructed'
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif
        if (Glob_ProcID==0) then
          write(*,'(1x,a)',advance='no') 'E='
          call writereal(6,Glob_CurrEnergy)
          write(*,'(5x,a,1x,i0)') 'prototype function is',wbfu
          do i=1,nfo
            write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') &
              nfru+i,':',IndSetBest(i,1),IndSetBest(i,2)
            call writerealarradv(6,ParSetBest(1:npt,i),npt)
          enddo
          write(*,*) 'Optimizing nonlinear parameters'
        endif

        do i=1,nfo
          x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
        enddo
        x_best(1:nv)=x(1:nv)
        NumOfEnergyEval=0
        NumOfGradEval=0

        select case (OptimizationType)
        case(0)
          !The selected candidate is already the best optimizer point.
        case(1)
          !The two Cartesian premultiplier indices are discrete basis metadata,
          !not DRMNG variables. Preserve the 2D G-method scan, but evaluate
          !each alternative through the Q suffix transaction. The final
          !winning pair is reconstructed explicitly because the nonlinear-
          !parameter generation cannot detect a metadata-only change.
          NumOfFailures=0
          call GenerateRndIntSeq(nfo,IndOptSequence)
          do i=1,nfo
            ii=IndOptSequence(i)
            j=Glob_Index(nfru+ii,1)
            j2=Glob_Index(nfru+ii,2)
            jbest=j
            jbest2=j2
            if (.not.Glob_IsIndexFixed(1)) then
              do jj=1,Glob_n
                if (jj/=j) then
                  Glob_Index(nfru+ii,1)=jj
                  call EvaluateQAppendedTrial(nfru,K,Evalue,ErrCode)
                  if (ErrCode/=Q_METHOD_SUCCESS) then
                    NumOfFailures=NumOfFailures+1
                    Glob_Index(nfru+ii,1)=jbest
                    if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                      if (Glob_ProcID==0) then
                        write(*,*) 'Error EC0128 in BasisEnlQ: number of failures in energy calculations'
                        write(*,*) 'during the optimization of the first index exceeded limit'
                      endif
                      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
                    endif
                  else if (Evalue<Glob_CurrEnergy) then
                    Glob_CurrEnergy=Evalue
                    jbest=jj
                  endif
                endif
              enddo
            endif
            Glob_Index(nfru+ii,1)=jbest
            if (.not.Glob_IsIndexFixed(2)) then
              do jj2=1,Glob_n
                if (jj2/=j2) then
                  Glob_Index(nfru+ii,2)=jj2
                  call EvaluateQAppendedTrial(nfru,K,Evalue,ErrCode)
                  if (ErrCode/=Q_METHOD_SUCCESS) then
                    NumOfFailures=NumOfFailures+1
                    Glob_Index(nfru+ii,2)=jbest2
                    if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                      if (Glob_ProcID==0) then
                        write(*,*) 'Error EC0129 in BasisEnlQ: number of failures in energy calculations'
                        write(*,*) 'during the optimization of the second index exceeded limit'
                      endif
                      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
                    endif
                  else if (Evalue<Glob_CurrEnergy) then
                    Glob_CurrEnergy=Evalue
                    jbest2=jj2
                  endif
                endif
              enddo
            endif
            Glob_Index(nfru+ii,2)=jbest2
          enddo
          call EvaluateQAppendedTrial(nfru,K,Glob_CurrEnergy,ErrCode)
          if (ErrCode/=Q_METHOD_SUCCESS) then
            if (Glob_ProcID==0) write(*,*) &
              'Error EC0130 in BasisEnlQ: optimized indices cannot be reconstructed'
            call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
          endif

          IV(1:LIV)=IV_init(1:LIV)
          V(1:LV)=V_init(1:LV)
          if (nfru>=nfo) then
            t=max(abs((E_init-Glob_CurrEnergy))/(abs(E_init)+abs(Glob_CurrEnergy)), &
                  10000*epsilon(Glob_CurrEnergy))
          else
            t=ONE
          endif
          D(1:nv)=t
          ExitNeeded=(MaxEnergyEval<=0)
          NumOfFailures=0
          E_best=Glob_CurrEnergy

          do while (.not.ExitNeeded)
            if (Glob_ProcID==0) call DRMNG(D,Glob_CurrEnergy,grad,IV,LIV,LV,nv,V,x)
            call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            select case (IV(1))
            case (1)
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              Evalue=EnergyQA(.true.,ErrCode)
              NumOfEnergyEval=NumOfEnergyEval+1
              if (ErrCode/=Q_METHOD_SUCCESS) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
              else
                Glob_CurrEnergy=Evalue
                if (Evalue<E_best) then
                  E_best=Evalue
                  x_best(1:nv)=x(1:nv)
                endif
              endif
            case (2)
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              call EnergyQB(Evalue,grad,.true.,ErrCode)
              NumOfGradEval=NumOfGradEval+1
              if (ErrCode/=Q_METHOD_SUCCESS) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
              else if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nv)=x(1:nv)
              endif
            case (3:10)
              ExitNeeded=.true.
            endselect
            if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
              if (Glob_ProcID==0) write(*,*) &
                'Warning WC0112 in BasisEnlQ: Q evaluation failures reached the limit'
            endif
            if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
          enddo
        endselect

        do i=1,nfo
          Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
        enddo
        IsEnergyBad=.false.
        Evalue=EnergyQAM(.true.,ErrCode)
        if (ErrCode==Q_METHOD_SUCCESS) then
          Glob_CurrEnergy=Evalue
        else
          IsEnergyBad=.true.
          Glob_CurrEnergy=E_init
          if (Glob_ProcID==0) write(*,*) &
            'Warning WC0113 in BasisEnlQ: optimized candidate is rejected after a failed solve'
        endif

        IsOverlapBad=.false.
        if ((ErrCode==Q_METHOD_SUCCESS).and.(OverlapThreshold>ZERO)) then
          ii=0
          do i=nfru+1,K
            do j=1,i-1
              if (abs(QCanonicalMatrixElement(Glob_S,i,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                Glob_CurrEnergy=E_init
                if (Glob_ProcID==0) then
                  if (ii==1) write(*,*) &
                    'Warning WC0110: overlap threshold exceeded; Q candidate is rejected'
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') &
                    ii,':',i,j,'    S='
                  call writerealadv(6,QCanonicalMatrixElement(Glob_S,i,j))
                endif
              endif
            enddo
          enddo
        endif

        IsAnyLinCoeffBad=.false.
        if ((ErrCode==Q_METHOD_SUCCESS).and.(LinCoeffThreshold>ZERO)) then
          ii=0
          do i=1,K
            if (abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              Glob_CurrEnergy=E_init
              if (Glob_ProcID==0) then
                if (ii==1) write(*,*) &
                  'Warning WC0111: linear-coefficient threshold exceeded; Q candidate is rejected'
                write(*,'(1x,i6,a1,i6,a6)',advance='no') ii,':',i,'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif
        AttemptToGetGoodFunc=AttemptToGetGoodFunc+1
      enddo

      if (IsOverlapBad.or.IsAnyLinCoeffBad.or.IsEnergyBad) then
        if (Glob_ProcID==0) write(*,*) &
          'Error EC0130 in BasisEnlQ: unable to construct an acceptable candidate block'
        call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
      endif

      if (Glob_ProcID==0) then
        write(*,'(1x,a,1x,i0,a,i0)') &
          'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
        write(*,*) 'E=',Glob_CurrEnergy
        do i=1,nfo
          write(*,'(1x,i6,a1)',advance='no') nfru+i,':'
          call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
        enddo
        if (wmu==1) then
          rgm1_counter=rgm1_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/ &
                  Glob_NonlinParam(j,wbfu+i-1)
                ms1=ms1+abs(t)
              enddo
            enddo
          endif
        endif
        if (wmu==2) then
          rgm2_counter=rgm2_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i))/ &
                  Glob_NonlinParam(j,wbfu+i-1)
                ms2=ms2+abs(t)
              enddo
            enddo
          endif
        endif
      endif

      Glob_CurrBasisSize=K
      do i=1,nfo
        Glob_History(nfru+i)%Energy=Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone=0
        Glob_History(nfru+i)%InitFuncAtLastStep=0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt=0
        Glob_FuncNum(nfru+i)=nfru+i
      enddo
      if (Glob_ProcID==0) call SaveResults(Sort='no')
    enddo

    call StoreMatricesInSwapFile()
    call ClearQWorkspace()

    deallocate(V_init)
    deallocate(V)
    deallocate(D)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x)
    deallocate(ParSetBest)
    deallocate(ParSet)
    deallocate(IndOptSequence)
    deallocate(IndSetBest)
    deallocate(IndSet)
    deallocate(Glob_WkGR)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_D)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      write(*,*) 'Random selection statistics:'
      write(*,'(1x,a,1x,i0,1x,a)') &
        'Method 1 of generating basis functions was used',rgm1_counter,'times'
      if (rgm1_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms1/(npt*rgm1_counter)
      write(*,'(1x,a,1x,i0,1x,a)') &
        'Method 2 of generating basis functions was used',rgm2_counter,'times'
      if (rgm2_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms2/(npt*rgm2_counter)
      write(*,*)
      write(*,*) 'Routine BasisEnlQ finished'
    endif
    ErrorCode=Q_METHOD_SUCCESS

  end subroutine BasisEnlQ

  subroutine OptCycleQ(K,FuncBegin,FuncEnd,NumOfFuncToOpt,NumOfFuncToShift, &
                       NumCycles,MaxEnergyEval,OverlapThreshold,LinCoeffThreshold, &
                       SavingFreq,ErrorCode)
!Subroutine OptCycleQ follows OptCycleG's scheduling, DRMNG, acceptance,
!failure limits, saving, and history semantics. G reverses the requested range
!before placing it at the end, so its optimizer block for CurrFunc is ordered
!from CurrFunc+NumActive-1 down to CurrFunc. Q puts those descending
!canonical indices in the active map without calling ReverseFuncOrder,
!PermuteFunctions, or any matrix permutation helper.
!
!Arguments:
    integer,intent(in)  :: K,FuncBegin,FuncEnd,NumOfFuncToOpt,NumOfFuncToShift
    integer,intent(in)  :: NumCycles,MaxEnergyEval,SavingFreq
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
    integer,intent(out) :: ErrorCode
!Local variables:
    integer i,j,a,ii,totsteps
    integer npt,nfo,nv,nvmax,cbs,CurrCycle,CurrFunc,CurrFuncBegin
    integer NumOfFailures,NumOfEnergyEval,NumOfGradEval,ErrCode
    integer ActiveI
    integer,allocatable :: ActiveFunction(:)
    real(wp) Evalue,E_best,E_prev,t
    real(wp),allocatable :: x(:),grad(:),x_init(:),x_best(:)
    logical IsSwapFileOK,ExitNeeded,LastIter
    logical IsOverlapBad,IsAnyLinCoeffBad
!Arrays used by DRMNG
    real(wp),allocatable :: D(:),V(:),V_init(:)
    integer,parameter :: LIV=60
    integer IV(LIV),IV_init(LIV),LV,ALG

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    cbs=Glob_CurrBasisSize
    if (K/=cbs) return
    if ((FuncBegin<1).or.(FuncEnd<FuncBegin).or.(FuncEnd>cbs)) return
    if ((NumOfFuncToOpt<1).or.(NumOfFuncToShift<1)) return
    if ((NumCycles<0).or.(MaxEnergyEval<0).or.(SavingFreq<1)) return

    Glob_GSEPSolutionMethod='Q'
    Glob_OverlapPenaltyAllowed=.false.
    npt=Glob_npt
    nvmax=NumOfFuncToOpt*npt
    Glob_HSLeadDim=cbs
    Glob_HSBuffLen=cbs*NumOfFuncToOpt
    Glob_nfa=cbs

!Checking if cyclic optimization is already completed for this basis size
    if ((Glob_History(cbs)%CyclesDone>=NumCycles).and. &
        (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd)) then
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,*) 'Routine OptCycleQ started'
        write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
        write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
          FuncBegin,'-',FuncEnd,'is already completed'
        write(*,*) 'Exiting OptCycleQ...'
        write(*,*) 'Routine OptCycleQ finished'
      endif
      ErrorCode=Q_METHOD_SUCCESS
      return
    endif

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine OptCycleQ started'
      write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
      write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin,'-',FuncEnd,'will be performed'
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval',MaxEnergyEval
    endif

!Allocate the canonical physical matrices and matrix-element work buffers.
!Unlike G, Q stores both diagonals directly in Glob_H and Glob_S and therefore
!does not allocate Glob_diagH or DSYGVX workspace.
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_D(2*npt,NumOfFuncToOpt,cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_WkGR(nvmax))

!Allocate arrays used by DRMNG and the canonical active-index schedule.
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2+1
    allocate(V(LV))
    allocate(V_init(LV))
    allocate(x(nvmax))
    allocate(x_init(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))
    allocate(ActiveFunction(NumOfFuncToOpt))

    call PrepareQWorkspace(cbs,cbs,NumOfFuncToOpt,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0150 in OptCycleQ: Q workspace cannot be allocated'
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif

!Call DIVSET to get default values in IV and V arrays. These settings are kept
!identical to OptCycleG so a G/Q comparison changes the eigensolver only.
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
    IV_init(1)=12

!Restore canonical physical matrices when possible. A missing swap file causes
!one full matrix-element assembly, after which every optimizer evaluation
!changes only the active rows and columns.
    Glob_H=ZERO
    Glob_S=ZERO
    Glob_diagS=ZERO
    Glob_c=ONE
    call ReadSwapFileAndDistributeData(IsSwapFileOK)
    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,*) &
        'Computing matrix elements and constructing fresh QR factors...'
      call ComputeMatElem(1,cbs)
    else
      if (Glob_ProcID==0) write(*,*) 'Constructing fresh QR factors...'
    endif
    Q_Workspace%MatricesAreCanonical=.true.
    call FactorizeQFresh(ErrCode)
    if (ErrCode==Q_METHOD_SUCCESS) call SolveQ(Glob_CurrEnergy,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,'(1x,a,1x,i0)') &
        'Error EC0151 in OptCycleQ: initial Q energy cannot be computed, status',ErrCode
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif
    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

!Normalize restart history exactly as in OptCycleG. Function numbers stay in
!their input positions, so no preliminary or restart permutation is required.
    if (Glob_History(cbs)%InitFuncAtLastStep<FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    if (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd) then
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
    endif

    totsteps=0
    do CurrCycle=Glob_History(cbs)%CyclesDone+1,NumCycles
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0,1x,a)') 'Cycle',CurrCycle,'began'
      endif
      CurrFuncBegin=Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      do CurrFunc=CurrFuncBegin,FuncEnd,NumOfFuncToShift
        totsteps=totsteps+1
        nfo=min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
        Glob_nfo=nfo
        Glob_nfru=cbs-nfo
        nv=nfo*npt

        !G reverses the requested range before moving it to the trailing block.
        !This descending map gives DRMNG the same block order while every
        !physical basis function remains at its canonical index.
        do a=1,nfo
          ActiveFunction(a)=CurrFunc+nfo-a
        enddo
        call SetQActiveFunctions(ActiveFunction(1:nfo),ErrCode)
        if (ErrCode==Q_METHOD_SUCCESS) call CaptureQMatrixParameters(ErrCode)
        if (ErrCode/=Q_METHOD_SUCCESS) then
          if (Glob_ProcID==0) write(*,*) &
            'Error EC0152 in OptCycleQ: active Q transaction cannot be initialized'
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif

        if (Glob_ProcID==0) then
          write(*,*)
          if (nfo>1) then
            write(*,'(1x,a,1x,i0,a,i0)') &
              'Optimizing functions',CurrFunc,'-',CurrFunc+nfo-1
          else
            write(*,'(1x,a,1x,i0)') 'Optimizing function',CurrFunc
          endif
          if (Glob_AreParamPrintedInCycleOptX) then
            write(*,*) 'Nonlinear parameters before optimization:'
            do a=1,nfo
              ActiveI=ActiveFunction(a)
              write(*,'(1x,i6,a1)',advance='no') Glob_FuncNum(ActiveI),':'
              call writerealarradv(6,Glob_NonlinParam(1:npt,ActiveI),npt)
            enddo
          endif
        endif

        IV(1:LIV)=IV_init(1:LIV)
        V(1:LV)=V_init(1:LV)
        do a=1,nfo
          ActiveI=ActiveFunction(a)
          x((a-1)*npt+1:a*npt)=Glob_NonlinParam(1:npt,ActiveI)
          x_init((a-1)*npt+1:a*npt)=Glob_NonlinParam(1:npt,ActiveI)
        enddo

        t=max(ONE/(cbs*cbs*sqrt(ONE*cbs)),10000*epsilon(Glob_CurrEnergy))
        do a=1,nfo
          do j=1,npt
            D(npt*(a-1)+j)=t
          enddo
        enddo

        ExitNeeded=.false.
        NumOfFailures=0
        NumOfEnergyEval=0
        NumOfGradEval=0
        if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        E_best=Glob_CurrEnergy
        x_best(1:nv)=x(1:nv)
        E_prev=Glob_CurrEnergy

        do while (.not.ExitNeeded)
          if (Glob_ProcID==0) call DRMNG(D,Glob_CurrEnergy,grad,IV,LIV,LV,nv,V,x)
          call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          select case (IV(1))
          case (1)
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do a=1,nfo
              ActiveI=ActiveFunction(a)
              Glob_NonlinParam(1:npt,ActiveI)=x((a-1)*npt+1:a*npt)
            enddo
            Evalue=EnergyQA(.true.,ErrCode)
            NumOfEnergyEval=NumOfEnergyEval+1
            if (ErrCode/=Q_METHOD_SUCCESS) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
            else
              Glob_CurrEnergy=Evalue
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nv)=x(1:nv)
              endif
            endif
          case (2)
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do a=1,nfo
              ActiveI=ActiveFunction(a)
              Glob_NonlinParam(1:npt,ActiveI)=x((a-1)*npt+1:a*npt)
            enddo
            call EnergyQB(Evalue,grad,.true.,ErrCode)
            NumOfGradEval=NumOfGradEval+1
            if (ErrCode/=Q_METHOD_SUCCESS) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
            else
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nv)=x(1:nv)
              endif
            endif
          case (3:10)
            ExitNeeded=.true.
          endselect
          if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
            if (Glob_ProcID==0) then
              write(*,'(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0150 in OptCycleQ: number of failed Q energy or gradient', &
                'calculations during nonlinear-parameter optimization reached', &
                'the limit of',Glob_MaxEnergyFailsAllowed
            endif
          endif
          if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        enddo

        !DRMNG's last requested point need not be its lowest-energy point.
        !Reassemble x_best transactionally so nonlinear parameters, physical
        !matrices, factors, energy, and coefficients all describe one state.
        do a=1,nfo
          ActiveI=ActiveFunction(a)
          Glob_NonlinParam(1:npt,ActiveI)=x_best((a-1)*npt+1:a*npt)
        enddo
        Evalue=EnergyQAM(.true.,ErrCode)
        if (ErrCode==Q_METHOD_SUCCESS) then
          Glob_CurrEnergy=Evalue
        else
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,a)') &
              'Warning WC0151 in OptCycleQ: failed to evaluate the best point.', &
              'The original nonlinear parameters will be restored.'
          endif
        endif

        !Check every unordered overlap pair touching the active set once. This
        !includes inactive functions with canonical indices greater than an
        !active function, which a simple lower-triangle prefix loop would miss.
        IsOverlapBad=.false.
        if ((ErrCode==Q_METHOD_SUCCESS).and.(OverlapThreshold>ZERO)) then
          ii=0
          do a=1,nfo
            ActiveI=ActiveFunction(a)
            do j=1,cbs
              if (j==ActiveI) cycle
              if ((Q_Workspace%ActivePosition(j)>0).and. &
                  (Q_Workspace%ActivePosition(j)<a)) cycle
              if (abs(QCanonicalMatrixElement(Glob_S,ActiveI,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                if (Glob_ProcID==0) then
                  if (ii==1) write(*,*) &
                    'Warning WC0152: overlap of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') &
                    ii,':',Glob_FuncNum(ActiveI),Glob_FuncNum(j),'    S='
                  call writerealadv(6,QCanonicalMatrixElement(Glob_S,ActiveI,j))
                endif
              endif
            enddo
          enddo
        endif

        IsAnyLinCoeffBad=.false.
        if ((ErrCode==Q_METHOD_SUCCESS).and.(LinCoeffThreshold>ZERO)) then
          ii=0
          do i=1,cbs
            if (abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              if (Glob_ProcID==0) then
                if (ii==1) write(*,*) &
                  'Warning WC0153: absolute value of linear parameters exceeds threshold. ', &
                  'Nonlinear parameters will be left unchanged'
                write(*,'(1x,i6,a1,i6,a6)',advance='no') &
                  ii,':',Glob_FuncNum(i),'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif

        if ((Glob_ProcID==0).and.(ErrCode==Q_METHOD_SUCCESS).and. &
            (.not.IsOverlapBad).and.(.not.IsAnyLinCoeffBad)) then
          write(*,'(1x,a,1x,i0,a,i0)') &
            'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
          write(*,*) 'E=',Glob_CurrEnergy
          if (Glob_AreParamPrintedInCycleOptX) then
            write(*,*) 'Nonlinear parameters after optimization:'
            do a=1,nfo
              ActiveI=ActiveFunction(a)
              write(*,'(1x,i6,a1)',advance='no') Glob_FuncNum(ActiveI),':'
              call writerealarradv(6,Glob_NonlinParam(1:npt,ActiveI),npt)
            enddo
          endif
        endif

        if ((ErrCode/=Q_METHOD_SUCCESS).or.IsOverlapBad.or.IsAnyLinCoeffBad) then
          do a=1,nfo
            ActiveI=ActiveFunction(a)
            Glob_NonlinParam(1:npt,ActiveI)=x_init((a-1)*npt+1:a*npt)
          enddo
          Evalue=EnergyQA(.true.,ErrCode)
          if (ErrCode==Q_METHOD_SUCCESS) then
            Glob_CurrEnergy=Evalue
          else
            !ApplyQTrial restores a provable matrix/factor generation on an
            !update failure. A solve failure leaves the restored physical
            !point represented, so retaining E_prev is safe for scheduling.
            if (Glob_ProcID==0) write(*,'(1x,a,1x,a)') &
              'Warning WC0154 in OptCycleQ: original-point energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_CurrEnergy=E_prev
          endif
        endif

        LastIter=(CurrFunc>FuncEnd-NumOfFuncToShift)
        Glob_History(cbs)%Energy=Glob_CurrEnergy
        if (LastIter) then
          Glob_History(cbs)%InitFuncAtLastStep=0
          Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
        else
          Glob_History(cbs)%InitFuncAtLastStep=CurrFunc
        endif

        if (Glob_ProcID==0) then
          if ((totsteps<=Glob_MinMandSavSteps).or. &
              (mod(totsteps,SavingFreq)==0).or. &
              (CurrFunc+NumOfFuncToShift>=FuncEnd)) then
            !Canonical order is already the user-visible order; Q never needs
            !SaveResults' sorting workspace or a temporary basis permutation.
            call SaveResults(Sort='no')
          endif
        endif
      enddo

      if (Glob_ProcID==0) then
        write(*,*)
        write(*,*) 'Cycle',CurrCycle,' finished'
      endif
      if (CurrCycle/=NumCycles) &
        Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    enddo

    call StoreMatricesInSwapFile()
    call ClearQWorkspace()

    deallocate(ActiveFunction)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x_init)
    deallocate(x)
    deallocate(V_init)
    deallocate(V)
    deallocate(D)
    deallocate(Glob_WkGR)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_diagS)
    deallocate(Glob_S)
    deallocate(Glob_H)

    ErrorCode=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) write(*,*) 'Routine OptCycleQ finished'

  end subroutine OptCycleQ

  subroutine FullOpt1Q(InitFunc,FinalFunc,MaxEnergyEval,OverlapThreshold,MaxOverlapPenalty, &
                       DataSaveMinTimeInterv,HessianSaveMinTimeInterv,HessFileName,ErrorCode)
!Subroutine FullOpt1Q simultaneously optimizes the nonlinear parameters of
!canonical functions InitFunc:FinalFunc. It follows FullOpt1G's DRMNG controls,
!smooth overlap penalty, Hessian restart/save policy, timed data saves, and
!history accounting. The active range is represented by an explicit map and is
!never moved to the end of the basis; consequently QR factors, H, S,
!coefficients, function identities, and Hessian blocks remain in one stable
!ordering throughout the operation.
!
!Arguments:
    integer,intent(in)     :: InitFunc,FinalFunc,MaxEnergyEval
    real(wp),intent(in)    :: OverlapThreshold,MaxOverlapPenalty
    real(4),intent(in)     :: DataSaveMinTimeInterv,HessianSaveMinTimeInterv
    character(*),intent(in) :: HessFileName
    integer,intent(out)    :: ErrorCode
!Local variables:
    integer i,j,npt,nfo,nfa,nv,ErrCode
    integer NumOfEnergyEval,NumOfGradEval,NumOfFailures
    integer NumOfEnergyEvalDuringFullOpt_Init
    integer,allocatable :: ActiveFunction(:)
    logical IsSwapFileOK,ExitNeeded,SaveHessian,IsHessFileOK,IsHessSaveSuccess
    real(wp) Evalue,CurrentEnergy,t
    real(wp) MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap
    real(4) TimeOfLastSave,TimeOfLastHessSave
    real(wp),allocatable :: x(:),grad(:),D(:),V(:)
    integer,parameter :: LIV=60
    integer IV(LIV),LV,ALG,IVLMAT

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    nfa=Glob_CurrBasisSize
    if ((InitFunc<1).or.(FinalFunc<InitFunc).or.(FinalFunc>nfa)) return
    if (MaxEnergyEval<0) return
    nfo=FinalFunc-InitFunc+1
    npt=Glob_npt
    nv=nfo*npt

    if (OverlapThreshold>=ONE) then
      Glob_OverlapPenaltyAllowed=.false.
    else
      Glob_OverlapPenaltyAllowed=.true.
      Glob_OverlapPenaltyThreshold2=OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty=MaxOverlapPenalty
    endif

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine FullOpt1Q started'
      write(*,*) 'Simultaneous optimization of nonlinear parameters of basis functions'
      write(*,*) InitFunc,'  through',FinalFunc,'  will be attempted'
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Overlap threshold is ',abs(OverlapThreshold)
        write(*,*) 'Max value of a pair overlap penalty is ',Glob_MaxOverlapPenalty
        write(*,*) 'Displayed optimizer objectives include the overlap penalty'
      else
        write(*,*) 'No constraints on overlaps will be imposed'
      endif
    endif

    Glob_GSEPSolutionMethod='Q'
    Glob_nfa=nfa
    Glob_nfo=nfo
    Glob_nfru=nfa-nfo
    Glob_HSLeadDim=nfa
    Glob_HSBuffLen=max(min(nfa*(nfa+1)/2,1000),30*nfa)

    allocate(Glob_H(nfa,nfa))
    allocate(Glob_S(nfa,nfa))
    allocate(Glob_diagS(nfa))
    allocate(Glob_D(2*npt,nfo,nfa))
    allocate(Glob_c(nfa))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_WkGR(nv))
    allocate(x(nv))
    allocate(grad(nv))
    allocate(ActiveFunction(nfo))
    LV=71+nv*(nv+13)/2+1
    if (Glob_ProcID==0) then
      allocate(D(nv))
      allocate(V(LV))
    endif

    call PrepareQWorkspace(nfa,nfa,nfo,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0160 in FullOpt1Q: Q workspace cannot be allocated'
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif
    do i=1,nfo
      ActiveFunction(i)=InitFunc+i-1
    enddo

    if ((HessFileName==' ').or.(HessFileName=='none').or. &
        (HessFileName=='NONE').or.(HessFileName=='None')) then
      SaveHessian=.false.
    else
      SaveHessian=.true.
    endif
    IsHessFileOK=.false.

!Restore or assemble the full canonical physical problem once, then install
!the active map and construct the initial factors. No temporary permutations
!or SaveResults sorting workspace are required for an interior active range.
    Glob_H=ZERO
    Glob_S=ZERO
    Glob_diagS=ZERO
    Glob_c=ONE
    call ReadSwapFileAndDistributeData(IsSwapFileOK)
    if (.not.IsSwapFileOK) call ComputeMatElem(1,nfa)
    Q_Workspace%MatricesAreCanonical=.true.
    call SetQActiveFunctions(ActiveFunction,ErrCode)
    if (ErrCode==Q_METHOD_SUCCESS) call CaptureQMatrixParameters(ErrCode)
    if (ErrCode==Q_METHOD_SUCCESS) call FactorizeQFresh(ErrCode)
    if (ErrCode==Q_METHOD_SUCCESS) Glob_CurrEnergy=EnergyQA(.false.,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,'(1x,a,1x,i0)') &
        'Error EC0161 in FullOpt1Q: initial Q energy cannot be computed, status',ErrCode
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif

    call GetQOverlapStatistics(MaxAbsOverlap,MinAbsOverlap, &
      AverageAbsOverlap,ErrCode)
    if (Glob_ProcID==0) then
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                           ',Glob_TotalOverlapPenalty
        write(*,*) 'Initial energy (including overlap penalty)',Glob_CurrEnergy
      else
        write(*,*) 'Initial energy                            ',Glob_CurrEnergy
      endif
      write(*,*) 'Maximal overlap                           ',MaxAbsOverlap
      write(*,*) 'Minimal overlap                           ',MinAbsOverlap
      write(*,*) 'Average abs value of overlap              ',AverageAbsOverlap
    endif

    call CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave=Glob_TimeSinceStart
    TimeOfLastHessSave=Glob_TimeSinceStart

!Initialize DRMNG and preserve its packed Hessian ordering as
![parameters of InitFunc, parameters of InitFunc+1, ...]. This is exactly the
!Q active-map order and therefore remains stable across restarts.
    ALG=2
    if (Glob_ProcID==0) then
      call DIVSET(ALG,IV,LIV,LV,V)
      IV(17)=1000000; IV(18)=1000000
      IV(19)=-1
      IV(20)=0; IV(22)=0; IV(23)=-1; IV(24)=0
      V(31)=0.0_wp
      V(32)=2*epsilon(V(32))
      V(37)=2*epsilon(V(37))
      V(35)=Glob_MaxScStepAllowedInOpt
      IV(1)=12
    endif
    do i=1,nfo
      x((i-1)*npt+1:i*npt)= &
        Glob_NonlinParam(1:npt,ActiveFunction(i))
    enddo

    if (Glob_ProcID==0) then
      if (SaveHessian) then
        IVLMAT=IV(42)
        call ReadHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessFileOK)
        if (IsHessFileOK) IV(25)=0
      endif
      if ((.not.Glob_FullOptSaveD).or.(.not.IsHessFileOK).or. &
          (.not.SaveHessian)) then
        t=max(ONE/(nfa*nfa*sqrt(ONE*nfa)),10000*epsilon(Glob_CurrEnergy))
        D(1:nv)=t
      endif
    endif

    ExitNeeded=(MaxEnergyEval<=0)
    NumOfFailures=0
    NumOfEnergyEval=0
    NumOfGradEval=0
    CurrentEnergy=Glob_CurrEnergy
    NumOfEnergyEvalDuringFullOpt_Init= &
      Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt

    do while (.not.ExitNeeded)
      if (Glob_ProcID==0) call DRMNG(D,CurrentEnergy,grad,IV,LIV,LV,nv,V,x)
      call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      select case (IV(1))
      case (1)
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,ActiveFunction(i))= &
            x((i-1)*npt+1:i*npt)
        enddo
        Evalue=EnergyQA(.true.,ErrCode)
        NumOfEnergyEval=NumOfEnergyEval+1
        if (ErrCode/=Q_METHOD_SUCCESS) then
          NumOfFailures=NumOfFailures+1
          IV(2)=1
        else
          CurrentEnergy=Evalue
        endif

        if ((Glob_ProcID==0).and.(ErrCode==Q_METHOD_SUCCESS)) then
          if (Evalue<Glob_CurrEnergy) then
            call CPU_TIME(Glob_TimeSinceStart)
            if (Glob_TimeSinceStart-TimeOfLastSave>DataSaveMinTimeInterv) then
              if (Glob_OverlapPenaltyAllowed) then
                Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
              else
                Glob_CurrEnergy=Evalue
              endif
              Glob_History(nfa)%Energy=Glob_CurrEnergy
              Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt= &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              call SaveResults(Sort='no')
              write(*,*) 'Data file has been updated'
              call GetQOverlapStatistics(MaxAbsOverlap,MinAbsOverlap, &
                AverageAbsOverlap,ErrCode)
              write(*,*) 'Some current statistics:'
              if (Glob_OverlapPenaltyAllowed) then
                write(*,*) 'Energy (without overlap penalty)  ', &
                  Evalue-Glob_TotalOverlapPenalty
                write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
                write(*,*) 'Energy (including overlap penalty)',Evalue
              else
                write(*,*) 'Energy                            ',Evalue
              endif
              write(*,*) 'Maximal overlap                   ',MaxAbsOverlap
              write(*,*) 'Minimal overlap                   ',MinAbsOverlap
              write(*,*) 'Average abs value of overlap      ',AverageAbsOverlap
              TimeOfLastSave=Glob_TimeSinceStart
            endif
            if ((Glob_TimeSinceStart-TimeOfLastHessSave> &
                 HessianSaveMinTimeInterv).and.SaveHessian) then
              call SaveHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessSaveSuccess)
              TimeOfLastHessSave=Glob_TimeSinceStart
            endif
          endif
        endif
      case (2)
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,ActiveFunction(i))= &
            x((i-1)*npt+1:i*npt)
        enddo
        call EnergyQB(Evalue,grad,.true.,ErrCode)
        NumOfGradEval=NumOfGradEval+1
        if (ErrCode/=Q_METHOD_SUCCESS) then
          NumOfFailures=NumOfFailures+1
          IV(2)=0
        endif
      case (3:10)
        ExitNeeded=.true.
      endselect

      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
        if (Glob_ProcID==0) write(*,*) &
          'Error EC0162 in FullOpt1Q: Q evaluation failures exceeded limit'
        call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
      endif
      if (NumOfEnergyEval>=MaxEnergyEval) then
        if (Glob_ProcID==0) then
          write(*,*) 'Warning WC0130 in FullOpt1Q: number of energy evaluations reached limit'
          write(*,*) 'Optimization is terminated'
        endif
        ExitNeeded=.true.
      endif
    enddo

!DRMNG returns its current accepted x. Re-evaluate it transactionally so the
!saved physical matrices, QR factors, coefficients, and parameters all belong
!to the same final generation.
    do i=1,nfo
      Glob_NonlinParam(1:npt,ActiveFunction(i))=x((i-1)*npt+1:i*npt)
    enddo
    Evalue=EnergyQA(.true.,ErrCode)
    if (ErrCode/=Q_METHOD_SUCCESS) then
      if (Glob_ProcID==0) write(*,'(1x,a,1x,i0)') &
        'Error EC0163 in FullOpt1Q: final Q energy cannot be computed, status',ErrCode
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif
    if (Glob_OverlapPenaltyAllowed) then
      Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
    else
      Glob_CurrEnergy=Evalue
    endif

    call GetQOverlapStatistics(MaxAbsOverlap,MinAbsOverlap, &
      AverageAbsOverlap,ErrCode)
    if (Glob_ProcID==0) then
      write(*,*)
      write(*,'(1x,a,1x,i0,a,i0)') &
        'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
      write(*,*) 'Final energy and overlap statistics:'
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Energy (without overlap penalty)  ', &
          Evalue-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
        write(*,*) 'Energy (including overlap penalty)',Evalue
      else
        write(*,*) 'Energy                            ',Evalue
      endif
      write(*,*) 'Maximal overlap                   ',MaxAbsOverlap
      write(*,*) 'Minimal overlap                   ',MinAbsOverlap
      write(*,*) 'Average abs value of overlap      ',AverageAbsOverlap
    endif

    Glob_History(nfa)%Energy=Glob_CurrEnergy
    Glob_History(nfa)%NumOfEnergyEvalDuringFullOpt= &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
    if (Glob_ProcID==0) call SaveResults(Sort='no')
    call StoreMatricesInSwapFile()
    Glob_OverlapPenaltyAllowed=.false.
    call ClearQWorkspace()

    deallocate(ActiveFunction)
    deallocate(grad)
    deallocate(x)
    if (Glob_ProcID==0) then
      deallocate(V)
      deallocate(D)
    endif
    deallocate(Glob_WkGR)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_diagS)
    deallocate(Glob_S)
    deallocate(Glob_H)

    ErrorCode=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) write(*,*) 'Routine FullOpt1Q finished'

  end subroutine FullOpt1Q

  subroutine BasisEnlG(Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval, &
                       OverlapThreshold,LinCoeffThreshold)
!Subroutine BasisEnlG enlarges the basis set from initial
!Kstart-1 functions to Kstop functions by means of trials of
!randomly selected candidates and then optimizing them. At each
!trial it generates a set of Kstep functions for NTrials times
!based on the existing distribution of nonlinear parameters.
!Only the set that lowers the energy the most is left.
!After that the nonlinear parameters
!are optimized according to the value of OptimizationType.
!MaxEnergyEval is the maximal number of the energy evaluations allowed
!for this optimization.
!Each newly accepted basis function is checked for pair linear
!dependency with other functions in the basis. If the absolute
!value of the overlap is greater than OverlapThreshold then
!such function is rejected and random trials and optimization
!take place again until either a good function/functions are
!generated or a certain limit of such repetitions is reached
!(the limit is given by the value of global constant
!Glob_BadOverlapOrLinCoeffLim). If the value of parameter OverlapThreshold
!is set to negative or zero then no such check is performed.
!A Similar check is performed for all linear parameters. If
!adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!by magnitude than LinCoeffThreshold then such new function/functions
!are rejected and the procedure is repeated until a good set is generated
!or a certain number of failures happens (Glob_BadOverlapOrLinCoeffLim).
!To avoid this check it is enough to set LinCoeffThreshold to a value
!equal or smaller than zero.
!It is assumed that the basis of Kstart-1 functions has already
!been selected and the corresponding nonlinear parameters are
!stored in proper global arrays. Subroutines EnergyGA, EnergyGAM, and EnergyGB
!are called to evaluate the energy and its gradient.

!Arguments:
    integer,intent(in)     :: Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
    integer      i,j,K,AttemptToGetGoodFunc,ii,jj,jbest,ii2,jj2,jbest2,j2,q
    integer      np,npt,nfo,nfa,nfru,nfrup1,nvmax,nv
    integer      OpenFileErr,ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
    logical      IsSwapFileOK,IsEnergyImproved,ExitNeeded
    logical      IsOverlapBad,IsAnyLinCoeffBad,IsEnergyBad
    integer      wbfu_t,wmu_t,wbfu,wmu,rgm1_counter,rgm2_counter,BlockSizeForDSYGVX
    real(wp)  ms1,ms2
    real(wp)  Evalue,E_init,E_best
    real(wp)  t
    real(wp),allocatable,dimension(:,:)   :: ParSet,ParSetBest
    integer,allocatable,dimension(:,:)    :: IndSet,IndSetBest
    integer,allocatable,dimension(:)      :: IndOptSequence
    real(wp),allocatable,dimension(:)     :: x,x_best,grad
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V,V_init
    integer,parameter    :: LIV=60
    integer                 IV(LIV), IV_init(LIV)
    integer                 LV
    integer                 ALG
!Allocatable work space
    real(wp),allocatable,dimension(:)            :: WorkBuffReal
    integer,allocatable,dimension(:)                :: WorkBuffInt
    type(Glob_HistoryStep),allocatable,dimension(:) :: TempHistory
    real(wp),allocatable,dimension(:,:)          :: TempParam
    integer,allocatable,dimension(:)                :: TempFunc
!====================================================
!These variables are used when a finite difference gradient is computed
!real(wp),allocatable,dimension(:)     ::    fx,fgrad
!real(wp)                                    deltax,Evalue1
!====================================================

    wbfu_t=0
    wmu_t=0

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*)  'Routine BasisEnlG started'
      write(*,'(1x,a,1x,i0)') 'Kstart =',Kstart
      write(*,'(1x,a,1x,i0)') 'Kstop =',Kstop
      write(*,'(1x,a,1x,i0)') 'Kstep =',Kstep
      write(*,'(1x,a,1x,i0)') 'OptimizationType =',OptimizationType
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval =',MaxEnergyEval
    endif

!Setting the values of some global variables
    Glob_GSEPSolutionMethod='G'
    Glob_OverlapPenaltyAllowed=.false.
    Glob_nfa=Kstart+Kstep
    Glob_nfo=Kstep
    Glob_HSLeadDim=Kstop
    Glob_HSBuffLen=Kstop*Kstep
    np=Glob_np
    npt=Glob_npt
    nfo=Glob_nfo
    nfa=Glob_nfa
    nvmax=Kstep*Glob_npt

!Reallocate arrays that contain the information
!about basis functions and optimization process.
    call ReallocateBasisFuncData(Kstop,Glob_CurrBasisSize)

!Allocate some global arrays
    allocate(Glob_H(Kstop,Kstop))
    allocate(Glob_S(Kstop,Kstop))
    allocate(Glob_diagH(Kstop))
    allocate(Glob_diagS(Kstop))
    allocate(Glob_c(Kstop))
    allocate(Glob_D(2*npt,Kstep,Kstop))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',Kstop,Kstop,Kstop,Kstop)
    Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*Kstop,8*Kstop)
    allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    allocate(Glob_IWorkForDSYGVX(5*Kstop))

!Allocate workspace for EnergyGB
    allocate(Glob_WkGR(Kstep*npt))

!Allocate workspace
    allocate(ParSet(npt,Kstep))
    allocate(ParSetBest(npt,Kstep))
    allocate(IndSet(Kstep,2))
    allocate(IndSetBest(Kstep,2))
    allocate(IndOptSequence(Kstep))
    allocate(x(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))

!Allocate arrays used by DRMNG
    nvmax=npt*Kstep
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2 + 1
    allocate(V(LV))
    allocate(V_init(LV))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0 !set summary print format
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
    IV_init(1)=12 !DIVSET has been called and some default values were changed

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Calculating the initial energy
    ErrCode=0
    if (Kstart>1) then
      if (IsSwapFileOK) then
        if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
        Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.false.,ErrCode)
      else
        if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
        Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.true.,ErrCode)
      endif
    else
      Glob_CurrEnergy=huge(Glob_CurrEnergy)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0125 in BasisEnlG: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

    rgm1_counter=0
    rgm2_counter=0
    ms1=ZERO
    ms2=ZERO
    K=Kstart-1

!Main loop begins here
    do while (K<Kstop)
      if (K+Kstep<=Kstop) then
        nfo=Kstep
        nfru=K
        K=K+Kstep
      else
        nfo=Kstop-K
        nfru=K
        K=Kstop
      endif
      call linalg_setparam(K)  !reset linalg flags to account for changes in the basis size
      Glob_nfa=K
      Glob_nfru=nfru
      Glob_nfo=nfo
      nfrup1=nfru+1
      nv=nfo*npt
      E_init=Glob_CurrEnergy
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0)') 'Current basis size is',Glob_CurrBasisSize
        if (nfo>1) then
          write(*,'(1x,a,1x,i0,a,i0)') 'Selecting functions',nfrup1,'-',K
        else
          write(*,'(1x,a,1x,i0)') 'Selecting function',K
        endif
      endif

      IsOverlapBad=.true.
      IsAnyLinCoeffBad=.true.
      IsEnergyBad=.false.
      AttemptToGetGoodFunc=1
      do while ((IsOverlapBad.or.IsAnyLinCoeffBad.or.IsEnergyBad).and. &
                (AttemptToGetGoodFunc<=Glob_BadOverlapOrLinCoeffLim))
        !Random selection
        NumOfFailures=0
        IsEnergyImproved=.false.
        wbfu=0
        wmu=0
        do i=1,NTrials
          if (Glob_ProcID==0) call GenerateTrialParam(nfo,ParSet,IndSet,wbfu_t,wmu_t)
          call MPI_BCAST(ParSet,npt*nfo,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_BCAST(IndSet,2*nfo,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          Glob_NonlinParam(1:npt,nfrup1:K)=ParSet(1:npt,1:nfo)
          Glob_Index(nfrup1:K,1)=IndSet(1:nfo,1)
          Glob_Index(nfrup1:K,2)=IndSet(1:nfo,2)
          Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
          if (ErrCode==0) then
            if (Evalue<Glob_CurrEnergy) then
              Glob_CurrEnergy=Evalue
              ParSetBest(1:npt,1:nfo)=ParSet(1:npt,1:nfo)
              IndSetBest(1:nfo,1)=IndSet(1:nfo,1)
              IndSetBest(1:nfo,2)=IndSet(1:nfo,2)
              IsEnergyImproved=.true.
              wbfu=wbfu_t
              wmu=wmu_t
            endif
          else
            NumOfFailures=NumOfFailures+1
          endif
        enddo
        if (NumOfFailures*ONE/NTrials>Glob_MaxFracOfTrialFailsAllowed) then
          if (Glob_ProcID==0) then
            write(*,*) 'Error EC0126 in BasisEnlG: the number of eigenvalue problem solution failures'
            write(*,*) 'in random selection process exceeded limit'
            write(*,'(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials,'%'
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        if (.not.(IsEnergyImproved)) then
          if (Glob_ProcID==0) then
            write(*,*) 'Error EC0127 in BasisEnlG: random selection did not result'
            write(*,*) 'in any energy improvement'
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        Glob_NonlinParam(1:npt,nfrup1:K)=ParSetBest(1:npt,1:nfo)
        Glob_Index(nfrup1:K,1)=IndSetBest(1:nfo,1)
        Glob_Index(nfrup1:K,2)=IndSetBest(1:nfo,2)
        Glob_CurrEnergy=EnergyGA(nfrup1,K,.true.,ErrCode)
        if (Glob_ProcID==0) then
          write (*,'(1x,a)',advance='no') 'E='
          call writereal(6,Glob_CurrEnergy)
          write (*,'(5x,a,1x,i0)') 'prototype function is',wbfu
          do i=1,nfo
            write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') nfru+i,':',IndSetBest(i,1),IndSetBest(i,2)
            call writerealarradv(6,ParSetBest(1:npt,i),npt)
          enddo
          write (*,*) 'Optimizing nonlinear parameters'
        endif

        !Optimization of the nonlinear parameters
        select case (OptimizationType)
        case(0) !No optimization
          do i=1,nfo
            x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          enddo
        case(1)

          !We first optimize indices (if they are not fixed)
          NumOfFailures=0
          !Generate a random sequence which will define the order in which
          !indices should be optimized (one index at a time)
          call GenerateRndIntSeq(nfo,IndOptSequence)
          !Loop where the indices are optimized. Note: there is some room for improvement
          !here it is programmed in a simple way when all matrix element of functions
          !nfrup1 through K are computed each time while it is not always necessary.
          do i=1,nfo
            ii=IndOptSequence(i) !!  ii should be the same for two indices
            j=Glob_Index(nfru+ii,1)
            j2=Glob_Index(nfru+ii,2)
            jbest=j
            jbest2=j2
            if (Glob_VectorCouplingScheme==2) then
              !Coupling scheme 2 requires the two indices to be equal. They therefore cannot
              !be varied independently and are varied together here. Note that if either of
              !the indices is fixed then both of them are effectively fixed, in which case
              !there is nothing to optimize.
              if (.not.(Glob_IsIndexFixed(1).or.Glob_IsIndexFixed(2))) then
                do jj=1,Glob_n
                  if (jj/=j) then
                    Glob_Index(nfru+ii,1)=jj
                    Glob_Index(nfru+ii,2)=jj
                    Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,1)=jbest
                      Glob_Index(nfru+ii,2)=jbest
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0128 in BasisEnlG: number of failures in energy calculations'
                          write(*,*) 'during the optimization of x,y-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest=jj
                      endif
                    endif
                  endif
                enddo
              endif
              jbest2=jbest
              Glob_Index(nfru+ii,1)=jbest
              Glob_Index(nfru+ii,2)=jbest2
            else
              !Coupling scheme 1 requires the two indices to be different, which is why the
              !value of the other index is excluded from the trial values below. Both index
              !orderings remain allowed, but a single-index move cannot turn one ordering
              !into the other. Scheme 0 imposes no restrictions on the indices.
              if (.not.Glob_IsIndexFixed(1)) then
                do jj=1,Glob_n
                  if ((jj/=j).and.((Glob_VectorCouplingScheme/=1).or.(jj/=j2))) then
                    Glob_Index(nfru+ii,1)=jj
                    Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,1)=jbest
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0128 in BasisEnlG: number of failures in energy calculations'
                          write(*,*) 'during the optimization of x-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest=jj
                      endif
                    endif
                  endif
                enddo
              endif
              Glob_Index(nfru+ii,1)=jbest
              if (.not.Glob_IsIndexFixed(2)) then
                do jj2=1,Glob_n
                  if ((jj2/=j2).and.((Glob_VectorCouplingScheme/=1).or.(jj2/=jbest))) then
                    Glob_Index(nfru+ii,2)=jj2
                    Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,2)=jbest2
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0129 in BasisEnlG: number of failures in energy calculations'
                          write(*,*) 'during the optimization of y-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest2=jj2
                      endif
                    endif
                  endif
                enddo
              endif
              Glob_Index(nfru+ii,2)=jbest2
            endif
          enddo
          
          !Now we optimize nonlinear parameters

          !Setting IV and V values as was in their initial copies
          IV(1:LIV)=IV_init(1:LIV)
          V(1:LV)=V_init(1:LV)

          nv=nfo*npt
          do i=1,nfo
            x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          enddo
          if (nfru>=nfo) then
            t=max(abs((E_init-Glob_CurrEnergy))/(abs(E_init)+abs(Glob_CurrEnergy)), &
                  10000*epsilon(Glob_CurrEnergy))
          else
            t=ONE
          endif
          !if (Glob_ProcID==0) write(*,*) 'scaling coeff=',t !remove later
          do i=1,nfo
            !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
            do j=1,npt
              !Make sure none of the D(i) will be zero or smaller than the threshold
              !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
              D(npt*(i-1)+j)=t
              !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
            enddo
          enddo

          ExitNeeded=.false.
          NumOfFailures=0
          NumOfEnergyEval=0
          NumOfGradEval=0
          if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
          E_best=Glob_CurrEnergy
          x_best(1:nfo*npt)=x(1:nfo*npt)

          do while (.not.(ExitNeeded))
            if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
            call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            select case (IV(1))
            case (1) !Only energy is needed
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
              NumOfEnergyEval=NumOfEnergyEval+1
              if (ErrCode/=0) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
              else
                Glob_CurrEnergy=Evalue
                if (Evalue<E_best) then
                  E_best=Evalue
                  x_best(1:nfo*npt)=x(1:nfo*npt)
                endif
              endif
            case (2) !Only gradient is needed
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              call EnergyGB(Evalue,grad,.true.,ErrCode)
              NumOfGradEval=NumOfGradEval+1
              if (ErrCode/=0) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
              else
                if (Evalue<E_best) then
                  E_best=Evalue
                  x_best(1:nfo*npt)=x(1:nfo*npt)
                endif
              endif
              !===================================
!                  The lines below need to be uncommented when finite
!                  difference gradient is shown
!                  allocate(fx(nvmax))
!                  allocate(fgrad(nvmax))
!                  allocate(fgrad_5_points(nvmax))
!                  write(*,*)'             Analytic              finite diff gradient       5point finite diff gradient'
!                  do j=1,nv
!                        fx(1:nv)=x(1:nv)
!                        deltax=x(j)*1.0Q-3
!                        fx(j)=x(j)+deltax
!                         do i=1,nfo
!                           Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                        enddo
!                        Evalue1=EnergyGA(nfrup1,K,.true.,ErrCode)
!                        fx(j)=x(j)+2*deltax
!                         do i=1,nfo
!                           Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                        enddo
!                        Evalue_plus_2h=EnergyGA(nfrup1,K,.true.,ErrCode)
!                        fx(j)=x(j)-deltax
!                         do i=1,nfo
!                            Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                        enddo
!                        Evalue=EnergyGA(nfrup1,K,.true.,ErrCode)
!                        fx(j)=x(j)-2*deltax
!                         do i=1,nfo
!                            Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                        enddo
!                        Evalue_minus_2h=EnergyGA(nfrup1,K,.true.,ErrCode)
!                        fgrad(j)=(Evalue1-Evalue)/(2*deltax)
!                        fgrad_5_points(j)=((Evalue_minus_2h-Evalue_plus_2h)/12+TWO*(Evalue1-Evalue)/THREE)/deltax
!                    enddo
!                    !f_2nd_der(j)=(4*(Evalue1+Evalue)/3-(Evalue_minus_2h+Evalue_plus_2h)/12-5*Glob_CurrEnergy/TWO)/(deltax*deltax)
!                    do j=1,nv
!                          write(*,*) j,'  ',grad(j),'  ' ,fgrad(j),'   ',fgrad_5_points(j)
!                    enddo
!                    write(*,*)
!                    do i=1,nfo
!                           Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
!                    enddo
!                  deallocate(fgrad)
!                  deallocate(fgrad_5_points)
!                  deallocate(fx)
!
              !===================================
            case (3:8) !Some kind of convergence has been reached
              ExitNeeded=.true.
            case (9:10)   !Function evaluation limit has been reached.
              !This never suppose to happen because we
              !count the number of function evaluations ourselves.
              ExitNeeded=.true.
            endselect
            if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
              if (Glob_ProcID==0) then
                write(*,'(1x,a,1x,a,1x,a,1x,i0)') &
                  'Warning WC0112 in BasisEnlG: number of failures in energy or gradient', &
                  'calculations during the optimization of nonlinear parameters', &
                  'reached the limit of',Glob_MaxEnergyFailsAllowed
              endif
              !call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
            endif
            if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
          enddo
        endselect !(OptimizationType)

        !Calculate the energy and the linear coefficients at the best point found
        do i=1,nfo
          Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
        enddo
        IsEnergyBad=.false.
        Evalue=EnergyGAM(nfrup1,K,.true.,ErrCode)
    !!We run EnergyGA again because EnergyGAM might give slightly different
    !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
        !Glob_CurrEnergy=EnergyGA(nfrup1,K,.false.,ErrCode)
        if (ErrCode==0) then
          Glob_CurrEnergy=Evalue
        else
          !Reject the generated basis function(s) and make another attempt to
          !generate them, just like it is done when the overlap or the linear
          !parameters are found to be bad
          IsEnergyBad=.true.
          Glob_CurrEnergy=E_init
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,a,1x,a)') &
              'Warning WC0113 in BasisEnlG: failed to evaluate energy after the optimization', &
              'of nonlinear parameters. Generated basis function(s) are rejected', &
              'and a new attempt to generate them will be made'
          endif
        endif

        !Checking if overlaps are OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad=.false.
        if (OverlapThreshold>ZERO) then
          ii=0
          do i=nfrup1,K
            do j=1,i-1
              if (abs(Glob_S(i,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                Glob_CurrEnergy=E_init
                if (Glob_ProcID==0) then
                  if (ii==1) then
                    write(*,*) 'Warning WC0110: overlap of the following functions exceeds threshold. ', &
                      'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  endif
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') ii,':',i,j,'    S='
                  call writerealadv(6,Glob_S(i,j))
                endif
              endif
            enddo
          enddo
        endif
        !Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad=.false.
        if (LinCoeffThreshold>ZERO) then
          ii=0
          do i=1,K
            if(abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              Glob_CurrEnergy=E_init
              if (Glob_ProcID==0) then
                if (ii==1) then
                  write(*,*) 'Warning WC0111: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                endif
                write(*,'(1x,i6,a1,i6,a6)',advance='no') ii,':',i,'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif
        AttemptToGetGoodFunc=AttemptToGetGoodFunc+1

      enddo !(IsOverlapBad.or.IsAnyLinCoeffBad).and. &
      !(AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)

      if (Glob_ProcID==0) then
        write (*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
        write (*,*) 'E=',Glob_CurrEnergy
        do i=1,nfo
          write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') nfru+i,':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
          call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
        enddo
      endif

      !getting statistics about the distribution of generated parameters
      if (Glob_ProcID==0) then
        if (wmu==1) then
          rgm1_counter=rgm1_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i)) &
                   /Glob_NonlinParam(j,wbfu+i-1)
                ms1=ms1+abs(t)
              enddo
            enddo
          endif
        endif
        if (wmu==2) then
          rgm2_counter=rgm2_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i)) &
                   /Glob_NonlinParam(j,wbfu+i-1)
                ms2=ms2+abs(t)
              enddo
            enddo
          endif
        endif
      endif
      Glob_CurrBasisSize=K
      do i=1,nfo
        Glob_History(nfru+i)%Energy=Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone=0
        Glob_History(nfru+i)%InitFuncAtLastStep=0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt=0
      enddo
      do i=1,nfo
        Glob_FuncNum(nfru+i)=nfru+i
      enddo

      if (Glob_ProcID==0) call SaveResults(Sort='no')
    enddo
!Main loop ends here

    call StoreMatricesInSwapFile()

!Deallocate arrays used by DRMNG
    deallocate(V_init)
    deallocate(V)
    deallocate(D)

!Deallocate workspace
    deallocate(IndOptSequence)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x)
    deallocate(IndSetBest)
    deallocate(IndSet)
    deallocate(ParSetBest)
    deallocate(ParSet)

!Deallocate workspace for EnergyGB
    deallocate(Glob_WkGR)

!Deallocate workspace for DSYGVX
    deallocate(Glob_IWorkForDSYGVX)
    deallocate(Glob_WorkForDSYGVX)

!Deallocate global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_D)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    deallocate(Glob_diagH)
    deallocate(GLob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      write(*,*) 'Random selection statistics:'
      write(*,'(1x,a,1x,i0,1x,a)') 'Method 1 of generating basis functions was used',rgm1_counter,'times'
      if (rgm1_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms1/(npt*rgm1_counter)
      write(*,'(1x,a,1x,i0,1x,a)') 'Method 2 of generating basis functions was used',rgm2_counter,'times'
      if (rgm2_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms2/(npt*rgm2_counter)
      write(*,*)
      write(*,*) 'Routine BasisEnlG finished'
    endif

  end subroutine BasisEnlG

  subroutine BasisEnlI(Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval, &
                       OverlapThreshold,LinCoeffThreshold)
!Subroutine BasisEnlI enlarges the basis set from initial
!Kstart-1 functions to Kstop functions by means of trials of
!randomly selected candidates and then optimizing them. At each
!trial it generates a set of Kstep functions for NTrials times
!based on the existing distribution of nonlinear parameters.
!Only the set that lowers the energy the most is left. Then
!the Z-index of the best candidate is optimized by trying to
!consequently change Z-indices. After that the nonlinear parameters
!are optimized according to the value of OptimizationType.
!MaxEnergyEval is the maximal number of the energy evaluations allowed
!for this optimization.
!Each newly accepted basis function is checked for pair linear
!dependency with other functions in the basis. If the absolute
!value of the overlap is greater than OverlapThreshold then
!such function is rejected and random trials and optimization
!take place again until either a good function/functions are
!generated or a certain limit of such repetitions is reached
!(the limit is given by the value of global constant
!Glob_BadOverlapOrLinCoeffLim). If the value of parameter OverlapThreshold
!is set to negative or zero then no such check is performed.
!A similar check is performed for all linear parameters. If
!adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!by magnitude than LinCoeffThreshold then such new function/functions
!are rejected and the procedure is repeated until a good set is generated
!or a certain number of failures happens (Glob_BadOverlapOrLinCoeffLim).
!To avoid this check it is enough to set LinCoeffThreshold to a value
!equal or smaller than zero.
!It is assumed that the basis of Kstart-1 functions has already
!been selected and the corresponding nonlinear parameters are
!stored in proper global arrays. Subroutines EnergyIA, EnergyIAM, and EnergyIB
!are called to evaluate the energy and its gradient.

!Arguments:
    integer,intent(in)     :: Kstart,Kstop,Kstep,NTrials,OptimizationType,MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
!Local variables:
    integer      i,j,K,AttemptToGetGoodFunc,ii,jj,jbest,ii2,jj2,jbest2,j2,q
    integer      np,npt,nfo,nfa,nfru,nfrup1,nvmax,nv
    integer      OpenFileErr,ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
    logical      IsSwapFileOK,IsEnergyImproved,ExitNeeded
    logical      IsOverlapBad,IsAnyLinCoeffBad,IsEnergyBad
    integer      wbfu_t,wmu_t,wbfu,wmu,rgm1_counter,rgm2_counter
    real(wp)  ms1,ms2
    real(wp)  Evalue,E_init,E_best
    real(wp)  t
    real(wp),allocatable,dimension(:,:)   :: ParSet,ParSetBest
    integer,allocatable,dimension(:,:)         :: IndSet,IndSetBest
    real(wp),allocatable,dimension(:)     :: x,x_best,grad
    real(wp),allocatable,dimension(:)     :: v_good
    integer,allocatable,dimension(:)         :: IndOptSequence
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V,V_init
    integer,parameter    :: LIV=60
    integer                 IV(LIV), IV_init(LIV)
    integer                 LV
    integer                 ALG
!Allocatable work space
    real(wp),allocatable,dimension(:)            :: WorkBuffReal
    integer,allocatable,dimension(:)                :: WorkBuffInt
    type(Glob_HistoryStep),allocatable,dimension(:) :: TempHistory
    real(wp),allocatable,dimension(:,:)          :: TempParam
    integer,allocatable,dimension(:,:)              :: TempInd
    integer,allocatable,dimension(:)                :: TempFunc
!====================================================
!These variables are used when a finite difference gradient is computed
    real(wp),allocatable,dimension(:)     ::    fx,fgrad,fgrad_5_points
    real(wp)                                    deltax,Evalue1,Evalue_plus_2h,Evalue_minus_2h
!====================================================

    wbfu_t=0
    wmu_t=0

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine BasisEnlI started'
      write(*,'(1x,a,1x,i0)') 'Kstart =',Kstart
      write(*,'(1x,a,1x,i0)') 'Kstop =',Kstop
      write(*,'(1x,a,1x,i0)') 'Kstep =',Kstep
      write(*,'(1x,a,1x,i0)') 'OptimizationType =',OptimizationType
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval =',MaxEnergyEval
    endif

!Setting the values of some global variables
    Glob_GSEPSolutionMethod='I'
    Glob_OverlapPenaltyAllowed=.false.
    Glob_nfa=Kstart+Kstep
    Glob_nfo=Kstep
    Glob_HSLeadDim=Kstop
    Glob_HSBuffLen=Kstop*Kstep
    np=Glob_np
    npt=Glob_npt
    nfo=Glob_nfo
    nfa=Glob_nfa
    nvmax=Kstep*Glob_npt

!Reallocate arrays that contain the information
!about basis functions and optimization process.
    call ReallocateBasisFuncData(Kstop,Glob_CurrBasisSize)

!Allocate some global arrays
    allocate(Glob_H(Kstop,Kstop))
    allocate(Glob_S(Kstop,Kstop))
    allocate(Glob_diagS(Kstop))
    allocate(Glob_invD(Kstop))
    allocate(Glob_c(Kstop))
    allocate(Glob_D(2*npt,Kstep,Kstop))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GSEPIIS, which is called
!inside EnergyIA, EnergyIAM, and EnergyIB
    allocate(Glob_WorkForGSEPIIS(Kstop))
    allocate(Glob_LastEigvector(Kstop))
    Glob_LastEigvector(1:Kstop)=ONE

!Allocate an array where a copy of the last eigenvector that was obtained in
!a successful energy evaluation is kept. Such a copy is needed because a failed
!inverse iteration process may leave Glob_LastEigvector in a state that makes
!it useless (or even harmful) as a starting vector for subsequent evaluations
    allocate(v_good(Kstop))
    v_good(1:Kstop)=ONE

!Allocate workspace for EnergyIB
    allocate(Glob_WkGR(Kstep*npt))

!Allocate workspace
    allocate(ParSet(npt,Kstep))
    allocate(ParSetBest(npt,Kstep))
    allocate(IndSet(Kstep,2))
    allocate(IndSetBest(Kstep,2))
    allocate(x(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))
    allocate(IndOptSequence(nfo))

!Allocate arrays used by DRMNG
    nvmax=npt*Kstep
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2 + 1
    allocate(V(LV))
    allocate(V_init(LV))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0 !set summary print format
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
    IV_init(1)=12 !DIVSET has been called and some default values were changed

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Calculating the initial energy
    ErrCode=0
    if (Kstart>1) then
      if (IsSwapFileOK) then
        if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
        Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.false.,ErrCode)
      else
        if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
        Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.true.,ErrCode)
      endif
    else
      Glob_CurrEnergy=huge(Glob_CurrEnergy)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0135 in BasisEnlI: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

    rgm1_counter=0
    rgm2_counter=0
    ms1=ZERO
    ms2=ZERO
    K=Kstart-1

!Main loop begins here
    do while (K<Kstop)
      if (K+Kstep<=Kstop) then
        nfo=Kstep
        nfru=K
        K=K+Kstep
      else
        nfo=Kstop-K
        nfru=K
        K=Kstop
      endif
      call linalg_setparam(K)  !reset linalg flags to account for changes in the basis size
      Glob_nfa=K
      Glob_nfru=nfru
      Glob_nfo=nfo
      nfrup1=nfru+1
      nv=nfo*npt
      E_init=Glob_CurrEnergy
      Glob_InvItTempCounter1=0
      Glob_InvItTempCounter2=0
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0)') 'Current basis size is',Glob_CurrBasisSize
        if (nfo>1) then
          write(*,'(1x,a,1x,i0,a,i0)') 'Selecting functions',nfrup1,'-',K
        else
          write(*,'(1x,a,1x,i0)') 'Selecting function',K
        endif
      endif
      IsOverlapBad=.true.
      IsAnyLinCoeffBad=.true.
      IsEnergyBad=.false.
      AttemptToGetGoodFunc=1
      do while ((IsOverlapBad.or.IsAnyLinCoeffBad.or.IsEnergyBad).and. &
                (AttemptToGetGoodFunc<=Glob_BadOverlapOrLinCoeffLim))
        !Random selection
        NumOfFailures=0
        IsEnergyImproved=.false.
        wbfu=0
        wmu=0
        do i=1,NTrials
          if (Glob_ProcID==0) call GenerateTrialParam(nfo,ParSet,IndSet,wbfu_t,wmu_t)
          call MPI_BCAST(ParSet,npt*nfo,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_BCAST(IndSet,2*nfo,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          Glob_NonlinParam(1:npt,nfrup1:K)=ParSet(1:npt,1:nfo)
          Glob_Index(nfrup1:K,1)=IndSet(1:nfo,1)
          Glob_Index(nfrup1:K,2)=IndSet(1:nfo,2)
          Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
          if (ErrCode==0) then
            v_good(1:K)=Glob_LastEigvector(1:K)
            if (Evalue<Glob_CurrEnergy) then
              Glob_CurrEnergy=Evalue
              ParSetBest(1:npt,1:nfo)=ParSet(1:npt,1:nfo)
              IndSetBest(1:nfo,1)=IndSet(1:nfo,1)
              IndSetBest(1:nfo,2)=IndSet(1:nfo,2)
              IsEnergyImproved=.true.
              wbfu=wbfu_t
              wmu=wmu_t
            endif
          else
            NumOfFailures=NumOfFailures+1
            !Restore the last good eigenvector as the vector left by the failed
            !inverse iteration process may be unusable as a starting vector
            Glob_LastEigvector(1:K)=v_good(1:K)
          endif
        enddo
        if (NumOfFailures*ONE/NTrials>Glob_MaxFracOfTrialFailsAllowed) then
          if (Glob_ProcID==0) then
            write(*,*) 'Error EC0136 in BasisEnlI: the number of eigenvalue problem solution failures'
            write(*,*) 'in random selection process exceeded limit'
            write(*,'(1x,a28,f7.3,a1)') 'The fraction of failures is ', &
              (100*NumOfFailures*ONE)/NTrials,'%'
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        if (.not.(IsEnergyImproved)) then
          if (Glob_ProcID==0) then
            write(*,*) 'Error EC0137 in BasisEnlI: random selection did not result'
            write(*,*) 'in any energy improvement'
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        Glob_NonlinParam(1:npt,nfrup1:K)=ParSetBest(1:npt,1:nfo)
        Glob_Index(nfrup1:K,1)=IndSetBest(1:nfo,1)
        Glob_Index(nfrup1:K,2)=IndSetBest(1:nfo,2)
        Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
        if (ErrCode==0) then
          Glob_CurrEnergy=Evalue
          v_good(1:K)=Glob_LastEigvector(1:K)
        else
          Glob_LastEigvector(1:K)=v_good(1:K)
        endif
        if (Glob_ProcID==0) then
          write (*,'(1x,a)',advance='no') 'E='
          call writereal(6,Glob_CurrEnergy)
          write (*,'(5x,a,1x,i0)') 'prototype function is',wbfu
          do i=1,nfo
            write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') nfru+i,':',IndSetBest(i,1),IndSetBest(i,2)
            call writerealarradv(6,ParSetBest(1:npt,i),npt)
          enddo
          write (*,*) 'Optimizing nonlinear parameters'
        endif

        !Optimization of the nonlinear parameters
        select case (OptimizationType)
        case(0) !No optimization
          do i=1,nfo
            x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          enddo
        case(1)

          !We first optimize indices (if they are not fixed)
          NumOfFailures=0
          !Generate a random sequence which will define the order in which
          !indices should be optimized (one index at a time)
          call GenerateRndIntSeq(nfo,IndOptSequence)
          !Loop where the indices are optimized. Note: there is some room for improvement
          !here it is programmed in a simple way when all matrix element of functions
          !nfrup1 through K are computed each time while it is not always necessary.
          do i=1,nfo
            ii=IndOptSequence(i) !!  ii should be the same for two indices
            j=Glob_Index(nfru+ii,1)
            j2=Glob_Index(nfru+ii,2)
            jbest=j
            jbest2=j2
            if (Glob_VectorCouplingScheme==2) then
              !Coupling scheme 2 requires the two indices to be equal. They therefore cannot
              !be varied independently and are varied together here. Note that if either of
              !the indices is fixed then both of them are effectively fixed, in which case
              !there is nothing to optimize.
              if (.not.(Glob_IsIndexFixed(1).or.Glob_IsIndexFixed(2))) then
                do jj=1,Glob_n
                  if (jj/=j) then
                    Glob_Index(nfru+ii,1)=jj
                    Glob_Index(nfru+ii,2)=jj
                    Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,1)=jbest
                      Glob_Index(nfru+ii,2)=jbest
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0138 in BasisEnlI: number of failures in energy calculations'
                          write(*,*) 'during the optimization of x,y-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest=jj
                      endif
                    endif
                  endif
                enddo
              endif
              jbest2=jbest
              Glob_Index(nfru+ii,1)=jbest
              Glob_Index(nfru+ii,2)=jbest2
            else
              !Coupling scheme 1 requires the two indices to be different, which is why the
              !value of the other index is excluded from the trial values below. Both index
              !orderings remain allowed, but a single-index move cannot turn one ordering
              !into the other. Scheme 0 imposes no restrictions on the indices.
              if (.not.Glob_IsIndexFixed(1)) then
                do jj=1,Glob_n
                  if ((jj/=j).and.((Glob_VectorCouplingScheme/=1).or.(jj/=j2))) then
                    Glob_Index(nfru+ii,1)=jj
                    Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,1)=jbest
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0138 in BasisEnlI: number of failures in energy calculations'
                          write(*,*) 'during the optimization of x-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest=jj
                      endif
                    endif
                  endif
                enddo
              endif
              Glob_Index(nfru+ii,1)=jbest
              if (.not.Glob_IsIndexFixed(2)) then
                do jj2=1,Glob_n
                  if ((jj2/=j2).and.((Glob_VectorCouplingScheme/=1).or.(jj2/=jbest))) then
                    Glob_Index(nfru+ii,2)=jj2
                    Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
                    if (ErrCode/=0) then
                      NumOfFailures=NumOfFailures+1
                      Glob_Index(nfru+ii,2)=jbest2
                      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
                        if (Glob_ProcID==0) then
                          write(*,*) 'Error EC0139 in BasisEnlI: number of failures in energy calculations'
                          write(*,*) 'during the optimization of y-indicies exceeded limit'
                        endif
                        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
                      endif
                    else
                      if (Evalue<Glob_CurrEnergy) then
                        Glob_CurrEnergy=Evalue
                        jbest2=jj2
                      endif
                    endif
                  endif
                enddo
              endif
              Glob_Index(nfru+ii,2)=jbest2
            endif
          enddo

          !Now we optimize nonlinear parameters

          !Setting IV and V values as was in their initial copies
          IV(1:LIV)=IV_init(1:LIV)
          V(1:LV)=V_init(1:LV)
          !  print*,'IA: this is IV(1) and IV(2) and IV(4))', IV(1),IV(2),IV(3),IV(4)
          nv=nfo*npt
          do i=1,nfo
            x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          enddo
          if (nfru>=nfo) then
            t=max(abs((E_init-Glob_CurrEnergy))/(abs(E_init)+abs(Glob_CurrEnergy)), &
                  10000*epsilon(Glob_CurrEnergy))
          else
            t=ONE
          endif
          !if (Glob_ProcID==0) write(*,*) 'scaling coeff=',t !remove later
          do i=1,nfo
            !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
            do j=1,npt
              !Make sure none of the D(i) will be zero or smaller than the threshold
              !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
              D(npt*(i-1)+j)=t
              !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
            enddo
          enddo

          ExitNeeded=.false.
          NumOfFailures=0
          NumOfEnergyEval=0
          NumOfGradEval=0
          if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
          E_best=Glob_CurrEnergy
          x_best(1:nfo*npt)=x(1:nfo*npt)

          do while (.not.(ExitNeeded))
            if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
            call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            select case (IV(1))
            case (1) !Only energy is needed
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
              NumOfEnergyEval=NumOfEnergyEval+1
              if (ErrCode/=0) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
                !Restore the last good eigenvector as the vector left by the failed
                !inverse iteration process may be unusable as a starting vector
                Glob_LastEigvector(1:K)=v_good(1:K)
              else
                Glob_CurrEnergy=Evalue
                if (Evalue<E_best) then
                  E_best=Evalue
                  x_best(1:nfo*npt)=x(1:nfo*npt)
                endif
              endif
            case (2) !Only gradient is needed
              call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
              do i=1,nfo
                Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
              enddo
              call EnergyIB(Evalue,grad,.true.,ErrCode)
              NumOfGradEval=NumOfGradEval+1
              if (ErrCode/=0) then
                NumOfFailures=NumOfFailures+1
                IV(2)=0
                !Restore the last good eigenvector as the vector left by the failed
                !inverse iteration process may be unusable as a starting vector
                Glob_LastEigvector(1:K)=v_good(1:K)
              else
                if (Evalue<E_best) then
                  E_best=Evalue
                  x_best(1:nfo*npt)=x(1:nfo*npt)
                endif
              endif
              !===================================
              !The lines below need to be uncommented when finite
              !difference gradient is shown
              !I compute also 2nd derivative, if +, then its minimum
!                  allocate(fx(nvmax))
!                  allocate(fgrad(nvmax))
!                  allocate(fgrad_5_points(nvmax))
!                  write(*,*)'             Analytic              finite diff gradient       5point finite diff gradient'
!                  do j=1,nv
!                      fx(1:nv)=x(1:nv)
!                      deltax=x(j)*1.0Q-3
!                      fx(j)=x(j)+deltax
!                       do i=1,nfo
!                          Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                      enddo
!                      Evalue1=EnergyIA(nfrup1,K,.true.,ErrCode)
!                      fx(j)=x(j)+2*deltax
!                       do i=1,nfo
!                          Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                      enddo
!                      Evalue_plus_2h=EnergyIA(nfrup1,K,.true.,ErrCode)
!                      fx(j)=x(j)-deltax
!                       do i=1,nfo
!                          Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                      enddo
!                      Evalue=EnergyIA(nfrup1,K,.true.,ErrCode)
!                      fx(j)=x(j)-2*deltax
!                       do i=1,nfo
!                           Glob_NonlinParam(1:npt,nfru+i)=fx((i-1)*npt+1:i*npt)
!                      enddo
!                      Evalue_minus_2h=EnergyIA(nfrup1,K,.true.,ErrCode)
!                      fgrad(j)=(Evalue1-Evalue)/(2*deltax)
!                      fgrad_5_points(j)=((Evalue_minus_2h-Evalue_plus_2h)/12+TWO*(Evalue1-Evalue)/THREE)/deltax
!                  enddo
!                  !f_2nd_der(j)=(4*(Evalue1+Evalue)/3-(Evalue_minus_2h+Evalue_plus_2h)/12-5*Glob_CurrEnergy/TWO)/(deltax*deltax)
!                  do j=1,nv
!                        write(*,*) j,'  ',grad(j),'  ' ,fgrad(j),'   ',fgrad_5_points(j)
!                  enddo
!                  write(*,*)
!                  do i=1,nfo
!                         Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
!                  enddo
!
!                  deallocate(fgrad)
!                  deallocate(fgrad_5_points)
!                  deallocate(fx)
              !===================================
            case (3:8) !Some kind of convergence has been reached
              ExitNeeded=.true.
            case (9:10) !Function evaluation limit has been reached.
              !This never suppose to happen because we
              !count the number of function evaluations ourselves.
              ExitNeeded=.true.
            endselect
            if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
              if (Glob_ProcID==0) then
                write(*,'(1x,a,1x,a,1x,a,1x,i0)') &
                  'Warning WC0117 in BasisEnlI: number of failures in energy or gradient', &
                  'calculations during the optimization of nonlinear parameters', &
                  'reached the limit of',Glob_MaxEnergyFailsAllowed
              endif
              !call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
            endif
            if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
          enddo
        endselect !(OptimizationType)

        !Calculate the energy and the linear coefficients at the best point found
        do i=1,nfo
          Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
        enddo
        IsEnergyBad=.false.
        Evalue=EnergyIAM(nfrup1,K,.true.,ErrCode)
        if (ErrCode==0) then
          Glob_CurrEnergy=Evalue
        else
          !Reject the generated basis function(s) and make another attempt to
          !generate them, just like it is done when the overlap or the linear
          !parameters are found to be bad
          IsEnergyBad=.true.
          Glob_CurrEnergy=E_init
          Glob_LastEigvector(1:K)=v_good(1:K)
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,a,1x,a)') &
              'Warning WC0118 in BasisEnlI: failed to evaluate energy after the optimization', &
              'of nonlinear parameters. Generated basis function(s) are rejected', &
              'and a new attempt to generate them will be made'
          endif
        endif

        !Checking if overlaps are OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad=.false.
        if (OverlapThreshold>ZERO) then
          ii=0
          do i=nfrup1,K
            do j=1,i-1
              if (abs(Glob_S(i,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                Glob_CurrEnergy=E_init
                if (Glob_ProcID==0) then
                  if (ii==1) then
                    write(*,*) 'Warning WC0115: overlap of the following functions exceeds threshold. ', &
                      'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                  endif
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') ii,':',i,j,'    S='
                  call writerealadv(6,Glob_S(i,j))
                endif
              endif
            enddo
          enddo
        endif
        !Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad=.false.
        if (LinCoeffThreshold>ZERO) then
          ii=0
          do i=1,K
            if(abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              Glob_CurrEnergy=E_init
              if (Glob_ProcID==0) then
                if (ii==1) then
                  write(*,*) 'Warning WC0116: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Generated basis function(s) are rejected and a new attempt to generate them will be made'
                endif
                write(*,'(1x,i6,a1,i6,a6)',advance='no') ii,':',i,'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif
        AttemptToGetGoodFunc=AttemptToGetGoodFunc+1

      enddo !(IsOverlapBad.or.IsAnyLinCoeffBad).and. &
      !(AttemptToGetGoodOverlap<=Glob_BasisEnlBadOverlapLim)

      if (Glob_ProcID==0) then
        write (*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
        write (*,*) 'E=',Glob_CurrEnergy
        do i=1,nfo
          write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') nfru+i,':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
          call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
        enddo
        write (*,*) 'Average number of iterations in GSEPIIS: ', &
          (Glob_InvItTempCounter2*ONE)/Glob_InvItTempCounter1
      endif

      !getting statistics about the distribution of generated parameters
      if (Glob_ProcID==0) then
        if (wmu==1) then
          rgm1_counter=rgm1_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i)) &
                   /Glob_NonlinParam(j,wbfu+i-1)
                ms1=ms1+abs(t)
              enddo
            enddo
          endif
        endif
        if (wmu==2) then
          rgm2_counter=rgm2_counter+1
          if (nfru>nfo) then
            do i=1,nfo
              do j=1,npt
                t=(Glob_NonlinParam(j,wbfu+i-1)-Glob_NonlinParam(j,nfru+i)) &
                   /Glob_NonlinParam(j,wbfu+i-1)
                ms2=ms2+abs(t)
              enddo
            enddo
          endif
        endif
      endif
      Glob_CurrBasisSize=K
      do i=1,nfo
        Glob_History(nfru+i)%Energy=Glob_CurrEnergy
        Glob_History(nfru+i)%CyclesDone=0
        Glob_History(nfru+i)%InitFuncAtLastStep=0
        Glob_History(nfru+i)%NumOfEnergyEvalDuringFullOpt=0
      enddo
      do i=1,nfo
        Glob_FuncNum(nfru+i)=nfru+i
      enddo

      if (Glob_ProcID==0) call SaveResults(Sort='no')
    enddo
!Main loop ends here

    call StoreMatricesInSwapFile()

!Deallocate arrays used by DRMNG
    deallocate(V_init)
    deallocate(V)
    deallocate(D)

!Deallocate workspace
    deallocate(IndOptSequence)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x)
    deallocate(IndSetBest)
    deallocate(IndSet)
    deallocate(ParSetBest)
    deallocate(ParSet)

!Deallocate workspace for EnergyIB
    deallocate(Glob_WkGR)

!Deallocate workspace for subroutine GSEPIIS, which is called
!inside EnergyIA, EnergyIAM, and EnergyIB
    deallocate(v_good)
    deallocate(Glob_LastEigvector)
    deallocate(Glob_WorkForGSEPIIS)

!Deallocate global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_D)
    deallocate(Glob_c)
    deallocate(Glob_invD)
    deallocate(Glob_diagS)
    deallocate(GLob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      write(*,*) 'Random selection statistics:'
      write(*,'(1x,a,1x,i0,1x,a)') 'Method 1 of generating basis functions was used',rgm1_counter,'times'
      if (rgm1_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms1/(npt*rgm1_counter)
      write(*,'(1x,a,1x,i0,1x,a)') 'Method 2 of generating basis functions was used',rgm2_counter,'times'
      if (rgm2_counter/=0) write(*,'(1x,a48,e13.6)') &
        'Average shift factor from prototype function is ',ms2/(npt*rgm2_counter)
      write(*,*)
      write(*,*) 'Routine BasisEnlI finished'
    endif

  end subroutine BasisEnlI

  subroutine OptCycleG(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                       NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold, SavingFreq)
!Subroutine OptCycleG performs an optimization of the basis set of K functions
!by means of cyclic optimization of a set of NumOfFuncToOpt functions at a time.
!After a step is made, the subroutine takes another set of NumOfFuncToOpt
!functions, whose numbers are shifted by NumFuncToShift (normally, one
!would probably want to use NumFuncToShift=NumOfFuncToOpt). This strategy
!is applied NumCycles times in a row to functions whose numbers range from
!FuncBegin to FuncEnd. MaxEnergyEval is the limit of the energy evaluations
!for each step. If it happens that at some step the basis functions being
!optimized become linearly dependent with other basis functions in the basis
!then such a change is rejected and this step is omitted. The pair linear
!dependency is checked by comparing the absolute value of the overlap
!with OverlapThreshold. If OverlapThreshold<=0 then no linear dependency
!check is performed. A Similar check is performed for all linear parameters.
!If adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions
!are rejected and the procedure is repeated until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avoid this check it is enough to set LinCoeffThreshold to a value
!equal or smaller than zero. Subroutines EnergyGA, EnergyGAM, and EnergyGB
!are called to evaluate the energy and its gradient.
!Parameter SavingFreq defines how often the results should be saved. If the
!value is equal to 1 than the results are saved after each cycle step. If the
!value is equal to 2, 3, ... k then the results are saved after every second,
!third, ... k-th cycle step.

!Arguments:
    integer,intent(in)     :: K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift
    integer,intent(in)     :: NumCycles, MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
    integer,intent(in)     :: SavingFreq
!Local variables:
    integer      i,j,m,ip,AttemptToGetGoodOverlap,ii,totsteps
    integer      np,npt,nfo,nfa,nfru,nfrup1,nv,nvmax,nfco,fbn,cbs,nfs
    integer      CurrCycle,CurrFunc,CurrFuncBegin
    integer      q,nr,OptIterCounter
    integer      ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
    integer      BlockSizeForDSYGVX
    real(wp)  Evalue,E_best,E_prev
    logical      IsSwapFileOK,ExitNeeded,LastIter
    logical      IsOverlapBad,IsAnyLinCoeffBad
    real(wp)  t
    real(wp),allocatable,dimension(:,:) :: NonlinParamTemp
    integer,allocatable,dimension(:)       :: FuncNumTemp
    real(wp),allocatable,dimension(:)   :: TempR
    real(wp),allocatable,dimension(:)   :: x,grad,x_init,x_best
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V,V_init
    integer,parameter    :: LIV=60
    integer                 IV(LIV),IV_init(LIV)
    integer                 LV
    integer                 ALG

!Setting the values of some global variables
    cbs=Glob_CurrBasisSize
    Glob_GSEPSolutionMethod='G'
    Glob_OverlapPenaltyAllowed=.false.
    np=Glob_np
    npt=Glob_npt
    nvmax=NumOfFuncToOpt*npt
    Glob_HSLeadDim=Glob_CurrBasisSize
    Glob_HSBuffLen=Glob_CurrBasisSize*NumOfFuncToOpt
    Glob_nfa=Glob_CurrBasisSize
    nfco=FuncEnd-FuncBegin+1
    fbn=FuncBegin+Glob_CurrBasisSize-FuncEnd

!Checking if cyclic optimization is already completed for this basis size
    if ((Glob_History(cbs)%CyclesDone>=NumCycles).and. &
        (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd)) then
      write(*,*)
      write(*,*) 'Routine OptCycleG started'
      write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
      write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin,'-',FuncEnd,'is already completed'
      write(*,*) 'Exiting OptCycleG...'
      write(*,*) 'Routine OptCycleG finished'
      return
    endif

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine OptCycleG started'
      write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
      write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin,'-',FuncEnd,'will be performed'
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval',MaxEnergyEval
    endif

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_D(2*npt,NumOfFuncToOpt,cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
    Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
    allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    allocate(Glob_IWorkForDSYGVX(5*cbs))

!Allocate workspace for EnergyGB
    allocate(Glob_WkGR(NumOfFuncToOpt*npt))

!Allocate workspace for SaveResults (we will use Sort='yes'
!option, which requires workspace)
    allocate(Glob_IntWorkArrForSaveResults(cbs))

!Allocate arrays used by DRMNG
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2 + 1
    allocate(V(LV))
    allocate(V_init(LV))

!Allocate workspace
    allocate(x(nvmax))
    allocate(x_init(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))
    allocate(NonlinParamTemp(npt,cbs-FuncBegin+1))
    allocate(FuncNumTemp(cbs-FuncBegin+1))
    allocate(TempR(cbs-FuncBegin+1))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0 !set summary print format
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
    IV_init(1)=12 !DIVSET has been called and some default values were changed

!Read swap file (if necessary) and distribute the data
    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Changing the order of the functions to be optimized and, if necessary,
!the corresponding matrix elements to reverse
    call ReverseFuncOrder(FuncBegin,FuncEnd)
    if (IsSwapFileOK) call ReverseMatElemOrder(FuncBegin,FuncEnd)

!Shifting the set of basis functions to be optimized to
!the very end and, if necessary, doing the permutation of matrix elements
!to reflect this change.
    call PermuteFunctions(FuncBegin,FuncEnd,FuncNumTemp,NonlinParamTemp)
    if (IsSwapFileOK) call PermuteMatrixElements(FuncBegin,FuncEnd,TempR)

!If the last optimized function number is greater than FuncEnd-1 or smaller than FuncBegin
!we need to change it so that the optimization begins from FuncBegin
    if (Glob_History(cbs)%InitFuncAtLastStep<FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    if (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd) then
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
    endif

!We need to make an initial permutation to shift the functions that were
!already optimized (if any) in the current optimization cycle
    ip=Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
    if (ip>0) then
      call PermuteFunctions(fbn,cbs-ip,FuncNumTemp,NonlinParamTemp)
      if (IsSwapFileOK) call PermuteMatrixElements(fbn,cbs-ip,TempR)
    endif

!Calculating the initial energy
    if (IsSwapFileOK) then
      !Getting initial energy
      if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyGA(1,cbs,.false.,ErrCode)
    else
      !Getting initial energy
      if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyGA(1,cbs,.true.,ErrCode)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0145 in OptCycleG: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

!Here comes main optimization cycle
    totsteps=0
    do CurrCycle=Glob_History(cbs)%CyclesDone+1,NumCycles
      !Doing cycle number CurrCycle
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0,1x,a)') 'Cycle',CurrCycle,'began'
      endif
      CurrFuncBegin=Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      q=FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
      OptIterCounter=0
      do CurrFunc=CurrFuncBegin,FuncEnd,NumOfFuncToShift
        !Note that CurrFunc counts functions as they were not shifted back
        !by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
        !initial input file.
        totsteps=totsteps+1
        nfo=min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
        Glob_nfo=nfo
        nv=nfo*npt
        nfru=cbs-nfo
        Glob_nfru=nfru
        nfrup1=nfru+1
        OptIterCounter=OptIterCounter+1
        if (Glob_ProcID==0) then
          write(*,*)
          if (nfo>1) then
            write(*,'(1x,a,1x,i0,a,i0)') 'Optimizing functions',CurrFunc,'-',CurrFunc+nfo-1
          else
            write(*,'(1x,a,1x,i0)') 'Optimizing function',CurrFunc
          endif
        endif

        !Permute basis functions so that the last nfo functions, which have just
        !been optimized are replaced with those that will be optimized at next step
        nfs=min(NumOfFuncToShift,nfo)
        if (OptIterCounter/=1) then
          i=cbs-nfo-nfs+1
          j=cbs-nfs
          call PermuteFunctions(i,j,FuncNumTemp,NonlinParamTemp)
          call PermuteMatrixElements(i,j,TempR)
          nr=NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)
          if (q-nfo>=nr*2) then
            !normal permutation, there is sufficient number of functions left
            m=cbs-nfo
            j=m-nr
            i=j-nr+1
          else
            !abnormal permutation, not sufficient number of functions left
            m=cbs-nfo
            j=m-nr
            i=cbs-q+1
          endif
          call PermuteFunctions2(i,j,m,FuncNumTemp,NonlinParamTemp)
          call PermuteMatrixElements2(i,j,m,TempR)
        endif

        if (Glob_ProcID==0) then
          if (Glob_AreParamPrintedInCycleOptX) then
            write (*,*) 'Nonlinear parameters before optimization:'
            do i=1,nfo
              write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') Glob_FuncNum(nfru+i),':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
              call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
            enddo
          endif
        endif

        !Setting IV and V values as was in their initial copies
        IV(1:LIV)=IV_init(1:LIV)
        V(1:LV)=V_init(1:LV)
        do i=1,nfo
          x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          x_init((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
        enddo

        t=max(ONE/(cbs*cbs*sqrt(ONE*cbs)),10000*epsilon(Glob_CurrEnergy))
        do i=1,nfo
          !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          do j=1,npt
            !Make sure none of the D(i) will be zero or smaller than the threshold
            !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j)=t
            !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          enddo
        enddo

        ExitNeeded=.false.
        NumOfFailures=0
        NumOfEnergyEval=0
        NumOfGradEval=0
        if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        E_best=Glob_CurrEnergy
        x_best(1:nfo*npt)=x(1:nfo*npt)
!Remember the energy that corresponds to the current (i.e. unchanged) values
!of the nonlinear parameters. It is needed in case the optimization of this
!function/set of functions has to be abandoned
        E_prev=Glob_CurrEnergy

        do while (.not.(ExitNeeded))
          if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
          call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          select case (IV(1))
          case (1) !Only energy is needed
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do i=1,nfo
              Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
            enddo
            Evalue=EnergyGA(nfrup1,cbs,.true.,ErrCode)
            NumOfEnergyEval=NumOfEnergyEval+1
            if (ErrCode/=0) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
            else
              Glob_CurrEnergy=Evalue
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nfo*npt)=x(1:nfo*npt)
              endif
            endif
          case (2) !Only gradient is needed
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do i=1,nfo
              Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
            enddo
            call EnergyGB(Evalue,grad,.true.,ErrCode)
            NumOfGradEval=NumOfGradEval+1
            if (ErrCode/=0) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
            else
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nfo*npt)=x(1:nfo*npt)
              endif
            endif
          case (3:8) !Some kind of convergence has been reached
            ExitNeeded=.true.
          case (9:10) !Function evaluation limit has been reached.
            !This is never supposed to happen because we
            !count the number of function evaluations ourselves.
            ExitNeeded=.true.
          endselect
          if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
            if (Glob_ProcID==0) then
              write(*,'(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0123 in OptCycleG: number of failures in energy or gradient', &
                'calculations during the optimization of nonlinear parameters', &
                'reached the limit of',Glob_MaxEnergyFailsAllowed
            endif
            !call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
          endif
          if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        enddo !while

        !Compute the energy and the linear coefficients at the best point found
        do i=1,nfo
          Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
        enddo
        Evalue=EnergyGAM(nfrup1,cbs,.true.,ErrCode)
    !!We run EnergyGA again because EnergyGAM might give slightly different
    !!energy than EnergyGA. EnergyGAM was needed to compute linear coefficients
        !Glob_CurrEnergy=EnergyGA(nfrup1,cbs,.false.,ErrCode)
        if (ErrCode==0) then
          Glob_CurrEnergy=Evalue
        else
          !Note that Glob_CurrEnergy is intentionally left unchanged here as
          !the value returned by a failed energy evaluation is meaningless
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,a,1x,a)') &
              'Warning WC0120 in OptCycleG: failed to evaluate energy after optimization', &
              'of nonlinear parameters. The values of the nonlinear parameters', &
              'will be left unchanged.'
          endif
        endif

        !Checking if overlap is OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad=.false.
        if (OverlapThreshold>ZERO) then
          ii=0
          do i=nfrup1,cbs
            do j=1,i-1
              if (abs(Glob_S(i,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                if (Glob_ProcID==0) then
                  if (ii==1) then
                    write(*,*) 'Warning WC0121: overlap of the following functions exceeds threshold. ', &
                      'Nonlinear parameters will be left unchanged'
                  endif
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') &
                    ii,':',Glob_FuncNum(i),Glob_FuncNum(j),'    S='
                  call writerealadv(6,Glob_S(i,j))
                endif
              endif
            enddo
          enddo
        endif
        !Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad=.false.
        if (LinCoeffThreshold>ZERO) then
          ii=0
          do i=1,cbs
            if(abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              if (Glob_ProcID==0) then
                if (ii==1) then
                  write(*,*) 'Warning WC0122: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                endif
                write(*,'(1x,i6,a1,i6,a6)',advance='no') ii,':',i,'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif

        if ((Glob_ProcID==0).and.(ErrCode==0).and.(.not.IsOverlapBad).and.(.not.IsAnyLinCoeffBad)) then
          write (*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
          write (*,*) 'E=',Glob_CurrEnergy
          if (Glob_AreParamPrintedInCycleOptX) then
            write (*,*) 'Nonlinear parameters after optimization:'
            do i=1,nfo
              write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') Glob_FuncNum(nfru+i),':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
              call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
            enddo
          endif
        endif

        if ((ErrCode/=0).or.IsOverlapBad.or.IsAnyLinCoeffBad) then
          !restore initial values of nonlinear parameters
          do i=1,nfo
            Glob_NonlinParam(1:npt,nfru+i)=x_init((i-1)*npt+1:i*npt)
          enddo
          Evalue=EnergyGA(nfrup1,cbs,.true.,ErrCode)
          if (ErrCode==0) then
            Glob_CurrEnergy=Evalue
          else
            !Leave this function (or set of functions) as it was and proceed to
            !the next one. The nonlinear parameters are the original ones and
            !the matrix elements that correspond to them have been computed,
            !so it is safe to continue
            if (Glob_ProcID==0) write(*,'(1x,a,1x,a)') &
              'Warning WC0124 in OptCycleG: energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_CurrEnergy=E_prev
          endif
        endif

        if (CurrFunc>FuncEnd-NumOfFuncToShift) then
          LastIter=.true.
        else
          LastIter=.false.
        endif

        Glob_History(cbs)%Energy=Glob_CurrEnergy
        if (LastIter) then
          Glob_History(cbs)%InitFuncAtLastStep=0
          Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
        else
          Glob_History(cbs)%InitFuncAtLastStep=CurrFunc
        endif

        if (Glob_ProcID==0) then
          if ((totsteps<=Glob_MinMandSavSteps).or.(mod(totsteps,SavingFreq)==0).or. &
              (CurrFunc+NumOfFuncToShift>=FuncEnd)) then
            call SaveResults(Sort='yes')
          endif
        endif

      enddo !end cycle CurrCycle

      if (Glob_ProcID==0) then
        write(*,*)
        write(*,*) 'Cycle',CurrCycle,' finished'
      endif
      if (CurrCycle/=NumCycles) then
        Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
        if (Glob_ProcID==0) write(*,'(1x,a47)',advance='no') &
          'Ordering basis functions and matrix elements...'
        call SortBasisFuncAndMatElem(fbn,cbs,FuncNumTemp,NonlinParamTemp,TempR)
        if (Glob_ProcID==0) write(*,*) 'done'
      endif

    enddo !End of main optimization cycle

    if (Glob_ProcID==0) write(*,'(1x,a53)',advance='no') &
      'Final ordering basis functions and matrix elements...'
    call SortBasisFuncAndMatElem(FuncBegin,cbs,FuncNumTemp,NonlinParamTemp,TempR)
    call ReverseFuncOrder(FuncBegin,cbs)
    call ReverseMatElemOrder(FuncBegin,cbs)
    if (Glob_ProcID==0) write(*,*) 'done'

    call StoreMatricesInSwapFile()

!deallocate workspace
    deallocate(TempR)
    deallocate(FuncNumTemp)
    deallocate(NonlinParamTemp)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x_init)
    deallocate(x)

!deallocate arrays used by DRMNG
    deallocate(V_init)
    deallocate(V)
    deallocate(D)

!deallocate workspace for SaveResults
    deallocate(Glob_IntWorkArrForSaveResults)

!deallocate workspace for EnergyGB
    deallocate(Glob_WkGR)

!Deallocate workspace for DSYGVX
    deallocate(Glob_IWorkForDSYGVX)
    deallocate(Glob_WorkForDSYGVX)

!Deallocate some global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_diagS)
    deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) write (*,*) 'Routine OptCycleG finished'

  end subroutine OptCycleG

  subroutine OptCycleI(K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift, &
                       NumCycles, MaxEnergyEval, OverlapThreshold, LinCoeffThreshold, SavingFreq)
!Subroutine OptCycleI performs an optimization of the basis set of K functions
!by means of cyclic optimization of a set of NumOfFuncToOpt functions at a time.
!After a step is made, the subroutine takes another set of NumOfFuncToOpt
!functions, whose numbers are shifted by NumFuncToShift (normally, one
!would probably want to use NumFuncToShift=NumOfFuncToOpt). This strategy
!is applied NumCycles times in a row to functions whose numbers range from
!FuncBegin to FuncEnd. MaxEnergyEval is the limit of the energy evaluations
!for each step. If it happens that at some step the basis functions being
!optimized become linearly dependent with other basis functions in the basis
!then such a change is rejected and this step is omitted. The pair linear
!dependency is checked by comparing the absolute value of the overlap
!with OverlapThreshold. If OverlapThreshold<=0 then no linear dependency
!check is performed. A similar check is performed for all linear parameters.
!If adding new function/functions makes any linear parameter
!of any function (not only those being added but any) is greater
!in magnitude than LinCoeffThreshold then such new function/functions
!are rejected and the procedure is repeated until a good set is generated
!or certain number of failures is reached (Glob_BadOverlapOrLinCoeffLim).
!To avoid this check it is enough to set LinCoeffThreshold to a value
!equal or smaller than zero. Subroutines EnergyIA, EnergyIAM, and EnergyIB
!are called to evaluate the energy and its gradient.
!Parameter SavingFreq defines how often the results should be saved. If the
!value is equal to 1 than the results are saved after each cycle step. If the
!value is equal to 2, 3, ... k then the results are saved after every second,
!third, ... k-th cycle step.

!Arguments:
    integer,intent(in)     :: K, FuncBegin, FuncEnd, NumOfFuncToOpt, NumOfFuncToShift
    integer,intent(in)     :: NumCycles, MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,LinCoeffThreshold
    integer,intent(in)     :: SavingFreq
!Local variables:
    integer      i,j,m,ip,AttemptToGetGoodOverlap,ii,totsteps
    integer      np,npt,nfo,nfa,nfru,nfrup1,nv,nvmax,nfco,fbn,cbs,nfs
    integer      CurrCycle,CurrFunc,CurrFuncBegin
    integer      q,nr,OptIterCounter
    integer      ErrCode,NumOfFailures,NumOfEnergyEval,NumOfGradEval
    real(wp)  Evalue,E_best,E_prev
    logical      IsSwapFileOK,ExitNeeded,LastIter
    logical      IsOverlapBad,IsAnyLinCoeffBad
    real(wp)  t
    real(wp),allocatable,dimension(:,:) :: NonlinParamTemp
    integer,allocatable,dimension(:)       :: FuncNumTemp
    real(wp),allocatable,dimension(:)   :: TempR
    real(wp),allocatable,dimension(:)   :: x,grad,x_init,x_best
    real(wp),allocatable,dimension(:)   :: v_good
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V,V_init
    integer,parameter    :: LIV=60
    integer                 IV(LIV),IV_init(LIV)
    integer                 LV
    integer                 ALG

!Setting the values of some global variables
    cbs=Glob_CurrBasisSize
    Glob_GSEPSolutionMethod='I'
    Glob_OverlapPenaltyAllowed=.false.
    np=Glob_np
    npt=Glob_npt
    nvmax=NumOfFuncToOpt*npt
    Glob_HSLeadDim=Glob_CurrBasisSize
    Glob_HSBuffLen=Glob_CurrBasisSize*NumOfFuncToOpt
    Glob_nfa=Glob_CurrBasisSize
    nfco=FuncEnd-FuncBegin+1
    fbn=FuncBegin+Glob_CurrBasisSize-FuncEnd

!Checking if cyclic optimization is already completed for this basis size
    if ((Glob_History(cbs)%CyclesDone>=NumCycles).and. &
        (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd)) then
      write(*,*)
      write(*,*) 'Routine OptCycleI started'
      write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
      write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin,'-',FuncEnd,'is already completed'
      write(*,*) 'Exiting OptCycleI...'
      write(*,*) 'Routine OptCycleI finished'
      return
    endif

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine OptCycleI started'
      write(*,'(1x,a,1x,i0)') 'Basis size is',cbs
      write(*,'(1x,a,1x,i0,a,i0,1x,a)') 'Cyclic optimization of basis functions', &
        FuncBegin,'-',FuncEnd,'will be performed'
      write(*,'(1x,a,1x,i0)') 'MaxEnergyEval',MaxEnergyEval
    endif

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_invD(cbs))
    allocate(Glob_D(2*npt,NumOfFuncToOpt,cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GSEPIIS, which is called
!inside EnergyIA, EnergyIAM, and EnergyIB
    allocate(Glob_WorkForGSEPIIS(cbs))
    allocate(Glob_LastEigvector(cbs))
    Glob_LastEigvector(1:cbs)=ONE

!Allocate an array where a copy of the last eigenvector that was obtained in
!a successful energy evaluation is kept. Such a copy is needed because a failed
!inverse iteration process may leave Glob_LastEigvector in a state that makes
!it useless (or even harmful) as a starting vector for subsequent evaluations
    allocate(v_good(cbs))

!Allocate workspace for EnergyIB
    allocate(Glob_WkGR(NumOfFuncToOpt*npt))

!Allocate workspace for SaveResults (we will use Sort='yes'
!option, which requires workspace)
    allocate(Glob_IntWorkArrForSaveResults(cbs))

!Allocate arrays used by DRMNG
    allocate(D(nvmax))
    LV=71+nvmax*(nvmax+13)/2 + 1
    allocate(V(LV))
    allocate(V_init(LV))

!Allocate workspace
    allocate(x(nvmax))
    allocate(x_init(nvmax))
    allocate(x_best(nvmax))
    allocate(grad(nvmax))
    allocate(NonlinParamTemp(npt,cbs-FuncBegin+1))
    allocate(FuncNumTemp(cbs-FuncBegin+1))
    allocate(TempR(cbs-FuncBegin+1))

!Setting some parameters for DRMNG
!We do it outside of the main loop so that no time
!is wasted for doing exactly the same operation over
!and over again

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    call DIVSET(ALG,IV_init,LIV,LV,V_init)
    IV_init(17)=1000000
    IV_init(18)=1000000
    IV_init(19)=0 !set summary print format
    IV_init(20)=0; IV_init(22)=0; IV_init(23)=-1; IV_init(24)=0
    V_init(31)=0.0_wp
    V_init(32)=2*epsilon(V_init(32))
    V_init(37)=2*epsilon(V_init(37))
!V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
!VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
!MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
    V_init(35)=Glob_MaxScStepAllowedInOpt*ONE
!V(35)=0.1*ONE
    IV_init(1)=12 !DIVSET has been called and some default values were changed

!Read swap file (if necessary) and distribute the data
    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Changing the order of the functions to be optimized and, if necessary,
!the corresponding matrix elements to reverse
    call ReverseFuncOrder(FuncBegin,FuncEnd)
    if (IsSwapFileOK) call ReverseMatElemOrder(FuncBegin,FuncEnd)

!Shifting the set of basis functions to be optimized to
!the very end and, if necessary, doing the permutation of matrix elements
!to reflect this change.
    call PermuteFunctions(FuncBegin,FuncEnd,FuncNumTemp,NonlinParamTemp)
    if (IsSwapFileOK) call PermuteMatrixElements(FuncBegin,FuncEnd,TempR)

!If the last optimized function number is greater than FuncEnd-1 or smaller than FuncBegin
!we need to change it so that the optimization begins from FuncBegin
    if (Glob_History(cbs)%InitFuncAtLastStep<FuncBegin) &
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
    if (Glob_History(cbs)%InitFuncAtLastStep>=FuncEnd) then
      Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
      Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
    endif

!We need to make an initial permutation to shift the functions that were
!already optimized (if any) in the current optimization cycle
    ip=Glob_History(cbs)%InitFuncAtLastStep-FuncBegin+NumOfFuncToShift
    if (ip>0) then
      call PermuteFunctions(fbn,cbs-ip,FuncNumTemp,NonlinParamTemp)
      if (IsSwapFileOK) call PermuteMatrixElements(fbn,cbs-ip,TempR)
    endif

!Calculating the initial energy
    if (IsSwapFileOK) then
      !Getting initial energy
      if (Glob_ProcID==0) write(*,*) 'Solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyIA(1,cbs,.false.,ErrCode)
    else
      !Getting initial energy
      if (Glob_ProcID==0) write(*,*) 'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyIA(1,cbs,.true.,ErrCode)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0150 in OptCycleI: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) write(*,*) 'Initial energy ',Glob_CurrEnergy

!Here comes main optimization cycle
    totsteps=0
    do CurrCycle=Glob_History(cbs)%CyclesDone+1,NumCycles
      !Doing cycle number CurrCycle
      if (Glob_ProcID==0) then
        write(*,*)
        write(*,'(1x,a,1x,i0,1x,a)') 'Cycle',CurrCycle,'began'
      endif
      CurrFuncBegin=Glob_History(cbs)%InitFuncAtLastStep+NumOfFuncToShift
      q=FuncEnd-Glob_History(cbs)%InitFuncAtLastStep-NumOfFuncToShift+1
      OptIterCounter=0
      do CurrFunc=CurrFuncBegin,FuncEnd,NumOfFuncToShift
        !Note that CurrFunc counts functions as they were not shifted back
        !by Glob_CurrBasisSize-FuncEnd, that is according to their numbering in the
        !initial input file.
        totsteps=totsteps+1
        nfo=min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
        Glob_nfo=nfo
        nv=nfo*npt
        nfru=cbs-nfo
        Glob_nfru=nfru
        nfrup1=nfru+1
        OptIterCounter=OptIterCounter+1
        Glob_InvItTempCounter1=0
        Glob_InvItTempCounter2=0
        if (Glob_ProcID==0) then
          write(*,*)
          if (nfo>1) then
            write(*,'(1x,a,1x,i0,a,i0)') 'Optimizing functions',CurrFunc,'-',CurrFunc+nfo-1
          else
            write(*,'(1x,a,1x,i0)') 'Optimizing function',CurrFunc
          endif
        endif

        !Permute basis functions so that the last nfo functions, which have just
        !been optimized are replaced with those that will be optimized at next step
        nfs=min(NumOfFuncToShift,nfo)
        if (OptIterCounter/=1) then
          i=cbs-nfo-nfs+1
          j=cbs-nfs
          call PermuteFunctions(i,j,FuncNumTemp,NonlinParamTemp)
          call PermuteMatrixElements(i,j,TempR)
          nr=NumOfFuncToShift*NumOfRowsToPermForUnitShift(OptIterCounter-1)
          if (q-nfo>=nr*2) then
            !normal permutation, there is sufficient number of functions left
            m=cbs-nfo
            j=m-nr
            i=j-nr+1
          else
            !abnormal permutation, not sufficient number of functions left
            m=cbs-nfo
            j=m-nr
            i=cbs-q+1
          endif
          call PermuteFunctions2(i,j,m,FuncNumTemp,NonlinParamTemp)
          call PermuteMatrixElements2(i,j,m,TempR)
          !call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
          Glob_CurrEnergy=EnergyIA(i,cbs,.false.,ErrCode)
          if (ErrCode/=0) then
            if (Glob_ProcID==0) write(*,*) &
              'Error EC0151 in OptCycleI: energy cannot be computed after permuting basis functions'
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
          endif
        endif

        if (Glob_ProcID==0) then
          if (Glob_AreParamPrintedInCycleOptX) then
            write (*,*) 'Nonlinear parameters before optimization:'
            do i=1,nfo
              write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') Glob_FuncNum(nfru+i),':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
              call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
            enddo
          endif
        endif

        !Setting IV and V values as was in their initial copies
        IV(1:LIV)=IV_init(1:LIV)
        V(1:LV)=V_init(1:LV)
        do i=1,nfo
          x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
          x_init((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,nfru+i)
        enddo

        t=max(ONE/(cbs*cbs*sqrt(ONE*cbs)),10000*epsilon(Glob_CurrEnergy))
        do i=1,nfo
          !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          do j=1,npt
            !Make sure none of the D(i) will be zero or smaller than the threshold
            !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j)=t
            !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          enddo
        enddo

        ExitNeeded=.false.
        NumOfFailures=0
        NumOfEnergyEval=0
        NumOfGradEval=0
        if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        E_best=Glob_CurrEnergy
        x_best(1:nfo*npt)=x(1:nfo*npt)
!Remember the energy and the eigenvector that correspond to the current
!(i.e. unchanged) values of the nonlinear parameters. They are needed in case
!the optimization of this function/set of functions has to be abandoned
        E_prev=Glob_CurrEnergy
        v_good(1:cbs)=Glob_LastEigvector(1:cbs)

        do while (.not.(ExitNeeded))
          if (Glob_ProcID==0) call DRMNG(D, Glob_CurrEnergy, grad, IV, LIV, LV, nv, V, x)
          call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          select case (IV(1))
          case (1) !Only energy is needed
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do i=1,nfo
              Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
            enddo
            Evalue=EnergyIA(nfrup1,cbs,.true.,ErrCode)
            NumOfEnergyEval=NumOfEnergyEval+1
            if (ErrCode/=0) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
              !Restore the last good eigenvector as the vector left by the failed
              !inverse iteration process may be unusable as a starting vector
              Glob_LastEigvector(1:cbs)=v_good(1:cbs)
            else
              Glob_CurrEnergy=Evalue
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nfo*npt)=x(1:nfo*npt)
              endif
            endif
          case (2) !Only gradient is needed
            call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
            do i=1,nfo
              Glob_NonlinParam(1:npt,nfru+i)=x((i-1)*npt+1:i*npt)
            enddo
            call EnergyIB(Evalue,grad,.true.,ErrCode)
            NumOfGradEval=NumOfGradEval+1
            if (ErrCode/=0) then
              NumOfFailures=NumOfFailures+1
              IV(2)=0
              !Restore the last good eigenvector as the vector left by the failed
              !inverse iteration process may be unusable as a starting vector
              Glob_LastEigvector(1:cbs)=v_good(1:cbs)
            else
              if (Evalue<E_best) then
                E_best=Evalue
                x_best(1:nfo*npt)=x(1:nfo*npt)
              endif
            endif
          case (3:8) !Some kind of convergence has been reached
            ExitNeeded=.true.
          case (9:10) !Function evaluation limit has been reached.
            !This is never supposed to happen because we
            !count the number of function evaluations ourselves.
            ExitNeeded=.true.
          endselect
          if (NumOfFailures==Glob_MaxEnergyFailsAllowed) then
            if (Glob_ProcID==0) then
              write(*,'(1x,a,1x,a,1x,a,1x,i0)') &
                'Warning WC0128 in OptCycleI: number of failures in energy or gradient', &
                'calculations during the optimization of nonlinear parameters', &
                'reached the limit of',Glob_MaxEnergyFailsAllowed
            endif
            !call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
          endif
          if (NumOfEnergyEval>=MaxEnergyEval) ExitNeeded=.true.
        enddo !while

        !Compute the energy and the linear coefficients at the best point found
        do i=1,nfo
          Glob_NonlinParam(1:npt,nfru+i)=x_best((i-1)*npt+1:i*npt)
        enddo
        Evalue=EnergyIAM(nfrup1,cbs,.true.,ErrCode)
        if (ErrCode==0) then
          Glob_CurrEnergy=Evalue
        else
          !Note that Glob_CurrEnergy is intentionally left unchanged here as
          !the value returned by a failed energy evaluation is meaningless
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,a,1x,a)') &
              'Warning WC0125 in OptCycleI: failed to evaluate energy after optimization', &
              'of nonlinear parameters. The values of the nonlinear parameters', &
              'will be left unchanged.'
          endif
          Glob_LastEigvector(1:cbs)=v_good(1:cbs)
        endif

        !Checking if overlap is OK (only in case OverlapThreshold>ZERO)
        IsOverlapBad=.false.
        if (OverlapThreshold>ZERO) then
          ii=0
          do i=nfrup1,cbs
            do j=1,i-1
              if (abs(Glob_S(i,j))>OverlapThreshold) then
                ii=ii+1
                IsOverlapBad=.true.
                if (Glob_ProcID==0) then
                  if (ii==1) then
                    write(*,*) 'Warning WC0126: overlap of the following functions exceeds threshold. ', &
                      'Nonlinear parameters will be left unchanged'
                  endif
                  write(*,'(1x,i6,a1,i6,i6,a6)',advance='no') &
                    ii,':',Glob_FuncNum(i),Glob_FuncNum(j),'    S='
                  call writerealadv(6,Glob_S(i,j))
                endif
              endif
            enddo
          enddo
        endif
        !Checking if linear coefficients are OK (only in case LinCoeffThreshold>ZERO)
        IsAnyLinCoeffBad=.false.
        if (LinCoeffThreshold>ZERO) then
          ii=0
          do i=1,cbs
            if(abs(Glob_c(i))>LinCoeffThreshold) then
              ii=ii+1
              IsAnyLinCoeffBad=.true.
              if (Glob_ProcID==0) then
                if (ii==1) then
                  write(*,*) 'Warning WC0127: absolute value of linear parameters of the following functions exceeds threshold. ', &
                    'Nonlinear parameters will be left unchanged'
                endif
                write(*,'(1x,i6,a1,i6,a6)',advance='no') ii,':',i,'    c='
                call writerealadv(6,Glob_c(i))
              endif
            endif
          enddo
        endif

        if ((Glob_ProcID==0).and.(ErrCode==0).and.(.not.IsOverlapBad).and.(.not.IsAnyLinCoeffBad)) then
          write (*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
          write (*,*) 'E=',Glob_CurrEnergy
          if (Glob_AreParamPrintedInCycleOptX) then
            write (*,*) 'Nonlinear parameters after optimization:'
            do i=1,nfo
              write(*,'(1x,i6,a1,i6,1x,i6)',advance='no') Glob_FuncNum(nfru+i),':',Glob_Index(nfru+i,1),Glob_Index(nfru+i,2)
              call writerealarradv(6,Glob_NonlinParam(1:npt,nfru+i),npt)
            enddo
            write (*,'(1x,a41,f8.4)') 'Average number of iterations in GSEPIIS: ', &
              (Glob_InvItTempCounter2*ONE)/Glob_InvItTempCounter1
          endif
        endif

        if ((ErrCode/=0).or.IsOverlapBad.or.IsAnyLinCoeffBad) then
          !restore initial values of nonlinear parameters
          do i=1,nfo
            Glob_NonlinParam(1:npt,nfru+i)=x_init((i-1)*npt+1:i*npt)
          enddo
          Evalue=EnergyIA(nfrup1,cbs,.true.,ErrCode)
          if (ErrCode/=0) then
            !The matrix elements are correct, so the failure can only come from
            !the inverse iteration process itself. Restore the last good
            !eigenvector and try once more. No matrix elements need be recomputed
            Glob_LastEigvector(1:cbs)=v_good(1:cbs)
            Evalue=EnergyIA(nfrup1,cbs,.false.,ErrCode)
          endif
          if (ErrCode==0) then
            Glob_CurrEnergy=Evalue
          else
            !Leave this function (or set of functions) as it was and proceed to
            !the next one. The nonlinear parameters are the original ones and
            !the matrix elements that correspond to them have been computed. The
            !factorization of H-Glob_ApproxEnergy*S is redone at the next step,
            !so it is safe to continue
            if (Glob_ProcID==0) write(*,'(1x,a,1x,a)') &
              'Warning WC0129 in OptCycleI: energy cannot be computed.', &
              'Proceeding to the next basis function'
            Glob_LastEigvector(1:cbs)=v_good(1:cbs)
            Glob_CurrEnergy=E_prev
          endif
        endif

        if (CurrFunc>FuncEnd-NumOfFuncToShift) then
          LastIter=.true.
        else
          LastIter=.false.
        endif

        Glob_History(cbs)%Energy=Glob_CurrEnergy
        if (LastIter) then
          Glob_History(cbs)%InitFuncAtLastStep=0
          Glob_History(cbs)%CyclesDone=Glob_History(cbs)%CyclesDone+1
        else
          Glob_History(cbs)%InitFuncAtLastStep=CurrFunc
        endif

        if (Glob_ProcID==0) then
          if ((totsteps<=Glob_MinMandSavSteps).or.(mod(totsteps,SavingFreq)==0).or. &
              (CurrFunc+NumOfFuncToShift>=FuncEnd)) then
            call SaveResults(Sort='yes')
          endif
        endif

      enddo !end cycle CurrCycle

      if (Glob_ProcID==0) then
        write(*,*)
        write(*,*) 'Cycle',CurrCycle,' finished'
      endif
      if (CurrCycle/=NumCycles) then
        Glob_History(cbs)%InitFuncAtLastStep=FuncBegin-NumOfFuncToShift
        if (Glob_ProcID==0) write(*,'(1x,a47)',advance='no') &
          'Ordering basis functions and matrix elements...'
        call SortBasisFuncAndMatElem(fbn,cbs,FuncNumTemp,NonlinParamTemp,TempR)
        !call EnergyIA because we need to refactorize Glob_H-Glob_ApproxEnergy*Glob_S
        Glob_CurrEnergy=EnergyIA(fbn,cbs,.false.,ErrCode)
        if (ErrCode/=0) then
          if (Glob_ProcID==0) write(*,*) &
            'Error EC0154 in OptCycleI: energy cannot be computed after sorting basis functions'
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        if (Glob_ProcID==0) write(*,*) 'done'
      endif

    enddo !End of main optimization cycle

    if (Glob_ProcID==0) write(*,'(1x,a53)',advance='no') &
      'Final ordering basis functions and matrix elements...'
    call SortBasisFuncAndMatElem(FuncBegin,cbs,FuncNumTemp,NonlinParamTemp,TempR)
    call ReverseFuncOrder(FuncBegin,cbs)
    call ReverseMatElemOrder(FuncBegin,cbs)
    if (Glob_ProcID==0) write(*,*) 'done'

    call StoreMatricesInSwapFile()

!deallocate workspace
    deallocate(TempR)
    deallocate(FuncNumTemp)
    deallocate(NonlinParamTemp)
    deallocate(grad)
    deallocate(x_best)
    deallocate(x_init)
    deallocate(x)

!deallocate arrays used by DRMNG
    deallocate(V_init)
    deallocate(V)
    deallocate(D)

!deallocate workspace for SaveResults
    deallocate(Glob_IntWorkArrForSaveResults)

!deallocate workspace for EnergyGB
    deallocate(Glob_WkGR)

!Deallocate workspace for GSEPIIS
    deallocate(Glob_WorkForGSEPIIS)
    deallocate(v_good)
    deallocate(Glob_LastEigvector)

!Deallocate some global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_invD)
    deallocate(Glob_diagS)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) write (*,*) 'Routine OptCycleI finished'

  end subroutine OptCycleI

  subroutine FullOpt1G(InitFunc,FinalFunc,MaxEnergyEval,OverlapThreshold,MaxOverlapPenalty,  &
                       DataSaveMinTimeInterv,HessianSaveMinTimeInterv,HessFileName)
!Subroutine FullOpt1G performs simultaneous optimization of
!nonlinear parameters of basis functions whose number
!ranges from InitFunct to FinalFunc. Optimization routine
!used in FullOpt1G is DRMNG. Parameter MaxEnergyEval defines
!the limit of the energy evaluations (just the energies, not
!the gradient) in the optimization process. The optimization uses penalty
!functions for the overlaps. If an overlap, Sij, exceeds the threshold,
!OverlapThreshold, then a penalty is applied for this (the penalty function is
!quadratic and smooth so that the actual objective function is differentiable,
!although the second derivative is not smooth). If the value of OverlapThreshold
!is set equal or higher than 1.0 then no penalty will be applied. Parameter
!MaxOverlapPenalty defines the magnitude of penalty. When MaxOverlapPenalty=1.0
!and Sij=1.0 then the penalty due to this pair overlap will be 1.0. The routine
!saves the hessian in the file whose name is defined by
!variable HessFileName (may include the full path if necessary).
!This save occurs periodically, with the interval approximately equal
!to the value of parameter HessianSaveMinTimeInterv.
!Similarly, the input/output file (which contains basis set) is
!updated periodically, with the interval approximately equal
!to the value of parameter DataSaveMinTimeInterv.
!If there is no need to save the hessian at all, one can set
!HessFileName=' ' or HessFileName='none' or HessFileName='NONE'.
!Any other name will be treated as a file name. It is important to
!mention that the hessian itself and the corresponding file
!may require a lot of memory/storage (approximately k*k/2 real values, where
!k is the dimensionality of the problem; k=CurrBasisSize*Glob_npt).
!Once the routine starts, it actually tries
!to locate file HessFileName and read the hessian from there. If this
!operation is successful, this hessian is used in the initialization of
!routine DRMNG. Subroutines EnergyGA,  and EnergyGB
!are called to evaluate the energy and its gradient.

!Arguments:
    integer,intent(in)     :: InitFunc,FinalFunc,MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,MaxOverlapPenalty
    real(4),intent(in)     :: DataSaveMinTimeInterv,HessianSaveMinTimeInterv
    character(Glob_FileNameLength),intent(in) :: HessFileName
!Local variables:
    integer      i,j,q
    integer      np,npt,nfo,nfa,nfru,nv
    integer      OpenFileErr,ErrCode,NumOfEnergyEval,NumOfGradEval,NumOfFailures
    integer      NumOfEnergyEvalDuringFullOpt_Init
    integer      BlockSizeForDSYGVX
    logical      IsSwapFileOK,ExitNeeded
    logical      SaveHessian,IsHessFileOK,IsHessSaveSuccess
    real(wp)  Evalue,CurrentEnergy
    real(4)      TimeOfLastSave, TimeOfLastHessSave
    real(wp)  t
    integer      tas,InitFuncNew
    real(wp)  MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap
    real(wp),allocatable,dimension(:,:)   :: NonlinParamTemp
    integer,allocatable,dimension(:)         :: FuncNumTemp
    real(wp),allocatable,dimension(:)     :: TempR
    real(wp),allocatable,dimension(:)     :: x,grad
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V
    integer,parameter    :: LIV=60
    integer                 IV(LIV)
    integer                 LV
    integer                 ALG
    integer                 IVLMAT
!!====================================================
!!These variables are used when a finite difference gradient is computed
    real(wp)                              deltax,Evalue1,Evalue_plus_2h,Evalue_minus_2h,res_5_point
!!====================================================

    IsHessFileOK=.false.

    if (OverlapThreshold>=ONE) then
      Glob_OverlapPenaltyAllowed=.false.
    else
      Glob_OverlapPenaltyAllowed=.true.
      Glob_OverlapPenaltyThreshold2=OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty=MaxOverlapPenalty
    endif
    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine FullOpt1G started'
      write(*,*) 'Simultaneous optimization of nonlinear parameters of basis functions'
      write(*,*) InitFunc,'  through',FinalFunc,'  will be attempted'
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Overlap threshold is ',abs(OverlapThreshold)
        write(*,*) 'Max value of a pair overlap penalty is ',Glob_MaxOverlapPenalty
        write(*,*) 'Warning! The energy value that will be shown during the optimization'
        write(*,*) 'may differ from the actual energy'
      else
        write(*,*) 'No constraints on overlaps will be imposed'
      endif
    endif

!Setting the values of some global variables
    Glob_GSEPSolutionMethod='G'
    Glob_nfa=Glob_CurrBasisSize
    Glob_nfo=FinalFunc-InitFunc+1
    Glob_nfru=Glob_CurrBasisSize-Glob_nfo
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    nfa=Glob_nfa
    nfo=Glob_nfo
    nfru=Glob_nfru
    nv=nfo*npt
    InitFuncNew=nfa+InitFunc-FinalFunc
    Glob_HSBuffLen=max(min(nfa*(nfa+1)/2,1000),30*nfa)

!Allocate some global arrays
    allocate(Glob_H(nfa,nfa))
    allocate(Glob_S(nfa,nfa))
    allocate(Glob_diagH(nfa))
    allocate(Glob_diagS(nfa))
    allocate(Glob_D(2*npt,nfo,nfa))
    allocate(Glob_c(nfa))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',nfa,nfa,nfa,nfa)
    Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*nfa,8*nfa)
    allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
    allocate(Glob_IWorkForDSYGVX(5*nfa))

!Allocate workspace for EnergyGB
    allocate(Glob_WkGR(nfo*npt))

!Allocate arrays used by DRMNG
    LV=71 + nv*(nv+13)/2 + 1
    if (Glob_ProcID==0) then
      allocate(D(nv))
      allocate(V(LV))
    endif

!Allocate workspace
    allocate(x(nv))
    allocate(grad(nv))

!Setting up a logical variable that determines
!whether the hessian should be saved from time to time
    if ((HessFileName==' ').or.(HessFileName=='none').or. &
        (HessFileName=='NONE').or.(HessFileName=='None')) then
      SaveHessian=.false.
    else
      SaveHessian=.true.
    endif

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Shifting basis functions InitFunc through FinalFunc to
!the very end and, if necessary, doing the permutation of
!matrix elements to reflect this change in function order.
    if (FinalFunc/=nfa) then
      !allocate space
      tas=FinalFunc-InitFunc+1
      allocate(NonlinParamTemp(1:npt,1:tas))
      allocate(FuncNumTemp(1:tas))
      allocate(TempR(1:tas))
      call PermuteFunctions(InitFunc,FinalFunc,FuncNumTemp,NonlinParamTemp)
      if (IsSwapFileOK) call PermuteMatrixElements(InitFunc,FinalFunc,TempR)
      deallocate(TempR)
      deallocate(FuncNumTemp)
      deallocate(NonlinParamTemp)
      !Allocate workspace for SaveResults (it uses Sort='yes'
      !option, which requires workspace)
      allocate(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
    endif

!Calculating the initial energy
    if (IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.false.,ErrCode)
    else
      if (Glob_ProcID==0) write(*,'(1x,a59)',advance='no') &
        'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyGA(1,Glob_CurrBasisSize,.true.,ErrCode)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0160 in FullOpt1G: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    if (Glob_ProcID==0) write(*,*) ' done'
    if (Glob_ProcID==0) then
      call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                           ',Glob_TotalOverlapPenalty
        write(*,*) 'Initial energy (including overlap penalty ',Glob_CurrEnergy
      else
        write(*,*) 'Initial energy                            ',Glob_CurrEnergy
      endif
      write(*,*)   'Maximal overlap                           ',MaxAbsOverlap
      write(*,*)   'Minimal overlap                           ',MinAbsOverlap
      write(*,*)   'Average abs value of overlap              ',AverageAbsOverlap
    endif

    call CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave=Glob_TimeSinceStart
    TimeOfLastHessSave=Glob_TimeSinceStart

!Setting parameters for DRMNG

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    if (Glob_ProcID==0) then
      call DIVSET(ALG,IV,LIV,LV,V)
      IV(17)=1000000; IV(18)=1000000
      IV(19)=-1 !set summary print format
      IV(20)=0; IV(22)=0; IV(23)=-1; IV(24)=0
      V(31)=0.0_wp
      V(32)=2*epsilon(V(32))
      V(37)=2*epsilon(V(37))
      !V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
      !VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
      !MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
      V(35)=Glob_MaxScStepAllowedInOpt
      IV(1)=12 !DIVSET has been called and some default values were changed
      nv=nfo*npt
    endif
    do i=1,nfo
      x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,InitFuncNew+i-1)
    enddo

    if (Glob_ProcID==0) then
      !If SaveHessian=.true. then try to read the Hessian
      !from the file
      if (SaveHessian) then
        IVLMAT=IV(42)
        call ReadHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessFileOK)
        if (IsHessFileOK) IV(25)=0
      endif
      if ((.not.Glob_FullOptSaveD).or.(.not.IsHessFileOK).or.(.not.SaveHessian)) then
!    !Set the scaling vector
!    do i=1,nfo
!      t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
!      do j=1,np
!        !Make sure none of the D(i) will be zero or smaller than the threshold
!        D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
!      enddo
!    enddo
        !Set the scaling vector
        t=max(ONE/(nfa*nfa*sqrt(ONE*nfa)),10000*epsilon(Glob_CurrEnergy))
        do i=1,nfo
          !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          do j=1,npt
            !Make sure none of the D(i) will be zero or smaller than the threshold
            !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j)=t
            !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          enddo
        enddo
      endif
    endif

    ExitNeeded=.false.
    NumOfFailures=0
    NumOfEnergyEval=0
    NumOfGradEval=0
    NumOfEnergyEvalDuringFullOpt_Init=Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt

    do while (.not.(ExitNeeded))
      if (Glob_ProcID==0) call DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
      call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      select case (IV(1))
      case (1) !Only energy is needed
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
        enddo
        Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
        NumOfEnergyEval=NumOfEnergyEval+1
        if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
          IV(2)=1
        else
          CurrentEnergy=Evalue
        endif
        if (Glob_ProcID==0) then
          if (Evalue<Glob_CurrEnergy) then
            !Save data and Hessian if necessary
            call CPU_TIME(Glob_TimeSinceStart)
            if (Glob_TimeSinceStart-TimeOfLastSave>DataSaveMinTimeInterv) then
              !Save the results if more than DataSaveMinTimeInterv seconds
              !have passed since the last save
              if (Glob_OverlapPenaltyAllowed) then
                Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
              else
                Glob_CurrEnergy=Evalue
              endif
              !Changing history
              Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
              Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              if (FinalFunc==nfa) then
                call SaveResults(Sort='no')
              else
                call SaveResults(Sort='yes')
              endif
              write(*,*) 'Data file has been updated'
              call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
              write(*,*) 'Some current statistics:'
              if (Glob_OverlapPenaltyAllowed) then
                write(*,*) 'Energy (without overlap penalty)  ',Evalue-Glob_TotalOverlapPenalty
                write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
                write(*,*) 'Energy (including overlap penalty)',Evalue
              else
                write(*,*) 'Energy                            ',Evalue
              endif
              write(*,*)   'Maximal overlap                   ',MaxAbsOverlap
              write(*,*)   'Minimal overlap                   ',MinAbsOverlap
              write(*,*)   'Average abs value of overlap      ',AverageAbsOverlap
              TimeOfLastSave=Glob_TimeSinceStart
            endif
            if ((Glob_TimeSinceStart-TimeOfLastHessSave>HessianSaveMinTimeInterv) &
                .and.(SaveHessian)) then
              !Save the Hessian if more than HessianSaveMinTimeInterv seconds
              !have passed since the last save
              call SaveHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessSaveSuccess)
              TimeOfLastHessSave=Glob_TimeSinceStart
            endif
          endif
        endif
      case (2) !Only gradient is needed
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
        enddo
        call EnergyGB(Evalue,grad,.true.,ErrCode)
        NumOfGradEval=NumOfGradEval+1
        if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
          IV(2)=0
        endif
!        !===================================
!        !The lines below need to be uncommented when finite
!    !difference gradient is shown
!    if (Glob_ProcID==0) write(*,*) '    analytic gradient       finite diff gradient      5points finite dif gradient'
!      do i=1,nfo
!        do j=1,npt
!                 deltax=Glob_NonlinParam(j,i+nfru)*1.0Q-6
!                 Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)+deltax
!                 Evalue1=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!                 Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)+deltax
!                 Evalue_plus_2h=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!                 Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)-3*deltax
!                 Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!                 Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)-deltax
!                 Evalue_minus_2h=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
!                 res_5_point= (Evalue_minus_2h-Evalue_plus_2h)/12+TWO*(Evalue1-Evalue)/THREE
!                !res_2nd_der= -(Evalue_minus_2h+Evalue_plus_2h)/12+FOUR*(Evalue1+Evalue)/THREE-FIVE*Glob_CurrEnergy/TWO
!                 if (Glob_ProcID==0) write(*,'(1x,i3,1x,e23.16,1x,e23.16,1x,e23.16)')&
!                   (i-1)*npt+j,grad((i-1)*npt+j),(Evalue1-Evalue)/(2*deltax), res_5_point/deltax
!                   Glob_NonlinParam(j,i+nfru)=Glob_NonlinParam(j,i+nfru)+2*deltax
!            enddo
!        if (Glob_ProcID==0) write(*,*)
!      enddo
!
!    !===================================
      case (3:8) !Some kind of convergence has been reached
        ExitNeeded=.true.
      case (9:10) !Function evaluation limit has been reached.
        !This never suppose to happen because we
        !count the number of function evaluations ourselves.
        ExitNeeded=.true.
      endselect
      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,a,1x,a)') &
            'Error EC0161 in FullOpt1G: number of failures in energy or gradient', &
            'calculations during the optimization of nonlinear parameters', &
            'exceeded limit'
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
      if (NumOfEnergyEval>=MaxEnergyEval) then
        if (Glob_ProcID==0) then
          write(*,*) 'Warning WC0130 in FullOpt1G: number of energy evaluations reached limit'
          write(*,*) 'Optimization is terminated'
        endif
        ExitNeeded=.true.
      endif
    enddo

!Calculate the energy at the best point found
    do i=1,nfo
      Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
    enddo
    Evalue=EnergyGA(InitFuncNew,nfa,.true.,ErrCode)
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0162 in FullOpt1G: failed to evaluate energy after the optimization'
        write(*,*) 'of nonlinear parameters'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    if (Glob_OverlapPenaltyAllowed) then
      Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
    else
      Glob_CurrEnergy=Evalue
    endif

!Printing the number of energy/gradient evaluations
!and the energy after the optimization
    if (Glob_ProcID==0) then
      write(*,*)
      write(*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
      write(*,*) 'Final energy and overlap statistics:'
      call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Energy (without overlap penalty)  ',Evalue-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
        write(*,*) 'Energy (including overlap penalty)',Evalue
      else
        write(*,*) 'Energy                            ',Evalue
      endif
      write(*,*)   'Maximal overlap                   ',MaxAbsOverlap
      write(*,*)   'Minimal overlap                   ',MinAbsOverlap
      write(*,*)   'Average abs value of overlap      ',AverageAbsOverlap
    endif

!Adding data to history
    Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
    Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval

!Shifting basis functions InitFunc through FinalFunc back from
!the very end to the middle, where they initially were, and, if
!necessary, doing the permutation of matrix elements to reflect
!this change in function order.
    if (FinalFunc/=nfa) then
      !allocate space
      tas=nfa-FinalFunc
      allocate(NonlinParamTemp(1:npt,1:tas))
      allocate(FuncNumTemp(1:tas))
      allocate(TempR(1:tas))
      call PermuteFunctions(InitFunc,InitFunc+tas-1,FuncNumTemp,NonlinParamTemp)
      if (Glob_UseSwapFile) call PermuteMatrixElements(InitFunc,InitFunc+tas-1,TempR)
      deallocate(TempR)
      deallocate(FuncNumTemp)
      deallocate(NonlinParamTemp)
      !deallocate workspace for SaveResults
      deallocate(Glob_IntWorkArrForSaveResults)
    endif

!saving results
    if (Glob_ProcID==0) call SaveResults(Sort='no')

    call StoreMatricesInSwapFile()

    if (Glob_OverlapPenaltyAllowed) Glob_OverlapPenaltyAllowed=.false.

!deallocate workspace
    deallocate(grad)
    deallocate(x)

!deallocate arrays used by DRMNG
    if (Glob_ProcID==0) then
      deallocate(V)
      deallocate(D)
    endif
!deallocate workspace for EnergyGB
    deallocate(Glob_WkGR)

!Deallocate workspace for DSYGVX
    deallocate(Glob_IWorkForDSYGVX)
    deallocate(Glob_WorkForDSYGVX)

!Deallocate some global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_diagS)
    deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) write (*,*) 'Routine FullOpt1G finished'

  end subroutine FullOpt1G

  subroutine FullOpt1I(InitFunc,FinalFunc,MaxEnergyEval,OverlapThreshold,MaxOverlapPenalty, &
                       DataSaveMinTimeInterv,HessianSaveMinTimeInterv,HessFileName)
!Subroutine FullOpt1I performs simultaneous optimization of
!nonlinear parameters of basis functions whose number
!ranges from InitFunct to FinalFunc. Optimization routine
!used in FullOpt1I is DRMNG. Parameter MaxEnergyEval defines
!the limit of the energy evaluations (just the energies, not
!the gradient) in the optimization process. The optimization uses penalty
!functions for the overlaps. If an overlap, Sij, exceeds the threshold,
!OverlapThreshold, then a penalty is applied for this (the penalty function is
!quadratic and smooth so that the actual objective function is differentiable,
!although the second derivative is not smooth). If the value of OverlapThreshold
!is set equal or higher than 1.0 then no penalty will be applied. Parameter
!MaxOverlapPenalty defines the magnitude of penalty. When MaxOverlapPenalty=1.0
!and Sij=1.0 then the penalty due to this pair overlap will be 1.0. The routine
!saves the hessian in the file whose name is defined by
!variable HessFileName (may include the full path if necessary).
!This save occurs periodically, with the interval approximately equal
!to the value of parameter HessianSaveMinTimeInterv.
!Similarly, the input/output file (which contains basis set) is
!updated periodically, with the interval approximately equal
!to the value of parameter DataSaveMinTimeInterv.
!If there is no need to save the hessian at all, one can set
!HessFileName=' ' or HessFileName='none' or HessFileName='NONE'.
!Any other name will be treated as a file name. It is important to
!mention that the hessian itself and the corresponding file
!may require a lot of memory/storage (approximately k*k/2 real values, where
!k is the dimensionality of the problem; k=CurrBasisSize*Glob_npt).
!Once the routine starts, it actually tries
!to locate file HessFileName and read the hessian from there. If this
!operation is successful, this hessian is used in the initialization of
!routine DRMNG. Subroutines EnergyIA,  and EnergyIB
!are called to evaluate the energy and its gradient.

!Arguments:
    integer,intent(in)     :: InitFunc,FinalFunc,MaxEnergyEval
    real(wp),intent(in) :: OverlapThreshold,MaxOverlapPenalty
    real(4),intent(in)     :: DataSaveMinTimeInterv,HessianSaveMinTimeInterv
    character(Glob_FileNameLength),intent(in) :: HessFileName
!Local variables:
    integer      i,j
    integer      np,npt,nfo,nfa,nfru,nv
    integer      OpenFileErr,ErrCode,NumOfEnergyEval,NumOfGradEval,NumOfFailures
    integer      NumOfEnergyEvalDuringFullOpt_Init
    logical      IsSwapFileOK,ExitNeeded
    logical      SaveHessian,IsHessFileOK,IsHessSaveSuccess
    real(wp)  Evalue,CurrentEnergy
    real(4)      TimeOfLastSave, TimeOfLastHessSave
    real(wp)  t
    integer      tas,InitFuncNew
    real(wp)  MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap
    real(wp),allocatable,dimension(:,:)   :: NonlinParamTemp
    integer,allocatable,dimension(:)         :: FuncNumTemp
    real(wp),allocatable,dimension(:)     :: TempR
    real(wp),allocatable,dimension(:)     :: x,grad
!Arrays used by DRMNG
    real(wp),allocatable,dimension(:)     :: D,V
    integer,parameter    :: LIV=60
    integer                 IV(LIV)
    integer                 LV
    integer                 ALG
    integer                 IVLMAT

    IsHessFileOK=.false.

    if (OverlapThreshold>=ONE) then
      Glob_OverlapPenaltyAllowed=.false.
    else
      Glob_OverlapPenaltyAllowed=.true.
      Glob_OverlapPenaltyThreshold2=OverlapThreshold*OverlapThreshold
      Glob_MaxOverlapPenalty=MaxOverlapPenalty
    endif
    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine FullOpt1I started'
      write(*,*) 'Simultaneous optimization of nonlinear parameters of basis functions'
      write(*,*) InitFunc,'  through',FinalFunc,'  will be attempted'
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Overlap threshold is ',abs(OverlapThreshold)
        write(*,*) 'Max value of a pair overlap penalty is ',Glob_MaxOverlapPenalty
        write(*,*) 'Warning! The energy value that will be shown during the optimization'
        write(*,*) 'may differ from the actual energy'
      else
        write(*,*) 'No constraints on overlaps will be imposed'
      endif
    endif

!Setting the values of some global variables
    Glob_GSEPSolutionMethod='I'
    Glob_nfa=Glob_CurrBasisSize
    Glob_nfo=FinalFunc-InitFunc+1
    Glob_nfru=Glob_CurrBasisSize-Glob_nfo
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    nfa=Glob_nfa
    nfo=Glob_nfo
    nfru=Glob_nfru
    nv=nfo*npt
    InitFuncNew=nfa+InitFunc-FinalFunc
    Glob_HSBuffLen=max(min(nfa*(nfa+1)/2,1000),30*nfa)

!Allocate some global arrays
    allocate(Glob_H(nfa,nfa))
    allocate(Glob_S(nfa,nfa))
    allocate(Glob_diagS(nfa))
    allocate(Glob_invD(nfa))
    allocate(Glob_D(2*npt,nfo,nfa))
    allocate(Glob_c(nfa))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))
    allocate(Glob_DkBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DkBuff2(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff1(2*npt,Glob_HSBuffLen))
    allocate(Glob_DlBuff2(2*npt,Glob_HSBuffLen))

!Allocate workspace for subroutine GSEPIIS, which is called
!inside EnergyIA and EnergyIB
    allocate(Glob_WorkForGSEPIIS(nfa))
    allocate(Glob_LastEigvector(nfa))
    Glob_LastEigvector(1:nfa)=ONE

!Allocate workspace for EnergyIB
    allocate(Glob_WkGR(nfo*npt))

!Allocate arrays used by DRMNG
    LV=71 + nv*(nv+13)/2 + 1
    if (Glob_ProcID==0) then
      allocate(D(nv))
      allocate(V(LV))
    endif

!Allocate workspace
    allocate(x(nv))
    allocate(grad(nv))

!Setting up a logical variable that determines
!whether the hessian should be saved from time to time
    if ((HessFileName==' ').or.(HessFileName=='none').or. &
        (HessFileName=='NONE').or.(HessFileName=='None')) then
      SaveHessian=.false.
    else
      SaveHessian=.true.
    endif

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

!Shifting basis functions InitFunc through FinalFunc to
!the very end and, if necessary, doing the permutation of
!matrix elements to reflect this change in function order.
    if (FinalFunc/=nfa) then
      !allocate space
      tas=FinalFunc-InitFunc+1
      allocate(NonlinParamTemp(1:npt,1:tas))
      allocate(FuncNumTemp(1:tas))
      allocate(TempR(1:tas))
      call PermuteFunctions(InitFunc,FinalFunc,FuncNumTemp,NonlinParamTemp)
      if (IsSwapFileOK) call PermuteMatrixElements(InitFunc,FinalFunc,TempR)
      deallocate(TempR)
      deallocate(FuncNumTemp)
      deallocate(NonlinParamTemp)
      !Allocate workspace for SaveResults (it uses Sort='yes'
      !option, which requires workspace)
      allocate(Glob_IntWorkArrForSaveResults(Glob_CurrBasisSize))
    endif

!Calculating the initial energy
    if (IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.false.,ErrCode)
    else
      if (Glob_ProcID==0) write(*,'(1x,a59)',advance='no') &
        'Computing matrix elements and solving eigenvalue problem...'
      Glob_CurrEnergy=EnergyIA(1,Glob_CurrBasisSize,.true.,ErrCode)
    endif
    if (ErrCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0165 in FullOpt1I: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    if (Glob_ProcID==0) write(*,*) ' done'
    if (Glob_ProcID==0) then
      call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Initial energy (without overlap penalty)  ', &
          Glob_CurrEnergy-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                           ',Glob_TotalOverlapPenalty
        write(*,*) 'Initial energy (including overlap penalty)',Glob_CurrEnergy
      else
        write(*,*) 'Initial energy                            ',Glob_CurrEnergy
      endif
      write(*,*)   'Maximal overlap                           ',MaxAbsOverlap
      write(*,*)   'Minimal overlap                           ',MinAbsOverlap
      write(*,*)   'Average abs value of overlap              ',AverageAbsOverlap
    endif

    call CPU_TIME(Glob_TimeSinceStart)
    TimeOfLastSave=Glob_TimeSinceStart
    TimeOfLastHessSave=Glob_TimeSinceStart

!Setting parameters for DRMNG

!Call DIVSET to get default values in IV and V arrays
!ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS
    ALG=2
    if (Glob_ProcID==0) then
      call DIVSET(ALG,IV,LIV,LV,V)
      IV(17)=1000000; IV(18)=1000000
      IV(19)=-1 !set summary print format
      IV(20)=0; IV(22)=0; IV(23)=-1; IV(24)=0
      V(31)=0.0_wp
      V(32)=2*epsilon(V(32))
      V(37)=2*epsilon(V(37))
      !V(35) GIVES THE MAXIMUM 2-NORM ALLOWED FOR D TIMES THE
      !VERY FIRST STEP THAT  DMNG ATTEMPTS.  THIS PARAMETER CAN
      !MARKEDLY AFFECT THE PERFORMANCE OF  DMNG.
      V(35)=Glob_MaxScStepAllowedInOpt
      IV(1)=12 !DIVSET has been called and some default values were changed
      nv=nfo*npt
    endif
    do i=1,nfo
      x((i-1)*npt+1:i*npt)=Glob_NonlinParam(1:npt,InitFuncNew+i-1)
    enddo

    if (Glob_ProcID==0) then
      !If SaveHessian=.true. then try to read the Hessian
      !from the file
      if (SaveHessian) then
        IVLMAT=IV(42)
        call ReadHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessFileOK)
        if (IsHessFileOK) IV(25)=0
      endif
      if ((.not.Glob_FullOptSaveD).or.(.not.IsHessFileOK).or.(.not.SaveHessian)) then
!    !Set the scaling vector
!    do i=1,nfo
!      t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
!      do j=1,np
!        !Make sure none of the D(i) will be zero or smaller than the threshold
!        D(npt*(i-1)+j)=1.0*ONE/max(abs(x(npt*(i-1)+j)),t)
!      enddo
!    enddo
        !Set the scaling vector
        t=max(ONE/(nfa*nfa*sqrt(ONE*nfa)),10000*epsilon(Glob_CurrEnergy))
        do i=1,nfo
          !t=maxval(abs(x(npt*(i-1)+1:npt*i-np)))/Glob_OptScalingThreshold
          do j=1,npt
            !Make sure none of the D(i) will be zero or smaller than the threshold
            !D(npt*(i-1)+j)=ONE/max(abs(x(npt*(i-1)+j)),t)
            D(npt*(i-1)+j)=t
            !write(*,*) 'i=',int(i,1),' j=',int(j,1),' D=',D(npt*(i-1)+j)
          enddo
        enddo
      endif
    endif

    ExitNeeded=.false.
    NumOfFailures=0
    NumOfEnergyEval=0
    NumOfGradEval=0
    NumOfEnergyEvalDuringFullOpt_Init=Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt

    do while (.not.(ExitNeeded))
      if (Glob_ProcID==0) call DRMNG(D, CurrentEnergy, grad, IV, LIV, LV, nv, V, x)
      call MPI_BCAST(IV,LIV,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      select case (IV(1))
      case (1) !Only energy is needed
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
        enddo
        Evalue=EnergyIA(InitFuncNew,nfa,.true.,ErrCode)
        NumOfEnergyEval=NumOfEnergyEval+1
        if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
          IV(2)=1
        else
          CurrentEnergy=Evalue
        endif
        if (Glob_ProcID==0) then
          if (Evalue<Glob_CurrEnergy) then
            !Save data and Hessian if necessary
            call CPU_TIME(Glob_TimeSinceStart)
            if (Glob_TimeSinceStart-TimeOfLastSave>DataSaveMinTimeInterv) then
              !Save the results if more than DataSaveMinTimeInterv seconds
              !have passed since the last save
              if (Glob_OverlapPenaltyAllowed) then
                Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
              else
                Glob_CurrEnergy=Evalue
              endif
              !Changing history
              Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
              Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
                NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval
              if (FinalFunc==nfa) then
                call SaveResults(Sort='no')
              else
                call SaveResults(Sort='yes')
              endif
              write(*,*) 'Data file has been updated'
              call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
              write(*,*) 'Some current statistics:'
              if (Glob_OverlapPenaltyAllowed) then
                write(*,*) 'Energy (without overlap penalty)  ',Evalue-Glob_TotalOverlapPenalty
                write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
                write(*,*) 'Energy (including overlap penalty)',Evalue
              else
                write(*,*) 'Energy                            ',Evalue
              endif
              write(*,*)   'Maximal overlap                   ',MaxAbsOverlap
              write(*,*)   'Minimal overlap                   ',MinAbsOverlap
              write(*,*)   'Average abs value of overlap      ',AverageAbsOverlap

              TimeOfLastSave=Glob_TimeSinceStart
            endif
            if ((Glob_TimeSinceStart-TimeOfLastHessSave>HessianSaveMinTimeInterv) &
                .and.(SaveHessian)) then
              !Save the Hessian if more than HessianSaveMinTimeInterv seconds
              !have passed since the last save
              call SaveHessianFile(V,IVLMAT,D,nv,HessFileName,IsHessSaveSuccess)
              TimeOfLastHessSave=Glob_TimeSinceStart
            endif
          endif
        endif
      case (2) !Only gradient is needed
        call MPI_BCAST(x,nv,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        do i=1,nfo
          Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
        enddo
        call EnergyIB(Evalue,grad,.true.,ErrCode)
        NumOfGradEval=NumOfGradEval+1
        if (ErrCode/=0) then
          NumOfFailures=NumOfFailures+1
          IV(2)=0
        endif
      case (3:8) !Some kind of convergence has been reached
        ExitNeeded=.true.
      case (9:10) !Function evaluation limit has been reached.
        !This never suppose to happen because we
        !count the number of function evaluations ourselves.
        ExitNeeded=.true.
      endselect
      if (NumOfFailures>Glob_MaxEnergyFailsAllowed) then
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,a,1x,a)') &
            'Error EC0166 in FullOpt1I: number of failures in energy or gradient', &
            'calculations during the optimization of nonlinear parameters', &
            'exceeded limit'
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
      if (NumOfEnergyEval>=MaxEnergyEval) then
        if (Glob_ProcID==0) then
          write(*,*) 'Warning WC0135 in FullOpt1I: number of energy evaluations reached limit'
          write(*,*) 'Optimization is terminated'
        endif
        ExitNeeded=.true.
      endif
    enddo

!Calculate the energy at the best point found
    do i=1,nfo
      Glob_NonlinParam(1:npt,InitFuncNew+i-1)=x((i-1)*npt+1:i*npt)
    enddo
    Evalue=EnergyIA(InitFuncNew,nfa,.true.,ErrCode)
    if (ErrCode/=0) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0167 in FullOpt1I: failed to evaluate energy after the optimization'
        write(*,*) 'of nonlinear parameters'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif
    if (Glob_OverlapPenaltyAllowed) then
      Glob_CurrEnergy=Evalue-Glob_TotalOverlapPenalty
    else
      Glob_CurrEnergy=Evalue
    endif

!Printing the number of energy/gradient evaluations
!and the energy after the optimization
    if (Glob_ProcID==0) then
      write(*,*)
      write(*,'(1x,a,1x,i0,a,i0)') 'Number of energy/gradient evaluations',NumOfEnergyEval,'/',NumOfGradEval
      write(*,*) 'Final energy and overlap statistics:'
      call GetOverlapStatistics(InitFuncNew,nfa,MaxAbsOverlap,MinAbsOverlap,AverageAbsOverlap)
      if (Glob_OverlapPenaltyAllowed) then
        write(*,*) 'Energy (without overlap penalty)  ',Evalue-Glob_TotalOverlapPenalty
        write(*,*) 'Overlap penalty                   ',Glob_TotalOverlapPenalty
        write(*,*) 'Energy (including overlap penalty)',Evalue
      else
        write(*,*) 'Energy                            ',Evalue
      endif
      write(*,*)   'Maximal overlap                   ',MaxAbsOverlap
      write(*,*)   'Minimal overlap                   ',MinAbsOverlap
      write(*,*)   'Average abs value of overlap      ',AverageAbsOverlap
    endif

!Adding data to history
    Glob_History(Glob_CurrBasisSize)%Energy=Glob_CurrEnergy
    Glob_History(Glob_CurrBasisSize)%NumOfEnergyEvalDuringFullOpt= &
      NumOfEnergyEvalDuringFullOpt_Init+NumOfEnergyEval

!Shifting basis functions InitFunc through FinalFunc back from
!the very end to the middle, where they initially were, and, if
!necessary, doing the permutation of matrix elements to reflect
!this change in function order.
    if (FinalFunc/=nfa) then
      !allocate space
      tas=nfa-FinalFunc
      allocate(NonlinParamTemp(1:npt,1:tas))
      allocate(FuncNumTemp(1:tas))
      allocate(TempR(1:tas))
      call PermuteFunctions(InitFunc,InitFunc+tas-1,FuncNumTemp,NonlinParamTemp)
      if (Glob_UseSwapFile) call PermuteMatrixElements(InitFunc,InitFunc+tas-1,TempR)
      deallocate(TempR)
      deallocate(FuncNumTemp)
      deallocate(NonlinParamTemp)
      !deallocate workspace for SaveResults
      deallocate(Glob_IntWorkArrForSaveResults)
    endif

!saving results
    if (Glob_ProcID==0) call SaveResults(Sort='no')

    call StoreMatricesInSwapFile()

    if (Glob_OverlapPenaltyAllowed) Glob_OverlapPenaltyAllowed=.false.

!deallocate workspace
    deallocate(grad)
    deallocate(x)

!deallocate arrays used by DRMNG
    if (Glob_ProcID==0) then
      deallocate(V)
      deallocate(D)
    endif
!deallocate workspace for EnergyIB
    deallocate(Glob_WkGR)

!Deallocate workspace for GSEPSII
    deallocate(Glob_LastEigvector)
    deallocate(Glob_WorkForGSEPIIS)

!Deallocate some global arrays
    deallocate(Glob_DlBuff2)
    deallocate(Glob_DlBuff1)
    deallocate(Glob_DkBuff2)
    deallocate(Glob_DkBuff1)
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_D)
    deallocate(Glob_invD)
    deallocate(Glob_diagS)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) write (*,*) 'Routine FullOpt1I finished'

  end subroutine FullOpt1I

  subroutine SolveEliminationGSEP(GSEPSolMethod,MatrixOrder,Evalue,ErrorCode)
!Subroutine SolveEliminationGSEP provides the common eigensolver boundary used
!by the elimination and separation BBOP routines. It constructs the initial
!factorization. After the basis change, Q cleanup drivers preserve this state
!through delete_symmetric or replace_symmetric; only G rebuilds and calls this
!routine a second time.
!
!The G path preserves the historical DSYGVX layout: it materializes the upper
!triangle and restores Hamiltonian and overlap diagonals before calling LAPACK.
!The Q path must not do that. StoreHS has already placed normalized, unshifted H
!and S in their canonical lower triangles, including their physical diagonals.
!qrlinalg receives that representation directly and owns only its factors.
!
!Arguments:
    character(1),intent(in) :: GSEPSolMethod
    integer,intent(in)      :: MatrixOrder
    real(wp),intent(out)    :: Evalue
    integer,intent(out)     :: ErrorCode
!Local variables:
    integer i,j,IFAIL(1),NumOfEigvalsFound
    real(wp) EVs(1)

    Evalue=huge(Evalue)
    ErrorCode=Q_METHOD_INVALID_ARGUMENT

    select case (GSEPSolMethod)
    case('G')
      do i=1,MatrixOrder
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,MatrixOrder
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo

      if (Glob_ProcID==0) then
        call DSYGVX(1,'V','I','U',MatrixOrder,Glob_H,Glob_HSLeadDim, &
          Glob_S,Glob_HSLeadDim,ZERO,ZERO,Glob_WhichEigenvalue, &
          Glob_WhichEigenvalue,Glob_AbsTolForDSYGVX,NumOfEigvalsFound, &
          EVs,Glob_c,MatrixOrder,Glob_WorkForDSYGVX, &
          Glob_LWorkForDSYGVX,Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
        Evalue=EVs(1)
      endif
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_c,MatrixOrder,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    case('Q')
      !There is no useful previous vector after a structural elimination or a
      !random separation. ONE gives inverse iteration a deterministic nonzero
      !starting vector and SolveQ returns an S-normalized coefficient vector.
      Glob_c(1:MatrixOrder)=ONE
      call PrepareQWorkspace(MatrixOrder,MatrixOrder,1,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) then
        Q_Workspace%MatricesAreCanonical=.true.
        call FactorizeQFresh(ErrorCode)
      endif
      if (ErrorCode==Q_METHOD_SUCCESS) call SolveQ(Evalue,ErrorCode)
    endselect

  end subroutine SolveEliminationGSEP

  subroutine DeleteQMaskedFunctions(RemoveMask,ErrorCode)
!Subroutine DeleteQMaskedFunctions removes every marked canonical basis
!function from both the qrlinalg state and the physical matrix representation.
!Deletions are submitted in descending canonical order, so an original index
!continues to identify the same row and column after every preceding deletion.
!For r removed functions this costs O(r*n**2), while the historical cleanup
!path recalculated O(n**2) matrix elements and constructed another O(n**3)
!factorization even though every survivor-survivor element was unchanged.
!
!The physical lower triangles are compacted only after all root-owned QR
!operations succeed. Their in-place ascending survivor copy is safe: each
!source row and column has an original index not smaller than its destination,
!and no write can destroy a matrix element needed by a later survivor. The
!upper triangles remain deliberately unspecified under the Q canonical-layout
!contract. Raw overlap norms and the current eigenvector are compacted by the
!same survivor map.
!
!A public qrlinalg deletion cannot fail after the complete metadata preflight,
!but the recovery path is still explicit. If a future library implementation
!introduces a failure, the untouched physical matrices reconstruct the old
!factorization before this routine reports the original error.
!
!Arguments:
    integer,intent(in)  :: RemoveMask(:)
    integer,intent(out) :: ErrorCode
!Local variables:
    integer i,j,OldI,OldJ,OldOrder,NewOrder,RootError,RecoveryError
    integer,allocatable :: Survivor(:)

    ErrorCode=Q_METHOD_INVALID_ARGUMENT
    OldOrder=Q_Workspace%MatrixOrder
    if (Glob_GSEPSolutionMethod/='Q') return
    if (.not.Q_Workspace%MatricesAreCanonical) return
    if (.not.Q_Workspace%FactorsMatchMatrices) return
    if (OldOrder<2) return
    if (size(RemoveMask)/=OldOrder) return
    NewOrder=count(RemoveMask==0)
    if ((NewOrder<1).or.(NewOrder>=OldOrder)) return
    if (.not.allocated(Glob_H)) return
    if (.not.allocated(Glob_S)) return
    if (.not.allocated(Glob_diagS)) return
    if (.not.allocated(Glob_c)) return
    if ((size(Glob_H,1)<OldOrder).or.(size(Glob_H,2)<OldOrder)) return
    if ((size(Glob_S,1)<OldOrder).or.(size(Glob_S,2)<OldOrder)) return
    if (size(Glob_diagS)<OldOrder) return
    if (size(Glob_c)<OldOrder) return

    allocate(Survivor(NewOrder),stat=RootError)
    if (RootError/=0) then
      ErrorCode=Q_METHOD_ALLOCATION_ERROR
      return
    endif
    j=0
    do i=1,OldOrder
      if (RemoveMask(i)==0) then
        j=j+1
        Survivor(j)=i
      endif
    enddo

    RootError=Q_METHOD_SUCCESS
    if (Glob_ProcID==0) then
      if (.not.Q_Workspace%Factors%is_valid()) RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%order()/=OldOrder) RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%get_capacity()/=Q_Workspace%Capacity) &
        RootError=Q_METHOD_INVALID_STATE
      if (Q_Workspace%Factors%get_shift()/=Glob_ApproxEnergy) &
        RootError=Q_METHOD_INVALID_STATE
      if (RootError==Q_METHOD_SUCCESS) then
        do i=OldOrder,1,-1
          if (RemoveMask(i)/=0) then
            call Q_Workspace%Factors%delete_symmetric(i,RootError)
            if (RootError/=Q_METHOD_SUCCESS) exit
          endif
        enddo
      endif
    endif
    call MPI_BCAST(RootError,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    if (RootError/=Q_METHOD_SUCCESS) then
      Q_Workspace%FactorsMatchMatrices=.false.
      call FactorizeQFresh(RecoveryError)
      if (RecoveryError/=Q_METHOD_SUCCESS) then
        ErrorCode=RecoveryError
      else
        ErrorCode=RootError
      endif
      deallocate(Survivor)
      return
    endif

    do i=1,NewOrder
      OldI=Survivor(i)
      Glob_c(i)=Glob_c(OldI)
      Glob_diagS(i)=Glob_diagS(OldI)
      do j=1,i
        OldJ=Survivor(j)
        Glob_H(i,j)=Glob_H(OldI,OldJ)
        Glob_S(i,j)=Glob_S(OldI,OldJ)
      enddo
    enddo

    Q_Workspace%MatrixOrder=NewOrder
    Q_Workspace%NumActive=0
    Q_Workspace%ActiveFunction=0
    Q_Workspace%ActivePosition=0
    Q_Workspace%MatrixParametersAreStored=.false.
    Q_Workspace%TrialIsReady=.false.
    Q_Workspace%TrialHasDerivatives=.false.
    Q_Workspace%AcceptedPointIsStored=.false.
    Q_Workspace%FactorsMatchMatrices=.true.
    Q_Workspace%LastEigenpairResidual=huge(ONE)
    Q_Workspace%LastFactorResidual=huge(ONE)
    ErrorCode=Q_METHOD_SUCCESS
    deallocate(Survivor)

  end subroutine DeleteQMaskedFunctions

  subroutine EliminateLittleContribFunc(LinCoeffThreshold,FileName,PrintInfoSpec,GSEPSolMethod)
!Subroutine EliminateLittleContribFunc eliminates basis
!functions whose contribution to the energy is small. More
!precisely, it eliminates functions that have linear coefficients
!whose absolute values are smaller than LinCoeffThreshold. It is
!important to note that this subroutine uses coefficients in
!front of normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1).
!The result is stored in file whose name is defined by parameter
!FileName. After saving the results the subroutine terminates
!the program.
!Parameter PrintInfoSpec specifies what information should be
!shown regarding the linear coefficients:
! PrintInfoSpec=0,1 : the subroutine does not print any specific info
!                   regarding linear coefficients.
! PrintInfoSpec=2   : the subroutine prints linear coefficients
!                     of all functions
! GSEPSolMethod     : optional solver selection. It defaults to G so existing
!                     callers retain their behavior; main passes Q explicitly.

!Arguments:
    real(wp),intent(in)                    :: LinCoeffThreshold
    character(Glob_FileNameLength),intent(in) :: FileName
    integer,intent(in)                        :: PrintInfoSpec
    character(1),intent(in),optional          :: GSEPSolMethod

!Local variables:
    integer      i,j
    integer      np,npt,cbs
    integer      OpenFileErr,ErrorCode
    logical      IsSwapFileOK
    integer      BlockSizeForDSYGVX
    real(wp)  Evalue
    real(wp)  Min_c,Max_c
    real(wp)  Aver_c
    real(wp),allocatable,dimension(:,:)   :: NonlinParamTemp
    integer,allocatable,dimension(:)      :: MaskArray
    integer,allocatable,dimension(:,:)    :: IndTemp
    character(Glob_FileNameLength)           :: ch_temp
    character(1)                              :: Method

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine EliminateLittleContribFunc started'
      write(*,*) 'Linear coefficient threshold value is',LinCoeffThreshold
    endif

!Setting the values of some global variables
    Method='G'
    if (present(GSEPSolMethod)) Method=GSEPSolMethod
    Glob_GSEPSolutionMethod=Method
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (Method=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (Method=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Reading data from swap file
    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
    else
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) ' done'
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
    endif

    call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0170 in EliminateLittleContribFunc: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*) 'Basis size before elimination',cbs
      write(*,*) 'Energy before elimination    ',Evalue
    endif

    allocate(NonlinParamTemp(1:npt,cbs))
    allocate(MaskArray(1:cbs))
    allocate(IndTemp(cbs,2))
    MaskArray=0

    if ((PrintInfoSpec>1).and.(Glob_ProcID==0)) then
      write(*,*) 'List of all linear coefficients:'
      do i=1,cbs
        write(*,'(i6,a4,f19.12)') i,'  c=',Glob_c(i)
      enddo
      write(*,*)
    endif

    Min_c=Glob_c(1)
    Max_c=Glob_c(1)
    Aver_c=ZERO
    j=0
    do i=1,cbs
      if (abs(Glob_c(i))>=LinCoeffThreshold) then
        NonlinParamTemp(1:npt,i-j)=Glob_NonlinParam(1:npt,i)
        IndTemp(i-j,1:2)=Glob_Index(i,1:2)
        IndTemp(i-j,1)=IndTemp(i,1)
        IndTemp(i-j,2)=IndTemp(i,2)
      else
        if ((j==0).and.(Glob_ProcID==0)) write(*,*) 'Little contributing function list:'
        j=j+1
        MaskArray(i)=1
        if (Glob_ProcID==0) &
          write(*,'(i6,a1,i6,a4,f19.12)') j,':',i,'  c=',Glob_c(i)
      endif
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Aver_c=Aver_c+abs(Glob_c(i))
    enddo
    Aver_c=Aver_c/cbs

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,'(1x,a41,e16.9)') 'Smallest by magnitude linear coeff.     =',Min_c
      write(*,'(1x,a41,e16.9)') 'Largest by magnitude linear coeff.      =',Max_c
      write(*,'(1x,a41,e16.9)') 'Average absolute value of linear coeff. =',Aver_c
      write(*,*)
    endif

    if (j==0) then
      if (Glob_ProcID==0) then
        write(*,*) 'There are no functions with the contribution'
        write(*,*) 'smaller than ',LinCoeffThreshold
        write(*,*) 'No output file have been written. Program will now stop'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    !An empty basis has no generalized eigenproblem and qrlinalg intentionally
    !has no valid order-zero state. Fail before changing the published basis
    !size or overwriting the input/output data with an unusable result.
    if (j==cbs) then
      if (Glob_ProcID==0) then
        write(*,*) 'All basis functions are below the coefficient threshold'
        write(*,*) 'No output file has been written. Program will now stop'
      endif
      call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
    endif

    if (Glob_ProcID==0) then
      write(*,*) 'Basis size before elimination',cbs
      write(*,*) 'Energy before elimination    ',Evalue
    endif

    Glob_NonlinParam(1:npt,1:cbs-j)=NonlinParamTemp(1:npt,1:cbs-j)
    Glob_Index(1:cbs-j,1:2)=IndTemp(1:cbs-j,1:2)

    if (Glob_ProcID==0) then
      if (Method=='Q') then
        write(*,*) 'Deleting selected rows and columns from the QR factorization...'
      else
        write(*,*) 'Computing matrix elements and solving eigenvalue problem with the'
        write(*,*) 'basis where little contributing functions are eliminated...'
      endif
    endif
    if (Method=='Q') then
      !Every survivor-survivor matrix element is unchanged. Preserve its
      !canonical value and delete the matching factor rows and columns instead
      !of evaluating the complete smaller matrix for a second time.
      call DeleteQMaskedFunctions(MaskArray,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) then
        Glob_CurrBasisSize=cbs-j
        cbs=Glob_CurrBasisSize
        call SolveQ(Evalue,ErrorCode)
      endif
    else
      Glob_CurrBasisSize=cbs-j
      cbs=Glob_CurrBasisSize
      call ComputeMatElem(1,cbs)
      call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    endif
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0171 in EliminateLittleContribFunc: energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) then
      write(*,*) 'Basis size after elimination ',cbs
      write(*,*) 'Energy after elimination     ',Evalue
    endif

    do i=1,cbs
      Glob_History(i)%Energy=ZERO
      Glob_History(i)%CyclesDone=0
      Glob_History(i)%InitFuncAtLastStep=0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
    enddo
    Glob_History(cbs)%Energy=Evalue

    Glob_CurrEnergy=Evalue
    Glob_LastEigvalTol=1.0E+35_wp
    Glob_BestEigvalTol=1.0E+35_wp
    Glob_WorstEigvalTol=1.0E-35_wp

    ch_temp=Glob_DataFileName
    Glob_DataFileName=FileName
    if (Glob_ProcID==0) call SaveResults(Sort='no')
    Glob_DataFileName=ch_temp

    deallocate(MaskArray)
    deallocate(IndTemp)
    deallocate(NonlinParamTemp)

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    if (Method=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

!Deallocate workspace for DSYGVX
    if (Method=='G') then
      deallocate(Glob_IWorkForDSYGVX)
      deallocate(Glob_WorkForDSYGVX)
    else
      call ClearQWorkspace()
    endif

    if (Glob_ProcID==0) then
      i=len_trim(FileName)
      write(*,*) 'New basis has been saved in file',FileName(1:i)
      write(*,*) 'Program will now stop'
    endif

    call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

  end subroutine EliminateLittleContribFunc

  subroutine EliminateLinDepFunc(LinDepThreshold,FileName,PrintInfoSpec,GSEPSolMethod)
!Subroutine EliminateLinDepFunc eliminates linearly dependent
!functions. It checks for pair linear dependency only. It
!removes those functions from the basis whose overlap (absolute value)
!with any of other basis function with a smaller number is greater
!than LinDepThreshold. It is important to note that this subroutine
!uses normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1).
!The result is stored in a file whose name is defined by parameter
!FileName. After saving the results the program is terminated.
!Parameter PrintInfoSpec specifies what information should be
!printed during linear dependency check:
! PrintInfoSpec=0 : the subroutine does not print any info
!                   regarding basis functions that are linearly dependent.
! PrintInfoSpec=1 : the subroutine prints the overlap values of
!                   linearly dependent functions.
! PrintInfoSpec=2 : same as the previous case, but in addition it also prints the
!                   nonlinear parameters of linearly dependent functions.
! GSEPSolMethod   : optional solver selection. It defaults to G so existing
!                   callers retain their behavior; main passes Q explicitly.

!Arguments:
    real(wp),intent(in)                    :: LinDepThreshold
    character(Glob_FileNameLength),intent(in) :: FileName
    integer,intent(in)                        :: PrintInfoSpec
    character(1),intent(in),optional          :: GSEPSolMethod

!Local variables:
    integer        i,j,k
    integer        np,npt,cbs
    integer        OpenFileErr,ErrorCode
    logical        IsSwapFileOK
    integer        BlockSizeForDSYGVX
    real(wp)    Evalue
    real(wp)    MaxOverlap,MinOverlap
    real(wp)    AverOverlap
    real(wp)    Min_c,Max_c
    real(wp)    Average_c
    integer,allocatable,dimension(:)    :: MaskArray
    character(Glob_FileNameLength)      :: ch_temp
    character(1)                        :: Method

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine EliminateLinDepFunc started'
      write(*,*) 'Overlap threshold value is',LinDepThreshold
    endif

!Setting the values of some global variables
    Method='G'
    if (present(GSEPSolMethod)) Method=GSEPSolMethod
    Glob_GSEPSolutionMethod=Method
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (Method=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (Method=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Allocate local workspace
    allocate(MaskArray(1:cbs))

!Reading data from swap file

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) ' done'
    endif

    if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
    call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0175 in EliminateLinDepFunc: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*) 'Pair linear dependency check:'
      write(*,*)
    endif
!Check overlap
    MaskArray(1:cbs)=0
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>LinDepThreshold) then
          MaskArray(j)=MaskArray(j)+1
          k=k+1
          if (Glob_ProcID==0) then
            if (PrintInfoSpec>0) write(*,'(i6,a1,i6,i6,a5,f17.14)') k,':',i,j,'   S=',Glob_S(j,i)
            if (PrintInfoSpec==2) then
              write(*,'(1x,i6,1x,i6,1x,i6)',advance='no') i,Glob_Index(i,1),Glob_Index(i,2)
              call writerealarradv(6,Glob_NonlinParam(1:Glob_npt,i),Glob_npt)
              write(*,*) '      c=',Glob_c(i)
              write(*,'(1x,i6,1x,i6,1x,i6)',advance='no') j,Glob_Index(j,1),Glob_Index(j,2)
              call writerealarradv(6,Glob_NonlinParam(1:Glob_npt,j),Glob_npt)
              write(*,*) '      c=',Glob_c(j)
            endif
            write(*,*)
          endif
        endif
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      !An order-one basis has no off-diagonal pair statistics.
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients:
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) '========== Summary before elimination: =========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size before elimination =',cbs
      write(*,*) 'Energy before elimination     =',Evalue
      write(*,*) '================================================'
      write(*,*)
    endif

    if (k==0) then
      if (Glob_ProcID==0) then
        write(*,*) 'No linearly dependent functions have been found'
        write(*,*) 'No file have been written. Program will now stop'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    j=0
    i=1
    do while (i+j<=cbs)
      if (MaskArray(i+j)>0) then
        j=j+1
      else
        Glob_NonlinParam(1:npt,i)=Glob_NonlinParam(1:npt,i+j)
        Glob_Index(i,1)=Glob_Index(i+j,1)
        Glob_Index(i,2)=Glob_Index(i+j,2)
        i=i+1
      endif
    enddo

    !k counts offending pairs for the diagnostic list. One later function can
    !overlap several earlier functions, but MaskArray removes that function
    !only once. j is the number of unique masked functions actually skipped by
    !the compaction loop and therefore defines the new basis order.
    if (Method=='Q') then
      if (Glob_ProcID==0) write(*,'(1x,a)',advance='no') &
        'Deleting selected rows and columns from QR factors...'
      call DeleteQMaskedFunctions(MaskArray,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) then
        Glob_CurrBasisSize=cbs-j
        cbs=Glob_CurrBasisSize
        call SolveQ(Evalue,ErrorCode)
      endif
    else
      Glob_CurrBasisSize=cbs-j
      cbs=Glob_CurrBasisSize
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) then
        write(*,*) ' done'
        write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      endif
      call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    endif
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0176 in EliminateLinDepFunc: energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!Check overlap
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      !Elimination can legitimately leave one surviving basis function.
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients:
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*)
      write(*,*) '========== Summary after elimination: =========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size after elimination  =',cbs
      write(*,*) 'Energy after elimination      =',Evalue
      write(*,*) '================================================'
      write(*,*)
    endif

    do i=1,cbs
      Glob_History(i)%Energy=ZERO
      Glob_History(i)%CyclesDone=0
      Glob_History(i)%InitFuncAtLastStep=0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
    enddo
    Glob_History(cbs)%Energy=Evalue

    Glob_CurrEnergy=Evalue
    Glob_LastEigvalTol=1.0E+35_wp
    Glob_BestEigvalTol=1.0E+35_wp
    Glob_WorstEigvalTol=1.0E-35_wp

    ch_temp=Glob_DataFileName
    Glob_DataFileName=FileName
    if (Glob_ProcID==0) call SaveResults(Sort='no')
    Glob_DataFileName=ch_temp

!deallocate local workspace
    deallocate(MaskArray)

!Deallocate workspace for DSYGVX
    if (Method=='G') then
      deallocate(Glob_IWorkForDSYGVX)
      deallocate(Glob_WorkForDSYGVX)
    else
      call ClearQWorkspace()
    endif

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    if (Method=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      i=len_trim(FileName)
      write(*,*) 'New basis has been saved in file',FileName(1:i)
      write(*,*) 'Program will now stop'
    endif

    call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

  end subroutine EliminateLinDepFunc

  subroutine SeparateLinDepFunc(LinDepThreshold,SeparationParam,FileName,PrintInfoSpec,GSEPSolMethod)
!Subroutine SeparateLinDepFunc does exactly the same thing as
!subroutine EliminateLinDepFunc does, but without throwing away
!linearly dependent functions. Instead, it changes the parameters
!of such functions randomly (the random shift is controlled by
!argument SeparationParam, so that a_new lies within
!interval [a_old*(1-SeparationParam),a_old*(1+SeparationParam)]).
!Optional GSEPSolMethod defaults to G for source compatibility. The Q path
!keeps the basis canonical and consumes only the lower H/S triangles.

!Arguments:
    real(wp),intent(in)                    :: LinDepThreshold
    real(wp),intent(in)                    :: SeparationParam
    character(Glob_FileNameLength),intent(in) :: FileName
    integer,intent(in)                        :: PrintInfoSpec
    character(1),intent(in),optional          :: GSEPSolMethod

!Local variables:
    integer        i,j,k,NumActive
    integer        np,npt,cbs
    integer        OpenFileErr,ErrorCode
    logical        IsSwapFileOK
    integer        BlockSizeForDSYGVX
    real(wp)    Evalue, r
    real(8)        r8
    real(wp)    MaxOverlap,MinOverlap
    real(wp)    AverOverlap
    real(wp)    Min_c,Max_c
    real(wp)    Average_c
    integer,allocatable,dimension(:)    :: MaskArray,ActiveFunction
    character(Glob_FileNameLength)      :: ch_temp
    character(1)                        :: Method

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine SeparateLinDepFunc started'
      write(*,*) 'Overlap threshold value is',LinDepThreshold
      write(*,*) 'Separation paramameter is ',SeparationParam
    endif

!Setting the values of some global variables
    Method='G'
    if (present(GSEPSolMethod)) Method=GSEPSolMethod
    Glob_GSEPSolutionMethod=Method
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (Method=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (Method=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Allocate local workspace
    allocate(MaskArray(1:cbs))

!Reading data from swap file

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) ' done'
    endif

    if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
    call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0180 in SeparateLinDepFunc: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*) 'Pair linear dependency check:'
      write(*,*)
    endif
!Check overlap
    MaskArray(1:cbs)=0
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>LinDepThreshold) then
          MaskArray(j)=MaskArray(j)+1
          k=k+1
          if (Glob_ProcID==0) then
            if (PrintInfoSpec>0) write(*,'(i6,a1,i6,i6,a5,f17.14)') k,':',i,j,'   S=',Glob_S(j,i)
            if (PrintInfoSpec==2) then
              write(*,'(1x,i6,1x,i6,1x,i6)',advance='no') i,Glob_Index(i,1),Glob_Index(i,2)
              call writerealarradv(6,Glob_NonlinParam(1:Glob_npt,i),Glob_npt)
              write(*,*) '      c=',Glob_c(i)
              write(*,'(1x,i6,1x,i6,1x,i6)',advance='no') j,Glob_Index(j,1),Glob_Index(j,2)
              call writerealarradv(6,Glob_NonlinParam(1:Glob_npt,j),Glob_npt)
              write(*,*) '      c=',Glob_c(j)
            endif
            write(*,*)
          endif
        endif
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients:
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) '========== Summary before separation: =========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size before separation =',cbs
      write(*,*) 'Energy before separation     =',Evalue
      write(*,*) '================================================'
      write(*,*)
    endif

    if (k==0) then
      if (Glob_ProcID==0) then
        write(*,*) 'No linearly dependent functions have been found'
        write(*,*) 'No file have been written. Program will now stop'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Method=='Q') then
      !The initial cleanup workspace holds one transaction column. Expand that
      !storage only to the number of unique functions selected by the overlap
      !mask, while retaining the already computed full QR factorization.
      NumActive=count(MaskArray>0)
      call EnsureQActiveCapacity(NumActive,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) then
        allocate(ActiveFunction(NumActive))
        j=0
        do i=1,cbs
          if (MaskArray(i)>0) then
            j=j+1
            ActiveFunction(j)=i
          endif
        enddo
        call SetQActiveFunctions(ActiveFunction,ErrorCode)
      endif
      if (ErrorCode==Q_METHOD_SUCCESS) call CaptureQMatrixParameters(ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) then
        if (Glob_ProcID==0) write(*,*) &
          'Error EC0181 in SeparateLinDepFunc: Q transaction cannot be prepared'
        call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
      endif
    endif

    j=0
    i=1
    do i=1,cbs
      if (MaskArray(i)>0) then
        do j=1,npt
          call random_number(r8)
          r=TWO*(r8-ONEHALF)*SeparationParam
          Glob_NonlinParam(j,i)=Glob_NonlinParam(j,i)*(1+r)
        enddo
      endif
    enddo
    call MPI_BCAST(Glob_NonlinParam,cbs*npt,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Method=='Q') then
      if (Glob_ProcID==0) write(*,'(1x,a)',advance='no') &
        'Computing selected matrix columns and updating QR...'
      call AssembleQTrial(.false.,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) call ApplyQTrial(ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) call SolveQ(Evalue,ErrorCode)
    else
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) then
        write(*,*) ' done'
        write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      endif
      call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    endif
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0181 in EliminateLinDepFunc: energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!Check overlap
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients:
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*)
      write(*,*) '========== Summary after separation: ==========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size after separation  =',cbs
      write(*,*) 'Energy after separation      =',Evalue
      write(*,*) '================================================'
      write(*,*)
    endif

    do i=1,cbs
      Glob_History(i)%Energy=ZERO
      Glob_History(i)%CyclesDone=0
      Glob_History(i)%InitFuncAtLastStep=0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
    enddo
    Glob_History(cbs)%Energy=Evalue

    Glob_CurrEnergy=Evalue
    Glob_LastEigvalTol=1.0E+35_wp
    Glob_BestEigvalTol=1.0E+35_wp
    Glob_WorstEigvalTol=1.0E-35_wp

    ch_temp=Glob_DataFileName
    Glob_DataFileName=FileName
    if (Glob_ProcID==0) call SaveResults(Sort='no')
    Glob_DataFileName=ch_temp

!deallocate local workspace
    if (allocated(ActiveFunction)) deallocate(ActiveFunction)
    deallocate(MaskArray)

!Deallocate workspace for DSYGVX
    if (Method=='G') then
      deallocate(Glob_IWorkForDSYGVX)
      deallocate(Glob_WorkForDSYGVX)
    else
      call ClearQWorkspace()
    endif

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    if (Method=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      i=len_trim(FileName)
      write(*,*) 'New basis has been saved in file',FileName(1:i)
      write(*,*) 'Program will now stop'
    endif

    call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

  end subroutine SeparateLinDepFunc

  subroutine SeparateFuncLargeCoeff(LCThreshold,SeparationParam,FileName,PrintInfoSpec,GSEPSolMethod)
!Subroutine SeparateFuncLargeCoeff changes linear parameters of
!those basis functions whose linear coefficients (more precisely their
!absolute values) exceed LCThreshold. It is important to note that this
!subroutine uses normalized functions (so that Glob_S is the overlap
!matrix of normalized functions, with Glob_S(i,i)=1). The change
!in nonlinear parameters is controlled by argument SeparationParam.
!Basically the nonlinear parameters of bad functions are chosen
!randomly from the interval [a_old*(1-SeparationParam), a_old*(1+SeparationParam)].
!After saving the results the program is terminated.
!Parameter PrintInfoSpec specifies what information should be
!shown:
! PrintInfoSpec=0 : the subroutine does not print any info
!                   about basis functions.
! PrintInfoSpec=1 : the subroutine prints the values of the linear and
!                   nonlinear parameters of bad functions.
! PrintInfoSpec=2 : same as the previous case, but in addition it also prints
!                   the linear parameters of all basis functions.
! GSEPSolMethod   : optional solver selection. It defaults to G so existing
!                   callers retain their behavior; main passes Q explicitly.

!Arguments:
    real(wp),intent(in)                    :: LCThreshold
    real(wp),intent(in)                    :: SeparationParam
    character(Glob_FileNameLength),intent(in) :: FileName
    integer,intent(in)                        :: PrintInfoSpec
    character(1),intent(in),optional          :: GSEPSolMethod

!Local variables:
    integer        i,j,k,NumActive
    integer        np,npt,cbs
    integer        OpenFileErr,ErrorCode
    logical        IsSwapFileOK
    integer        BlockSizeForDSYGVX
    real(wp)    Evalue, r
    real(8)        r8
    real(wp)    MaxOverlap,MinOverlap
    real(wp)    AverOverlap
    real(wp)    Min_c,Max_c
    real(wp)    Average_c
    integer,allocatable,dimension(:) :: ActiveFunction
    character(Glob_FileNameLength)  :: ch_temp
    character(1)                    :: Method

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine SeparateFuncLargeCoeff started'
      write(*,*) 'Linear coefficient threshold value is',LCThreshold
      write(*,*) 'Separation parameter is              ',SeparationParam
    endif

!Setting the values of some global variables
    Method='G'
    if (present(GSEPSolMethod)) Method=GSEPSolMethod
    Glob_GSEPSolutionMethod=Method
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize

!Allocate some global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (Method=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (Method=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Reading data from swap file
    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) ' done'
    endif

    if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
    call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) &
        'Error EC0185 in SeparateFuncLargeCoeff: initial energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*)
      if (PrintInfoSpec>=2) then
        write(*,*) 'Linear coefficients of all basis functions before separation:'
        do i=1,cbs
          write(*,'(1x,i6,a4,f19.12)')  i,'  c=',Glob_c(i)
        enddo
      endif
      write(*,*)
    endif

!Check overlap
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients
    k=0
    do i=1,cbs
      if (abs(Glob_c(i))>LCThreshold) then
        k=k+1
        if ((Glob_ProcID==0).and.(PrintInfoSpec>=1)) then
          if (k==1) write(*,*) 'Functions whose linear coefficients exceed threshold'
          write(*,'(1x,i5,a3,i6,a4,f19.12)') k,':  ',i,'  c=',Glob_c(i)
          write(*,*) Glob_NonlinParam(1:Glob_npt,i)
        endif
      endif
    enddo

    if (k==0) then
      if (Glob_ProcID==0) then
        write(*,*) 'There are no functions whose linear coefficients exceed threshold'
        write(*,*) 'No file have been written. Program will now stop'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

    if (Method=='Q') then
      !The coefficient scan already provides the exact active count. Reserve
      !only those transaction columns and preserve the initial factorization.
      NumActive=k
      call EnsureQActiveCapacity(NumActive,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) then
        allocate(ActiveFunction(NumActive))
        k=0
        do i=1,cbs
          if (abs(Glob_c(i))>LCThreshold) then
            k=k+1
            ActiveFunction(k)=i
          endif
        enddo
        call SetQActiveFunctions(ActiveFunction,ErrorCode)
      endif
      if (ErrorCode==Q_METHOD_SUCCESS) call CaptureQMatrixParameters(ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) then
        if (Glob_ProcID==0) write(*,*) &
          'Error EC0186 in SeparateFuncLargeCoeff: Q transaction cannot be prepared'
        call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
      endif
    endif
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*) '========== Summary before separation: =========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size before separation =',cbs
      write(*,*) 'Energy before separation     =',Evalue
      write(*,*) '================================================'
    endif

!now change the nonlinear parameters of bad functions
    do i=1,cbs
      if (abs(Glob_c(i))>LCThreshold) then
        do j=1,npt
          call random_number(r8)
          r=(r8-ONEHALF)*2*SeparationParam
          Glob_NonlinParam(j,i)=Glob_NonlinParam(j,i)*(1+r)
        enddo
      endif
    enddo
    call MPI_BCAST(Glob_NonlinParam,cbs*Glob_npt,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    if (Method=='Q') then
      if (Glob_ProcID==0) write(*,'(1x,a)',advance='no') &
        'Computing selected matrix columns and updating QR...'
      call AssembleQTrial(.false.,ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) call ApplyQTrial(ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) call SolveQ(Evalue,ErrorCode)
    else
      if (Glob_ProcID==0) write(*,'(1x,a28)',advance='no') 'Computing matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) then
        write(*,*) ' done'
        write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      endif
      call SolveEliminationGSEP(Method,cbs,Evalue,ErrorCode)
    endif
    if (ErrorCode/=0) then
      if (Glob_ProcID==0) write(*,*) 'Error EC0186 in SeparateFuncLargeCoeff: energy cannot be computed'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!Check overlap
    MaxOverlap=ZERO
    MinOverlap=huge(MinOverlap)/2
    AverOverlap=ZERO
    k=0
    do i=1,cbs
      do j=i+1,cbs
        if (abs(Glob_S(j,i))>abs(MaxOverlap)) MaxOverlap=Glob_S(j,i)
        if (abs(Glob_S(j,i))<abs(MinOverlap)) MinOverlap=Glob_S(j,i)
        AverOverlap=AverOverlap+abs(Glob_S(j,i))
      enddo
    enddo
    if (cbs>1) then
      AverOverlap=AverOverlap/(cbs*(cbs-1)/TWO)
    else
      MaxOverlap=ZERO
      MinOverlap=ZERO
      AverOverlap=ZERO
    endif

!Check linear coefficients
    Min_c=huge(Min_c)/2
    Max_c=ZERO
    Average_c=ZERO
    do i=1,cbs
      if (abs(Glob_c(i))>abs(Max_c)) Max_c=Glob_c(i)
      if (abs(Glob_c(i))<abs(Min_c)) Min_c=Glob_c(i)
      Average_c=Average_c+abs(Glob_c(i))
    enddo
    Average_c=Average_c/cbs

    if (Glob_ProcID==0) then
      write(*,*) ' done'
      write(*,*) '========== Summary after separation: ==========='
      write(*,*) 'Maximal by magnitude overlap         =',MaxOverlap
      write(*,*) 'Minimal by magnitude overlap         =',MinOverlap
      write(*,*) 'Average absolute value of overlap    =',AverOverlap
      write(*,*) 'Maximal by magnitude lin. coeff.     =',Max_c
      write(*,*) 'Minimal by magnitude lin. coeff.     =',Min_c
      write(*,*) 'Average absolute value of lin. coeff.=',Average_c
      write(*,*) 'Basis size after separation',cbs
      write(*,*) 'Energy after separation    ',Evalue
      write(*,*) '================================================'
    endif

    do i=1,cbs
      Glob_History(i)%Energy=ZERO
      Glob_History(i)%CyclesDone=0
      Glob_History(i)%InitFuncAtLastStep=0
      Glob_History(i)%NumOfEnergyEvalDuringFullOpt=0
    enddo
    Glob_History(cbs)%Energy=Evalue

    Glob_CurrEnergy=Evalue
    Glob_LastEigvalTol=1.0E+35_wp
    Glob_BestEigvalTol=1.0E+35_wp
    Glob_WorstEigvalTol=1.0E-35_wp

    ch_temp=Glob_DataFileName
    Glob_DataFileName=FileName
    if (Glob_ProcID==0) call SaveResults(Sort='no')
    Glob_DataFileName=ch_temp

    if (allocated(ActiveFunction)) deallocate(ActiveFunction)

!Deallocate workspace for DSYGVX
    if (Method=='G') then
      deallocate(Glob_IWorkForDSYGVX)
      deallocate(Glob_WorkForDSYGVX)
    else
      call ClearQWorkspace()
    endif

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    deallocate(Glob_diagS)
    if (Method=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) then
      i=len_trim(FileName)
      write(*,*) 'New basis has been saved in file',FileName(1:i)
      write(*,*) 'Program will now stop'
    endif

    call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

  end subroutine SeparateFuncLargeCoeff

  subroutine SaveHSWF(FileName1,FileName2,FileName3,FileName4,GSEPSolMethod)
!Subroutine SaveHSWF computes the Hamiltonian and overlap matrices, as well
!as the eigenvector corresponding to the normalized wave function and save
!them into files. Depending on the argument GSEPsolMethod, it can use
!either LAPACK subroutine DSYGVX or the inverse iteration method to solve GSEP.
!Input parameters:
!  FileName1 - the name of the file where the Hamiltonian matrix will be saved.
!              If FileName1 is equal to 'none','NONE', or 'None' then nothing
!              is saved.
!  FileName2 - the name of the file where the overlap matrix will be saved.
!              If FileName2 is equal to 'none','NONE', or 'None' then nothing
!              is saved.
!  FileName3 - the name of the file where the eigenvector (linear variational
!              parameters) will be saved. If FileName3 is equal to 'none','NONE',
!              or 'None' then nothing is saved.
!  FileName4 - the name of the file where the entire wave function (both
!              linear and nonlinear variational parameters) will be saved. If
!              FileName4 is equal to 'none','NONE',or 'None' then nothing is saved.
!  GSEPsolMethod - can be 'G', 'I', or 'Q'. It defines the method used to solve GSEP
!The format of the files that contains the Hamiltonian and overlap is such that each
!matrix element is placed in a separate line and is preceded by two integer indices, e.g.
!    23  78   0.123456789E+02
!The file that contains the eigenvector also contains one entry per line, preceded by its
!index, e.g.
!    57  0.123456789E+02
!The file that contains the wave function consists of a short header (containing the
!number of particles, masses, charges, and Young operator) and then
!lines containing the index i, linear coefficient of basis function i, and, after a
!colon, parameters of that function:
!   57  0.123456789E+02  :  0.987654321E+02 .........

!Parameters:
    character(Glob_FileNameLength),intent(in) :: FileName1
    character(Glob_FileNameLength),intent(in) :: FileName2
    character(Glob_FileNameLength),intent(in) :: FileName3
    character(Glob_FileNameLength),intent(in) :: FileName4
    character(1)        ::    GSEPSolMethod

!Local variables:
    integer        i,j
    integer        n,np,npt,cbs
    integer        ErrorCode
    logical        IsHNeeded,IsSNeeded,IsEVNeeded,IsWFNeeded
    logical        IsSwapFileOK
    integer        BlockSizeForDSYGVX
    integer        NumOfEigvecs,NumOfEigvalsFound
    real(wp)    Evalue
    real(wp),allocatable,dimension(:)      :: Eigvals
    real(wp),allocatable,dimension(:,:)    :: Eigvecs
    integer,allocatable,dimension(:)          :: IFAIL
    integer        NumOfIterations

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine SaveHSWF started'
      write(*,*) 'Number of basis functions',Glob_CurrBasisSize
      write(*,*) 'GSEP solution method ',GSEPsolMethod
    endif
    if ((GSEPsolMethod/='G').and.(GSEPsolMethod/='I').and. &
        (GSEPsolMethod/='Q')) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0190 in SaveHSWF: wrong GSEP solution method'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!Setting the values of some global and local variables
    Glob_GSEPSolutionMethod=GSEPsolMethod
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize
    if (GSEPsolMethod=='G') NumOfEigvecs=min(cbs,Glob_WhichEigenvalue+10)
    if (GSEPsolMethod=='I') NumOfEigvecs=1
    if (GSEPsolMethod=='Q') NumOfEigvecs=1

!Setting logical variables that determine if everything (H, S, eigenvector, wave function)
!needs to be saved
    if ((FileName1==' ').or.(FileName1=='none').or. &
        (FileName1=='NONE').or.(FileName1=='None')) then
      IsHNeeded=.false.
    else
      IsHNeeded=.true.
    endif
    if ((FileName2==' ').or.(FileName2=='none').or. &
        (FileName2=='NONE').or.(FileName2=='None')) then
      IsSNeeded=.false.
    else
      IsSNeeded=.true.
    endif
    if ((FileName3==' ').or.(FileName3=='none').or. &
        (FileName3=='NONE').or.(FileName3=='None')) then
      IsEVNeeded=.false.
    else
      IsEVNeeded=.true.
    endif
    if ((FileName4==' ').or.(FileName4=='none').or. &
        (FileName4=='NONE').or.(FileName4=='None')) then
      IsWFNeeded=.false.
    else
      IsWFNeeded=.true.
    endif

!Allocate global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (GSEPsolMethod=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    if (GSEPsolMethod=='I') allocate(Glob_invD(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (GSEPsolMethod=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Allocate workspace for subroutine GSEPIIS
    if (GSEPsolMethod=='I') then
      allocate(Glob_WorkForGSEPIIS(cbs))
      allocate(Glob_LastEigvector(cbs))
      Glob_LastEigvector(1:cbs)=ONE
    endif

!Allocate local arrays
    if (GSEPsolMethod=='G') then
      allocate(Eigvals(NumOfEigvecs))
      allocate(Eigvecs(cbs,NumOfEigvecs))
      allocate(IFAIL(cbs))
    endif

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a52)',advance='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) 'done'
    endif

    if (GSEPSolMethod=='G') then
      do i=1,cbs
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,cbs
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo

      !Saving matrices H and S:
      if (Glob_ProcID==0) then
        if (IsHNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving the Hamiltonian matrix...'
          open(2,file=FileName1)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              call writerealadv(2,Glob_H(i,j))
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
        if (IsSNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving the overlap matrix...'
          open(2,file=FileName2)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              call writerealadv(2,Glob_S(i,j))
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
      endif

      if ((IsEVNeeded).or.(IsWFNeeded)) then
        if (Glob_ProcID==0) then
          write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
          call DSYGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
                      ZERO,ZERO,1,NumOfEigvecs,Glob_AbsTolForDSYGVX, &
                      NumOfEigvalsFound,Eigvals,Eigvecs,cbs,Glob_WorkForDSYGVX,Glob_LWorkForDSYGVX, &
                      Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
          !SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!         VL, VU, IL, IU, ABSTOL,
!         M, W, Z, LDZ, WORK, LWORK,
!         IWORK, IFAIL, INFO )
        endif
        call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        if (ErrorCode/=0) then
          if (Glob_ProcID==0) then
            write(*,*) 'failed'
            write(*,*) &
              'Error EC0191 in SaveHSWF: routine DSYGVX failed with error code',ErrorCode
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif

        !sending the eigenvalue and the eigenvector to all processes
        if (Glob_ProcID==0) then
          Evalue=Eigvals(Glob_WhichEigenvalue)
          Glob_c(1:cbs)=Eigvecs(1:cbs,Glob_WhichEigenvalue)
        endif
        call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(Glob_c,cbs,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        Glob_CurrEnergy=Evalue

        if (Glob_ProcID==0) then
          write(*,*) 'done'
          write(*,*) 'Energy: ',Evalue
        endif

      endif

    endif !if (GSEPSolMethod=='G')

    if (GSEPSolMethod=='I') then
      !Saving matrices H and S:
      if (Glob_ProcID==0) then
        if (IsHNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving Hamiltonian matrix...'
          open(2,file=FileName1)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              if (i==j) call writerealadv(2,Glob_H(i,j)+Glob_ApproxEnergy)
              if (i>j) call writerealadv(2,Glob_H(i,j)+Glob_ApproxEnergy*Glob_S(i,j))
              if (i<j) call writerealadv(2,Glob_H(j,i)+Glob_ApproxEnergy*Glob_S(j,i))
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
        if (IsSNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving overlap matrix...'
          open(2,file=FileName2)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              if (i==j) then
                call writerealadv(2,ONE)
              else
                call writerealadv(2,Glob_S(i,j))
              endif
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
      endif

      if ((IsEVNeeded).or.(IsWFNeeded)) then
        if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
        if (cbs==1) then
          Glob_CurrEnergy=Glob_diagH(1)
          NumOfIterations=1
          ErrorCode=0
        else
          call GSEPIIS(1,cbs,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim, &
                       Glob_ApproxEnergy,Glob_LastEigvector,Glob_WorkForGSEPIIS,Glob_EigvalTol, &
                       Evalue,Glob_c,Glob_LastEigvalTol,Glob_MaxIterForGSEPIIS, &
                       0,NumOfIterations,ErrorCode)
          !GSEPIIS(k,n,M,nM,invD,B,nB, &
          !        apprlambda,v,w,Tol, &
          !        lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
          if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol
          if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol
          call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_BCAST(Glob_c,cbs,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        endif
        Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
        Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
        Glob_CurrEnergy=Evalue
        if (ErrorCode/=0) then
          if (Glob_ProcID==0) then
            write(*,*) 'failed'
            write(*,*) 'Error EC0192 in SaveHSWF: the energy cannot be computed'
          endif
          call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
        endif
        !print the energy
        if (Glob_ProcID==0) then
          write(*,*) 'done'
          write(*,*) 'Energy: ',Evalue
        endif
      endif
    endif

    if (GSEPSolMethod=='Q') then
      !Q physical matrices already use the canonical normalized lower
      !triangle. Materialize conceptual symmetry only while writing text; do
      !not overwrite the unused upper triangles in memory.
      if (Glob_ProcID==0) then
        if (IsHNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving Hamiltonian matrix...'
          open(2,file=FileName1)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              call writerealadv(2,QCanonicalMatrixElement(Glob_H,i,j))
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
        if (IsSNeeded) then
          write(*,'(1x,a)',advance='no') 'Saving overlap matrix...'
          open(2,file=FileName2)
          do i=1,cbs
            do j=1,cbs
              write(2,'(1x,i6,1x,i6,1x)',advance='no') i,j
              call writerealadv(2,QCanonicalMatrixElement(Glob_S,i,j))
            enddo
          enddo
          close(2)
          write(*,*) 'done'
        endif
      endif

      if (IsEVNeeded.or.IsWFNeeded) then
        if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') &
          'Solving eigenvalue problem...'
        Glob_c=ONE
        call PrepareQWorkspace(cbs,cbs,1,ErrorCode)
        Q_Workspace%MatricesAreCanonical=.true.
        if (ErrorCode==Q_METHOD_SUCCESS) call FactorizeQFresh(ErrorCode)
        if (ErrorCode==Q_METHOD_SUCCESS) call SolveQ(Evalue,ErrorCode)
        if (ErrorCode/=Q_METHOD_SUCCESS) then
          if (Glob_ProcID==0) then
            write(*,*) 'failed'
            write(*,'(1x,a,1x,i0)') &
              'Error EC0193 in SaveHSWF: Q energy cannot be computed, status',ErrorCode
          endif
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif
        Glob_CurrEnergy=Evalue
        if (Glob_ProcID==0) then
          write(*,*) 'done'
          write(*,*) 'Energy: ',Evalue
        endif
      endif
    endif

!Saving the eigenvector
    if ((IsEVNeeded).and.(Glob_ProcID==0)) then
      write(*,'(1x,a)',advance='no') 'Saving eigenvector...'
      open(2,file=FileName3)
      do i=1,cbs
        write(2,'(1x,i6,1x)',advance='no') i
        call writerealadv(2,Glob_c(i))
      enddo
      close(2)
      write(*,*) 'done'
    endif

!Saving the wave function
    if ((IsWFNeeded).and.(Glob_ProcID==0)) then
      write(*,'(1x,a)',advance='no') 'Saving wave function...'
      open(2,file=FileName4)
      write(2,'(1x,a)') 'RG_2D WAVE FUNCTION FILE'
      write(2,'(1x,a9,1x,i6)') 'PARTICLES',Glob_n+1
      write(2,'(1x,a6)',advance='no') 'MASSES'
      call writerealarradv(2,Glob_Mass,Glob_n+1)
      write(2,'(1x,a7)',advance='no') 'CHARGES'
      call writereal(2,Glob_PseudoCharge0)
      call writerealarradv(2,Glob_PseudoCharge,Glob_n)
      j=len_trim(Glob_YOperatorString)
      write(2,'(1x,a8)',advance='no') 'SYMMETRY'
      call writestringadv(2,Glob_YOperatorString,j)
      write(2,'(1x,a10,1x,i6)') 'BASIS_SIZE',Glob_CurrBasisSize
      write(2,'(1x,a14)',advance='no') 'CURRENT_ENERGY'
      call writerealadv(2,Glob_CurrEnergy)
      write(2,*) '=============================='
      do i=1,cbs
        write(2,'(1x,i6,1x)',advance='no') i
        call writereal(2,Glob_c(i))
        write(2,'(1x,a1,1x)',advance='no') ':'
        write(2,'(i6,1x,i6,1x)',advance='no') Glob_Index(i,1),Glob_Index(i,2)
        call writerealarradv(2,Glob_NonlinParam(1:Glob_npt,i),Glob_npt)
      enddo
      close(2)
      write(*,*) 'done'
    endif

    if (GSEPSolMethod=='G') then
      deallocate(IFAIL)
      deallocate(Eigvecs)
      deallocate(Eigvals)
    endif

    if (GSEPsolMethod=='I') then
      deallocate(Glob_LastEigvector)
      deallocate(Glob_WorkForGSEPIIS)
    endif
    if (GSEPSolMethod=='Q') call ClearQWorkspace()

!deallocate workspace for DSYGVX
    if (GSEPSolMethod=='G') then
      deallocate(Glob_WorkForDSYGVX)
      deallocate(Glob_IWorkForDSYGVX)
    endif

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    if (GSEPSolMethod=='I') deallocate(Glob_invD)
    deallocate(Glob_diagS)
    if (GSEPSolMethod=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if (Glob_ProcID==0) write (*,*) 'Routine SaveHSWF finished'

  end subroutine SaveHSWF

  subroutine ExpectationValues(SymmAdaptMethod,FileName1,FileName2,FileName3,FileName4,GSEPSolMethod)
!ExpectationValues computes expectation values in the basis of
!Glob_CurrBasisSize functions. Depending on the argument GSEPsolMethod,
!it can use either LAPACK subroutine DSYGVX or the inverse iteration method to
!solve GSEP.
!Input parameters:
!  SymmAdaptMethod  - defines how the expectation values should be
!calculated. If SymmAdaptMethod=1 then Y^{\dagger}Y operator is applied
!to the ket. This option is faster but requires manual symmetrization
!of the obtained expectation values of such operators as r_{ij} in the
!case of some Hamiltonian symmetries (like charge reversal symmetry).
!If SymmAdaptMethod=2 then Y operator is applied to both bra and ket.
!Since the square of the number of terms in Y operator may be significantly
!larger than the number of terms in Y^{\dagger}Y operator this option is slower.
!But in this case all the expectation values are correct and ready to use
!immediately.
!  FileName1 - the name of the file that defines grid for correlation functions.
!              If FileName1 is equal to 'none','NONE', or 'None' then correlation
!              functions are not computed.
!  FileName2 - the name of the file where correlation functions will be stored
!  FileName3 - the name of the file that defines grid for particle densities.
!              If FileName3 is equal to 'none','NONE', or 'None' then particle
!              densities are not computed.
!  FileName4 - the name of the file where particle densities will be stored
!  GSEPsolMethod - can be 'G', 'I', or 'Q'. It defines the method used to solve GSEP

!Parameters:
    integer,intent(in)  ::    SymmAdaptMethod
    character(Glob_FileNameLength),intent(in) :: FileName1
    character(Glob_FileNameLength),intent(in) :: FileName2
    character(Glob_FileNameLength),intent(in) :: FileName3
    character(Glob_FileNameLength),intent(in) :: FileName4
    character(1)        ::    GSEPSolMethod

!Local variables:
    integer        i,j,k,kk,counter,a,b,c,d,a1,b1
    integer        n,np,npt,cbs
    integer        OpenFileErr,ErrorCode
    logical        IsSwapFileOK
    integer        BlockSizeForDSYGVX
    integer        NumOfEigvecs,NumOfEigvalsFound
    real(wp)    Evalue
    real(wp),allocatable,dimension(:)      :: Eigvals
    real(wp),allocatable,dimension(:,:)    :: Eigvecs
    integer,allocatable,dimension(:)          :: IFAIL
    integer        NumOfExpcVals,NumOfIterations
    real(wp)    factor
    real(wp)    temp1,temp2
    real(wp),allocatable,dimension(:,:)    ::  IdentityPerm
    real(wp)    beta,mu
    logical        AreCorrFuncNeeded,ArePartDensNeeded
    logical        IsFile1OK,IsFile3OK

!Local variables used to store temporary data
!associated with certain expectation values
!Local variables used to store temporary data
!associated with certain expectation values
    integer                                    :: NumCFGridPoints
    real(wp),allocatable,dimension(:,:)     :: CFGrid
    real(wp),allocatable,dimension(:,:)     :: CFkl
    real(wp),allocatable,dimension(:,:)     :: CF
    integer                                    :: NumDensGridPoints
    real(wp),allocatable,dimension(:,:)     :: DensGrid
    real(wp),allocatable,dimension(:,:)     :: Denskl
    real(wp),allocatable,dimension(:,:)     :: Dens
    integer                                    :: NumOfCFAndDensExpVals
    real(wp),allocatable,dimension(:)       :: CFDMEkl_s
    real(wp),allocatable,dimension(:)       :: MEkl,MEkl_s,MEkl1,MEkl_s1,MEkl2,MEkl_s2
    real(wp)                                :: Hkl,Skl,Tkl,Vkl
    real(wp)                                :: Hkl1,Skl1,Tkl1,Vkl1
    real(wp)                                :: Hkl2,Skl2,Tkl2,Vkl2
    real(wp)                                :: MVkl,Darwinkl,drach_Darwinkl,OOkl
    real(wp)                                :: MVkl1,drach_MVkl1,Darwinkl1,drach_Darwinkl1,OOkl1
    real(wp)                                :: MVkl2,drach_MVkl2,Darwinkl2,drach_Darwinkl2,OOkl2
    real(wp)                                :: H,S,T,V,MV,drach_MV1,drach_MV2,Darwin,drach_Darwin,OO
    real(wp),allocatable,dimension(:,:)     :: rm2kl,rmkl,rkl,r2kl,deltarkl,drach_deltarkl,prvalkl
    real(wp),allocatable,dimension(:,:)     :: rm2kl1,rmkl1,rkl1,r2kl1,deltarkl1,drach_deltarkl1,prvalkl1
    real(wp),allocatable,dimension(:,:)     :: rm2kl2,rmkl2,rkl2,r2kl2,deltarkl2,drach_deltarkl2,prvalkl2
    real(wp),allocatable,dimension(:,:)     :: rm2,rm,r,r2,deltar,drach_deltar,prval
    real(wp),allocatable,dimension(:,:,:,:) :: rmrmkl,rmrmkl1,rmrmkl2

! spin-dependent stuff
    integer :: nFactorial, spinDependentValuesNeeded
    real(wp) :: Skk
    real(wp), allocatable, dimension(:, :, :) :: SziME, ketYMatrix, drach_SSFMatrix, &
                                                    drach_AnihMatrix, SSFMatrix, AnihMatrix
    real(wp), allocatable, dimension(:, :, :) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient
    real(wp), allocatable, dimension(:, :) :: SSFmassChargeCoefficient, SSNCmassChargeCoefficient, AnihMassChargeCoefficient
    real(wp), allocatable, dimension(:, :, :, :) :: SiSjME, SSNCspinME
    real(wp), allocatable, dimension(:) :: SO1kl, SO2kl, SO1, SO2, AMM1, AMM2, AMM1Fin, AMM2Fin, &
                                              AMM1kl, AMM2kl, AMM1Finkl, AMM2Finkl, SSNC, SSNCkl, &
                                              drach_SSF, drach_SSFe, drach_Anih, SSF, SSFe, Anih
    real(wp), allocatable, dimension(:) :: parityFactor, diagS
    real(wp), allocatable, dimension(:) :: spinFreeME
    integer :: positronPosition, numberOfSpinFunctions

! One can set this flag to zero to disable everything introduced by DT
    spinDependentValuesNeeded = 1

    if (Glob_ProcID==0) then
      write(*,*)
      write(*,*) 'Routine ExpectationValues started'
      write(*,*) 'Number of basis functions',Glob_CurrBasisSize
      write(*,*) 'GSEP solution method ',GSEPsolMethod
    endif
    if ((GSEPsolMethod/='G').and.(GSEPsolMethod/='I').and. &
        (GSEPsolMethod/='Q')) then
      if (Glob_ProcID==0) then
        write(*,*) 'Error EC0195 in ExpectationValues: wrong GSEP solution method'
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    endif

!Setting the values of some global and local variables
    Glob_GSEPSolutionMethod=GSEPsolMethod
    Glob_OverlapPenaltyAllowed=.false.
    Glob_HSLeadDim=Glob_CurrBasisSize
    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    Glob_HSBuffLen=max(min(Glob_CurrBasisSize*(Glob_CurrBasisSize+1)/2,1000),30*Glob_CurrBasisSize)
    cbs=Glob_CurrBasisSize
    if (GSEPsolMethod=='G') NumOfEigvecs=min(cbs,Glob_WhichEigenvalue+10)
    if (GSEPsolMethod=='I') NumOfEigvecs=1
    if (GSEPsolMethod=='Q') NumOfEigvecs=1

!Setting logical variables that determine whether correlation
!functions and particle densities need to be computed
    if ((FileName1==' ').or.(FileName1=='none').or. &
        (FileName1=='NONE').or.(FileName1=='None')) then
      AreCorrFuncNeeded=.false.
    else
      AreCorrFuncNeeded=.true.
    endif
    if ((FileName3==' ').or.(FileName3=='none').or. &
        (FileName3=='NONE').or.(FileName3=='None')) then
      ArePartDensNeeded=.false.
    else
      ArePartDensNeeded=.true.
    endif

!If correlation functions and/or particle densities are needed
!then we open files FileName1 and FileName3 that contain grids for
!correlation functions and particle densities respectively

!Here we determine the number of grid points for correlation
!function calculation
    if (AreCorrFuncNeeded) then
      if (Glob_ProcID==0) then
        IsFile1OK=.true.
        open(1,file=FileName1,status='old',iostat=OpenFileErr)
        if (OpenFileErr==0) then
          NumCFGridPoints=0
          do while (OpenFileErr==0)
            read (1,*,iostat=OpenFileErr) temp1,temp2
            NumCFGridPoints=NumCFGridPoints+1
          enddo
          NumCFGridPoints=NumCFGridPoints-1
          if ((OpenFileErr>0).or.(NumCFGridPoints==0)) then
            !Improper data in the input file
            write(*,*) 'Error EC0196 in ExpectationValues: improper data in file ',FileName1
            IsFile1OK=.false.
          endif
          close(1)
        else
          IsFile1OK=.false.
        endif
      endif
      call MPI_BCAST(IsFile1OK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(NumCFGridPoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      !stop if there were problems with correlation function grid file
      if (.not.IsFile1OK) then
        if (Glob_ProcID==0) then
          write(*,*) 'Error EC0197 in ExpectationValues: cannot open CF grid file ',FileName1
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
    endif

!Here we determine the number of grid points for particle density
!calculations
    if (ArePartDensNeeded) then
      if (Glob_ProcID==0) then
        IsFile3OK=.true.
        open(1,file=FileName3,status='old',iostat=OpenFileErr)
        if (OpenFileErr==0) then
          NumDensGridPoints=0
          do while (OpenFileErr==0)
            read (1,*,iostat=OpenFileErr) temp1,temp2
            NumDensGridPoints=NumDensGridPoints+1
          enddo
          NumDensGridPoints=NumDensGridPoints-1
          if ((OpenFileErr>0).or.(NumDensGridPoints==0)) then
            !Improper data in the input file
            write(*,*) 'Error EC0198 in ExpectationValues: improper data in file ',FileName3
            IsFile3OK=.false.
          endif
          close(1)
        else
          IsFile3OK=.false.
        endif
      endif
      call MPI_BCAST(IsFile3OK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(NumDensGridPoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      !stop if there were problems with particle density grid file
      if (.not.IsFile3OK) then
        if (Glob_ProcID==0) then
          write(*,*) 'Error EC0199 in ExpectationValues: cannot open density grid file ',FileName3
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
    endif

!Allocate arrays that will be used to store grid points and the
!function values for correlation function and particle density
!calculations

    NumOfCFAndDensExpVals=0

    if (AreCorrFuncNeeded) then
      allocate(CFGrid(2,NumCFGridPoints))
      NumOfCFAndDensExpVals=NumOfCFAndDensExpVals+NumCFGridPoints*n*(n+1)/2
      allocate(CFkl(n*(n+1)/2,NumCFGridPoints))
      allocate(CF(n*(n+1)/2,NumCFGridPoints))
    else
      !allocate just one point to have a valid pointer
      allocate(CFGrid(2,1))
      allocate(CFkl(1,1))
    endif

    if (ArePartDensNeeded) then
      allocate(DensGrid(2,NumDensGridPoints))
      NumOfCFAndDensExpVals=NumOfCFAndDensExpVals+NumDensGridPoints*(n+1)
      allocate(Denskl(n+1,NumDensGridPoints))
      allocate(Dens(n+1,NumDensGridPoints))
    else
      !allocate just one point to have a valid pointer
      allocate(DensGrid(2,1))
      allocate(Denskl(1,1))
    endif

    if((AreCorrFuncNeeded).or.(ArePartDensNeeded)) then
      allocate(CFDMEkl_s(NumOfCFAndDensExpVals))
    endif

!Now we open files FileName1 and FileName3 again, but this time
!we read the data from them into arrays CFGrid and DensGrid
!It is assumed that the data in FileName1 and FileName3 consists
!of two columns that define points in cylindrical coordinates.
!The first column is $\xhi_\rho$ and the second one is $\xi_z$
    if (AreCorrFuncNeeded) then
      if (Glob_ProcID==0) then
        open(1,file=FileName1,status='old')
        do i=1,NumCFGridPoints
          read (1,*) CFGrid(1,i),CFGrid(2,i)
        enddo
        close(1)
      endif
      call MPI_BCAST(CFGrid,2*NumCFGridPoints,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
    if (ArePartDensNeeded) then
      if (Glob_ProcID==0) then
        open(1,file=FileName3,status='old')
        do i=1,NumDensGridPoints
          read (1,*) DensGrid(1,i),DensGrid(2,i)
        enddo
        close(1)
      endif
      call MPI_BCAST(DensGrid,2*NumDensGridPoints,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

!Allocate global arrays
    allocate(Glob_H(cbs,cbs))
    allocate(Glob_S(cbs,cbs))
    if (GSEPsolMethod=='G') allocate(Glob_diagH(cbs))
    allocate(Glob_diagS(cbs))
    if (GSEPsolMethod=='I') allocate(Glob_invD(cbs))
    allocate(Glob_c(cbs))
    allocate(Glob_HklBuff1(Glob_HSBuffLen))
    allocate(Glob_HklBuff2(Glob_HSBuffLen))
    allocate(Glob_SklBuff1(Glob_HSBuffLen))
    allocate(Glob_SklBuff2(Glob_HSBuffLen))

!Allocate workspace for DSYGVX
    if (GSEPsolMethod=='G') then
      BlockSizeForDSYGVX=ILAENV(1,'DSYTRD','VIU',cbs,cbs,cbs,cbs)
      Glob_LWorkForDSYGVX=max((BlockSizeForDSYGVX+3)*cbs,8*cbs)
      allocate(Glob_WorkForDSYGVX(Glob_LWorkForDSYGVX))
      allocate(Glob_IWorkForDSYGVX(5*cbs))
    endif

!Allocate workspace for subroutine GSEPIIS
    if (GSEPsolMethod=='I') then
      allocate(Glob_WorkForGSEPIIS(cbs))
      allocate(Glob_LastEigvector(cbs))
      Glob_LastEigvector(1:cbs)=ONE
    endif

!Allocate local arrays
    if (GSEPsolMethod=='G') then
      allocate(Eigvals(NumOfEigvecs))
      allocate(Eigvecs(cbs,NumOfEigvecs))
      allocate(IFAIL(cbs))
    endif

    allocate(IdentityPerm(n,n))
    IdentityPerm(1:n,1:n)=ZERO
    do i=1,n
      IdentityPerm(i,i)=ONE
    enddo

!rm2          ME(1:n*(n+1)/2)
!rm           ME(n*(n+1)/2+1:2*n*(n+1)/2)
!r            ME(2*n*(n+1)/2+1:3*n*(n+1)/2)
!r2           ME(3*n*(n+1)/2+1:4*n*(n+1)/2)
!deltar       ME(4*n*(n+1)/2+1:5*n*(n+1)/2)
!drach_deltar ME(5*n*(n+1)/2+1:6*n*(n+1)/2)
!prval        ME(6*n*(n+1)/2+1:7*n*(n+1)/2)
!H            ME(7*n*(n+1)/2+1)
!S            ME(7*n*(n+1)/2+2)
!T            ME(7*n*(n+1)/2+3)
!V            ME(7*n*(n+1)/2+4)
!MV           ME(7*n*(n+1)/2+5)
!drach_MVkl   ME(7*n*(n+1)/2+6)
!Darwin       ME(7*n*(n+1)/2+7)
!drach_Darwin ME(7*n*(n+1)/2+8)
!OO           ME(7*n*(n+1)/2+9)
!rmrmkl       ME( 7*n*(n+1)/2+10 + (3*n**4+10*n**3+9*n**2+2*n)/24 : 6*n*(n+1)/2+9 + (3*n**4+10*n**3+9*n**2+2*n)/12 )
    NumOfExpcVals = 9 + 7*n*(n+1)/2 + (3*n**4+10*n**3+9*n**2+2*n)/12

    allocate(MEkl(NumOfExpcVals))
    allocate(MEkl_s(NumOfExpcVals))

    allocate(MEkl1(NumOfExpcVals))
    allocate(MEkl_s1(NumOfExpcVals))
    allocate(MEkl2(NumOfExpcVals))
    allocate(MEkl_s2(NumOfExpcVals))

    allocate(rm2kl(n,n))
    allocate(rm2(n,n))

    allocate(rm2kl1(n,n))
    allocate(rm2kl2(n,n))

    allocate(rmkl(n,n))
    allocate(rm(n,n))

    allocate(rmkl1(n,n))
    allocate(rmkl2(n,n))

    allocate(rkl(n,n))
    allocate(r(n,n))

    allocate(rkl1(n,n))
    allocate(rkl2(n,n))

    allocate(r2kl(n,n))
    allocate(r2(n,n))

    allocate(r2kl1(n,n))
    allocate(r2kl2(n,n))

    allocate(deltarkl(n,n))
    allocate(deltar(n,n))

    allocate(deltarkl1(n,n))
    allocate(deltarkl2(n,n))

    allocate(drach_deltarkl(n,n))
    allocate(drach_deltar(n,n))

    allocate(drach_deltarkl1(n,n))
    allocate(drach_deltarkl2(n,n))

    allocate(prvalkl(n,n))
    allocate(prval(n,n))

    allocate(prvalkl1(n,n))
    allocate(prvalkl2(n,n))

    allocate(rmrmkl(n,n,n,n))

    allocate(rmrmkl1(n,n,n,n))
    allocate(rmrmkl2(n,n,n,n))

    call ReadSwapFileAndDistributeData(IsSwapFileOK)

    if (.not.IsSwapFileOK) then
      if (Glob_ProcID==0) write(*,'(1x,a52)',advance='no') &
        'Computing Hamiltonian and overlap matrix elements...'
      call ComputeMatElem(1,cbs)
      if (Glob_ProcID==0) write(*,*) 'done'
    endif

    if (GSEPSolMethod=='G') then
      do i=1,cbs
        do j=1,i-1
          Glob_H(j,i)=Glob_H(i,j)
        enddo
        Glob_H(i,i)=Glob_diagH(i)
      enddo
      do i=1,cbs
        do j=1,i-1
          Glob_S(j,i)=Glob_S(i,j)
        enddo
        Glob_S(i,i)=ONE
      enddo

      if (Glob_ProcID==0) then
        write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
        call DSYGVX(1,'V','I','U',cbs,Glob_H,Glob_HSLeadDim,Glob_S,Glob_HSLeadDim,  &
                    ZERO,ZERO,1,NumOfEigvecs,Glob_AbsTolForDSYGVX, &
                    NumOfEigvalsFound,Eigvals,Eigvecs,cbs,Glob_WorkForDSYGVX,Glob_LWorkForDSYGVX, &
                    Glob_IWorkForDSYGVX,IFAIL,ErrorCode)
        !SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!       VL, VU, IL, IU, ABSTOL,
!       M, W, Z, LDZ, WORK, LWORK,
!       IWORK, IFAIL, INFO )
      endif
      call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (ErrorCode/=0) then
        if (Glob_ProcID==0) then
          write(*,*) 'failed'
          write(*,*) &
            'Error EC0200 in ExpectationValues: routine DSYGVX failed with error code',ErrorCode
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif

      !sending the eigenvalue and the eigenvector to all processes
      if (Glob_ProcID==0) then
        Evalue=Eigvals(Glob_WhichEigenvalue)
        Glob_c(1:cbs)=Eigvecs(1:cbs,Glob_WhichEigenvalue)
      endif
      call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_BCAST(Glob_c,cbs,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      Glob_CurrEnergy=Evalue

      !print the lower part of the spectrum
      if (Glob_ProcID==0) then
        write(*,*) 'done'
        write(*,*) 'Energy: ',Evalue
        write(*,*)
        write(*,*) 'Lowest eigenvalues:'
        do i=1,NumOfEigvalsFound
          write(*,*) i,' ',Eigvals(i)
        enddo
        write(*,*)
      endif
    endif !if (GSEPSolMethod=='G')

    if (GSEPSolMethod=='I') then
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') 'Solving eigenvalue problem...'
      if (cbs==1) then
        Glob_CurrEnergy=Glob_diagH(1)
        NumOfIterations=1
        ErrorCode=0
      else
        call GSEPIIS(1,cbs,Glob_H,Glob_HSLeadDim,Glob_invD,Glob_S,Glob_HSLeadDim, &
                     Glob_ApproxEnergy,Glob_LastEigvector,Glob_WorkForGSEPIIS,Glob_EigvalTol, &
                     Evalue,Glob_c,Glob_LastEigvalTol,Glob_MaxIterForGSEPIIS, &
                     0,NumOfIterations,ErrorCode)
        !GSEPIIS(k,n,M,nM,invD,B,nB, &
        !        apprlambda,v,w,Tol, &
        !        lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
        if (Glob_LastEigvalTol>Glob_WorstEigvalTol) Glob_WorstEigvalTol=Glob_LastEigvalTol
        if (Glob_LastEigvalTol>Glob_BestEigvalTol) Glob_BestEigvalTol=Glob_LastEigvalTol
        call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(Evalue,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(Glob_c,cbs,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      endif
      Glob_InvItTempCounter1=Glob_InvItTempCounter1+1
      Glob_InvItTempCounter2=Glob_InvItTempCounter2+NumOfIterations
      Glob_CurrEnergy=Evalue
      if (ErrorCode/=0) then
        if (Glob_ProcID==0) then
          write(*,*) 'failed'
          write(*,*) 'Error EC0201 in ExpectationValues: the energy cannot be computed'
        endif
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
      endif
      !print the energy
      if (Glob_ProcID==0) then
        write(*,*) 'done'
        write(*,*) 'Energy: ',Evalue
      endif
    endif

    if (GSEPSolMethod=='Q') then
      if (Glob_ProcID==0) write(*,'(1x,a29)',advance='no') &
        'Solving eigenvalue problem...'
      Glob_c=ONE
      call PrepareQWorkspace(cbs,cbs,1,ErrorCode)
      Q_Workspace%MatricesAreCanonical=.true.
      if (ErrorCode==Q_METHOD_SUCCESS) call FactorizeQFresh(ErrorCode)
      if (ErrorCode==Q_METHOD_SUCCESS) call SolveQ(Evalue,ErrorCode)
      if (ErrorCode/=Q_METHOD_SUCCESS) then
        if (Glob_ProcID==0) then
          write(*,*) 'failed'
          write(*,'(1x,a,1x,i0)') &
            'Error EC0202 in ExpectationValues: Q energy cannot be computed, status',ErrorCode
        endif
        call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
      endif
      Glob_CurrEnergy=Evalue
      if (Glob_ProcID==0) then
        write(*,*) 'done'
        write(*,*) 'Energy: ',Evalue
      endif
    endif

    if (Glob_ProcID==0) write(*,'(1x,a31)',advance='no') 'Computing expectation values...'

! in case we want to evaluate the spin-dependent operators mean values,
! we switch from spin-free to a regular formalism here
! note that we change global arrays Glob_YHYMatr and Glob_YHYCoeff here
    if (spinDependentValuesNeeded == 1) then

      nFactorial = 1
      do i = 2, n
        nFactorial = nFactorial * i
      enddo

      allocate(SOmassChargeCoefficient(n, n, 4))
      allocate(SSFmassChargeCoefficient(n, n))
      allocate(SSNCmassChargeCoefficient(n, n))
      allocate(AnihMassChargeCoefficient(n, n))
      allocate(AMMmassChargeCoefficient(n, n, 4))
      allocate(AMMFinmassChargeCoefficient(n, n, 4))
      allocate(parityFactor(nFactorial))

      allocate(ketYMatrix(1 : n, 1 : n, nFactorial))
      allocate(spinFreeME(nFactorial))
      allocate(SziME(n, 2, nFactorial))
      allocate(SiSjME(n, n, 2, nFactorial))
      allocate(SSNCspinME(n, n, 2, nFactorial))
      !allocate(SiSjCoeff(n, n, 2, nFactorial))

      call spinPreCalc(n, nFactorial, SziME, parityFactor, SSFmassChargeCoefficient,SSNCmassChargeCoefficient, &
            SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, AnihMassChargeCoefficient, ketYMatrix, &
                       Glob_YOperatorString, positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME, SSNCspinME)

      ! changing global variables here, care
      deallocate(Glob_YHYMatr, Glob_YHYCoeff)

      Glob_NumYHYTerms = nFactorial

      if (positronPosition > 0) then ! we do have a positron

        do i = 1, nFactorial
          if (nint(ketYMatrix(positronPosition, positronPosition, i)) == 0) then
            Glob_NumYHYTerms = Glob_NumYHYTerms - 1 ! exclude positron permutations
          endif
        enddo

      endif

      allocate(Glob_YHYMatr(n, n, Glob_NumYHYTerms))
      allocate(Glob_YHYCoeff(Glob_NumYHYTerms))

      k = 1
      do i = 1, nFactorial
        if (positronPosition > 0) then
          if (nint(ketYMatrix(positronPosition, positronPosition, i)) == 0) cycle
        endif

        Glob_YHYMatr(:, :, k) = ketYMatrix(:, :, i)
        Glob_YHYCoeff(k) = parityFactor(i) * spinFreeME(i)

        ! do kk = 1, numberOfSpinFunctions
        !   do j = 1, n
        !     do c = 1, n
        !       SiSjCoeff(j, c, kk, k) = parityFactor(i) * SiSjME(j, c, kk, i)
        !     enddo
        !   enddo
        ! enddo

        ! skipping positronic permutations
        do kk = 1, numberOfSpinFunctions
          do j = 1, n
            do c = 1, n
              SiSjME(j, c, kk, k) = parityFactor(i) * SiSjME(j, c, kk, i)
            enddo
          enddo
        enddo

        k = k + 1

      enddo

      allocate(diagS(cbs))
      diagS = ZERO

      ! we should recalculate mean values of a unity operator here (it should be proportional to the old values)
      Skk = ZERO
      do i = 1, cbs
        do a = 1, Glob_NumYHYTerms
          call OverlapMatrixElement_RG_2D(Glob_Index(i,1), Glob_Index(i,2), Glob_NonlinParam(1 : npt, i), &
                                       Glob_YHYMatr(1 : n, 1 : n, a), Skk)
          diagS(i) = diagS(i) + Glob_YHYCoeff(a) * Skk
        enddo ! Permutations from S_n
        Glob_diagS(i) = diagS(i)
      enddo

      allocate(drach_SSFMatrix(n, n, numberOfSpinFunctions))
      drach_SSFMatrix = ZERO
      allocate(SSFMatrix(n, n, numberOfSpinFunctions))
      SSFMatrix = ZERO

      allocate(drach_AnihMatrix(n, n, numberOfSpinFunctions))
      drach_AnihMatrix = ZERO
      allocate(AnihMatrix(n, n, numberOfSpinFunctions))
      AnihMatrix = ZERO

    endif

!main loop
    MEkl_s(1:NumOfExpcVals)=ZERO
    if (AreCorrFuncNeeded.or.ArePartDensNeeded) CFDMEkl_s(1:NumOfCFAndDensExpVals)=ZERO
    counter=0
    do i=1,cbs
      do j=1,i  !j=1,i
        counter=counter+1
        if (mod(counter,Glob_NumOfProcs)==Glob_ProcID) then
          if (i==j) then
            factor=Glob_c(i)*Glob_c(j)/sqrt(Glob_diagS(i)*Glob_diagS(i))  !1
          else
            factor=TWO*Glob_c(i)*Glob_c(j)/sqrt(Glob_diagS(i)*Glob_diagS(j))  !2
          endif
          if (SymmAdaptMethod==1) then
            do k=1,Glob_NumYHYTerms
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call MatrixElementsAll_RG_2D(Glob_Index(i,1),Glob_Index(i,2),Glob_Index(j,1),Glob_Index(j,2),  &
                                                Glob_NonlinParam(1:npt,i),Glob_NonlinParam(1:npt,j),                     &
                                                IdentityPerm,Glob_YHYMatr(1:n,1:n,k),Hkl,Skl,Tkl,Vkl,                    &
                                                rm2kl,rmkl,rkl,r2kl,deltarkl,drach_deltarkl,MVkl,drach_MVkl1,drach_MVkl2,   &
                                                Darwinkl,drach_Darwinkl,OOkl,rmrmkl,prvalkl,NumCFGridPoints,CFGrid,CFkl, &
                                                NumDensGridPoints,DensGrid,Denskl,AreCorrFuncNeeded,ArePartDensNeeded)
              c=0
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=rm2kl(b,a)!-rm2kl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=rmkl(b,a)!-rmkl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=rkl(b,a)!-rkl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=r2kl(b,a)!-r2kl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=deltarkl(b,a)!-deltarkl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=drach_deltarkl(b,a)!-drach_deltarkl1(b,a)
                enddo
              enddo
              do a=1,n
                do b=a,n
                  c=c+1; MEkl(c)=prvalkl(b,a)!-prvalkl1(b,a)
                enddo
              enddo
              c=c+1; MEkl(c)=Hkl!2-Hkl1
              c=c+1; MEkl(c)=Skl!2-Skl1
              c=c+1; MEkl(c)=Tkl!2-Tkl1
              c=c+1; MEkl(c)=Vkl!2-Vkl1
              c=c+1; MEkl(c)=MVkl!2-MVkl1
              c=c+1; MEkl(c)=drach_MVkl1!2-drach_MVkl1
              c=c+1; MEkl(c)=drach_MVkl2
              c=c+1; MEkl(c)=Darwinkl!2-Darwinkl1
              c=c+1; MEkl(c)=drach_Darwinkl!2-drach_Darwinkl1
              c=c+1; MEkl(c)=OOkl!2-OOkl1
              do a=1,n
                do b=a,n
                  do a1=a,n
                    do b1=a1,n
                      c=c+1; MEkl(c)=rmrmkl(a,b,a1,b1)!-rmrmkl1(a,b,a1,b1)
                    enddo
                  enddo
                enddo
              enddo

              do a=1,NumOfExpcVals
                MEkl_s(a)=MEkl_s(a)+factor*Glob_YHYCoeff(k)*MEkl(a)
              enddo

              if (spinDependentValuesNeeded == 1) then
                ! SSF term is special: it needs SiSj mean value with it
                ! not the spin-free value like the other terms here
                ! we build it from drachmanized deltas for each pair b, a
                do c = 1, numberOfSpinFunctions
                  do a = 1, n
                    do b = a + 1, n
                      drach_SSFMatrix(a, b, c) = drach_SSFMatrix(a, b, c) + factor * SiSjME(a, b, c, k) * drach_deltarkl(b, a)
                    enddo
                  enddo

                  do a = 1, n
                    do b = a + 1, n
                      SSFMatrix(a, b, c) = SSFMatrix(a, b, c) + factor * SiSjME(a, b, c, k) * deltarkl(b, a)
                    enddo
                  enddo

                  if (positronPosition > 0) then
                    ! with positron we also need annihilation correction
                    do a = 1, n
                      do b = a + 1, n
                        drach_AnihMatrix(a, b, c) = drach_AnihMatrix(a, b, c) + &
                                                    factor * (SiSjME(a, b, c, k) + THREE / FOUR * Glob_YHYCoeff(k)) &
                                                    * drach_deltarkl(b, a)
                      enddo
                    enddo

                    do a = 1, n
                      do b = a + 1, n
                        AnihMatrix(a, b, c) = AnihMatrix(a, b, c) + &
                                              factor * (SiSjME(a, b, c, k) + THREE / FOUR * Glob_YHYCoeff(k)) &
                                              * deltarkl(b, a)
                      enddo
                    enddo
                  endif

                enddo
              endif

              c=0
              if (AreCorrFuncNeeded) then
                do a=1,NumCFGridPoints
                  do b=1,n*(n+1)/2
                    c=c+1
                    CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YHYCoeff(k)*CFkl(b,a)
                  enddo
                enddo
              endif
              if (ArePartDensNeeded) then
                do a=1,NumDensGridPoints
                  do b=1,n+1
                    c=c+1
                    CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YHYCoeff(k)*Denskl(b,a)
                  enddo
                enddo
              endif

            enddo !k=1,Glob_NumYHYTerms
          endif !SymmAdaptMethod==1

          if (SymmAdaptMethod==2) then

            do k=1,Glob_NumYTerms
              do kk=1,Glob_NumYTerms
                call MatrixElementsAll_RG_2D(Glob_Index(i,1),Glob_Index(i,2),Glob_Index(j,1), Glob_Index(j,2),  &
                                                  Glob_NonlinParam(1:npt,i),Glob_NonlinParam(1:npt,j),                     &
                                                  Glob_YMatr(1:n,1:n,k),Glob_YMatr(1:n,1:n,kk),Hkl,Skl,Tkl,Vkl,            &
                                                  rm2kl,rmkl,rkl,r2kl,deltarkl,drach_deltarkl,MVkl,drach_MVkl1, drach_MVkl2,  &
                                                  Darwinkl,drach_Darwinkl,OOkl,rmrmkl,prvalkl,NumCFGridPoints,CFGrid,CFkl, &
                                                  NumDensGridPoints,DensGrid,Denskl,AreCorrFuncNeeded,ArePartDensNeeded)
                c=0
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=rm2kl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=rmkl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=rkl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=r2kl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=deltarkl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=drach_deltarkl(b,a)
                  enddo
                enddo
                do a=1,n
                  do b=a,n
                    c=c+1; MEkl(c)=prvalkl(b,a)
                  enddo
                enddo
                c=c+1; MEkl(c)=Hkl
                c=c+1; MEkl(c)=Skl
                c=c+1; MEkl(c)=Tkl
                c=c+1; MEkl(c)=Vkl
                c=c+1; MEkl(c)=MVkl
                c=c+1; MEkl(c)=drach_MVkl1
                c=c+1; MEkl(c)=drach_MVkl2
                c=c+1; MEkl(c)=Darwinkl
                c=c+1; MEkl(c)=drach_Darwinkl
                c=c+1; MEkl(c)=OOkl
                do a=1,n
                  do b=a,n
                    do a1=a,n
                      do b1=a1,n
                        c=c+1; MEkl(c)=rmrmkl(a,b,a1,b1)
                      enddo
                    enddo
                  enddo
                enddo

                do a=1,NumOfExpcVals
                  MEkl_s(a)=MEkl_s(a)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*MEkl(a)
                enddo

                c=0
                if (AreCorrFuncNeeded) then
                  do a=1,NumCFGridPoints
                    do b=1,n*(n+1)/2
                      c=c+1
                      CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*CFkl(b,a)
                    enddo
                  enddo
                endif
                if (ArePartDensNeeded) then
                  do a=1,NumDensGridPoints
                    do b=1,n+1
                      c=c+1
                      CFDMEkl_s(c)=CFDMEkl_s(c)+factor*Glob_YCoeff(k)*Glob_YCoeff(kk)*Denskl(b,a)
                    enddo
                  enddo
                endif

              enddo !kk=1,Glob_NumYTerms
            enddo !k=1,Glob_NumYTerms
          endif !SymmAdaptMethod==2

        endif
      enddo
    enddo

!Combining the results of all processes
    do a=1,NumOfExpcVals
      temp1=MEkl_s(a)
      call MPI_ALLREDUCE(temp1,temp2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      MEkl_s(a)=temp2
    enddo
    k=0
    if (AreCorrFuncNeeded) then
      k=NumCFGridPoints*n*(n+1)/2
      call MPI_ALLREDUCE(CFDMEkl_s,CF,k,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
    if (ArePartDensNeeded) then
      kk=NumDensGridPoints*(n+1)
      call MPI_ALLREDUCE(CFDMEkl_s(k+1:k+kk),Dens,kk,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

    if (spinDependentValuesNeeded == 1) then
      do c = 1, numberOfSpinFunctions
        do a = 1, n
          do b = a + 1, n
            temp1 = drach_SSFMatrix(a, b, c)
            call MPI_ALLREDUCE(temp1,temp2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
            drach_SSFMatrix(a, b, c) = temp2
          enddo
        enddo

        do a = 1, n
          do b = a + 1, n
            temp1 = SSFMatrix(a, b, c)
            call MPI_ALLREDUCE(temp1,temp2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
            SSFMatrix(a, b, c) = temp2
          enddo
        enddo

        if (positronPosition > 0) then
          do a = 1, n
            do b = a + 1, n
              temp1 = drach_AnihMatrix(a, b, c)
              call MPI_ALLREDUCE(temp1,temp2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
              drach_AnihMatrix(a, b, c) = temp2
            enddo
          enddo

          do a = 1, n
            do b = a + 1, n
              temp1 = AnihMatrix(a, b, c)
              call MPI_ALLREDUCE(temp1,temp2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
              AnihMatrix(a, b, c) = temp2
            enddo
          enddo
        endif

      enddo
    endif

!Extracting expectation values from arrays MEkl_s and MEkl_e
    c=0
    do a=1,n
      do b=a,n
        c=c+1
        rm2(b,a)=MEkl_s(c); rm2(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        rm(b,a)=MEkl_s(c); rm(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        r(b,a)=MEkl_s(c); r(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        r2(b,a)=MEkl_s(c); r2(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        deltar(b,a)=MEkl_s(c); deltar(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        drach_deltar(b,a)=MEkl_s(c); drach_deltar(a,b)=MEkl_s(c)
      enddo
    enddo
    do a=1,n
      do b=a,n
        c=c+1
        prval(b,a)=MEkl_s(c); prval(a,b)=MEkl_s(c)
      enddo
    enddo
    c=c+1; H=MEkl_s(c)
    c=c+1; S=MEkl_s(c)
    c=c+1; T=MEkl_s(c)
    c=c+1; V=MEkl_s(c)
    c=c+1; MV=MEkl_s(c)
    c=c+1; drach_MV1=MEkl_s(c)
    c=c+1; drach_MV2=MEkl_s(c)
    c=c+1; Darwin=MEkl_s(c)
    c=c+1; drach_Darwin=MEkl_s(c)
    c=c+1; OO=MEkl_s(c)
    rmrmkl(1:n,1:n,1:n,1:n)=ZERO
    do a=1,n
      do b=a,n
        do a1=a,n
          do b1=a1,n
            c=c+1; temp1=MEkl_s(c)
            rmrmkl(a,b,a1,b1)=temp1
            rmrmkl(a,b,b1,a1)=temp1
            rmrmkl(b,a,a1,b1)=temp1
            rmrmkl(b,a,b1,a1)=temp1
            rmrmkl(a1,b1,a,b)=temp1
            rmrmkl(a1,b1,b,a)=temp1
            rmrmkl(b1,a1,a,b)=temp1
            rmrmkl(b1,a1,b,a)=temp1
          enddo
        enddo
      enddo
    enddo

! spin-dependent expectation values calculation involves sum over all electronic permutations
    if (spinDependentValuesNeeded == 1) then

      ! we already have everything needed for SSF term calculation
      allocate(drach_SSF(numberOfSpinFunctions))
      drach_SSF = ZERO
      do k = 1, numberOfSpinFunctions
        do i = 1, n
          do j = i + 1, n
            drach_SSF(k) = drach_SSF(k) + drach_SSFMatrix(i, j, k) * SSFmassChargeCoefficient(i, j)
          enddo
        enddo
      enddo

      allocate(SSF(numberOfSpinFunctions))
      SSF = ZERO
      do k = 1, numberOfSpinFunctions
        do i = 1, n
          do j = i + 1, n
            SSF(k) = SSF(k) + SSFMatrix(i, j, k) * SSFmassChargeCoefficient(i, j)
          enddo
        enddo
      enddo

      ! we also calculate "electronic" SSF term as a test
      if (positronPosition > 0) then
        allocate(drach_SSFe(numberOfSpinFunctions))
        drach_SSFe = ZERO
        do k = 1, numberOfSpinFunctions
          do i = 1, n
            if (i == positronPosition) cycle
            do j = i + 1, n
              if (j == positronPosition) cycle
              drach_SSFe(k) = drach_SSFe(k) + drach_SSFMatrix(i, j, k) * SSFmassChargeCoefficient(i, j)
            enddo
          enddo
        enddo

        allocate(SSFe(numberOfSpinFunctions))
        SSFe = ZERO
        do k = 1, numberOfSpinFunctions
          do i = 1, n
            if (i == positronPosition) cycle
            do j = i + 1, n
              if (j == positronPosition) cycle
              SSFe(k) = SSFe(k) + SSFMatrix(i, j, k) * SSFmassChargeCoefficient(i, j)
            enddo
          enddo
        enddo

        allocate(drach_Anih(numberOfSpinFunctions))
        drach_Anih = ZERO
        do k = 1, numberOfSpinFunctions
          do i = 1, n
            do j = i + 1, n
              if (i /= positronPosition .and. j /= positronPosition) cycle
              drach_Anih(k) = drach_Anih(k) + drach_AnihMatrix(i, j, k) * AnihMassChargeCoefficient(i, j)
            enddo
          enddo
        enddo

        allocate(Anih(numberOfSpinFunctions))
        Anih = ZERO
        do k = 1, numberOfSpinFunctions
          do i = 1, n
            do j = i + 1, n
              if (i /= positronPosition .and. j /= positronPosition) cycle
              Anih(k) = Anih(k) + AnihMatrix(i, j, k) * AnihMassChargeCoefficient(i, j)
            enddo
          enddo
        enddo

      endif

      allocate(SO1(numberOfSpinFunctions), SO2(numberOfSpinFunctions), SSNC(numberOfSpinFunctions), &
               AMM1(numberOfSpinFunctions), AMM2(numberOfSpinFunctions), &
               AMM1Fin(numberOfSpinFunctions), AMM2Fin(numberOfSpinFunctions),&
               SO1kl(numberOfSpinFunctions), SO2kl(numberOfSpinFunctions), SSNCkl(numberOfSpinFunctions), &
               AMM1kl(numberOfSpinFunctions), AMM2kl(numberOfSpinFunctions),&
               AMM1Finkl(numberOfSpinFunctions), AMM2Finkl(numberOfSpinFunctions))

      SSNC = ZERO

      SO1 = ZERO
      SO2 = ZERO

      AMM1 = ZERO
      AMM2 = ZERO

      AMM1Fin = ZERO
      AMM2Fin = ZERO

      counter = 0
      do i = 1, cbs
        do j = 1, i
          counter = counter + 1
          if (mod(counter, Glob_NumOfProcs) == Glob_ProcID) then
            if (i == j) then
              factor = Glob_c(i) * Glob_c(j) / sqrt(diagS(i)*diagS(i))  !1
            else
              factor = TWO * Glob_c(i) * Glob_c(j) / sqrt(diagS(i)*diagS(j))  !2
            endif

            do a = 1, nFactorial ! Permutations from S_n introduced by A operator
              if (positronPosition > 0) then
                if (nint(ketYMatrix(positronPosition, positronPosition, a)) == 0) cycle
              endif ! we only need electronic permutations

              call spinDependentMatrixElements(Glob_Index(i,1),Glob_Index(j,1),Glob_Index(i,2), Glob_Index(j,2),     &
                                               Glob_NonlinParam(1 : npt, i), Glob_NonlinParam(1 : npt, j), &
                                               ketYMatrix(1 : n, 1 : n, a), &
                                     SziME(1 : n, 1 : numberOfSpinFunctions, a), SSNCspinME(1:n, 1:n, 1:numberOfSpinFunctions, a), &
                                               SSNCmassChargeCoefficient, SOmassChargeCoefficient, AMMmassChargeCoefficient, &
                                               AMMFinmassChargeCoefficient,&
                                               SSNCkl, SO1kl, SO2kl, AMM1kl, AMM2kl,  AMM1Finkl, AMM2Finkl, numberOfSpinFunctions)

              SSNC = SSNC + parityFactor(a) * factor * SSNCkl

              SO1 = SO1 + parityFactor(a) * factor * SO1kl
              SO2 = SO2 + parityFactor(a) * factor * SO2kl
              ! final value of SO and SSNC should be multiplied by the appropriate angular factors C_J
              ! (they are shown in spinData.txt)
              AMM1 = AMM1 + parityFactor(a) * factor * AMM1kl
              AMM2 = AMM2 + parityFactor(a) * factor * AMM2kl

              AMM1Fin = AMM1Fin + parityFactor(a) * factor * AMM1Finkl
              AMM2Fin = AMM2Fin + parityFactor(a) * factor * AMM2Finkl

            enddo ! Permutations from S_n

          endif ! ProcID check
        enddo
      enddo

      !Combining the results of all processes

      do k = 1, numberOfSpinFunctions

        temp1 = SO1(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        SO1(k) = temp2

        temp1 = SO2(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        SO2(k) = temp2

        temp1 = AMM1(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        AMM1(k) = temp2

        temp1 = AMM2(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        AMM2(k) = temp2

        temp1 = AMM1Fin(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        AMM1Fin(k) = temp2

        temp1 = AMM2Fin(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        AMM2Fin(k) = temp2

        temp1 = SSNC(k)
        call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        SSNC(k) = temp2

      enddo

    endif !spinDependentValuesNeeded

!Printing results
    if (Glob_ProcID==0) then

      !Opening an additional file where selected expectation values will be saved
      open(2,file=Glob_ExpValFileName,status='replace')

      write(*,*) 'done'
      write(*,*)
      write(*,*) 'Expectation values:'
      write(*,*)
      write(*,*) '                      H=',H
      write(*,*) '                      S=',S
      write(*,*) '                      T=',T
      write(*,*) '                      V=',V
      write(*,*) '                     MV=',MV
      if (Glob_ArePseudoParticleMassesTheSame) then
        write(*,*) '               drach_MV=',drach_MV1
      else
        write(*,*) '             drach_MV_l=',drach_MV1
        write(*,*) '             drach_MV_h=',drach_MV2
      endif
      write(*,*) '                 Darwin=',Darwin
      write(*,*) '           drach_Darwin=',drach_Darwin
      write(*,*) '                     OO=',OO

      if (spinDependentValuesNeeded == 1) then
        if (numberOfSpinFunctions == 1) then
          write(*,*) '                    SO1=',SO1(1)
          write(*,*) '                    SO2=',SO2(1)
          write(*,*) '                   AMM1=',AMM1(1)
          write(*,*) '                   AMM2=',AMM2(1)
          write(*,*) '                AMM1Fin=',AMM1Fin(1)
          write(*,*) '                AMM2Fin=',AMM2Fin(1)
          write(*,*) '              drach_SSF=',drach_SSF(1)
          write(*,*) '                    SSF=',SSF(1)
          write(*,*) '                    SSNC=',SSNC(1)
          if (positronPosition > 0) then
            write(*,*) '             drach_SSFe=',drach_SSFe(1)
            write(*,*) '             drach_Anih=',drach_Anih(1)
            write(*,*) '                   SSFe=',SSFe(1)
            write(*,*) '                   Anih=',Anih(1)
          endif
        elseif (numberOfSpinFunctions == 2) then
          write(*,*) '                  SO1_h=',SO1(1)
          write(*,*) '                  SO2_h=',SO2(1)
          write(*,*) '                 AMM1_h=',AMM1(1)
          write(*,*) '                 AMM2_h=',AMM2(1)
          write(*,*) '              AMM1Fin_h=',AMM1Fin(1)
          write(*,*) '              AMM2Fin_h=',AMM2Fin(1)
          write(*,*) '            drach_SSF_h=',drach_SSF(1)
          write(*,*) '           drach_SSFe_h=',drach_SSFe(1)
          write(*,*) '           drach_Anih_h=',drach_Anih(1)
          write(*,*) '                  SSF_h=',SSF(1)
          write(*,*) '                 SSFe_h=',SSFe(1)
          write(*,*) '                 Anih_h=',Anih(1)
          write(*,*) '                  SSNC_h=',SSNC(1)

          write(*,*) '                  SO1_l=',SO1(2)
          write(*,*) '                  SO2_l=',SO2(2)
          write(*,*) '                 AMM1_l=',AMM1(2)
          write(*,*) '                 AMM2_l=',AMM2(2)
          write(*,*) '              AMM1Fin_l=',AMM1Fin(2)
          write(*,*) '              AMM2Fin_l=',AMM2Fin(2)
          write(*,*) '            drach_SSF_l=',drach_SSF(2)
          write(*,*) '           drach_SSFe_l=',drach_SSFe(2)
          write(*,*) '           drach_Anih_l=',drach_Anih(2)
          write(*,*) '                  SSF_l=',SSF(2)
          write(*,*) '                 SSFe_l=',SSFe(2)
          write(*,*) '                 Anih_l=',Anih(2)
          write(*,*) '                  SSNC_l=',SSNC(2)

        endif
      endif
      write(*,*)

      write(*,*) '           (alpha^2)*MV=',MV*(Glob_FineStructConst**2)
      if (Glob_ArePseudoParticleMassesTheSame) then
        write(*,*) '     (alpha^2)*drach_MV=',drach_MV1*(Glob_FineStructConst**2)
      else
        write(*,*) '   (alpha^2)*drach_MV_l=',drach_MV1*(Glob_FineStructConst**2)
        write(*,*) '   (alpha^2)*drach_MV_h=',drach_MV2*(Glob_FineStructConst**2)
      endif
      write(*,*) '       (alpha^2)*Darwin=',Darwin*(Glob_FineStructConst**2)
      write(*,*) ' (alpha^2)*drach_Darwin=',drach_Darwin*(Glob_FineStructConst**2)
      write(*,*) '           (alpha^2)*OO=',OO*(Glob_FineStructConst**2)

      if (spinDependentValuesNeeded == 1) then
        if (numberOfSpinFunctions == 1) then
          write(*,*) '        (alpha^2)*SO1=',SO1(1)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SO2=',SO2(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM1=',AMM1(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM2=',AMM2(1)*(Glob_FineStructConst**2)
          write(*,*) '  (alpha^2)*drach_SSF=',drach_SSF(1)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SSF=',SSF(1)*(Glob_FineStructConst**2)
          write(*,*) '(alpha^2)*SSNC = ', SSNC(1)*(Glob_FineStructConst**2)
          if (positronPosition > 0) then
            write(*,*) ' (alpha^2)*drach_SSFe=',drach_SSFe(1)*(Glob_FineStructConst**2)
            write(*,*) ' (alpha^2)*drach_Anih=',drach_Anih(1)*(Glob_FineStructConst**2)
            write(*,*) '       (alpha^2)*SSFe=',SSFe(1)*(Glob_FineStructConst**2)
            write(*,*) '       (alpha^2)*Anih=',Anih(1)*(Glob_FineStructConst**2)
          endif
        elseif (numberOfSpinFunctions == 2) then
          write(*,*) '        (alpha^2)*SO1_h=',SO1(1)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SO2_h=',SO2(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM1_h=',AMM1(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM2_h=',AMM2(1)*(Glob_FineStructConst**2)
          write(*,*) '  (alpha^2)*drach_SSF_h=',drach_SSF(1)*(Glob_FineStructConst**2)
          write(*,*) ' (alpha^2)*drach_SSFe_h=',drach_SSFe(1)*(Glob_FineStructConst**2)
          write(*,*) ' (alpha^2)*drach_Anih_h=',drach_Anih(1)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SSF_h=',SSF(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*SSFe_h=',SSFe(1)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*Anih_h=',Anih(1)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SSNC_h=',SSNC(1)*(Glob_FineStructConst**2)

          write(*,*) '        (alpha^2)*SO1_l=',SO1(2)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SO2_l=',SO2(2)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM1_l=',AMM1(2)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*AMM2_l=',AMM2(2)*(Glob_FineStructConst**2)
          write(*,*) '  (alpha^2)*drach_SSF_l=',drach_SSF(2)*(Glob_FineStructConst**2)
          write(*,*) ' (alpha^2)*drach_SSFe_l=',drach_SSFe(2)*(Glob_FineStructConst**2)
          write(*,*) ' (alpha^2)*drach_Anih_l=',drach_Anih(2)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SSF_h=',SSF(2)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*SSFe_h=',SSFe(2)*(Glob_FineStructConst**2)
          write(*,*) '       (alpha^2)*Anih_h=',Anih(2)*(Glob_FineStructConst**2)
          write(*,*) '        (alpha^2)*SSNC_l=',SSNC(2)*(Glob_FineStructConst**2)

        endif
      endif

      write(*,*)

      write(2,'(a)',advance='no') '                  basis '
      write(2,*) cbs
      write(2,'(a)',advance='no') '                 Energy '
      call writerealadv(2,Evalue)
      write(2,'(a)',advance='no') '                      H '
      call writerealadv(2,H)
      write(2,'(a)',advance='no') '                      S '
      call writerealadv(2,S)
      write(2,'(a)',advance='no') '                      T '
      call writerealadv(2,T)
      write(2,'(a)',advance='no') '                      V '
      call writerealadv(2,V)
      write(2,'(a)',advance='no') '                     MV '
      call writerealadv(2,MV)
      if (Glob_ArePseudoParticleMassesTheSame) then
        write(2,'(a)',advance='no') '               drach_MV '
        call writerealadv(2,drach_MV1)
      else
        write(2,'(a)',advance='no') '             drach_MV_l '
        call writerealadv(2,drach_MV1)
        write(2,'(a)',advance='no') '             drach_MV_h '
        call writerealadv(2,drach_MV2)
      endif
      write(2,'(a)',advance='no') '                 Darwin '
      call writerealadv(2,Darwin)
      write(2,'(a)',advance='no') '           drach_Darwin '
      call writerealadv(2,drach_Darwin)
      write(2,'(a)',advance='no') '                     OO '
      call writerealadv(2,OO)
      write(2,'(a)',advance='no') '           (alpha^2)*MV '
      call writerealadv(2,MV*(Glob_FineStructConst**2))
      if (Glob_ArePseudoParticleMassesTheSame) then
        write(2,'(a)',advance='no') '     (alpha^2)*drach_MV '
        call writerealadv(2,drach_MV1*(Glob_FineStructConst**2))
      else
        write(2,'(a)',advance='no') '   (alpha^2)*drach_MV_l '
        call writerealadv(2,drach_MV1*(Glob_FineStructConst**2))
        write(2,'(a)',advance='no') '   (alpha^2)*drach_MV_h '
        call writerealadv(2,drach_MV2*(Glob_FineStructConst**2))
      endif
      write(2,'(a)',advance='no') '       (alpha^2)*Darwin '
      call writerealadv(2,Darwin*(Glob_FineStructConst**2))
      write(2,'(a)',advance='no') ' (alpha^2)*drach_Darwin '
      call writerealadv(2,drach_Darwin*(Glob_FineStructConst**2))
      write(2,'(a)',advance='no') '           (alpha^2)*OO '
      call writerealadv(2,OO*(Glob_FineStructConst**2))

      if ((Glob_NumOfIdentPartSets/=Glob_n+1).and.(SymmAdaptMethod==1)) then
        write(*,*) '(Warning! These values do not account for indistinguishability of'
        write(*,*) 'identical particles and other possible symmetries of the system)'
        write(*,*)
      endif
      do i=1,n
        write(*,'(1x,a22,i1)',advance='no') '                1/r^2_',i
        write(*,*) '=',rm2(i,i)
        do j=i+1,n
          write(*,'(1x,a21,i1,i1)',advance='no') '               1/r^2_',i,j
          write(*,*) '=',rm2(i,j)
        enddo
      enddo
      do i=1,n
        write(*,'(1x,a22,i1)',advance='no') '                  1/r_',i
        write(*,*) '=',rm(i,i)
        do j=i+1,n
          write(*,'(1x,a21,i1,i1)',advance='no') '                 1/r_',i,j
          write(*,*)'=',rm(i,j)
        enddo
      enddo
      write(*,*)
      do i=1,n
        write(*,'(1x,a22,i1)',advance='no') '                    r_',i
        write(*,*) '=',r(i,i)
        do j=i+1,n
          write(*,'(1x,a21,i1,i1)',advance='no') '                   r_',i,j
          write(*,*) '=',r(i,j)
        enddo
      enddo
      write(*,*)
      do i=1,n
        write(*,'(1x,a22,i1)',advance='no') '                  r^2_',i
        write(*,*) '=',r2(i,i)
        do j=i+1,n
          write(*,'(1x,a21,i1,i1)',advance='no') '                 r^2_',i,j
          write(*,*) '=',r2(i,j)
        enddo
      enddo
      write(*,*)
      do i=1,n
        write(*,'(1x,a21,i1,a1)',advance='no') '            delta(r_',i,')'
        write(*,*) '=',deltar(i,i)
        do j=i+1,n
          write(*,'(1x,a20,i1,i1,a1)',advance='no') '            delta(r_',i,j,')'
          write(*,*) '=',deltar(i,j)
        enddo
      enddo
      write(*,*)
      do i=1,n
        write(*,'(1x,a21,i1,a1)',advance='no') '       drach_delta(r_',i,')'
        write(*,*) '=',drach_deltar(i,i)
        do j=i+1,n
          write(*,'(1x,a20,i1,i1,a1)',advance='no') '      drach_delta(r_',i,j,')'
          write(*,*) '=',drach_deltar(i,j)
        enddo
      enddo
      do i=1,n
        write(*,'(1x,a21,i1,a1)',advance='no') '             prval(r_',i,')'
        write(*,*) '=',prval(i,i)
        do j=i+1,n
          write(*,'(1x,a20,i1,i1,a1)',advance='no') '            prval(r_',i,j,')'
          write(*,*) '=',prval(i,j)
        enddo
      enddo
      write(*,*)
      do i=1,n
        do j=i,n
          do a=i,n
            do b=a,n
              write(*,'(4x,a,i1)',advance='no') '1/(r_',i
              if (i/=j) write(*,'(i1)',advance='no') j
              write(*,'(a)',advance='no') '*'
              write(*,'(a,i1)',advance='no') 'r_',a
              if (a/=b) write(*,'(i1)',advance='no') b
              write(*,'(a)',advance='no') ')'
              write(*,*) '=',rmrmkl(i,j,a,b)
            enddo
          enddo
        enddo
      enddo
      write(*,*)
      if (Glob_NumOfIdentPartSets/=Glob_n+1) then
        write(*,*) 'Based on the particle mass and charge values it was determined'
        write(*,*) 'that the system has the following sets of identical particles:'
        do i=1,Glob_NumOfIdentPartSets
          j=Glob_NumOfPartInIdentPartSet(i)
          write(*,'(1x,a3,i2,a13)',advance='no') 'set',i,' :  particles'
          write(*,*) Glob_IdentPartList(1:j,i)
        enddo
        write(*,*)
        write(*,*) 'Properly symmetrized expectation values of two-particle quantities'
        write(*,*) 'that account for permutational symmetry of the above mentioned sets'
        write(*,*) 'of identical particles are:'
        write(*,*) '(Warning! An additional symmetrization might be necessary if the'
        write(*,*) 'Young operator contains other types of symmetries)'
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a6,i1,a3)',advance='no') '1/r^2_',a,' = '
            else
              write(*,'(a6,i1,i1,a3)',advance='no') '1/r^2_',a,b,' = '
            endif
            beta=beta+rm2(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,1x)',advance='no') '               1/r^2_',a,b
          if (a==b) write(2,'(a,i1,1x)',advance='no')    '                1/r^2_',a
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a4,i1,a3)',advance='no') '1/r_',a,' = '
            else
              write(*,'(a4,i1,i1,a3)',advance='no') '1/r_',a,b,' = '
            endif
            beta=beta+rm(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,1x)',advance='no') '                 1/r_',a,b
          if (a==b) write(2,'(a,i1,1x)',advance='no')    '                  1/r_',a
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a2,i1,a3)',advance='no') 'r_',a,' = '
            else
              write(*,'(a2,i1,i1,a3)',advance='no') 'r_',a,b,' = '
            endif
            beta=beta+r(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,1x)',advance='no') '                   r_',a,b
          if (a==b) write(2,'(a,i1,1x)',advance='no')    '                    r_',a
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a4,i1,a3)',advance='no') 'r^2_',a,' = '
            else
              write(*,'(a4,i1,i1,a3)',advance='no') 'r^2_',a,b,' = '
            endif
            beta=beta+r2(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,1x)',advance='no') '                 r^2_',a,b
          if (a==b) write(2,'(a,i1,1x)',advance='no')    '                  r^2_',a
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a8,i1,a4)',advance='no') 'delta(r_',a,') = '
            else
              write(*,'(a8,i1,i1,a4)',advance='no') 'delta(r_',a,b,') = '
            endif
            beta=beta+deltar(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,a1,1x)',advance='no') '            delta(r_',a,b,')'
          if (a==b) write(2,'(a,i1,a1,1x)',advance='no')    '             delta(r_',a,')'
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a14,i1,a4)',advance='no') 'drach_delta(r_',a,') = '
            else
              write(*,'(a14,i1,i1,a4)',advance='no') 'drach_delta(r_',a,b,') = '
            endif
            beta=beta+drach_deltar(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,a1,1x)',advance='no') '      drach_delta(r_',a,b,')'
          if (a==b) write(2,'(a,i1,a1,1x)',advance='no')    '       drach_delta(r_',a,')'
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
        do i=1,Glob_NumOfNoneqvPairSets
          beta=ZERO
          mu=ZERO
          k=Glob_NumOfPairsInEqvPairSet(i)
          write(*,'(1x)',advance='no')
          do j=1,k
            a=Glob_EqvPairList(1,j,i)
            b=Glob_EqvPairList(2,j,i)
            if (a==b) then
              write(*,'(a8,i1,a4)',advance='no') 'prval(r_',a,') = '
            else
              write(*,'(a8,i1,i1,a4)',advance='no') 'prval(r_',a,b,') = '
            endif
            beta=beta+prval(a,b)
          enddo
          call writerealadv(6,beta/k)
          !write to file
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a/=b) write(2,'(a,i1,i1,a1,1x)',advance='no') '            prval(r_',a,b,')'
          if (a==b) write(2,'(a,i1,a1,1x)',advance='no')    '             prval(r_',a,')'
          call writerealadv(2,beta/k)
        enddo
        write(*,*)
      endif

      close(2)

      !Saving correlation functions
      if (AreCorrFuncNeeded) then
        open(1,file=FileName2,status='replace')
        !first we print titles of all data columns
        write(1,'(6x,a7,21x,a4)',advance='no') '#xi_rho','xi_z'
        do i=1,Glob_NumOfNoneqvPairSets
          a=Glob_EqvPairList(1,1,i)
          b=Glob_EqvPairList(2,1,i)
          if (a>b) then
            c=a; a=b; b=c
          endif
          if (a==b) then
            write(1,'(21x,a1,i1,1x)',advance='no') 'g',a
          else
            write(1,'(21x,a1,i1,i1)',advance='no') 'g',a,b
          endif
        enddo
        write(1,*)
        !then we print the data columns themselves
        do kk=1,NumCFGridPoints
          write(1,'(2(1x,e23.16))',advance='no') CFGrid(1,kk),CFGrid(2,kk)
          do i=1,Glob_NumOfNoneqvPairSets
            mu=ZERO
            k=Glob_NumOfPairsInEqvPairSet(i)
            do j=1,k
              a=Glob_EqvPairList(1,j,i)
              b=Glob_EqvPairList(2,j,i)
              if (a>b) then
                c=a; a=b; b=c
              endif
              mu=mu+CF(b-(a-1)*(a-2*n)/2,kk)
            enddo
            write(1,'(1x,e23.16)',advance='no') mu/k
          enddo
          write(1,*)
        enddo
        close(1)
        i=len_trim(FileName2)
        write(*,*) 'Correlation functions have been stored in file',FileName2(1:i)
        write(*,*)
      endif

      !Saving particle densities
      if (ArePartDensNeeded) then
        open(1,file=FileName4,status='replace')
        !first we print titles of all data columns
        write(1,'(6x,a7,21x,a4)',advance='no') '#xi_rho','xi_z'
        do i=1,Glob_NumOfIdentPartSets
          write(1,'(20x,a3,i1)',advance='no') 'rho',Glob_IdentPartList(1,i)
        enddo
        write(1,*)
        !then we print the data columns themselves
        do kk=1,NumDensGridPoints
          write(1,'(2(1x,e23.16))',advance='no') DensGrid(1,kk),DensGrid(2,kk)
          do i=1,Glob_NumOfIdentPartSets
            mu=ZERO
            k=Glob_NumOfPartInIdentPartSet(i)
            do j=1,k
              mu=mu+Dens(Glob_IdentPartList(j,i),kk)
            enddo
            write(1,'(1x,e23.16)',advance='no') mu/k
          enddo
          write(1,*)
        enddo
        close(1)
        i=len_trim(FileName4)
        write(*,*) 'Particle densities have been stored in file',FileName4(1:i)
        write(*,*)
      endif

    endif

!deallocate local arrays

    deallocate(rm2kl)
    deallocate(rm2)

    deallocate(rm2kl1)
    deallocate(rm2kl2)

    deallocate(rmkl)
    deallocate(rm)

    deallocate(rmkl1)
    deallocate(rmkl2)

    deallocate(rkl)
    deallocate(r)

    deallocate(rkl1)
    deallocate(rkl2)

    deallocate(r2kl)
    deallocate(r2)

    deallocate(r2kl1)
    deallocate(r2kl2)

    deallocate(deltarkl)
    deallocate(deltar)

    deallocate(deltarkl1)
    deallocate(deltarkl2)

    deallocate(drach_deltarkl)
    deallocate(drach_deltar)
    deallocate(drach_deltarkl1)
    deallocate(drach_deltarkl2)

    deallocate(prvalkl)
    deallocate(prval)

    deallocate(prvalkl1)
    deallocate(prvalkl2)

    deallocate(MEkl)
    deallocate(MEkl_s)

    deallocate(MEkl1)
    deallocate(MEkl_s1)
    deallocate(MEkl2)
    deallocate(MEkl_s2)

    deallocate(IdentityPerm)

    if (GSEPSolMethod=='G') then
      deallocate(IFAIL)
      deallocate(Eigvecs)
      deallocate(Eigvals)
    endif

    if (GSEPsolMethod=='I') then
      deallocate(Glob_LastEigvector)
      deallocate(Glob_WorkForGSEPIIS)
    endif
    if (GSEPSolMethod=='Q') call ClearQWorkspace()

!deallocate workspace for DSYGVX
    if (GSEPSolMethod=='G') then
      deallocate(Glob_WorkForDSYGVX)
      deallocate(Glob_IWorkForDSYGVX)
    endif

!deallocate global arrays
    deallocate(Glob_SklBuff2)
    deallocate(Glob_SklBuff1)
    deallocate(Glob_HklBuff2)
    deallocate(Glob_HklBuff1)
    deallocate(Glob_c)
    if (GSEPSolMethod=='I') deallocate(Glob_invD)
    deallocate(Glob_diagS)
    if (GSEPSolMethod=='G') deallocate(Glob_diagH)
    deallocate(Glob_S)
    deallocate(Glob_H)

    if((AreCorrFuncNeeded).or.(ArePartDensNeeded)) then
      deallocate(CFDMEkl_s)
    endif

    deallocate(DensGrid)
    deallocate(Denskl)
    if (ArePartDensNeeded) deallocate(Dens)
    deallocate(CFGrid)
    deallocate(CFkl)
    if (AreCorrFuncNeeded) deallocate(CF)

    if (Glob_ProcID==0) write (*,*) 'Routine ExpectationValues finished'

  end subroutine ExpectationValues

end module workproc
