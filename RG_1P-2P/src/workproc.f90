module workproc
!This module contains basic work SUBROUTINEs
  USE matelem
  IMPLICIT NONE

CONTAINS

  SUBROUTINE READwf1wf2()
! SUBROUTINE READwf1wf2 READs data (number of particle,
! mass, charge, nonlinear variational
! parameters and other information) from the wave function
! FILE whose name is specIFied by global variable
! Glob_wf FILE FILEs.

!Local variables:
    INTEGER        :: OPENFILEErr
    INTEGER        :: READInt,READErr
    INTEGER        :: particle_n1,particle_n2
    REAL(wp)    :: Mass1(Glob_AllowedNumOfParticles),Mass2(Glob_AllowedNumOfParticles)
    REAL(wp)    :: PseudoCharge1(Glob_AllowedNumOfParticles),PseudoCharge2(Glob_AllowedNumOfParticles)
    REAL(wp)    :: READReal1,READReal2
    INTEGER        :: WorkInt(MAX(MAX(Glob_YOperatorStringLength,20),Glob_FILENameLength))
    INTEGER        :: WorkInt1(MAX(MAX(Glob_YOperatorStringLength,20),Glob_FILENameLength))
    INTEGER        :: WorkInt2(MAX(MAX(Glob_YOperatorStringLength,20),Glob_FILENameLength))
    INTEGER        :: i,j,k,l,Line1,Line2
    CHARACTER(70)  :: readchar
    logical        :: ErrorInDataFILE !,IsDRMCStep
!==============================================================================================

    ErrorInDataFILE=.false.

!-------------------------------------------------------------------
! opening the wave function FILEs of initial and final states
!-------------------------------------------------------------------
! Glob_WFFILE1 initial state's FILE
! Glob_WFFILE2 final state's FILE

    IF (Glob_ProcID==0) THEN

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(4x,A)') 'Reading Wave functions       '
      WRITE(*,'(A)') '================================================='
      WRITE(*,*)

!   initial state
      OPEN(1,FILE=Glob_WFFILE1,STATUS='old',IOSTAT=OPENFILEErr)

      IF (OPENFILEErr/=0) THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0210 in opening wave function of inital state, '
        WRITE(*,*) '   "',trim(adjustl(Glob_WFFILE1)), '"  FILE not found !!!'
        WRITE(*,*)
        ErrorInDataFILE=.true.
      ENDIF

!   final state
      OPEN(2,FILE=Glob_WFFILE2,STATUS='old',IOSTAT=OPENFILEErr)

      IF (OPENFILEErr/=0) THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0211 in opening wave function of final state, '
        WRITE(*,*) ' "',trim(adjustl(Glob_WFFILE2)), '" FILE not found !!!'
        WRITE(*,*)
        ErrorInDataFILE=.true.
      ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInDataFILE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFILE) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

!-------------------------------------------------------------------
! reading Header
!-------------------------------------------------------------------

    IF (Glob_ProcID==0) Line1=0
    IF (Glob_ProcID==0) Line2=0

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:30)
      WRITE(*,*) '  Reading " ',readchar(1:5) ,'" of the initial state from :    ', trim(adjustl(Glob_WFFILE1))
      Line1=Line1+1

!   final state
      READ(2,*) readchar(1:30)
      WRITE(*,*) '  Reading " ',readchar(1:5) ,'" of the final   state from :    ', trim(adjustl(Glob_WFFILE2))
      Line2=Line2+1
      WRITE (*,*)
      WRITE(*,*)

    ENDIF

!-------------------------------------------------------------------
! reading the number of particle from the wave function FILEs
!-------------------------------------------------------------------

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:9),READInt
      Line1=Line1+1

!   particle_n1 is the number of pseuDOparticles in Glob_WFFILE1
      particle_n1=READInt-1

      IF ((particle_n1<1).or.(readchar(1:9)/='PARTICLES')) THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0212 in reading the number of particles of the initial state, line: ',Line1
        WRITE(*,*)
        ErrorInDataFILE=.true.
      ENDIF

!   final state
      READ(2,*) readchar(1:9),READInt
      Line2=Line2+1
!   particle_n2 is the number of pseuDOparticles in Glob_WFFILE2
      particle_n2=READInt-1

      IF ((particle_n2<1).or.(readchar(1:9)/='PARTICLES')) THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0213 in reading the number of particles of the final state, line: ',Line2
        WRITE(*,*)
        ErrorInDataFILE=.true.
      ENDIF

!   comparing the number of particle in Glob_WFFILE1 and Glob_WFFILE2 FILEs
      IF (particle_n1/=particle_n2) THEN

        WRITE(*,*)
        WRITE(*,*) '  The number of particles in initial and final states is not the same !!!'
        WRITE(*,*) '  The program will STOP now.'
        WRITE(*,*)
        ErrorInDataFILE=.true.

      ELSE

        IF (particle_n1/=Glob_AllowedNumOfPseudoParticles) then
          write(*,*) ' '
          write (*,*) 'The version of the code you are running was compiled for the case'
          write (*,*) 'when the number of particles in the system is equal to', &
            Glob_AllowedNumOfParticles
          write (*,*) 'while the number of particles specified in the wave function files is',particle_n1+1
          write (*,*) 'Please make appropriate changes. Program will now stop.'
          write(*,*) ' '
          ErrorInDataFile=.true.
        EndIF
        Glob_n=particle_n1

      ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInDataFILE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFILE) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    CALL MPI_BCAST(Glob_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    Glob_2Raised3n2=TWO**((3*Glob_n)/TWO)
    Glob_PiRaised3n2=Glob_Pi**((3*Glob_n)/TWO)
    Glob_np=Glob_n*(Glob_n+1)/2
    Glob_npt=Glob_np

!-------------------------------------------------------------------
! reading the masses of particles from the wave function FILEs
!-------------------------------------------------------------------
! Mass1 masses which are READ from the Glob_WFFILE1
! Mass2 masses which are READ from the Glob_WFFILE2

    ALLOCATE(Glob_Mass(Glob_n+1))

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:6),Mass1(1:Glob_n+1)
      Line1=Line1+1

      IF (readchar(1:6)/='MASSES') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0214 in reading masses of particles of initial state, line: ',Line1
        WRITE(*,*)
      ENDIF

!   final state
      READ(2,*) readchar(1:6),Mass2(1:Glob_n+1)
      Line2=Line2+1

      IF (readchar(1:6)/='MASSES') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0215 in reading masses of particles of final state, line: ',Line2
        WRITE(*,*)
      ENDIF

!   comparison of the masses of two FILEs
      DO i=1,Glob_n+1

        IF (Mass1(i)/=Mass2(i))THEN

          IF(Mass1(i) >= 1.00E20 .AND. Mass2(i) >= 1.00E20) THEN

            WRITE(*,*) '        ************************************************* '
            WRITE(*,*) '        *                   WARNING!!!                  * '
            WRITE(*,'(4X,A,I1,A)') '     *        The mass of the ',i,'-th particle          * '
            WRITE(*,*) '        *  in initial and final states is not the same. * '
            WRITE(*,*) '        * But, both of them are equal to or larger than * '
            WRITE(*,*) '        *                   1.00E20.                    * '
            WRITE(*,*) '        *       It is assumed that they are same.       * '
            WRITE(*,*) '        * The largest value is using in the calculation.* '
            WRITE(*,*) '        ************************************************* '
            WRITE(*,*)

          ELSE

            WRITE(*,*)
            WRITE(*,*) '  The mass of the',i,'-th particle in initial and final states is not the same !!!'
            WRITE(*,*)
            ErrorInDataFILE=.true.

          ENDIF

        ENDIF

      ENDDO

      Glob_Mass=Mass1
      ! WRITE(*,'(1x,a6)',advance='no') readchar(1:6)
      ! CALL writerealarradv(6,Glob_Mass,Glob_n+1)

    ENDIF

    CALL MPI_BCAST(ErrorInDataFILE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFILE) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    CALL MPI_BCAST(Glob_Mass,Glob_n+1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!-------------------------------------------------------------------
! reading the charges of particles from the wave function FILEs
!-------------------------------------------------------------------

! PseudoCharge1 charges of the nucleus and electrons which are READ from the Glob_WFFILE1
! PseudoCharge2 charges of the nucleus and electrons which are READ from the Glob_WFFILE2

    ALLOCATE(Glob_PseudoCharge(Glob_n))

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:7),PseudoCharge1(1:Glob_n+1)
      Line1=Line1+1

      IF (readchar(1:7)/='CHARGES') THEN

        WRITE(*,*)
        WRITE(*,*) '  Error EC0216 in reading charges of particles of initial state, line: ',Line1
        WRITE(*,*)

      ENDIF

!   final state
      READ(2,*) readchar(1:7),PseudoCharge2(1:Glob_n+1)
      Line2=Line2+1

      IF (readchar(1:7)/='CHARGES') THEN

        WRITE(*,*)
        WRITE(*,*) '  Error EC0217 in reading charges of particles of final state, line: ',Line2
        WRITE(*,*)

      ENDIF

!   comparing charges of the two FILEs

      DO i=1,Glob_n

        IF (PseudoCharge1(i)/=PseudoCharge2(i))THEN

          WRITE(*,*)
          WRITE(*,*) '  The charge of the',i+1,' in initial and final states is not the same !!! '
          WRITE(*,*)
          ErrorInDataFILE=.true.

        ENDIF

      ENDDO

      IF(ErrorInDataFILE .eqv. .FALSE.) THEN

        Glob_PseudoCharge0=PseudoCharge1(1)
        DO i=1,Glob_n
          Glob_PseudoCharge(i)=PseudoCharge1(i+1)
        ENDDO

      ENDIF

      ! WRITE(*,'(1x,a7)',advance='no') readchar(1:7)
      ! CALL writeREAL(6,Glob_PseudoCharge0)
      ! CALL writerealarradv(6,Glob_PseudoCharge,Glob_n)
      WRITE(6,'(A10,7X,A10,25X,A6)') 'PARTICLE','MASS','CHARGE'
      DO i=1,Glob_n+1
        WRITE(6,'(I5,7X)',advance='no') i
        CALL writereal(6,Glob_Mass(i))
        WRITE(6,'(7X)',advance='no')
        CALL writerealadv(6,PseudoCharge1(i))
      ENDDO

    ENDIF

    CALL MPI_BCAST(ErrorInDataFILE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFILE) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    CALL MPI_BCAST(Glob_PseudoCharge0,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_PseudoCharge,Glob_n,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!-------------------------------------------------------------------
! reading Young operators
!-------------------------------------------------------------------

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:9),Glob_YOperatorString1
      Line1=Line1+1

      IF (readchar(1:9)/='SYMMETRY') THEN

        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) readchar(1:9)
        WRITE(*,*) '  Error EC0218 in reading symmetry of initial state, line: ',Line1
        WRITE(*,*)

      ENDIF

      DO i=1,Glob_YOperatorStringLength
        WorkInt(i)=ICHAR(Glob_YOperatorString1(i:i))
      ENDDO

      WorkInt1=WorkInt

!   final state
      READ(2,*) readchar(1:9),Glob_YOperatorString2
      Line2=Line2+1

      IF (readchar(1:9)/='SYMMETRY') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0219 in reading symmetry of final state, line: ',Line2
        WRITE(*,*)
      ENDIF

      DO i=1,Glob_YOperatorStringLength
        WorkInt(i)=ICHAR(Glob_YOperatorString2(i:i))
      ENDDO

      WorkInt2=WorkInt

      IF (Glob_ProcID==0) then
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,'(2X,A,A)') 'The Young Operator for L=1:  ', TRIM(Glob_YOperatorString1)
        WRITE(*,'(2X,A,A)') 'The Young Operator for L=2:  ', TRIM(Glob_YOperatorString2)
      ENDIF

! !   comparing Young operators of the two FILEs

      ! IF(WorkInt1(i)/=WorkInt2(i)) THEN
      ! WRITE(*,*)
      ! WRITE(*,*) 'the Young operators in initial and final states are not the same !!!'
      ! WRITE(*,*)
      ! ErrorInDataFILE=.true.
      ! ENDIF

    ENDIF

    CALL MPI_BCAST(ErrorInDataFILE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFILE) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    CALL MPI_BCAST(WorkInt1,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(WorkInt2,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    DO i=1,Glob_YOperatorStringLength

      Glob_YOperatorString1(i:i)=char(WorkInt1(i))
      Glob_YOperatorString2(i:i)=char(WorkInt2(i))

    ENDDO

!-------------------------------------------------------------------
! reading number of the basis functions for each state
!-------------------------------------------------------------------

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:10),Glob_CurrBasisSize1
      Line1=Line1+1

      IF (readchar(1:10)/='BASIS_SIZE') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0220 in reading BASIS_SIZE of initial state, line: ',Line1
        WRITE(*,*)
      ENDIF

      WRITE(*,*)
      WRITE(*,'(2X,A,I6)') 'INITIAL STATE BASIS SIZE:  ', Glob_CurrBasisSize1

!    final state
      READ(2,*) readchar(1:10),Glob_CurrBasisSize2
      Line2=Line2+1

      IF (readchar(1:10)/='BASIS_SIZE') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0221 in reading BASIS_SIZE of final state, line: ',Line2
        WRITE(*,*)
      ENDIF

      WRITE(*,'(2X,A,I6)') 'FINAL   STATE BASIS SIZE:  ', Glob_CurrBasisSize2

    ENDIF

    CALL MPI_BCAST(Glob_CurrBasisSize1,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_CurrBasisSize2,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

!-------------------------------------------------------------------
! reading energy of each state
!-------------------------------------------------------------------

    IF (Glob_ProcID==0) THEN

!   initial state
      READ(1,*) readchar(1:14),Glob_CurrEnergy1
      Line1=Line1+1

      IF (readchar(1:14)/='CURRENT_ENERGY') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0222 in reading CURRENT_ENERGY of initial state, line: ',Line1
        WRITE(*,*)
      ENDIF

      WRITE(*,*)
      WRITE(*,'(2X,A)',advance='no')  'INITIAL STATE ENERGY:  '
      CALL WRITErealadv(6,Glob_CurrEnergy1)

!   final state
      READ(2,*) readchar(1:14),Glob_CurrEnergy2
      Line2=Line2+1

      IF (readchar(1:14)/='CURRENT_ENERGY') THEN
        WRITE(*,*)
        WRITE(*,*) '  Error EC0223 in reading CURRENT_ENERGY of final state, line: ',Line2
        WRITE(*,*)
      ENDIF

      WRITE(*,'(2X,A)',advance='no')  'FINAL   STATE ENERGY:  '
      CALL WRITErealadv(6,Glob_CurrEnergy2)

    ENDIF

    CALL MPI_BCAST(Glob_CurrEnergy1,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_CurrEnergy2,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    IF (Glob_ProcID==0) THEN
      READ(1,*) readchar
      READ(2,*) readchar
    ENDIF

! reading the parameters of the basis functions
    ALLOCATE(Glob_c1(Glob_CurrBasisSize1))
    ALLOCATE(Glob_FuncNum1(Glob_CurrBasisSize1))
    ALLOCATE(Glob_NonlinParam1(Glob_npt,Glob_CurrBasisSize1))
    ALLOCATE(Glob_Index1(Glob_CurrBasisSize1))

    ALLOCATE(Glob_c2(Glob_CurrBasisSize2))
    ALLOCATE(Glob_FuncNum2(Glob_CurrBasisSize2))
    ALLOCATE(Glob_NonlinParam2(Glob_npt,Glob_CurrBasisSize2))
    allocate(Glob_Index2(Glob_CurrBasisSize2,2))

    IF (Glob_ProcID==0) THEN

!   initial state
      DO i=1,Glob_CurrBasisSize1
        READ(1,*) Glob_FuncNum1(i),Glob_c1(i),readchar(1:2),Glob_Index1(i),Glob_NonlinParam1(1:Glob_npt,i)
      ENDDO

!   final state
      DO i=1,Glob_CurrBasisSize2
        READ(2,*) Glob_FuncNum2(i),Glob_c2(i),readchar(1:2),Glob_Index2(i,1),Glob_Index2(i,2),Glob_NonlinParam2(1:Glob_npt,i)
      ENDDO

    ENDIF

    CALL MPI_BCAST(Glob_c1,Glob_CurrBasisSize1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_c2,Glob_CurrBasisSize2,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    CALL MPI_BCAST(Glob_NonlinParam1,Glob_npt*Glob_CurrBasisSize1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_NonlinParam2,Glob_npt*Glob_CurrBasisSize2,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    CALL MPI_BCAST(Glob_FuncNum1,Glob_CurrBasisSize1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_FuncNum2,Glob_CurrBasisSize2,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    CALL MPI_BCAST(Glob_Index1,Glob_CurrBasisSize1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_Index2,2*Glob_CurrBasisSize2,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

  END SUBROUTINE READwf1wf2

  SUBROUTINE DataInitForAYoungOp(YOpInput,AreYsIdentical)
!SUBROUTINE DataInitForAYoungOp initializes Young operators and on the output it gives
!Glob_NumYTerms0, Glob_YCoeff0, Glob_YMatr0, and Glob_NumYTerms1, Glob_YCoeff1, Glob_YMatr1.
!YOpInput = 1 or 2 in the cases L=1 or L=2, respectively.

!==============================================================================================
!Local variables:
    INTEGER, intent(in)                   :: YOpInput
    logical, intent(out)                  :: AreYsIdentical
    INTEGER                               :: n,npart
    INTEGER                               :: i,j,k,p,q,t,s,w,ii,jj,kk
    CHARACTER(1)                          :: c1
    INTEGER                               :: StrLen,NumFactY
    INTEGER                               :: TotNumOfYTerms,CurrNumOfTerms,TotNumOfYHYTerms
    INTEGER                               :: L,R,FirstLPos,LastRPos
    INTEGER                               :: Coeff,Cf3
    logical                               :: AreTermsIdentical
    INTEGER,allocatable,dimension(:)      :: TempSymCoeff,TempSymCoeff1,NumTermsInYOpFact
    INTEGER,allocatable,dimension(:,:,:)  :: TempSymMatr,TempSymMatr1
    INTEGER,allocatable,dimension(:,:)    :: Matr1,Matr2,Matr3,Matr4

    CHARACTER(Glob_YOperatorStringLength),allocatable,dimension(:) :: YOpStr,YHOpStr

!==============================================================================================

    IF (Glob_ProcID==0) Then

      WRITE(*,*)
      WRITE(*,'(2X,A,I0)') 'Initializing Young operator for L= ',YOpInput
      WRITE(*,*)

    EndIF

    SELECT CASE (YOpInput)

    CASE(1)
      Glob_YOperatorString=Glob_YOperatorString1

    CASE(2)
      Glob_YOperatorString=Glob_YOperatorString2

    ENDSELECT

    n=Glob_n
    npart=n+1

! Constructing the Young operator based on the content of
! a string variable Glob_YOperatorString
! First we throw away spaces and multiplication signs
! from Glob_YOperatorString

    StrLen=len_trim(Glob_YOperatorString)

    do i=1,StrLen

      c1=Glob_YOperatorString(i:i)

      IF ((c1==' ').or.(c1=='*')) then

        do j=i,StrLen-1
          Glob_YOperatorString(j:j)=Glob_YOperatorString(j+1:j+1)
        enddo

        Glob_YOperatorString(j:j)=' '

      ENDIF

    enddo

    StrLen=len_trim(Glob_YOperatorString)

!Checking for wrong symbols in Glob_YOperatorString
    do i=1,StrLen

      c1=Glob_YOperatorString(i:i)

      IF ((c1/='1').and.(c1/='2').and.(c1/='3').and.(c1/='4').and.(c1/='5').and. &
          (c1/='6').and.(c1/='7').and.(c1/='8').and.(c1/='9').and.(c1/='0').and. &
          (c1/='P').and.(c1/='+').and.(c1/='-').and.(c1/='*').and.(c1/=')').and. &
          (c1/='(')) THEN

        IF (Glob_ProcID==0) THEN
          WRITE(*,*) '  Error EC0115 in ProgramDataInit: the Young operator expression contains wrong symbols'
          WRITE(*,*) '  contains wrong symbols. Invalid character: "', c1, '" at position ', i
        ENDIF

        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

      ENDIF

    enddo

!Checking if the number of left and right brackets is the same
!and counting how many brackets there are
    L=0
    R=0

    do i=1,StrLen

      IF (Glob_YOperatorString(i:i)==')') R=R+1
      IF (Glob_YOperatorString(i:i)=='(') L=L+1

    enddo

    IF (R/=L) THEN

      IF(Glob_ProcID==0)WRITE(*,*)  &
        '  Error EC0116 in ProgramDataInit: the numer of left and right brackets in the Young operator is different'
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    ENDIF

! NumFactY is the number of factors in the Young operator,
! FirstLPos is the position of the first left bracket,
! LastRPos is the position of the last right bracket,

    IF (R/=0) then

      FirstLPos=scan(Glob_YOperatorString(1:StrLen),'(')
      LastRPos=scan(Glob_YOperatorString(1:StrLen),')',back=.true.)
      NumFactY=R
      i=0

      do j=1,R

        k=0
        i=i+1
        c1=Glob_YOperatorString(i:i)

        do while (c1/='(')
          IF (c1/='*') k=1
          i=i+1
          c1=Glob_YOperatorString(i:i)
        enddo

        IF (k==1) NumFactY=NumFactY+1

        do while (c1/=')')
          i=i+1
          c1=Glob_YOperatorString(i:i)
        enddo

      enddo

      IF (Glob_YOperatorString(StrLen:StrLen)/=')') NumFactY=NumFactY+1

    ELSE

      NumFactY=1

    ENDIF

    SELECT CASE (YOpInput)

    CASE(1)
      Glob_NumFactY1=NumFactY
    CASE(2)
      Glob_NumFactY2=NumFactY

    ENDSELECT

! Splitting Glob_YOperatorString into an array of smaller
! strings, YOpStr. Each column of this array will contain just
! one factor, with no brackets. A '+' or a '-' sign is added in
! front of the first term in a factor if needed. Multiplication
! signs are dropped.

    allocate(YOpStr(NumFactY))
    YOpStr=' '

    IF (R==0) then

      c1=Glob_YOperatorString(1:1)

      IF ((c1/='+').or.(c1/='-')) then
        YOpStr(1)(1:1)='+'
        YOpStr(1)(2:StrLen+1)=Glob_YOperatorString(1:StrLen)
      ELSE
        YOpStr(1)(1:StrLen)=Glob_YOperatorString(1:StrLen)
      ENDIF

    ELSE

      i=1
      k=1
      p=i
      q=0
      c1=Glob_YOperatorString(i:i)

      IF ((c1/='(').and.(c1/='+').and.(c1/='-')) then
        q=1
        YOpStr(k)(1:1)='+'
      ENDIF

      do while (Glob_YOperatorString(i:i)/='(')
        i=i+1
      enddo

      IF (i>1) then

        YOpStr(k)(p+q:i-1+q)=Glob_YOperatorString(p:i-1)
        !IF (YOpStr(k)(i-1+q:i-1+q)=='*') YOpStr(k)(i-1+q:i-1+q)=' '
        k=k+1

      ENDIF

      do j=1,R

        i=i+1
        p=i
        q=0
        c1=Glob_YOperatorString(i:i)

        IF ((c1/=')').and.(c1/='+').and.(c1/='-')) then
          q=1
          YOpStr(k)(1:1)='+'
        ENDIF

        do while (Glob_YOperatorString(i:i)/=')')
          i=i+1
        enddo

        YOpStr(k)(1+q:i+q-p)=Glob_YOperatorString(p:i-1)
        k=k+1
        i=i+1
        c1=Glob_YOperatorString(i:i)

        IF (c1=='*') then

          i=i+1
          c1=Glob_YOperatorString(i:i)
        ENDIF

        IF ((c1/='(').and.(c1/=' '))  then

          p=i
          YOpStr(k)(1:1)='+'

          do while ((c1/='(').and.(c1/=' '))
            i=i+1
            c1=Glob_YOperatorString(i:i)
          enddo

          YOpStr(k)(2:i+1-p)=Glob_YOperatorString(p:i-1)
          k=k+1

        ENDIF

      enddo

    ENDIF

! Checking of similarity of Young operators. If these are the same then AreYsIdentical=.true.
    AreYsIdentical=.false.

    IF ((YOpInput==2).and.(Glob_NumFactY1==NumFactY)) then

      IF(all(Glob_YOpStr1==YOpStr)) THEN

        AreYsIdentical=.true.
        deallocate(YOpStr)

        IF (Glob_ProcID==0) THEN
          WRITE(*,*) '  The Young operators are the same for L=1 and L=2 states '
        ENDIF

        return

      ELSE
        IF (Glob_ProcID==0) THEN
          WRITE(*,*) '        **************************************** '
          WRITE(*,*) '        *              WARNING!!!              * '
          WRITE(*,*) '        *         The Young operators          * '
          WRITE(*,*) '        *          are not the same            * '
          WRITE(*,*) '        *       for L=1 and L=2 states         * '
          WRITE(*,*) '        *     MAKE SURE THEY ARE CORRECT !!!   * '
          WRITE(*,*) '        **************************************** '
          WRITE(*,*)
        ENDIF
      end IF

    ENDIF

!!!
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

        IF (c1=='P') then

          k=0 !k counts the number of Permutations in the current term
          t=0

          do while ((c1/='+').and.(c1/='-').and.(c1/=' '))
            IF (c1=='P') k=k+1
            t=t+1
            c1=YOpStr(s)(j+t:j+t)
          enddo

          do t=1,k
            YHOpStr(i)(j+3*(k-t):j+3*(k-t)+2)=YOpStr(s)(j+3*(t-1):j+3*(t-1)+2)
          enddo

          j=j+3*k

        ELSE

          YHOpStr(i)(j:j)=c1
          j=j+1

        ENDIF

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
        IF ((YOpStr(k)(i:i)=='+').or.(YOpStr(k)(i:i)=='-')) j=j+1
      enddo

      NumTermsInYOpFact(k)=j
      TotNumOfYTerms=TotNumOfYTerms*j
      TotNumOfYHYTerms=TotNumOfYTerms*TotNumOfYTerms
    enddo

    IF (Glob_ProcID==0) THEN
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y operator:     ',TotNumOfYTerms
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y^{+}Y operator:',TotNumOfYHYTerms
    ENDIF

    SELECT CASE (YOpInput)

    CASE(1)
      allocate(Glob_YOpStr1(NumFactY))
      Glob_YOpStr1=YOpStr
      allocate(Glob_NumTermsInYOpFact1(NumFactY))
      Glob_NumTermsInYOpFact1=NumTermsInYOpFact

    CASE(2)
      allocate(Glob_NumTermsInYOpFact2(NumFactY))
      Glob_NumTermsInYOpFact2=NumTermsInYOpFact

    ENDSELECT

!Multiplying all factors in YOpStr and placing actual matrices
!and coefficients in arrays Glob_YMatr and Glob_YCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.

    allocate(Matr1(1:n,1:n))
    allocate(Matr2(1:n,1:n))
    allocate(Matr3(1:n,1:n))
    allocate(Matr4(1:n,1:n))

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

        IF (i-p>1) then
          read(YOpStr(j)(p:i-1),*) Coeff
        ELSE
          IF (YOpStr(j)(i-1:i-1)=='+') then
            Coeff=1
          ELSE
            Coeff=-1
          ENDIF
        ENDIF

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

        IF (j/=NumFactY) then

          IF (k==1) then

            Matr3(1:n,1:n)=Matr1(1:n,1:n)
            Cf3=Coeff

          ELSE
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
          ENDIF

        ELSE

          TempSymMatr(1:n,1:n,k)=Matr1(1:n,1:n)
          TempSymCoeff(k)=Coeff

        ENDIF

      enddo

      IF (j/=NumFactY) then

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

      ENDIF

      ! mark the identical terms (adding their coefficients
      ! and setting all of them but one to zero)
      t=CurrNumOfTerms

      do i=1,CurrNumOfTerms

        IF (TempSymCoeff(i)==0) cycle
        do s=i+1,CurrNumOfTerms

          IF (TempSymCoeff(s)==0) cycle

          IF (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
            TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
            IF (TempSymCoeff(i)==0) t=t-1
            TempSymCoeff(s)=0
            t=t-1
          ENDIF

        enddo
      enddo

      !reallocate arrays containing symmetry terms
      !to allow for multiplication by the next factor
      IF (j/=1) then

        allocate(TempSymCoeff1(t))
        allocate(TempSymMatr1(n,n,t))

        s=0

        do i=1,CurrNumOfTerms
          IF (TempSymCoeff(i)/=0) then

            s=s+1
            TempSymCoeff1(s)=TempSymCoeff(i)
            TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)

          ENDIF
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
      ENDIF

    enddo

    Glob_NumYTerms=t

    SELECT CASE (YOpInput)

    CASE(1)
      Glob_NumYTerms1=Glob_NumYTerms
      allocate(Glob_YCoeff1(Glob_NumYTerms))
      allocate(Glob_YMatr1(n,n,Glob_NumYTerms))

    CASE(2)
      Glob_NumYTerms2=Glob_NumYTerms
      allocate(Glob_YCoeff2(Glob_NumYTerms))
      allocate(Glob_YMatr2(n,n,Glob_NumYTerms))

    ENDSELECT

    s=0

    do i=1,CurrNumOfTerms
      IF (TempSymCoeff(i)/=0) then

        s=s+1

        SELECT CASE (YOpInput)

        CASE(1)
          Glob_YCoeff1(s)=TempSymCoeff(i)
          Glob_YMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        CASE(2)
          Glob_YCoeff2(s)=TempSymCoeff(i)
          Glob_YMatr2(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)

        ENDSELECT

      ENDIF
    enddo

    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)

    IF (Glob_ProcID==0) THEN
      write(*,'(1x,a,1x,i9)') 'Total number of terms in simplified Y operator:        ',Glob_NumYTerms
    ENDIF

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
! TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
! TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)

    SELECT CASE (YOpInput)

    CASE(1)
      TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff1(1:Glob_NumYTerms)
      TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr1(1:n,1:n,1:Glob_NumYTerms)
    CASE(2)
      TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff2(1:Glob_NumYTerms)
      TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr2(1:n,1:n,1:Glob_NumYTerms)

    ENDSELECT

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

        IF (i-p>1) then
          read(YHOpStr(j)(p:i-1),*) Coeff
        ELSE
          IF (YHOpStr(j)(i-1:i-1)=='+') then
            Coeff=1
          ELSE
            Coeff=-1
          ENDIF
        ENDIF

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

        IF (k==1) then

          Matr3(1:n,1:n)=Matr1(1:n,1:n)
          Cf3=Coeff

        ELSE

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

        ENDIF

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

      ! mark the identical terms (adding their coefficients
      ! and setting all of them but one to zero)
      t=CurrNumOfTerms
      do i=1,CurrNumOfTerms

        IF (TempSymCoeff(i)==0) cycle

        do s=i+1,CurrNumOfTerms

          IF (TempSymCoeff(s)==0) cycle

          IF (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then

            TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
            IF (TempSymCoeff(i)==0) t=t-1
            TempSymCoeff(s)=0
            t=t-1

          ENDIF

        enddo

      enddo

      !reallocate arrays containing symmetry terms
      !to allow for multiplication by the next factor
      IF (j/=1) then

        allocate(TempSymCoeff1(t))
        allocate(TempSymMatr1(n,n,t))
        s=0

        do i=1,CurrNumOfTerms
          IF (TempSymCoeff(i)/=0) then

            s=s+1
            TempSymCoeff1(s)=TempSymCoeff(i)
            TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)

          ENDIF
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

      ENDIF
    enddo

    Glob_NumYHYTerms=t

    SELECT CASE (YOpInput)

    CASE(1)
      Glob_NumYHYTerms1=Glob_NumYHYTerms
      allocate(Glob_YHYCoeff1(Glob_NumYHYTerms1))
      allocate(Glob_YHYMatr1(n,n,Glob_NumYHYTerms1))
    CASE(2)
      Glob_NumYHYTerms2=Glob_NumYHYTerms
      allocate(Glob_YHYCoeff2(Glob_NumYHYTerms2))
      allocate(Glob_YHYMatr2(n,n,Glob_NumYHYTerms2))

    ENDSELECT

    s=0
    do i=1,CurrNumOfTerms
      IF (TempSymCoeff(i)/=0) then

        s=s+1

        SELECT CASE (YOpInput)

        CASE(1)
          Glob_YHYCoeff1(s)=TempSymCoeff(i)
          Glob_YHYMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        CASE(2)
          Glob_YHYCoeff2(s)=TempSymCoeff(i)
          Glob_YHYMatr2(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)

        ENDSELECT

      ENDIF
    enddo

    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)

    IF (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in simplified Y^{+}Y operator:   ',Glob_NumYHYTerms
    ENDIF

    deallocate(Matr4)
    deallocate(Matr3)
    deallocate(Matr2)
    deallocate(Matr1)
    deallocate(NumTermsInYOpFact)
    deallocate(YHOpStr)
    deallocate(YOpStr)

  END SUBROUTINE DataInitForAYoungOp

  SUBROUTINE ProgramDataInit()
!SUBROUTINE ProgramDataInit initializes some data needed for
!calculations. It should be CALLed at the start of the program,
!right after reading input/output FILE.

!==============================================================================================
    INTEGER                               :: n,npart
    INTEGER                               :: i,j,k,p,q,t,s,w
    INTEGER                               :: pi,pj,pt,ps
    logical                               :: AreYsIdentical
    REAL(wp)                           :: mk,mi,m0
    INTEGER,allocatable,dimension(:)      :: IdentParticleSet
    INTEGER,allocatable,dimension(:,:)    :: IdentPseuDOPartPairSet
!==============================================================================================
!==============================================================================================

    IF (Glob_ProcID==0)THEN

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(4x,A)') 'Initializing program data'
      WRITE(*,'(A)') '================================================='

    ENDIF

    n=Glob_n
    npart=n+1

!Constructing Glob_MassMatrix
    ALLOCATE(Glob_MassMatrix(n,n))
    Glob_MassMatrix(1:n,1:n)=ONEHALF/Glob_Mass(1)

    DO i=1,n
      Glob_MassMatrix(i,i)=Glob_MassMatrix(i,i)+ONEHALF/Glob_Mass(i+1)
    ENDDO

!Determine the components of vector Glob_bvc
!that is used in evaluation of particle densities
    ALLOCATE(Glob_bvc(n,npart))
    Glob_MassTotal=sum(Glob_Mass(1:npart))

    DO i=1,npart
      Glob_bvc(1:n,i)=-Glob_Mass(2:n+1)/Glob_MassTotal
    ENDDO

    DO i=2,npart
      Glob_bvc(i-1,i)=Glob_bvc(i-1,i)+ONE
    ENDDO

!Determine the mass and the index of the lightest particle
!(reference particle excluded).
!and its index
    k=0
    mk=2*Glob_MassTotal

    DO i=1,n

      IF (Glob_Mass(i+1)<mk) THEN
        k=i
        mk=Glob_Mass(i+1)
      ENDIF

    ENDDO

    m0=Glob_Mass(1)

!alpha = sqrt( 0.5 * (m0^3 + m_k^3)/(m0*m_k*(m0 + m_k)^2) )
    Glob_dmva2 = (m0**3 + mk**3)/(TWO*m0*mk*(m0+mk)**2)
!Glob_dmvB(i,i) = (beta^2 + gamma_i^2)/(alpha^2 * M_ii) - M_ii
    Glob_dmvB(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO

    DO i=1,n

      mi=Glob_Mass(i+1)
      Glob_dmvB(i,i)=( (m0**3+mi**3)*mk*(m0+mk)**2 - (m0**3+mk**3)*mi*(m0+mi)**2 ) / ( TWO*(m0+mi)*(m0**3+mk**3)*m0*mi**2 )

    ENDDO

    Glob_dmvM(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO
    Glob_dmvM(1:n,1:n)=Glob_MassMatrix(1:n,1:n)
    Glob_dmvMB=Glob_dmvM+Glob_dmvB

!-----------------------------------------------------------------------------------
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
!
!-----------------------------------------------------------------------------------

    ALLOCATE(Glob_Transposit(n,n,npart,npart))
!First set all of them to be unit matrices
    Glob_Transposit(1:n,1:n,1:npart,1:npart)=0

    DO i=1,npart
      DO j=1,npart

        DO k=1,n
          Glob_Transposit(k,k,i,j)=1
        ENDDO

      ENDDO
    ENDDO

!Now continue depending on type of transposition (P1i or Pij)
    DO i=2,npart
      Glob_Transposit(1:n,i-1,1,i)=-1
    ENDDO

    DO i=2,npart
      DO j=i+1,npart

        Glob_Transposit(i-1,i-1,i,j)=0
        Glob_Transposit(j-1,j-1,i,j)=0
        Glob_Transposit(j-1,i-1,i,j)=1
        Glob_Transposit(i-1,j-1,i,j)=1
        Glob_Transposit(i-1,i-1,j,i)=0
        Glob_Transposit(j-1,j-1,j,i)=0
        Glob_Transposit(j-1,i-1,j,i)=1
        Glob_Transposit(i-1,j-1,j,i)=1

      ENDDO
    ENDDO

!Calculate Glob_YCoeff and Glob_YMatr for both case L=0 and L=1.
!The order of CALLing is important, it must be 0,1.

    DO i=1,2
      CALL DataInitForAYoungOp(i,AreYsIdentical)
    ENDDO

    IF(AreYsIdentical) THEN

      Glob_NumFactY2=Glob_NumFactY1
      ALLOCATE(Glob_NumTermsInYOpFact2(Glob_NumFactY2))
      Glob_NumTermsInYOpFact2=Glob_NumTermsInYOpFact1

      Glob_NumYTerms2=Glob_NumYTerms1
      Glob_NumYHYTerms2=Glob_NumYHYTerms1

      ALLOCATE(Glob_YCoeff2(Glob_NumYTerms2))
      ALLOCATE(Glob_YMatr2(n,n,Glob_NumYTerms2))

      Glob_YCoeff2=Glob_YCoeff1
      Glob_YMatr2=Glob_YMatr1

      ALLOCATE(Glob_YHYCoeff2(Glob_NumYTerms2))
      ALLOCATE(Glob_YHYMatr2(n,n,Glob_NumYTerms2))

      Glob_YHYCoeff2=Glob_YHYCoeff1
      Glob_YHYMatr2=Glob_YHYMatr1

    ENDIF

!Now we determine which particles are identical. This determination
!is based on the input values of masses and charges only. The information
!about the sets of identical particles may be needed for proper
!symmetrization of expectation values of operators that involve
!two-particle quantities (such as interparticle distances).
!The set of particles to which particle i belongs is labelled by IdentParticleSet(i)
!IF IdentParticleSet(i)=IdentParticleSet(j) THEN it means that
!particles i and j are identical. The largest value in IdentParticleSet
!gives the total number of identical particle sets

    ALLOCATE(IdentParticleSet(npart))
    IdentParticleSet(1)=1
    k=1

    DO i=2,npart

      s=0
      j=0

      DO while ((j<i-1).and.(s==0))

        j=j+1

        IF (j>1) THEN

          IF ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge(j-1)==Glob_PseudoCharge(i-1))) THEN

            IdentParticleSet(i)=IdentParticleSet(j)
            s=1

          ENDIF

        ELSE

          !j=1 CASE
          IF ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge0==Glob_PseudoCharge(i-1))) THEN
            IdentParticleSet(i)=IdentParticleSet(j)
            s=1
          ENDIF

        ENDIF

      ENDDO

      IF (s==0) THEN
        k=k+1
        IdentParticleSet(i)=k
      ENDIF

    ENDDO

    Glob_NumOfIdentPartSets=MAXval(IdentParticleSet(1:npart))

!Below we determine which pairs of pseuDOparticles are identical.
!The information about this is stored in array IdentPseuDOPartPairSet(1:n,1:n)
!Diagonal elements DO not actually designate pairs of pseuDOparticles but
!rather a single pseuDOparticle, which corresponds to a certain pair of particles.
!IF IdentPseuDOPartPairSet(i,j)=IdentPseuDOPartPairSet(k,l) THEN
!it means these pairs ij and kl should be equivalent.
!The largest value of array IdentPseuDOPartPairSet gives the number
!of nonequivalent pairs.

    ALLOCATE(IdentPseuDOPartPairSet(1:n,1:n))

    IdentPseuDOPartPairSet(1:n,1:n)=0
    k=0

    DO i=1,n
      DO j=i,n

        IF (i==j) THEN
          pi=1; pj=j+1
        ELSE
          pi=i+1; pj=j+1
        ENDIF

        w=0

        DO s=1,i

          IF (s==i) THEN
            q=j-1
          ELSE
            q=n
          ENDIF

          DO t=s,q

            IF (w==1) cycle

            IF (s==t) THEN
              ps=1; pt=t+1
            ELSE
              ps=s+1; pt=t+1
            ENDIF

            IF ((IdentParticleSet(ps)==IdentParticleSet(pi)).and. &
                (IdentParticleSet(pt)==IdentParticleSet(pj))) THEN
              w=1
              IdentPseuDOPartPairSet(i,j)=IdentPseuDOPartPairSet(s,t)
            ENDIF

          ENDDO

        ENDDO

        IF (w==0) THEN
          k=k+1
          IdentPseuDOPartPairSet(i,j)=k
        ENDIF

      ENDDO
    ENDDO

    Glob_NumOfNoneqvPairSets=MAXval(IdentPseuDOPartPairSet(1:n,1:n))

!Now we create arrays Glob_NumOfPartInIdentPartSet
!and Glob_IdentPartList

    ALLOCATE(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
    ALLOCATE(Glob_IdentPartList(npart,Glob_NumOfIdentPartSets))
    Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)=0
    Glob_IdentPartList(1:npart,1:Glob_NumOfIdentPartSets)=0

    DO i=1,npart
      k=IdentParticleSet(i)
      Glob_NumOfPartInIdentPartSet(k)=Glob_NumOfPartInIdentPartSet(k)+1
      Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k),k)=i
    ENDDO

!Create arrays Glob_NumOfPairsInEqvPairSet and
!Glob_EqvPairList
    ALLOCATE(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
    ALLOCATE(Glob_EqvPairList(2,n*(n+1)/2,Glob_NumOfNoneqvPairSets))
    Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets)=0
    Glob_EqvPairList(1:2,1:n*(n+1)/2,1:Glob_NumOfNoneqvPairSets)=0

    DO i=1,n
      DO j=i,n

        k=IdentPseuDOPartPairSet(i,j)
        Glob_NumOfPairsInEqvPairSet(k)=Glob_NumOfPairsInEqvPairSet(k)+1
        Glob_EqvPairList(1,Glob_NumOfPairsInEqvPairSet(k),k)=i
        Glob_EqvPairList(2,Glob_NumOfPairsInEqvPairSet(k),k)=j

      ENDDO
    ENDDO

    DEALLOCATE(IdentPseuDOPartPairSet)
    DEALLOCATE(IdentParticleSet)

  END SUBROUTINE ProgramDataInit

  SUBROUTINE ComputeTranDipoL1L2()

! SUBROUTINE ComputeTranDipoL1L2 computes expectation value of
!                      < L=1 | D{z} | L=2 >
!
!
!                                            c1_k * c2_l
!  <psi1|D{z}|psi2> = SUM{k,l}---------------------------------------*
!                                Sqrt( <psi1|psi1> * <psi2|psi2> )
!
!                        m_i         < Y1*f_k1 | z_i | Y2*f_l2 >
!  SUM{i}(q_i - Q_tot * -----)*-----------------------------------------------
!                        m0    Sqrt( <f_k1|f_k1> ) * Sqrt( <f_l2|f_l2> )
!
!
!
!        <psi{L1}|psi{L1}>   and   <psi{L2}|psi{L12}>
!
!        computed by SUBROUTINE ElementS1 ElementS2
!
!
!
!                              m_i         < Y1*f_k1 | z_i | Y2*f_l2 >
!        SUM{i}(q_i - Q_tot * -----)*-----------------------------------------------
!                              m0    Sqrt( <f_k1|f_k1> ) * Sqrt( <f_l2|f_l2> )
!
!        computed by SUBROUTINE MatrixElemenTranDipolMomen
!
!
! Expressions:
!
!        TDkl_1                        =        output of the MatrixElemenTranDipolMomen SUBROUTINE
!        TDkl                    =        transition dipole momentum after applying Glob_YCoeff and normalization factor
!        k,l                                =        Glob_CurrBasisSize1,Glob_CurrBasisSize2
!        Glob_ExpVals         =   (Sum_{k,l} c1_k * TDkl * Glob_c2(l))/sqrt(diagS1*diagS2)
!
!
! Verbosity levels:
!   0  -  Final results only
!   1  -  Progress info, dimensions, timing, sanity checks
!   2  -  Overlap arrays, intermediate sums, per-row contributions
!   3  -  Full element-by-element detail (very large output)

    IMPLICIT NONE
!==============================================================================================
! local variables
    INTEGER        ::  count1, count2, count_sec, time_TranDipoL, rate
    INTEGER        ::  n,np,npt
    INTEGER        ::  k,l,i,j,indx
    REAL(wp)    ::  ExpVal1, ExpVal2, ExpValLoc1, ExpValLoc2
    REAL(wp)    ::  temp0
    REAL(wp)    ::  Skk, Sll
    REAL(wp)    ::  TranDipolLength_kl_element, TranDipolLength_kl_Loc, TranDipolLength_kl
    REAL(wp)    ::  TranDipolVelocity_kl_element, TranDipolVelocity_kl_Loc, TranDipolVelocity_kl

    REAL(wp),allocatable,dimension(:)  ::  diagS1, diagS2, temp1

!==============================================================================================
    CALL SYSTEM_CLOCK(count1, rate)

    n=Glob_n
    np=Glob_np
    npt=Glob_npt

    IF(Glob_ProcID==0) THEN

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(4x,A)') 'Transition Dipole Moment calculation for '
      WRITE(*,'(6x,A)') '< Po(y) | D_x | Pe(z) > started '
      WRITE(*,'(A)') '================================================='

      IF(Verbose >= 1) THEN
        WRITE(*,*)
        WRITE(*,'(A)')          '   Run configuration:  '
        WRITE(*,'(A)')          '  ----------------------------'
        WRITE(*,'(4X,A,I10)')    'Number of particles (n)  = ', n
        WRITE(*,'(4X,A,I10)')    'Nonlinear parameters (np)= ', np
        WRITE(*,'(4X,A,I10)')    'Total parameters (npt)   = ', npt
        WRITE(*,'(4X,A,I10)')    'Basis size, first state  = ', Glob_CurrBasisSize1
        WRITE(*,'(4X,A,I10)')    'Basis size, final state  = ', Glob_CurrBasisSize2
        WRITE(*,'(4X,A,I10)')    'Total (k,l) pairs        = ', Glob_CurrBasisSize1 * Glob_CurrBasisSize2
        WRITE(*,'(4X,A,I10)')    'Y-coefficient terms(1st) = ', Glob_NumYTerms1
        WRITE(*,'(4X,A,I10)')    'Y-coefficient terms(2nd) = ', Glob_NumYTerms2
        WRITE(*,*)
      ENDIF

    ENDIF

!=============================================================================
! computing diagonal elements of overlap matrix L=1
!   <psi1|psi1>
!=============================================================================
    ALLOCATE(diagS1(Glob_CurrBasisSize1))
    ALLOCATE(temp1(Glob_CurrBasisSize1))

    diagS1=ZERO
    temp1=ZERO

    IF(Glob_ProcID==0) THEN

      WRITE(*,*)
      WRITE(*,'(2x,A)') 'Section 1:'
      WRITE(*,'(4x,A)') 'Diagonal overlap (P^o)'
      WRITE(*,'(2x,A)') '----------------------------'

      IF(Verbose >= 1)THEN
        WRITE(*,'(4X,A,I4)')    'Number of permutation terms (YHY): ',Glob_NumYHYTerms1
      ENDIF

    ENDIF

    CALL SYSTEM_CLOCK(count_sec)

    DO k=1+Glob_ProcID,Glob_CurrBasisSize1,Glob_NumOfProcs

      temp0=ZERO

      DO i=1,Glob_NumYHYTerms1
        CALL OverLapElement_S_Po(Glob_Index1(k), Glob_NonlinParam1(1:np,k),Glob_YHYMatr1(1:n,1:n,i), Skk)
        temp0=temp0+Glob_YHYCoeff1(i)*Skk

        IF(Verbose >= 3) THEN
          WRITE(*,'(8X,A,I4,A,I4,A,ES14.6,A,ES14.6)') &
            'k=', k, '  YHY term i=', i, '   Skk=', Skk, '   coeff=', Glob_YHYCoeff1(i)
        ENDIF

      ENDDO

      temp1(k)=temp0

      IF(Verbose >= 3) THEN
        WRITE(*,'(6X,A,I4,A,ES16.8)') '>> k=', k, '   diagS1(k) [local] = ', temp0
      ENDIF

    ENDDO

    CALL MPI_ALLREDUCE(temp1,diagS1,Glob_CurrBasisSize1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    DEALLOCATE(temp1)

    CALL SYSTEM_CLOCK(count2)

    IF(Glob_ProcID==0) THEN

      WRITE(*,'(4X,A)') 'Diagonal overlap S1 ...  done!'

      IF(Verbose >= 1) THEN
        WRITE(*,*)
        WRITE(*,'(4X,A,ES16.8)') 'min( diagS1 ) = ', MINVAL(diagS1)
        WRITE(*,'(4X,A,ES16.8)') 'max( diagS1 ) = ', MAXVAL(diagS1)
        WRITE(*,'(4X,A,F10.3,A)') 'Wall time     = ', REAL(count2 - count_sec) / REAL(rate), ' seconds'
      ENDIF

      IF(Verbose >= 2) THEN
        WRITE(*,'(4X,A)') 'Full diagS1 array:'
        DO k = 1, Glob_CurrBasisSize1
          WRITE(*,'(8X,A,I4,A,ES22.14)') 'diagS1(', k, ') = ', diagS1(k)
        ENDDO
        WRITE(*,'(4X,A)') 'YHY coefficients L=1:'
        DO i = 1, Glob_NumYHYTerms1
          WRITE(*,'(8X,A,I4,A,ES22.14)') 'YHYCoeff1(', i, ') = ', Glob_YHYCoeff1(i)
        ENDDO
      ENDIF

      WRITE(*,*)

    ENDIF

!=============================================================================
! computing diagonal elements of overlap matrix P^e
!   <psi2|psi2>
!=============================================================================

    IF(Glob_ProcID==0) THEN

      WRITE(*,*)
      WRITE(*,'(2x,A)') 'Section 2:'
      WRITE(*,'(4x,A)') 'Diagonal overlap (P^e)'
      WRITE(*,'(2x,A)') '----------------------------'

      IF(Verbose >= 1)THEN
        WRITE(*,'(4X,A,I4)')    'Number of permutation terms (YHY): ',Glob_NumYHYTerms2
      ENDIF

    ENDIF

    ALLOCATE(diagS2(Glob_CurrBasisSize2))
    ALLOCATE(temp1(Glob_CurrBasisSize2))

    diagS2=ZERO
    temp1=ZERO

    CALL SYSTEM_CLOCK(count_sec)

    DO l=1+Glob_ProcID,Glob_CurrBasisSize2,Glob_NumOfProcs

      temp0=ZERO

      DO i=1,Glob_NumYHYTerms2

        CALL OverLapElement_S_Pe(Glob_Index2(l,1), Glob_Index2(l,2), Glob_NonlinParam2(1:npt, l), Glob_YHYMatr2(1:n, 1:n, i), Sll)
        temp0=temp0+Glob_YHYCoeff2(i)*Sll

        IF(Verbose >= 3) THEN
          WRITE(*,'(8X,A,I4,A,I4,A,ES14.6,A,ES14.6)') &
            'l=', l, '  YHY term i=', i, '   Sll=', Sll, '   coeff=', Glob_YHYCoeff2(i)
        ENDIF

      ENDDO

      temp1(l)=temp0

      IF(Verbose >= 3) THEN
        WRITE(*,'(6X,A,I4,A,ES16.8)') '>> l=', l, '   diagS2(l) [local] = ', temp0
      ENDIF

    ENDDO

    CALL MPI_ALLREDUCE(temp1,diagS2,Glob_CurrBasisSize2,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    DEALLOCATE(temp1)

    CALL SYSTEM_CLOCK(count2)

    IF(Glob_ProcID==0) THEN

      WRITE(*,'(4X,A)') 'Diagonal overlap S2 ...  done!'

      IF(Verbose >= 1) THEN
        WRITE(*,*)
        WRITE(*,'(4X,A,ES16.8)') 'min( diagS2 ) = ', MINVAL(diagS2)
        WRITE(*,'(4X,A,ES16.8)') 'max( diagS2 ) = ', MAXVAL(diagS2)
        WRITE(*,'(4X,A,F10.3,A)') 'Wall time     = ', REAL(count2 - count_sec) / REAL(rate), ' seconds'
      ENDIF

      IF(Verbose >= 2) THEN

        WRITE(*,'(4X,A)') 'Full diagS2 array:'
        DO l = 1, Glob_CurrBasisSize2
          WRITE(*,'(8X,A,I4,A,ES22.14)') 'diagS2(', l, ') = ', diagS2(l)
        ENDDO

        WRITE(*,'(4X,A)') 'YHY coefficients L=2:'
        DO i = 1, Glob_NumYHYTerms2
          WRITE(*,'(8X,A,I4,A,ES22.14)') 'YHYCoeff2(', i, ') = ', Glob_YHYCoeff2(i)
        ENDDO

      ENDIF

      WRITE(*,*)

    ENDIF

!=============================================================================
! main loop
!=============================================================================

    IF(Glob_ProcID==0) THEN
      WRITE(*,*)
      WRITE(*,'(A)') '   Section 3:'
      WRITE(*,'(A)') '    Oscillator Strength'
      WRITE(*,'(A)') '  ----------------------------'
      WRITE(*,*)
    ENDIF

    indx = ZERO
    ExpValLoc1 = ZERO
    ExpValLoc2 = ZERO
    TranDipolLength_kl=ZERO
    TranDipolVelocity_kl = ZERO

    CALL SYSTEM_CLOCK(count_sec)

    DO k=1,Glob_CurrBasisSize1
      DO l=1,Glob_CurrBasisSize2

        indx=indx+1

        IF(mod(indx,Glob_NumOfProcs)==Glob_ProcID) THEN

          TranDipolLength_kl_Loc   = ZERO
          TranDipolVelocity_kl_Loc = ZERO
          TranDipolLength_kl   = ZERO
          TranDipolVelocity_kl = ZERO

          IF(Verbose >= 3) THEN
            WRITE(*,*)
            WRITE(*,'(6X,A)') '- - - - - - - - - - - - - - - - - - - - - - - - -'
            WRITE(*,'(6X,A,I4,A,I4)')  'Pair  k=', k, '   l=', l
            WRITE(*,'(6X,A,ES14.6)')   'c1(k)  = ', Glob_c1(k)
            WRITE(*,'(6X,A,ES14.6)')   'c2(l)  = ', Glob_c2(l)
            WRITE(*,'(6X,A)')          'Nonlinear params L=1:'
            WRITE(*,'(8X,10ES13.5)')   Glob_NonlinParam1(1:np,k)
            WRITE(*,'(6X,A)')          'Nonlinear params L=2:'
            WRITE(*,'(8X,10ES13.5)')   Glob_NonlinParam2(1:npt,l)
          ENDIF

          DO i=1,Glob_NumYTerms1
            DO j=1,Glob_NumYTerms2

              CALL MatElemTranDipoleMoment_RG_1P_2P(                            &
                Glob_Index1(k),             &
                Glob_NonlinParam1(1:np,k),  &
                Glob_YMatr1(1:n,1:n,i),     &
                Glob_Index2(l,1),           &
                Glob_Index2(l,2),           &
                Glob_NonlinParam2(1:np,l),  &
                Glob_YMatr2(1:n,1:n,j),     &
                TranDipolLength_kl_element, &
                TranDipolVelocity_kl_element)

              IF(Verbose >= 3) THEN
                WRITE(*,'(10X,A,I3,A,I3,A,ES14.6,A,ES14.6)') &
                  'Y(i=', i, ', j=', j, ')   Len_elem=', TranDipolLength_kl_element, &
                  '   Vel_elem=', TranDipolVelocity_kl_element
                WRITE(*,'(10X,A,ES14.6,A,ES14.6)') &
                  '                     YCoeff1(i)=', Glob_YCoeff1(i), &
                  '   YCoeff2(j)=', Glob_YCoeff2(j)
              ENDIF

              TranDipolLength_kl_Loc = TranDipolLength_kl_Loc + &
                                       Glob_YCoeff1(i) * TranDipolLength_kl_element * Glob_YCoeff2(j)

              TranDipolVelocity_kl_Loc = TranDipolVelocity_kl_Loc + &
                                         Glob_YCoeff1(i) * TranDipolVelocity_kl_element * Glob_YCoeff2(j)

            ENDDO
          ENDDO

          IF(diagS1(k)*diagS2(l) <= ZERO) THEN
            WRITE(*,'(A)')            '  *** ERROR: negative overlap product ***'
            WRITE(*,'(4X,A,I4,A,I4)') '  at  k=',k,'   l=',l
            WRITE(*,'(4X,A,ES22.10)') '  diagS1(k) = ', diagS1(k)
            WRITE(*,'(4X,A,ES22.10)') '  diagS2(l) = ', diagS2(l)
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
          ENDIF

          TranDipolLength_kl = TranDipolLength_kl_Loc / sqrt(diagS1(k)*diagS2(l))
          ExpValLoc1 = ExpValLoc1 + Glob_c1(k) * TranDipolLength_kl * Glob_c2(l)

          TranDipolVelocity_kl = TranDipolVelocity_kl_Loc / sqrt(diagS1(k)*diagS2(l))
          ExpValLoc2 = ExpValLoc2 + Glob_c1(k) * TranDipolVelocity_kl * Glob_c2(l)

          IF(Verbose >= 3)THEN
            WRITE(*,'(8X,A,ES14.6)') 'S1(kk)                 = ', diagS1(k)
            WRITE(*,'(8X,A,ES14.6)') 'S2(ll)                 = ', diagS2(l)
            WRITE(*,'(8X,A,ES14.6)') 'TranDipolLength_kl     = ', TranDipolLength_kl
            WRITE(*,'(8X,A,ES14.6)') 'TranDipolVelocity_kl   = ', TranDipolVelocity_kl
            WRITE(*,'(8X,A,ES14.6)') 'Running ExpValLoc1     = ', ExpValLoc1
            WRITE(*,'(8X,A,ES14.6)') 'Running ExpValLoc2     = ', ExpValLoc2
          ENDIF

        ENDIF
      ENDDO

      IF(Verbose >= 2 .AND. Glob_ProcID==0) THEN
        WRITE(*,'(4X,A,I4,A,ES16.8,A,ES16.8)') &
          'After row k=', k, '   ExpValLoc1=', ExpValLoc1, '   ExpValLoc2=', ExpValLoc2
      ENDIF

    ENDDO

    IF(Verbose >= 2 .AND. Glob_ProcID==0)THEN
      WRITE(*,*)
      WRITE(*,'(4X,A)') '--- Local sums before MPI reduce ---'
      WRITE(*,'(4X,A,ES16.8)') 'ExpValLoc1 (Length)   = ', ExpValLoc1
      WRITE(*,'(4X,A,ES16.8)') 'ExpValLoc2 (Velocity) = ', ExpValLoc2
    ENDIF

    CALL MPI_ALLREDUCE(ExpValLoc1,ExpVal1,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    Glob_ExpVals1 = ExpVal1

    CALL MPI_ALLREDUCE(ExpValLoc2,ExpVal2,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    Glob_ExpVals2 = ExpVal2

    DEALLOCATE(diagS1)
    DEALLOCATE(diagS2)

    CALL SYSTEM_CLOCK(count2, rate)
    time_TranDipoL = count2 - count1

    IF(Glob_ProcID==0)THEN

      WRITE(*,'(A)') '    Transition Dipole Moment calculation done!'

      IF(Verbose >= 1)THEN
        WRITE(*,'(4X,A,F10.3,A)') '  Wall time     = ', REAL(count2 - count_sec) / REAL(rate), ' seconds'
        WRITE(*,'(4X,A,F10.3,A)') 'Total wall time = ', REAL(time_TranDipoL) / REAL(rate), ' seconds'
        WRITE(*,*)
      ENDIF

    ENDIF

  END SUBROUTINE ComputeTranDipoL1L2

end module workproc
