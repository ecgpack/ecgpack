module workproc
!This module contains basic work subroutines
  use matelem
  implicit none

contains

  subroutine Readwf0wf1()
!Subroutine Readwf0wf1 reads data (number of particle,
!mass, charge, nonlinear variational
!parameters and other information) from the wave function
!file whose name is specIFied by global variable
!Glob_wf File files.

!Local variables:
    integer        :: OpenFileErr
    integer        :: ReadInt,ReadErr
    integer        :: particle_n0,particle_n1
    real(wp)    :: Mass0(Glob_AllowedNumOfParticles),Mass1(Glob_AllowedNumOfParticles)
    real(wp)    :: PseudoCharge0(Glob_AllowedNumOfParticles),PseudoCharge1(Glob_AllowedNumOfParticles)
    real(wp)    :: PseudoCharge00,PseudoCharge01
    real(wp)    :: ReadReal0,ReadReal1
    real(wp)    :: RepulsionScalingParamPlus0,RepulsionScalingParamPlus1
    real(wp)    :: RepulsionScalingParamMinus0,RepulsionScalingParamMinus1
    real(wp)    :: RepulsionScalingParam0,RepulsionScalingParam1
    real(wp)    :: AttractionScalingParam0,AttractionScalingParam1
    logical        :: AttrScalParamSupplied0,AttrScalParamSupplied1
    integer        :: WorkInt(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
    integer        :: WorkInt0(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
    integer        :: WorkInt1(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
    integer        :: i,j,k,l,Line0,Line1
    character(70)  :: ReadChar
    logical        :: ErrorInDataFile !,IsDRMCStep

    ErrorInDataFile=.false.
    Glob_PiRaised3n2=Glob_Pi**((3*Glob_n)/TWO)

! opening the wave function files of initial and final states
! Glob_WFfile0 initial state's file
! Glob_WFfile1 final state's file
    IF (Glob_ProcID==0) then

!        initial state
      open(1,file=Glob_WFfile0,status='old',iostat=OpenFileErr)
      IF (OpenFileErr/=0) then
        write(*,*) ' '
        write (*,*) ' Error EC0210 in opening wave function of inital state, '
        write (*,*) ' "',trim(adjustl(Glob_WFfile0)), '"  file not found !!!'
        write(*,*) ' '
        ErrorInDataFile=.true.
      EndIF

!        final state
      open(2,file=Glob_WFfile1,status='old',iostat=OpenFileErr)
      IF (OpenFileErr/=0) then
        write(*,*) ' '
        write (*,*) ' Error EC0211 in opening wave function of final state, '
        write (*,*) ' "',trim(adjustl(Glob_WFfile1)), '" file not found !!!'
        write(*,*) ' '
        ErrorInDataFile=.true.
      EndIF

    EndIF
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop

    IF (Glob_ProcID==0) Line0=0
    IF (Glob_ProcID==0) Line1=0

! Reading Header
    IF (Glob_ProcID==0) then
!        initial state
      read(1,*) ReadChar(1:30)
      write(*,*) 'Reading " ',ReadChar(1:5) ,'" of the initial state from :    ', trim(adjustl(Glob_WFfile0))
      Line0=Line0+1

!        final state
      read(2,*) ReadChar(1:30)
      write(*,*) 'Reading " ',ReadChar(1:5) ,'" of the final   state from :    ', trim(adjustl(Glob_WFfile1))
      Line1=Line1+1
      write (*,*)

    EndIF

!Reading the number of particle from the wave function files
    IF (Glob_ProcID==0) then

!        initial state
      read(1,*) ReadChar(1:9),ReadInt
      Line0=Line0+1
!        particle_n0 is the number of pseudoparticles in Glob_WFfile0
      particle_n0=ReadInt-1
      IF ((particle_n0<1).or.(ReadChar(1:9)/='PARTICLES')) then
        write(*,*) ' '
        write(*,*) 'Error EC0212 in reading the number of particles of the initial state, line: ',Line0
        write(*,*) ' '
        ErrorInDataFile=.true.
      EndIF

!        final state
      read(2,*) ReadChar(1:9),ReadInt
      Line1=Line1+1
!        particle_n1 is the number of pseudoparticles in Glob_WFfile1
      particle_n1=ReadInt-1
      IF ((particle_n1<1).or.(ReadChar(1:9)/='PARTICLES')) then
        write(*,*) ' '
        write(*,*) 'Error EC0213 in reading the number of particles of the final state, line: ',Line1
        write(*,*) ' '
        ErrorInDataFile=.true.
      EndIF

!     comparing the number of particle in Glob_WFfile0 and Glob_WFfile1 files
      IF (particle_n0/=particle_n1) then
        write(*,*) ' '
        write(*,*) 'the number of particles in initial and final states is not the same !!!'
        write(*,*) ' '
        ErrorInDataFile=.true.
      Else
        IF (particle_n0/=Glob_AllowedNumOfPseudoParticles) then
          write(*,*) ' '
          write (*,*) 'The version of the code you are running was compiled for the case'
          write (*,*) 'when the number of particles in the system is equal to', &
            Glob_AllowedNumOfParticles
          write (*,*) 'while the number of particles specIFied in the wave function files is',particle_n0+1
          write (*,*) 'Please make appropriate changes. Program will now stop.'
          write(*,*) ' '
          ErrorInDataFile=.true.
        EndIF
        Glob_n=particle_n0
      EndIF

    EndIF
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    call MPI_BCAST(Glob_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    Glob_2Raised3n2=TWO**((3*Glob_n)/TWO)
    Glob_np=Glob_n*(Glob_n+1)/2
    Glob_npt=Glob_np

! Reading the masses of particles from the wave function files
! Mass0 masses which are read from the Glob_WFfile0
! Mass1 masses which are read from the Glob_WFfile1
    allocate(Glob_Mass(Glob_n+1))
    IF (Glob_ProcID==0) then

!        initial state
      read(1,*) ReadChar(1:6),Mass0(1:Glob_n+1)
      Line0=Line0+1
      IF (ReadChar(1:6)/='MASSES') then
        write(*,*) ' '
        write(*,*) 'Error EC0214 in reading masses of particles of initial state, line: ',Line0
        write(*,*) ' '
      EndIF
!        final state
      read(2,*) ReadChar(1:6),Mass1(1:Glob_n+1)
      Line1=Line1+1
      IF (ReadChar(1:6)/='MASSES') then
        write(*,*) ' '
        write(*,*) 'Error EC0215 in reading masses of particles of final state, line: ',Line1
        write(*,*) ' '
      EndIF

!        comparison of the masses of two files
      Do i=1,Glob_n+1
        IF (Mass0(i)/=Mass1(i))then
          write(*,*) ' '
          write(*,*) 'the mass of the',i,'th particle in initial and final states is not the same !!! '
          write(*,*) ' '
          ErrorInDataFile=.true.
        EndIF
      EndDo
!        Glob_Mass=Mass0=Mass1
      Glob_Mass=Mass0
      write(*,'(1x,a6)',advance='no') ReadChar(1:6)
      call writerealarradv(6,Glob_Mass,Glob_n+1)

    EndIF
    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    call MPI_BCAST(Glob_Mass,Glob_n+1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

! Reading the charges of particles from the wave function files
! PseudoCharge0   charges of the electrons which are read from the Glob_WFfile0
! PseudoCharge1   charges of the electrons which are read from the Glob_WFfile1
! PseudoCharge00 charge  of the nucleus   which are read from the Glob_WFfile0
! PseudoCharge01 charge  of the nucleus   which are read from the Glob_WFfile0
    allocate(Glob_PseudoCharge(Glob_n))
    IF (Glob_ProcID==0) then

!        initial state
      read(1,*) ReadChar(1:7),PseudoCharge00,PseudoCharge0(1:Glob_n)
      Line0=Line0+1
      IF (ReadChar(1:7)/='CHARGES') then
        write(*,*) ' '
        write(*,*) 'Error EC0216 in reading charges of particles of initial state, line: ',Line0
        write(*,*) ' '
      EndIF

!        final state
      read(2,*) ReadChar(1:7),PseudoCharge01,PseudoCharge1(1:Glob_n)
      Line1=Line1+1
      IF (ReadChar(1:7)/='CHARGES') then
        write(*,*) ' '
        write(*,*) 'Error EC0217 in reading charges of particles of final state, line: ',Line1
        write(*,*) ' '
      EndIF

!        comparing charges of the two files
      IF (PseudoCharge00/=PseudoCharge01)then
        write(*,*) ' '
        write(*,*) 'the charge of first particle in initial and final states is not the same !!!'
        write(*,*) ' '
        ErrorInDataFile=.true.
      EndIF

      Do i=1,Glob_n
        IF (PseudoCharge0(i)/=PseudoCharge1(i))then
          write(*,*) ' '
          write(*,*) 'the charge of the',i+1,' in initial and final states is not the same !!! '
          write(*,*) ' '
          ErrorInDataFile=.true.
        EndIF
      EndDo
!        Glob_PseudoCharge0=PseudoCharge00=PseudoCharge01
      Glob_PseudoCharge0=PseudoCharge00
!        Glob_PseudoCharge=PseudoCharge0=PseudoCharge1
      Glob_PseudoCharge=PseudoCharge0
      write(*,'(1x,a7)',advance='no') ReadChar(1:7)
      call writereal(6,Glob_PseudoCharge0)
      call writerealarradv(6,Glob_PseudoCharge,Glob_n)
    EndIF

    call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    IF (ErrorInDataFile) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
    call MPI_BCAST(Glob_PseudoCharge0,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_PseudoCharge,Glob_n,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

! Reading Young operators
    IF (Glob_ProcID==0) then

!        initial state
      read(1,*) ReadChar(1:9),Glob_YOperatorString0
      Line0=Line0+1
      IF (ReadChar(1:9)/='SYMMETRY') then
        write(*,*) ' '
        write(*,*) ReadChar(1:9)
        write(*,*) 'Error EC0218 in reading symmetry of initial state, line: ',Line0
        write(*,*) ' '
      EndIF
      Do i=1,Glob_YOperatorStringLength
        WorkInt(i)=ichar(Glob_YOperatorString0(i:i))
      EndDo
      WorkInt0=WorkInt

!        final state
      read(2,*) ReadChar(1:9),Glob_YOperatorString1
      Line1=Line1+1
      IF (ReadChar(1:9)/='SYMMETRY') then
        write(*,*) ' '
        write(*,*) 'Error EC0219 in reading symmetry of final state, line: ',Line1
        write(*,*) ' '
      EndIF
      Do i=1,Glob_YOperatorStringLength
        WorkInt(i)=ichar(Glob_YOperatorString1(i:i))
      EndDo
      WorkInt1=WorkInt

    EndIF
    call MPI_BCAST(WorkInt0,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(WorkInt1,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    Do i=1,Glob_YOperatorStringLength
      Glob_YOperatorString0(i:i)=char(WorkInt0(i))
      Glob_YOperatorString1(i:i)=char(WorkInt1(i))
    EndDo

! Reading number of the basis functions for each state
    IF (Glob_ProcID==0) then

!        initial state
      read(1,*) ReadChar(1:10),Glob_CurrBasisSize0
      Line0=Line0+1
      IF (ReadChar(1:10)/='BASIS_SIZE') then
        write(*,*) ' '
        write(*,*) 'Error EC0220 in reading BASIS_SIZE of initial state, line: ',Line0
        write(*,*) ' '
      EndIF
      write(*,*) 'INITIAL STATE BASIS SIZE:  ',Glob_CurrBasisSize0

!         final state
      read(2,*) ReadChar(1:10),Glob_CurrBasisSize1
      Line1=Line1+1
      IF (ReadChar(1:10)/='BASIS_SIZE') then
        write(*,*) ' '
        write(*,*) 'Error EC0221 in reading BASIS_SIZE of initial state, line: ',Line1
        write(*,*) ' '
      EndIF
      write(*,*) 'FINAL STATE BASIS SIZE:    ',Glob_CurrBasisSize1
    EndIF
    call MPI_BCAST(Glob_CurrBasisSize0,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_CurrBasisSize1,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

! Reading energy of each state
    IF (Glob_ProcID==0) then
!        initial state
      read(1,*) ReadChar(1:14),Glob_CurrEnergy0
      Line0=Line0+1
      IF (ReadChar(1:14)/='CURRENT_ENERGY') then
        write(*,*) ' '
        write(*,*) 'Error EC0222 in reading CURRENT_ENERGY of initial state, line: ',Line0
        write(*,*) ' '
      EndIF
      write(*,'(a33)',advance='no')  ' INITIAL STATE CURRENT_ENERGY:    '
      call writerealadv(6,Glob_CurrEnergy0)

!        final state
      read(2,*) ReadChar(1:14),Glob_CurrEnergy1
      Line1=Line1+1
      IF (ReadChar(1:14)/='CURRENT_ENERGY') then
        write(*,*) ' '
        write(*,*) 'Error EC0223 in reading CURRENT_ENERGY of final state, line: ',Line1
        write(*,*) ' '
      EndIF
      write(*,'(a33)',advance='no')  ' FINAL STATE  CURRENT_ENERGY:    '
      call writerealadv(6,Glob_CurrEnergy1)
    EndIF
    call MPI_BCAST(Glob_CurrEnergy0,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_CurrEnergy1,1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    IF (Glob_ProcID==0) then
      read(1,*) ReadChar
      read(2,*) ReadChar
    EndIF

! Reading the parameters of the basis functions
    allocate(Glob_c0(Glob_CurrBasisSize0))
    allocate(Glob_FuncNum0(Glob_CurrBasisSize0))
    allocate(Glob_NonlinParam0(Glob_npt,Glob_CurrBasisSize0))
    allocate(Glob_c1(Glob_CurrBasisSize1))
    allocate(Glob_FuncNum1(Glob_CurrBasisSize1))
    allocate(Glob_NonlinParam1(Glob_npt,Glob_CurrBasisSize1))
    allocate(Glob_Index(Glob_CurrBasisSize1,2))

    IF (Glob_ProcID==0) then
      Do i=1,Glob_CurrBasisSize0
        read(1,*) Glob_FuncNum0(i),Glob_c0(i),ReadChar(1:2),Glob_NonlinParam0(1:Glob_npt,i)
      EndDo
      Do i=1,Glob_CurrBasisSize1
        read(2,*) Glob_FuncNum1(i),Glob_c1(i),ReadChar(1:2),Glob_Index(i,1),Glob_Index(i,2),Glob_NonlinParam1(1:Glob_npt,i)
      EndDo
    EndIF

    call MPI_BCAST(Glob_c0,Glob_CurrBasisSize0,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_c1,Glob_CurrBasisSize1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_NonlinParam0,Glob_npt*Glob_CurrBasisSize0,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_NonlinParam1,Glob_npt*Glob_CurrBasisSize1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_FuncNum0,Glob_CurrBasisSize0,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_FuncNum1,Glob_CurrBasisSize1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    call MPI_BCAST(Glob_Index,2*Glob_CurrBasisSize1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

  end subroutine Readwf0wf1

  subroutine DataInitForAYoungOp(YOpInput)
!Subroutine DataInitForAYoungOp initializes Young operators and on the output it gives
!Glob_NumYTerms0, Glob_YCoeff0, Glob_YMatr0, and Glob_NumYTerms1, Glob_YCoeff1, Glob_YMatr1.
!YOpInput = 0 or 1 in the cases L=0 or L=2, respectively.

    integer, intent(in)                   :: YOpInput

    integer                               :: n,npart
    integer                               :: i,j,k,p,q,t,s,w,ii,jj,kk
    character(1)                          :: c1
    integer                               :: StrLen,NumFactY
    integer                               :: TotNumOfYTerms,CurrNumOfTerms,TotNumOfYHYTerms
    integer                               :: L,R,FirstLPos,LastRPos
    integer,allocatable,dimension(:)      :: TempSymCoeff,TempSymCoeff1,NumTermsInYOpFact
    integer,allocatable,dimension(:,:,:)  :: TempSymMatr,TempSymMatr1
    integer,allocatable,dimension(:,:)    :: Matr1,Matr2,Matr3,Matr4
    integer                               :: Coeff,Cf3
    logical                               :: AreTermsIdentical
    character(Glob_YOperatorStringLength),allocatable,dimension(:) :: YOpStr,YHOpStr

    if (Glob_ProcID==0) Then
      write(*,*) ' '
      write(*,*) 'Initializing Young operator for L=',YOpInput
    EndIF
    select case (YOpInput)
    case(0)
      Glob_YOperatorString=Glob_YOperatorString0
    case(1)
      Glob_YOperatorString=Glob_YOperatorString1
    endselect

    n=Glob_n
    npart=n+1

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

    select case (YOpInput)
    case(0)
      Glob_NumFactY0=NumFactY
    case(1)
      Glob_NumFactY1=NumFactY
    endselect

!Splitting Glob_YOperatorString into an array of smaller
!strings, YOpStr. Each column of this array will contain just
!one factor, with no brackets. A '+' or a '-' sign is added in
!front of the first term in a factor if needed. Multiplication
!signs are dropped.
    allocate(YOpStr(NumFactY))
    YOpStr=' '
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
  !!!
      TotNumOfYHYTerms=TotNumOfYTerms*TotNumOfYTerms
    enddo
    if (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y operator:     ',TotNumOfYTerms
      write(*,'(1x,a,1x,i9)') 'Total number of terms in nonsimplified Y^{+}Y operator:',TotNumOfYHYTerms
    endif

    select case (YOpInput)
    case(0)
      allocate(Glob_YOpStr0(NumFactY))
      Glob_YOpStr0=YOpStr
      allocate(Glob_NumTermsInYOpFact0(NumFactY))
      Glob_NumTermsInYOpFact0=NumTermsInYOpFact
    case(1)
      allocate(Glob_YOpStr1(NumFactY))
      Glob_YOpStr1=YOpStr
      allocate(Glob_NumTermsInYOpFact1(NumFactY))
      Glob_NumTermsInYOpFact1=NumTermsInYOpFact
    endselect

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
    select case (YOpInput)
    case(0)
      Glob_NumYTerms0=Glob_NumYTerms
      allocate(Glob_YCoeff0(Glob_NumYTerms))
      allocate(Glob_YMatr0(n,n,Glob_NumYTerms))
    case(1)
      Glob_NumYTerms1=Glob_NumYTerms
      allocate(Glob_YCoeff1(Glob_NumYTerms))
      allocate(Glob_YMatr1(n,n,Glob_NumYTerms))
    endselect
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        select case (YOpInput)
        case(0)
          Glob_YCoeff0(s)=TempSymCoeff(i)
          Glob_YMatr0(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        case(1)
          Glob_YCoeff1(s)=TempSymCoeff(i)
          Glob_YMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        endselect
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
! TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff(1:Glob_NumYTerms)
! TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr(1:n,1:n,1:Glob_NumYTerms)

    select case (YOpInput)
    case(0)
      TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff0(1:Glob_NumYTerms)
      TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr0(1:n,1:n,1:Glob_NumYTerms)
    case(1)
      TempSymCoeff(1:Glob_NumYTerms)=Glob_YCoeff1(1:Glob_NumYTerms)
      TempSymMatr(1:n,1:n,1:Glob_NumYTerms)=Glob_YMatr1(1:n,1:n,1:Glob_NumYTerms)
    endselect

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
    select case (YOpInput)
    case(0)
      Glob_NumYHYTerms0=Glob_NumYHYTerms
      allocate(Glob_YHYCoeff0(Glob_NumYHYTerms0))
      allocate(Glob_YHYMatr0(n,n,Glob_NumYHYTerms0))
    case(1)
      Glob_NumYHYTerms1=Glob_NumYHYTerms
      allocate(Glob_YHYCoeff1(Glob_NumYHYTerms1))
      allocate(Glob_YHYMatr1(n,n,Glob_NumYHYTerms1))
    endselect

    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        select case (YOpInput)
        case(0)
          Glob_YHYCoeff0(s)=TempSymCoeff(i)
          Glob_YHYMatr0(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        case(1)
          Glob_YHYCoeff1(s)=TempSymCoeff(i)
          Glob_YHYMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
        endselect
      endif
    enddo
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    if (Glob_ProcID==0) then
      write(*,'(1x,a,1x,i9)') 'Total number of terms in simplified Y^{+}Y operator:   ',Glob_NumYHYTerms
    endif

    deallocate(Matr4)
    deallocate(Matr3)
    deallocate(Matr2)
    deallocate(Matr1)
    deallocate(NumTermsInYOpFact)
    deallocate(YHOpStr)
    deallocate(YOpStr)

  end subroutine DataInitForAYoungOp

  subroutine ProgramDataInit()
!Subroutine ProgramDataInit initializes some data needed for
!calculations. It should be called at the start of the program,
!right after reading input/output file.

    integer                               :: n,npart
    integer                               :: i,j,k,p,q,t,s,w
    integer                               :: pi,pj,pt,ps
    logical                               :: AreYsIdentical
    real(wp)                           :: mk,mi,m0
    integer,allocatable,dimension(:)      :: IdentParticleSet
    integer,allocatable,dimension(:,:)    :: IdentPseudoPartPairSet

    IF (Glob_ProcID==0)then
      write(*,*)''
      write(*,*) 'Initializing program data ... '
    EndIF
    n=Glob_n
    npart=n+1

!Constructing Glob_MassMatrix
    allocate(Glob_MassMatrix(n,n))
    Glob_MassMatrix(1:n,1:n)=ONEHALF/Glob_Mass(1)
    Do i=1,n
      Glob_MassMatrix(i,i)=Glob_MassMatrix(i,i)+ONEHALF/Glob_Mass(i+1)
    EndDo

!Determine the components of vector Glob_bvc
!that is used in evaluation of particle densities
    allocate(Glob_bvc(n,npart))
    Glob_MassTotal=sum(Glob_Mass(1:npart))
    Do i=1,npart
      Glob_bvc(1:n,i)=-Glob_Mass(2:n+1)/Glob_MassTotal
    EndDo
    Do i=2,npart
      Glob_bvc(i-1,i)=Glob_bvc(i-1,i)+ONE
    EndDo

!Determine the mass and the index of the lightest particle
!(reference particle excluded).
!and its index
    k=0
    mk=2*Glob_MassTotal
    Do i=1,n
      IF (Glob_Mass(i+1)<mk) then
        k=i
        mk=Glob_Mass(i+1)
      EndIF
    EndDo

    m0=Glob_Mass(1)
!alpha = sqrt( 0.5 * (m0^3 + m_k^3)/(m0*m_k*(m0 + m_k)^2) )
    Glob_dmva2 = (m0**3 + mk**3)/(TWO*m0*mk*(m0+mk)**2)
!Glob_dmvB(i,i) = (beta^2 + gamma_i^2)/(alpha^2 * M_ii) - M_ii
    Glob_dmvB(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO
    Do i=1,n
      mi=Glob_Mass(i+1)
      Glob_dmvB(i,i)=( (m0**3+mi**3)*mk*(m0+mk)**2 - (m0**3+mk**3)*mi*(m0+mi)**2 ) / ( TWO*(m0+mi)*(m0**3+mk**3)*m0*mi**2 )
    EndDo
    Glob_dmvM(1:Glob_AllowedNumOfPseudoParticles,1:Glob_AllowedNumOfPseudoParticles)=ZERO
    Glob_dmvM(1:n,1:n)=Glob_MassMatrix(1:n,1:n)
    Glob_dmvMB=Glob_dmvM+Glob_dmvB

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
    Do i=1,npart
      Do j=1,npart
        Do k=1,n
          Glob_Transposit(k,k,i,j)=1
        EndDo
      EndDo
    EndDo
!Now continue depending on type of transposition (P1i or Pij)
    Do i=2,npart
      Glob_Transposit(1:n,i-1,1,i)=-1
    EndDo
    Do i=2,npart
      Do j=i+1,npart
        Glob_Transposit(i-1,i-1,i,j)=0
        Glob_Transposit(j-1,j-1,i,j)=0
        Glob_Transposit(j-1,i-1,i,j)=1
        Glob_Transposit(i-1,j-1,i,j)=1
        Glob_Transposit(i-1,i-1,j,i)=0
        Glob_Transposit(j-1,j-1,j,i)=0
        Glob_Transposit(j-1,i-1,j,i)=1
        Glob_Transposit(i-1,j-1,j,i)=1
      EndDo
    EndDo

!Calculate Glob_YCoeff and Glob_YMatr for both case L=0 and L=1.
!The order of calling is important, it must be 0,1.
    Do i=0,1
      call DataInitForAYoungOp(i)
    EndDo

!Now we determine which particles are identical. This determination
!is based on the input values of masses and charges only. The information
!about the sets of identical particles may be needed for proper
!symmetrization of expectation values of operators that involve
!two-particle quantities (such as interparticle distances).
!The set of particles to which particle i belongs is labelled by IdentParticleSet(i)
!IF IdentParticleSet(i)=IdentParticleSet(j) then it means that
!particles i and j are identical. The largest value in IdentParticleSet
!gives the total number of identical particle sets
    allocate(IdentParticleSet(npart))
    IdentParticleSet(1)=1
    k=1
    Do i=2,npart
      s=0
      j=0
      do while ((j<i-1).and.(s==0))
        j=j+1
        IF (j>1) then
          IF ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge(j-1)==Glob_PseudoCharge(i-1))) then
            IdentParticleSet(i)=IdentParticleSet(j)
            s=1
          EndIF
        Else
          !j=1 case
          IF ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge0==Glob_PseudoCharge(i-1))) then
            IdentParticleSet(i)=IdentParticleSet(j)
            s=1
          EndIF
        EndIF
      EndDo
      IF (s==0) then
        k=k+1
        IdentParticleSet(i)=k
      EndIF
    EndDo
    Glob_NumOfIdentPartSets=maxval(IdentParticleSet(1:npart))

!Below we determine which pairs of pseudoparticles are identical.
!The information about this is stored in array IdentPseudoPartPairSet(1:n,1:n)
!Diagonal elements do not actually designate pairs of pseudoparticles but
!rather a single pseudoparticle, which corresponds to a certain pair of particles.
!IF IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(k,l) then
!it means these pairs ij and kl should be equivalent.
!The largest value of array IdentPseudoPartPairSet gives the number
!of nonequivalent pairs.
    allocate(IdentPseudoPartPairSet(1:n,1:n))
    IdentPseudoPartPairSet(1:n,1:n)=0
    k=0
    Do i=1,n
      Do j=i,n
        IF (i==j) then
          pi=1; pj=j+1
        Else
          pi=i+1; pj=j+1
        EndIF
        w=0
        do s=1,i
          IF (s==i) then
            q=j-1
          Else
            q=n
          EndIF
          Do t=s,q
            IF (w==1) cycle
            IF (s==t) then
              ps=1; pt=t+1
            Else
              ps=s+1; pt=t+1
            EndIF
            IF ((IdentParticleSet(ps)==IdentParticleSet(pi)).and. &
                (IdentParticleSet(pt)==IdentParticleSet(pj))) then
              w=1
              IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(s,t)
            EndIF
          EndDo
        EndDo
        IF (w==0) then
          k=k+1
          IdentPseudoPartPairSet(i,j)=k
        EndIF
      EndDo
    EndDo
    Glob_NumOfNoneqvPairSets=maxval(IdentPseudoPartPairSet(1:n,1:n))

!Now we create arrays Glob_NumOfPartInIdentPartSet
!and Glob_IdentPartList
    allocate(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
    allocate(Glob_IdentPartList(npart,Glob_NumOfIdentPartSets))
    Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)=0
    Glob_IdentPartList(1:npart,1:Glob_NumOfIdentPartSets)=0
    Do i=1,npart
      k=IdentParticleSet(i)
      Glob_NumOfPartInIdentPartSet(k)=Glob_NumOfPartInIdentPartSet(k)+1
      Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k),k)=i
    EndDo

!Create arrays Glob_NumOfPairsInEqvPairSet and
!Glob_EqvPairList
    allocate(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
    allocate(Glob_EqvPairList(2,n*(n+1)/2,Glob_NumOfNoneqvPairSets))
    Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets)=0
    Glob_EqvPairList(1:2,1:n*(n+1)/2,1:Glob_NumOfNoneqvPairSets)=0
    Do i=1,n
      Do j=i,n
        k=IdentPseudoPartPairSet(i,j)
        Glob_NumOfPairsInEqvPairSet(k)=Glob_NumOfPairsInEqvPairSet(k)+1
        Glob_EqvPairList(1,Glob_NumOfPairsInEqvPairSet(k),k)=i
        Glob_EqvPairList(2,Glob_NumOfPairsInEqvPairSet(k),k)=j
      EndDo
    EndDo

    deallocate(IdentPseudoPartPairSet)
    deallocate(IdentParticleSet)

  end subroutine ProgramDataInit

  subroutine ComputeSpinDep()

!local variables
    integer :: i, j, n, a, ptr, k, npt, counter
    integer :: nFactorial
    real(wp) :: Skk, temp1, temp2
    real(wp), allocatable, dimension(:, :, :) :: ketYMatrix
    real(wp), allocatable, dimension(:, :) :: SiPlusME, SiMinusME, SziME
    real(wp), allocatable, dimension(:, :, :) :: SSNCspinME
    real(wp), allocatable, dimension(:, :) :: SSNCmassChargeCoefficient
    real(wp), allocatable, dimension(:) :: parityFactor, diagS_0, diagS_1, diagS_test_0, diagS_test_1
    real(wp), allocatable, dimension(:, :) :: spinFreeME, SOspinME, spinCoeff
    real(wp) :: SSNCkl, SSNC, factor

    n = Glob_n
    npt = Glob_npt
    nFactorial = 1
    do i = 2, n
      nFactorial = nFactorial * i
    enddo

    allocate(SSNCmassChargeCoefficient(n, n))
    allocate(parityFactor(nFactorial))

    allocate(ketYMatrix(1 : n, 1 : n, nFactorial))
    allocate(spinFreeME(nFactorial, 2))
    allocate(SziME(n, nFactorial))
    allocate(SiMinusME(n, nFactorial))
    allocate(SiPlusME(n, nFactorial))
    allocate(SSNCspinME(n, n, nFactorial))
    allocate(spinCoeff(nFactorial, 2))

    call spinPreCalc(n, nFactorial, parityFactor, SSNCmassChargeCoefficient, &
                     ketYMatrix, Glob_YOperatorString0, Glob_YOperatorString1, &
                     SSNCspinME, SiMinusME, SiPlusME, SziME, spinFreeME)

    do i = 1, nFactorial
      spinCoeff(i,:) = parityFactor(i) * spinFreeME(i,:)
    enddo

! we should recalculate mean values of a unity operator here (it should be proportional to the old values)
    allocate(diagS_1(Glob_CurrBasisSize1), diagS_0(Glob_CurrBasisSize0))
    diagS_1 = ZERO
    diagS_0 = ZERO
    Skk = ZERO
    do i = 1, Glob_CurrBasisSize0
      do ptr = 1, nFactorial
        call OverlapMatrixElement_RG_0S(Glob_NonlinParam0(1:Glob_np, i), ketYMatrix(1 : n, 1 : n, ptr), Skk)
        diagS_0(i) = diagS_0(i) + spinCoeff(ptr,1) * Skk
      enddo ! Permutations from S_n
    enddo

    Skk = ZERO
    do i = 1, Glob_CurrBasisSize1
      do ptr = 1, nFactorial
        call OverlapMatrixElement_RG_2D(Glob_Index(i,1), Glob_Index(i,2), Glob_NonlinParam1(1 : npt, i), &
                                     ketYMatrix(1 : n, 1 : n, ptr), Skk)
        diagS_1(i) = diagS_1(i) + spinCoeff(ptr,2) * Skk
      enddo ! Permutations from S_n
    enddo

    SSNC = ZERO
    counter = 0
    do i = 1, Glob_CurrBasisSize1
      do j = 1, Glob_CurrBasisSize0
        counter = counter + 1
        if (mod(counter, Glob_NumOfProcs) == Glob_ProcID) then
          factor = Glob_c1(i) * Glob_c0(j) / sqrt(diagS_1(i)*diagS_0(j))
          do a = 1, nFactorial ! Permutations from S_n introduced by A operator

            call spinDependentMatrixElements(Glob_Index(i,1),Glob_Index(i,2), Glob_NonlinParam1(1 : npt, i), &
                                             Glob_NonlinParam0(1 : npt, j), ketYMatrix(1 : n, 1 : n, a), SSNCspinME(1:n, 1:n, a), &
                                             SSNCmassChargeCoefficient, SSNCkl)

            SSNC = SSNC + parityFactor(a) * factor * SSNCkl
          enddo ! Permutations from S_n
        endif ! ProcID check
      enddo !Glob_CurrBasisSize0
    enddo !Glob_CurrBasisSize1

!Combining the results of all processes

    temp1 = SSNC
    call MPI_ALLREDUCE(temp1, temp2, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    SSNC = temp2

!Printing results
    if (Glob_ProcID==0) then

      !Opening an additional file where selected expectation values will be saved
      open(2,file=Glob_ExpValFileName,status='replace')

      write(*,*) '                    SSNC=', SSNC
      write(*,*) '          (alpha^2)*SSNC=', SSNC*(Glob_FineStructConst**2)

      write(2,*) '                    SSNC=', SSNC
      write(2,*) '          (alpha^2)*SSNC=', SSNC*(Glob_FineStructConst**2)

    endif

  end subroutine ComputeSpinDep

end module workproc
