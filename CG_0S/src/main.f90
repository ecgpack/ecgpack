program main

  use workproc
  implicit none

!Local variables
  integer      i,iw,Kstart,Kstop,Kstep,OpenFileErr,OptimizationType

!Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

  if (Glob_ProcID==0) then
    write (*,'(1x,a,1x,a,1x,a)') 'Program',Glob_BasisType,'started'
    write (*,'(1x,a,1x,i0)') 'Working precision (wp):',wp
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  call ReadIOFile()
  if (Glob_IsOptCycleScripted) call ReadBlackList()
  call ProgramDataInit()

!Seed the random number generators. The base seed is selected by the optional
!RANDOM_SEED line of the data file (Glob_RndSeedInput, set by ReadIOFile): the
!integer given there, the integer held by the environment variable ECG_RND_SEED,
!or, with SYSTEM_CLOCK (the default), the system clock. The clock count is
!combined with the time of day from date_and_time because some compilers (e.g.
!nvfortran) count system_clock from the start of the program, which alone would
!give the same seed in every run. If ECG_RND_SEED is not set or does not hold an
!integer, the clock is used instead.
!Rank 0 determines the base seed, prints it, and broadcasts it, so that a run can
!be repeated by setting RANDOM_SEED to the printed value (with the same number of
!MPI processes). Every rank offsets the base seed by its rank,
!so that the ranks drawing their own random numbers get different sequences. The
!rank seed is computed with 64-bit integers and reduced modulo huge(0) so that
!large base seeds cannot overflow. It initializes both of the code's own random
!number generators, urand (uniform) and drnor (normal), whose sequences do not
!depend on the compiler (the intrinsic random_number is not used). drnor_start
!needs a nonzero seed, hence the substitutions of 519 for 0.
  block
    integer,parameter :: ik8=selected_int_kind(18)
    integer       RndSeedBase,RndSeedErr,RndSeedStat,RndTime(8),k
    integer(ik8)  RndClock
    character(32) RndSeedEnv,RndSeedSource
    real(8)       RndSeedEcho

    if (Glob_ProcID==0) then
      RndSeedErr=1
      if (Glob_RndSeedInput=='ECG_RND_SEED') then
        RndSeedSource='ECG_RND_SEED'
        call get_environment_variable('ECG_RND_SEED',RndSeedEnv,status=RndSeedStat)
        if (RndSeedStat==0) read(RndSeedEnv,*,iostat=RndSeedErr) RndSeedBase
        if (RndSeedErr/=0) then
          write(*,*) 'Warning: environment variable ECG_RND_SEED is not set or not an integer.'
          write(*,*) 'The random seed is taken from the system clock instead.'
        endif
      elseif (Glob_RndSeedInput/='SYSTEM_CLOCK') then
        RndSeedSource='data file'
        read(Glob_RndSeedInput,*,iostat=RndSeedErr) RndSeedBase
      endif
      if (RndSeedErr/=0) then
        RndSeedSource='SYSTEM_CLOCK'
        call system_clock(count=RndClock)
        call date_and_time(values=RndTime)
        RndClock=RndClock+RndTime(8)+1000_ik8*(RndTime(7)+60_ik8*(RndTime(6)+60_ik8*(RndTime(5)+24_ik8*RndTime(3))))
        RndSeedBase=int(mod(RndClock,int(huge(0),ik8)))
      endif
      if (RndSeedBase==0) RndSeedBase=519
      write(*,'(1x,a,1x,i0,1x,a,a,a)') 'Random seed:',RndSeedBase,'(from ',trim(RndSeedSource),')'
      write(*,*)
    endif
    call MPI_BCAST(RndSeedBase,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)

    k=int(mod(int(RndSeedBase,ik8)+7919_ik8*Glob_ProcID,int(huge(0),ik8)))
    if (k==0) k=519
    call urand_start(k)
    RndSeedEcho=drnor_start(k)
  end block

!Empty swap file as it may contain some garbage left
!after last run (if there was a failure)
  if ((Glob_ProcID==0).and.Glob_UseSwapFile) then
    open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
    if (OpenFileErr==0) then
      write(1) 'Swap file is empty'
      close(1)
    endif
  endif

  OptimizationType=1

  do i=1,Glob_NumOfBBOPSteps
    call linalg_setparam(Glob_CurrBasisSize)

    Glob_CurrBBOPStep=i
    Glob_ApproxEnergy=Glob_CurrEnergy*Glob_InvItParameter

    select case (Glob_BBOP(i)%Action)

    case('BASIS_ENL')
      Kstart=Glob_BBOP(i)%A
      Kstop=Glob_BBOP(i)%B
      Kstep=Glob_BBOP(i)%C
      if (Kstop>Glob_CurrBasisSize) then
        if (Kstart<=Glob_CurrBasisSize+1) then
          Kstart=Glob_CurrBasisSize+1
          select case (Glob_BBOP(i)%GSEPSolutionMethod)
          case('G')
            call BasisEnlG(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R)
          case('I')
            call BasisEnlI(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R)
          endselect
        else
          if (Glob_ProcID==0) then
            write(*,'(1x,a,1x,i0)') 'Error EC0001 in main: incorrect BBOP step',i
            write(*,*) 'One or more parameters in BASIS_ENL are incorrect'
          endif
        endif
      endif

    case('OPT_CYCLE')
      if (Glob_ProcID==0) then
        if ((Glob_BBOP(i)%C>Glob_BBOP(i)%A).or.(Glob_BBOP(i)%B>Glob_BBOP(i)%C)) then
          write(*,'(1x,a,1x,i0)') 'Error EC0002 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in OPT_CYCLE are incorrect'
        endif
      endif
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
          (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
        if (Glob_History(Glob_CurrBasisSize)%CyclesDone<Glob_BBOP(i)%F) then
          select case (Glob_BBOP(i)%GSEPSolutionMethod)
          case('G')
            call OptCycleG(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          case('I')
            call OptCycleI(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          endselect
        endif
      endif

    case('FULL_OPT1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
          (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call FullOpt1G(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q, &
                         Glob_BBOP(i)%R,real(Glob_BBOP(i)%E,wp),real(Glob_BBOP(i)%F,wp), &
                         Glob_BBOP(i)%FileName1)
        case('I')
          call FullOpt1I(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q, &
                         Glob_BBOP(i)%R,real(Glob_BBOP(i)%E,wp),real(Glob_BBOP(i)%F,wp), &
                         Glob_BBOP(i)%FileName1)
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0003 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in FULL_OPT1 are incorrect'
        endif
      endif

    case('ELIM_LCFN')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call EliminateLittleContribFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1, &
                                          Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LCFN'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0004 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in ELIM_LCFN are incorrect'
        endif
      endif

    case('ELIM_LND1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call EliminateLinDepFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0005 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in ELIM_LND1 are incorrect'
        endif
      endif

    case('SEPR_LND1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call SeparateLinDepFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1, &
                                  Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0006 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in SEPR_LND1 are incorrect'
        endif
      endif

    case('SEPR_FLCF')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call SeparateFuncLargeCoeff(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1, &
                                      Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0007 in main: incorrect BBOP step',i
          write(*,*) 'One or more parameters in SEPR_FLCF are incorrect'
        endif
      endif

    case('EXPC_VALS')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call ExpectationValuesG(1,Glob_FileNameNone,Glob_FileNameNone,Glob_FileNameNone, &
                                  Glob_FileNameNone)
        case('I')
          if (Glob_ProcID==0) write(*,*) &
            'Sorry, GSEP soluton method I has not been implemented in EXPC_VALS'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0008 in main: incorrect BBOP step',i
          write(*,*) 'Second parameter in EXPC_VALS is incorrect'
        endif
      endif

    case('DENSITIES')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call ExpectationValuesG(1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
                                  Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4)
        case('I')
          if (Glob_ProcID==0) write(*,*) &
            'Sorry, GSEP soluton method I has not been implemented in DENSITIES'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,'(1x,a,1x,i0)') 'Error EC0009 in main: incorrect BBOP step',i
          write(*,*) 'Second parameters in DENSITIES is incorrect'
        endif
      endif

    case('SAVE_FILE')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        if (Glob_ProcID==0) then
          iw=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(*,*)
          write(*,*) 'Saving basis in file ',Glob_BBOP(i)%FileName1(1:iw),'...'
          call SaveResults(Filename=Glob_BBOP(i)%FileName1,Sort='no')
          write(*,*) ' done'
          write(*,*)
        endif
      endif

    endselect
  enddo

  if (Glob_ProcID==0) then
    if (Glob_UseSwapFile) then
      !Empty swap file as it may contain the data saved after
      !the last basis building and optimization program step
      open(1,file=Glob_SwapFileName,form='binary',status='replace',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        write(1) 'Swap file is empty'
        close(1)
        write(*,*)
        write(*,*) 'Swap file has been cleaned up'
      endif
    endif
    write(*,*) ' '
    write(*,*) 'Basis Building and Optimization Program is completed'
    write(*,*) 'Program stopped'
  endif

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main
