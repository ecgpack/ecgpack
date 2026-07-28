program main

  use workproc
  implicit none

!Local variables
  character(len=100) :: line
  character(len=20)  :: mode, transition
  integer :: ios
  logical :: spin_calc, scalar_calc
  spin_calc = .false.
  scalar_calc = .false.

  !Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

  if (Glob_ProcID==0) then
    write (*,*) 'Program RG_1P-1P started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  open(unit=10, file='inp.txt', status='old', action='read')
  read(10,'(A)',iostat=ios) line
  if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop 'Error reading inp.txt'

  ! Reads first two whitespace-separated strings
  read(line,*,iostat=ios) mode, transition
  if (ios /= 0) call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop 'Cannot parse input line'

  select case (trim(mode))
  case ('SPIN')
      spin_calc = .true.
  case ('SCALAR')
      scalar_calc = .true.
  case ('ALL')
      spin_calc = .true.
      scalar_calc = .true.
  case default
      if (Glob_ProcID==0) then
        write(*,*) 'Unknown mode of calculation: ', trim(mode)
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
  end select

  select case (trim(transition))
  case ('3P_1P')
      Glob_selectTransition = 2
  case ('XP_XP')
      Glob_selectTransition = 1
  case default
      if (Glob_ProcID==0) then
        write(*,*) 'Unknown transition: ', trim(transition)
      endif
      call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
  end select

  if (Glob_selectTransition == 2 .and. scalar_calc) then
    if (Glob_ProcID==0) then
      write(*,*) 'Scalar calculation for XP -> XP transition is not supported'
    endif
    call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !stop
  end if

  call Readwf0wf1()
  call ProgramDataInit()
  if (Glob_ProcID==0) then
    write(*,*) ' '
    write(*,*)'Young string L=0: ',Glob_YOperatorString0
    write(*,*)'Young string L=1: ',Glob_YOperatorString1
  endif


  if (Glob_selectTransition == 1) then
    if (spin_calc) then
      if (Glob_ProcID==0) then
         write(*,*) 'Calculating XP -> XP transition'
         write(*,*) 'Calculating spin-dependent matrix elements'
      endif
      call ComputeSpinDep()
    endif

    if (scalar_calc) then
      if (Glob_ProcID==0) then
        write(*,*) 'Calculating XP -> XP transition'
        write(*,*) 'Calculating scalar matrix elements'
      endif
      call ComputeScalar()
    endif
  else if (Glob_selectTransition == 2) then
    if (Glob_ProcID==0) then
      write(*,*) 'Calculating XP -> XP transition'
      write(*,*) 'Calculating spin-dependent matrix elements'
    endif
    call ComputeSpinDep()
  endif

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main

