program main

  use workproc
  implicit none
  !Local variables
  character(len=100) :: line
  character(len=20)  :: transition
  integer :: ios

  !Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

  if (Glob_ProcID==0) then
    write (*,*) 'Program RG_0S-2P started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  open(unit=10, file='inp.txt', status='old', action='read')
  read(10,'(A)',iostat=ios) line
  if (ios /= 0) stop 'Error reading inp.txt'

  ! Reads first two whitespace-separated strings
  read(line,*,iostat=ios) transition
  if (ios /= 0) stop 'Cannot parse input line'

  select case (trim(transition))
  case ('3P_1S')
      Glob_selectTransition = 1
  case ('3P_3S')
      Glob_selectTransition = 2
  case ('4P_2S')
      Glob_selectTransition = 3
  case default
      if (Glob_ProcID==0) then
        write(*,*) 'Unknown transition: ', trim(transition)
      endif
      stop
  end select

  call Readwf0wf1()
  call ProgramDataInit()
  if (Glob_ProcID==0) then
    write(*,*) ' '
    write(*,*)'Young string L=0: ',Glob_YOperatorString0
    write(*,*)'Young string L=1: ',Glob_YOperatorString1
  endif

  if (Glob_ProcId==0) then 
    if (Glob_selectTransition == 1)  write(*,*) 'Calculating 3P -> 1S transition'
    if (Glob_selectTransition == 2)  write(*,*) 'Calculating 3P -> 3S transition'
    if (Glob_selectTransition == 3)  write(*,*) 'Calculating 4P -> 2S transition'
  endif
        
  call ComputeSpinDep()

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main

