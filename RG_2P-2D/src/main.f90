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
    write (*,*) 'Program RG_2P-2D started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  open(unit=10, file='inp.txt', status='old', action='read')
  read(10,'(A)',iostat=ios) line
  if (ios /= 0) stop 'Error reading inp.txt'

  ! Reads first two whitespace-separated strings
  read(line,*,iostat=ios) transition
  if (ios /= 0) stop 'Cannot parse input line'

  !selectTransition = 1 -- calculate 3P_2 -> 1D_2 matelem  <3P | H_SO | 1D>
  !selectTransition = 2 -- calculate 3P_1 -> 3D_1 matelem <3P | H_SO + H_SSNC | 3D>
  !selectTransition = 3 -- calculate 1P_1 -> 3D_1 matelem <1P | H_SO | 3D>
  !selectTransition = 4 -- calculate 4P_5/2 -> 2D_5/2 matelem <4P | H_SO + H_SSNC | 2D>

  select case (trim(transition))
  case ('3P_1D')
      Glob_selectTransition = 1
  case ('3P_3D')
      Glob_selectTransition = 2
  case ('1P_3D')
      Glob_selectTransition = 3
  case ('4P_2D')
      Glob_selectTransition = 4
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

  call ComputeSpinDep()




  call MPI_FINALIZE(Glob_MPIErrCode)


end program main

