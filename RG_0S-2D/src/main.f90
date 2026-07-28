program main

  use workproc
  implicit none

!Local variables
  integer        :: i,j
  character(70)  :: ReadChar
  real(wp)    :: LineStrength
  real(wp)    :: OscillatorStrength
  real(wp)    :: DeltaE
  real(wp)    :: DeltaE_cm

!Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

  if (Glob_ProcID==0) then
    write (*,*) 'Program RG_0S-2D started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  call Readwf0wf1()
  call ProgramDataInit()
  if (Glob_ProcID==0) then
    write(*,*) ' '
    write(*,*)'Young string L=0: ',Glob_YOperatorString0
    write(*,*)'Young string L=1: ',Glob_YOperatorString1
  endif
  call ComputeSpinDep()

  DeltaE = Glob_CurrEnergy0-Glob_CurrEnergy1

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main

