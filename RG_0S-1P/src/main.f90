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
    write (*,*) 'Program RG_0S-1P started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write (*,*)
  endif

  call Readwf0wf1()
  call ProgramDataInit()

  if (Glob_ProcID==0) then
    write(*,*) ' '
    write(*,*)'Symmetry L=0: ',Glob_YOperatorString0
    write(*,*)'Symmetry L=1: ',Glob_YOperatorString1
  endif

  DeltaE = Glob_CurrEnergy0-Glob_CurrEnergy1

  call ComputeTranDipoL0L1()

  if (Glob_ProcID==0) then
    write(*,*) '  '
    write(*,*) '  '

    write(*,'(a30)',advance='no') ' E{initial} - E{final}  =     '
    call writereal(6,abs(DeltaE))
    write(*,'(a10)') '   Hartree'

    write(*,*) '  '
    write(*,*) '  '
    write(*,*) 'Transition Dipole Moment Integral in Length Form:  '
    write(*,*) '  '
    write(*,'(a35)') ' <L=0|x|L=1>    =       0.0'
    write(*,'(a35)') ' <L=0|y|L=1>    =       0.0'
    write(*,'(a30)',advance='no') ' <L=0|z|L=1>    =     '
    call writerealadv(6,Glob_ExpVals1)

    write(*,*) '  '
    write(*,*) '  '
    write(*,*) 'Transition Dipole Moment Integral in Velocity Form:  '
    write(*,*) '  '
    write(*,'(a30)') ' <L=0|P(x)|L=1>    =     '
    write(*,'(a30)') ' <L=0|P(y)|L=1>    =     '
    write(*,'(a30)',advance='no') ' <L=0|P(z)|L=1>    = -i *'
    call writerealadv(6,Glob_ExpVals2)

    write(*,*) '  '
    write(*,*) '  '
    write(*,*) 'Data Reader and Matrix Calculator Program is completed'
    write(*,*) 'Program has stopped'
    write(*,*) ' '
  endif

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main

