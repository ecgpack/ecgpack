
program main
  use workproc
  implicit none

!Local variables
  integer        :: i,j
  character(70)  :: ReadChar
  real(wp)    :: LineStrength
  real(wp)    :: OscillatorStrength
  real(wp)    :: DeltaE

!Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

  if (Glob_ProcID==0) then
    write (*,*) 'Program RG_1P-2P started'
    write (*,'(1x,a,1x,i0)') 'Number of parallel MPI processes running:',Glob_NumOfProcs
    write(*,*)
    write(*,'(2X,A)')'Transition under consideration:'
    write(*,'(4X,A)')'P^o, y component  -->  P^e, z component'
    write(*,'(2X,A)')'Calculating the x component of the transition dipole moment integral'
    write(*,*)
  endif

  call READwf1wf2()
  call ProgramDataInit()

  DeltaE = Glob_CurrEnergy1-Glob_CurrEnergy2
  call ComputeTranDipoL1L2()

  if (Glob_ProcID==0) then
    write(*,*)
    write(*,*)
    write(*,'(A)') '   Energy difference:     '
    write(*,'(A)') '  ------------------------'
    write(*,*)
    write(*,'(4X,A)',advance='no') '    E(P) - E(D)  = '
    call writereal(6,abs(DeltaE))
    write(*,'(A)') '   Hartree'

    write(*,*)
    write(*,*)
    write(*,'(A)') '   Length Gauge:          '
    write(*,'(A)') '  ------------------------'
    write(*,*)
    write(*,'(4X,A)',advance='no') ' < Po | x | Pe >   = '
    call writerealadv(6,Glob_ExpVals1)
    write(*,'(4X,A)') ' < Po | y | Pe >   =   0.0'
    write(*,'(4X,A)') ' < Po | z | Pe >   =   0.0'

    write(*,*)
    write(*,*)
    write(*,'(A)') '   Velocity Gauge:        '
    write(*,'(A)') '  ------------------------'
    write(*,*)
    write(*,'(4X,A)',advance='no') '< Po | P(x) | Pe >  ='
    call writerealadv(6,Glob_ExpVals2)
    write(*,'(4X,A)') '< Po | P(y) | Pe >  =  0.0'
    write(*,'(4X,A)') '< Po | P(z) | Pe >  =  0.0'

    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,'(A)') ' Program completed successfully.      '
    write(*,*)

  endif

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main
