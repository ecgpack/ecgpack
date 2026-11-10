! Compare the working-kind DBESK implementation against independent mpmath
! reference values for K_0(x), K_1(x), EXP(x)*K_0(x), and EXP(x)*K_1(x).
!
! The routines tested are those of the ECGPACK module bessel
! (RG_0S/src/bessel.f90 and bessel_coeff.f90, identical in RG_1P, RG_2D, and
! RG_2P); the working kind is the one selected by RG_0S/src/wp_def_<PREC>.f90.
!
! By default the program reads bessel_k01_mpmath_reference.dat from the
! current directory.  A different path may be supplied as the first command-
! line argument.
program test_bessel_k01_wp
  use globvars, only: wp
  use bessel, only: DBESK
  implicit none

  integer :: Failures,IOStat,KODE,NZ,ReferenceUnit,Tests
  character(512) :: InputLine,ReferenceFile
  real(wp) :: Computed(4),Reference(4),Tolerance,X,Y(2)

  ReferenceFile = 'bessel_k01_mpmath_reference.dat'
  if (command_argument_count().ge.1) then
    call get_command_argument(1,ReferenceFile)
  endif

  open(newunit=ReferenceUnit,file=trim(ReferenceFile),status='old', &
       action='read',iostat=IOStat)
  if (IOStat.ne.0) then
    error stop 'Unable to open the mpmath Bessel reference-data file'
  endif

  ! This allows for accumulated roundoff in the series and recurrence paths,
  ! while still requiring accuracy appropriate to the selected working kind.
  Tolerance = 1000.0_wp*epsilon(1.0_wp)
  Failures = 0
  Tests = 0

  write(*,'(A)') 'Modified Bessel K_0/K_1 comparison with mpmath'
  write(*,'(A,I0)') 'Working-kind decimal precision: ',precision(1.0_wp)
  write(*,'(A,ES14.6E3)') 'Relative-error tolerance:   ',Tolerance
  write(*,'(/,A)') '       x        quantity            computed'// &
                   '                              reference'// &
                   '                             relative error  result'

  do
    read(ReferenceUnit,'(A)',iostat=IOStat) InputLine
    if (IOStat.lt.0) exit
    if (IOStat.gt.0) error stop 'Error while reading Bessel reference data'
    InputLine = adjustl(InputLine)
    if (len_trim(InputLine).eq.0) cycle
    if (InputLine(1:1).eq.'#') cycle

    read(InputLine,*,iostat=IOStat) X,Reference
    if (IOStat.ne.0) error stop 'Malformed line in Bessel reference data'

    KODE = 1
    call DBESK(X,0.0_wp,KODE,2,Y,NZ)
    Computed(1:2) = Y

    KODE = 2
    call DBESK(X,0.0_wp,KODE,2,Y,NZ)
    if (NZ.ne.0) error stop 'Unexpected underflow in scaled DBESK test'
    Computed(3:4) = Y

    call CheckValue('K_0(x)       ',X,Computed(1),Reference(1), &
                    Tolerance,Tests,Failures)
    call CheckValue('K_1(x)       ',X,Computed(2),Reference(2), &
                    Tolerance,Tests,Failures)
    call CheckValue('exp(x)*K_0(x)',X,Computed(3),Reference(3), &
                    Tolerance,Tests,Failures)
    call CheckValue('exp(x)*K_1(x)',X,Computed(4),Reference(4), &
                    Tolerance,Tests,Failures)
  enddo
  close(ReferenceUnit)

  write(*,'(/,A,I0,A,I0,A)') 'Summary: ',Tests-Failures,' passed, ', &
                              Failures,' failed'
  if (Failures.ne.0) error stop 'Modified Bessel K_0/K_1 test failed'

contains

  subroutine CheckValue(Name,Argument,Value,Expected,AllowedError, &
                        NumberOfTests,NumberOfFailures)
    character(*), intent(in) :: Name
    real(wp), intent(in) :: AllowedError,Argument,Expected,Value
    integer, intent(inout) :: NumberOfFailures,NumberOfTests

    character(4) :: Result
    real(wp) :: RelativeError

    ! Very large X can put the unscaled function below the working kind's
    ! exponent range.  In that case a zero or subnormal result is correct;
    ! the corresponding scaled-function test still measures relative error.
    if (abs(Expected).lt.tiny(1.0_wp)) then
      if (abs(Value).lt.tiny(1.0_wp)) then
        RelativeError = 0.0_wp
      else
        RelativeError = huge(1.0_wp)
      endif
    else
      RelativeError = abs((Value-Expected)/Expected)
    endif

    NumberOfTests = NumberOfTests+1
    if (RelativeError.le.AllowedError) then
      Result = 'PASS'
    else
      Result = 'FAIL'
      NumberOfFailures = NumberOfFailures+1
    endif

    write(*,'(ES11.3E3,2X,A15,2X,ES44.35E4,2X,ES44.35E4,2X,'// &
             'ES13.5E3,2X,A4)') Argument,Name,Value,Expected, &
                                RelativeError,Result
  end subroutine CheckValue

end program test_bessel_k01_wp
