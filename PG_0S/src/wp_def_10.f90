MODULE wp_def
  ! This module contains basic constants and basic input/output subroutines
  ! which depend on the choice of precision, compiler, and MPI implementation

  ! include 'mpif.h'
  USE mpi

  ! Without this, an undeclared name in the write routines below would be
  ! implicitly typed - DEFAULT REAL, i.e. single precision - in the one module
  ! whose entire job is to define the working precision.
  IMPLICIT NONE

  ! Define kind parameter for real type
  INTEGER, PARAMETER :: wp = 10

  ! This is data type identifier for MPI corresponding to real type of kind wp
  INTEGER, PARAMETER :: MPI_WP = MPI_REAL16

  ! This is the number of particles in the system that should be set by the user.
  ! For reasons related to performance, it is made a fixed (compile time) parameter.
  ! ReadIOFile requires the data file to specify EXACTLY this many particles
  ! (PARTICLES record), so it is an exact count, not an upper bound. BH (bh/ECG5): 8;
  INTEGER, PARAMETER :: Glob_AllowedNumOfParticles=8

CONTAINS

  ! Subroutines writereal/writerealadv (one value) and writerealarr/
  ! writerealarradv (arrays) write REAL(wp) values. All real output of the
  ! program goes through them, so the edit descriptor here is the only
  ! place where the output width is fixed. ES29.20 gives 21 significant
  ! digits, an exact decimal round trip for the 64-bit mantissa of x87
  ! extended precision. The exponent width is chosen per value: without an
  ! explicit width Fortran drops the 'E' of a three-digit exponent, and at
  ! this precision exponents reach four digits, so |r| >= 1e99 or
  ! |r| < 1e-98 uses ES29.20E4. Both forms are 29 columns wide.

  FUNCTION IsWideExponent(r)
    ! .TRUE. when r needs more than two exponent digits, so that writereal
    ! below keeps the 'E' instead of letting Fortran drop it. The bounds
    ! carry a decade of margin: a value just under 1e100 can be ROUNDED UP
    ! to 1.0...E+100 by the formatting itself and would then lose its 'E'
    ! on the two-digit branch.
    LOGICAL  :: IsWideExponent
    REAL(wp) :: r
    IsWideExponent = (r /= 0.0_wp) .AND. &
                     ((ABS(r) >= 1.0e99_wp) .OR. (ABS(r) < 1.0e-98_wp))
  END FUNCTION IsWideExponent

  SUBROUTINE writereal(u, r)
    INTEGER  :: u  ! i/o unit
    REAL(wp) :: r  ! real number that needs to be written
    IF (IsWideExponent(r)) THEN
      WRITE(u, '(1x,ES29.20E4)', ADVANCE='no') r
    ELSE
      WRITE(u, '(1x,ES29.20)', ADVANCE='no') r
    ENDIF
  END SUBROUTINE writereal

  SUBROUTINE writerealadv(u, r)
    INTEGER  :: u  ! i/o unit
    REAL(wp) :: r  ! real number that needs to be written
    IF (IsWideExponent(r)) THEN
      WRITE(u, '(1x,ES29.20E4)') r
    ELSE
      WRITE(u, '(1x,ES29.20)') r
    ENDIF
  END SUBROUTINE writerealadv

  SUBROUTINE writerealarr(u, r, k)
    INTEGER  :: u     ! i/o unit
    INTEGER  :: k     ! the number of elements to write (writing begins with element 1)
    REAL(wp) :: r(k)  ! real array that needs to be written
    ! k is declared BEFORE r(k) on purpose. A variable used in a specification
    ! expression must already have its type established at that point (F2018
    ! 10.1.11), and with IMPLICIT NONE in force a later "integer k" does not
    ! establish it. gfortran accepts the reverse order; stricter front ends
    ! need not.
    INTEGER :: i
    DO i = 1, k
      CALL writereal(u, r(i))
    ENDDO
  END SUBROUTINE writerealarr

  SUBROUTINE writerealarradv(u, r, k)
    INTEGER  :: u     ! i/o unit
    INTEGER  :: k     ! the number of elements to write (writing begins with element 1)
    REAL(wp) :: r(k)  ! real array that needs to be written
    ! k is declared BEFORE r(k) on purpose. A variable used in a specification
    ! expression must already have its type established at that point (F2018
    ! 10.1.11), and with IMPLICIT NONE in force a later "integer k" does not
    ! establish it. gfortran accepts the reverse order; stricter front ends
    ! need not.
    INTEGER :: i
    DO i = 1, k-1
      CALL writereal(u, r(i))
    ENDDO
    CALL writerealadv(u, r(k))
  END SUBROUTINE writerealarradv

  ! Subroutines writestring and writestringadv realize nonadvanced (advanced)
  ! output of strings

  SUBROUTINE writestring(u, s, k)
    INTEGER      :: u  ! i/o unit
    CHARACTER(*) :: s  ! string that needs to be written
    INTEGER      :: k  ! the length of the string (k first characters) to write.
    INTEGER      :: i
    WRITE(u, '(1x)', ADVANCE='no')
    DO i = 1, k
      WRITE(u, '(a1)', ADVANCE='no') s(i:i)
    ENDDO
  END SUBROUTINE writestring

  SUBROUTINE writestringadv(u, s, k)
    INTEGER      :: u  ! i/o unit
    CHARACTER(*) :: s  ! string that needs to be written
    INTEGER      :: k  ! the length of the string (k first characters) to write.
    INTEGER      :: i
    WRITE(u, '(1x)', ADVANCE='no')
    DO i = 1, k-1
      WRITE(u, '(a1)', ADVANCE='no') s(i:i)
    ENDDO
    WRITE(u, '(a1)') s(k:k)
  END SUBROUTINE writestringadv

END MODULE wp_def
