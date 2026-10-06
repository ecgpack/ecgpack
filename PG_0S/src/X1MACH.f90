!=====================================================================
!  THIRD-PARTY SOURCE
!  SLATEC machine constants I1MACH and D1MACH
!=====================================================================
!
!  These are the machine constant enquiry functions I1MACH and D1MACH
!  of the SLATEC Common Mathematical Library, obtained from Netlib:
!
!      https://www.netlib.org/slatec/
!
!  They originate in the PORT framework of Fox, Hall and Schryer at
!  Bell Laboratories:
!
!      P. A. Fox, A. D. Hall and N. L. Schryer, "Framework for a
!      portable library", ACM Trans. Math. Software 4, 2 (June 1978),
!      pp. 177-188.
!
!  The copy here is a Fortran 90 adaptation rather than the original
!  SLATEC source: the function bodies obtain the machine constants
!  from Fortran 90 intrinsic enquiry functions (RADIX, DIGITS,
!  MINEXPONENT, MAXEXPONENT, HUGE, TINY, EPSILON, BIT_SIZE) instead of
!  the per-architecture DATA tables of the original fixed-form
!  routines.  The SLATEC prologue documentation has been retained, and
!  records that adaptation in its revision history:
!
!      790101  DATE WRITTEN
!      960329  Modified for Fortran 90 (BE after suggestions by EHG)
!
!---------------------------------------------------------------------
!  License
!---------------------------------------------------------------------
!
!  The SLATEC Common Mathematical Library is in the public domain, so
!  no license conditions attach to this file.  From the SLATEC guide
!  published at https://www.netlib.org/slatec/guide:
!
!      "The Library is in the public domain and distributed by the
!       Energy Science and Technology Software Center."
!
!---------------------------------------------------------------------
!  Modifications made for ECGPACK
!---------------------------------------------------------------------
!
!  The real kind has been made selectable: DOUBLE PRECISION
!  declarations were replaced by REAL(wp) and the functions obtain the
!  kind parameter wp from module wp_def, so that the constants
!  returned describe the working precision actually in use, whether
!  double (fp64), extended (fp80) or quadruple (fp128).
!
!  See THIRD-PARTY-NOTICES.md in the root of the ECGPACK repository.
!=====================================================================
!Fortran-90 versions of routines I1MACH and D1MACH
!where double precision is replaced with the quadruple one.

!DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
!get the definition of wp parameter from module wp_def
        use wp_def
        IMPLICIT NONE
        INTEGER :: I
        REAL :: X
        REAL(wp) :: XX
!***BEGIN PROLOGUE  I1MACH
!***PURPOSE  Return integer machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0 .LE. X(I) .LT. B for I=1,...,T,
!                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by EHG).
!   980727  Modified value of I1MACH(6) (BE after suggestion by EHG).
!***END PROLOGUE  I1MACH
!
        X  = 1.0
        XX = 1.0_wp

        SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          I1MACH = 0 ! Punch unit is no longer used
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = 4            ! Characters per integer is hopefully no
          ! longer used.
          ! If it is used it has to be set manually.
          ! The value 4 is correct on IEEE-machines.
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X)
        CASE (11)
          I1MACH = DIGITS(X)
        CASE (12)
          I1MACH = MINEXPONENT(X)
        CASE (13)
          I1MACH = MAXEXPONENT(X)
        CASE (14)
          I1MACH = DIGITS(XX)
        CASE (15)
          I1MACH = MINEXPONENT(XX)
        CASE (16)
          I1MACH = MAXEXPONENT(XX)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
9000      FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
          STOP
        END SELECT
        RETURN
      END

!DECK D1MACH
      FUNCTION D1MACH (I)
!get the definition of wp parameter from module wp_def
        use wp_def
        IMPLICIT NONE
        real(wp) :: D1MACH
        INTEGER :: I
        REAL(wp) :: B, X
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)
!***END PROLOGUE  D1MACH
!
        X = 1.0_wp
        B = RADIX(X)
        SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)               ! the smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          D1MACH = EPSILON(X)/B          ! the smallest relative spacing.
        CASE (4)
          D1MACH = EPSILON(X)            ! the largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
9000      FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
          STOP
        END SELECT
        RETURN
      END
