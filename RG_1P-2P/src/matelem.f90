module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 and L=2 Gaussians.

  use globvars
  implicit none

contains

  SUBROUTINE OverLapElement_S_Po(m_k, vechLk, P, Skk)
!==============================================================================================
! This subroutine calculates the diagonal overlap for L=1 correlated Gaussians:
!
! fk = z{m_k} * exp[-r'(Lk*Lk')r]
!
!                  <P*fk|P*fk>               ||L_kk||^3       v_k' * inv_tAkk * tv_k
! <psi_k|psi_k> = -------------- = 2^(3n/2) ------------- * --------------------------
!                    <fk|fk>                 |tAkk|^{3/2}      v_k' * inv_Akk * v_k
!
! where tAkk = Ak + P'*Ak*P,  tv_k = P'*v_k
!
!----------------------------------------------------------------------------------------------
! Input:
!   m_k      :: integer that determines which z-component is in the premultiplier
!   vechLk   :: Array of length n(n+1)/2 of exponential parameters
!   P        :: The symmetry permutation matrix of size n x n
!
! Output:
!   Skk      :: <psi_k|psi_k>
!==============================================================================================

    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! ARGUMENTS
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN)      :: m_k
    REAL(wp), INTENT(IN)  :: vechLk(Glob_np)
    REAL(wp), INTENT(IN)  :: P(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: Skk

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !---------------------------------------------------------------------------
    INTEGER                  :: n
    INTEGER                  :: i, j, k, indx

    REAL(wp)              :: Lk(Glob_n, Glob_n)
    REAL(wp)              :: tAk(Glob_n, Glob_n), tAkk(Glob_n, Glob_n)
    REAL(wp)              :: inv_tAkk(Glob_n, Glob_n), inv_Akk(Glob_n, Glob_n)
    REAL(wp)              :: W1(Glob_n, Glob_n)
    REAL(wp)              :: tvk(Glob_n)
    REAL(wp)              :: temp1
    REAL(wp)              :: det_Lk, det_tAkk, tau3

    !---------------------------------------------------------------------------
    ! 1. Build lower-triangular Lk from vechLk
    !---------------------------------------------------------------------------
    indx = 0
    DO i = 1, Glob_n
      DO j = i, Glob_n
        indx = indx + 1
        Lk(i, j) = ZERO
        Lk(j, i) = vechLk(indx)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 2. Compute inv_Akk using Cholesky factor Lk
    !    inv_Akk = (1/2) * Lk^{-T} * Lk^{-1} = (2*Ak)^{-1} = Akk^{-1}
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    DO i = 1, Glob_n
      W1(i,i) = ONE / Lk(i,i)
      DO j = i + 1, Glob_n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - Lk(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / Lk(j,j)
      END DO
    END DO

    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_Akk(i,j) = ONEHALF * temp1
        inv_Akk(j,i) = ONEHALF * temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 3. Compute Ak = Lk * Lk'
    !---------------------------------------------------------------------------
    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = 1, i
          temp1 = temp1 + Lk(i, k) * Lk(j,k)
        END DO
        tAk(i,j) = temp1
        tAk(j,i) = temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 4. Apply symmetry permutation and form tAkk = Ak + P' * Ak * P
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    ! W1 = P' * Ak
    DO i = 1, Glob_n
      DO j = 1, Glob_n
        temp1 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + P(k,j) * tAk(k,i)
        END DO
        W1(j,i) = temp1
      END DO
    END DO

    ! tAkk = Ak + P' * Ak * P
    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + W1(i, k) * P(k,j)
        END DO
        tAkk(i,j) = tAk(i,j) + temp1
        tAkk(j,i) = tAkk(i,j)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 5. Determinant of Lk (product of diagonal elements)
    !---------------------------------------------------------------------------
    det_Lk = ONE
    DO i = 1, Glob_n
      det_Lk = det_Lk * Lk(i,i)
    END DO

    !---------------------------------------------------------------------------
    ! 6. Cholesky factorization of tAkk -> W1, and determinant
    !---------------------------------------------------------------------------
    det_tAkk = ONE

    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = tAkk(i,j)
        DO k = i - 1, 1, -1
          temp1 = temp1 - W1(i, k) * W1(j,k)
        END DO

        IF (i == j) THEN
          W1(i,i) = SQRT(temp1)

          IF (temp1 <= 0.0_wp) THEN
            WRITE(*,*) 'ERROR: tAkk not SPD at i=',i
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
          END IF

          det_tAkk = det_tAkk * temp1
        ELSE
          W1(j,i) = temp1 / W1(i,i)
          W1(i,j) = ZERO
        END IF

      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 7. Invert tAkk using its Cholesky factor (stored in W1)
    !---------------------------------------------------------------------------
    DO i = 1, Glob_n
      W1(i,i) = ONE / W1(i,i)
      DO j = i + 1, Glob_n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - W1(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / W1(j,j)
      END DO
    END DO

    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_tAkk(i,j) = temp1
        inv_tAkk(j,i) = temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 8. Compute tau3 = v_k' * inv_tAkk * tv_k
    !    where v_k = e_{m_k} (unit vector) and tv_k = P(m_k, :)
    !---------------------------------------------------------------------------
    DO i = 1, Glob_n
      tvk(i) = P(m_k, i)
    END DO

    tau3 = ZERO
    DO i = 1, Glob_n
      tau3 = tau3 + inv_tAkk(m_k, i) * tvk(i)
    END DO

    !---------------------------------------------------------------------------
    ! 9. Final overlap calculation
    !    Skk = 2^{3n/2} * ||L_kk||^3 / |tAkk|^{3/2} * tau3 / inv_Akk(m_k, m_k)
    !---------------------------------------------------------------------------
    temp1 = (det_Lk * det_Lk) / det_tAkk
    temp1 = ABS(temp1) * SQRT(ABS(temp1))

    Skk = Glob_2Raised3n2 * temp1 * tau3 / inv_Akk(m_k, m_k)

  END SUBROUTINE OverLapElement_S_Po

  SUBROUTINE OverLapElement_S_Pe(ml_1, ml_2, vechLl, P, Sll)
!
! The subroutine was copied from "RG_2D" code.
! This subroutine computes normalized overlap matrix element for
! a real P^e (L=1) correlated Gaussian, where the ket is symmetry
! transformed by permutation matrix P.
!
!       —                                  —
!      |                                    |
! fl = | x(ml_1)*y(ml_2) - y(ml_1)*x(ml_2)  | exp[-r'(Ll*Ll')r]
!      |                                    |
!       —                                  —
!
!        ml_1 and ml_2 are integers between 1 and Glob_n
!         (Glob_n is the number of pseudoparticles).
!        Symmetry adaption is applied to the ket using permutation matrix P
!
!
!
!                  <fl | P fl>                   || L_l ||^3
! <psi_l|psi_l> = ------------- =  2^(3*n/2) * ----------------
!                    <fl|fl>                     det(tAll)^(3/2)
!
!                                      (tgamma1 * tgamma2 - tgamma5 * tgamma6)
!                                  * --------------------------------------------
!                                        (gamma1 * gamma2 - gamma5 * gamma6)
!
! with
!
!   Al      = Ll * Ll'
!   tAll    = Al + P' * Al * P
!
!   gamma1  = v_l' * inv_All  * v_l
!   gamma2  = w_l' * inv_All  * w_l
!   gamma5  = v_l' * inv_All  * w_l
!   gamma6  = w_l' * inv_All  * v_l
!
!   tgamma1 = v_l' * inv_tAll * tv_l
!   tgamma2 = w_l' * inv_tAll * tw_l
!   tgamma5 = v_l' * inv_tAll * tw_l
!   tgamma6 = w_l' * inv_tAll * tv_l
!
!   tv_l    = P' * v_l
!   tw_l    = P' * w_l
!
!----------------------------------------------------------------------------------------------
! Input:
!
!  ml_1, ml_2     :: integers that determine which pseudoparticles carry l=1 momentum
!  vechLl         :: Array of length (n(n+1)/2) of exponential parameters.
!  P              :: The symmetry permutation matrix of size n x n
!
! Output:
!  Sll            :: Overlap matrix element (normalized)
!
!===============================================================================================

    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! ARGUMENTS
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN)      :: ml_1, ml_2
    REAL(wp), INTENT(IN)  :: vechLl(Glob_np)
    REAL(wp), INTENT(IN)  :: P(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: Sll

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !---------------------------------------------------------------------------
    INTEGER                  :: i, j, k, indx
    INTEGER                  :: n
    ! Vectors for Bra/Ket components
    REAL(wp)              :: tvl(Glob_n), twl(Glob_n)
    REAL(wp)              :: vl_inv_tAll(Glob_n), wl_inv_tAll(Glob_n)
    ! Matrices
    REAL(wp)              :: Ll(Glob_n, Glob_n)
    REAL(wp)              :: Al(Glob_n, Glob_n), tAll(Glob_n, Glob_n)
    REAL(wp)              :: inv_tAll(Glob_n, Glob_n), inv_All(Glob_n, Glob_n)
    REAL(wp)              :: W1(Glob_n, Glob_n)
    ! Scalars
    REAL(wp)              :: temp1, temp2
    REAL(wp)              :: det_tAll, det_Ll
    REAL(wp)              :: tgamma1, tgamma2, tgamma5, tgamma6
    REAL(wp)              :: gamma1, gamma2, gamma5, gamma6
    REAL(wp)              :: tm, m

    !---------------------------------------------------------------------------
    ! 1. Construct Lower-Triangular Matrix Ll from Packed Vector
    !---------------------------------------------------------------------------
    indx = 0
    DO i = 1, Glob_n
      DO j = i, Glob_n
        indx = indx + 1
        Ll(i,j) = ZERO
        Ll(j,i) = vechLl(indx)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 2. Compute Inverse of All (Unpermuted)
    !    We first compute Ll^-1, then construct inv_All = 0.5 * (L L')^-1
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    ! Invert Lower Triangular Matrix Ll -> Store in W1
    DO i = 1, Glob_n
      W1(i,i) = ONE / Ll(i,i)

      DO j = i + 1, Glob_n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - Ll(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / Ll(j,j)
      END DO

    END DO

    ! Construct inv_All = 0.5 * W1' * W1
    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_All(i,j) = ONEHALF * temp1
        inv_All(j,i) = ONEHALF * temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 3. Compute Matrix Al = Ll * Ll'
    !---------------------------------------------------------------------------
    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = ZERO
        DO k = 1, i
          temp1 = temp1 + Ll(i, k) * Ll(j,k)
        END DO
        Al(i,j) = temp1
        Al(j,i) = temp1

      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 4. Apply symmetry permutation and form tAll = Al + P' * Al * P
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    ! W1 = P' * Al
    DO i = 1, Glob_n
      DO j = 1, Glob_n
        temp1 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + P(k,j) * Al(k,i)
        END DO
        W1(j,i) = temp1
      END DO
    END DO

    ! tAll = Al + P' * Al * P
    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + W1(i, k) * P(k,j)
        END DO
        tAll(i,j) = Al(i,j) + temp1
        tAll(j,i) = tAll(i,j)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 5. Determinants and Inversion of tAll
    !---------------------------------------------------------------------------

    ! Determinant of Ll is the product of diagonal elements
    det_Ll = ONE
    DO i = 1, Glob_n
      det_Ll = det_Ll * Ll(i,i)
    END DO

    ! Cholesky Decomposition of tAll -> Store in W1
    ! Also accumulate the determinant of tAll
    det_tAll = ONE

    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = tAll(i,j)
        DO k = i - 1, 1, -1
          temp1 = temp1 - W1(i, k) * W1(j,k)
        END DO

        IF (i == j) THEN
          W1(i,i) = SQRT(temp1)

          IF (temp1 <= 0.0_wp) THEN
            WRITE(*,*) 'ERROR: tAll not SPD at i=',i
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
          END IF

          det_tAll = det_tAll * temp1
        ELSE
          W1(j,i) = temp1 / W1(i,i)
          W1(i,j) = ZERO
        END IF

      END DO
    END DO

    ! Invert tAll using its Cholesky factors (W1)
    DO i = 1, Glob_n
      W1(i,i) = ONE / W1(i,i)
      DO j = i + 1, Glob_n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - W1(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / W1(j,j)
      END DO
    END DO

    ! Construct the full inverse matrix inv_tAll
    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_tAll(i,j) = temp1
        inv_tAll(j,i) = temp1

      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 6. Compute Gamma Terms (Vector-Matrix-Vector products)
    !---------------------------------------------------------------------------

    ! Extract transformed vectors from Permutation Matrix
    DO i = 1, Glob_n
      tvl(i) = P(ml_1, i)
      twl(i) = P(ml_2, i)
    END DO

    ! Calculate Matrix-Vector Products: vector' * inv_tAll
    DO i=1,Glob_n
      vl_inv_tAll(i) = inv_tAll(ml_1,i)
      wl_inv_tAll(i) = inv_tAll(ml_2,i)
    END DO

    ! Calculate Numerator Gammas (tgamma)
    tgamma1 = ZERO
    tgamma2 = ZERO
    tgamma5 = ZERO
    tgamma6 = ZERO

    DO i = 1, Glob_n
      ! tgamma1 = vl' * inv_tAll * tvl
      tgamma1 = tgamma1 + vl_inv_tAll(i) * tvl(i)

      ! tgamma2 = wl' * inv_tAll * twl
      tgamma2 = tgamma2 + wl_inv_tAll(i) * twl(i)

      ! tgamma5 = vl' * inv_tAll * twl
      tgamma5 = tgamma5 + vl_inv_tAll(i) * twl(i)

      ! tgamma6 = wl' * inv_tAll * tvl
      tgamma6 = tgamma6 + wl_inv_tAll(i) * tvl(i)
    END DO

    ! Calculate Denominator Gammas (from inv_All direct access)
    gamma1 = inv_All(ml_1, ml_1)
    gamma2 = inv_All(ml_2, ml_2)
    gamma5 = inv_All(ml_1, ml_2)
    gamma6 = inv_All(ml_2, ml_1)

    ! Combine Gammas
    tm = tgamma1 * tgamma2 - tgamma5 * tgamma6
    m =  gamma1 *  gamma2 -  gamma5 *  gamma6

    !---------------------------------------------------------------------------
    ! 7. Final Overlap Calculation
    !---------------------------------------------------------------------------

    ! Determinant Ratio: (|Ll| * |Ll|) / |tAll|
    temp1 = (det_Ll * det_Ll) / det_tAll

    ! Apply power 3/2 (mathematically equivalent to |L|^3 / |A|^1.5)
    temp1 = ABS(temp1) * SQRT(ABS(temp1))

    ! Final Result
    Sll = Glob_2Raised3n2 * temp1 * tm / m

  END SUBROUTINE OverLapElement_S_Pe

  SUBROUTINE MatElemTranDipoleMoment_RG_1P_2P(mk, vechLk, Pk, ml_1, ml_2, vechLl, Pl, &
                                          TranDipolLength_kl, TranDipolVelocity_kl)
!
! This subroutine computes the symmetry-adapted, normalized transition
! dipole matrix element between a real P^o (L=1) bra and a real P^e (L=1)
! ket, using explicitly correlated Gaussians.
!
! Bra  (P^o):
!
!   fk = z_{mk} exp[-r'(Lk*Lk')r]
!
!        mk is an integer between 1 and n (n = number of pseudoparticles).
!        Symmetry adaption is applied to the bra using permutation matrix Pk.
!
!   Pk fk = (Pk z_{mk}) exp[-r' (Pk'*Lk*Lk'*Pk) r]
!
!
! Ket  (P^e):
!
!        —                                  —
!       |                                    |
!  fl = | x(ml_1)*y(ml_2) - y(ml_1)*x(ml_2)  | exp[-r'( Ll' * Ll )r]
!       |                                    |
!        —                                  —
!
!        ml_1 and ml_2 are integers between 1 and n.
!        Symmetry adaption is applied to the ket using permutation matrix Pl.
!
!
!            —                                                        —
!           |                                                          |
!   Pl fl = | Pl x_{ml_1})*(Pl y_{ml_2}) - (Pl y_{ml_1})*(Pl x_{ml_2}) | * exp[-r' (Pl'*Ll*Ll'*Pl) r]
!           |                                                          |
!            —                                                        —
!
!
!
! This subroutine computes the following expressions, which are parts of
! the transition dipole moment calculation.
!
!
!====================================================================
! Transition dipole integral in Length form  [Eq.(99)]
!====================================================================
!
!                      ---
!                      \                  m_i         < Pk*fk | x_i | Pl*fl >
! TranDipolLength_kl = /  (q_i - Q_tot * -----) * --------------------------------------
!                      ---                m_0       Sqrt( <fk|fk> ) * Sqrt( <fl|fl> )
!                       i
!
!
!
!             2^(3*n/2)      |det_Lk|^1.5 * |det_Ll|^1.5
!          = ----------- * --------------------------------
!                2                 |det_tAkl|^1.5
!
!                                        1
!          * -------------------------------------------------------------
!              Sqrt(vk'*inv_Akk*vk) * Sqrt(gamma1*gamma2 - gamma5*gamma6)
!
!
!             ---
!             \                  m_i
!          *  /  (q_i - Q_tot * -----) *( tgamma5 * tgamma6 - tgamma1 * tgamma2 )
!             ---                m_0
!              i
!
!
!====================================================================
! Transition dipole integral in Velocity form  [Eq.(134)]
!====================================================================
!
!
! TranDipolVelocity_kl =
!
!               ---
!               \     q_0   q_i          < Pk*fk | px_i | Pl*fl >
!               /   (----- - ---) * --------------------------------------
!               ---   m_0   m_i       Sqrt( <fk|fk> ) * Sqrt( <fl|fl> )
!                i
!
!                                 |det_Lk|^1.5 * |det_Ll|^1.5
!           =   i * 2^(3*n/2) * ---------------------------------
!                                       |det_tAkl|^1.5
!
!                                          1
!           * --------------------------------------------------------------
!              Sqrt(vk'*inv_Akk*vk) * Sqrt(gamma1*gamma2 - gamma5*gamma6)
!
!
!              ---
!              \     q_0   q_i
!           *  /   (----- - ---) (tgamma1 * tgamma2_v - tgamma5 * tgamma6_v)
!              ---   m_0   m_i
!               i
!
! where the imaginary factor i is dropped (it cancels in |<...>|^2).
!
!
!====================================================================
! Gamma definitions
!====================================================================
!
! Normalization gammas (from inv_All, the P^e norm):
!   gamma1 = v_l' * inv_All * v_l = inv_All(ml_1, ml_1)
!   gamma2 = w_l' * inv_All * w_l = inv_All(ml_2, ml_2)
!   gamma5 = v_l' * inv_All * w_l = inv_All(ml_1, ml_2)
!   gamma6 = w_l' * inv_All * v_l = inv_All(ml_2, ml_1)
!
! Length-gauge tilde gammas (from inv_tAkl):
!   tgamma1 = tvk' * inv_tAkl * tvl       (i-independent)
!   tgamma5 = tvk' * inv_tAkl * twl       (i-independent)
!   tgamma2 = v_i' * inv_tAkl * twl       (i-dependent, = inv_tAkl_twl(i))
!   tgamma6 = v_i' * inv_tAkl * tvl       (i-dependent, = inv_tAkl_tvl(i))
!
! Velocity-gauge tilde gammas:
!   tgamma1   = tvk' * inv_tAkl * tvl     (same as length)
!   tgamma5   = tvk' * inv_tAkl * twl     (same as length)
!   tgamma2_v = u_i' * inv_tAkl * twl     (i-dependent, u_i = tAk * v_i)
!   tgamma6_v = u_i' * inv_tAkl * tvl     (i-dependent, u_i = tAk * v_i)
!
! Transformed vectors:
!   tvk = Pk' * v_k   ->  tvk(j) = Pk(mk, j)
!   tvl = Pl' * v_l   ->  tvl(j) = Pl(ml_1, j)
!   twl = Pl' * w_l   ->  twl(j) = Pl(ml_2, j)
!
!----------------------------------------------------------------------------------------------
! Input:
!   mk              :: integer that determines the P^o bra pseudoparticle (z-component)
!   ml_1, ml_2      :: integers that determine the P^e ket pseudoparticles (xy - yx angular part)
!   vechLk, vechLl  :: Arrays of length n(n+1)/2 of exponential parameters for bra and ket
!   Pk, Pl          :: Symmetry permutation matrices of size n x n for bra and ket
!
! Output:
!   TranDipolLength_kl   :: Normalized length-gauge matrix element
!   TranDipolVelocity_kl :: Normalized velocity-gauge matrix element
!
!==============================================================================================

    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! ARGUMENTS
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN)      :: mk, ml_1, ml_2
    REAL(wp), INTENT(IN)  :: vechLk(Glob_np), vechLl(Glob_np)
    REAL(wp), INTENT(IN)  :: Pk(Glob_n, Glob_n), Pl(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: TranDipolLength_kl, TranDipolVelocity_kl

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !---------------------------------------------------------------------------
    INTEGER                  :: Qtotal
    INTEGER                  :: i, j, k, indx
    REAL(wp)              :: Lk(Glob_n, Glob_n), Ll(Glob_n, Glob_n)
    REAL(wp)              :: tAk(Glob_n, Glob_n), tAl(Glob_n, Glob_n), tAkl(Glob_n, Glob_n)
    REAL(wp)              :: W1(Glob_n, Glob_n), W2(Glob_n, Glob_n)

    REAL(wp)              :: inv_Akk(Glob_n, Glob_n), inv_All(Glob_n, Glob_n), inv_tAkl(Glob_n, Glob_n)
    REAL(wp)              :: tvk(Glob_np)                ! Transformed bra vector: Pk' * e_{mk}
    REAL(wp)              :: tvl(Glob_np), twl(Glob_np)  ! Transformed ket vectors: Pl' * e_{ml_1}, Pl' * e_{ml_2}

    REAL(wp)              :: tvk_inv_tAkl(Glob_np)       ! tvk' * inv_tAkl (row vector)
    REAL(wp)              :: inv_tAkl_twl(Glob_np)       ! inv_tAkl * twl   (column vector)
    REAL(wp)              :: inv_tAkl_tvl(Glob_np)       ! inv_tAkl * tvl   (column vector)
    REAL(wp)              :: det_Lk, det_Ll, det_tAkl
    REAL(wp)              :: temp01, temp02, temp03, temp1, temp2
    REAL(wp)              :: charge_mass0

    REAL(wp)              :: tgamma1, tgamma5, tgamma2, tgamma6     ! Length-gauge tilde gammas
    REAL(wp)              :: tgamma2_v, tgamma6_v                   ! Velocity-gauge tilde gammas (u_i-based)
    REAL(wp)              :: gamma1, gamma2, gamma5, gamma6         ! Normalization gammas (P^e norm)
    REAL(wp)              :: tm1, tm3, tm_num, m_den, m1, m3
    REAL(wp)              :: tm1_v, tm3_v, tm_num_v

!---------------------------------------------------------------------------
!  Build lower-triangular matrices Lk and Ll from packed vector storage
!---------------------------------------------------------------------------
    indx = 0

    DO i = 1, Glob_n
      DO j = i, Glob_n
        indx = indx + 1
        Lk(i,j) = ZERO
        Lk(j,i) = vechLk(indx)
        Ll(i,j) = ZERO
        Ll(j,i) = vechLl(indx)
      END DO
    END DO

!---------------------------------------------------------------------------
!  Compute A matrices: Ak = Lk * Lk',  Al = Ll * Ll'
!---------------------------------------------------------------------------
    DO i = 1, Glob_n
      DO j = i, Glob_n
        ! Compute Ak
        temp1 = ZERO
        DO k = 1, i
          temp1 = temp1 + Lk(i, k) * Lk(j,k)
        END DO
        tAk(i,j) = temp1
        tAk(j,i) = temp1

        ! Compute Al
        temp1 = ZERO
        DO k = 1, i
          temp1 = temp1 + Ll(i, k) * Ll(j,k)
        END DO
        tAl(i,j) = temp1
        tAl(j,i) = temp1
      END DO
    END DO

!---------------------------------------------------------------------------
! Apply symmetry permutations to form
!      tAk <- Pk' * Ak * Pk    (permuted bra exponent)
!      tAl <- Pl' * Al * Pl    (permuted ket exponent)
!      tAkl = tAk + tAl        (combined exponent matrix)
!---------------------------------------------------------------------------
    W1(:,:) = ZERO
    W2(:,:) = ZERO

    ! W1 = Pl' * Al  (intermediate for ket)
    ! W2 = Ak * Pk   (intermediate for bra)
    DO i = 1, Glob_n
      DO j = 1, Glob_n

        temp1 = ZERO
        temp2 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + Pl(k,j) * tAl(k,i)
          temp2 = temp2 + tAk(j,k) * Pk(k,i)
        END DO
        W1(j,i) = temp1
        W2(j,i) = temp2

      END DO
    END DO

    ! Complete the permutation and sum:
    !   tAl <- Pl' * Al * Pl = W1 * Pl
    !   tAk <- Pk' * Ak * Pk = Pk' * W2
    !   tAkl = tAk + tAl
    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = ZERO
        temp2 = ZERO
        DO k = 1, Glob_n
          temp1 = temp1 + W1(j,k) * Pl(k,i)
          temp2 = temp2 + Pk(k,j) * W2(k,i)
        END DO
        tAl(j,i)  = temp1
        tAl(i,j)  = temp1
        tAk(j,i)  = temp2
        tAk(i,j)  = temp2
        tAkl(j,i) = temp1 + temp2
        tAkl(i,j) = temp1 + temp2

      END DO
    END DO

!---------------------------------------------------------------------------
!  Determinants of Lk, Ll and Cholesky factorization / inversion of tAkl
!---------------------------------------------------------------------------

    ! Determinants of L matrices (product of diagonal elements)
    det_Lk = ONE
    det_Ll = ONE

    DO i = 1, Glob_n
      det_Lk = det_Lk * Lk(i,i)
      det_Ll = det_Ll * Ll(i,i)
    END DO

    ! Cholesky decomposition of tAkl -> store factor in W1
    ! Also accumulate det(tAkl) as the product of diagonal^2 terms
    det_tAkl = ONE
    W1(:,:) = ZERO

    DO i = 1, Glob_n
      DO j = i, Glob_n
        temp1 = tAkl(i,j)

        DO k = i - 1, 1, -1
          temp1 = temp1 - W1(i, k) * W1(j,k)
        END DO

        IF (i == j) THEN
          IF (temp1 <= 0.0_wp) THEN
            WRITE(*,*) 'ERROR: tAkl not SPD at i=',i,' temp1=',temp1,' mk=',mk,' ml1=',ml_1,' ml2=',ml_2
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
          END IF
          W1(i,i) = SQRT(temp1)
          det_tAkl = det_tAkl * temp1
        ELSE
          W1(j,i) = temp1 / W1(i,i)
          W1(i,j) = ZERO
        END IF

      END DO
    END DO

    ! Invert the Cholesky factor (lower triangular) in-place
    DO i = 1, Glob_n

      W1(i,i) = ONE / W1(i,i)
      DO j = i + 1, Glob_n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - W1(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / W1(j,j)
      END DO

    END DO

    ! Construct full inverse: inv_tAkl = W1' * W1
    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_tAkl(i,j) = temp1
        inv_tAkl(j,i) = temp1

      END DO
    END DO

!---------------------------------------------------------------------------
! Compute transformed vectors (bra and ket)
!    tvk(j) = Pk(mk, j)    i.e.   mk-th row of Pk = Pk' * e_{mk}
!    tvl(j) = Pl(ml_1, j)  i.e. ml_1-th row of Pl = Pl' * e_{ml_1}
!    twl(j) = Pl(ml_2, j)  i.e. ml_2-th row of Pl = Pl' * e_{ml_2}
!---------------------------------------------------------------------------
    DO i = 1, Glob_n
      tvk(i) = Pk(mk, i)     ! Bra: permuted P^o vector
      tvl(i) = Pl(ml_1, i)   ! Ket: permuted P^e vector v_l
      twl(i) = Pl(ml_2, i)   ! Ket: permuted P^e vector w_l
    END DO

!---------------------------------------------------------------------------
! Compute inverse of unpermuted A matrices:
!      inv_Akk = 0.5 * (Lk Lk')^{-1}   (for P^o normalization)
!      inv_All = 0.5 * (Ll Ll')^{-1}   (for P^e normalization)
!---------------------------------------------------------------------------
    W1(:,:) = ZERO
    W2(:,:) = ZERO

    ! Invert lower-triangular Lk -> W1,  Ll -> W2
    DO i = 1, Glob_n

      W1(i,i) = ONE / Lk(i,i)
      W2(i,i) = ONE / Ll(i,i)

      DO j = i + 1, Glob_n
        temp1 = ZERO
        temp2 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - Lk(j,k) * W1(k,i)
          temp2 = temp2 - Ll(j,k) * W2(k,i)
        END DO
        W1(j,i) = temp1 / Lk(j,j)
        W2(j,i) = temp2 / Ll(j,j)
      END DO

    END DO

    ! inv_Akk = 0.5 * W1' * W1,  inv_All = 0.5 * W2' * W2
    DO i = 1, Glob_n
      DO j = i, Glob_n

        temp1 = ZERO
        temp2 = ZERO
        DO k = j, Glob_n
          temp1 = temp1 + W1(k,i) * W1(k,j)
          temp2 = temp2 + W2(k,i) * W2(k,j)
        END DO
        inv_Akk(i,j) = ONEHALF * temp1
        inv_Akk(j,i) = ONEHALF * temp1
        inv_All(i,j) = ONEHALF * temp2
        inv_All(j,i) = ONEHALF * temp2

      END DO
    END DO

!---------------------------------------------------------------------------
! Compute P^e normalization denominator
!    From Eq.(32):  norm ~ (gamma1*gamma2 - gamma5*gamma6)
!    gamma1 = v_l' * inv_All * v_l = inv_All(ml_1, ml_1)
!    gamma2 = w_l' * inv_All * w_l = inv_All(ml_2, ml_2)
!    gamma5 = v_l' * inv_All * w_l = inv_All(ml_1, ml_2)
!    gamma6 = w_l' * inv_All * v_l = inv_All(ml_2, ml_1)
!---------------------------------------------------------------------------
    gamma1 = inv_All(ml_1, ml_1)
    gamma2 = inv_All(ml_2, ml_2)
    gamma5 = inv_All(ml_1, ml_2)
    gamma6 = inv_All(ml_2, ml_1)

    m1 = gamma1 * gamma2
    m3 = gamma5 * gamma6

    ! P^e norm denominator: sqrt(gamma1*gamma2 - gamma5*gamma6)
    m_den = SQRT(m1 - m3)

    IF(Verbose==2 .OR. Verbose==3) THEN

      IF (m_den /= m_den) THEN
        WRITE(*,*)
        WRITE(*,*) 'WARNING: m_den is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      END IF

      IF(Verbose==3) THEN
        WRITE(*,*)
        WRITE(*,'(12X,A31,ES12.4)') 'Denominator (m_den): ', m_den
      END IF

    END IF

!---------------------------------------------------------------------------
! Compute constant prefactor
!    Length:   2^(3n/2) / 2  *  |Lk|^{3/2} |Ll|^{3/2} / |tAkl|^{3/2}
!                            *  1 / sqrt(vk' inv_Akk vk)
!    Velocity: prefactor is 2x the length prefactor (factor i dropped)
!---------------------------------------------------------------------------

    ! Determinant ratio: |Lk|*|Ll| / |tAkl|
    temp1 = ABS(det_Ll * det_Lk) / det_tAkl

    IF(Verbose==2 .OR. Verbose==3) THEN
      IF(Verbose==3) &
        WRITE(*,'(12X,A32,ES12.4)') '|Lk|^1.5*|Ll|^1.5 / |tAkl|^1.5: ', temp1
      IF (temp1 /= temp1) THEN
        WRITE(*,*) 'WARNING: |Lk|^1.5*|Ll|^1.5 / |tAkl|^1.5 is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      ENDIF
    END IF

    ! Raise to power 3/2: (|Lk|*|Ll|)^{3/2} / |tAkl|^{3/2}
    temp1 = temp1 * SQRT(temp1)

    ! Multiply by 2^(3n/2) / 2
    temp1 = temp1 * Glob_2Raised3n2
    temp1 = temp1 / TWO

    IF(Verbose==2 .OR. Verbose==3) THEN
      IF(Verbose==3) &
        WRITE(*,'(12X,A31,ES12.4)') 'SQRT(inv_Akk(mk,mk)): ', SQRT(inv_Akk(mk, mk))

      IF (SQRT(inv_Akk(mk, mk)) /= SQRT(inv_Akk(mk, mk))) THEN
        WRITE(*,*) 'WARNING: SQRT(inv_Akk(mk, mk)) is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      ENDIF
    END IF

    ! Divide by sqrt(vk' * inv_Akk * vk) = sqrt(inv_Akk(mk, mk))
    ! This is the P^o normalization factor
    temp1 = temp1 / SQRT(inv_Akk(mk, mk))

    ! Store the common prefactor divided by the P^e norm denominator
    ! Length gauge uses this directly; velocity gauge is multiplied by 2
    TranDipolLength_kl = temp1 / m_den
    TranDipolVelocity_kl = TranDipolLength_kl * TWO

!---------------------------------------------------------------------------
!  Pre-compute i-independent quantities for the dipole summation
!---------------------------------------------------------------------------

    ! Compute tvk' * inv_tAkl (row vector, used for tgamma1 and tgamma5)
    DO i = 1, Glob_n
      temp1 = ZERO

      DO j = 1, Glob_n
        temp1 = temp1 + tvk(j) * inv_tAkl(j,i)
      END DO

      tvk_inv_tAkl(i) = temp1
    END DO

    ! tgamma1 = tvk' * inv_tAkl * tvl   (i-independent)
    ! tgamma5 = tvk' * inv_tAkl * twl   (i-independent)
    tgamma1 = ZERO
    tgamma5 = ZERO

    DO j = 1, Glob_n
      tgamma1 = tgamma1 + tvk_inv_tAkl(j) * tvl(j)
      tgamma5 = tgamma5 + tvk_inv_tAkl(j) * twl(j)
    END DO

    ! Pre-compute inv_tAkl * twl and inv_tAkl * tvl (column vectors)
    ! These provide tgamma2(i) and tgamma6(i) for the length gauge:
    !   tgamma2 = v_i' * inv_tAkl * twl = inv_tAkl_twl(i)
    !   tgamma6 = v_i' * inv_tAkl * tvl = inv_tAkl_tvl(i)

    DO i = 1, Glob_n

      temp1 = ZERO
      temp2 = ZERO

      DO j = 1, Glob_n
        temp1 = temp1 + inv_tAkl(i,j) * twl(j)
        temp2 = temp2 + inv_tAkl(i,j) * tvl(j)
      END DO

      inv_tAkl_twl(i) = temp1
      inv_tAkl_tvl(i) = temp2

    END DO

!---------------------------------------------------------------------------
! Dipole moment summation loop over pseudoparticles i = 1..Glob_n
!     Length gauge:   sum_i (q_i - Q_tot*m_i/m_0) *
!                           (tgamma5*tgamma6 - tgamma1*tgamma2)
!
!     Velocity gauge: sum_i (q_0/m_0 - q_i/m_i) *
!                           (tgamma1*tgamma2_v - tgamma5*tgamma6_v)
!---------------------------------------------------------------------------
    Qtotal = Glob_PseudoCharge0

    DO i = 1, Glob_n
      Qtotal = Qtotal + Glob_PseudoCharge(i)
    END DO

    charge_mass0 = Glob_PseudoCharge0 / Glob_Mass(1)

    temp01 = ZERO
    temp03 = ZERO

    DO i = 1, Glob_n

      !------------------------------------------
      ! Length-gauge gammas (i-dependent part):
      !   tgamma2 = v_i' * inv_tAkl * twl = inv_tAkl_twl(i)
      !   tgamma6 = v_i' * inv_tAkl * tvl = inv_tAkl_tvl(i)
      !------------------------------------------
      tgamma2 = inv_tAkl_twl(i)
      tgamma6 = inv_tAkl_tvl(i)

      ! Length-gauge numerator: tgamma5*tgamma6 - tgamma1*tgamma2
      tm1 = tgamma1 * tgamma2
      tm3 = tgamma5 * tgamma6
      tm_num = tm3 - tm1

      IF (tm_num /= tm_num) THEN
        WRITE(*,*) 'WARNING: tm_num is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      END IF

      ! Length-gauge charge-mass factor: (q_i - Q_tot * m_i / m_total)
      temp02 = Glob_PseudoCharge(i) - Qtotal * Glob_Mass(i + 1) / Glob_MassTotal

      ! Accumulate length-gauge sum
      temp01 = temp01 + temp02 * tm_num

      !------------------------------------------
      ! Velocity-gauge gammas (i-dependent part):
      !   u_i = tAk * v_i  ->  u_i(j) = tAk(j,i)
      !   tgamma2_v = u_i' * inv_tAkl * twl
      !   tgamma6_v = u_i' * inv_tAkl * tvl
      !------------------------------------------
      tgamma2_v = ZERO
      tgamma6_v = ZERO

      DO j = 1, Glob_n
        tgamma2_v = tgamma2_v + tAk(j,i) * inv_tAkl_twl(j)
        tgamma6_v = tgamma6_v + tAk(j,i) * inv_tAkl_tvl(j)
      END DO

      ! Velocity-gauge numerator: tgamma1*tgamma2_v - tgamma5*tgamma6_v
      tm1_v = tgamma1 * tgamma2_v
      tm3_v = tgamma5 * tgamma6_v
      tm_num_v = tm1_v - tm3_v

      IF (tm_num_v /= tm_num_v) THEN
        WRITE(*,*) 'WARNING: tm_num_v is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      END IF

      ! Velocity-gauge charge-mass factor: (q_0/m_0 - q_i/m_i)
      temp02 = charge_mass0 - Glob_PseudoCharge(i)/Glob_Mass(i+1)
      temp03 = temp03 + temp02 * tm_num_v

      IF(Verbose==2 .OR. Verbose==3) THEN

        IF(Verbose==3) THEN
          WRITE(*,'(29X,A3,I1,A10,ES12.4)') 'i= ',i,', tm_num: ', tm_num
          WRITE(*,'(35X,A10,ES12.4)') 'tm_num_v: ',tm_num_v
        ENDIF

      END IF

    END DO

    !---------------------------------------------------------------------
    ! Final results:
    !     Length:   prefactor * sum
    !     Velocity: prefactor * sum    (imaginary unit i is dropped)
    !---------------------------------------------------------------------

    TranDipolLength_kl   = TranDipolLength_kl * temp01

    IF(Verbose==2 .OR. Verbose==3) THEN

      IF (TranDipolLength_kl /= TranDipolLength_kl) THEN
        WRITE(*,*) 'WARNING: TranDipolLength_kl is NaN'
        call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
      END IF

      IF(Verbose==3) THEN
        WRITE(*,'(12X,A31,ES12.4)') 'TranDipolLength_kl: ', TranDipolLength_kl
      END IF

    END IF

    TranDipolVelocity_kl = TranDipolVelocity_kl * temp03

  END SUBROUTINE MatElemTranDipoleMoment_RG_1P_2P

end module matelem

