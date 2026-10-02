!DEC$ DECLARE

MODULE matelem
  ! Module matelem contains procedures of computing
  ! matrix elements.

  USE globvars
  USE wp_def
  USE mpi


  IMPLICIT NONE


  PRIVATE :: compute_g3, compute_products_1, compute_products_2, &
             compute_products_3, compute_products_4, compute_products_5, &
             compute_dVk, &
             compute_series

CONTAINS


  SUBROUTINE MatrixElementsOpt(mk, PreLk, PreAk, PreMAk, ml, vechLl, SymMatrixBuf, &
      SklBuf, HklBuf, dSkBuf, dSlBuf, dHkBuf, dHlBuf, &
      gradflag)
    ! Symmetry adapted matrix elements in a basis of correlated Gaussians
    ! with r1 premultipliers,  fk = r1^mk * exp[-r'(Lk*Lk' (x) I3) r].
    ! The symmetry adaption is applied to the ket:
    ! Al --> SymMatrix'*Ll*Ll'*SymMatrix.
    ! Input:
    !   PreLk, PreAk, PreMAk :: bra side, precomputed by PrecomputeMatrices_MA_L_A:
    !                           Lk, Ak = Lk*Lk' and M*Ak (no permutation needed)
    !   vechLl         :: n(n+1)/2 exponent parameters of the ket (Ll lower triangular)
    !   mk, ml         :: premultiplier powers
    !   SymMatrixBuf   :: the symmetry permutation matrices
    !   gradflag       :: 0 no gradient; 1 dSk, dHk only; 2 also dSl, dHl
    ! Output (normalized):
    !   Skl, Hkl       :: overlap and Hamiltonian matrix elements
    !   dSk, dHk       :: derivatives with respect to vech[Lk]
    !   dSl, dHl       :: derivatives with respect to vech[Ll]

    ! Parameters
    ! Use Glob_n from globvars as the number of pseudoparticles.
    !================================================================
    IMPLICIT NONE


    INTEGER, INTENT(IN) :: mk, ml, gradflag

    REAL(wp), INTENT(IN) :: PreLk(Glob_n, Glob_n), PreAk(Glob_n, Glob_n), PreMAk(Glob_n, Glob_n)
    REAL(wp), INTENT(IN) :: vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrixBuf(Glob_NumYHYTerms, Glob_n, Glob_n)

    REAL(wp), INTENT(OUT) :: SklBuf(Glob_NumYHYTerms), HklBuf(Glob_NumYHYTerms)
    REAL(wp), INTENT(OUT) :: dSkBuf(Glob_NumYHYTerms, Glob_npt), &
                                dSlBuf(Glob_NumYHYTerms, Glob_npt), &
                                dHkBuf(Glob_NumYHYTerms, Glob_npt), &
                                dHlBuf(Glob_NumYHYTerms, Glob_npt)

    ! Local variables

    REAL(wp) :: SymMatrix(Glob_n, Glob_n)
    REAL(wp) :: Skl, Hkl
    REAL(wp) :: dSk(Glob_npt), dSl(Glob_npt), dHk(Glob_npt), dHl(Glob_npt)

    INTEGER :: m                                            ! matrix index
    INTEGER :: i, j, k, indx, info, jj, ii, kk, mko2, mlo2
    INTEGER :: p, q, t

    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: detLk, detLl, detAkl
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: SUM0, SUM1, SUM2, SUM3
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: trT, tl, tk
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: Ttmp
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: a
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: Tkl, Vkl
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: g1, g2, g3

    REAL(wp) :: mkml, mkpml
    REAL(wp) :: coabij
    REAL(wp) :: Vtmp(Glob_n, Glob_n)

    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Lk, Ll, Ak, Al, invAkl, PLl, invL
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Temp
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: R, b, c, coab
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ddSk, ddSl, ddTk, ddTl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ddSk_over_Skl, ddSl_over_Skl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: dRk, dRl, ddVk, ddVl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: invLk, invLl, invAk, invAl, chAl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ALk, ALl, AAk, AAl, AlLl, MAk, MAl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Jkk, Jll, Jklk, Jkll
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: AAlMAl, AAkMAk, AAkMAl, AAlMAk
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Tllk, Tkkl, Tklk, Tlkl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: TlMkJk, TkMlJl, TlMlJk, TkMkJl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: JkMlk, JlMkl, JlMlk, JkMkl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: da, db, dc
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: RLk, RLl, JRLk, JRLl, RJLk, RJLl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_npt)       :: dTk, dTl, dVk, dVl

    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Per-pair tracing (Verbose >= 4 floods: one line per pair on every rank)
    IF (Verbose >= 4) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,i0,a,i0)') 'MatrixElementsOpt: powers ', mk, ' ', ml, &
                                            '  gradflag ', gradflag, '  rank ', Glob_ProcID
    ENDIF

    ! Bra side: Lk, Ak = Lk*Lk' and MAk = M*Ak are precomputed once per
    ! function; they are the same for every symmetry term
    DO j = 1, Glob_n
      DO i = 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          Lk(m, i, j) = PreLk(i, j)
          Ak(m, i, j) = PreAk(i, j)
          MAk(m, i, j) = PreMAk(i, j)
        ENDDO
      ENDDO
    ENDDO

    ! Build the lower triangular matrix Ll from the elements in vechLl
    Ll = ZERO
    indx = 0
    DO j = 1, Glob_n
      DO i = j, Glob_n
        indx = indx + 1
        DO m = 1, Glob_NumYHYTerms
          Ll(m, i, j) = vechLl(indx)
        ENDDO
      ENDDO
    ENDDO

    ! Apply the symmetry transformation SymMatrix to Ll  PLl=SymMatrix'*Ll
    PLl = ZERO
    DO j = 1, Glob_n
      DO k = j, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            PLl(m, i, j) = PLl(m, i, j) + SymMatrixBuf(m, k, i) * Ll(m, k, j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Build AL=P'*LL*LL'*P (Ak is precomputed)
    Al = ZERO
    DO j = 1, Glob_n
      DO k = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            Al(m, i, j) = Al(m, i, j) + PLl(m, i, k) * PLl(m, j, k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Compute Akl the upper triangle of Akl=Ak+Al is in invL
    ! (we store the upper part of Akl for convience in the
    ! Cholesky decomposition code below)
    DO j = 1, Glob_n
      DO i = 1, j - 1
        DO m = 1, Glob_NumYHYTerms
          invL(m, j, i) = ZERO
        ENDDO
      ENDDO
      DO i = j, Glob_n
        DO m = 1, Glob_NumYHYTerms
          invL(m, j, i) = Ak(m, i, j) + Al(m, i, j)
        ENDDO
      ENDDO
    ENDDO

    ! compute det(LK) and det(LL)
    detLk = ONE
    detLl = ONE
    DO i = 1, Glob_n
      DO m = 1, Glob_NumYHYTerms
        detLk(m) = detLk(m) * Lk(m, i, i)
        detLl(m) = detLl(m) * Ll(m, i, i)
      ENDDO
    ENDDO

    ! Compute the Cholesky decomposition of Akl (invL)
    ! This code is adapted from Numerical Recipes  choldc.for
    ! det(Akl) is computed during the decomposition and stored in
    ! detAkl
    !
    ! need to do the decomp on Al too to get invAl since
    ! Al is P'Ll*(P'Ll)' i.e. we lost the original chol decomp
    ! we'll put it in chAl
    chAl = Al
    detAkl = ONE
    DO i = 1, Glob_n
      DO j = i, Glob_n
        DO m = 1, Glob_NumYHYTerms
          SUM0(m) = invL(m, i, j)
          SUM1(m) = chAl(m, i, j)
        ENDDO
        DO k = i - 1, 1, -1
          DO m = 1, Glob_NumYHYTerms
            SUM0(m) = SUM0(m) - invL(m, i, k) * invL(m, j, k)
            SUM1(m) = SUM1(m) - chAl(m, i, k) * chAl(m, j, k)
          ENDDO
        ENDDO
        IF (i == j) THEN
          DO m = 1, Glob_NumYHYTerms
            invL(m, i, i) = SQRT(SUM0(m))
            chAl(m, i, i) = SQRT(SUM1(m))  ! don't need the det for Al
            detAkl(m) = detAkl(m) * SUM0(m)
          ENDDO
        ELSE
          DO m = 1, Glob_NumYHYTerms
            invL(m, j, i) = SUM0(m) / invL(m, i, i)
            chAl(m, j, i) = SUM1(m) / chAl(m, i, i)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    ! invL now contains the Cholesky decomp,L, of Akl, where Akl=L*L'
    ! detAkl now contains det(Akl)
    ! chAl now contains the Cholesky decomp of Al
    !
    ! Compute the inverse of Akl, Ak, Al, Lk, Ll
    !
    ! Start with invL, invLk, invLl, chAl
    invLk = Lk
    invLl = Ll
    ! [invert the lower triangle matrices]
    DO i = 1, Glob_n
      DO m = 1, Glob_NumYHYTerms
        invL(m, i, i) = ONE / invL(m, i, i)
        invLk(m, i, i) = ONE / invLk(m, i, i)
        invLl(m, i, i) = ONE / invLl(m, i, i)
        chAl(m, i, i) = ONE / chAl(m, i, i)
      ENDDO
      DO j = i + 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          SUM0(m) = ZERO
          SUM1(m) = ZERO
          SUM2(m) = ZERO
          SUM3(m) = ZERO
        ENDDO
        DO k = i, j - 1
          DO m = 1, Glob_NumYHYTerms
            SUM0(m) = SUM0(m) - invL(m, j, k) * invL(m, k, i)
            SUM1(m) = SUM1(m) - invLk(m, j, k) * invLk(m, k, i)
            SUM2(m) = SUM2(m) - invLl(m, j, k) * invLl(m, k, i)
            SUM3(m) = SUM3(m) - chAl(m, j, k) * chAl(m, k, i)
          ENDDO
        ENDDO
        DO m = 1, Glob_NumYHYTerms
          invL(m, j, i) = SUM0(m) / invL(m, j, j)
          invL(m, i, j) = ZERO  ! zero out the upper triangle for convenience
          invLk(m, j, i) = SUM1(m) / invLk(m, j, j)
          invLk(m, i, j) = ZERO
          invLl(m, j, i) = SUM2(m) / invLl(m, j, j)
          invLl(m, i, j) = ZERO
          chAl(m, j, i) = SUM3(m) / chAl(m, j, j)
          chAl(m, i, j) = ZERO
        ENDDO
      ENDDO
    ENDDO

    ! We now have invLk, invLl (this is not invP'Ll), inv(chAl) invL
    !
    ! Now do invAkl=invL'*invL
    !       invAk=invLk'*invLk
    !       invAl=chAl'*chAl    NOTE that this has the SymMatrix term in it
    invAkl = ZERO
    invAk = ZERO
    invAl = ZERO
    DO j = 1, Glob_n
      DO k = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            invAkl(m, i, j) = invAkl(m, i, j)+ invL(m, k, i) * invL(m, k, j)
            invAk(m, i, j) = invAk(m, i, j) + invLk(m, k, i) * invLk(m, k, j)
            invAl(m, i, j) = invAl(m, i, j) + chAl(m, k, i) * chAl(m, k, j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! OK looks like we got most of the important terms
    ! Lk, Ll, PLl, Ak, Al=P'Ll(P'Ll)', invAkl, invAk, invAl
    ! invLk, invLl [not inv(P'Ll)]


    !********************************************************
    ! Compute intermediate matrix results LOTS OF THEM!
    !********************************************************

    ! ALk=invAkl*Lk
    ! ALl=invAkl*PLl
    ! AlLl=invAl*PLl need this because of sym
    ! AAk=invAkl*Ak
    ! AAl=invAkl*Al

    CALL compute_products_1(Lk, Ak, Al, invAkl, PLl, invAl, ALk, ALl, AAk, AAl, AlLl)


    ! Now a bunch of rank one mats and do it a column at a time!
    ! Jkk = invAk*J11*invAk*Lk
    ! Jll = invAl*J11*invAl*Ll
    ! Jklk = invAkl*J11*invAkl*Lk
    ! Jkll = invAkl*J11*invAkl*Ll
    DO j = 1, Glob_n
      DO i = 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          Jkk(m, i, j) = invAk(m, i, 1) * invLk(m, j, 1)
          Jll(m, i, j) = invAl(m, i, 1) * AlLl(m, 1, j)  ! using AlLl to get sym right
          Jklk(m, i, j) = invAkl(m, i, 1) * ALk(m, 1, j)
          Jkll(m, i, j) = invAkl(m, i, 1) * ALl(m, 1, j)
        ENDDO
      ENDDO
    ENDDO

    ! MAl=M*Al (MAk=M*Ak is precomputed)

    CALL compute_products_2(Al, MAl)

    ! AAlMAl=invAkl*Al*M*Al
    ! AAkMAk=invAkl*Ak*M*Ak
    ! AAkMAl=invAkl*Ak*M*Al
    ! AAlMAk=invAkl*Al*M*Ak

    CALL compute_products_3(AAk, AAl, MAk, MAl, AAlMAl, AAkMAk, AAkMAl, AAlMAk)


    ! A few kinetic energy terms (scalers)
    ! trT=tr[invAkl*Ak*M*Al]
    ! tl=(invAkl*Al*M*Al*invAkl)(1,1)
    ! tk=(invAkl*Ak*M*Ak*invAkl)(1,1)
    trT = ZERO
    tl = ZERO
    tk = ZERO
    DO i = 1, Glob_n
      DO m = 1, Glob_NumYHYTerms
        trT(m) = trT(m) + AAkMAl(m, i, i)
        tl(m) = tl(m) + AAlMAl(m, 1, i) * invAkl(m, i, 1)
        tk(m) = tk(m) + AAkMAk(m, 1, i) * invAkl(m, i, 1)
      ENDDO
    ENDDO

    ! More kinetic energy (gradient) matrix terms
    IF (gradflag > 0) THEN
      ! Tllk=invAkl*Al*M*Al*invAkl*Lk
      ! Tkkl=invAkl*Ak*M*Ak*invAkl*Ll
      ! Tklk=invAkl*Ak*M*Al*invAkl*Lk
      ! Tlkl=invAkl*Al*M*Ak*invAkl*Ll

      CALL compute_products_4(ALk, ALl, AAlMAl, AAkMAk, AAkMAl, AAlMAk, &
          Tllk, Tkkl, Tklk, Tlkl)

      ! Still more kinetic energy (gradient) matrix terms
      ! TlMkJk = invAkl*Al*M*Al*invAkl*Lk
      ! TkMlJl = invAkl*Ak*M*Ak*invAkl*Ll
      ! TlMlJk = invAkl*Ak*M*Al*invAkl*Lk
      ! TkMkJl = invAkl*Al*M*Ak*invAkl*Ll

      CALL compute_products_5(AAlMAl, AAkMAk, AAkMAl, AAlMAk, Jklk, Jkll, &
          TlMkJk, TkMlJl, TlMlJk, TkMkJl)

      ! yet more kinetic terms ... nice rank 1 ones
      ! JkMlk = invAkl*J11*invAkl*Ak*M*Al*invAkl*Lk
      ! JlMkl = invAkl*J11*invAkl*Al*M*Ak*invAkl*Ll
      ! JlMlk = invAkl*J11*invAkl*Al*M*Al*invAkl*Lk
      ! JkMkl = invAkl*J11*invAkl*Ak*M*Ak*invAkl*Ll
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            JkMlk(m, i, j) = invAkl(m, i, 1) * Tklk(m, 1, j)
            JlMkl(m, i, j) = invAkl(m, i, 1) * Tlkl(m, 1, j)
            JlMlk(m, i, j) = invAkl(m, i, 1) * Tllk(m, 1, j)
            JkMkl(m, i, j) = invAkl(m, i, 1) * Tkkl(m, 1, j)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    ! Thats enough for now ... 6 more i,j dependent potential terms later
    ! We have enough to compute Skl dSkl Tkl dTkl, so, on to it ...

    !********************************************************
    ! Compute Skl
    !
    ! Skl = gm1(mk,ml) 2^(3n/2) *
    !    sqrt( (invAkl(1,1)/invAk(1,1))^mk *
    !          (invAkl(1,1)/invAl(1,1))^ml *
    !          (||Lk||||Ll||/|Akl|)^3  )
    !********************************************************
    mko2 = mk / 2
    mlo2 = ml / 2
    DO m = 1, Glob_NumYHYTerms
      SklBuf(m) = (invAkl(m, 1, 1) / invAk(m, 1, 1))**(mko2) * Glob_gm1(mko2, mlo2) &
          *((invAkl(m, 1, 1) / invAl(m, 1, 1))**(mlo2) * Glob_gm1(mko2, mlo2)) &
          *SQRT((ONE * (2**(3 * Glob_n))) * (ABS(detLk(m) * detLl(m) / detAkl(m))**3))
    ENDDO


    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in the matrix gm1
    ! in module constants


    !********************************************************
    ! gradient terms dSk dSl
    !********************************************************
    mkpml = ONE * (mk + ml)

    IF (gradflag > 0) THEN
      ! dSk
      DO j = 1, Glob_n
        DO i = j, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddSk(m, i, j) = (ONE + ONEHALF) * invLk(m, j, i) - THREE * ALk(m, i, j)
          ENDDO
        ENDDO
      ENDDO
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddSk(m, i, j) = ddSk(m, i, j) + mk / invAk(m, 1, 1) * Jkk(m, i, j)
            ddSk(m, i, j) = ddSk(m, i, j) - mkpml / invAkl(m, 1, 1) * Jklk(m, i, j)
            ddSk_over_Skl(m, i, j) = ddSk(m, i, j)
            ddSk(m, i, j) = SklBuf(m) * ddSk(m, i, j)
          ENDDO
        ENDDO
      ENDDO
      ! load it into dSk (a vector)
      indx = 0
      DO j = 1, Glob_n
        DO i = j, Glob_n
          indx = indx + 1
          DO m = 1, Glob_NumYHYTerms
            dSkBuf(m, indx) = ddSk(m, i, j)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    IF (gradflag > 1) THEN
      ! dSl
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            Temp(m, i, j) = ml / invAl(m, 1, 1) * Jll(m, i, j)
            Temp(m, i, j) = Temp(m, i, j) - mkpml / invAkl(m, 1, 1) * Jkll(m, i, j)
            Temp(m, i, j) = Temp(m, i, j) - THREE * ALl(m, i, j)
          ENDDO
        ENDDO
      ENDDO

      ! Multiply by SymMatrix
      ddSl = ZERO
      DO k = 1, Glob_n
        DO j = 1, Glob_n
          DO i = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              ddSl(m, i, j) = ddSl(m, i, j) + SymMatrixBuf(m, i, k) * Temp(m, k, j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, Glob_n
        DO i = j, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddSl(m, i, j) = ddSl(m, i, j) + (ONE + ONEHALF) * invLl(m, j, i)
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddSl_over_Skl(m, i, j) = ddSl(m, i, j)
            ddSl(m, i, j) = SklBuf(m) * ddSl(m, i, j)
          ENDDO
        ENDDO
      ENDDO

      ! Load it into dSl (a vector)
      indx = 0
      DO j = 1, Glob_n
        DO i = j, Glob_n
          indx = indx + 1
          DO m = 1, Glob_NumYHYTerms
            dSlBuf(m, indx) = ddSl(m, i, j)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    !********************************************************
    ! Compute Tkl
    !
    ! Tkl = 2 Skl [ 3 tr(invAkl*Ak*M*Al)
    !       + invAkl(1,1)^-1 { mk*ml*M11/(mk+ml+1)
    !                       - mk (invAklAlMAlinvAkl)(1,1)
    !                       - ml (invAklAkMAkinvAkl)(1,1) } ]
    !********************************************************
    mkml = ONE * mk * ml

    DO m = 1, Glob_NumYHYTerms
      Tkl(m) = mkml * Glob_MassMatrix(1, 1) / (mkpml + ONE)
      Tkl(m) = Tkl(m) - mk * tl(m)
      Tkl(m) = Tkl(m) - ml * tk(m)
      Ttmp(m) = Tkl(m)  ! need this in the gradient
      Tkl(m) = Tkl(m) / invAkl(m, 1, 1)
      Tkl(m) = Tkl(m) + THREE * trT(m)
      Tkl(m) = Tkl(m) * TWO * SklBuf(m)
    ENDDO

    !********************************************************
    ! gradient terms dTk dTl
    !********************************************************
    IF (gradflag > 0) THEN
      DO m = 1, Glob_NumYHYTerms
        tk(m) = invAkl(m, 1, 1)  ! drop this here since we don't need tk anymore
      ENDDO
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddTk(m, i, j) = mk * (JlMlk(m, i, j) + TlMlJk(m, i, j)) &
                - ml * (JkMlk(m, i, j) + TlMkJk(m, i, j))
            ddTk(m, i, j) = TWO * ddTk(m, i, j) / tk(m)
            ddTk(m, i, j) = ddTk(m, i, j) + TWO * Ttmp(m) / (tk(m) * tk(m)) * Jklk(m, i, j)

            ddTk(m, i, j) = ddTk(m, i, j) + SIX * Tllk(m, i, j)
            ddTk(m, i, j) = TWO * SklBuf(m) * ddTk(m, i, j)
            ddTk(m, i, j) = ddTk(m, i, j) + Tkl(m) * ddSk_over_Skl(m, i, j)
          ENDDO
        ENDDO
      ENDDO
      ! Load it into dTk (a vector)
      indx = 0
      DO j = 1, Glob_n
        DO i = j, Glob_n
          indx = indx + 1
          DO m = 1, Glob_NumYHYTerms
            dTk(m, indx) = ddTk(m, i, j)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (gradflag > 1) THEN
      ! dTl
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            Temp(m, i, j) = ml * (JkMkl(m, i, j) + TkMkJl(m, i, j)) &
                - mk * (JlMkl(m, i, j) + TkMlJl(m, i, j))
            Temp(m, i, j) = TWO * (Temp(m, i, j) / tk(m) + Ttmp(m) / (tk(m) * tk(m)) * Jkll(m, i, j))

            Temp(m, i, j) = Temp(m, i, j) + SIX * Tkkl(m, i, j)
            Temp(m, i, j) = TWO * SklBuf(m) * Temp(m, i, j)
          ENDDO
        ENDDO
      ENDDO
      ! Multiply by SymMatrix
      ddTl = ZERO
      DO k = 1, Glob_n
        DO j = 1, Glob_n
          DO i = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              ddTl(m, i, j) = ddTl(m, i, j) + SymMatrixBuf(m, i, k) * Temp(m, k, j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ddTl(m, i, j) = ddTl(m, i, j) + Tkl(m) * ddSl_over_Skl(m, i, j)
          ENDDO
        ENDDO
      ENDDO
      ! load it into dTl (a vector)
      indx = 0
      DO j = 1, Glob_n
        DO i = j, Glob_n
          indx = indx + 1
          DO m = 1, Glob_NumYHYTerms
            dTl(m, indx) = ddTl(m, i, j)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    !********************************************************************
    ! Build intermediate results for potential energy  matrix element
    ! a = invAkl(1,1)
    ! b(i,i) = invAkl(i,i)
    ! b(i,j) = invAkl(i,i) + invAkl(j,j) - 2 invAkl(i,j)
    ! c(i,i) = invAkl(1,i)^2
    ! c(i,j) = invAkl(1,i)^2 + invAkl(1,j)^2 - 2 invAkl(1,i) invAkl(1,j)
    !********************************************************************

    ! a
    DO m = 1, Glob_NumYHYTerms
      a(m) = invAkl(m, 1, 1)
    ENDDO

    ! b
    b = ONE
    DO i = 1, Glob_n
      DO m = 1, Glob_NumYHYTerms
        b(m, i, i) = invAkl(m, i, i)
      ENDDO
    ENDDO
    DO j = 1, Glob_n - 1
      DO i = j + 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          b(m, i, j) = invAkl(m, i, i) + invAkl(m, j, j) - TWO * invAkl(m, i, j)
        ENDDO
      ENDDO
    ENDDO

    ! c
    c = ONE
    DO i = 1, Glob_n
      DO m = 1, Glob_NumYHYTerms
        c(m, i, i) = invAkl(m, 1, i) * invAkl(m, 1, i)
      ENDDO
    ENDDO
    DO j = 1, Glob_n - 1
      DO i = j + 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          c(m, i, j) = (invAkl(m, i, 1) - invAkl(m, j, 1)) * (invAkl(m, i, 1) - invAkl(m, j, 1))
        ENDDO
      ENDDO
    ENDDO

    !**************************************************
    ! Potential energy EVEN POWERS ONLY!!!!!!
    !**************************************************

    ! We will compute all of the 1/rij terms at once
    ! and place them in R, then sum over R*ChargeMatrix to
    ! get Vkl

    ! coab=c/ab  a matrix
    ! compute only lower triangle -SGI
    DO j = 1, Glob_n
      DO i = j, Glob_n
        DO m = 1, Glob_NumYHYTerms
          coab(m, i, j) = c(m, i, j) / (a(m) * b(m, i, j))
        ENDDO
      ENDDO
    ENDDO

    p = (mk + ml) / 2
    R = ZERO

    DO j = 1, Glob_n
      DO k = j, Glob_n
        CALL compute_g3(p, coab(1, k, j), R(1, k, j))
      ENDDO
    ENDDO

    DO m = 1, Glob_NumYHYTerms
      g1(m) = SklBuf(m) * Glob_gm2(p)
    ENDDO
    DO j = 1, Glob_n
      DO k = j, Glob_n
        DO m = 1, Glob_NumYHYTerms
          R(m, k, j) = (R(m, k, j) + ONE) * g1(m) / SQRT(b(m, k, j))
        ENDDO
      ENDDO
    ENDDO

    ! Now compute the potential energy
    Vkl = ZERO
    DO j = 1, Glob_n
      DO i = j, Glob_n
        DO m = 1, Glob_NumYHYTerms
          Vkl(m) = Vkl(m) + R(m, i, j) * Glob_ChargeMatrix(i, j)
        ENDDO
      ENDDO
    ENDDO

    DO m = 1, Glob_NumYHYTerms
      Hkl = Tkl(m) + Vkl(m)
      HklBuf(m) = Hkl
    ENDDO

    !***************************************************************
    ! potential energy gradient elements EVEN powers of r only
    !***************************************************************
    ! Each term of the potential gradient is computed individually and
    ! accumulated in dRk. We sum over these terms multiplied by
    ! the ChargeMatrix element for that term. i=j terms are fundimentaly
    ! different than i/=j terms so we split the loop
    IF (gradflag > 0) THEN
      !********
      ! dVK
      !********

      CALL compute_dVk(SklBuf, invAkl, a, b, coab, p, R, ddSk_over_Skl, ALk, Jklk, dVk)

    ENDIF

    IF (gradflag > 1) THEN
      !********
      ! dVl
      !********

      ! da is independent of i and j
      DO jj = 1, Glob_n
        DO ii = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            da(m, ii, jj) = -TWO * Jkll(m, ii, jj)
          ENDDO
        ENDDO
      ENDDO

      ! diagonal 1/rij terms
      ddVl = ZERO
      DO i = 1, Glob_n
        ! Build the matrices need for this component
        ! RL = invAkl*Jij*invAkl*Ll
        ! RJLl = invAkl*Jij*invAklJ11*invAkl*Ll
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              RLl(m, ii, jj) = invAkl(m, ii, i) * ALl(m, i, jj)
              RJLl(m, ii, jj) = invAkl(m, ii, i) * Jkll(m, i, jj)
            ENDDO
          ENDDO
        ENDDO

        ! JRLl=invAkl*J11*invAkl*Jij*invAkl*Ll
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              JRLl(m, ii, jj) = invAkl(m, ii, 1) * RLl(m, 1, jj)
            ENDDO
          ENDDO
        ENDDO

        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              db(m, ii, jj) = -TWO * RLl(m, ii, jj)
              dc(m, ii, jj) = -TWO * (JRLl(m, ii, jj) + RJLl(m, ii, jj))
            ENDDO
          ENDDO
        ENDDO

        ! now we have the pices to build dRl
        DO m = 1, Glob_NumYHYTerms
          coabij = coab(m, i, i)
          SUM0(m) = ZERO
          g1(m) = ONE - coabij
          g2(m) = ONE
        ENDDO
        DO ii = 1, p
          DO m = 1, Glob_NumYHYTerms
            SUM0(m) = SUM0(m) + Glob_gm3(ii) * ii * g2(m)
            g2(m) = g2(m) * g1(m)
          ENDDO
        ENDDO
        DO m = 1, Glob_NumYHYTerms
          SUM0(m) = SUM0(m) / SQRT(b(m, i, i))
        ENDDO
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              coabij = coab(m, i, i)
              dRl(m, ii, jj) = SUM0(m) &
                  * (coabij * (da(m, ii, jj) * b(m, i, i) &
                  + a(m) * db(m, ii, jj)) - dc(m, ii, jj)) &
                  / (a(m) * b(m, i, i))
              dRl(m, ii, jj) = dRl(m, ii, jj) * SklBuf(m) * Glob_gm2(p) &
                  - R(m, i, i) / (TWO * b(m, i, i)) * db(m, ii, jj)
            ENDDO
          ENDDO
        ENDDO

        ! Multiply by SymMatrix
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              Temp(m, ii, jj) = dRl(m, ii, jj)
              dRl(m, ii, jj) = ZERO
            ENDDO
          ENDDO
        ENDDO
        DO kk = 1, Glob_n
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                dRl(m, ii, jj) = dRl(m, ii, jj) + SymMatrixBuf(m, ii, kk) * Temp(m, kk, jj)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              dRl(m, ii, jj) = dRl(m, ii, jj) + R(m, i, i) * ddSl_over_Skl(m, ii, jj)
              ddVl(m, ii, jj) = ddVl(m, ii, jj) + dRl(m, ii, jj) * Glob_ChargeMatrix(i, i)
            ENDDO
          ENDDO
        ENDDO

      ENDDO

      ! Now do the off diagonal terms
      DO j = 1, Glob_n - 1
        DO i = j + 1, Glob_n
          ! Build the matrices need for this component
          ! RLl = invAkl*Jij*invAkl*Ll
          ! RJLl = invAkl*Jij*invAklJ11*invAkl*Ll
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                RLl(m, ii, jj) = (invAkl(m, ii, i) - invAkl(m, ii, j)) &
                    * (ALl(m, i, jj) - ALl(m, j, jj))
                RJLl(m, ii, jj) = (invAkl(m, ii, i) - invAkl(m, ii, j)) &
                    * (Jkll(m, i, jj) - Jkll(m, j, jj))
              ENDDO
            ENDDO
          ENDDO
          ! JRLl=invAkl*J11*invAkl*Jij*invAkl*Ll
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                JRLl(m, ii, jj) = invAkl(m, ii, 1) * RLl(m, 1, jj)
              ENDDO
            ENDDO
          ENDDO
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                db(m, ii, jj) = -TWO * RLl(m, ii, jj)
                dc(m, ii, jj) = -TWO * (JRLl(m, ii, jj) + RJLl(m, ii, jj))
              ENDDO
            ENDDO
          ENDDO

          ! Now we have the pices to build dRl
          DO m = 1, Glob_NumYHYTerms
            coabij = coab(m, i, j)
            SUM0(m) = ZERO
            g1(m) = ONE - coabij
            g2(m) = ONE
          ENDDO
          DO ii = 1, p
            DO m = 1, Glob_NumYHYTerms
              SUM0(m) = SUM0(m) + Glob_gm3(ii) * ii * g2(m)
              g2(m) = g2(m) * g1(m)
            ENDDO
          ENDDO
          DO m = 1, Glob_NumYHYTerms
            SUM0(m) = SUM0(m) / SQRT(b(m, i, j))
          ENDDO
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                coabij = coab(m, i, j)
                dRl(m, ii, jj) = SUM0(m) &
                    * (coabij * (da(m, ii, jj) * b(m, i, j) &
                    + a(m) * db(m, ii, jj)) - dc(m, ii, jj)) &
                    / (a(m) * b(m, i, j))
                dRl(m, ii, jj) = dRl(m, ii, jj) * SklBuf(m) * Glob_gm2(p) &
                    - R(m, i, j) / (TWO * b(m, i, j)) * db(m, ii, jj)
              ENDDO
            ENDDO
          ENDDO

          ! Multiply by SymMatrix
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                Temp(m, ii, jj) = dRl(m, ii, jj)
                dRl(m, ii, jj) = ZERO
              ENDDO
            ENDDO
          ENDDO
          DO kk = 1, Glob_n
            DO jj = 1, Glob_n
              DO ii = 1, Glob_n
                DO m = 1, Glob_NumYHYTerms
                  dRl(m, ii, jj) = dRl(m, ii, jj) + SymMatrixBuf(m, ii, kk) * Temp(m, kk, jj)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO jj = 1, Glob_n
            DO ii = 1, Glob_n
              DO m = 1, Glob_NumYHYTerms
                dRl(m, ii, jj) = dRl(m, ii, jj) + R(m, i, j) * ddSl_over_Skl(m, ii, jj)
                ddVl(m, ii, jj) = ddVl(m, ii, jj) + dRl(m, ii, jj) * Glob_ChargeMatrix(i, j)
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO
      ! load it into dVk (a vector)
      indx = 0
      DO j = 1, Glob_n
        DO i = j, Glob_n
          indx = indx + 1
          DO m = 1, Glob_NumYHYTerms
            dVl(m, indx) = ddVl(m, i, j)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    SELECT CASE (gradflag)
    CASE (0)
      DO i = 1, Glob_npt
        DO m = 1, Glob_NumYHYTerms
          dHkBuf(m, i) = ZERO
          dHlBuf(m, i) = ZERO
          dSkBuf(m, i) = ZERO
          dSlBuf(m, i) = ZERO
        ENDDO
      ENDDO
    CASE (1)
      DO i = 1, Glob_npt
        DO m = 1, Glob_NumYHYTerms
          dHkBuf(m, i) = dTk(m, i) + dVk(m, i)
          dHlBuf(m, i) = ZERO
          dSlBuf(m, i) = ZERO
        ENDDO
      ENDDO
    CASE (2)
      DO i = 1, Glob_npt
        DO m = 1, Glob_NumYHYTerms
          dHkBuf(m, i) = dTk(m, i) + dVk(m, i)
          dHlBuf(m, i) = dTl(m, i) + dVl(m, i)
        ENDDO
      ENDDO
    END SELECT

  END SUBROUTINE MatrixElementsOpt


  SUBROUTINE PrecomputeMatrices_MA_L_A(np, Nmax, NonlinParam, MassMatrix, Lk, Ak, MAk)
    ! Precomputes, once per basis function, the matrices that depend on that
    ! function alone, so that they need not be rebuilt for every pair :
    !   Lk  :: lower triangular matrix filled from vech[Lk]
    !   Ak  :: Lk*Lk'
    !   MAk :: MassMatrix*Ak (optional). Needs no permutation: the symmetry
    !          projection acts on the ket only, so this is the bra side.
    ! The loops follow the order used in MatrixElementsOpt (vech filled
    ! column by column, sums over k = 1..n).
    ! Input:
    !   np          :: number of nonlinear parameters per function (Glob_npt)
    !   Nmax        :: number of basis functions
    !   NonlinParam :: nonlinear parameters of all functions (np x Nmax)
    !   MassMatrix  :: mass matrix (n x n)
    ! Output:
    !   Lk, Ak, MAk :: (n x n x Nmax)

    IMPLICIT NONE

    INTEGER, PARAMETER :: nn = Glob_n

    INTEGER, INTENT(IN) :: np, Nmax
    REAL(wp), INTENT(IN) :: NonlinParam(np, Nmax), MassMatrix(nn, nn)
    REAL(wp), INTENT(OUT) :: Lk(nn, nn, Nmax), Ak(nn, nn, Nmax)
    REAL(wp), INTENT(OUT), OPTIONAL :: MAk(nn, nn, Nmax)

    INTEGER :: f, i, j, k, indx
    REAL(wp) :: temp1

    DO f = 1, Nmax

      ! Lk from vech[Lk]; the upper triangle is zero
      Lk(1:nn, 1:nn, f) = ZERO
      indx = 0
      DO j = 1, nn
        DO i = j, nn
          indx = indx + 1
          Lk(i, j, f) = NonlinParam(indx, f)
        ENDDO
      ENDDO

      ! Ak = Lk*Lk'
      DO j = 1, nn
        DO i = 1, nn
          temp1 = ZERO
          DO k = 1, nn
            temp1 = temp1 + Lk(i, k, f) * Lk(j, k, f)
          ENDDO
          Ak(i, j, f) = temp1
        ENDDO
      ENDDO

      ! MAk = MassMatrix*Ak
      IF (PRESENT(MAk)) THEN
        DO j = 1, nn
          DO i = 1, nn
            temp1 = ZERO
            DO k = 1, nn
              temp1 = temp1 + MassMatrix(i, k) * Ak(k, j, f)
            ENDDO
            MAk(i, j, f) = temp1
          ENDDO
        ENDDO
      ENDIF

    ENDDO

  END SUBROUTINE PrecomputeMatrices_MA_L_A




  SUBROUTINE compute_g3(p, coab, g3)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p

    REAL(wp), DIMENSION(Glob_NumYHYTerms), INTENT(IN)  :: coab
    REAL(wp), DIMENSION(Glob_NumYHYTerms), INTENT(OUT) :: g3

    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: g1, g2

    INTEGER :: m, i

    ! Tuning parameters:
    INTEGER, PARAMETER :: p3x = 9
    INTEGER, PARAMETER :: p4x = 16


    SELECT CASE (p)

    CASE (:0)

      ! Empty series. This branch is not strictly required today because
      ! the only caller zeroes R before the call loop (see "R = ZERO"
      ! above the compute_g3 loop), but relying on that leaves g3
      ! (INTENT(OUT)) unassigned - the same trap that bit compute_series.

      DO m = 1, Glob_NumYHYTerms
        g3(m) = ZERO
      ENDDO


    CASE (1)

      ! Unwind

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g3(m) = Glob_gm3(1) * g1(m)
      ENDDO


    CASE (2)

      ! Unwind + Horner's Rule

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g3(m) = (Glob_gm3(1) + Glob_gm3(2) * g1(m)) * g1(m)
      ENDDO


    CASE (3)

      ! Unwind + Horner's Rule

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g3(m) = (Glob_gm3(1) + (Glob_gm3(2) + Glob_gm3(3) * g1(m)) * g1(m)) * g1(m)
      ENDDO


    CASE (4:p3x-1)

      ! Unroll 2X

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g2(m) = ONE
        g3(m) = ZERO
      ENDDO

      DO i = 1, p-1, 2
        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+1) * g2(m)
        ENDDO
      ENDDO

      IF (IAND(p, 1) == 1) THEN
        i = p

        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)
        ENDDO
      ENDIF


    CASE (p3x:p4x-1)

      ! Unroll 3X

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g2(m) = ONE
        g3(m) = ZERO
      ENDDO

      DO i = 1, p-2, 3
        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+1) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+2) * g2(m)
        ENDDO
      ENDDO

      DO i = p - MOD(p, 3) + 1, p  ! Let p = 3*q + r, 0 <= r < 3; then p-mod(p,3) = 3*q
        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)
        ENDDO
      ENDDO


    CASE (p4x:)

      ! Unroll 4X

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE - coab(m)
        g2(m) = ONE
        g3(m) = ZERO
      ENDDO

      DO i = 1, p-3, 4
        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+1) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+2) * g2(m)

          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i+3) * g2(m)
        ENDDO
      ENDDO

      ! do i=p-mod(p,4)+1,p
      DO i = p - IAND(p, 3) + 1, p
        DO m = 1, Glob_NumYHYTerms
          g2(m) = g2(m) * g1(m)
          g3(m) = g3(m) + Glob_gm3(i) * g2(m)
        ENDDO
      ENDDO

    END SELECT

  END SUBROUTINE compute_g3


  SUBROUTINE compute_products_1(Lk, Ak, Al, invAkl, PLl, invAl, &
                                ALk, ALl, AAk, AAl, AlLl)
    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Lk, Ak, Al, invAkl, PLl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: invAl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ALk, ALl, AAk, AAl, AlLl

    INTEGER :: m, i, j, k


    ALk = ZERO
    ALl = ZERO
    AAk = ZERO
    AAl = ZERO
    AlLl = ZERO

    ! CDIR$ UNROLL (2)
    DO k = 1, Glob_n
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            ALk(m, i, j) = ALk(m, i, j) + invAkl(m, i, k) * Lk(m, k, j)
            ALl(m, i, j) = ALl(m, i, j) + invAkl(m, i, k) * PLl(m, k, j)
            AlLl(m, i, j) = AlLl(m, i, j) + invAl(m, i, k) * PLl(m, k, j)
            AAk(m, i, j) = AAk(m, i, j) + invAkl(m, i, k) * Ak(m, k, j)
            AAl(m, i, j) = AAl(m, i, j) + invAkl(m, i, k) * Al(m, k, j)  ! could do I-invAkl*Ak
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_products_1


  SUBROUTINE compute_products_2(Al, MAl)
    ! MAl = M*Al for the ket; the bra's M*Ak is precomputed by
    ! PrecomputeMatrices_MA_L_A

    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Al
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: MAl

    INTEGER :: m, i, j, k


    MAl = ZERO

    ! CDIR$ UNROLL (2)
    DO k = 1, Glob_n
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            MAl(m, i, j) = MAl(m, i, j) + Glob_MassMatrix(i, k) * Al(m, k, j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_products_2


  SUBROUTINE compute_products_3(AAk, AAl, MAk, MAl, &
                                AAlMAl, AAkMAk, AAkMAl, AAlMAk)

    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: AAk, AAl, MAk, MAl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: AAlMAl, AAkMAk, AAkMAl, AAlMAk

    INTEGER :: m, i, j, k


    AAlMAl = ZERO
    AAkMAk = ZERO
    AAkMAl = ZERO
    AAlMAk = ZERO

    ! CDIR$ UNROLL (2)
    DO k = 1, Glob_n
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            AAlMAl(m, i, j) = AAlMAl(m, i, j) + AAl(m, i, k) * MAl(m, k, j)
            AAkMAk(m, i, j) = AAkMAk(m, i, j) + AAk(m, i, k) * MAk(m, k, j)
            AAkMAl(m, i, j) = AAkMAl(m, i, j) + AAk(m, i, k) * MAl(m, k, j)
            AAlMAk(m, i, j) = AAlMAk(m, i, j) + AAl(m, i, k) * MAk(m, k, j)  ! MAk-invAkl*Ak...
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_products_3


  SUBROUTINE compute_products_4(ALk, ALl, AAlMAl, AAkMAk, AAkMAl, AAlMAk, &
                                Tllk, Tkkl, Tklk, Tlkl)

    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ALk, ALl
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: AAlMAl, AAkMAk, AAkMAl, AAlMAk
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Tllk, Tkkl, Tklk, Tlkl

    INTEGER :: m, i, j, k


    Tllk = ZERO
    Tkkl = ZERO
    Tklk = ZERO
    Tlkl = ZERO

    ! CDIR$ UNROLL (2)
    DO k = 1, Glob_n
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            Tllk(m, i, j) = Tllk(m, i, j) + AAlMAl(m, i, k) * ALk(m, k, j)
            Tkkl(m, i, j) = Tkkl(m, i, j) + AAkMAk(m, i, k) * ALl(m, k, j)
            Tklk(m, i, j) = Tklk(m, i, j) + AAkMAl(m, i, k) * ALk(m, k, j)
            Tlkl(m, i, j) = Tlkl(m, i, j) + AAlMAk(m, i, k) * ALl(m, k, j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_products_4


  SUBROUTINE compute_products_5(AAlMAl, AAkMAk, AAkMAl, AAlMAk, Jklk, Jkll, &
                                TlMkJk, TkMlJl, TlMlJk, TkMkJl)

    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Jklk, Jkll
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: AAlMAl, AAkMAk, AAkMAl, AAlMAk
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: TlMkJk, TkMlJl, TlMlJk, TkMkJl

    INTEGER :: m, i, j, k


    TlMkJk = ZERO
    TkMlJl = ZERO
    TlMlJk = ZERO
    TkMkJl = ZERO

    ! DIR$ UNROLL(2)
    DO k = 1, Glob_n
      DO j = 1, Glob_n
        DO i = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            TlMkJk(m, i, j) = TlMkJk(m, i, j) + AAlMAk(m, i, k) * Jklk(m, k, j)
            TkMlJl(m, i, j) = TkMlJl(m, i, j) + AAkMAl(m, i, k) * Jkll(m, k, j)
            TlMlJk(m, i, j) = TlMlJk(m, i, j) + AAlMAl(m, i, k) * Jklk(m, k, j)
            TkMkJl(m, i, j) = TkMkJl(m, i, j) + AAkMAk(m, i, k) * Jkll(m, k, j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_products_5


  SUBROUTINE compute_dVk(SklBuf, invAkl, a, b, coab, p, R, ddSk_over_Skl, ALk, Jklk, dVk)

    IMPLICIT NONE


    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: SklBuf  ! input

    REAL(wp), DIMENSION(Glob_NumYHYTerms)                 :: SUM0    ! local
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: invAkl  ! input
    REAL(wp), DIMENSION(Glob_NumYHYTerms)                 :: a       ! input

    REAL(wp)                                              :: coabij      ! local
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: R, b, coab  ! input

    INTEGER :: p  ! input

    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_npt)       :: dVk            ! output
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ddSk_over_Skl  ! input

    REAL(wp)                                              :: dRk   ! local
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ddVk  ! local

    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: ALk   ! input
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: Jklk  ! input

    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: da, db, dc       ! local
    REAL(wp), DIMENSION(Glob_NumYHYTerms, Glob_n, Glob_n) :: RLk, JRLk, RJLk  ! local

    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: g1, g2            ! local
    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: factor1, factor2

    INTEGER :: m, i, j, ii, jj, indx


    ! da is independent of i and j
    DO jj = 1, Glob_n
      DO ii = 1, Glob_n
        DO m = 1, Glob_NumYHYTerms
          da(m, ii, jj) = -TWO * Jklk(m, ii, jj)
        ENDDO
      ENDDO
    ENDDO


    ! diagonal 1/rij terms
    ddVk = ZERO

    DO i = 1, Glob_n

      ! Build the matrices need for this component
      ! RLk= nvAkl*Jij*invAkl*Lk
      ! RJLk=invAkl*Jij*invAklJ11*invAkl*Lk
      DO jj = 1, Glob_n
        DO ii = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            RLk(m, ii, jj) = invAkl(m, ii, i) * ALk(m, i, jj)
            RJLk(m, ii, jj) = invAkl(m, ii, i) * Jklk(m, i, jj)
          ENDDO
        ENDDO
      ENDDO


      ! JRLk=invAkl*J11*invAkl*Jij*invAkl*Lk
      DO jj = 1, Glob_n
        DO ii = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            JRLk(m, ii, jj) = invAkl(m, ii, 1) * RLk(m, 1, jj)
          ENDDO
        ENDDO
      ENDDO


      DO jj = 1, Glob_n
        DO ii = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            db(m, ii, jj) = -TWO * RLk(m, ii, jj)
            dc(m, ii, jj) = -TWO * (JRLk(m, ii, jj) + RJLk(m, ii, jj))
          ENDDO
        ENDDO
      ENDDO


      ! Now we have the pices to build dRk
      CALL compute_series(p, coab(1, i, i), SUM0)

      DO m = 1, Glob_NumYHYTerms
        SUM0(m) = SUM0(m) / SQRT(b(m, i, i))
      ENDDO

      DO m = 1, Glob_NumYHYTerms
        factor1(m) = SUM0(m) / (a(m) * b(m, i, i)) * SklBuf(m) * Glob_gm2(p)
        factor2(m) = R(m, i, i) / (TWO * b(m, i, i))
      ENDDO

      DO jj = 1, Glob_n
        DO ii = 1, Glob_n
          DO m = 1, Glob_NumYHYTerms
            coabij = coab(m, i, i)

            dRk = coabij * (da(m, ii, jj) * b(m, i, i) + a(m) * db(m, ii, jj)) - dc(m, ii, jj)
            dRk = dRk * factor1(m) - factor2(m) * db(m, ii, jj) + R(m, i, i) * ddSk_over_Skl(m, ii, jj)

            ddVk(m, ii, jj) = ddVk(m, ii, jj) + dRk * Glob_ChargeMatrix(i, i)
          ENDDO
        ENDDO
      ENDDO

    ENDDO


    ! Now do the off diagonal terms
    DO j = 1, Glob_n - 1
      DO i = j + 1, Glob_n

        ! Build the matrices need for this component
        ! RLk=invAkl*Jij*invAkl*Lk
        ! RJLk=invAkl*Jij*invAklJ11*invAkl*Lk
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              RLk(m, ii, jj) = (invAkl(m, ii, i) - invAkl(m, ii, j)) * &
                               (ALk(m, i, jj) - ALk(m, j, jj))

              RJLk(m, ii, jj) = (invAkl(m, ii, i) - invAkl(m, ii, j)) * &
                               (Jklk(m, i, jj) - Jklk(m, j, jj))
            ENDDO
          ENDDO
        ENDDO


        ! JRLk=invAkl*J11*invAkl*Jij*invAkl*Lk
        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              JRLk(m, ii, jj) = invAkl(m, ii, 1) * RLk(m, 1, jj)
            ENDDO
          ENDDO
        ENDDO


        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              db(m, ii, jj) = -TWO * RLk(m, ii, jj)
              dc(m, ii, jj) = -TWO * (JRLk(m, ii, jj) + RJLk(m, ii, jj))
            ENDDO
          ENDDO
        ENDDO


        ! Now we have the pices to build dRk
        CALL compute_series(p, coab(1, i, j), SUM0)

        DO m = 1, Glob_NumYHYTerms
          SUM0(m) = SUM0(m) / SQRT(b(m, i, j))
        ENDDO

        DO m = 1, Glob_NumYHYTerms
          factor1(m) = SUM0(m) / (a(m) * b(m, i, j)) * SklBuf(m) * Glob_gm2(p)
          factor2(m) = R(m, i, j) / (TWO * b(m, i, j))
        ENDDO

        DO jj = 1, Glob_n
          DO ii = 1, Glob_n
            DO m = 1, Glob_NumYHYTerms
              coabij = coab(m, i, j)

              dRk = coabij * (da(m, ii, jj) * b(m, i, j) + a(m) * db(m, ii, jj)) - dc(m, ii, jj)
              dRk = dRk * factor1(m) - factor2(m) * db(m, ii, jj) + R(m, i, j) * ddSk_over_Skl(m, ii, jj)

              ddVk(m, ii, jj) = ddVk(m, ii, jj) + dRk * Glob_ChargeMatrix(i, j)
            ENDDO
          ENDDO
        ENDDO

      ENDDO
    ENDDO


    ! Load it into dVk (a vector)
    indx = 0

    DO j = 1, Glob_n
      DO i = j, Glob_n
        indx = indx + 1

        DO m = 1, Glob_NumYHYTerms
          dVk(m, indx) = ddVk(m, i, j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE compute_dVk


  SUBROUTINE compute_series(p, coab, sum)

    IMPLICIT NONE

    INTEGER, INTENT(IN)                                   :: p
    REAL(wp), DIMENSION(Glob_NumYHYTerms), INTENT(IN)  :: coab
    REAL(wp), DIMENSION(Glob_NumYHYTerms), INTENT(OUT) :: sum

    REAL(wp), DIMENSION(Glob_NumYHYTerms) :: g1, g2

    INTEGER :: m, i

    SELECT CASE (p)

    CASE (:0)

      ! For p=0 the series is empty, so its value is zero. Without this
      ! branch sum (INTENT(OUT)) is never assigned and the caller keeps
      ! whatever its actual argument held on entry - in compute_dVk that
      ! is the local SUM0, which a debug build fills with NaN. The NaN
      ! then reaches DMNG and trips -ffpe-trap=invalid inside DV2NRM.
      DO m = 1, Glob_NumYHYTerms
        sum(m) = ZERO
      ENDDO

    CASE (1)

      DO m = 1, Glob_NumYHYTerms
        sum(m) = Glob_gm3(1)
      ENDDO

    CASE (2)

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE-coab(m)
        sum(m) = Glob_gm3(1)+Glob_gm3(2)*2*g1(m)
      ENDDO

    CASE (3)

      DO m = 1, Glob_NumYHYTerms
        g1(m) = ONE-coab(m)
        sum(m) = Glob_gm3(1)+(Glob_gm3(2)*2+Glob_gm3(3)*3*g1(m))*g1(m)
      ENDDO

    CASE (4:)

      DO m = 1, Glob_NumYHYTerms
        sum(m) = ZERO
        g1(m) = ONE-coab(m)
        g2(m) = ONE
      ENDDO

      DO i = 1, p-1, 2
        DO m = 1, Glob_NumYHYTerms
          sum(m) = sum(m)+Glob_gm3(i)*(i)*g2(m)
          g2(m) = g2(m)*g1(m)
          sum(m) = sum(m)+Glob_gm3(i+1)*(i+1)*g2(m)
          g2(m) = g2(m)*g1(m)
        ENDDO
      ENDDO

      IF (IAND(p, 1) == 1) THEN
        i = p
        DO m = 1, Glob_NumYHYTerms
          sum(m) = sum(m)+Glob_gm3(i)*i*g2(m)
        ENDDO
      ENDIF

    END SELECT

  END SUBROUTINE compute_series


  SUBROUTINE ExpValuesMatElem(mk, PreLk, PreAk, PreMAk, ml, vechLl, SymMatrix, &
            Skl, Tkl, Vkl, rkl_m1, rkl, rkl_2, deltarkl, &
            MassVelocity_kl, Darwin_kl, Darwin1_kl, OrbitOrbit_kl, &
            NumCFGridPoints, CFGrid, CFkl, AreCorrFuncNeeded)

    IMPLICIT NONE

    ! Symmetry adapted matrix elements of the operators needed for the
    ! expectation values (overlap, kinetic and potential energy, 1/r, r,
    ! r^2, delta, mass-velocity, Darwin, orbit-orbit, correlation function).
    ! Input:
    !   PreLk, PreAk, PreMAk :: bra side, precomputed by PrecomputeMatrices_MA_L_A:
    !                           Lk, Ak = Lk*Lk' and M*Ak (no permutation needed)
    !   vechLl         :: n(n+1)/2 exponent parameters of the ket (Ll lower triangular)
    !   mk, ml         :: premultiplier powers
    !   SymMatrix      :: the symmetry permutation matrix
    !   NumCFGridPoints   :: number of grid points for the correlation function
    !   CFGrid            :: grid points where the correlation function is computed
    !   AreCorrFuncNeeded :: whether CFkl is computed (CFGrid and CFkl are
    !                        not referenced otherwise)
    ! Output (normalized):
    !   Skl, Tkl, Vkl  :: overlap, kinetic and potential energy terms,
    !                     followed by the remaining operator terms
    !   CFkl           :: matrix elements of the nucleus-nucleus correlation
    !                     function (pair 0-1, i.e. r_1), the only one available
    !                     in PG_0S - the g_1 of the reference codes

    ! Parameters
    INTEGER :: n

    INTEGER, INTENT(IN)      :: mk, ml
    REAL(wp), INTENT(IN)  :: PreLk(Glob_n, Glob_n), PreAk(Glob_n, Glob_n), PreMAk(Glob_n, Glob_n)
    REAL(wp), INTENT(IN)  :: vechLl(Glob_npt)
    REAL(wp), INTENT(IN)  :: SymMatrix(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: Skl, Tkl, Vkl, rkl_m1(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: rkl(Glob_n, Glob_n), rkl_2(Glob_n, Glob_n), deltarkl(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: MassVelocity_kl, Darwin_kl, Darwin1_kl, OrbitOrbit_kl
    INTEGER, INTENT(IN)   :: NumCFGridPoints
    REAL(wp), INTENT(IN)  :: CFGrid(NumCFGridPoints)
    REAL(wp), INTENT(OUT) :: CFkl(NumCFGridPoints)
    LOGICAL, INTENT(IN)   :: AreCorrFuncNeeded

    ! Local variables
    INTEGER     :: i, j, k, indx, info, jj, ii, kk, mko2, mlo2
    REAL(wp) :: detLk, detLl, detAkl
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3
    REAL(wp) :: Lk(Glob_n, Glob_n), Ll(Glob_n, Glob_n), Ak(Glob_n, Glob_n), Al(Glob_n, Glob_n)
    REAL(wp) :: invAkl(Glob_n, Glob_n), PLl(Glob_n, Glob_n), invL(Glob_n, Glob_n)
    REAL(wp) :: Temp(Glob_n, Glob_n)
    REAL(wp) :: trT, tl, tk, mkml, mkpml, Ttmp
    REAL(wp) :: a, coabij
    REAL(wp) :: R(Glob_n, Glob_n), b(Glob_n, Glob_n), c(Glob_n, Glob_n), coab(Glob_n, Glob_n), Vtmp(Glob_n, Glob_n)
    INTEGER     :: p, q, t
    REAL(wp) :: invLk(Glob_n, Glob_n), invLl(Glob_n, Glob_n), invAk(Glob_n, Glob_n), invAl(Glob_n, Glob_n), chAl(Glob_n, Glob_n)
    REAL(wp) :: r1kl, r1kl_2
    REAL(wp) :: ALk(Glob_n, Glob_n), ALl(Glob_n, Glob_n), AAk(Glob_n, Glob_n), AAl(Glob_n, Glob_n), AlLl(Glob_n, Glob_n)
    REAL(wp) :: MAk(Glob_n, Glob_n), MAl(Glob_n, Glob_n)
    REAL(wp) :: Jkk(Glob_n, Glob_n), Jll(Glob_n, Glob_n), Jklk(Glob_n, Glob_n), Jkll(Glob_n, Glob_n)
    REAL(wp) :: AAlMAl(Glob_n, Glob_n), AAkMAk(Glob_n, Glob_n), AAkMAl(Glob_n, Glob_n), AAlMAk(Glob_n, Glob_n)
    REAL(wp) :: Tllk(Glob_n, Glob_n), Tkkl(Glob_n, Glob_n), Tklk(Glob_n, Glob_n), Tlkl(Glob_n, Glob_n)
    REAL(wp) :: TlMkJk(Glob_n, Glob_n), TkMlJl(Glob_n, Glob_n), TlMlJk(Glob_n, Glob_n), TkMkJl(Glob_n, Glob_n)
    REAL(wp) :: JkMlk(Glob_n, Glob_n), JlMkl(Glob_n, Glob_n), JlMlk(Glob_n, Glob_n), JkMkl(Glob_n, Glob_n)
    REAL(wp) :: da(Glob_n, Glob_n), db(Glob_n, Glob_n), dc(Glob_n, Glob_n)
    REAL(wp) :: RLk(Glob_n, Glob_n), RLl(Glob_n, Glob_n)
    REAL(wp) :: JRLk(Glob_n, Glob_n), JRLl(Glob_n, Glob_n), RJLk(Glob_n, Glob_n), RJLl(Glob_n, Glob_n)
    REAL(wp) :: g1, g2, g3

    ! DK
    REAL(wp) :: cvel  ! speed of light

    ! DK

    n = Glob_n

    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Bra side: Lk, Ak = Lk*Lk' and MAk = M*Ak are precomputed once per function
    Lk = PreLk
    Ak = PreAk
    MAk = PreMAk

    ! Build the lower triangular matrix Ll from the elements in vechLl
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Ll(i, j) = vechLl(indx)
      ENDDO
    ENDDO

    ! Apply the symmetry transformation SymMatrix to Ll  PLl=SymMatrix'*Ll
    PLl = ZERO
    DO j = 1, n
      DO k = j, n
        DO i = 1, n
          PLl(i, j) = PLl(i, j)+SymMatrix(k, i)*Ll(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! Build AL=P'*LL*LL'*P (Ak is precomputed)
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Al(i, j) = Al(i, j)+PLl(i, k)*PLl(j, k)
        ENDDO
      ENDDO
    ENDDO

    ! Compute Akl the upper triangle of Akl=Ak+Al is in invL
    ! (we store the upper part of Akl for convience in the
    ! Cholesky decomposition code below)
    DO j = 1, n
      DO i = 1, j-1
        invL(j, i) = ZERO
      ENDDO
      DO i = j, n
        invL(j, i) = Ak(i, j)+Al(i, j)
      ENDDO
    ENDDO

    ! compute det(LK) and det(LL)
    detLk = ONE
    detLl = ONE
    DO i = 1, n
      detLk = detLk*Lk(i, i)
      detLl = detLl*Ll(i, i)
    ENDDO

    ! Compute the Cholesky decomposition of Akl (invL)
    ! This code is adapted from Numerical Recipes  choldc.for
    ! det(Akl) is computed during the decomposition and stored in
    ! detAkl
    !
    ! need to DO the decomp on Al too to get invAl since
    ! Al is P'Ll*(P'Ll)' i.e. we lost the original chol decomp
    ! we'll put it in chAl
    chAl = Al
    detAkl = ONE
    DO i = 1, n
      DO j = i, n
        SUM0 = invL(i, j)
        SUM1 = chAl(i, j)
        DO k = i-1, 1, -1
          SUM0 = SUM0-invL(i, k)*invL(j, k)
          SUM1 = SUM1-chAl(i, k)*chAl(j, k)
        ENDDO
        IF (i == j) THEN
          invL(i, i) = SQRT(SUM0)
          chAl(i, i) = SQRT(SUM1)  ! don't need the det for Al
          detAkl = detAkl*SUM0
        ELSE
          invL(j, i) = SUM0/invL(i, i)
          chAl(j, i) = SUM1/chAl(i, i)
        ENDIF
      ENDDO
    ENDDO

    ! invL now CONTAINS the Cholesky decomp,L, of Akl, where Akl=L*L'
    ! detAkl now CONTAINS det(Akl)
    ! chAl now CONTAINS the Cholesky decomp of Al
    !
    ! Compute the inverse of Akl, Ak, Al, Lk, Ll
    !
    ! Start with invL, invLk, invLl, chAl
    invLk = Lk
    invLl = Ll
    ! [invert the lower triangle matrices]
    DO i = 1, n
      invL(i, i) = ONE/invL(i, i)
      invLk(i, i) = ONE/invLk(i, i)
      invLl(i, i) = ONE/invLl(i, i)
      chAl(i, i) = ONE/chAl(i, i)
      DO j = i+1, n
        SUM0 = ZERO
        SUM1 = ZERO
        SUM2 = ZERO
        SUM3 = ZERO
        DO k = i, j-1
          SUM0 = SUM0-invL(j, k)*invL(k, i)
          SUM1 = SUM1-invLk(j, k)*invLk(k, i)
          SUM2 = SUM2-invLl(j, k)*invLl(k, i)
          SUM3 = SUM3-chAl(j, k)*chAl(k, i)
        ENDDO
        invL(j, i) = SUM0/invL(j, j)
        invL(i, j) = ZERO  ! zero out the upper triangle for convenience
        invLk(j, i) = SUM1/invLk(j, j)
        invLk(i, j) = ZERO
        invLl(j, i) = SUM2/invLl(j, j)
        invLl(i, j) = ZERO
        chAl(j, i) = SUM3/chAl(j, j)
        chAl(i, j) = ZERO
      ENDDO
    ENDDO

    ! We now have invLk, invLl (this is not invP'Ll), inv(chAl) invL
    !
    ! Now DO invAkl=invL'*invL
    !       invAk=invLk'*invLk
    !       invAl=chAl'*chAl    NOTE that this has the SymMatrix term in it
    invAkl = ZERO
    invAk = ZERO
    invAl = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          invAkl(i, j) = invAkl(i, j)+ invL(k, i)*invL(k, j)
          invAk(i, j) = invAk(i, j)+invLk(k, i)*invLk(k, j)
          invAl(i, j) = invAl(i, j)+chAl(k, i)*chAl(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! OK looks like we got most of the important terms
    ! Lk, Ll, PLl, Ak, Al=P'Ll(P'Ll)', invAkl, invAk, invAl
    ! invLk, invLl [not inv(P'Ll)]


    !********************************************************
    ! Compute intermediate matrix results LOTS OF THEM!
    !********************************************************

    ! ALk=invAkl*Lk
    ! ALl=invAkl*PLl
    ! AlLl=invAl*PLl need this because of sym
    ! AAk=invAkl*Ak
    ! AAl=invAkl*Al
    ALk = ZERO
    ALl = ZERO
    AAk = ZERO
    AAl = ZERO
    AlLl = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          ALk(i, j) = ALk(i, j) +invAkl(i, k)*Lk(k, j)
          ALl(i, j) = ALl(i, j) +invAkl(i, k)*PLl(k, j)
          AlLl(i, j) = AlLl(i, j)+invAl(i, k)*PLl(k, j)
          AAk(i, j) = AAk(i, j) +invAkl(i, k)*Ak(k, j)
          AAl(i, j) = AAl(i, j) +invAkl(i, k)*Al(k, j)  ! could DO I-invAkl*Ak
        ENDDO
      ENDDO
    ENDDO

    DO j = 1, n
      Jkk(1:n, j) = invAk(1:n, 1)*invLk(j, 1)
      Jll(1:n, j) = invAl(1:n, 1)*AlLl(1, j)  ! using AlLl to get sym right
      Jklk(1:n, j) = invAkl(1:n, 1)*ALk(1, j)
      Jkll(1:n, j) = invAkl(1:n, 1)*ALl(1, j)
    ENDDO

    ! MAl=M*Al (MAk=M*Ak is precomputed)
    MAl = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          MAl(i, j) = MAl(i, j)+Glob_MassMatrix(i, k)*Al(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! AAlMAl=invAkl*Al*M*Al
    ! AAkMAk=invAkl*Ak*M*Ak
    ! AAkMAl=invAkl*Ak*M*Al
    ! AAlMAk=invAkl*Al*M*Ak
    AAlMAl = ZERO
    AAkMAk = ZERO
    AAkMAl = ZERO
    AAlMAk = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          AAlMAl(i, j) = AAlMAl(i, j)+AAl(i, k)*MAl(k, j)
          AAkMAk(i, j) = AAkMAk(i, j)+AAk(i, k)*MAk(k, j)
          AAkMAl(i, j) = AAkMAl(i, j)+AAk(i, k)*MAl(k, j)
          AAlMAk(i, j) = AAlMAk(i, j)+AAl(i, k)*MAk(k, j)  ! MAk-invAkl*Ak...
        ENDDO
      ENDDO
    ENDDO


    ! A few kinetic energy terms (scalers)
    ! trT=tr[invAkl*Ak*M*Al]
    ! tl=(invAkl*Al*M*Al*invAkl)(1,1)
    ! tk=(invAkl*Ak*M*Ak*invAkl)(1,1)
    trT = ZERO
    tl = ZERO
    tk = ZERO
    DO i = 1, n
      trT = trT+AAkMAl(i, i)
      tl = tl+AAlMAl(1, i)*invAkl(i, 1)
      tk = tk+AAkMAk(1, i)*invAkl(i, 1)
    ENDDO

    ! Thats enough for now ... 6 more i,j dependent potential terms later
    ! We have enough to compute Skl dSkl Tkl dTkl, so, on to it ...

    !********************************************************
    ! Compute Skl
    !
    ! Skl = Glob_gm1(mk,ml) 2^(3n/2) *
    !    SQRT( (invAkl(1,1)/invAk(1,1))^mk *
    !          (invAkl(1,1)/invAl(1,1))^ml *
    !          (||Lk||||Ll||/|Akl|)^3  )
    !********************************************************
    mko2 = mk/2
    mlo2 = ml/2
    Skl = (invAkl(1, 1)/invAk(1, 1))**(mko2)*Glob_gm1(mko2, mlo2) &
        *((invAkl(1, 1)/invAl(1, 1))**(mlo2) * Glob_gm1(mko2, mlo2)) &
        *SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))


    ! Glob_gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             SQRT( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in the matrix Glob_gm1
    ! in module constants

    r1kl = Skl*((mk+ml+2)/TWO)*Glob_gm2((mk+ml)/2)*SQRT(invAkl(1, 1))
    r1kl_2 = Skl*((mk+ml+3)/TWO)*invAkl(1, 1)


    !********************************************************
    ! Compute Tkl
    !
    ! Tkl = 2 Skl [ 3 tr(invAkl*Ak*M*Al)
    !       + invAkl(1,1)^-1 { mk*ml*M11/(mk+ml+1)
    !                       - mk (invAklAlMAlinvAkl)(1,1)
    !                       - ml (invAklAkMAkinvAkl)(1,1) } ]
    !********************************************************
    mkpml = ONE*(mk+ml)
    mkml = ONE*mk*ml

    Tkl = mkml*Glob_MassMatrix(1, 1)/(mkpml+ONE)
    Tkl = Tkl-mk*tl
    Tkl = Tkl-ml*tk
    Ttmp = Tkl  ! need this in the gradient
    Tkl = Tkl/invAkl(1, 1)
    Tkl = Tkl+THREE*trT
    Tkl = Tkl*TWO*Skl

    ! a
    a = invAkl(1, 1)

    ! b
    b = ONE
    DO i = 1, n
      b(i, i) = invAkl(i, i)
    ENDDO
    DO j = 1, n-1
      DO i = j+1, n
        b(i, j) = invAkl(i, i)+invAkl(j, j)-TWO*invAkl(i, j)
      ENDDO
    ENDDO

    ! c
    c = ONE
    DO i = 1, n
      c(i, i) = invAkl(1, i)*invAkl(1, i)
    ENDDO
    DO j = 1, n-1
      DO i = j+1, n
        c(i, j) = (invAkl(i, 1)-invAkl(j, 1))*(invAkl(i, 1)-invAkl(j, 1))
      ENDDO
    ENDDO

    !**************************************************
    ! Potential energy EVEN POWERS ONLY!!!!!!
    !**************************************************

    ! We will compute all of the 1/rij terms at once
    ! and place them in R, THEN sum over R*Glob_ChargeMatrix to
    ! get Vkl

    ! coab=c/ab  a matrix
    coab = c/(a*b)

    p = (mk+ml)/2
    R = ZERO

    DO j = 1, n
      DO k = j, n
        g1 = ONE-coab(k, j)
        g2 = ONE
        g3 = ZERO
        DO i = 1, p
          g2 = g2*g1
          g3 = g3+Glob_gm3(i)*g2
        ENDDO
        R(k, j) = g3
      ENDDO
    ENDDO

    g1 = Skl*Glob_gm2(p)
    DO j = 1, n
      DO k = j, n
        R(k, j) = (R(k, j)+ONE)*g1/SQRT(b(k, j))
      ENDDO
    ENDDO


    ! Now compute the potential energy
    Vkl = ZERO
    DO j = 1, n
      DO i = j, n
        Vkl = Vkl+R(i, j)*Glob_ChargeMatrix(i, j)
      ENDDO
    ENDDO

    ! Distances
    !!!
    DO j = 1, n
      DO i = j, n
        rkl_m1(i, j) = R(i, j)
        rkl_m1(j, i) = R(i, j)
      ENDDO
    ENDDO

    rkl = ZERO
    rkl_2 = ZERO
    DO j = 1, n
      DO k = j, n
        g1 = ONE-coab(k, j)
        g2 = ONE
        g3 = ZERO
        DO i = 1, p
          g2 = g2*g1
          g3 = g3+Glob_gm3(i)*i*g2
        ENDDO
        rkl(k, j) = R(k, j)*(p*c(k, j)/a + b(k, j)) - Skl*(Glob_gm2(p)/SQRT(b(k, j)))*g3*(c(k, j)/a)
        rkl(j, k) = rkl(k, j)
      ENDDO
    ENDDO
    ! rkl(1,1)=Skl*(p+ONE)*Glob_gm2(p)*SQRT(invAkl(1,1))
    DO j = 1, n
      DO k = j, n
        rkl_2(k, j) = Skl*(p*c(k, j)/a + 3*b(k, j)/TWO)
        rkl_2(j, k) = rkl_2(k, j)
      ENDDO
    ENDDO
    ! rkl_2(1,1)=Skl*(p+THREE/TWO)*invAkl(1,1)

    ! Evaluation of the nucleus-nucleus correlation function (pair 0-1)
    IF (AreCorrFuncNeeded) THEN
      g1 = 1/invAkl(1, 1)
      DO i = 1, NumCFGridPoints
        g2 = g1*CFGrid(i)*CFGrid(i)
        IF (g2 == ZERO) THEN
          IF (p == 0) THEN
            CFkl(i) = Skl*g1*SQRT(g1)*ONEHALF/Glob_Pi
          ELSE
            CFkl(i) = ZERO
          ENDIF
        ELSE
          CFkl(i) = EXP(-g2+p*LOG(g2)-Glob_lngamma(p))*Skl*g1*SQRT(g1)*ONEHALF/Glob_Pi
        ENDIF
      ENDDO
    ENDIF

    !!!
    ! Computing expectation values of delta functions
    DO i = 1, n
      IF (i == 1) THEN
        IF (p == 0) THEN
          deltarkl(1, 1) = delta_R1(n, Skl, invAkl, detAkl)
        ELSE
          deltarkl(1, 1) = ZERO
        ENDIF
      ELSE
        deltarkl(i, i) = delta_Ri(n, mk, ml, i, Ak, Al, Skl, invAkl, detAkl)
      ENDIF
      DO j = 1, i-1
        deltarkl(i, j) = delta_Rij(n, mk, ml, j, i, Ak, Al, Skl, invAkl, detAkl)
        deltarkl(j, i) = deltarkl(i, j)
      ENDDO
    ENDDO


    !!!
    ! cvel = 137.0359895_wp
    cvel = ONE
    CALL Darwin(n, i, j, mk, ml, cvel, Skl, Ak, Al, invAkl, Darwin_kl)

    Darwin1_kl = Darwin_kl

    OrbitOrbit_kl = At_Orb_Orb1(Ak, Al, Skl, invAkl, rkl_m1)


    ! compute mass-velocity correction
    !
    ! MassVelocity_kl=-0.125/c^2{
    !               1/ M_1^3 <\nabla' J \nabla \phi_k|\nabla' J \nabla \phi_l>
    !               + sum_{i=1}^{n} 1/M(i+1)^3  <\nabla' Jii \nabla \phi_k|\nabla' Jii\nabla \phi_l>
    ! where (J)_kl=(1)_kl, (Jii)_kl=\delta_{ki}\delta{il}


    !!!
    ! cvel=137.0359895_wp
    cvel = ONE

    MassVelocity_kl = ZERO

    MassVelocity_kl = -0.125_wp/(cvel*cvel*Glob_Mass(1)**3)*&
            nJnnJn_kl(n, mk, ml, Skl, Ak, Al, invAkl)

    DO i = 1, Glob_n
      MassVelocity_kl = MassVelocity_kl-&
           nJiinnJiin_kl(n, i, mk, ml, Skl, Ak, Al, invAkl)/&
           (8.0_wp*Glob_Mass(i+1)**3*cvel*cvel)

    ENDDO
  !
  !    MassVelocity_kl= -ONE/(2.0_wp*Glob_Mass(1)*cvel*cvel)*mass_kl (n,mk,ml,Skl,Ak,Al,invAkl,Vkl)-&
  !         ONE/(2.0_wp*cvel*cvel*Glob_Mass(2))*mass_kl (n,mk,ml,Skl,Ak,Al,invAkl,Vkl)


  !    MassVelocity_kl = &
  !         -nJnnJn_kl(n,mk,ml,Skl,Ak,Al,invAkl)/(8.0_wp*cvel*cvel)*&
  !         ( ONE/Glob_Mass(1)**3 +ONE/Glob_Mass(1+1)**3 )


  CONTAINS


    !!!

    REAL(wp) FUNCTION At_Orb_Orb1(Ak, Al, Skl, invAkl, rkl_m1)
      IMPLICIT NONE
      REAL(wp), DIMENSION(Glob_n, Glob_n), INTENT(IN) :: Ak, Al, invAkl, rkl_m1
      REAL(wp), INTENT(IN)                            :: Skl

      REAL(wp), DIMENSION(Glob_n)         :: O_O1
      REAL(wp), DIMENSION(Glob_n, Glob_n) :: O_O2, O_O3

      CALL At_Orb_Orb_pomoc(Ak, Al, Skl, invAkl, rkl_m1, O_O1, O_O2, O_O3)
      At_Orb_Orb1 = (-SUM(O_O1)-SUM(O_O2)+SUM(O_O3))/(2.0_wp)
    END FUNCTION At_Orb_Orb1


    ! The charge products are taken from Glob_ChargeMatrix (diagonal: q_i*q_0,
    ! strict lower triangle: q_i*q_j) instead of Glob_PseudoCharge directly,
    ! so that the REPULSION/ATTRACTION_SCALING_PARAM records of the data file
    ! act on them as in the PG_0S reference (scaled charge products).
    SUBROUTINE At_Orb_Orb_pomoc(Ak, Al, Skl, invAkl, rkl_m1, O_O1, O_O2, O_O3)
      IMPLICIT NONE
      REAL(wp), DIMENSION(Glob_n, Glob_n), INTENT(IN)  :: Ak, Al, invAkl, rkl_m1
      REAL(wp), INTENT(IN)                             :: Skl
      REAL(wp), DIMENSION(Glob_n), INTENT(OUT)         :: O_O1
      REAL(wp), DIMENSION(Glob_n, Glob_n), INTENT(OUT) :: O_O2, O_O3

      REAL(wp), DIMENSION(Glob_n, Glob_n) :: Bmac, Cmac, Kmac
      INTEGER                                :: i, j, ii, jj


      O_O1 = ZERO
      O_O2 = ZERO
      O_O3 = ZERO
      DO ii = 1, Glob_n
        ! suma I
        ! I-1
        O_O1(ii) = -12.0_wp*Al(ii, ii)*rkl_m1(ii, ii)
        ! I-2
        DO i = 1, Glob_n
          DO j = 1, Glob_n
            Bmac(i, j) = 8*Al(i, ii)*Al(ii, j)  ! MONIKA (i,ii)*(ii,j)
          ENDDO
        ENDDO
        ! O_O1(ii)=O_O1(ii)+8.0_wp*At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,0)
        ! I-3
        ! Bmac=ZERO;
        Bmac(ii, 1:Glob_n) = Bmac(ii, 1:Glob_n)+4*Al(ii, ii)*(Ak(ii, 1:Glob_n)+5.0_wp*Al(ii, 1:Glob_n))
        ! O_O1(ii)=O_O1(ii)+At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,0)
        O_O1(ii) = O_O1(ii)+ME_rXr_over_rij(Bmac, ii, ii, invAkl, rkl_m1(ii, ii))
        ! I-4
        Bmac = ZERO; Bmac(ii, 1:Glob_n) = Al(ii, 1:Glob_n)
        DO i = 1, Glob_n
          DO j = 1, Glob_n
            Cmac(i, j) = (Ak(i, ii)+Al(i, ii))*Al(ii, j)
          ENDDO
        ENDDO
        ! O_O1(ii)=O_O1(ii)-8.0_wp*At_rijm1rBrrCr(Glob_n,Skl,invAkl,Bmac,Cmac,ii,0)
        O_O1(ii) = O_O1(ii)-8.0_wp*ME_rXr_rYr_over_rij(Bmac, Cmac, ii, ii, invAkl, rkl_m1(ii, ii))
        O_O1(ii) = O_O1(ii)*Glob_ChargeMatrix(ii, ii)/Glob_Mass(1)/Glob_Mass(ii+1)

        DO jj = 1, Glob_n
          ! startuje suma II
          IF (jj .NE. ii) THEN
            ! II-1
            O_O2(ii, jj) = -12.0_wp*Al(ii, jj)*rkl_m1(ii, ii)
            ! II-2
            DO i = 1, Glob_n
              DO j = 1, Glob_n
                Bmac(i, j) = 8*Al(i, ii)*Al(jj, j)  ! MONIKA (i,ii)*(ii,j)
              ENDDO
            ENDDO
            ! O_O2(ii,jj)=O_O2(ii,jj)+8.0_wp*At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,0)
            ! II-3
            ! Bmac=ZERO;
            Bmac(ii, 1:Glob_n) = Bmac(ii, 1:Glob_n)+4*Al(ii, jj)*(Ak(ii, 1:Glob_n)+2.0_wp*Al(ii, 1:Glob_n))
            ! O_O2(ii,jj)=O_O2(ii,jj)+At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,0)
            ! II-4
            ! Bmac=ZERO;
            Bmac(ii, 1:Glob_n) = Bmac(ii, 1:Glob_n)+12*Al(ii, ii)*Al(jj, 1:Glob_n)
            ! O_O2(ii,jj)=O_O2(ii,jj)+At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,0)
            O_O2(ii, jj) = O_O2(ii, jj)+ME_rXr_over_rij(Bmac, ii, ii, invAkl, rkl_m1(ii, ii))
            ! II-5
            Bmac = ZERO;
            Bmac(ii, 1:Glob_n) = Al(jj, 1:Glob_n)
            DO i = 1, Glob_n
              DO j = 1, Glob_n
                Cmac(i, j) = (Ak(i, ii)+Al(i, ii))*Al(ii, j)
              ENDDO
            ENDDO
            ! O_O2(ii,jj)=O_O2(ii,jj)-8.0_wp*At_rijm1rBrrCr(Glob_n,Skl,invAkl,Bmac,Cmac,ii,0)
            O_O2(ii, jj) = O_O2(ii, jj)-8.0_wp*ME_rXr_rYr_over_rij(Bmac, Cmac, ii, ii, invAkl, rkl_m1(ii, ii))
            O_O2(ii, jj) = O_O2(ii, jj)*Glob_ChargeMatrix(ii, ii)/Glob_Mass(1)/Glob_Mass(ii+1)
          ENDIF
          ! startuje III suma
          IF (jj .GT. ii) THEN
            ! III-1
            O_O3(ii, jj) = -12.0_wp*Al(jj, ii)*rkl_m1(ii, jj)
            ! III-2
            DO i = 1, Glob_n
              DO j = 1, Glob_n
                Bmac(i, j) = 8*Al(i, ii)*Al(jj, j)
              ENDDO
            ENDDO
            ! O_O3(ii,jj) = O_O3(ii,jj)+8.0_wp*At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,jj)
            ! III-3
            ! Bmac=ZERO;
            Bmac(ii, 1:Glob_n) = Bmac(ii, 1:Glob_n)-4*Al(jj, ii)*(Ak(jj, 1:Glob_n)+4.0_wp*Al(jj, 1:Glob_n))
            Bmac(jj, 1:Glob_n) = Bmac(jj, 1:Glob_n)+4*Al(jj, ii)*(Ak(jj, 1:Glob_n)+4.0_wp*Al(jj, 1:Glob_n))
            ! O_O3(ii,jj) = O_O3(ii,jj)+At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,jj)
            ! III-4
            ! Bmac=ZERO;
            Bmac(ii, 1:Glob_n) = Bmac(ii, 1:Glob_n)-4*Al(jj, jj)*Al(ii, 1:Glob_n)
            Bmac(jj, 1:Glob_n) = Bmac(jj, 1:Glob_n)+4*Al(jj, jj)*Al(ii, 1:Glob_n)
            ! O_O3(ii,jj) = O_O3(ii,jj)+At_rijm1rBr(Glob_n,Skl,invAkl,Bmac,ii,jj)
            O_O3(ii, jj) = O_O3(ii, jj)+ME_rXr_over_rij(Bmac, ii, jj, invAkl, rkl_m1(ii, jj))
            ! III-5
            Bmac = ZERO; Cmac = ZERO
            Bmac(ii, 1:Glob_n) = Al(jj, 1:Glob_n)
            Bmac(jj, 1:Glob_n) = -Al(jj, 1:Glob_n)
            DO i = 1, Glob_n
              DO j = 1, Glob_n
                Cmac(i, j) = (Al(i, jj)+Ak(i, jj))*Al(ii, j)
              ENDDO
            ENDDO
            ! O_O3(ii,jj) = O_O3(ii,jj) + 8.0_wp*At_rijm1rBrrCr(Glob_n,Skl,invAkl,Bmac,Cmac,ii,jj)
            O_O3(ii, jj) = O_O3(ii, jj) + 8.0_wp*ME_rXr_rYr_over_rij(Bmac, Cmac, ii, jj, invAkl, rkl_m1(ii, jj))
            O_O3(ii, jj) = O_O3(ii, jj)*Glob_ChargeMatrix(jj, ii)/Glob_Mass(ii+1)/Glob_Mass(jj+1)
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE At_Orb_Orb_pomoc


    FUNCTION ME_rXr_over_rij(X, i, j, invCkl, ME_1_over_rij)
      ! FUNCTION ME_rXr_over_rij computes the following
      ! matrix element with simple Gaussians phi_k and phi_l:
      !<phi_k| r'Xr/r_{ij} |phi_l>
      ! Here X is an arbitrary (i.e. nonsymmetric) complex matrix.
      ! Index i can be equal to j. In the latter CASE
      !<phi_k| r'Xr/r_{i} |phi_l> is computed
      ! Input:
      !            X  :: Glob_n x Glob_n matrix
      !           i,j :: indices denoting i and j.
      !        invCkl :: Glob_n x Glob_n matrix where the inverse of Ck+tCl is stored
      ! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element

      IMPLICIT NONE

      REAL(wp) :: ME_rXr_over_rij

      ! Arguments:
      REAL(wp) :: X(Glob_n, Glob_n), invCkl(Glob_n, Glob_n), ME_1_over_rij
      INTEGER     :: i, j
      ! Local variables
      REAL(wp) :: TrCJ
      REAL(wp) :: CXi(Glob_n), CXj(Glob_n)
      REAL(wp) :: TrCX, TrCXCJ
      INTEGER     :: m, p, q

      IF (i == j) THEN
        TrCJ = invCkl(i, i)
      ELSE
        TrCJ = invCkl(i, i)+invCkl(j, j)-invCkl(j, i)-invCkl(j, i)
      ENDIF

      TrCX = ZERO
      DO m = 1, Glob_n
        DO p = 1, Glob_n
          TrCX = TrCX+invCkl(m, p)*X(p, m)
        ENDDO
      ENDDO

      DO m = 1, Glob_n
        CXi(m) = ZERO
        DO p = 1, Glob_n
          CXi(m) = CXi(m)+invCkl(i, p)*X(p, m)
        ENDDO
      ENDDO
      IF (j /= i) THEN
        DO m = 1, Glob_n
          CXj(m) = ZERO
          DO p = 1, Glob_n
            CXj(m) = CXj(m)+invCkl(j, p)*X(p, m)
          ENDDO
        ENDDO
      ENDIF

      IF (i == j) THEN
        TrCXCJ = ZERO
        DO m = 1, Glob_n
          TrCXCJ = TrCXCJ+CXi(m)*invCkl(m, i)
        ENDDO
      ELSE
        TrCXCJ = ZERO
        DO m = 1, Glob_n
          TrCXCJ = TrCXCJ+(CXi(m)-CXj(m))*(invCkl(m, i)-invCkl(m, j))
        ENDDO
      ENDIF

      ME_rXr_over_rij = ME_1_over_rij*(3*TrCX-TrCXCJ/TrCJ)/2

    END FUNCTION ME_rXr_over_rij


    FUNCTION ME_rXr_rYr_over_rij(X, Y, i, j, invCkl, ME_1_over_rij)
      ! FUNCTION ME_rXr_rYr_over_rij computes the following
      ! matrix element with simple Gaussians phi_k and phi_l:
      !<phi_k| (r'Xr)(r'Yr)/r_{ij} |phi_l>
      ! Here X and Y are arbitrary (i.e. nonsymmetric) matrices.
      ! Index i can be equal to j. In the latter CASE
      !<phi_k| (r'Xr)(r'Yr)/r_{i} |phi_l> is computed
      ! Input:
      !            X  :: Glob_n x Glob_n matrix
      !            Y  :: Glob_n x Glob_n matrix
      !           i,j :: indices denoting i and j.
      !        invCkl :: Glob_n x Glob_n matrix where the inverse of Ck+Cl is stored
      ! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
      IMPLICIT NONE

      REAL(wp) :: ME_rXr_rYr_over_rij
      ! Arguments:
      REAL(wp) :: X(Glob_n, Glob_n), Y(Glob_n, Glob_n), invCkl(Glob_n, Glob_n), ME_1_over_rij
      INTEGER     :: i, j
      ! Local variables
      REAL(wp) :: TrCJ
      REAL(wp) :: Ys(Glob_n, Glob_n), Xs(Glob_n, Glob_n), CX(Glob_n, Glob_n), CY(Glob_n, Glob_n), CXCYi(Glob_n), CXCYj(Glob_n)
      REAL(wp) :: TrCX, TrCY, TrCXCY, TrCXCJ, TrCYCJ, TrCXCYCJ
      INTEGER     :: m, p, q

      IF (i == j) THEN
        TrCJ = invCkl(i, i)
      ELSE
        TrCJ = invCkl(i, i)+invCkl(j, j)-invCkl(j, i)-invCkl(j, i)
      ENDIF

      DO p = 1, Glob_n
        Ys(p, p) = Y(p, p)
        DO q = p+1, Glob_n
          Ys(p, q) = (Y(p, q)+Y(q, p))/2
          Ys(q, p) = Ys(p, q)
        ENDDO
      ENDDO
      DO p = 1, Glob_n
        Xs(p, p) = X(p, p)
        DO q = p+1, Glob_n
          Xs(p, q) = (X(p, q)+X(q, p))/2
          Xs(q, p) = Xs(p, q)
        ENDDO
      ENDDO

      DO p = 1, Glob_n
        DO q = 1, Glob_n
          CY(p, q) = ZERO
          DO m = 1, Glob_n
            CY(p, q) = CY(p, q)+invCkl(p, m)*Ys(m, q)
          ENDDO
        ENDDO
      ENDDO
      DO p = 1, Glob_n
        DO q = 1, Glob_n
          CX(p, q) = ZERO
          DO m = 1, Glob_n
            CX(p, q) = CX(p, q)+invCkl(p, m)*Xs(m, q)
          ENDDO
        ENDDO
      ENDDO

      DO m = 1, Glob_n
        CXCYi(m) = ZERO
        DO p = 1, Glob_n
          CXCYi(m) = CXCYi(m)+CX(i, p)*CY(p, m)
        ENDDO
      ENDDO
      IF (j /= i) THEN
        DO m = 1, Glob_n
          CXCYj(m) = ZERO
          DO p = 1, Glob_n
            CXCYj(m) = CXCYj(m)+CX(j, p)*CY(p, m)
          ENDDO
        ENDDO
      ENDIF

      TrCY = ZERO
      TrCX = ZERO
      DO m = 1, Glob_n
        TrCY = TrCY+CY(m, m)
        TrCX = TrCX+CX(m, m)
      ENDDO
      TrCXCY = ZERO
      DO m = 1, Glob_n
        DO p = 1, Glob_n
          TrCXCY = TrCXCY+CX(m, p)*CY(p, m)
        ENDDO
      ENDDO

      IF (i == j) THEN
        TrCXCJ = ZERO
        DO m = 1, Glob_n
          TrCXCJ = TrCXCJ+CX(i, m)*invCkl(m, i)
        ENDDO
        TrCYCJ = ZERO
        DO m = 1, Glob_n
          TrCYCJ = TrCYCJ+CY(i, m)*invCkl(m, i)
        ENDDO
        TrCXCYCJ = ZERO
        DO m = 1, Glob_n
          TrCXCYCJ = TrCXCYCJ+CXCYi(m)*invCkl(m, i)
        ENDDO
      ELSE
        TrCXCJ = ZERO
        DO m = 1, Glob_n
          TrCXCJ = TrCXCJ+(CX(i, m)-CX(j, m))*(invCkl(m, i)-invCkl(m, j))
        ENDDO
        TrCYCJ = ZERO
        DO m = 1, Glob_n
          TrCYCJ = TrCYCJ+(CY(i, m)-CY(j, m))*(invCkl(m, i)-invCkl(m, j))
        ENDDO
        TrCXCYCJ = ZERO
        DO m = 1, Glob_n
          TrCXCYCJ = TrCXCYCJ+(CXCYi(m)-CXCYj(m))*(invCkl(m, i)-invCkl(m, j))
        ENDDO
      ENDIF

      ME_rXr_rYr_over_rij = ME_1_over_rij*(&
           9*TrCY*TrCX + 6*TrCXCY &
           - (3*(TrCY*TrCXCJ + TrCX*TrCYCJ) + 4*TrCXCYCJ)/TrCJ &
           + 3*TrCXCJ*TrCYCJ/(TrCJ*TrCJ) &
          )/4

    END FUNCTION ME_rXr_rYr_over_rij


    !
    ! Begin DK
    !
    ! Procedures needed to claculate relativistic correction
    !


    SUBROUTINE r1m2(n, mk, ml, Skl, invAkl, r1m2_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-2}|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl
      REAL(wp), INTENT(OUT)                     :: r1m2_kl

      REAL(wp) :: pak1, pak2
      pak1 = 2.0_wp/REAL(mk+ml+1)
      pak2 = Skl/invAkl(1, 1)
      r1m2_kl = pak1*pak2
    END SUBROUTINE r1m2


    SUBROUTINE r1m4(n, mk, ml, Skl, invAkl, r1m4_kl)
      ! Calculate matrix element:
      !    <\phi_k|r_1^{-4}|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl
      REAL(wp), INTENT(OUT)                     :: r1m4_kl

      r1m4_kl = 4.0_wp/((mk+ml)**2-1)/invAkl(1, 1)**2*Skl
    END SUBROUTINE r1m4


    SUBROUTINE rBr(n, mk, ml, Skl, B, invAkl, rBr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r'Br|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B
      REAL(wp), INTENT(OUT)                     :: rBr_kl

      INTEGER     :: i, j
      REAL(wp) :: ABA11, trAB, pak

      pak = REAL((mk+ml), wp)/invAkl(1, 1)

      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          trAB = trAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO
      rBr_kl = 0.5_wp*Skl*(3.0_wp*trAB+ pak*ABA11)
    END SUBROUTINE rBr


    SUBROUTINE r1m2_rBr(n, mk, ml, Skl, B, invAkl, r1m2_rBr_kl)
      ! Calculate matrix element:
      !  <\phi_k|r_1^{-2} r'Br|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B
      REAL(wp), INTENT(OUT)                     :: r1m2_rBr_kl

      INTEGER     :: i, j
      REAL(wp) :: ABA11, trAB

      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          trAB = trAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO
      r1m2_rBr_kl = ONE/(mk+ml+ONE)*&
           Skl/invAkl(1, 1)*(3.0_wp*trAB+(mk+ml-2)*ABA11/invAkl(1, 1))

    END SUBROUTINE r1m2_rBr

    SUBROUTINE nower1m2_rBr(n, mk, ml, Skl, B, invAkl, nower1m2_rBr_kl)
      ! Calculate matrix element:
      ! <\phi_k|r_1^{-2} r'Br|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B
      REAL(wp), INTENT(OUT)                     :: nower1m2_rBr_kl

      INTEGER     :: i, j
      REAL(wp) :: ABA11, trAB
      REAL(wp) :: pak1, pak2, pak3

      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          TrAB = TrAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO
      pak1 = Skl/invAkl(1, 1)
      pak2 = REAL((mk+ml-2), wp)/REAL((mk+ml+1), wp)
      pak3 = ABA11/invAkl(1, 1)

      nower1m2_rBr_kl = ONE/(mk+ml+ONE)*&
           pak1*3.0_wp*TrAB + pak2*pak1*pak3

    END SUBROUTINE nower1m2_rBr


    SUBROUTINE r1m4_rBr(n, mk, ml, Skl, B, invAkl, r1m4_rBr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-4} r'Br|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B
      REAL(wp), INTENT(OUT)                     :: r1m4_rBr_kl


      INTEGER     :: i, j
      REAL(wp) :: ABA11, trAB


      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          trAB = trAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO
      r1m4_rBr_kl = 2.0_wp/((mk+ml)**2-ONE)*Skl/invAkl(1, 1)**2&
           *(3.0_wp*trAB+(mk+ml-4)*ABA11/invAkl(1, 1))

    END SUBROUTINE r1m4_rBr


    SUBROUTINE rBrrCr(n, mk, ml, Skl, B, C, invAkl, rBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|(r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: rBrrCr_kl


      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11

      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      rBrrCr_kl = 0.25_wp*Skl* &
           (9.0_wp*trAB*trAC+6.0_wp*trABAC&
           +3.0_wp*(mk+ml)/invAkl(1, 1)*(trAB*ACA11+trAC*ABA11)&
           +(mk+ml-2)*(mk+ml)/invAkl(1, 1)**2*ABA11*ACA11&
           +4.0_wp*(mk+ml)/invAkl(1, 1)*ABACA11)
      !   write(*,*)'rBrrCr_kl=',rBrrCr_kl
    END SUBROUTINE rBrrCr

    SUBROUTINE nowerBrrCr(n, mk, ml, Skl, B, C, invAkl, nowerBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|(r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: nowerBrrCr_kl


      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11
      REAL(wp) :: pak1, pak2, pak3

      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      pak1 = ACA11/invAkl(1, 1)
      pak2 = ABA11/invAkl(1, 1)
      pak3 = ABACA11/invAkl(1, 1)

      nowerBrrCr_kl = 0.25_wp*Skl* &
           (9.0_wp*trAB*trAC+6.0_wp*trABAC&
           +3.0_wp*(mk+ml)*(trAB*pak1+trAC*pak2)&
           +(mk+ml-2)*(mk+ml)*pak1*pak2&
           +4.0_wp*(mk+ml)*pak3)

    END SUBROUTINE nowerBrrCr


    SUBROUTINE r1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m2_rBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-2} (r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: r1m2_rBrrCr_kl

      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11

      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      r1m2_rBrrCr_kl = 0.5_wp/(mk+ml+1)*Skl/invAkl(1, 1)* &
           (9.0_wp*trAB*trAC+6.0_wp*trABAC&
           +3.0_wp*(mk+ml-2)/invAkl(1, 1)*(trAB*ACA11+trAC*ABA11)&
           +(mk+ml-4)*(mk+ml-2)/invAkl(1, 1)**2*ABA11*ACA11&
           +4.0_wp*(mk+ml-2)/invAkl(1, 1)*ABACA11)

    END SUBROUTINE r1m2_rBrrCr

    SUBROUTINE nower1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m2_rBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-2} (r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: nower1m2_rBrrCr_kl

      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11

      REAL(wp) :: pak1, pak2, pak3, pak4, pak5, pak6, pak7, pak8, pak9, pak10, pak11
      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      pak1 = ACA11/invAkl(1, 1)
      pak2 = ABA11/invAkl(1, 1)
      pak3 = ABACA11/invAkl(1, 1)
      pak4 = Skl/invAkl(1, 1)
      pak5 = REAL((mk+ml-2), wp)/REAL((mk+ml+1), wp)
      pak6 = REAL((mk+ml-4)*(mk+ml-2), wp)/REAL((mk+ml+1), wp)
      pak7 = REAL((mk+ml-2), wp)/REAL((mk+ml+1), wp)
      pak8 = 0.5_wp/REAL((mk+ml+1), wp)
      pak9 = pak6*pak1*pak2
      pak10 = 4.0_wp*pak7*pak3
      pak11 = 3.0_wp*pak5

      nower1m2_rBrrCr_kl = pak8*pak4*(9.0_wp*TrAB*TrAC + 6.0_wp*TrABAC)+&
                0.5_wp*pak4*(pak11*(TrAB*pak1 + TrAC*pak2) + pak9 + pak10)

    END SUBROUTINE nower1m2_rBrrCr


    SUBROUTINE r1m4_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m4_rBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-4} (r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: r1m4_rBrrCr_kl

      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11

      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      r1m4_rBrrCr_kl = ONE/((mk+ml)**2-1)*Skl/invAkl(1, 1)**2* &
           (9.0_wp*trAB*trAC+6.0_wp*trABAC&
           +3.0_wp*(mk+ml-4)/invAkl(1, 1)*(trAB*ACA11+trAC*ABA11)&
           +(mk+ml-6)*(mk+ml-4)/invAkl(1, 1)**2*ABA11*ACA11&
           +4.0_wp*(mk+ml-4)/invAkl(1, 1)*ABACA11)

    END SUBROUTINE r1m4_rBrrCr

    SUBROUTINE nower1m4_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m4_rBrrCr_kl)
      ! Calculate matrix element:
      !     <\phi_k|r_1^{-4} (r'Br)(r'Cr)|\phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, B, C
      REAL(wp), INTENT(OUT)                     :: nower1m4_rBrrCr_kl


      INTEGER     :: i, j, k, l
      REAL(wp) :: trAB, trAc, trABAC
      REAL(wp) :: ABA11, ACA11, ABACA11
      REAL(wp) :: pak1, pak2, pak3, pak4, pak5, pak6, pak7, pak8, pak9, pak10

      trAB = ZERO
      trAC = ZERO
      ABA11 = ZERO
      ACA11 = ZERO
      trABAC = ZERO
      ABACA11 = ZERO

      DO i = 1, n
        DO j = 1, n
          trAB = trAB+invAkl(i, j)*B(j, i)
          trAC = trAC+invAkl(i, j)*C(j, i)
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          DO k = 1, n
            DO l = 1, n
              trABAC = trABAC+invAkl(i, j)*B(j, k)*invAkl(k, l)*C(l, i)
              ABACA11 = ABACA11+invAkl(1, i)*B(i, j)*invAkl(j, k)*C(k, l)&
                             *invAkl(l, 1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      pak1 = Skl/invAkl(1, 1)**2
      pak2 = ACA11/invAkl(1, 1)
      pak3 = ABA11/invAkl(1, 1)
      pak4 = ABACA11/invAkl(1, 1)
      pak5 = REAL((mk+ml-6), wp)/REAL((mk+ml+1), wp)
      pak6 = REAL((mk+ml-4), wp)/REAL((mk+ml-1), wp)
      pak7 = REAL(mk, wp)/REAL((mk+ml-1), wp)
      pak8 = REAL(ml, wp)/REAL((mk+ml+1), wp)
      pak9 = pak7*pak8*pak1
      pak10 = pak8*pak1*pak6

      nower1m4_rBrrCr_kl = pak9* 9.0_wp*TrAB*TrAC +&
           pak9*6.0_wp*TrABAC+&
           mk*3.0_wp*pak10*(TrAB*pak2 + TrAC*pak3)+&
           pak1*REAL(ml*mk, wp)*pak5*pak6*pak2*pak3+&
           pak10*4.0_wp*REAL(mk, wp)*pak4

    END SUBROUTINE nower1m4_rBrrCr


    REAL(wp) FUNCTION nJnnJn_kl (n, mk, ml, Skl, Ak, Al, invAkl)
      ! Calculate matrix element:
      !     <\nabla' J \nabla \phi_k|\nabla' J \nabla \phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, Ak, Al


      INTEGER                          :: i, j, k, l
      REAL(wp)                      :: trJAk, trJAl, ABA11, TrAB, sum111
      REAL(wp)                      :: ACA11, TrAC
      REAL(wp), DIMENSION(1:n, 1:n) :: B, C
      REAL(wp)                      :: r1m2_kl, r1m4_kl
      REAL(wp)                      :: rBr_kl, r1m4_rBr_kl, r1m2_rBr_kl
      REAL(wp)                      :: rBrrCr_kl, r1m2_rBrrCr_kl, r1m4_rBrrCr_kl
      REAL(wp), DIMENSION(1:16)     :: Wyn, Wod, wynw, Stare
      REAL(wp)                      :: pak1, pak2, pak1a, pak1b, pak2a, pak2b, pak
      REAL(wp)                      :: pak3, pak4, pak5, pak6, pak7
      REAL(wp)                      :: nowerBrrCr_kl, nower1m2_rBr_kl, nower1m2_rBrrCr_kl, nower1m4_rBrrCr_kl

      trJAk = ZERO
      trJAl = ZERO
      DO i = 1, n
        DO j = 1, n
          trJAk = trJAk+Ak(i, j)
          trJAl = trJAl+Al(i, j)
        ENDDO
      ENDDO

      Wyn = ZERO
      Wod = ZERO
      Stare = ZERO
      Wyn(1) = 36.0_wp*trJAk*trJAl*Skl
      Wod(1) = 36.0_wp*Ak(1, 1)*Al(1, 1)*Skl

      !    write(*,*)'1  ', wyn(1)

      DO i = 1, n
        DO j = 1, n
          B(i, j) = ZERO
          C(i, j) = ZERO
          DO k = 1, n
            DO l = 1, n
              B(i, j) = B(i, j)+Al(i, k)*Al(l, j)
              C(i, j) = C(i, j)+Ak(i, k)*Ak(l, j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      CALL rBr(n, mk, ml, Skl, B, invAkl, rBr_kl)
      Wyn(3) = -24.0_wp*trJAk*rBr_kl

      !    write(*,*)'3',wyn(3)

      CALL rBr(n, mk, ml, Skl, C, invAkl, rBr_kl)
      Wyn(2) = -24.0_wp*trJAl*rBr_kl

      Wod(2) = -12.0_wp*Ak(1, 1)*Ak(1, 1)*Al(1, 1)*Skl*&
           REAL((mk+ml+3), wp)*invAkl(1, 1)

      !     write(*,*)'2',wyn(2)


      ! begin Wyn(4)********************************************
      CALL rBrrCr(n, mk, ml, Skl, C, B, invAkl, rBrrCr_kl)
      Stare(4) = +16.0_wp*rBrrCr_kl

      CALL nowerBrrCr(n, mk, ml, Skl, C, B, invAkl, nowerBrrCr_kl)
      Wyn(4) = 16.0_wp*nowerBrrCr_kl


      !    write(*,*) '4',wyn(4),Stare(4)
      ! end Wyn(4)*******************************************************
      ! begin Wyn(5) i Wyn(6)********************************************
      CALL r1m2(n, mk, ml, Skl, invAkl, r1m2_kl)
      Stare(5) = -6.0_wp*(mk+1)*mk*trJAl*r1m2_kl
      Stare(6) = -6.0_wp*(ml+1)*ml*trJAk*r1m2_kl

      pak1 = REAL((mk+1)*mk, wp)/REAL((mk+ml+1), wp)
      pak2 = Skl/invAkl(1, 1)

      Wyn(5) = -12.0_wp*TrJAl*pak1*pak2

      pak1 = REAL((ml+1)*ml, wp)/REAL((mk+ml+1), wp)
      pak2 = Skl/invAkl(1, 1)

      Wyn(6) = -12.0_wp*TrJAk*pak1*pak2

      !    write(*,*)'5',wyn(5),stare(5)
      !    write(*,*) '6',wyn(6),stare(6)
      ! end Wyn(5) i Wyn(6)*******************************************************
      ! DO Wyn(7)********************************************************
      CALL r1m2_rBr(n, mk, ml, Skl, B, invAkl, r1m2_rBr_kl)
      Stare(7) = +4.0_wp*(mk+1)*mk*r1m2_rBr_kl

      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          trAB = trAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO

      pak1a = REAL(mk*(mk+1), wp)/(REAL((mk+ml+1), wp))
      pak1b = Skl/invAkl(1, 1)
      pak2a = REAL((mk+ml-2), wp)/(REAL((mk+ml+1), wp))
      pak2b = ABA11/invAkl(1, 1)

      Wyn(7) = 12.0_wp*pak1a*pak1b*TrAB +&
           4.0_wp*REAL(mk*(mk+1), wp)*pak2a*pak1b*pak2b


      Wod(7) = 4.0_wp*Skl*Al(1, 1)*Al(1, 1)*REAL(mk*(mk+1), wp)
      !    write(*,*) '7',wyn(7),stare(7),wod(7)
      ! end Wyn(7)********************************************************

      ! begin Wyn(8)*******************************************************
      CALL r1m2_rBr(n, mk, ml, Skl, C, invAkl, r1m2_rBr_kl)
      Stare(8) = +4.0_wp*(ml+1)*ml*r1m2_rBr_kl

      ACA11 = ZERO
      trAC = ZERO
      DO i = 1, n
        DO j = 1, n
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          trAC = trAC+invAkl(i, j)*C(j, i)
        ENDDO
      ENDDO

      pak1a = REAL(ml*(ml+1), wp)/(REAL((mk+ml+1), wp))
      pak1b = Skl/invAkl(1, 1)
      pak2a = REAL((mk+ml-2), wp)/(REAL((mk+ml+1), wp))
      pak2b = ACA11/invAkl(1, 1)

      Wyn(8) = 12.0_wp*pak1a*pak1b*TrAC+&
           4.0_wp*REAL(ml*(ml+1), wp)*pak2a*pak1b*pak2b

      Wod(8) = 4.0_wp*Skl*Ak(1, 1)*Ak(1, 1)*REAL(ml*(ml+1), wp)

      !    write(*,*) '8',wyn(8),stare(8),wod(8)
      ! end Wyn(8)*********************************************************
      ! begin Wyn(9) i Wyn(10)*********************************************
      DO i = 1, n
        DO j = 1, n
          B(i, j) = ZERO
          C(i, j) = ZERO
          B(i, 1) = B(i, 1) + Ak(i, j)
          B(1, i) = B(1, i) + Ak(j, i)
          C(i, 1) = C(i, 1) + Al(i, j)
          C(1, i) = C(1, i) + Al(j, i)
        ENDDO
      ENDDO

      CALL r1m2_rBr(n, mk, ml, Skl, B, invAkl, r1m2_rBr_kl)
      Stare(9) = +12.0_wp*mk*trJAl*r1m2_rBr_kl

      CALL nower1m2_rBr(n, mk, ml, Skl, B, invAkl, nower1m2_rBr_kl)
      Wyn(9) = 12.0_wp*REAL(mk, wp)*TrJAl*nower1m2_rBr_kl
      Wod(9) = 24.0_wp*REAL(mk, wp)*Ak(1, 1)*Al(1, 1)*Skl

      !    write(*,*) '9',wyn(9),Stare(9),wod(9)

      CALL r1m2_rBr(n, mk, ml, Skl, C, invAkl, r1m2_rBr_kl)
      Stare(10) = +12.0_wp*ml*trJAk*r1m2_rBr_kl

      CALL nower1m2_rBr(n, mk, ml, Skl, C, invAkl, nower1m2_rBr_kl)
      Wyn(10) = 12.0_wp*ml*TrJAk*nower1m2_rBr_kl
      Wod(10) = 24.0_wp*REAL(ml, wp)*Al(1, 1)*Ak(1, 1)*Skl

      !    write(*,*) '10',wyn(10),Stare(10),wod(10)
      ! end Wyn(9) i Wyn(10)***********************************************
      ! begin Wyn(11)******************************************************
      DO i = 1, n
        DO j = 1, n
          B(i, j) = ZERO
          C(i, j) = ZERO
          B(i, 1) = B(i, 1)+Ak(i, j)
          B(1, i) = B(1, i)+Ak(j, i)
          DO k = 1, n
            DO l = 1, n
              C(i, j) = C(i, j)+Al(i, k)*Al(l, j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL r1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m2_rBrrCr_kl)
      Stare(11) = -8.0_wp*mk*r1m2_rBrrCr_kl

      CALL nower1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m2_rBrrCr_kl)
      Wyn(11) = -8.0_wp*REAL(mk, wp)*nower1m2_rBrrCr_kl
      Wod(11) = -8.0_wp*Al(1, 1)*Al(1, 1)*Ak(1, 1)*&
           Skl*REAL(mk*(mk+ml+3), wp)*invAkl(1, 1)

      !    write(*,*) '11',wyn(11),Stare(11),wod(11)

      ! end Wyn(11) *****************************************
      ! begin Wyn(12) **********************************************
      DO i = 1, n
        DO j = 1, n
          B(i, j) = ZERO
          C(i, j) = ZERO
          B(i, 1) = B(i, 1)+Al(i, j)
          B(1, i) = B(1, i)+Al(j, i)
          DO k = 1, n
            DO l = 1, n
              C(i, j) = C(i, j)+Ak(i, k)*Ak(l, j)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL r1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m2_rBrrCr_kl)
      Stare(12) = -8.0_wp*ml*r1m2_rBrrCr_kl

      CALL nower1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m2_rBrrCr_kl)
      Wyn(12) = -8.0_wp*ml*nower1m2_rBrrCr_kl

      Wod(12) = -8.0_wp*Al(1, 1)*Ak(1, 1)*Ak(1, 1)*Skl*REAL(ml*(mk+ml+3), wp)*invAkl(1, 1)

      !    write(*,*) '12',wyn(12),Stare(12)

      ! begin Wyn(13) **********************************************

      CALL r1m4(n, mk, ml, Skl, invAkl, r1m4_kl)
      Stare(13) = mk*ml*(mk+1)*(ml+1)*r1m4_kl

      pak1 = REAL(mk*(ml+1), wp)/REAL((mk+ml-1), wp)
      pak1a = pak1/invAkl(1, 1)
      pak2 = REAL(ml*(mk+1), wp)/REAL((mk+ml+1), wp)
      pak2a = pak2/invAkl(1, 1)
      pak = pak1a*pak2a

      Wyn(13) = 4.0_wp*Skl*pak
      Wod(13) = Wyn(13)
      !    write(*,*) '13',wyn(13),Stare(13),Wod(13)

      ! begin Wyn(14) Wyn(15)******************************************************************
      DO i = 1, n
        DO j = 1, n
          B(i, j) = ZERO
          C(i, j) = ZERO
          B(i, 1) = B(i, 1)+Ak(i, j)
          B(1, i) = B(1, i)+Ak(j, i)
          C(i, 1) = C(i, 1)+Al(i, j)
          C(1, i) = C(1, i)+Al(j, i)
        ENDDO
      ENDDO

      CALL r1m4_rBr(n, mk, ml, Skl, C, invAkl, r1m4_rBr_kl)
      Stare(14) = -2.0_wp*(mk+1)*mk*ml*r1m4_rBr_kl

      CALL r1m4_rBr(n, mk, ml, Skl, B, invAkl, r1m4_rBr_kl)
      Stare(15) = -2.0_wp*mk*(ml+1)*ml*r1m4_rBr_kl


      ACA11 = ZERO
      trAC = ZERO
      DO i = 1, n
        DO j = 1, n
          ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
          trAC = trAC+invAkl(i, j)*C(j, i)
        ENDDO
      ENDDO

      pak1 = ACA11/invAkl(1, 1)
      pak2 = Skl/invAkl(1, 1)**2.0_wp
      pak3 = REAL((mk+ml-4), wp)/REAL((mk+ml-1), wp)
      pak4 = REAL((mk+1), wp)/REAL((mk+ml+1), wp)
      pak5 = REAL(mk, wp)/REAL((mk+ml-1), wp)
      pak6 = REAL(mk, wp)*REAL(ml, wp)*pak4
      pak7 = REAL(ml, wp)*pak4

      Wyn(14) = -4.0_wp*pak5*pak2*3.0_wp*TrAC*pak7 - &
           4.0_wp*pak2*pak3*pak1*pak6


      ABA11 = ZERO
      trAB = ZERO
      DO i = 1, n
        DO j = 1, n
          ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
          trAB = trAB+invAkl(i, j)*B(j, i)
        ENDDO
      ENDDO

      pak1 = ABA11/invAkl(1, 1)
      pak2 = Skl/invAkl(1, 1)**2
      pak3 = REAL((mk+ml-4), wp)/REAL((mk+ml-1), wp)
      pak4 = REAL((ml+1), wp)/REAL((mk+ml+1), wp)
      pak5 = REAL(mk, wp)/REAL((mk+ml-1), wp)
      pak6 = REAL(mk, wp)*REAL(ml, wp)*pak4
      pak7 = REAL(ml, wp)*pak4

      Wyn(15) = -4.0_wp*pak5*pak2*3.0_wp*TrAB*pak7 - &
           4.0_wp*pak2*pak3*pak1*pak6

      !    write(*,*) '14',wyn(14),Stare(14),wod(14)
      !    write(*,*) '15',wyn(15),stare(15),wod(15)

      ! Wyn(16)***************************************************************
      CALL r1m4_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m4_rBrrCr_kl)
      Stare(16) = +4.0_wp*mk*ml*r1m4_rBrrCr_kl

      CALL nower1m4_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m4_rBrrCr_kl)
      Wyn(16) = +4.0_wp*nower1m4_rBrrCr_kl

      Wod(16) = 16.0_wp*REAL(mk*ml, wp)*Ak(1, 1)*Al(1, 1)*Skl

      !    write(*,*) '16',wyn(16),stare(16),wod(16)
      ! end Wyn(16)***********************************************************

      sum111 = ZERO

      DO i = 1, 16
        sum111 = sum111 + Wyn(i)
      ENDDO

      ! r^2
      !    nJnnJn_kl = Wod(2)+Wod(3)+Wod(11)+Wod(12)
      ! r^-2
      !
      !    nJnnJn_kl = Wod(5)+Wod(6)+Wod(14)+Wod(15)
      ! r^-4
      !       nJnnJn_kl = Wod(13)
      ! r^4
      !    nJnnJn_kl = Wod(4)
      ! s
      !    nJnnJn_kl = Wod(1) + Wod(7) + Wod(8) + Wod(9) + Wod(10) + Wod(16)

      nJnnJn_kl = sum111
    END FUNCTION nJnnJn_kl


    REAL(wp) FUNCTION nJiinnJiin_kl (n, ii, mk, ml, Skl, Ak, Al, invAkl)
      ! Calculate matrix element:
      !     <\nabla' Jii \nabla \phi_k|\nabla' Jii \nabla \phi_l>
      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, ii, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl, Ak, Al


      INTEGER                          :: i, j, k, l
      REAL(wp), DIMENSION(1:n, 1:n) :: B, C
      REAL(wp)                      :: r1m2_kl, r1m4_kl, sum111
      REAL(wp)                      :: rBr_kl, r1m2_rBr_kl, r1m4_rBr_kl
      REAL(wp)                      :: rBrrCr_kl, r1m2_rBrrCr_kl, r1m4_rBrrCr_kl
      REAL(wp), DIMENSION(1:16)     :: Wyn, Stare

      REAL(wp) :: pak1, pak2, pak3, pak4, pak5, pak6, pak7, pak8, pak9
      REAL(wp) :: nowerBrrCr_kl, nower1m2_rBr_kl, nower1m2_rBrrCr_kl, nower1m4_rBrrCr_kl
      REAL(wp) :: ACA11, TrAB, ABA11, TrAC, TrJAk, TrJAl

      Wyn = ZERO
      Stare = ZERO

      Wyn(1) = 36.0_wp*Ak(ii, ii)*Al(ii, ii)*Skl
      !    write(*,*) '1 ii',wyn(1)
      DO i = 1, n
        DO j = 1, n
          B(i, j) = Ak(i, ii)*Ak(ii, j)
          C(i, j) = Al(i, ii)*Al(ii, j)
        ENDDO
      ENDDO
      CALL rBr(n, mk, ml, Skl, B, invAkl, rBr_kl)
      Wyn(2) = -24.0_wp*Al(ii, ii)*rBr_kl
      !    write(*,*) '2 ii',wyn(2)
      CALL rBr(n, mk, ml, Skl, C, invAkl, rBr_kl)
      Wyn(3) = -24.0_wp*Ak(ii, ii)*rBr_kl
      !    write(*,*) '3 ii',wyn(3)

      ! begin Wyn(4)*********************************************
      CALL rBrrCr(n, mk, ml, Skl, B, C, invAkl, rBrrCr_kl)
      Stare(4) = +16.0_wp*rBrrCr_kl

      CALL nowerBrrCr(n, mk, ml, Skl, C, B, invAkl, nowerBrrCr_kl)
      Wyn(4) = 16.0_wp*nowerBrrCr_kl

      !    write(*,*) '4 ii',wyn(4)
      !*********************************************************
      !
      ! begin Wyn(5) i Wyn(6)************************************
      IF (ii .EQ. 1) THEN
        CALL r1m2(n, mk, ml, Skl, invAkl, r1m2_kl)
        Stare(5) = -6.0_wp*Al(1, 1)*mk*(mk+1)*r1m2_kl
        Stare(6) = -6.0_wp*Ak(1, 1)*ml*(ml+1)*r1m2_kl

        pak1 = REAL((mk+1)*mk, wp)/REAL((mk+ml+1), wp)
        pak2 = Skl/invAkl(1, 1)

        Wyn(5) = -12.0_wp*Al(1, 1)*pak1*pak2


        pak1 = REAL((ml+1)*ml, wp)/REAL((mk+ml+1), wp)
        pak2 = Skl/invAkl(1, 1)

        Wyn(6) = -12.0_wp*Ak(1, 1)*pak1*pak2


        !    write(*,*) '5 ii',wyn(5),stare(5)
        !    write(*,*) '6 ii',wyn(6),stare(6)

        ! begin Wyn(9) i Wyn(10)**************************************
        DO i = 1, n
          DO j = 1, n
            B(i, j) = ZERO
            C(i, j) = ZERO
          ENDDO
          B(1, i) = Ak(1, i)
          B(i, 1) = B(i, 1)+Ak(i, 1)
          C(1, i) = Al(1, i)
          C(i, 1) = C(i, 1)+Al(i, 1)
        ENDDO
        CALL r1m2_rBr(n, mk, ml, Skl, B, invAkl, r1m2_rBr_kl)
        Stare(9) = +12.0_wp*Al(1, 1)*mk*r1m2_rBr_kl

        CALL r1m2_rBr(n, mk, ml, Skl, C, invAkl, r1m2_rBr_kl)
        Stare(10) = +12.0_wp*Ak(1, 1)*ml*r1m2_rBr_kl

        CALL nower1m2_rBr(n, mk, ml, Skl, B, invAkl, nower1m2_rBr_kl)
        Wyn(9) = 12.0_wp*REAL(mk, wp)*Al(1, 1)*nower1m2_rBr_kl

        CALL nower1m2_rBr(n, mk, ml, Skl, C, invAkl, nower1m2_rBr_kl)
        Wyn(10) = 12.0_wp*ml*Ak(1, 1)*nower1m2_rBr_kl

        !       write(*,*) '9 ii', Wyn(9), stare(9)
        !       write(*,*) '10 ii',wyn(10),stare(10)

        ! begin Wyn(7) i Wyn(8)**************************************

        DO i = 1, n
          DO j = 1, n
            B(i, j) = Ak(i, ii)*Ak(ii, j)
            C(i, j) = Al(i, ii)*Al(ii, j)
          ENDDO
        ENDDO

        CALL r1m2_rBr(n, mk, ml, Skl, B, invAkl, r1m2_rBr_kl)
        Stare(8) = +4.0_wp*ml*(ml+1)*r1m2_rBr_kl
        CALL r1m2_rBr(n, mk, ml, Skl, C, invAkl, r1m2_rBr_kl)
        Stare(7) = +4.0_wp*mk*(mk+1)*r1m2_rBr_kl

        ABA11 = ZERO
        trAB = ZERO
        DO i = 1, n
          DO j = 1, n
            ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
            trAB = trAB+invAkl(i, j)*B(j, i)
          ENDDO
        ENDDO
        pak1 = REAL(ml*(ml+1), wp)/(REAL((mk+ml+1), wp))
        pak2 = Skl/invAkl(1, 1)
        pak3 = REAL((mk+ml-2), wp)/(REAL((mk+ml+1), wp))
        pak4 = ABA11/invAkl(1, 1)

        Wyn(8) = 12.0_wp*pak1*pak2*TrAB+&
             4.0_wp*REAL(ml*(ml+1), wp)*pak3*pak2*pak4

        ACA11 = ZERO
        trAC = ZERO
        DO i = 1, n
          DO j = 1, n
            ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
            trAC = trAC+invAkl(i, j)*C(j, i)
          ENDDO
        ENDDO
        pak1 = REAL(mk*(mk+1), wp)/(REAL((mk+ml+1), wp))
        pak2 = Skl/invAkl(1, 1)
        pak3 = REAL((mk+ml-2), wp)/(REAL((mk+ml+1), wp))
        pak4 = ACA11/invAkl(1, 1)

        Wyn(7) = 12.0_wp*pak1*pak2*TrAC +&
              4.0_wp*REAL(mk*(mk+1), wp)*pak3*pak2*pak4


        !    write(*,*) '7 ii',wyn(7),stare(7)
        !    write(*,*) '8 ii',wyn(8),stare(8)
        ! begin Wyn(11)********************************


        DO i = 1, n
          DO j = 1, n
            B(i, j) = ZERO
          ENDDO
          B(1, i) = Ak(1, i)
          B(i, 1) = B(i, 1)+Ak(i, 1)
        ENDDO
        CALL r1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m2_rBrrCr_kl)
        Stare(11) = -8.0_wp*mk*r1m2_rBrrCr_kl

        CALL nower1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m2_rBrrCr_kl)
        Wyn(11) = -8.0_wp*REAL(mk, wp)*nower1m2_rBrrCr_kl


        !       write(*,*) '11 ii',wyn(11),stare(11)

        ! begin Wyn(12)********************************

        DO i = 1, n
          DO j = 1, n
            B(i, j) = ZERO
            C(i, j) = Ak(i, ii)*Ak(ii, j)
          ENDDO
          B(1, i) = Al(1, i)
          B(i, 1) = B(i, 1)+Al(i, 1)
        ENDDO
        CALL r1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m2_rBrrCr_kl)
        Stare(12) = -8.0_wp*ml*r1m2_rBrrCr_kl

        CALL nower1m2_rBrrCr(n, mk, ml, Skl, B, C, invAkl, nower1m2_rBrrCr_kl)
        Wyn(12) = -8.0_wp*ml*nower1m2_rBrrCr_kl

        !       write(*,*) '12 ii',wyn(12),stare(12)

        ! begin Wyn(13) **********************************************

        CALL r1m4(n, mk, ml, Skl, invAkl, r1m4_kl)
        Stare(13) = mk*ml*(mk+1)*(ml+1)*r1m4_kl

        pak1 = REAL(mk*(ml+1), wp)/REAL((mk+ml-1), wp)
        pak2 = pak1/invAkl(1, 1)
        pak3 = REAL(ml*(mk+1), wp)/REAL((mk+ml+1), wp)
        pak4 = pak3/invAkl(1, 1)
        pak5 = pak2*pak4

        Wyn(13) = 4.0_wp*Skl*pak5

        !   write(*,*) '13 ii',wyn(13),Stare(13)

        ! begin Wyn(14) Wyn(15)*********************************************************
        DO i = 1, n
          DO j = 1, n
            B(i, j) = ZERO
            C(i, j) = ZERO
          ENDDO
          B(1, i) = Ak(1, i)
          B(i, 1) = B(i, 1)+Ak(i, 1)
          C(1, i) = Al(1, i)
          C(i, 1) = C(i, 1)+Al(i, 1)
        ENDDO

        CALL r1m4_rBr(n, mk, ml, Skl, C, invAkl, r1m4_rBr_kl)
        Stare(14) = -2.0_wp*(mk+1)*mk*ml*r1m4_rBr_kl

        CALL r1m4_rBr(n, mk, ml, Skl, B, invAkl, r1m4_rBr_kl)
        Stare(15) = -2.0_wp*mk*(ml+1)*ml*r1m4_rBr_kl

        ACA11 = ZERO
        trAC = ZERO
        DO i = 1, n
          DO j = 1, n
            ACA11 = ACA11+invAkl(1, i)*C(i, j)*invAkl(j, 1)
            trAC = trAC+invAkl(i, j)*C(j, i)
          ENDDO
        ENDDO

        pak1 = ACA11/invAkl(1, 1)
        pak2 = Skl/invAkl(1, 1)**2
        pak3 = REAL((mk+ml-4), wp)/REAL((mk+ml-1), wp)
        pak4 = REAL((mk+1), wp)/REAL((mk+ml+1), wp)
        pak5 = REAL(mk, wp)/REAL((mk+ml-1), wp)
        pak6 = REAL(mk, wp)*REAL(ml, wp)*pak4
        pak7 = REAL(ml, wp)*pak4

        Wyn(14) = -4.0_wp*pak5*pak2*3.0_wp*TrAC*pak7 - &
             4.0_wp*pak2*pak3*pak1*pak6


        ABA11 = ZERO
        trAB = ZERO
        DO i = 1, n
          DO j = 1, n
            ABA11 = ABA11+invAkl(1, i)*B(i, j)*invAkl(j, 1)
            trAB = trAB+invAkl(i, j)*B(j, i)
          ENDDO
        ENDDO

        pak1 = ABA11/invAkl(1, 1)
        pak2 = Skl/invAkl(1, 1)**2
        pak3 = REAL((mk+ml-4), wp)/REAL((mk+ml-1), wp)
        pak4 = REAL((ml+1), wp)/REAL((mk+ml+1), wp)
        pak5 = REAL(mk, wp)/REAL((mk+ml-1), wp)
        pak6 = REAL(mk, wp)*REAL(ml, wp)*pak4
        pak7 = REAL(ml, wp)*pak4

        Wyn(15) = -4.0_wp*pak5*pak2*3.0_wp*TrAB*pak7 - &
             4.0_wp*pak2*pak3*pak1*pak6


        !    write(*,*) '14 ii',wyn(14),Stare(14)
        !    write(*,*) '15 ii',wyn(15),stare(15)

        ! Wyn(16)***************************************************************

        CALL r1m4_rBrrCr(n, mk, ml, Skl, B, C, invAkl, r1m4_rBrrCr_kl)
        Stare(16) = +4.0_wp*mk*ml*r1m4_rBrrCr_kl
        Wyn(16) = stare(16)
        !       write(*,*) '16 ii',wyn(16),stare(16)

      ENDIF


      sum111 = ZERO
      DO i = 1, 16
        sum111 = sum111 + Wyn(i)
      ENDDO

      nJiinnJiin_kl = sum111

    END FUNCTION nJiinnJiin_kl


    !
    ! end of mass-velocity correction
    !


    !

    !**************************************************************************
    !
    ! Procedures needed to claculate relativistic correction: Darvin correction
    !
    !**************************************************************************

    SUBROUTINE suma_macierzy(n, Ak, Al, Msum)

      IMPLICIT NONE
      INTEGER                          :: n, i, j
      REAL(wp), DIMENSION(1:n, 1:n) :: Ak, Al, Msum


      DO i = 1, n
        DO j = 1, n
          Msum(i, j) = Ak(i, j) + Al(i, j)
        ENDDO
      ENDDO

    END SUBROUTINE suma_macierzy


    SUBROUTINE matrix_Dkl(n, i, j, Ak, Al, Dkl)
      ! Calculate matrix Dkl (n-1)x(n-1)
      ! formed from Akl by adding the jth row to the ith one,
      ! and adding jth column to the ith column, and crossing out
      ! the jth column and row

      IMPLICIT NONE
      INTEGER, INTENT(IN)                               :: n, i, j
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n)      :: Ak, Al
      REAL(wp), INTENT(OUT), DIMENSION(1:n-1, 1:n-1) :: Dkl

      REAL(wp), DIMENSION(1:n, 1:n) :: Akl
      INTEGER                          :: ii, jj

      CALL suma_macierzy(n, Ak, Al, Akl)


      Dkl(1:j-1, 1:j-1) = Akl(1:j-1, 1:j-1)
      Dkl(j:n-1, 1:j-1) = Akl(j+1:n, 1:j-1)
      Dkl(1:j-1, j:n-1) = Akl(1:j-1, j+1:n)
      Dkl(j:n-1, j:n-1) = Akl(j+1:n, j+1:n)

      Dkl(1:j-1, i) = Dkl(1:j-1, i)+Akl(1:j-1, j)
      Dkl(i, 1:j-1) = Dkl(i, 1:j-1)+Akl(j, 1:j-1)
      Dkl(j:n-1, i) = Dkl(j:n-1, i)+Akl(j+1:n, j)
      Dkl(i, j:n-1) = Dkl(i, j:n-1)+Akl(j, j+1:n)
      Dkl(i, i) = Dkl(i, i)+Akl(j, j)


    END SUBROUTINE matrix_Dkl

    SUBROUTINE matrix_Dkl_ii(n, i, Ak, Al, Dkl_ii)
      ! Calculate matrix Dkl_ii
      ! Dkl_ii is formed from Akl by crossing aut i-th column and i-th row

      IMPLICIT NONE
      INTEGER, INTENT(IN)                               :: n, i
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n)      :: Ak, Al
      REAL(wp), INTENT(OUT), DIMENSION(1:n-1, 1:n-1) :: Dkl_ii

      REAL(wp), DIMENSION(1:n, 1:n) :: Akl

      CALL suma_macierzy(n, Ak, Al, Akl)

      Dkl_ii(1:i-1, 1:i-1) = Akl(1:i-1, 1:i-1)
      Dkl_ii(i:n-1, 1:i-1) = Akl(i+1:n, 1:i-1)
      Dkl_ii(1:i-1, i:n-1) = Akl(1:i-1, i+1:n)
      Dkl_ii(i:n-1, i:n-1) = Akl(i+1:n, i+1:n)

    END SUBROUTINE matrix_Dkl_ii

    SUBROUTINE choldc(a, n, np, wek)

      IMPLICIT NONE
      INTEGER                          :: n, np
      REAL(wp), DIMENSION(1:n, 1:n) :: a
      REAL(wp), DIMENSION(1:n)      :: wek

      INTEGER     :: i, j, k
      REAL(wp) :: sum

      DO i = 1, n
        DO j = i, n
          sum = a(i, j)
          DO k = i-1, 1, -1
            sum = sum-a(i, k)*a(j, k)
          ENDDO
          IF (i .EQ. j) THEN
            ! IF(sum.LE.0.) PAUSE 'Press to play'
            wek(i) = SQRT(sum)
          ELSE
            a(j, i) = sum/wek(i)
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE choldc


    SUBROUTINE dett(wek, n, detD)

      IMPLICIT NONE
      INTEGER                     :: n, i, j, k
      REAL(wp)                 :: detD
      REAL(wp), DIMENSION(1:n) :: wek

      detD = ONE
      DO i = 1, n
        detD = detD*wek(i)
      ENDDO
      detD = detD*detD

    END SUBROUTINE dett


    SUBROUTINE invD(a, n, wek, invDkl)

      IMPLICIT NONE
      INTEGER                          :: n, i, j, k
      REAL(wp), DIMENSION(1:n)      :: wek
      REAL(wp), DIMENSION(1:n, 1:n) :: a, invDkl
      REAL(wp)                      :: sum

      DO i = 1, n
        a(i, i) = wek(i)
      ENDDO

      DO i = 1, n
        a(i, i) = ONE/a(i, i)
        DO j = i+1, n
          sum = ZERO
          DO k = i, j-1
            sum = sum-a(j, k)*a(k, i)
          ENDDO
          a(j, i) = sum/a(j, j)
          a(i, j) = ZERO
        ENDDO
      ENDDO
      DO i = 1, n
        DO j = 1, n
          invDkl(i, j) = ZERO
        ENDDO
      ENDDO
      DO j = 1, n
        DO k = 1, n
          DO i = 1, n
            invDkl(i, j) = invDkl(i, j)+ a(k, i)*a(k, j)
          ENDDO
        ENDDO
      ENDDO
    END SUBROUTINE invD


    REAL(wp) FUNCTION delta_Rij (n, mk, ml, i, j, Ak, Al, Skl, invAkl, detAkl)
      ! Calculate matrix element:
      !    < phi_k | delta (r_ij) |  phi_l>

      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, i, j, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl, detAkl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: Ak, Al, invAkl

      INTEGER                              :: p, ii, jj
      REAL(wp), DIMENSION(1:n-1, 1:n-1) :: Dkl, workDkl, inv_Dkl
      REAL(wp), DIMENSION(1:n)          :: wek
      REAL(wp)                          :: det_Dkl
      REAL(wp)                          :: aa, dd, pak1, pak2, pak3

      p = (mk+ml)/2

      DO ii = 1, n-1
        DO jj = 1, n-1
          Dkl(ii, jj) = ZERO
        ENDDO
      ENDDO

      CALL matrix_Dkl(n, i, j, Ak, Al, Dkl)

      DO ii = 1, n-1
        workDkl(1:n-1, ii) = Dkl(1:n-1, ii)
      ENDDO

      CALL choldc(workDkl, n-1, n-1, wek)

      CALL dett(wek, n-1, det_Dkl)


      CALL invD(workDkl, n-1, wek, inv_Dkl)

      dd = inv_Dkl(1, 1)
      aa = invAkl(1, 1)
      pak1 = dd/aa
      pak2 = detAkl/det_Dkl
      pak3 = ONE/Glob_Pi

      delta_Rij = Skl*(pak2)**(3.0_wp/2.0_wp)*&
                  pak3**(3.0_wp/2.0_wp)*(pak1)**p

    END FUNCTION delta_Rij


    REAL(wp) FUNCTION delta_Ri (n, mk, ml, i, Ak, Al, Skl, invAkl, detAkl)
      ! Calculate matrix element:
      !    < phi_k | delta (r_i) |  phi_l>

      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, i, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl, detAkl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: Ak, Al, invAkl

      INTEGER                              :: p, ii, jj
      REAL(wp), DIMENSION(1:n-1, 1:n-1) :: Dkl_ii, workDkl_ii, inv_Dkl_ii
      REAL(wp), DIMENSION(1:n)          :: wek
      REAL(wp)                          :: det_Dkl_ii
      REAL(wp)                          :: aa, dd, pak1, pak2, pak3

      p = (mk+ml)/2

      CALL matrix_Dkl_ii(n, i, Ak, Al, Dkl_ii)
      DO ii = 1, n-1
        workDkl_ii(1:n-1, ii) = Dkl_ii(1:n-1, ii)
      ENDDO

      CALL choldc(workDkl_ii, n-1, n-1, wek)

      CALL dett(wek, n-1, det_Dkl_ii)

      CALL invD(workDkl_ii, n-1, wek, inv_Dkl_ii)


      dd = inv_Dkl_ii(1, 1)
      aa = invAkl(1, 1)
      pak1 = dd/aa
      pak2 = ONE/Glob_Pi
      pak3 = detAkl/det_Dkl_ii

      delta_Ri = Skl*(pak3)**(3.0_wp/2.0_wp)*&
                 pak2**(3.0_wp/2.0_wp)*(pak1)**p

    END FUNCTION delta_Ri


    REAL(wp) FUNCTION delta_R1 (n, Skl, invAkl, detAkl)
      ! Calculate matrix element:
      !    < phi_k | delta (r_1) |  phi_l>

      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n
      REAL(wp), INTENT(IN)                      :: Skl, detAkl
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: invAkl

      REAL(wp) :: aa, gmp32

      ! gmp32 = 1/(2*Glob_Pi)/Gamma(3/2)

      gmp32 = 0.17958712212516656168908198362769275528_wp
      aa = invAkl(1, 1)

      delta_R1 = Skl*gmp32/(aa**(3.0_wp/2.0_wp))

    END FUNCTION delta_R1


    !#####################################################
    !###    general SUBROUTINE for Darwin correction   ###
    !#####################################################

    ! The charge products are taken from Glob_ChargeMatrix (diagonal: q_i*q_0,
    ! strict lower triangle: q_i*q_j) instead of Glob_PseudoCharge directly,
    ! so that the REPULSION/ATTRACTION_SCALING_PARAM records of the data file
    ! act on them as in the PG_0S reference (scaled charge products).
    SUBROUTINE Darwin(n, i, j, mk, ml, cvel, Skl, Ak, Al, invAkl, Darwin_kl)
      ! Calculate matrix element Darwin correction:
      !      darwin_kl = 1.0d0_wp/(8.0d0_wp*cvel*cvel)*
      !                  ( sums )

      IMPLICIT NONE
      INTEGER, INTENT(IN)                          :: n, mk, ml
      REAL(wp), INTENT(IN)                      :: Skl, cvel
      REAL(wp), INTENT(IN), DIMENSION(1:n, 1:n) :: Ak, Al, invAkl
      REAL(wp), INTENT(OUT)                     :: darwin_kl

      INTEGER     :: p, i, j, k, l
      REAL(wp) :: sum1, sum2, sum3, Const, ilo

      Const = -Glob_Pi/(2.0_wp*cvel*cvel)
      p = (mk+ml)/2
      sum1 = ZERO
      sum2 = ZERO
      sum3 = ZERO


      IF (p .EQ. 0) THEN
        sum1 = Glob_ChargeMatrix(1, 1)*delta_R1(n, Skl, invAkl, detAkl)*&
               (ONE/(Glob_Mass(1)*Glob_Mass(1)) + ONE/(Glob_Mass(2)*Glob_Mass(2)))

      ELSE
        sum1 = ZERO
      ENDIF

      DO i = 2, n
        sum2 = sum2+Glob_ChargeMatrix(i, i)*&
             delta_Ri(n, mk, ml, i, Ak, Al, Skl, invAkl, detAkl)*&
             (ONE/(Glob_Mass(1)*Glob_Mass(1)) + ONE/(Glob_Mass(i+1)*Glob_Mass(i+1)))
      ENDDO

      DO i = 1, n
        ilo = ONE/(Glob_Mass(i+1)*Glob_Mass(i+1))
        DO j = 1, i-1
          sum3 = sum3 + ilo * Glob_ChargeMatrix(i, j)*&
               delta_Rij (n, mk, ml, j, i, Ak, Al, Skl, invAkl, detAkl)
        ENDDO
        DO j = i+1, n
          sum3 = sum3 + ilo * Glob_ChargeMatrix(j, i)*&
               delta_Rij (n, mk, ml, i, j, Ak, Al, Skl, invAkl, detAkl)
        ENDDO
      ENDDO


      Darwin_kl = Const*(sum1 + sum2 + sum3)

    END SUBROUTINE Darwin

    !
    ! end of Darvin correction
    !


    ! sorting SUBROUTINE

    SUBROUTINE sort(n, par, Input, output)
      ! sorting Input vector ans ding its terms witht given order
      ! IF par=0 no sorting is done

      IMPLICIT NONE

      INTEGER(2), INTENT(IN)                  :: n, par
      REAL(wp), DIMENSION(1:n), INTENT(IN) :: Input
      REAL(wp), INTENT(OUT)                :: output


      REAL(wp), DIMENSION(1:n) :: Dod, Uje
      LOGICAL(1), DIMENSION(1:n)  :: maska
      INTEGER(2)                  :: i, j, k, ipos
      REAL(wp)                 :: wynu, wynd, wynk


      output = 0

      IF (par == 0) THEN
        DO i = 1, n
          output = output+Input(i)
        ENDDO
        RETURN
      ENDIF

      Uje = ZERO
      Dod = ZERO
      wynd = ZERO
      wynu = ZERO
      j = 0
      k = 0

      DO i = 1, n
        IF (Input(i) < ZERO) THEN
          j = j+1
          Uje(j) = Input(i)
        ELSE
          k = k+1
          Dod(k) = Input(i)
        ENDIF
      ENDDO

      maska = .FALSE.
      maska(1:k) = .TRUE.

      DO i = 1, k
        ipos = MINLOC(Dod, 1, maska)
        maska(ipos) = .FALSE.
        wynd = wynd+Dod(ipos)
      ENDDO

      maska = .FALSE.
      maska(1:j) = .TRUE.

      DO i = 1, j
        ipos = MINLOC(ABS(Uje), 1, maska)
        maska(ipos) = .FALSE.
        wynu = wynu+Uje(ipos)
      ENDDO

      output = wynd+wynu


    END SUBROUTINE sort

    !
    ! end DK
    !


  END SUBROUTINE ExpValuesMatElem


END MODULE matelem

