!DEC$ DECLARE

MODULE matelem
  ! Module matelem contains the procedures that compute the matrix
  ! elements of the PG_1P basis
  !
  !   f_k = x_i * r_1^m_k * exp[-r'(L_k L_k' (x) I3) r]
  !
  ! (all-particle explicitly correlated Gaussians with an x_i angular
  ! prefactor and an even power of the internuclear distance r_1, for
  ! N = 1 rotationally excited states of diatomic systems):
  !   MatrixElements   - overlap, Hamiltonian and their derivatives with
  !                      respect to vech L_k and vech L_l
  !   MatrixElementsXi - the same for one prefactor x_j of the symmetry
  !                      transformed ket (called by MatrixElements)
  !   ExpValuesMatElem - the nucleus-nucleus correlation function and the
  !                      interparticle distances <r_ij>, <r_ij^2>
  !   ExpValuesXi      - the same for one prefactor x_j of the symmetry
  !                      transformed ket (called by ExpValuesMatElem)
  !   KetNormRatios    - normalization ratios of the ket terms when the
  !                      symmetry operator turns x_il into a sum of x_j
  !   GramTermPair, GramTermCorr - the Gram-determinant term of the
  !                      Coulomb integrals and its gradient (see
  !                      MatrixElementsXi); from 5 pseudoparticles on, in
  !                      double precision, GramTermSetup, GramTermCorrPairs
  !                      instead
  !   PowerZeroGradCorr - correction of the Coulomb gradient for mk = ml = 0
  !   CoulombCorrSetup, CoulombCorrDerivs, AddCorrToGrad - their helpers
  !   MatrixElements0  - overlap and Hamiltonian of the N = 0 basis
  !                      r_1^m_k * exp[...] (kept from the original code,
  !                      not called by this program)
  !   MatElemTraMom, TMMatElem - kernels of the N = 0 -> N = 1 transition
  !                      moments <psi(N=0)|x_i|psi(N=1)> (kept from the
  !                      original code, not called by this program)
  ! The Coulomb terms take the charge products from Glob_ChargeMatrix, into
  ! which ProgramDataInit has folded the optional REPULSION_SCALING_PARAM
  ! and ATTRACTION_SCALING_PARAM factors of the data file.

  USE globvars
  USE wp_def
  USE mpi

  IMPLICIT NONE

  ! The power sums of all pseudoparticle pairs are advanced together
  ! (CoulombPowerSums, GramTermCorrPairs) in double precision only: in x87
  ! extended (wp = 10) and software quad (wp = 16) arithmetic this is
  ! slower than the original loops over one pair at a time.
  LOGICAL, PARAMETER :: PairedPowerSums = (wp == 8)

CONTAINS


  SUBROUTINE MatrixElements(mk, ik, vechLk, ml, il, vechLl, SymMatrix, &
                            Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! Symmetry adapted matrix elements in the basis of correlated Gaussians
    ! with x_i r_1^m premultipliers for N = 1 rotationally excited states,
    !   f_k = x_i * r_1^mk * exp[-r'(Lk*Lk' (x) I3) r].
    ! The symmetry operator acts on the whole ket, P f_l(r) = f_l(SymMatrix*r):
    ! the exponent becomes Al = SymMatrix'*Ll*Ll'*SymMatrix and the
    ! prefactor x_il becomes sum_j SymMatrix(il, j) x_j (r_1 is unchanged,
    ! row 1 of SymMatrix is +-e_1). MatrixElementsXi computes the matrix
    ! elements for one prefactor x_j of the transformed ket:
    !  - row il of SymMatrix has one entry +-1 (every permutation of the
    !    particles 2 .. N, and x_1 under the exchange of the two nuclei):
    !    one call with j = its column, the sign applied to all outputs;
    !  - otherwise (x_il, il > 1, under the exchange of the two nuclei:
    !    x_il --> x_il - x_1) the terms are summed. MatrixElementsXi
    !    normalizes the ket with its own x_j, so every term is rescaled
    !    to the norm of the true ket (KetNormRatios).
    ! For the identity term this is one call with j = il, as before.
    ! GramTermCorr (GramTermCorrPairs from 5 pseudoparticles on in double
    ! precision, PairedPowerSums) then adds
    ! the Gram-determinant term of the Coulomb integrals where it can be
    ! present (see MatrixElementsXi), and for mk = ml = 0 PowerZeroGradCorr
    ! corrects the Coulomb gradients.
    ! Input:
    !   mk, ml         :: premultiplier powers (even) of the bra and the ket
    !   ik, il         :: pseudoparticle indices i of the x_i prefactors of
    !                     the bra and the ket (1 .. n)
    !   vechLk, vechLl :: n(n+1)/2 exponent parameters (Lk, Ll lower triangular)
    !   SymMatrix      :: the symmetry permutation matrix
    !   gradflag       :: 0 no gradient; 1 dSk, dHk only; 2 also dSl, dHl
    ! Output (normalized):
    !   Skl, Hkl       :: overlap and Hamiltonian matrix elements
    !   dSk, dHk       :: derivatives with respect to vech[Lk]
    !   dSl, dHl       :: derivatives with respect to vech[Ll]
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, il, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: Skl, Hkl
    REAL(wp), INTENT(OUT) :: dSk(Glob_npt), dSl(Glob_npt), dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: j, t, NumKetTerms
    INTEGER :: KetIndex(n)                        ! x-index j of each ket term
    REAL(wp) :: KetWeight(n)                      ! SymMatrix(il, j) of each ket term
    REAL(wp) :: KetRatio(n)                       ! sqrt(N_j/N_ket), see KetNormRatios
    REAL(wp) :: dlnKetRatio(Glob_npt, n)          ! its log-derivative, vech[Ll]
    REAL(wp) :: f, St, Ht
    LOGICAL :: IsGramPossible
    REAL(wp) :: dSkt(Glob_npt), dSlt(Glob_npt), dHkt(Glob_npt), dHlt(Glob_npt)


    ! The terms of the transformed prefactor: row il of SymMatrix
    NumKetTerms = 0
    DO j = 1, n
      IF (SymMatrix(il, j) /= ZERO) THEN
        NumKetTerms = NumKetTerms+1
        KetIndex(NumKetTerms) = j
        KetWeight(NumKetTerms) = SymMatrix(il, j)
      ENDIF
    ENDDO

    ! One term: x_il --> +-x_j
    IF ((NumKetTerms == 1) .AND. (ABS(KetWeight(1)) == ONE)) THEN
      CALL MatrixElementsXi(mk, ik, vechLk, ml, KetIndex(1), vechLl, SymMatrix, &
                            Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
      IsGramPossible = (mk+ml > 0) .AND. (ik /= 1) .AND. (KetIndex(1) /= 1) .AND. (n > 2)
      IF (IsGramPossible) THEN
        IF ((n >= 5) .AND. PairedPowerSums) THEN
          CALL GramTermCorrPairs(mk, ik, vechLk, ml, KetIndex(1), vechLl, SymMatrix, &
                                 Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
        ELSE
          CALL GramTermCorr(mk, ik, vechLk, ml, KetIndex(1), vechLl, SymMatrix, &
                            Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
        ENDIF
      ENDIF
      IF ((mk+ml == 0) .AND. (gradflag > 0)) CALL PowerZeroGradCorr(ik, KetIndex(1), vechLk, vechLl, &
                                                                SymMatrix, Skl, dHk, dHl, gradflag)
      IF (KetWeight(1) < ZERO) THEN
        Skl = -Skl
        Hkl = -Hkl
        dSk = -dSk
        dSl = -dSl
        dHk = -dHk
        dHl = -dHl
      ENDIF
      RETURN
    ENDIF

    ! Several terms: x_il --> sum_t KetWeight(t) x_KetIndex(t)
    CALL KetNormRatios(ml, il, vechLl, SymMatrix, NumKetTerms, KetIndex, &
                       KetRatio, dlnKetRatio, (gradflag > 1))

    Skl = ZERO
    Hkl = ZERO
    dSk = ZERO
    dSl = ZERO
    dHk = ZERO
    dHl = ZERO
    DO t = 1, NumKetTerms
      CALL MatrixElementsXi(mk, ik, vechLk, ml, KetIndex(t), vechLl, SymMatrix, &
                            St, Ht, dSkt, dSlt, dHkt, dHlt, gradflag)
      IsGramPossible = (mk+ml > 0) .AND. (ik /= 1) .AND. (KetIndex(t) /= 1) .AND. (n > 2)
      IF (IsGramPossible) THEN
        IF ((n >= 5) .AND. PairedPowerSums) THEN
          CALL GramTermCorrPairs(mk, ik, vechLk, ml, KetIndex(t), vechLl, SymMatrix, &
                                 St, Ht, dSkt, dSlt, dHkt, dHlt, gradflag)
        ELSE
          CALL GramTermCorr(mk, ik, vechLk, ml, KetIndex(t), vechLl, SymMatrix, &
                            St, Ht, dSkt, dSlt, dHkt, dHlt, gradflag)
        ENDIF
      ENDIF
      IF ((mk+ml == 0) .AND. (gradflag > 0)) CALL PowerZeroGradCorr(ik, KetIndex(t), vechLk, vechLl, &
                                                                SymMatrix, St, dHkt, dHlt, gradflag)
      f = KetWeight(t)*KetRatio(t)
      Skl = Skl+f*St
      Hkl = Hkl+f*Ht
      dSk = dSk+f*dSkt
      dHk = dHk+f*dHkt
      IF (gradflag > 1) THEN
        dSl = dSl+f*(dSlt+St*dlnKetRatio(1:Glob_npt, t))
        dHl = dHl+f*(dHlt+Ht*dlnKetRatio(1:Glob_npt, t))
      ENDIF
    ENDDO

  END SUBROUTINE MatrixElements





  SUBROUTINE MatrixElementsXi(mk, ik, vechLk, ml, jl, vechLl, SymMatrix, &
                              Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! Symmetry adapted matrix elements in the basis of correlated Gaussians
    ! with x_i r_1^m premultipliers for N = 1 rotationally excited states,
    !   f_k = x_i * r_1^mk * exp[-r'(Lk*Lk' (x) I3) r],
    ! for one term x_jl of the prefactor of the symmetry transformed ket
    ! (see MatrixElements): the ket is x_jl * r_1^ml * exp[-r'(Al (x) I3) r],
    ! Al = SymMatrix'*Ll*Ll'*SymMatrix, normalized with its own x_jl.
    !
    ! The Coulomb terms follow eqs. (33)-(35) of Kirnosov, Sharkey,
    ! Adamowicz, J. Phys. B 48, 195101 (2015), summed over all pairs of
    ! eq. (34). Eq. (31) of the paper leaves out a term: the inverse of the
    ! rank-2 update A + u J11 + t^2 J_pg has, besides the terms of eq. (31),
    ! the term u t^2 G_ij/det, where G_ij is the 3x3 Gram determinant of
    ! (e1, v, e_i) and (e1, v, e_j) in the metric invAkl (v = e_p for 1/r_p,
    ! e_p - e_g for 1/r_pg). G_ij vanishes when e_i or e_j lies in
    ! span{e1, v} - always for n = 2, and whenever the bra or the ket
    ! prefactor is x_1 - and the formulas of the paper are then exact.
    ! Otherwise MatrixElements adds the term with GramTermCorr (see
    ! GramTermPair); it is not computed here, so that the code of this
    ! routine and its results where the term is absent stay as they were.
    ! Input:
    !   mk, ml         :: premultiplier powers (even) of the bra and the ket
    !   ik             :: pseudoparticle index i of the x_i prefactor of the bra
    !   jl             :: pseudoparticle index j of the x_j term of the
    !                     transformed ket (1 .. n)
    !   vechLk, vechLl :: n(n+1)/2 exponent parameters (Lk, Ll lower triangular)
    !   SymMatrix      :: the symmetry permutation matrix
    !   gradflag       :: 0 no gradient; 1 dSk, dHk only; 2 also dSl, dHl
    ! Output (normalized):
    !   Skl, Hkl       :: overlap and Hamiltonian matrix elements
    !   dSk, dHk       :: derivatives with respect to vech[Lk]
    !   dSl, dHl       :: derivatives with respect to vech[Ll]
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, jl, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: Skl, Hkl
    REAL(wp), INTENT(OUT) :: dSk(Glob_npt), dSl(Glob_npt), dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: i, j, k, m, p, q, t, ii, jj, kk, indx, Msum, MMsum, M1, M2, M3, mko2, mlo2
    INTEGER :: parti, partj
    REAL(wp) :: mkml, mkpml, temp1, temp2, temp3, temp4, temp5, temp6
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3, chi1, A11, A1i, A1j, Aij, oA11
    REAL(wp) :: detLk, detLl, detAkl, a, coabij, g1, g2, g3, g4, d, h
    REAL(wp) :: tau1, tau2, tau3, Tkl, Vkl, hoa
    REAL(wp) :: XMij, XM11, XM1j, XM1i, AAkMij, AAkM1j, AAkMj1, AAkMi1, AAkM11
    REAL(wp) :: AAlMji, AAlMj1, AAlM1i, AAlMi1, AAlM11
    REAL(wp) :: Lk(n, n), Ll(n, n), Ak(n, n), Al(n, n)
    REAL(wp) :: invAkl(n, n), PLl(n, n), invL(n, n), Akl(n, n)
    REAL(wp) :: Temp(n, n), R(n, n), b(n, n), c(n, n), coab(n, n), Vtmp(n, n)
    REAL(wp) :: fob(n, n), g5(n, n)
    REAL(wp) :: dkhoa(n, n), dkb(n, n), dkfob(n, n), dkcoab(n, n)
    REAL(wp) :: dlhoa(n, n), dlb(n, n), dlfob(n, n), dlcoab(n, n)
    REAL(wp) :: dTk(Glob_npt), dTl(Glob_npt), dVk(Glob_npt), dVl(Glob_npt)
    REAL(wp) :: ddSk(n, n), ddSl(n, n), ddTk(n, n), ddTl(n, n)
    REAL(wp) :: ddSk_over_Skl(n, n), ddSl_over_Skl(n, n), Templ(n, n), TemplT(n, n)
    REAL(wp) :: dRk(n, n), dRl(n, n), ddVk(n, n), ddVl(n, n)
    REAL(wp) :: invLk(n, n), invLl(n, n), invAk(n, n), invAl(n, n), chAl(n, n)
    REAL(wp) :: ALk(n, n), ALl(n, n), AAk(n, n), AAl(n, n), AlLl(n, n)
    REAL(wp) :: dkA11(n, n), dkA1j(n, n), dkA1i(n, n), dkAij(n, n), dkApq(n, n)
    REAL(wp) :: dlA11(n, n), dlA1j(n, n), dlA1i(n, n), dlAij(n, n), dlApq(n, n)
    REAL(wp) :: dkA(n, n, n, n), dlA(n, n, n, n)
    REAL(wp) :: dkAk11(n, n), dkAk1i(n, n), dkAkii(n, n)  ! ,dkAkij(n,n)
    REAL(wp) :: dlAl11(n, n), dlAl1j(n, n), dlAljj(n, n)  ! ,dlAlij(n,n)
    REAL(wp) :: MAk(n, n), MAl(n, n), Jkk(n, n), Jll(n, n), Jklk(n, n), Jkll(n, n)
    REAL(wp) :: MAlALk(n, n), AAlM(n, n)
    REAL(wp) :: MAlALl(n, n), MAkALl(n, n), MLl(n, n), AAkMAlALl(n, n)
    REAL(wp) :: AAlMAkALl(n, n), AAkMLl(n, n)
    REAL(wp) :: dkAAlM1i(n, n), dkAAlM11(n, n), dkAAlMi1(n, n), dkAAlMj1(n, n)
    REAL(wp) :: dkAAlMji(n, n)
    REAL(wp) :: dkAAkMij(n, n), dkAAkM1j(n, n), dkAAkM11(n, n), dkAAkMi1(n, n)
    REAL(wp) :: dkAAkMj1(n, n)
    REAL(wp) :: AAlMAkA(n, n), AAlMAkALk(n, n), AAkMAlA(n, n), AAkMAlALk(n, n)
    REAL(wp) :: dkXij(n, n), dkX11(n, n), dkX1j(n, n), dkX1i(n, n), dktau1(n, n)
    REAL(wp) :: dlAAlM1i(n, n), dlAAlM11(n, n), dlAAlMi1(n, n), dlAAlMj1(n, n)
    REAL(wp) :: dlAAlMji(n, n)
    REAL(wp) :: dlAAkMij(n, n), dlAAkM1j(n, n), dlAAkM11(n, n), dlAAkMi1(n, n)
    REAL(wp) :: dlAAkMj1(n, n)
    REAL(wp) :: dlXij(n, n), dlX11(n, n), dlX1j(n, n), dlX1i(n, n), dltau1(n, n)
    REAL(wp) :: AAkM(n, n), MAkALk(n, n), MLk(n, n), AAlMLk(n, n)
    REAL(wp) :: AAlMAl(n, n), AAkMAk(n, n), AAkMAl(n, n), AAlMAk(n, n)
    REAL(wp) :: da(n, n), db(n, n), dc(n, n), dtau1lM(n, n), dtau3lM(n, n)
    REAL(wp) :: dSkvech(Glob_npt), dSlvech(Glob_npt)
    REAL(wp) :: dtau1k(Glob_npt), dtau1l(Glob_npt)  ! ,dtau3k(NParam),dtau3l(NParam)


    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Build the lower triangular matrices Lk and Ll from
    ! the elements in vechLk and vechLl
    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
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

    ! Build AK=LK*Lk' and AL=P'*LL*LL'*P
    Ak = ZERO
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Ak(i, j) = Ak(i, j)+Lk(i, k)*Lk(j, k)
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
      DO i = 1, n
        Akl(i, j) = Ak(i, j)+Al(i, j)
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
    ! need to do the decomp on Al too to get invAl since
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
    ! Now do invAkl=invL'*invL
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

    !********************************************************
    ! Compute Skl
    !********************************************************

    parti = ik
    partj = jl

    mko2 = mk/2
    mlo2 = ml/2

    mkpml = ONE*(mk+ml)

    Skl = ZERO

    A11 = invAkl(1, 1)
    oA11 = ONE/A11
    A1i = invAkl(1, parti)
    A1j = invAkl(1, partj)
    Aij = invAkl(parti, partj)

    Skl = Skl+(invAk(parti, parti)+(TWO*mk)/THREE*invAk(1, parti)*invAk(1, parti)/invAk(1, 1))
    Skl = Skl*(invAl(partj, partj)+(TWO*ml)/THREE*invAl(1, partj)*invAl(1, partj)/invAl(1, 1))
    Skl = TWO/SQRT(Skl)
    Skl = Skl*(Aij+mkpml*ONETHIRD*A1i*A1j*oA11)
    Skl = Skl*SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))
    Skl = Skl*((A11/invAk(1, 1))**(mko2))*((A11/invAl(1, 1))**(mlo2))
    Skl = Skl*Glob_gm1(mko2, mlo2)*Glob_gm1(mko2, mlo2)

    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in Glob_gm1 (module globvars)
    !
    ! gj1(mk,ml) = (3+mk+ml)/sqrt(9+6*(mk+ml)+4*mk*ml)

    !********************************************************
    ! gradient terms dSk dSl
    !********************************************************
    IF (gradflag > 0) THEN

      ! ALk=invAkl*Lk
      ALk = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ALk(i, j) = ALk(i, j) +invAkl(i, k)*Lk(k, j)
          ENDDO
        ENDDO
      ENDDO

      ! d(Aij)/dLk
      DO j = 1, n
        dkA11(1:n, j) = -TWO*invAkl(1:n, 1)*ALk(1, j)
        dkA1j(1:n, j) = -invAkl(1:n, partj)*ALk(1, j)-invAkl(1:n, 1)*ALk(partj, j)
        dkA1i(1:n, j) = -invAkl(1:n, 1)*ALk(parti, j)-invAkl(1:n, parti)*ALk(1, j)
        dkAij(1:n, j) = -invAkl(1:n, partj)*ALk(parti, j)-invAkl(1:n, parti)*ALk(partj, j)

        dkAk11(1:n, j) = -TWO*invAk(1:n, 1)*invLk(j, 1)
        dkAk1i(1:n, j) = -invAk(1:n, 1)*invLk(j, parti)-invAk(1:n, parti)*invLk(j, 1)
        dkAkii(1:n, j) = -TWO*invAk(1:n, parti)*invLk(j, parti)
      ENDDO

      DO p = 1, n
        DO q = 1, n
          DO j = 1, n
            dkApq(1:n, j) = -invAkl(1:n, q)*ALk(p, j)-invAkl(1:n, p)*ALk(q, j)
          ENDDO
          dkA(p, q, :, :) = dkApq
        ENDDO
      ENDDO

      ! dSk
      DO j = 1, n
        DO i = j, n
          ddSk(i, j) = (ONE+ONEHALF)*invLk(j, i)-THREE*ALk(i, j)
        ENDDO
      ENDDO
      ddSk = ddSk+(mko2+mlo2)*oA11*dkA11-mko2/invAk(1, 1)*dkAk11
      ddSk = ddSk+ONE/(Aij+mkpml*ONETHIRD*&
          A1i*A1j*oA11)*&
          (dkAij+mkpml*ONETHIRD*A1i*oA11*&
          (A1j/A1i*dkA1i+dkA1j-&
          A1j*oA11*dkA11))
      ddSk = ddSk-ONEHALF/(invAk(parti, parti)+TWO*mk*ONETHIRD*&
          invAk(1, parti)*invAk(1, parti)/invAk(1, 1))*&
          (dkAkii+TWO*mk*ONETHIRD*invAk(1, parti)/invAk(1, 1)*&
          (TWO*dkAk1i-invAk(1, parti)/invAk(1, 1)*dkAk11))

      ddSk_over_Skl = ddSk
      ddSk = Skl*ddSk
      ! load it into dSk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dSk(indx) = ddSk(i, j)
        ENDDO
      ENDDO
    ENDIF

    IF (gradflag > 1) THEN

      ! ALl=invAkl*P'*Ll
      ! AlLl=invAl*P'*Ll need this because of sym
      ! AAl=invAkl*Al
      ALl = ZERO
      AlLl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ALl(i, j) = ALl(i, j) +invAkl(i, k)*PLl(k, j)
            AlLl(i, j) = AlLl(i, j)+invAl(i, k)*PLl(k, j)
          ENDDO
        ENDDO
      ENDDO

      ! Jll = invAl*J11*invAl*Ll
      ! Jkll = invAkl*J11*invAkl*Ll
      DO j = 1, n
        Jll(1:n, j) = invAl(1:n, 1)*AlLl(1, j)  ! using AlLl to get sym right
        Jkll(1:n, j) = invAkl(1:n, 1)*ALl(1, j)
      ENDDO

      ! d(Aij)/dLl
      DO j = 1, n
        dlA11(1:n, j) = -TWO*invAkl(1:n, 1)*ALl(1, j)
        dlA1j(1:n, j) = -invAkl(1:n, partj)*ALl(1, j)-invAkl(1:n, 1)*ALl(partj, j)
        dlA1i(1:n, j) = -invAkl(1:n, 1)*ALl(parti, j)-invAkl(1:n, parti)*ALl(1, j)
        dlAij(1:n, j) = -invAkl(1:n, partj)*ALl(parti, j)-invAkl(1:n, parti)*ALl(partj, j)

        ! dlAlij(1:n,j)=-invAl(1:n,partj)*AlLl(parti,j)-invAl(1:n,parti)*AlLl(partj,j)
        dlAl11(1:n, j) = -TWO*invAl(1:n, 1)*AlLl(1, j)
        dlAl1j(1:n, j) = -invAl(1:n, partj)*AlLl(1, j)-invAl(1:n, 1)*AlLl(partj, j)
        dlAljj(1:n, j) = -TWO*invAl(1:n, partj)*AlLl(partj, j)
      ENDDO

      DO p = 1, n
        DO q = 1, n
          DO j = 1, n
            dlApq(1:n, j) = -invAkl(1:n, q)*ALl(p, j)-invAkl(1:n, p)*ALl(q, j)
          ENDDO
          dlA(p, q, :, :) = dlApq
        ENDDO
      ENDDO

      ! dSl

      ! Multiply by SymMatrix
      ddSl = ZERO
      Templ = (ONE+ONEHALF)*AlLl-THREE*ALl
      Templ = Templ+(mko2+mlo2)*oA11*dlA11-mlo2/invAl(1, 1)*dlAl11
      Templ = Templ+ONE/(Aij+mkpml*ONETHIRD*&
          A1i*A1j*oA11)*&
          (dlAij+mkpml*ONETHIRD*A1i*oA11*&
          (A1j/A1i*dlA1i+dlA1j-&
          A1j*oA11*dlA11))
      Templ = Templ-ONEHALF/(invAl(partj, partj)+TWO*ml*ONETHIRD*&
          invAl(1, partj)*invAl(1, partj)/invAl(1, 1))*&
          (dlAljj+TWO*ml*ONETHIRD*invAl(1, partj)/invAl(1, 1)*&
          (TWO*dlAl1j-invAl(1, partj)/invAl(1, 1)*dlAl11))
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ddSl(i, j) = ddSl(i, j)+SymMatrix(i, k)*Templ(k, j)
          ENDDO
        ENDDO
      ENDDO


      ddSl_over_Skl = ddSl
      ddSl = Skl*ddSl
      ! Load it into dSl (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dSl(indx) = ddSl(i, j)
        ENDDO
      ENDDO
    ENDIF

    !*********************
    ! Kinetic Specific Terms
    !**********************

    ! AAk=invAkl*Ak
    ! AAl=invAkl*Al
    ! MAk=M*Ak
    ! MAl=M*Al
    AAk = ZERO
    AAl = ZERO
    MAk = ZERO
    MAl = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          AAk(i, j) = AAk(i, j)+invAkl(i, k)*Ak(k, j)
          AAl(i, j) = AAl(i, j)+invAkl(i, k)*Al(k, j)
          MAk(i, j) = MAk(i, j)+Glob_MassMatrix(i, k)*Ak(k, j)
          MAl(i, j) = MAl(i, j)+Glob_MassMatrix(i, k)*Al(k, j)
        ENDDO
      ENDDO
    ENDDO

    AAkMAl = ZERO
    AAlMAk = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          AAkMAl(i, j) = AAkMAl(i, j)+AAk(i, k)*MAl(k, j)
          AAlMAk(i, j) = AAlMAk(i, j)+AAl(i, k)*MAk(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! A few kinetic energy terms (scalers)
    ! tau1=tr[invAkl*Ak*M*Al]

    tau1 = ZERO
    XMij = ZERO
    XM11 = ZERO
    XM1i = ZERO
    XM1j = ZERO
    AAkMij = ZERO
    AAkM1j = ZERO
    AAkMj1 = ZERO
    AAkMi1 = ZERO
    AAkM11 = ZERO

    AAlMji = ZERO
    AAlMj1 = ZERO
    AAlM1i = ZERO
    AAlMi1 = ZERO
    AAlM11 = ZERO

    DO i = 1, n
      tau1 = tau1+AAkMAl(i, i)
      XMij = XMij+AAkMAl(parti, i)*invAkl(i, partj)+AAkMAl(partj, i)*invAkl(i, parti)
      XM11 = XM11+TWO*AAkMAl(1, i)*invAkl(i, 1)
      XM1i = XM1i+AAkMAl(parti, i)*invAkl(i, 1)+AAkMAl(1, i)*invAkl(i, parti)
      XM1j = XM1j+AAkMAl(1, i)*invAkl(i, partj)+AAkMAl(partj, i)*invAkl(i, 1)
      AAkMij = AAkMij+AAk(parti, i)*Glob_MassMatrix(i, partj)
      AAkM1j = AAkM1j+AAk(1, i)*Glob_MassMatrix(i, partj)
      AAkMj1 = AAkMj1+AAk(partj, i)*Glob_MassMatrix(i, 1)
      AAkMi1 = AAkMi1+AAk(parti, i)*Glob_MassMatrix(i, 1)
      AAkM11 = AAkM11+AAk(1, i)*Glob_MassMatrix(i, 1)
      AAlMji = AAlMji+AAl(partj, i)*Glob_MassMatrix(i, parti)
      AAlM1i = AAlM1i+AAl(1, i)*Glob_MassMatrix(i, parti)
      AAlMi1 = AAlMi1+AAl(parti, i)*Glob_MassMatrix(i, 1)
      AAlMj1 = AAlMj1+AAl(partj, i)*Glob_MassMatrix(i, 1)
      AAlM11 = AAlM11+AAl(1, i)*Glob_MassMatrix(i, 1)
    ENDDO

    !********************************************************
    ! Compute Tkl
    !********************************************************
    mkml = mk*ml

    AAlM = ZERO
    AAkM = ZERO
    AAlMAkA = ZERO
    AAkMAlA = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          AAlM(i, j) = AAlM(i, j)+AAl(i, k)*Glob_MassMatrix(k, j)
          AAkM(i, j) = AAkM(i, j)+AAk(i, k)*Glob_MassMatrix(k, j)
          AAlMAkA(i, j) = AAlMAkA(i, j)+AAlMAk(i, k)*invAkl(k, j)
          AAkMAlA(i, j) = AAkMAlA(i, j)+AAkMAl(i, k)*invAkl(k, j)
        ENDDO
      ENDDO
    ENDDO

    Tkl = TWO*Glob_MassMatrix(parti, partj)
    Tkl = Tkl+TWO*ONETHIRD*oA11*(ml*A1j*Glob_MassMatrix(1, parti)+&
        mk*A1i*Glob_MassMatrix(1, partj))
    Tkl = Tkl-TWO*(AAlMji+AAkMij+(mkpml*ONETHIRD*oA11)*(A1j*AAlM1i+A1i*AAkM1j))
    Tkl = Tkl+mkml*(Glob_MassMatrix(1, 1)/((mkpml+ONE)*A11*ONEHALF)*(Aij+&
        (mkpml-2)*A1i*A1j*oA11*ONETHIRD) )
    Tkl = Tkl-FOUR/(ONE+mkpml)*(mkpml+1)/(SIX*A11*A11)*&
        (A11*(THREE*Aij*(mk*AAlM(1, 1)+ml*AAkM(1, 1))+&
        A1j*(mk*AAlM(parti, 1)+ml*AAkM(parti, 1)))+&
        A1i*(A11*(mk*AAlM(partj, 1)+ml*AAkM(partj, 1))+&
        A1j*(mk*AAlM(1, 1)+ml*AAkM(1, 1))*(mkpml-TWO)))
    Tkl = Tkl+FOUR/(SIX*A11*A11)*&
        (mkpml*(mkpml-2)*ONEHALF*A1i*A1j*(AAkMAlA(1, 1)+AAlMAkA(1, 1))+&
        mkpml*A11*(A1j*(AAkMAlA(1, parti)+AAlMAkA(1, parti))+&
        A1i*(AAkMAlA(1, partj)+AAlMAkA(1, partj))+&
        THREE*ONEHALF*(Aij*(AAlMAkA(1, 1)+AAkMAlA(1, 1))+&
        TWO*A1j*A1i*tau1))+THREE*invAkl(1, 1)*invAkl(1, 1)*&
        (AAkMAlA(parti, partj)+AAlMAkA(parti, partj)+THREE*Aij*tau1))
    Tkl = Tkl*Skl/(Aij+mkpml*A1i*A1j*oA11*ONETHIRD)

    IF (gradflag > 0) THEN

      MAlALk = ZERO
      MAkALk = ZERO
      MLk = ZERO
      AAlMAkALk = ZERO
      AAkMAlALk = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            MAlALk(i, j) = MAlALk(i, j)+MAl(i, k)*ALk(k, j)
            MAkALk(i, j) = MAkALk(i, j)+MAk(i, k)*ALk(k, j)
            MLk(i, j) = MLk(i, j)+Glob_MassMatrix(i, k)*Lk(k, j)
            AAlMAkALk(i, j) = AAlMAkALk(i, j)+AAlMAk(i, k)*ALk(k, j)
            AAkMAlALk(i, j) = AAkMAlALk(i, j)+AAkMAl(i, k)*ALk(k, j)
          ENDDO
        ENDDO
      ENDDO

      AAlMLk = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            AAlMLk(i, j) = AAlMLk(i, j)+AAl(i, k)*MLk(k, j)
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, n
        dkAAlMji(1:n, j) = -AAlM(1:n, parti)*ALk(partj, j)-invAkl(1:n, partj)*MAlALk(parti, j)
        dkAAlM1i(1:n, j) = -AAlM(1:n, parti)*ALk(1, j)-invAkl(1:n, 1)*MAlALk(parti, j)
        dkAAlM11(1:n, j) = -AAlM(1:n, 1)*ALk(1, j)-invAkl(1:n, 1)*MAlALk(1, j)
        dkAAlMi1(1:n, j) = -AAlM(1:n, 1)*ALk(parti, j)-invAkl(1:n, parti)*MAlALk(1, j)
        dkAAlMj1(1:n, j) = -AAlM(1:n, 1)*ALk(partj, j)-invAkl(1:n, partj)*MAlALk(1, j)

        dkAAkMij(1:n, j) = -AAkM(1:n, partj)*ALk(parti, j)-invAkl(1:n, parti)*MAkALk(partj, j)+&
            invAkl(1:n, parti)*MLk(partj, j)+Glob_MassMatrix(1:n, partj)*ALk(parti, j)
        dkAAkM1j(1:n, j) = -AAkM(1:n, partj)*ALk(1, j)-invAkl(1:n, 1)*MAkALk(partj, j)+&
            invAkl(1:n, 1)*MLk(partj, j)+Glob_MassMatrix(1:n, partj)*ALk(1, j)
        dkAAkM11(1:n, j) = -AAkM(1:n, 1)*ALk(1, j)-invAkl(1:n, 1)*MAkALk(1, j)+&
            invAkl(1:n, 1)*MLk(1, j)+Glob_MassMatrix(1:n, 1)*ALk(1, j)
        dkAAkMi1(1:n, j) = -AAkM(1:n, 1)*ALk(parti, j)-invAkl(1:n, parti)*MAkALk(1, j)+&
            invAkl(1:n, parti)*MLk(1, j)+Glob_MassMatrix(1:n, 1)*ALk(parti, j)
        dkAAkMj1(1:n, j) = -AAkM(1:n, 1)*ALk(partj, j)-invAkl(1:n, partj)*MAkALk(1, j)+&
            invAkl(1:n, partj)*MLk(1, j)+Glob_MassMatrix(1:n, 1)*ALk(partj, j)

        dkXij(1:n, j) = -invAkl(1:n, parti)*AAlMAkALk(partj, j)-invAkl(1:n, partj)*AAlMAkALk(parti, j)-&
            AAkMAlA(1:n, parti)*ALk(partj, j)-AAkMAlA(1:n, partj)*ALk(parti, j)+&
            invAkl(1:n, parti)*AAlMLk(partj, j)+invAkl(1:n, partj)*AAlMLk(parti, j)+&
            AAlM(parti, 1:n)*ALk(partj, j)+AAlM(partj, 1:n)*ALk(parti, j)-&
            AAlMAkA(1:n, parti)*ALk(partj, j)-AAlMAkA(1:n, partj)*ALk(parti, j)-&
            invAkl(1:n, parti)*AAkMAlALk(partj, j)-invAkl(1:n, partj)*AAkMAlALk(parti, j)
        dkX11(1:n, j) = -invAkl(1:n, 1)*AAlMAkALk(1, j)-invAkl(1:n, 1)*AAlMAkALk(1, j)-&
            AAkMAlA(1:n, 1)*ALk(1, j)-AAkMAlA(1:n, 1)*ALk(1, j)+&
            invAkl(1:n, 1)*AAlMLk(1, j)+invAkl(1:n, 1)*AAlMLk(1, j)+&
            AAlM(1, 1:n)*ALk(1, j)+AAlM(1, 1:n)*ALk(1, j)-&
            AAlMAkA(1:n, 1)*ALk(1, j)-AAlMAkA(1:n, 1)*ALk(1, j)-&
            invAkl(1:n, 1)*AAkMAlALk(1, j)-invAkl(1:n, 1)*AAkMAlALk(1, j)
        dkX1j(1:n, j) = -invAkl(1:n, 1)*AAlMAkALk(partj, j)-invAkl(1:n, partj)*AAlMAkALk(1, j)-&
            AAkMAlA(1:n, 1)*ALk(partj, j)-AAkMAlA(1:n, partj)*ALk(1, j)+&
            invAkl(1:n, 1)*AAlMLk(partj, j)+invAkl(1:n, partj)*AAlMLk(1, j)+&
            AAlM(1, 1:n)*ALk(partj, j)+AAlM(partj, 1:n)*ALk(1, j)-&
            AAlMAkA(1:n, 1)*ALk(partj, j)-AAlMAkA(1:n, partj)*ALk(1, j)-&
            invAkl(1:n, 1)*AAkMAlALk(partj, j)-invAkl(1:n, partj)*AAkMAlALk(1, j)
        dkX1i(1:n, j) = -invAkl(1:n, parti)*AAlMAkALk(1, j)-invAkl(1:n, 1)*AAlMAkALk(parti, j)-&
            AAkMAlA(1:n, parti)*ALk(1, j)-AAkMAlA(1:n, 1)*ALk(parti, j)+&
            invAkl(1:n, parti)*AAlMLk(1, j)+invAkl(1:n, 1)*AAlMLk(parti, j)+&
            AAlM(parti, 1:n)*ALk(1, j)+AAlM(1, 1:n)*ALk(parti, j)-&
            AAlMAkA(1:n, parti)*ALk(1, j)-AAlMAkA(1:n, 1)*ALk(parti, j)-&
            invAkl(1:n, parti)*AAkMAlALk(1, j)-invAkl(1:n, 1)*AAkMAlALk(parti, j)
      ENDDO

      dktau1 = -AAlMAkALk-AAkMAlALk+AAlMLk+MAlALk

      ddTk = (-TWO/(THREE*A11*A11)*(ml*A1j*Glob_MassMatrix(1, parti)+mk*A1i*Glob_MassMatrix(1, partj))*dkA11+&
          TWO/(THREE*A11)*(ml*Glob_MassMatrix(1, parti)*dkA1j+mk*Glob_MassMatrix(1, partj)*dkA1i))
      ddTk = ddTk-TWO*( ( dkAAlMji+dkAAkMij - &
          mkpml/(THREE*A11*A11)*(A1j*AAlM1i+A1i*AAkM1j)*dkA11+&
          mkpml/(THREE*A11)*(dkA1j*AAlM1i+dkA1i*AAkM1j+A1j*dkAAlM1i+A1i*dkAAkM1j) )  )
      ddTk = ddTk+mkml*( Glob_MassMatrix(1, 1)*&
          ( TWO*(dkAij - (mkpml-TWO)/(THREE*A11*A11)*A1i*A1j*dkA11 +&
          (mkpml-TWO)/(THREE*A11)*(A1i*dkA1j+dkA1i*A1j))/(A11*(mkpml+ONE)) - &
          TWO*(Aij+(mkpml-TWO)/THREE*A1i*A1j/A11)/(A11*A11*(mkpml+ONE))*dkA11 )  )
      ddTk = ddTk-TWO*( ( -dkA11*TWO/(THREE*A11*A11*A11)*&
          (THREE*mk*A11*Aij*AAlM11+THREE*ml*A11*Aij*AAkM11+&
          mk*A11*A1j*AAlMi1+ml*A11*A1j*AAkMi1+&
          mk*A1i*A11*AAlMj1+ml*A1i*A11*AAkMj1+&
          (mkpml-TWO)*mk*A1i*A1j*AAlM11+(mkpml-TWO)*ml*A1i*A1j*AAkM11) + &
          ONE/(THREE*A11*A11)*&
          (THREE*mk*(dkA11*Aij*AAlM11+A11*dkAij*AAlM11+A11*Aij*dkAAlM11)+&
          THREE*ml*(dkA11*Aij*AAkM11+A11*dkAij*AAkM11+A11*Aij*dkAAkM11)+&
          mk*(dkA11*A1j*AAlMi1+A11*dkA1j*AAlMi1+A11*A1j*dkAAlMi1)+&
          ml*(dkA11*A1j*AAkMi1+A11*dkA1j*AAkMi1+A11*A1j*dkAAkMi1)+&
          mk*(dkA1i*A11*AAlMj1+A1i*dkA11*AAlMj1+A1i*A11*dkAAlMj1)+&
          ml*(dkA1i*A11*AAkMj1+A1i*dkA11*AAkMj1+A1i*A11*dkAAkMj1)+&
          (mkpml-TWO)*mk*(dkA1i*A1j*AAlM11+A1i*dkA1j*AAlM11+A1i*A1j*dkAAlM11)+&
          (mkpml-TWO)*ml*(dkA1i*A1j*AAkM11+A1i*dkA1j*AAkM11+A1i*A1j*dkAAkM11)  ) ))
      ddTk = ddTk+FOUR*( -dkA11/(THREE*A11*A11*A11)*(&
          (mkpml-TWO)*mkpml/TWO*A1i*A1j*XM11 + &
          mkpml*(A11*A1j*XM1i+A11*A1i*XM1j+THREE/TWO*A11*XM11*Aij+THREE*A11*A1i*A1j*tau1)+&
          THREE*(A11*A11*XMij+THREE*A11*A11*Aij*tau1)) +&
          ONE/(SIX*A11*A11)*( &
          (mkpml-TWO)*mkpml/TWO*(dkA1i*A1j*XM11+A1i*dkA1j*XM11+A1i*A1j*dkX11) +&
          mkpml*(dkA11*A1j*XM1i+A11*dkA1j*XM1i+A11*A1j*dkX1i+&
          dkA11*A1i*XM1j+A11*dkA1i*XM1j+A11*A1i*dkX1j+&
          THREE/TWO*(dkA11*XM11*Aij+A11*dkX11*Aij+A11*XM11*dkAij)+&
          THREE*(dkA11*A1i*A1j*tau1+A11*dkA1i*A1j*tau1+&
          A11*A1i*dkA1j*tau1+A11*A1i*A1j*dktau1)) +&
          THREE*(TWO*dkA11*A11*XMij+A11*A11*dkXij+&
          THREE*(TWO*dkA11*A11*Aij*tau1+A11*A11*dkAij*tau1+&
          A11*A11*Aij*dktau1)) ))
      ddTk = ddTk*Skl/(Aij+mkpml*A1i*A1j*oA11*ONETHIRD)
      ddTk = ddTk+(ddSk_over_Skl-(dkAij+mkpml/(THREE*A11)*&
          (dkA1i*A1j+A1i*dkA1j-A1i*A1j*dkA11/A11))/(Aij+mkpml/THREE*A1i*A1j/A11))*Tkl

      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dTk(indx) = ddTk(i, j)
        ENDDO
      ENDDO

    ENDIF

    IF (gradflag > 1) THEN

      MAlALl = ZERO
      MAkALl = ZERO
      MLl = ZERO
      AAlMAkALl = ZERO
      AAkMAlALl = ZERO
      AAkMLl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            MAlALl(i, j) = MAlALl(i, j)+MAl(i, k)*ALl(k, j)
            MAkALl(i, j) = MAkALl(i, j)+MAk(i, k)*ALl(k, j)
            MLl(i, j) = MLl(i, j)+Glob_MassMatrix(i, k)*PLl(k, j)
            AAlMAkALl(i, j) = AAlMAkALl(i, j)+AAlMAk(i, k)*ALl(k, j)
            AAkMAlALl(i, j) = AAkMAlALl(i, j)+AAkMAl(i, k)*ALl(k, j)
            AAkMLl(i, j) = AAkMLl(i, j)+AAkM(i, k)*PLl(k, j)
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, n
        dlAAkMij(1:n, j) = -AAkM(1:n, partj)*ALl(parti, j)-invAkl(1:n, parti)*MAkALl(partj, j)
        dlAAkM1j(1:n, j) = -AAkM(1:n, partj)*ALl(1, j)-invAkl(1:n, 1)*MAkALl(partj, j)
        dlAAkM11(1:n, j) = -AAkM(1:n, 1)*ALl(1, j)-invAkl(1:n, 1)*MAkALl(1, j)
        dlAAkMi1(1:n, j) = -AAkM(1:n, 1)*ALl(parti, j)-invAkl(1:n, parti)*MAkALl(1, j)
        dlAAkMj1(1:n, j) = -AAkM(1:n, 1)*ALl(partj, j)-invAkl(1:n, partj)*MAkALl(1, j)

        dlAAlMji(1:n, j) = -AAlM(1:n, parti)*ALl(partj, j)-invAkl(1:n, partj)*MAlALl(parti, j)+&
            invAkl(1:n, partj)*MLl(parti, j)+Glob_MassMatrix(1:n, parti)*ALl(partj, j)
        dlAAlM1i(1:n, j) = -AAlM(1:n, parti)*ALl(1, j)-invAkl(1:n, 1)*MAlALl(parti, j)+&
            invAkl(1:n, 1)*MLl(parti, j)+Glob_MassMatrix(1:n, parti)*ALl(1, j)
        dlAAlM11(1:n, j) = -AAlM(1:n, 1)*ALl(1, j)-invAkl(1:n, 1)*MAlALl(1, j)+&
            invAkl(1:n, 1)*MLl(1, j)+Glob_MassMatrix(1:n, 1)*ALl(1, j)
        dlAAlMi1(1:n, j) = -AAlM(1:n, 1)*ALl(parti, j)-invAkl(1:n, parti)*MAlALl(1, j)+&
            invAkl(1:n, parti)*MLl(1, j)+Glob_MassMatrix(1:n, 1)*ALl(parti, j)
        dlAAlMj1(1:n, j) = -AAlM(1:n, 1)*ALl(partj, j)-invAkl(1:n, partj)*MAlALl(1, j)+&
            invAkl(1:n, partj)*MLl(1, j)+Glob_MassMatrix(1:n, 1)*ALl(partj, j)

        dlXij(1:n, j) = -invAkl(1:n, parti)*AAlMAkALl(partj, j)-invAkl(1:n, partj)*AAlMAkALl(parti, j)-&
            AAkMAlA(1:n, parti)*ALl(partj, j)-AAkMAlA(1:n, partj)*ALl(parti, j)+&
            invAkl(1:n, parti)*AAkMLl(partj, j)+invAkl(1:n, partj)*AAkMLl(parti, j)+&
            AAkM(parti, 1:n)*ALl(partj, j)+AAkM(partj, 1:n)*ALl(parti, j)-&
            AAlMAkA(1:n, parti)*ALl(partj, j)-AAlMAkA(1:n, partj)*ALl(parti, j)-&
            invAkl(1:n, parti)*AAkMAlALl(partj, j)-invAkl(1:n, partj)*AAkMAlALl(parti, j)
        dlX11(1:n, j) = -invAkl(1:n, 1)*AAlMAkALl(1, j)-invAkl(1:n, 1)*AAlMAkALl(1, j)-&
            AAkMAlA(1:n, 1)*ALl(1, j)-AAkMAlA(1:n, 1)*ALl(1, j)+&
            invAkl(1:n, 1)*AAkMLl(1, j)+invAkl(1:n, 1)*AAkMLl(1, j)+&
            AAkM(1, 1:n)*ALl(1, j)+AAkM(1, 1:n)*ALl(1, j)-&
            AAlMAkA(1:n, 1)*ALl(1, j)-AAlMAkA(1:n, 1)*ALl(1, j)-&
            invAkl(1:n, 1)*AAkMAlALl(1, j)-invAkl(1:n, 1)*AAkMAlALl(1, j)
        dlX1j(1:n, j) = -invAkl(1:n, 1)*AAlMAkALl(partj, j)-invAkl(1:n, partj)*AAlMAkALl(1, j)-&
            AAkMAlA(1:n, 1)*ALl(partj, j)-AAkMAlA(1:n, partj)*ALl(1, j)+&
            invAkl(1:n, 1)*AAkMLl(partj, j)+invAkl(1:n, partj)*AAkMLl(1, j)+&
            AAkM(1, 1:n)*ALl(partj, j)+AAkM(partj, 1:n)*ALl(1, j)-&
            AAlMAkA(1:n, 1)*ALl(partj, j)-AAlMAkA(1:n, partj)*ALl(1, j)-&
            invAkl(1:n, 1)*AAkMAlALl(partj, j)-invAkl(1:n, partj)*AAkMAlALl(1, j)
        dlX1i(1:n, j) = -invAkl(1:n, parti)*AAlMAkALl(1, j)-invAkl(1:n, 1)*AAlMAkALl(parti, j)-&
            AAkMAlA(1:n, parti)*ALl(1, j)-AAkMAlA(1:n, 1)*ALl(parti, j)+&
            invAkl(1:n, parti)*AAkMLl(1, j)+invAkl(1:n, 1)*AAkMLl(parti, j)+&
            AAkM(parti, 1:n)*ALl(1, j)+AAkM(1, 1:n)*ALl(parti, j)-&
            AAlMAkA(1:n, parti)*ALl(1, j)-AAlMAkA(1:n, 1)*ALl(parti, j)-&
            invAkl(1:n, parti)*AAkMAlALl(1, j)-invAkl(1:n, 1)*AAkMAlALl(parti, j)
      ENDDO

      dltau1 = -AAlMAkALl-AAkMAlALl+AAkMLl+MAkALl

      ! The kinetic part is assembled with respect to PLl in TemplT, and
      ! Templ keeps d(ln Skl)/dPLl from above: the overlap term of the
      ! kinetic gradient below and the potential gradient (dRl) need it
      ! before the SymMatrix transformation (ddSl_over_Skl is already
      ! transformed)
      TemplT = &
          (-TWO/(THREE*A11*A11)*(ml*A1j*Glob_MassMatrix(1, parti)+mk*A1i*Glob_MassMatrix(1, partj))*dlA11+&
          TWO/(THREE*A11)*(ml*Glob_MassMatrix(1, parti)*dlA1j+mk*Glob_MassMatrix(1, partj)*dlA1i))

      TemplT = TemplT-TWO*( dlAAlMji+dlAAkMij-mkpml/(THREE*A11*A11)*(A1j*AAlM1i+A1i*AAkM1j)*dlA11+&
          mkpml/(THREE*A11)*(dlA1j*AAlM1i+dlA1i*AAkM1j+A1j*dlAAlM1i+A1i*dlAAkM1j) )
      TemplT = TemplT+mkml*Glob_MassMatrix(1, 1)*&
          ( TWO*(dlAij - (mkpml-TWO)/(THREE*A11*A11)*A1i*A1j*dlA11 +&
          (mkpml-TWO)/(THREE*A11)*(A1i*dlA1j+dlA1i*A1j))/(A11*(mkpml+ONE)) - &
          TWO*(Aij+(mkpml-TWO)/THREE*A1i*A1j/A11)/(A11*A11*(mkpml+ONE))*dlA11 )
      TemplT = TemplT-TWO*( -dlA11*TWO/(THREE*A11*A11*A11)*&
          (THREE*mk*A11*Aij*AAlM11+THREE*ml*A11*Aij*AAkM11+&
          mk*A11*A1j*AAlMi1+ml*A11*A1j*AAkMi1+&
          mk*A1i*A11*AAlMj1+ml*A1i*A11*AAkMj1+&
          (mkpml-TWO)*mk*A1i*A1j*AAlM11+(mkpml-TWO)*ml*A1i*A1j*AAkM11) + &
          ONE/(THREE*A11*A11)*&
          (THREE*mk*(dlA11*Aij*AAlM11+A11*dlAij*AAlM11+A11*Aij*dlAAlM11)+&
          THREE*ml*(dlA11*Aij*AAkM11+A11*dlAij*AAkM11+A11*Aij*dlAAkM11)+&
          mk*(dlA11*A1j*AAlMi1+A11*dlA1j*AAlMi1+A11*A1j*dlAAlMi1)+&
          ml*(dlA11*A1j*AAkMi1+A11*dlA1j*AAkMi1+A11*A1j*dlAAkMi1)+&
          mk*(dlA1i*A11*AAlMj1+A1i*dlA11*AAlMj1+A1i*A11*dlAAlMj1)+&
          ml*(dlA1i*A11*AAkMj1+A1i*dlA11*AAkMj1+A1i*A11*dlAAkMj1)+&
          (mkpml-TWO)*mk*(dlA1i*A1j*AAlM11+A1i*dlA1j*AAlM11+A1i*A1j*dlAAlM11)+&
          (mkpml-TWO)*ml*(dlA1i*A1j*AAkM11+A1i*dlA1j*AAkM11+A1i*A1j*dlAAkM11)  ))
      TemplT = TemplT+FOUR*( -dlA11/(THREE*A11*A11*A11)*(&
          (mkpml-TWO)*mkpml/TWO*A1i*A1j*XM11 + &
          mkpml*(A11*A1j*XM1i+A11*A1i*XM1j+THREE/TWO*A11*XM11*Aij+THREE*A11*A1i*A1j*tau1)+&
          THREE*(A11*A11*XMij+THREE*A11*A11*Aij*tau1)) +&
          ONE/(SIX*A11*A11)*( &
          (mkpml-TWO)*mkpml/TWO*(dlA1i*A1j*XM11+A1i*dlA1j*XM11+A1i*A1j*dlX11) +&
          mkpml*(dlA11*A1j*XM1i+A11*dlA1j*XM1i+A11*A1j*dlX1i+&
          dlA11*A1i*XM1j+A11*dlA1i*XM1j+A11*A1i*dlX1j+&
          THREE/TWO*(dlA11*XM11*Aij+A11*dlX11*Aij+A11*XM11*dlAij)+&
          THREE*(dlA11*A1i*A1j*tau1+A11*dlA1i*A1j*tau1+&
          A11*A1i*dlA1j*tau1+A11*A1i*A1j*dltau1)) +&
          THREE*(TWO*dlA11*A11*XMij+A11*A11*dlXij+&
          THREE*(TWO*dlA11*A11*Aij*tau1+A11*A11*dlAij*tau1+&
          A11*A11*Aij*dltau1))  ))

      TemplT = TemplT*Skl/(Aij+mkpml*A1i*A1j/(THREE*A11))
      TemplT = TemplT+(Templ-(dlAij+mkpml/(THREE*A11)*(dlA1i*A1j+&
          A1i*dlA1j-A1i*A1j*dlA11/A11))/(Aij+mkpml/THREE*A1i*A1j/A11))*Tkl

      ddTl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ddTl(i, j) = ddTl(i, j)+SymMatrix(i, k)*TemplT(k, j)
          ENDDO
        ENDDO
      ENDDO
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dTl(indx) = ddTl(i, j)
        ENDDO
      ENDDO


    ENDIF

    !********************************************************************
    ! Build intermediate results for potential energy  matrix element
    ! a = invAkl(1,1)
    ! b(p,p)=invAkl(p,p)
    ! b(p,q)=invAkl(p,p)+invAkl(q,q)-2*invAkl(p,q)
    ! c(p,p)=invAkl(p,p)**2
    ! c(p,q)=(invAkl(1,p)-invAkl(1,q))**2
    ! d = invAkl(parti,partj)
    ! f(p,p)=invAkl(parti,p)*invAkl(partj,p)
    ! f(p,q)=(invAkl(parti,p)-invAkl(parti,q))*(invAkl(partj,p)-invAkl(partj,q))
    ! h = invAkl(1,parti)*invAkl(1,partj)
    ! coab(p,q)=ONE-c(p,q)/(a*b(p,q))
    !********************************************************************

    ! b
    b = ONE
    DO p = 1, n
      b(p, p) = invAkl(p, p)
    ENDDO
    DO q = 1, n-1
      DO p = q+1, n
        b(p, q) = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
      ENDDO
    ENDDO

    ! coab
    coab = ONE
    DO i = 1, n
      coab(i, i) = ONE-invAkl(1, i)*invAkl(1, i)/(invAkl(1, 1)*invAkl(i, i))
    ENDDO
    DO j = 1, n-1
      DO i = j+1, n
        coab(i, j) = ONE-(invAkl(i, 1)-invAkl(j, 1))*(invAkl(i, 1)-invAkl(j, 1))/&
            invAkl(1, 1)/b(i, j)
      ENDDO
    ENDDO

    ! f/b
    fob = ONE
    DO p = 1, n
      fob(p, p) = invAkl(parti, p)*invAkl(partj, p)/b(p, p)
    ENDDO
    DO q = 1, n-1
      DO p = q+1, n
        fob(p, q) = (invAkl(parti, p)-invAkl(parti, q))*&
            (invAkl(partj, p)-invAkl(partj, q))/b(p, q)
      ENDDO
    ENDDO

    ! h/a
    hoa = invAkl(1, parti)*invAkl(1, partj)/invAkl(1, 1)

    !**************************************************
    ! Potential energy EVEN POWERS ONLY!!!!!!
    !**************************************************

    ! We will compute all of the 1/rij terms at once
    ! and place them in R, then sum over R*ChargeMatrix to
    ! get Vkl

    m = (mk+ml)/2
    R = ZERO

    R(1, 1) = Skl/(SQRT(b(1, 1)))*Glob_gm2(m)*(ONE-fob(1, 1)/&
        (THREE*invAkl(parti, partj)+TWO*m*hoa))

    ! All pairs of eq. (34) of Kirnosov, Sharkey, Adamowicz, J. Phys. B 48,
    ! 195101 (2015): R(k, j), j < k, is the 1/r_jk term and R(k, k) the
    ! 1/r_k term; R(1, 1) is computed above.
    ! From 6 pseudoparticles on, in double precision (PairedPowerSums),
    ! CoulombPowerSums forms the same sums for all pairs together (n and
    ! wp are constants: only one branch is compiled).
    IF ((n >= 6) .AND. PairedPowerSums) THEN

      CALL CoulombPowerSums(m, coab, fob, invAkl(parti, partj), hoa, R)

    ELSE

      DO k = 2, n
        DO j = 1, k
          g1 = coab(k, j)
          g2 = ONE
          g3 = ZERO
          DO i = 1, m
            g2 = g2*g1
            g3 = g3+Glob_gm3(i)*g2*((TWO*i+THREE)*invAkl(parti, partj)-&
                (TWO*i+ONE)*fob(k, j)+TWO*(m-i)*hoa)
          ENDDO
          R(k, j) = g3
        ENDDO
      ENDDO

    ENDIF

    g4 = Skl*Glob_gm2(m)/(THREE*invAkl(parti, partj)+TWO*m*hoa)
    DO k = 2, n
      DO j = 1, k
        R(k, j) = (R(k, j)+THREE*invAkl(parti, partj)-fob(k, j)+TWO*m*hoa&
            )*g4/SQRT(b(k, j))
      ENDDO
    ENDDO

    ! Now compute the potential energy
    Vkl = ZERO
    DO j = 1, n
      DO i = j, n
        Vkl = Vkl+R(i, j)*Glob_ChargeMatrix(i, j)
      ENDDO
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

      ! dhoa is independent of i and j
      dkhoa = hoa*(dkA1i/A1i+dkA1j/A1j-dkA11/A11)

      ! diagonal 1/rij terms
      ddVk = ZERO

      ! p=q=1
      dRk = R(1, 1)*(ddSk_over_Skl-dkA11/(TWO*invAkl(1, 1))-&
          (THREE*dkAij+TWO*m*dkhoa)/(THREE*Aij+TWO*m*hoa)+&
          (THREE*dkAij+(TWO*m-ONE)*dkhoa)/(THREE*Aij+(TWO*m-ONE)*hoa))
      ddVk = ddVk+dRk*Glob_ChargeMatrix(1, 1)
      ! p=q != 1
      DO p = 2, n
        dkfob = fob(p, p)*(dkA(parti, p, :, :)/invAkl(parti, p)+dkA(partj, p, :, :)/&
            invAkl(partj, p)-dkA(p, p, :, :)/invAkl(p, p))
        dkb = dkA(p, p, :, :)
        dkcoab = (invAkl(1, p)*invAkl(1, p)/invAkl(1, 1)/invAkl(p, p))*&
            (dkA(1, 1, :, :)/invAkl(1, 1)+dkA(p, p, :, :)/invAkl(p, p)-&
            TWO*dkA(1, p, :, :)/invAkl(1, p))

        q = p
        g2 = ONE
        g5 = ZERO
        DO k = 1, (m-1)
          g2 = g2*coab(p, q)
          g5 = g5+g2*Glob_gm3(k)*(k/coab(p, q)*((TWO*k+THREE)*Aij-(TWO*k+ONE)*fob(p, q)+&
              TWO*(m-k)*hoa)*dkcoab+(TWO*k+THREE)*dkAij-(TWO*k+ONE)*dkfob+TWO*&
              (m-k)*dkhoa)
        ENDDO
        dRk = R(p, q)*(ddSk_over_Skl-dkb/(TWO*b(p, q))-&
            (THREE*dkAij+TWO*m*dkhoa)/(THREE*Aij+TWO*m*hoa))+&
            Skl/(SQRT(b(p, q))*(THREE*Aij+TWO*m*hoa))*Glob_gm2(m)*(THREE*dkAij-&
            dkfob+TWO*m*dkhoa+Glob_gm3(m)*g2*coab(p, q)*(m/coab(p, q)*((TWO*&
            m+THREE)*Aij-(TWO*m+ONE)*fob(p, q))*dkcoab+(TWO*m+THREE)*&
            dkAij-(TWO*m+ONE)*dkfob)+g5)
        ddVk = ddVk+dRk*Glob_ChargeMatrix(p, p)
      ENDDO

      ! Now do the off diagonal terms
      DO q = 1, n-1
        DO p = q+1, n

          dkb = dkA(p, p, :, :)+dkA(q, q, :, :)-TWO*dkA(p, q, :, :)

          ! product forms: finite also where a factor of fob or coab is zero (for example k = l with a
          ! symmetry term that swaps p and q, invAkl(1,p) = invAkl(1,q))
          dkfob = ((dkA(parti, p, :, :)-dkA(parti, q, :, :))*(invAkl(partj, p)-invAkl(partj, q))+&
              (invAkl(parti, p)-invAkl(parti, q))*(dkA(partj, p, :, :)-dkA(partj, q, :, :))-&
              fob(p, q)*dkb)/b(p, q)

          dkcoab = ((invAkl(1, p)-invAkl(1, q))*(invAkl(1, p)-invAkl(1, q))*&
              (dkA(1, 1, :, :)/invAkl(1, 1)+dkb/b(p, q))-&
              TWO*(invAkl(1, p)-invAkl(1, q))*(dkA(1, p, :, :)-dkA(1, q, :, :)))/&
              (invAkl(1, 1)*b(p, q))

          g2 = ONE
          g5 = ZERO
          DO k = 1, (m-1)
            g2 = g2*coab(p, q)
            g5 = g5+g2*Glob_gm3(k)*(k/coab(p, q)*((TWO*k+THREE)*Aij-(TWO*k+ONE)*fob(p, q)+&
                TWO*(m-k)*hoa)*dkcoab+(TWO*k+THREE)*dkAij-(TWO*k+ONE)*dkfob+TWO*&
                (m-k)*dkhoa)
          ENDDO
          dRk = R(p, q)*(ddSk_over_Skl-dkb/(TWO*b(p, q))-&
              (THREE*dkAij+TWO*m*dkhoa)/(THREE*Aij+TWO*m*hoa))+&
              Skl/(SQRT(b(p, q))*(THREE*Aij+TWO*m*hoa))*Glob_gm2(m)*(THREE*dkAij-&
              dkfob+TWO*m*dkhoa+Glob_gm3(m)*g2*coab(p, q)*(m/coab(p, q)*((TWO*&
              m+THREE)*Aij-(TWO*m+ONE)*fob(p, q))*dkcoab+(TWO*m+THREE)*&
              dkAij-(TWO*m+ONE)*dkfob)+g5)
          ddVk = ddVk+dRk*Glob_ChargeMatrix(p, q)

        ENDDO
      ENDDO

      ! Load it into dVk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dVk(indx) = ddVk(i, j)
        ENDDO
      ENDDO

    ENDIF

    IF (gradflag > 1) THEN

      !********
      ! dVl
      !********

      ! dhoa is independent of i and j
      dlhoa = hoa*(dlA1i/A1i+dlA1j/A1j-dlA11/A11)

      ! diagonal 1/rij terms
      ddVl = ZERO

      ! p=q=1
      dRl = R(1, 1)*(Templ-dlA11/(TWO*invAkl(1, 1))-&
          (THREE*dlAij+TWO*m*dlhoa)/(THREE*Aij+TWO*m*hoa)+&
          (THREE*dlAij+(TWO*m-ONE)*dlhoa)/(THREE*Aij+(TWO*m-ONE)*hoa))
      ddVl = ddVl+dRl*Glob_ChargeMatrix(1, 1)

      ! p=q != 1
      DO p = 2, n
        dlfob = fob(p, p)*(dlA(parti, p, :, :)/invAkl(parti, p)+dlA(partj, p, :, :)/&
            invAkl(partj, p)-dlA(p, p, :, :)/invAkl(p, p))
        dlb = dlA(p, p, :, :)
        dlcoab = (invAkl(1, p)*invAkl(1, p)/invAkl(1, 1)/invAkl(p, p))*&
            (dlA(1, 1, :, :)/invAkl(1, 1)+dlA(p, p, :, :)/invAkl(p, p)-&
            TWO*dlA(1, p, :, :)/invAkl(1, p))

        q = p
        g2 = ONE
        g5 = ZERO
        DO k = 1, (m-1)
          g2 = g2*coab(p, q)
          g5 = g5+g2*Glob_gm3(k)*(k/coab(p, q)*((TWO*k+THREE)*Aij-(TWO*k+ONE)*fob(p, q)+&
              TWO*(m-k)*hoa)*dlcoab+(TWO*k+THREE)*dlAij-(TWO*k+ONE)*dlfob+TWO*&
              (m-k)*dlhoa)
        ENDDO
        dRl = R(p, q)*(Templ-dlb/(TWO*b(p, q))-&
            (THREE*dlAij+TWO*m*dlhoa)/(THREE*Aij+TWO*m*hoa))+&
            Skl/(SQRT(b(p, q))*(THREE*Aij+TWO*m*hoa))*Glob_gm2(m)*(THREE*dlAij-&
            dlfob+TWO*m*dlhoa+Glob_gm3(m)*g2*coab(p, q)*(m/coab(p, q)*((TWO*&
            m+THREE)*Aij-(TWO*m+ONE)*fob(p, q))*dlcoab+(TWO*m+THREE)*&
            dlAij-(TWO*m+ONE)*dlfob)+g5)
        ddVl = ddVl+dRl*Glob_ChargeMatrix(p, p)
      ENDDO

      ! Now do the off diagonal terms
      DO q = 1, n-1
        DO p = q+1, n

          dlb = dlA(p, p, :, :)+dlA(q, q, :, :)-TWO*dlA(p, q, :, :)

          ! product forms: finite also where a factor of fob or coab is zero (for example k = l with a
          ! symmetry term that swaps p and q, invAkl(1,p) = invAkl(1,q))
          dlfob = ((dlA(parti, p, :, :)-dlA(parti, q, :, :))*(invAkl(partj, p)-invAkl(partj, q))+&
              (invAkl(parti, p)-invAkl(parti, q))*(dlA(partj, p, :, :)-dlA(partj, q, :, :))-&
              fob(p, q)*dlb)/b(p, q)

          dlcoab = ((invAkl(1, p)-invAkl(1, q))*(invAkl(1, p)-invAkl(1, q))*&
              (dlA(1, 1, :, :)/invAkl(1, 1)+dlb/b(p, q))-&
              TWO*(invAkl(1, p)-invAkl(1, q))*(dlA(1, p, :, :)-dlA(1, q, :, :)))/&
              (invAkl(1, 1)*b(p, q))

          g2 = ONE
          g5 = ZERO
          DO k = 1, (m-1)
            g2 = g2*coab(p, q)
            g5 = g5+g2*Glob_gm3(k)*(k/coab(p, q)*((TWO*k+THREE)*Aij-(TWO*k+ONE)*fob(p, q)+&
                TWO*(m-k)*hoa)*dlcoab+(TWO*k+THREE)*dlAij-(TWO*k+ONE)*dlfob+TWO*&
                (m-k)*dlhoa)
          ENDDO
          dRl = R(p, q)*(Templ-dlb/(TWO*b(p, q))-&
              (THREE*dlAij+TWO*m*dlhoa)/(THREE*Aij+TWO*m*hoa))+&
              Skl/(SQRT(b(p, q))*(THREE*Aij+TWO*m*hoa))*Glob_gm2(m)*(THREE*dlAij-&
              dlfob+TWO*m*dlhoa+Glob_gm3(m)*g2*coab(p, q)*(m/coab(p, q)*((TWO*&
              m+THREE)*Aij-(TWO*m+ONE)*fob(p, q))*dlcoab+(TWO*m+THREE)*&
              dlAij-(TWO*m+ONE)*dlfob)+g5)
          ddVl = ddVl+dRl*Glob_ChargeMatrix(p, q)

        ENDDO
      ENDDO

      Templ = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            Templ(i, j) = Templ(i, j)+SymMatrix(i, k)*ddVl(k, j)
          ENDDO
        ENDDO
      ENDDO

      !!Load it into dVl (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dVl(indx) = Templ(i, j)
        ENDDO
      ENDDO

    ENDIF


    Hkl = Tkl+Vkl
    SELECT CASE (gradflag)
    CASE (0)
      dHk = ZERO
      dHl = ZERO
      dSk = ZERO
      dSl = ZERO
    CASE (1)
      dHk = dTk+dVk
      dHl = ZERO
      dSl = ZERO
    CASE (2)
      dHk = dTk+dVk
      dHl = dTl+dVl
    END SELECT


  END SUBROUTINE MatrixElementsXi





  SUBROUTINE CoulombPowerSums(m, coab, fob, Aij, hoa, R)
    ! The power sums of the Coulomb integrals of MatrixElementsXi for all
    ! pairs (k, j) of eq. (34), 2 <= k <= n, j <= k:
    !   R(k, j) = sum_{i=1}^{m} gm3(i) coab(k,j)^i ((2i+3) Aij
    !             - (2i+1) fob(k,j) + 2(m-i) hoa),
    ! used from 6 pseudoparticles on. The sums of the pairs are independent,
    ! so they are advanced together, one element of the Pair* vectors per
    ! pair: every element gets exactly the operations of the sum done alone
    ! in MatrixElementsXi, but the chains no longer wait for each other.
    ! Aij and hoa are copied into every element, so that the products with
    ! them stay inside the expression as in the sum done alone. With fewer
    ! pseudoparticles the compiler sees through the copies, or compiles
    ! MatrixElementsXi differently, so there the sums stay in the kernel.
    ! Input:
    !   m               :: (mk+ml)/2
    !   coab, fob       :: as in MatrixElementsXi
    !   Aij             :: invAkl(parti, partj)
    !   hoa             :: as in MatrixElementsXi
    ! Input/Output:
    !   R               :: R(k, j), k >= 2, set to the power sums
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles
    INTEGER, PARAMETER :: NumPairs = n*(n+1)/2-1  ! the pairs (k, j), 2 <= k <= n, j <= k

    INTEGER, INTENT(IN) :: m
    REAL(wp), INTENT(IN) :: coab(n, n), fob(n, n), Aij, hoa
    REAL(wp), INTENT(INOUT) :: R(n, n)

    ! Local variables
    INTEGER :: i, j, k, t
    REAL(wp) :: PairG1(NumPairs)   ! coab of each pair
    REAL(wp) :: PairG2(NumPairs)   ! coab**i of each pair
    REAL(wp) :: PairG3(NumPairs)   ! power sum of each pair
    REAL(wp) :: PairFob(NumPairs)  ! fob of each pair
    REAL(wp) :: PairAij(NumPairs)  ! Aij, the same for all pairs
    REAL(wp) :: PairHoa(NumPairs)  ! hoa, the same for all pairs


    t = 0
    DO k = 2, n
      DO j = 1, k
        t = t+1
        PairG1(t) = coab(k, j)
        PairFob(t) = fob(k, j)
        PairAij(t) = Aij
        PairHoa(t) = hoa
        PairG2(t) = ONE
        PairG3(t) = ZERO
      ENDDO
    ENDDO

    DO i = 1, m
      DO t = 1, NumPairs
        PairG2(t) = PairG2(t)*PairG1(t)
        PairG3(t) = PairG3(t)+Glob_gm3(i)*PairG2(t)*((TWO*i+THREE)*PairAij(t)-&
            (TWO*i+ONE)*PairFob(t)+TWO*(m-i)*PairHoa(t))
      ENDDO
    ENDDO

    t = 0
    DO k = 2, n
      DO j = 1, k
        t = t+1
        R(k, j) = PairG3(t)
      ENDDO
    ENDDO

  END SUBROUTINE CoulombPowerSums





  SUBROUTINE GramTermPair(invAkl, parti, partj, m, p, q, IsOn, Gn, Q0, Q1)
    ! The Gram-determinant term of the Coulomb integral R(p, q) of
    ! MatrixElementsXi for the bra prefactor x_i (i = parti), the ket term
    ! x_j (j = partj) and m = (mk+ml)/2. With a = invAkl(1,1),
    ! Aij = invAkl(i,j), A1i = invAkl(1,i), v = e_p (q = p) or e_p - e_q,
    !   beta = v'*invAkl*v,  c = (invAkl*v)(1),  d_i = (invAkl*v)(i),
    !   coab = 1 - c^2/(a*beta), hoa = A1i*A1j/a, fob = d_i*d_j/beta,
    !   e    = c*(A1i*d_j + d_i*A1j)/(a*beta),
    !   Gn   = coab*Aij - hoa - fob + e = G_ij/(a*beta),
    !   Q0   = sum_{t=0}^{m-1} (2t+1) gm3(t) coab^t,  Q1 = dQ0/dcoab,
    ! G_ij being the 3x3 Gram determinant of (e1, v, e_i) and (e1, v, e_j)
    ! in the metric invAkl. The bracket of eq. (35) of the paper gets
    ! -Gn*Q0. The term is absent (IsOn false) when m = 0 or e_i or e_j
    ! lies in span{e1, v}: e_i lies in it when i = 1, or i = p and v is
    ! e_p or e_p - e_1 (q = p or q = 1).
    ! Input:
    !   invAkl          :: inverse of Akl
    !   parti, partj    :: x-indices of the bra and the ket term
    !   m               :: (mk+ml)/2
    !   p, q            :: the pair (q <= p)
    ! Output:
    !   IsOn            :: the term is present
    !   Gn, Q0, Q1      :: as above (if IsOn)
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    REAL(wp), INTENT(IN) :: invAkl(n, n)
    INTEGER, INTENT(IN) :: parti, partj, m, p, q
    LOGICAL, INTENT(OUT) :: IsOn
    REAL(wp), INTENT(OUT) :: Gn, Q0, Q1

    ! Local variables
    INTEGER :: t
    LOGICAL :: IsSpanI, IsSpanJ
    REAL(wp) :: a, Aij, A1i, A1j, beta, c, di, dj, coab, w


    IsSpanI = (parti == 1) .OR. ((parti == p) .AND. ((q == p) .OR. (q == 1)))
    IsSpanJ = (partj == 1) .OR. ((partj == p) .AND. ((q == p) .OR. (q == 1)))
    IsOn = (m > 0) .AND. (.NOT. IsSpanI) .AND. (.NOT. IsSpanJ)
    IF (.NOT. IsOn) RETURN

    a = invAkl(1, 1)
    Aij = invAkl(parti, partj)
    A1i = invAkl(1, parti)
    A1j = invAkl(1, partj)

    IF (p == q) THEN
      beta = invAkl(p, p)
      c = invAkl(1, p)
      di = invAkl(parti, p)
      dj = invAkl(partj, p)
    ELSE
      beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
      c = invAkl(1, p)-invAkl(1, q)
      di = invAkl(parti, p)-invAkl(parti, q)
      dj = invAkl(partj, p)-invAkl(partj, q)
    ENDIF

    coab = ONE-c*c/(a*beta)
    Gn = coab*Aij-A1i*A1j/a-di*dj/beta+c*(A1i*dj+di*A1j)/(a*beta)

    ! Q0 and Q1 = dQ0/dcoab; w = coab^(t-1) at the top of the loop
    w = ONE
    Q0 = ONE
    Q1 = ZERO
    DO t = 1, m-1
      Q1 = Q1+(TWO*t+ONE)*Glob_gm3(t)*t*w
      w = w*coab
      Q0 = Q0+(TWO*t+ONE)*Glob_gm3(t)*w
    ENDDO

  END SUBROUTINE GramTermPair





  SUBROUTINE GramTermSetup(invAkl, parti, partj, m, p, q, IsOn, Gn, coab)
    ! The Gram-determinant term of the Coulomb integral R(p, q) of
    ! MatrixElementsXi for the bra prefactor x_i (i = parti), the ket term
    ! x_j (j = partj) and m = (mk+ml)/2. With a = invAkl(1,1),
    ! Aij = invAkl(i,j), A1i = invAkl(1,i), v = e_p (q = p) or e_p - e_q,
    !   beta = v'*invAkl*v,  c = (invAkl*v)(1),  d_i = (invAkl*v)(i),
    !   coab = 1 - c^2/(a*beta), hoa = A1i*A1j/a, fob = d_i*d_j/beta,
    !   e    = c*(A1i*d_j + d_i*A1j)/(a*beta),
    !   Gn   = coab*Aij - hoa - fob + e = G_ij/(a*beta),
    ! G_ij being the 3x3 Gram determinant of (e1, v, e_i) and (e1, v, e_j)
    ! in the metric invAkl. The bracket of eq. (35) of the paper gets
    ! -Gn*Q0 with the power sum
    !   Q0   = sum_{t=0}^{m-1} (2t+1) gm3(t) coab^t,  Q1 = dQ0/dcoab,
    ! which GramTermCorrPairs forms for all pairs together. The term is absent
    ! (IsOn false) when m = 0 or e_i or e_j lies in span{e1, v}: e_i lies
    ! in it when i = 1, or i = p and v is e_p or e_p - e_1 (q = p or q = 1).
    ! Input:
    !   invAkl          :: inverse of Akl
    !   parti, partj    :: x-indices of the bra and the ket term
    !   m               :: (mk+ml)/2
    !   p, q            :: the pair (q <= p)
    ! Output:
    !   IsOn            :: the term is present
    !   Gn, coab        :: as above (if IsOn)
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    REAL(wp), INTENT(IN) :: invAkl(n, n)
    INTEGER, INTENT(IN) :: parti, partj, m, p, q
    LOGICAL, INTENT(OUT) :: IsOn
    REAL(wp), INTENT(OUT) :: Gn, coab

    ! Local variables
    LOGICAL :: IsSpanI, IsSpanJ
    REAL(wp) :: a, Aij, A1i, A1j, beta, c, di, dj


    IsSpanI = (parti == 1) .OR. ((parti == p) .AND. ((q == p) .OR. (q == 1)))
    IsSpanJ = (partj == 1) .OR. ((partj == p) .AND. ((q == p) .OR. (q == 1)))
    IsOn = (m > 0) .AND. (.NOT. IsSpanI) .AND. (.NOT. IsSpanJ)
    IF (.NOT. IsOn) RETURN

    a = invAkl(1, 1)
    Aij = invAkl(parti, partj)
    A1i = invAkl(1, parti)
    A1j = invAkl(1, partj)

    IF (p == q) THEN
      beta = invAkl(p, p)
      c = invAkl(1, p)
      di = invAkl(parti, p)
      dj = invAkl(partj, p)
    ELSE
      beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
      c = invAkl(1, p)-invAkl(1, q)
      di = invAkl(parti, p)-invAkl(parti, q)
      dj = invAkl(partj, p)-invAkl(partj, q)
    ENDIF

    coab = ONE-c*c/(a*beta)
    Gn = coab*Aij-A1i*A1j/a-di*dj/beta+c*(A1i*dj+di*A1j)/(a*beta)

  END SUBROUTINE GramTermSetup





  SUBROUTINE CoulombCorrSetup(vechLk, vechLl, SymMatrix, Lk, PLl, invAkl)
    ! Lk, PLl = SymMatrix'*Ll and invAkl = inv(Lk*Lk' + PLl*PLl') as in
    ! MatrixElementsXi (Cholesky factor of Akl), for the corrections
    ! GramTermCorr and PowerZeroGradCorr, which recompute them so that
    ! MatrixElementsXi stays as it was.
    ! Input:
    !   vechLk, vechLl  :: n(n+1)/2 exponent parameters
    !   SymMatrix       :: the symmetry permutation matrix
    ! Output:
    !   Lk, PLl, invAkl
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: Lk(n, n), PLl(n, n), invAkl(n, n)

    ! Local variables
    INTEGER :: i, j, k, indx
    REAL(wp) :: SUM0
    REAL(wp) :: Ll(n, n), invL(n, n)


    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
        Ll(i, j) = vechLl(indx)
      ENDDO
    ENDDO

    PLl = ZERO
    DO j = 1, n
      DO k = j, n
        DO i = 1, n
          PLl(i, j) = PLl(i, j)+SymMatrix(k, i)*Ll(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! Akl in invL, then its Cholesky factor
    invL = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          invL(i, j) = invL(i, j)+Lk(i, k)*Lk(j, k)+PLl(i, k)*PLl(j, k)
        ENDDO
      ENDDO
    ENDDO

    DO i = 1, n
      DO j = i, n
        SUM0 = invL(i, j)
        DO k = i-1, 1, -1
          SUM0 = SUM0-invL(i, k)*invL(j, k)
        ENDDO
        IF (i == j) THEN
          invL(i, i) = SQRT(SUM0)
        ELSE
          invL(j, i) = SUM0/invL(i, i)
        ENDIF
      ENDDO
    ENDDO

    DO i = 1, n
      invL(i, i) = ONE/invL(i, i)
      DO j = i+1, n
        SUM0 = ZERO
        DO k = i, j-1
          SUM0 = SUM0-invL(j, k)*invL(k, i)
        ENDDO
        invL(j, i) = SUM0/invL(j, j)
        invL(i, j) = ZERO
      ENDDO
    ENDDO

    invAkl = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          invAkl(i, j) = invAkl(i, j)+invL(k, i)*invL(k, j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CoulombCorrSetup





  SUBROUTINE CoulombCorrDerivs(invAkl, Lx, dA)
    ! dA(p, q, :, :) = d invAkl(p, q)/d Lx for Lx = Lk or PLl, as dkA and
    ! dlA of MatrixElementsXi: with ALx = invAkl*Lx,
    !   dA(p, q, a, b) = -invAkl(a, q)*ALx(p, b) - invAkl(a, p)*ALx(q, b).
    ! Input:
    !   invAkl          :: inverse of Akl
    !   Lx              :: Lk or PLl
    ! Output:
    !   dA
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    REAL(wp), INTENT(IN) :: invAkl(n, n), Lx(n, n)
    REAL(wp), INTENT(OUT) :: dA(n, n, n, n)

    ! Local variables
    INTEGER :: i, j, k, p, q
    REAL(wp) :: ALx(n, n)


    ALx = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          ALx(i, j) = ALx(i, j)+invAkl(i, k)*Lx(k, j)
        ENDDO
      ENDDO
    ENDDO

    DO p = 1, n
      DO q = 1, n
        DO j = 1, n
          dA(p, q, 1:n, j) = -invAkl(1:n, q)*ALx(p, j)-invAkl(1:n, p)*ALx(q, j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CoulombCorrDerivs





  SUBROUTINE GramTermCorr(mk, ik, vechLk, ml, jl, vechLl, SymMatrix, &
                          Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! Adds to the output of MatrixElementsXi (one ket term x_jl) the
    ! Gram-determinant term of its Coulomb integrals (see GramTermPair):
    !   R(p,q) += Cpq*Skl,  Cpq = -Gn*Q0*gm2(m)/((3 Aij + 2 m hoa) sqrt(beta)),
    ! for every pair where the term is present, so Hkl gets
    ! Skl * sum charge(p,q)*Cpq, and the gradients get
    ! sum charge(p,q)*(Cpq*dS + Skl*dCpq), with
    !   dCpq  = Cpq*(-d(3 Aij + 2 m hoa)/(3 Aij + 2 m hoa) - dbeta/(2 beta))
    !           + gm2(m)/((3 Aij + 2 m hoa) sqrt(beta)) * d(-Gn*Q0),
    !   dG    = sum of the cofactors of the 3x3 Gram matrix times the
    !           derivatives of its elements,
    !   dGn   = dG/(a*beta) - Gn*(da/a + dbeta/beta),
    !   dcoab = -2 c dc/(a*beta) + (1 - coab)*(da/a + dbeta/beta),
    !   d(-Gn*Q0) = -(dGn*Q0 + Gn*Q1*dcoab).
    ! The derivatives of invAkl are those of MatrixElementsXi: with respect
    ! to Lk for dHk, with respect to PLl = SymMatrix'*Ll (then times
    ! SymMatrix) for dHl. Nothing is divided by c, d_i or A1i, which may
    ! vanish. invAkl is recomputed (CoulombCorrSetup) so that
    ! MatrixElementsXi, and its results where the term is absent, stay
    ! exactly as they were; this routine is called only when the term can
    ! be present (mk+ml > 0, ik /= 1, jl /= 1, n > 2).
    ! Input:
    !   mk, ik, vechLk  :: the bra
    !   ml, jl, vechLl  :: the ket term (power, x-index in the transformed
    !                      ket, exponent parameters)
    !   SymMatrix       :: the symmetry permutation matrix
    !   Skl, dSk, dSl   :: overlap and its gradients from MatrixElementsXi
    !   gradflag        :: as in MatrixElementsXi
    ! Input/Output:
    !   Hkl, dHk, dHl   :: Hamiltonian and its gradients
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, jl, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(IN) :: Skl, dSk(Glob_npt), dSl(Glob_npt)
    REAL(wp), INTENT(INOUT) :: Hkl, dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: i, j, k, p, q, m, side, parti, partj, indx
    LOGICAL :: IsOn
    REAL(wp) :: Gn, Q0, Q1, Cpq, chpq, g4den
    REAL(wp) :: a, Aij, A1i, A1j, beta, c, di, dj, kappa
    REAL(wp) :: Lk(n, n), PLl(n, n), invAkl(n, n), dA(n, n, n, n)
    REAL(wp) :: dbeta(n, n), dc(n, n), ddi(n, n), ddj(n, n), dhoa(n, n)
    REAL(wp) :: dG(n, n), dGn(n, n), dcoab(n, n), dCpq(n, n), dCsum(n, n), dCl(n, n)


    m = (mk+ml)/2
    parti = ik
    partj = jl

    CALL CoulombCorrSetup(vechLk, vechLl, SymMatrix, Lk, PLl, invAkl)

    a = invAkl(1, 1)
    Aij = invAkl(parti, partj)
    A1i = invAkl(1, parti)
    A1j = invAkl(1, partj)
    g4den = THREE*Aij+TWO*m*A1i*A1j/a

    ! The energy
    DO p = 2, n
      DO q = 1, p
        CALL GramTermPair(invAkl, parti, partj, m, p, q, IsOn, Gn, Q0, Q1)
        IF (.NOT. IsOn) CYCLE
        IF (p == q) THEN
          beta = invAkl(p, p)
        ELSE
          beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
        ENDIF
        Hkl = Hkl-Glob_ChargeMatrix(p, q)*Gn*Q0*Glob_gm2(m)/(g4den*SQRT(beta))*Skl
      ENDDO
    ENDDO

    IF (gradflag == 0) RETURN

    ! The gradients: side 1 with respect to Lk, side 2 with respect to PLl
    DO side = 1, MIN(gradflag, 2)

      IF (side == 1) THEN
        CALL CoulombCorrDerivs(invAkl, Lk, dA)
      ELSE
        CALL CoulombCorrDerivs(invAkl, PLl, dA)
      ENDIF

      dhoa = (dA(1, parti, :, :)*A1j+A1i*dA(1, partj, :, :))/a-A1i*A1j*dA(1, 1, :, :)/(a*a)

      ! sum over the pairs of charge*Skl*dCpq in dCsum and of charge*Cpq in
      ! chpq (its product with dS is added to the vector below)
      dCsum = ZERO
      chpq = ZERO
      DO p = 2, n
        DO q = 1, p
          CALL GramTermPair(invAkl, parti, partj, m, p, q, IsOn, Gn, Q0, Q1)
          IF (.NOT. IsOn) CYCLE

          IF (p == q) THEN
            beta = invAkl(p, p)
            c = invAkl(1, p)
            di = invAkl(parti, p)
            dj = invAkl(partj, p)
            dbeta = dA(p, p, :, :)
            dc = dA(1, p, :, :)
            ddi = dA(parti, p, :, :)
            ddj = dA(partj, p, :, :)
          ELSE
            beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
            c = invAkl(1, p)-invAkl(1, q)
            di = invAkl(parti, p)-invAkl(parti, q)
            dj = invAkl(partj, p)-invAkl(partj, q)
            dbeta = dA(p, p, :, :)+dA(q, q, :, :)-TWO*dA(p, q, :, :)
            dc = dA(1, p, :, :)-dA(1, q, :, :)
            ddi = dA(parti, p, :, :)-dA(parti, q, :, :)
            ddj = dA(partj, p, :, :)-dA(partj, q, :, :)
          ENDIF

          kappa = c*c/(a*beta)
          Cpq = -Gn*Q0*Glob_gm2(m)/(g4den*SQRT(beta))

          dG = (beta*Aij-di*dj)*dA(1, 1, :, :)+(a*Aij-A1i*A1j)*dbeta+&
              (a*beta-c*c)*dA(parti, partj, :, :)+(A1i*dj+di*A1j-TWO*c*Aij)*dc+&
              (c*A1j-a*dj)*ddi+(c*A1i-a*di)*ddj+&
              (c*dj-beta*A1j)*dA(1, parti, :, :)+(c*di-beta*A1i)*dA(1, partj, :, :)

          dGn = dG/(a*beta)-Gn*(dA(1, 1, :, :)/a+dbeta/beta)

          dcoab = -TWO*c*dc/(a*beta)+kappa*(dA(1, 1, :, :)/a+dbeta/beta)

          dCpq = Cpq*(-(THREE*dA(parti, partj, :, :)+TWO*m*dhoa)/g4den-dbeta/(TWO*beta))+&
              Glob_gm2(m)/(g4den*SQRT(beta))*(-(dGn*Q0+Gn*Q1*dcoab))

          dCsum = dCsum+Glob_ChargeMatrix(p, q)*Skl*dCpq
          chpq = chpq+Glob_ChargeMatrix(p, q)*Cpq
        ENDDO
      ENDDO

      CALL AddCorrToGrad(side, SymMatrix, chpq, dCsum, dSk, dSl, dHk, dHl)

    ENDDO

  END SUBROUTINE GramTermCorr





  SUBROUTINE GramTermCorrPairs(mk, ik, vechLk, ml, jl, vechLl, SymMatrix, &
                               Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! GramTermCorr for 5 or more pseudoparticles: the same operations, but
    ! the power sums Q0, Q1 of all pairs with the term are formed together
    ! (one element of the On* vectors per pair) and once per call, for the
    ! energy and both gradient sides. With fewer pseudoparticles the
    ! compiler would fuse some products of the gradients differently, so
    ! there GramTermCorr is used.
    ! Adds to the output of MatrixElementsXi (one ket term x_jl) the
    ! Gram-determinant term of its Coulomb integrals (see GramTermSetup):
    !   R(p,q) += Cpq*Skl,  Cpq = -Gn*Q0*gm2(m)/((3 Aij + 2 m hoa) sqrt(beta)),
    ! for every pair where the term is present, so Hkl gets
    ! Skl * sum charge(p,q)*Cpq, and the gradients get
    ! sum charge(p,q)*(Cpq*dS + Skl*dCpq), with
    !   dCpq  = Cpq*(-d(3 Aij + 2 m hoa)/(3 Aij + 2 m hoa) - dbeta/(2 beta))
    !           + gm2(m)/((3 Aij + 2 m hoa) sqrt(beta)) * d(-Gn*Q0),
    !   dG    = sum of the cofactors of the 3x3 Gram matrix times the
    !           derivatives of its elements,
    !   dGn   = dG/(a*beta) - Gn*(da/a + dbeta/beta),
    !   dcoab = -2 c dc/(a*beta) + (1 - coab)*(da/a + dbeta/beta),
    !   d(-Gn*Q0) = -(dGn*Q0 + Gn*Q1*dcoab).
    ! The derivatives of invAkl are those of MatrixElementsXi: with respect
    ! to Lk for dHk, with respect to PLl = SymMatrix'*Ll (then times
    ! SymMatrix) for dHl. Nothing is divided by c, d_i or A1i, which may
    ! vanish. invAkl is recomputed (CoulombCorrSetup) so that
    ! MatrixElementsXi, and its results where the term is absent, stay
    ! exactly as they were; this routine is called only when the term can
    ! be present (mk+ml > 0, ik /= 1, jl /= 1, n > 2).
    ! Input:
    !   mk, ik, vechLk  :: the bra
    !   ml, jl, vechLl  :: the ket term (power, x-index in the transformed
    !                      ket, exponent parameters)
    !   SymMatrix       :: the symmetry permutation matrix
    !   Skl, dSk, dSl   :: overlap and its gradients from MatrixElementsXi
    !   gradflag        :: as in MatrixElementsXi
    ! Input/Output:
    !   Hkl, dHk, dHl   :: Hamiltonian and its gradients
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, jl, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(IN) :: Skl, dSk(Glob_npt), dSl(Glob_npt)
    REAL(wp), INTENT(INOUT) :: Hkl, dHk(Glob_npt), dHl(Glob_npt)

    ! Parameters
    INTEGER, PARAMETER :: NumPairs = n*(n+1)/2-1  ! the pairs (p, q), 2 <= p <= n, q <= p

    ! Local variables
    INTEGER :: i, j, k, p, q, m, side, parti, partj, indx, t, NumOn
    LOGICAL :: IsOn
    REAL(wp) :: Gn, Q0, Q1, Cpq, chpq, g4den, coab
    REAL(wp) :: a, Aij, A1i, A1j, beta, c, di, dj, kappa
    REAL(wp) :: Lk(n, n), PLl(n, n), invAkl(n, n), dA(n, n, n, n)
    REAL(wp) :: dbeta(n, n), dc(n, n), ddi(n, n), ddj(n, n), dhoa(n, n)
    REAL(wp) :: dG(n, n), dGn(n, n), dcoab(n, n), dCpq(n, n), dCsum(n, n), dCl(n, n)
    INTEGER :: OnP(NumPairs)      ! p of each pair with the term
    INTEGER :: OnQ(NumPairs)      ! q of each pair with the term
    REAL(wp) :: OnGn(NumPairs)    ! Gn of each pair with the term
    REAL(wp) :: OnCoab(NumPairs)  ! coab of each pair with the term
    REAL(wp) :: OnW(NumPairs)     ! coab**(t-1), then coab**t, in the power sums
    REAL(wp) :: OnQ0(NumPairs)    ! Q0 of each pair with the term
    REAL(wp) :: OnQ1(NumPairs)    ! Q1 = dQ0/dcoab of each pair with the term


    m = (mk+ml)/2
    parti = ik
    partj = jl

    CALL CoulombCorrSetup(vechLk, vechLl, SymMatrix, Lk, PLl, invAkl)

    a = invAkl(1, 1)
    Aij = invAkl(parti, partj)
    A1i = invAkl(1, parti)
    A1j = invAkl(1, partj)
    g4den = THREE*Aij+TWO*m*A1i*A1j/a

    ! The pairs where the term is present, in the order of the sums below
    NumOn = 0
    DO p = 2, n
      DO q = 1, p
        CALL GramTermSetup(invAkl, parti, partj, m, p, q, IsOn, Gn, coab)
        IF (.NOT. IsOn) CYCLE
        NumOn = NumOn+1
        OnP(NumOn) = p
        OnQ(NumOn) = q
        OnGn(NumOn) = Gn
        OnCoab(NumOn) = coab
      ENDDO
    ENDDO

    ! The power sums Q0 and, for the gradients, Q1 of these pairs, advanced
    ! together (each element gets the operations of a sum done alone);
    ! OnW = coab^(t-1) at the top of the loop
    OnW(1:NumOn) = ONE
    OnQ0(1:NumOn) = ONE
    OnQ1(1:NumOn) = ZERO
    IF (gradflag == 0) THEN
      DO t = 1, m-1
        DO i = 1, NumOn
          OnW(i) = OnW(i)*OnCoab(i)
          OnQ0(i) = OnQ0(i)+(TWO*t+ONE)*Glob_gm3(t)*OnW(i)
        ENDDO
      ENDDO
    ELSE
      DO t = 1, m-1
        DO i = 1, NumOn
          OnQ1(i) = OnQ1(i)+(TWO*t+ONE)*Glob_gm3(t)*t*OnW(i)
          OnW(i) = OnW(i)*OnCoab(i)
          OnQ0(i) = OnQ0(i)+(TWO*t+ONE)*Glob_gm3(t)*OnW(i)
        ENDDO
      ENDDO
    ENDIF

    ! The energy
    DO i = 1, NumOn
      p = OnP(i)
      q = OnQ(i)
      IF (p == q) THEN
        beta = invAkl(p, p)
      ELSE
        beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
      ENDIF
      Hkl = Hkl-Glob_ChargeMatrix(p, q)*OnGn(i)*OnQ0(i)*Glob_gm2(m)/(g4den*SQRT(beta))*Skl
    ENDDO

    IF (gradflag == 0) RETURN

    ! The gradients: side 1 with respect to Lk, side 2 with respect to PLl
    DO side = 1, MIN(gradflag, 2)

      IF (side == 1) THEN
        CALL CoulombCorrDerivs(invAkl, Lk, dA)
      ELSE
        CALL CoulombCorrDerivs(invAkl, PLl, dA)
      ENDIF

      dhoa = (dA(1, parti, :, :)*A1j+A1i*dA(1, partj, :, :))/a-A1i*A1j*dA(1, 1, :, :)/(a*a)

      ! sum over the pairs of charge*Skl*dCpq in dCsum and of charge*Cpq in
      ! chpq (its product with dS is added to the vector below)
      dCsum = ZERO
      chpq = ZERO
      DO i = 1, NumOn
        p = OnP(i)
        q = OnQ(i)
        Gn = OnGn(i)
        Q0 = OnQ0(i)
        Q1 = OnQ1(i)

        IF (p == q) THEN
          beta = invAkl(p, p)
          c = invAkl(1, p)
          di = invAkl(parti, p)
          dj = invAkl(partj, p)
          dbeta = dA(p, p, :, :)
          dc = dA(1, p, :, :)
          ddi = dA(parti, p, :, :)
          ddj = dA(partj, p, :, :)
        ELSE
          beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
          c = invAkl(1, p)-invAkl(1, q)
          di = invAkl(parti, p)-invAkl(parti, q)
          dj = invAkl(partj, p)-invAkl(partj, q)
          dbeta = dA(p, p, :, :)+dA(q, q, :, :)-TWO*dA(p, q, :, :)
          dc = dA(1, p, :, :)-dA(1, q, :, :)
          ddi = dA(parti, p, :, :)-dA(parti, q, :, :)
          ddj = dA(partj, p, :, :)-dA(partj, q, :, :)
        ENDIF

        kappa = c*c/(a*beta)
        Cpq = -Gn*Q0*Glob_gm2(m)/(g4den*SQRT(beta))

        dG = (beta*Aij-di*dj)*dA(1, 1, :, :)+(a*Aij-A1i*A1j)*dbeta+&
            (a*beta-c*c)*dA(parti, partj, :, :)+(A1i*dj+di*A1j-TWO*c*Aij)*dc+&
            (c*A1j-a*dj)*ddi+(c*A1i-a*di)*ddj+&
            (c*dj-beta*A1j)*dA(1, parti, :, :)+(c*di-beta*A1i)*dA(1, partj, :, :)

        dGn = dG/(a*beta)-Gn*(dA(1, 1, :, :)/a+dbeta/beta)

        dcoab = -TWO*c*dc/(a*beta)+kappa*(dA(1, 1, :, :)/a+dbeta/beta)

        dCpq = Cpq*(-(THREE*dA(parti, partj, :, :)+TWO*m*dhoa)/g4den-dbeta/(TWO*beta))+&
            Glob_gm2(m)/(g4den*SQRT(beta))*(-(dGn*Q0+Gn*Q1*dcoab))

        dCsum = dCsum+Glob_ChargeMatrix(p, q)*Skl*dCpq
        chpq = chpq+Glob_ChargeMatrix(p, q)*Cpq
      ENDDO

      CALL AddCorrToGrad(side, SymMatrix, chpq, dCsum, dSk, dSl, dHk, dHl)

    ENDDO

  END SUBROUTINE GramTermCorrPairs





  SUBROUTINE PowerZeroGradCorr(ik, jl, vechLk, vechLl, SymMatrix, Skl, dHk, dHl, gradflag)
    ! For mk = ml = 0 (m = 0) the Coulomb gradients of MatrixElementsXi
    ! (dRk and dRl of the pairs (p, q), p >= 2) add the term t = m of the
    ! bracket of eq. (35) separately from the term t = 0, which for m = 0
    ! is the same term, 3 Aij - fob: its derivative is counted twice, the
    ! second time as
    !   charge(p,q) * Skl*gm2(0)/(3 Aij sqrt(beta)) * coab*(3 dAij - dfob)
    ! (Glob_gm3(0) = 1, coab^0 = 1). This routine subtracts it; the energy
    ! is not affected. As in GramTermCorr, invAkl and its derivatives are
    ! recomputed, so that MatrixElementsXi and its results for m > 0 stay
    ! exactly as they were; dfob is taken without dividing by d_i or d_j:
    !   dfob = (dd_i*d_j + d_i*dd_j)/beta - fob*dbeta/beta.
    ! Input:
    !   ik, jl          :: x-indices of the bra and the ket term
    !   vechLk, vechLl  :: n(n+1)/2 exponent parameters
    !   SymMatrix       :: the symmetry permutation matrix
    !   Skl             :: the overlap from MatrixElementsXi
    !   gradflag        :: 1 dHk, 2 dHk and dHl (nothing for 0)
    ! Input/Output:
    !   dHk, dHl        :: Hamiltonian gradients
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: ik, jl, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(IN) :: Skl
    REAL(wp), INTENT(INOUT) :: dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: p, q, side, parti, partj
    REAL(wp) :: Aij, beta, c, di, dj, coab, fob
    REAL(wp) :: Lk(n, n), PLl(n, n), invAkl(n, n), dA(n, n, n, n)
    REAL(wp) :: dbeta(n, n), ddi(n, n), ddj(n, n), dfob(n, n), dCsum(n, n)
    REAL(wp) :: dSdummy(Glob_npt)


    IF (gradflag == 0) RETURN

    parti = ik
    partj = jl
    dSdummy = ZERO

    CALL CoulombCorrSetup(vechLk, vechLl, SymMatrix, Lk, PLl, invAkl)
    Aij = invAkl(parti, partj)

    DO side = 1, MIN(gradflag, 2)

      IF (side == 1) THEN
        CALL CoulombCorrDerivs(invAkl, Lk, dA)
      ELSE
        CALL CoulombCorrDerivs(invAkl, PLl, dA)
      ENDIF

      dCsum = ZERO
      DO p = 2, n
        DO q = 1, p

          IF (p == q) THEN
            beta = invAkl(p, p)
            c = invAkl(1, p)
            di = invAkl(parti, p)
            dj = invAkl(partj, p)
            dbeta = dA(p, p, :, :)
            ddi = dA(parti, p, :, :)
            ddj = dA(partj, p, :, :)
          ELSE
            beta = invAkl(p, p)+invAkl(q, q)-TWO*invAkl(p, q)
            c = invAkl(1, p)-invAkl(1, q)
            di = invAkl(parti, p)-invAkl(parti, q)
            dj = invAkl(partj, p)-invAkl(partj, q)
            dbeta = dA(p, p, :, :)+dA(q, q, :, :)-TWO*dA(p, q, :, :)
            ddi = dA(parti, p, :, :)-dA(parti, q, :, :)
            ddj = dA(partj, p, :, :)-dA(partj, q, :, :)
          ENDIF

          coab = ONE-c*c/(invAkl(1, 1)*beta)
          fob = di*dj/beta
          dfob = (ddi*dj+di*ddj)/beta-fob*dbeta/beta

          dCsum = dCsum-Glob_ChargeMatrix(p, q)*Skl*Glob_gm2(0)/(THREE*Aij*SQRT(beta))*&
              coab*(THREE*dA(parti, partj, :, :)-dfob)
        ENDDO
      ENDDO

      CALL AddCorrToGrad(side, SymMatrix, ZERO, dCsum, dSdummy, dSdummy, dHk, dHl)

    ENDDO

  END SUBROUTINE PowerZeroGradCorr





  SUBROUTINE AddCorrToGrad(side, SymMatrix, chpq, dCsum, dSk, dSl, dHk, dHl)
    ! Adds a correction of the Coulomb gradient, chpq*dS + dCsum, to dHk
    ! (side 1, dCsum with respect to Lk) or to dHl (side 2, dCsum with
    ! respect to PLl, taken to Ll with SymMatrix as in MatrixElementsXi);
    ! the lower triangle is loaded into the vector as vech.
    ! Input:
    !   side            :: 1 bra, 2 ket
    !   SymMatrix       :: the symmetry permutation matrix
    !   chpq            :: factor of the overlap gradient
    !   dCsum           :: the correction as a matrix
    !   dSk, dSl        :: overlap gradients from MatrixElementsXi
    ! Input/Output:
    !   dHk, dHl        :: Hamiltonian gradients
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: side
    REAL(wp), INTENT(IN) :: SymMatrix(n, n), chpq, dCsum(n, n)
    REAL(wp), INTENT(IN) :: dSk(Glob_npt), dSl(Glob_npt)
    REAL(wp), INTENT(INOUT) :: dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: i, j, k, indx
    REAL(wp) :: dCl(n, n)


    IF (side == 1) THEN
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dHk(indx) = dHk(indx)+chpq*dSk(indx)+dCsum(i, j)
        ENDDO
      ENDDO
    ELSE
      ! d/dLl = SymMatrix * d/dPLl
      dCl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            dCl(i, j) = dCl(i, j)+SymMatrix(i, k)*dCsum(k, j)
          ENDDO
        ENDDO
      ENDDO
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dHl(indx) = dHl(indx)+chpq*dSl(indx)+dCl(i, j)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE AddCorrToGrad





  SUBROUTINE KetNormRatios(ml, il, vechLl, SymMatrix, NumKetTerms, KetIndex, &
                           KetRatio, dlnKetRatio, IsGradNeeded)
    ! When the symmetry operator turns the ket prefactor x_il into
    ! sum_t w_t x_j(t) (w = row il of SymMatrix, j(t) = KetIndex(t)),
    ! MatrixElementsXi and ExpValuesXi normalize term t with the norm N_j(t)
    ! of x_j(t) r_1^ml exp[-r'(Al (x) I3) r], Al = PLl*PLl', PLl =
    ! SymMatrix'*Ll, whereas the ket has the norm N_w of (w'x) r_1^ml
    ! exp[...], which is the norm of f_l itself. With B = inv(Al), the
    ! norm of (u'x) r_1^ml exp[...] is, up to a factor common to all u,
    !   N_u = u'B u + (2 ml/3) (e1'B u)^2/B(1,1),
    ! and with P = PLl
    !   d(u'B v)/dP(a,b) = -[(B u)_a (P'B v)_b + (B v)_a (P'B u)_b].
    ! The derivatives with respect to PLl are taken to Ll with SymMatrix,
    ! as in MatrixElementsXi.
    ! Input:
    !   ml            :: premultiplier power of the ket
    !   il            :: x-index of the ket before the transformation
    !   vechLl        :: n(n+1)/2 exponent parameters of the ket
    !   SymMatrix     :: the symmetry permutation matrix
    !   NumKetTerms   :: number of terms
    !   KetIndex      :: x-index j(t) of each term
    !   IsGradNeeded  :: compute dlnKetRatio
    ! Output:
    !   KetRatio(t)   :: sqrt(N_j(t)/N_w)
    !   dlnKetRatio(:, t) :: d ln KetRatio(t)/d vech[Ll] (if IsGradNeeded)
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: ml, il, NumKetTerms
    REAL(wp), INTENT(IN) :: vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    INTEGER, INTENT(IN) :: KetIndex(n)
    REAL(wp), INTENT(OUT) :: KetRatio(n), dlnKetRatio(Glob_npt, n)
    LOGICAL, INTENT(IN) :: IsGradNeeded

    ! Local variables
    INTEGER :: i, j, k, t, indx
    REAL(wp) :: SUM1, Nw, Nt
    REAL(wp) :: Ll(n, n), PLl(n, n), Al(n, n), chAl(n, n), B(n, n)
    REAL(wp) :: u(n), dNw(n, n), dNt(n, n), dlnP(n, n), dln(n, n)


    ! PLl = SymMatrix'*Ll and Al = PLl*PLl'
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Ll(i, j) = vechLl(indx)
      ENDDO
    ENDDO

    PLl = ZERO
    DO j = 1, n
      DO k = j, n
        DO i = 1, n
          PLl(i, j) = PLl(i, j)+SymMatrix(k, i)*Ll(k, j)
        ENDDO
      ENDDO
    ENDDO

    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Al(i, j) = Al(i, j)+PLl(i, k)*PLl(j, k)
        ENDDO
      ENDDO
    ENDDO

    ! B = inv(Al) from the Cholesky factor of Al, as in MatrixElementsXi
    chAl = Al
    DO i = 1, n
      DO j = i, n
        SUM1 = chAl(i, j)
        DO k = i-1, 1, -1
          SUM1 = SUM1-chAl(i, k)*chAl(j, k)
        ENDDO
        IF (i == j) THEN
          chAl(i, i) = SQRT(SUM1)
        ELSE
          chAl(j, i) = SUM1/chAl(i, i)
        ENDIF
      ENDDO
    ENDDO

    DO i = 1, n
      chAl(i, i) = ONE/chAl(i, i)
      DO j = i+1, n
        SUM1 = ZERO
        DO k = i, j-1
          SUM1 = SUM1-chAl(j, k)*chAl(k, i)
        ENDDO
        chAl(j, i) = SUM1/chAl(j, j)
        chAl(i, j) = ZERO
      ENDDO
    ENDDO

    B = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          B(i, j) = B(i, j)+chAl(k, i)*chAl(k, j)
        ENDDO
      ENDDO
    ENDDO

    ! Norm of the ket, u = w
    u = SymMatrix(il, 1:n)
    CALL NormOfKetTerm(u, Nw, dNw)

    DO t = 1, NumKetTerms

      u = ZERO
      u(KetIndex(t)) = ONE
      CALL NormOfKetTerm(u, Nt, dNt)

      KetRatio(t) = SQRT(Nt/Nw)

      IF (IsGradNeeded) THEN

        dlnP = ONEHALF*(dNt/Nt-dNw/Nw)

        ! d/dLl = SymMatrix * d/dPLl, then load it into a vector
        dln = ZERO
        DO k = 1, n
          DO j = 1, n
            DO i = 1, n
              dln(i, j) = dln(i, j)+SymMatrix(i, k)*dlnP(k, j)
            ENDDO
          ENDDO
        ENDDO
        indx = 0
        DO j = 1, n
          DO i = j, n
            indx = indx+1
            dlnKetRatio(indx, t) = dln(i, j)
          ENDDO
        ENDDO

      ENDIF

    ENDDO

  CONTAINS

    SUBROUTINE NormOfKetTerm(uu, Nu, dNu)
      ! N_u and dN_u/dPLl for the prefactor u'x (see KetNormRatios)
      REAL(wp), INTENT(IN) :: uu(n)
      REAL(wp), INTENT(OUT) :: Nu, dNu(n, n)

      REAL(wp) :: Bu(n), PBu(n), Be1(n), PBe1(n), quu, q1u, q11, c2
      INTEGER :: aa, bb

      Bu = MATMUL(B, uu)
      PBu = MATMUL(TRANSPOSE(PLl), Bu)
      Be1 = B(1:n, 1)
      PBe1 = MATMUL(TRANSPOSE(PLl), Be1)

      quu = DOT_PRODUCT(uu, Bu)
      q1u = Bu(1)
      q11 = B(1, 1)
      c2 = TWO*ml/THREE

      Nu = quu+c2*q1u*q1u/q11

      ! d(u'B v)/dP(a,b) = -[(Bu)_a (P'Bv)_b + (Bv)_a (P'Bu)_b]
      DO bb = 1, n
        DO aa = 1, n
          dNu(aa, bb) = -TWO*Bu(aa)*PBu(bb)&
              +c2*(TWO*q1u/q11*(-(Be1(aa)*PBu(bb)+Bu(aa)*PBe1(bb)))&
              -q1u*q1u/(q11*q11)*(-TWO*Be1(aa)*PBe1(bb)))
        ENDDO
      ENDDO

    END SUBROUTINE NormOfKetTerm

  END SUBROUTINE KetNormRatios





  SUBROUTINE ExpValuesMatElem(mk, ik, vechLk, ml, il, vechLl, SymMatrix, &
                              NumCFGridPoints, CFGridX, ycf, CFkl, AreCorrFuncNeeded, &
                              ArePairDistNeeded, rkl, rkl_2)
    ! Matrix elements of the expectation-value operators between two basis
    ! functions x_ik r_1^mk exp[...] and x_il r_1^ml exp[...]: the
    ! nucleus-nucleus correlation function on one row of the (x,y) grid,
    ! and the interparticle distances <r_ij> and <r_ij^2>. The symmetry
    ! operator turns the ket prefactor x_il into sum_j SymMatrix(il, j) x_j,
    ! handled as in MatrixElements: one call of ExpValuesXi with the sign
    ! when row il of SymMatrix has one entry +-1, otherwise the terms
    ! rescaled with KetNormRatios and summed.
    ! Input:
    !   mk, ml           :: premultiplier powers of the bra and the ket
    !   ik, il           :: x-indices of the prefactors of the bra and the ket
    !   vechLk, vechLl   :: n(n+1)/2 exponent parameters
    !   SymMatrix        :: the symmetry permutation matrix (ket)
    !   NumCFGridPoints  :: number of x grid points of the correlation function
    !   CFGridX          :: the x grid points
    !   ycf              :: the y grid point of this call
    !   AreCorrFuncNeeded  :: compute the correlation function at all
    !   ArePairDistNeeded  :: compute rkl and rkl_2 (they do not depend on
    !                         the grid, so once per pair is enough)
    ! Output (unnormalized), only those asked for:
    !   CFkl(1:NumCFGridPoints) :: correlation function at (CFGridX(i), ycf)
    !   rkl, rkl_2       :: <r_ij> and <r_ij^2>; the diagonal (i,i) is the
    !                       distance of pseudoparticle i from particle 0
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, il
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    INTEGER, INTENT(IN) :: NumCFGridPoints
    REAL(wp), INTENT(IN) :: CFGridX(NumCFGridPoints)
    REAL(wp), INTENT(IN) :: ycf
    LOGICAL, INTENT(IN) :: AreCorrFuncNeeded, ArePairDistNeeded
    REAL(wp), INTENT(OUT) :: CFkl(NumCFGridPoints)
    REAL(wp), INTENT(OUT) :: rkl(n, n), rkl_2(n, n)

    ! Local variables
    INTEGER :: j, t, NumKetTerms
    INTEGER :: KetIndex(n)                        ! x-index j of each ket term
    REAL(wp) :: KetWeight(n)                      ! SymMatrix(il, j) of each ket term
    REAL(wp) :: KetRatio(n)                       ! sqrt(N_j/N_ket), see KetNormRatios
    REAL(wp) :: dlnKetRatio(Glob_npt, n)          ! not used here
    REAL(wp) :: f
    REAL(wp) :: CFt(NumCFGridPoints), rt(n, n), r2t(n, n)


    ! The terms of the transformed prefactor: row il of SymMatrix
    NumKetTerms = 0
    DO j = 1, n
      IF (SymMatrix(il, j) /= ZERO) THEN
        NumKetTerms = NumKetTerms+1
        KetIndex(NumKetTerms) = j
        KetWeight(NumKetTerms) = SymMatrix(il, j)
      ENDIF
    ENDDO

    ! One term: x_il --> +-x_j
    IF ((NumKetTerms == 1) .AND. (ABS(KetWeight(1)) == ONE)) THEN
      CALL ExpValuesXi(mk, ik, vechLk, ml, KetIndex(1), vechLl, SymMatrix, &
                       NumCFGridPoints, CFGridX, ycf, CFkl, AreCorrFuncNeeded, &
                       ArePairDistNeeded, rkl, rkl_2)
      IF (KetWeight(1) < ZERO) THEN
        IF (AreCorrFuncNeeded) CFkl = -CFkl
        IF (ArePairDistNeeded) THEN
          rkl = -rkl
          rkl_2 = -rkl_2
        ENDIF
      ENDIF
      RETURN
    ENDIF

    ! Several terms: x_il --> sum_t KetWeight(t) x_KetIndex(t)
    CALL KetNormRatios(ml, il, vechLl, SymMatrix, NumKetTerms, KetIndex, &
                       KetRatio, dlnKetRatio, .FALSE.)

    IF (AreCorrFuncNeeded) CFkl = ZERO
    IF (ArePairDistNeeded) THEN
      rkl = ZERO
      rkl_2 = ZERO
    ENDIF
    DO t = 1, NumKetTerms
      CALL ExpValuesXi(mk, ik, vechLk, ml, KetIndex(t), vechLl, SymMatrix, &
                       NumCFGridPoints, CFGridX, ycf, CFt, AreCorrFuncNeeded, &
                       ArePairDistNeeded, rt, r2t)
      f = KetWeight(t)*KetRatio(t)
      IF (AreCorrFuncNeeded) CFkl = CFkl+f*CFt
      IF (ArePairDistNeeded) THEN
        rkl = rkl+f*rt
        rkl_2 = rkl_2+f*r2t
      ENDIF
    ENDDO

  END SUBROUTINE ExpValuesMatElem





  SUBROUTINE ExpValuesXi(mk, ik, vechLk, ml, jl, vechLl, SymMatrix, &
                         NumCFGridPoints, CFGridX, ycf, CFkl, AreCorrFuncNeeded, &
                         ArePairDistNeeded, rkl, rkl_2)
    ! The matrix elements of ExpValuesMatElem for one term x_jl of the
    ! prefactor of the symmetry transformed ket, x_jl r_1^ml exp[-r'(Al (x)
    ! I3) r], normalized with its own x_jl.
    ! For ik = jl = 1 the formulas are those of the x_1 r_1^m basis
    ! (Sharkey, Kirnosov and Adamowicz, J. Chem. Phys. 139, 164119 (2013)),
    ! exactly as in the original code. Otherwise, with p = (mk+ml)/2,
    ! a = invAkl(1,1), Aij = invAkl(i,j), A1i = invAkl(1,i) (i = ik,
    ! j = jl), hoa = A1i*A1j/a, P0 = Aij + (2p/3) hoa and, for the pair
    ! vector v (v = e_k, or e_k - e_l), b = v'*invAkl*v, cv = (invAkl*v)(1),
    ! d_i = (invAkl*v)(i), kap = cv^2/(a*b), w = 1 - kap, fob = d_i*d_j/b,
    ! e = cv*(A1i*d_j + d_i*A1j)/(a*b):
    !   correlation function (the density of r_1 = (x, y, 0), x along the
    !   axis of the prefactors; E[x_i x_j | r_1] = (A1i A1j/a^2) x^2
    !   + (Aij - hoa)/2):
    !     g(x,y) = Skl/(pi P0) [(A1i A1j/a^2) x^2 + (Aij - hoa)/2]
    !              * (r^2/a)^p exp(-r^2/a) a^(-3/2)/Gamma(p+3/2)
    !   <r^2> = Skl [b (3/2 + (p-1) kap) + b (kap Aij + fob + (2p/3) e)/P0]
    !   <r>   = Skl gm2(p) sqrt(b)/(3 P0)
    !           * { sum_{t=0}^{p-1} gm3(t) w^t [(1 + (p-t) kap) T_t - D_t]
    !               + gm3(p) w^p [(2p+3) Aij - (2p-1) fob] }
    !     T_t = (2t+3) Aij - (2t+1)(w Aij + e) + (2p+1) hoa
    !     D_t = -(2t+3) fob - (2t+1)(kap w Aij - w fob - w e - 2 kap fob)
    !           + (2p+1)(kap hoa - e)
    ! <r> is -d/ds of the Coulomb integral of MatrixElementsXi (with the
    ! Gram term) for invAkl --> inv(Akl + s v v'); for x_1 all three reduce
    ! to the x_1 formulas.
    ! Input and output: as ExpValuesMatElem, with jl in place of il.
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ik, ml, jl
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    INTEGER, INTENT(IN) :: NumCFGridPoints
    REAL(wp), INTENT(IN) :: CFGridX(NumCFGridPoints)
    REAL(wp), INTENT(IN) :: ycf
    LOGICAL, INTENT(IN) :: AreCorrFuncNeeded, ArePairDistNeeded
    REAL(wp), INTENT(OUT) :: CFkl(NumCFGridPoints)
    REAL(wp), INTENT(OUT) :: rkl(n, n), rkl_2(n, n)

    ! Local variables
    INTEGER :: i, j, k, indx, info, jj, ii, kk, mko2, mlo2
    REAL(wp) :: detLk, detLl, detAkl, Skl
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3
    REAL(wp) :: Lk(n, n), Ll(n, n), Ak(n, n), Al(n, n), invAkl(n, n), invL(n, n), PLl(n, n)
    INTEGER :: p, q, t
    REAL(wp) :: invLk(n, n), invLl(n, n), invAk(n, n), invAl(n, n), chAl(n, n)
    REAL(wp) :: g1, g2, chi1, g3
    REAL(wp) :: a, coabij, r1kl
    REAL(wp) :: Temp(n, n), R(n, n), b(n, n), c(n, n), coab(n, n), Vtmp(n, n)

    ! Prefactors x_i (bra) and x_j (ket term) other than x_1 x_1
    LOGICAL :: IsX1X1
    INTEGER :: parti, partj
    REAL(wp) :: Aij, A1i, A1j, hoa, P0, SklOverP0
    REAL(wp) :: cv, di, dj, kap, w, fob, e, Tt, Dt

    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Build the lower triangular matrices Lk and Ll from
    ! the elements in vechLk and vechLl
    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
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
    ! Build AK=LK*LL' and AL=P'*LL*LL'*P
    Ak = ZERO
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Ak(i, j) = Ak(i, j)+Lk(i, k)*Lk(j, k)
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
    ! need to do the decomp on Al too to get invAl since
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
    ! Now do invAkl=invL'*invL
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
    ! Compute Skl
    !
    ! Skl = 2 * (3+mk+ml)/sqrt(9+6*(mk+ml)+4*mk*ml) *
    !      gm1(mk,ml) 2^(3n/2) *
    !      sqrt( (invAkl(1,1)/invAk(1,1))^(mk+1)  *
    !            (invAkl(1,1)/invAl(1,1))^(ml+2)  *
    !            (||Lk||||Ll||/|Akl|)^3            )
    !********************************************************
    mko2 = mk/2
    mlo2 = ml/2

    parti = ik
    partj = jl
    IsX1X1 = (parti == 1) .AND. (partj == 1)

    IF (IsX1X1) THEN

      chi1 = ZERO
      chi1 = chi1+(NINE+SIX*(mk+ml)+FOUR*mk*ml)
      chi1 = ((THREE+mk+ml)/SQRT(chi1))
      Skl = ZERO
      Skl = Skl+SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))
      Skl = Skl*((invAkl(1, 1)/invAk(1, 1))**(mko2))*((invAkl(1, 1)/invAl(1, 1))**(mlo2))
      Skl = Skl*Glob_gm1(mko2, mlo2)*Glob_gm1(mko2, mlo2)
      Skl = Skl*SQRT((invAkl(1, 1)/invAk(1, 1)))*SQRT((invAkl(1, 1)/invAl(1, 1)))
      Skl = Skl*TWO*chi1

    ELSE

      ! The overlap of MatrixElementsXi, kept as Skl/P0 and P0 so that the
      ! formulas below do not divide by P0
      A1i = invAkl(1, parti)
      A1j = invAkl(1, partj)
      Aij = invAkl(parti, partj)
      hoa = A1i*A1j/invAkl(1, 1)
      P0 = Aij+(mk+ml)*ONETHIRD*hoa

      SklOverP0 = (invAk(parti, parti)+(TWO*mk)/THREE*invAk(1, parti)*invAk(1, parti)/invAk(1, 1))
      SklOverP0 = SklOverP0*(invAl(partj, partj)+(TWO*ml)/THREE*invAl(1, partj)*invAl(1, partj)/invAl(1, 1))
      SklOverP0 = TWO/SQRT(SklOverP0)
      SklOverP0 = SklOverP0*SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))
      SklOverP0 = SklOverP0*((invAkl(1, 1)/invAk(1, 1))**(mko2))*((invAkl(1, 1)/invAl(1, 1))**(mlo2))
      SklOverP0 = SklOverP0*Glob_gm1(mko2, mlo2)*Glob_gm1(mko2, mlo2)
      Skl = SklOverP0*P0

    ENDIF

    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in Glob_gm1 (module globvars)

    ! Computing nucleus-nucleus correlation function if needed

    p = (mk+ml)/2

    IF (AreCorrFuncNeeded .AND. IsX1X1) THEN
      g1 = 1/invAkl(1, 1)
      DO i = 1, NumCFGridPoints
        g2 = g1*(CFGridX(i)*CFGridX(i)+ycf*ycf)
        IF (g2 == ZERO) THEN
          CFkl(i) = ZERO
        ELSE
          CFkl(i) = Skl*(ONEHALF/Glob_Pi)*(SIX/(THREE+TWO*p))*CFGridX(i)*CFGridX(i)
          CFkl(i) = CFkl(i)*EXP(-Glob_lngamma(p)-g2+p*LOG(g2)+FIVE*ONEHALF*LOG(g1))
        ENDIF
      ENDDO
    ENDIF

    ! The general prefactors (see the header); at r = 0 only the term
    ! (Aij - hoa)/2 with p = 0 is left
    IF (AreCorrFuncNeeded .AND. (.NOT. IsX1X1)) THEN
      g1 = 1/invAkl(1, 1)
      DO i = 1, NumCFGridPoints
        g2 = g1*(CFGridX(i)*CFGridX(i)+ycf*ycf)
        g3 = A1i*A1j*g1*g1*CFGridX(i)*CFGridX(i)+ONEHALF*(Aij-hoa)
        IF (g2 == ZERO) THEN
          IF (p == 0) THEN
            CFkl(i) = SklOverP0/Glob_Pi*g3*EXP(-Glob_lngamma(p)+THREE*ONEHALF*LOG(g1))
          ELSE
            CFkl(i) = ZERO
          ENDIF
        ELSE
          CFkl(i) = SklOverP0/Glob_Pi*g3
          CFkl(i) = CFkl(i)*EXP(-Glob_lngamma(p)-g2+p*LOG(g2)+THREE*ONEHALF*LOG(g1))
        ENDIF
      ENDDO
    ENDIF

    ! Computing average interparticle distances

    !********************************************************************
    ! Build intermediate results for <r_ij>
    ! a = invAkl(1,1)
    ! b(i,i) = invAkl(i,i)
    ! b(i,j) = invAkl(i,i) + invAkl(j,j) - 2 invAkl(i,j)
    ! c(i,i) = invAkl(1,i)^2
    ! c(i,j) = invAkl(1,i)^2 + invAkl(1,j)^2 - 2 invAkl(1,i) invAkl(1,j)
    !********************************************************************
    IF (ArePairDistNeeded) THEN
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

      ! coab=c/ab  a matrix
      coab = c/(a*b)

      IF (IsX1X1) THEN

        rkl = ZERO
        DO j = 1, n
          DO i = j, n
            g1 = ONE-coab(i, j)
            g2 = ONE
            g3 = ZERO
            DO k = 1, p
              g2 = g2*g1
              g3 = g3+Glob_gm3(k)*g2*(p-k+ONE)*(p-k+TWO)*(c(i, j)-THREE*a*b(i, j)/(TWO*k-ONE))
            ENDDO
            rkl(i, j) = Skl*Glob_gm2(p)*((p+ONE)*(p+TWO)*(THREE*a*b(i, j)+c(i, j))+g3)
            rkl(i, j) = rkl(i, j)/(FOUR*(p+THREE*ONEHALF)*SQRT(b(i, j))*a)
            rkl(j, i) = rkl(i, j)
          ENDDO
        ENDDO

        rkl_2 = ZERO
        DO j = 1, n
          DO i = j, n
            rkl_2(i, j) = Skl*(THREE*a*b(i, j)+TWO*(ONE+p)*c(i, j))/(TWO*a)
            rkl_2(j, i) = rkl_2(i, j)
          ENDDO
        ENDDO

      ELSE

        ! The general prefactors (see the header); coab here is kap
        rkl = ZERO
        rkl_2 = ZERO
        DO j = 1, n
          DO i = j, n

            IF (i == j) THEN
              cv = invAkl(1, i)
              di = invAkl(parti, i)
              dj = invAkl(partj, i)
            ELSE
              cv = invAkl(1, i)-invAkl(1, j)
              di = invAkl(parti, i)-invAkl(parti, j)
              dj = invAkl(partj, i)-invAkl(partj, j)
            ENDIF
            kap = coab(i, j)
            w = ONE-kap
            fob = di*dj/b(i, j)
            e = cv*(A1i*dj+di*A1j)/(a*b(i, j))

            ! <r>: the terms t = 0 .. p-1, then t = p
            g2 = ONE
            g3 = ZERO
            DO t = 0, p-1
              Tt = (TWO*t+THREE)*Aij-(TWO*t+ONE)*(w*Aij+e)+(TWO*p+ONE)*hoa
              Dt = -(TWO*t+THREE)*fob-(TWO*t+ONE)*(kap*w*Aij-w*fob-w*e-TWO*kap*fob)+&
                  (TWO*p+ONE)*(kap*hoa-e)
              g3 = g3+Glob_gm3(t)*g2*((ONE+(p-t)*kap)*Tt-Dt)
              g2 = g2*w
            ENDDO
            g3 = g3+Glob_gm3(p)*g2*((TWO*p+THREE)*Aij-(TWO*p-ONE)*fob)
            rkl(i, j) = SklOverP0*Glob_gm2(p)*SQRT(b(i, j))*g3/THREE
            rkl(j, i) = rkl(i, j)

            ! <r^2>
            rkl_2(i, j) = Skl*b(i, j)*(THREE*ONEHALF+(p-1)*kap)+&
                SklOverP0*b(i, j)*(kap*Aij+fob+TWO*p*e/THREE)
            rkl_2(j, i) = rkl_2(i, j)

          ENDDO
        ENDDO

      ENDIF

    ENDIF

  END SUBROUTINE ExpValuesXi





  SUBROUTINE MatrixElements0(mk, vechLk, ml, vechLl, SymMatrix, &
                             Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! Symmetry adapted matrix elements in the N = 0 basis of correlated
    ! Gaussians with r_1 premultipliers,
    !   f_k = r_1^mk * exp[-r'(Lk*Lk' (x) I3) r],
    ! i.e. the PG_0S basis. Kept from the original code, where it served the
    ! N = 0 -> N = 1 transition moments; this program does not call it.
    ! The symmetry adaption is applied to the ket:
    ! Al --> SymMatrix'*Ll*Ll'*SymMatrix.
    ! Input:
    !   mk, ml         :: premultiplier powers (even)
    !   vechLk, vechLl :: n(n+1)/2 exponent parameters (Lk, Ll lower triangular)
    !   SymMatrix      :: the symmetry permutation matrix
    !   gradflag       :: 0 no gradient; 1 dSk, dHk only; 2 also dSl, dHl
    ! Output (normalized):
    !   Skl, Hkl       :: overlap and Hamiltonian matrix elements
    !   dSk, dHk       :: derivatives with respect to vech[Lk]
    !   dSl, dHl       :: derivatives with respect to vech[Ll]
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ml, gradflag
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: Skl, Hkl
    REAL(wp), INTENT(OUT) :: dSk(Glob_npt), dSl(Glob_npt), dHk(Glob_npt), dHl(Glob_npt)

    ! Local variables
    INTEGER :: i, j, k, indx, info, jj, ii, kk, mko2, mlo2
    REAL(wp) :: detLk, detLl, detAkl
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3
    REAL(wp) :: Lk(n, n), Ll(n, n), Ak(n, n), Al(n, n), invAkl(n, n), PLl(n, n), invL(n, n)
    REAL(wp) :: Temp(n, n)
    REAL(wp) :: trT, tl, tk, mkml, mkpml, Ttmp
    REAL(wp) :: a, coabij
    REAL(wp) :: R(n, n), b(n, n), c(n, n), coab(n, n), Vtmp(n, n)
    INTEGER :: p, q, t
    REAL(wp) :: Tkl, Vkl
    REAL(wp) :: dTk(Glob_npt), dTl(Glob_npt), dVk(Glob_npt), dVl(Glob_npt)
    REAL(wp) :: ddSk(n, n), ddSl(n, n), ddTk(n, n), ddTl(n, n)
    REAL(wp) :: ddSk_over_Skl(n, n), ddSl_over_Skl(n, n)
    REAL(wp) :: dRk(n, n), dRl(n, n), ddVk(n, n), ddVl(n, n)
    REAL(wp) :: invLk(n, n), invLl(n, n), invAk(n, n), invAl(n, n), chAl(n, n)

    REAL(wp) :: ALk(n, n), ALl(n, n), AAk(n, n), AAl(n, n), AlLl(n, n), MAk(n, n), MAl(n, n)
    REAL(wp) :: Jkk(n, n), Jll(n, n), Jklk(n, n), Jkll(n, n)
    REAL(wp) :: AAlMAl(n, n), AAkMAk(n, n), AAkMAl(n, n), AAlMAk(n, n)
    REAL(wp) :: Tllk(n, n), Tkkl(n, n), Tklk(n, n), Tlkl(n, n)
    REAL(wp) :: TlMkJk(n, n), TkMlJl(n, n), TlMlJk(n, n), TkMkJl(n, n)
    REAL(wp) :: JkMlk(n, n), JlMkl(n, n), JlMlk(n, n), JkMkl(n, n)
    REAL(wp) :: da(n, n), db(n, n), dc(n, n)
    REAL(wp) :: RLk(n, n), RLl(n, n), JRLk(n, n), JRLl(n, n), RJLk(n, n), RJLl(n, n)
    REAL(wp) :: g1, g2, g3

    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Build the lower triangular matrices Lk and Ll from
    ! the elements in vechLk and vechLl
    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
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

    ! Build AK=LK*LL' and AL=P'*LL*LL'*P
    Ak = ZERO
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Ak(i, j) = Ak(i, j)+Lk(i, k)*Lk(j, k)
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
    ! need to do the decomp on Al too to get invAl since
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
    ! Now do invAkl=invL'*invL
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
          AAl(i, j) = AAl(i, j) +invAkl(i, k)*Al(k, j)  ! could do I-invAkl*Ak
        ENDDO
      ENDDO
    ENDDO

    ! Now a bunch of rank one mats and do it a column at a time!
    ! Jkk = invAk*J11*invAk*Lk
    ! Jll = invAl*J11*invAl*Ll
    ! Jklk = invAkl*J11*invAkl*Lk
    ! Jkll = invAkl*J11*invAkl*Ll
    DO j = 1, n
      Jkk(1:n, j) = invAk(1:n, 1)*invLk(j, 1)
      Jll(1:n, j) = invAl(1:n, 1)*AlLl(1, j)  ! using AlLl to get sym right
      Jklk(1:n, j) = invAkl(1:n, 1)*ALk(1, j)
      Jkll(1:n, j) = invAkl(1:n, 1)*ALl(1, j)
    ENDDO

    ! MAk=M*Ak
    ! MAl=M*Al
    MAk = ZERO
    MAl = ZERO
    DO k = 1, n
      DO j = 1, n
        DO i = 1, n
          MAk(i, j) = MAk(i, j)+Glob_MassMatrix(i, k)*Ak(k, j)
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

    ! More kinetic energy (gradient) matrix terms
    IF (gradflag > 0) THEN
      ! Tllk=invAkl*Al*M*Al*invAkl*Lk
      ! Tkkl=invAkl*Ak*M*Ak*invAkl*Ll
      ! Tklk=invAkl*Ak*M*Al*invAkl*Lk
      ! Tlkl=invAkl*Al*M*Ak*invAkl*Ll
      Tllk = ZERO
      Tkkl = ZERO
      Tklk = ZERO
      Tlkl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            Tllk(i, j) = Tllk(i, j)+AAlMAl(i, k)*ALk(k, j)
            Tkkl(i, j) = Tkkl(i, j)+AAkMAk(i, k)*ALl(k, j)
            Tklk(i, j) = Tklk(i, j)+AAkMAl(i, k)*ALk(k, j)
            Tlkl(i, j) = Tlkl(i, j)+AAlMAk(i, k)*ALl(k, j)
          ENDDO
        ENDDO
      ENDDO

      ! Still more kinetic energy (gradient) matrix terms
      ! TlMkJk = invAkl*Al*M*Al*invAkl*Lk
      ! TkMlJl = invAkl*Ak*M*Ak*invAkl*Ll
      ! TlMlJk = invAkl*Ak*M*Al*invAkl*Lk
      ! TkMkJl = invAkl*Al*M*Ak*invAkl*Ll
      TlMkJk = ZERO
      TkMlJl = ZERO
      TlMlJk = ZERO
      TkMkJl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            TlMkJk(i, j) = TlMkJk(i, j)+AAlMAk(i, k)*Jklk(k, j)
            TkMlJl(i, j) = TkMlJl(i, j)+AAkMAl(i, k)*Jkll(k, j)
            TlMlJk(i, j) = TlMlJk(i, j)+AAlMAl(i, k)*Jklk(k, j)
            TkMkJl(i, j) = TkMkJl(i, j)+AAkMAk(i, k)*Jkll(k, j)
          ENDDO
        ENDDO
      ENDDO

      ! yet more kinetic terms ... nice rank 1 ones
      ! JkMlk = invAkl*J11*invAkl*Ak*M*Al*invAkl*Lk
      ! JlMkl = invAkl*J11*invAkl*Al*M*Ak*invAkl*Ll
      ! JlMlk = invAkl*J11*invAkl*Al*M*Al*invAkl*Lk
      ! JkMkl = invAkl*J11*invAkl*Ak*M*Ak*invAkl*Ll
      DO j = 1, n
        JkMlk(1:n, j) = invAkl(1:n, 1)*Tklk(1, j)
        JlMkl(1:n, j) = invAkl(1:n, 1)*Tlkl(1, j)
        JlMlk(1:n, j) = invAkl(1:n, 1)*Tllk(1, j)
        JkMkl(1:n, j) = invAkl(1:n, 1)*Tkkl(1, j)
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
    !!$  Skl = ( invAkl(1,1)/invAk(1,1) )**mk
    !!$  Skl = Skl * ( ( invAkl(1,1)/invAl(1,1) )**ml)
    !!$  Skl = Skl * (DABS( detLk*detLl/detAkl )**3)
    !!$  Skl = DSQRT(Skl)
    !!$  Skl = Skl * DSQRT(DBLE(2**(3*n)))
    !!$  Skl = Skl * gm1(mk,ml)
    mko2 = mk/2
    mlo2 = ml/2
    Skl = (invAkl(1, 1)/invAk(1, 1))**(mko2)*Glob_gm1(mko2, mlo2)   &
        *((invAkl(1, 1)/invAl(1, 1))**(mlo2) * Glob_gm1(mko2, mlo2))  &
        *SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))


    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in Glob_gm1 (module globvars)


    !********************************************************
    ! gradient terms dSk dSl
    !********************************************************
    mkpml = ONE*(mk+ml)

    IF (gradflag > 0) THEN
      ! dSk
      DO j = 1, n
        DO i = j, n
          ddSk(i, j) = (ONE+ONEHALF)*invLk(j, i)-THREE*ALk(i, j)
        ENDDO
      ENDDO
      ddSk = ddSk+mk/invAk(1, 1)*Jkk
      ddSk = ddSk-mkpml/invAkl(1, 1)*Jklk
      ddSk_over_Skl = ddSk
      ddSk = Skl*ddSk
      ! load it into dSk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dSk(indx) = ddSk(i, j)
        ENDDO
      ENDDO

    ENDIF

    IF (gradflag > 1) THEN
      ! dSl
      Temp = ml/invAl(1, 1)*Jll
      Temp = Temp-mkpml/invAkl(1, 1)*Jkll
      Temp = Temp-THREE*ALl

      ! Multiply by SymMatrix
      ddSl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ddSl(i, j) = ddSl(i, j)+SymMatrix(i, k)*Temp(k, j)
          ENDDO
        ENDDO
      ENDDO

      DO j = 1, n
        DO i = j, n
          ddSl(i, j) = ddSl(i, j) + (ONE+ONEHALF) * invLl(j, i)
        ENDDO
      ENDDO

      ddSl_over_Skl = ddSl
      ddSl = Skl*ddSl
      ! Load it into dSl (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dSl(indx) = ddSl(i, j)
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
    mkml = ONE*mk*ml

    Tkl = mkml*Glob_MassMatrix(1, 1)/(mkpml+ONE)
    Tkl = Tkl-mk*tl
    Tkl = Tkl-ml*tk
    Ttmp = Tkl  ! need this in the gradient
    Tkl = Tkl/invAkl(1, 1)
    Tkl = Tkl+THREE*trT
    Tkl = Tkl*TWO*Skl

    !********************************************************
    ! gradient terms dTk dTl
    !********************************************************
    IF (gradflag > 0) THEN
      tk = invAkl(1, 1)  ! drop this here since we don't need tk anymore
      ddTk = mk*(JlMlk+TlMlJk)-ml*(JkMlk+TlMkJk)
      ddTk = TWO*ddTk/tk
      ddTk = ddTk+TWO*Ttmp/(tk*tk)*Jklk

      ddTk = ddTk+SIX*Tllk
      ddTk = TWO*Skl*ddTk
      ddTk = ddTk+Tkl*ddSk_over_Skl
      ! Load it into dTk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dTk(indx) = ddTk(i, j)
        ENDDO
      ENDDO
    ENDIF

    IF (gradflag > 1) THEN
      ! dTl
      Temp = ml*(JkMkl+TkMkJl)-mk*(JlMkl+TkMlJl)
      Temp = TWO*(Temp/tk+Ttmp/(tk*tk)*Jkll)

      Temp = Temp+SIX*Tkkl
      Temp = TWO*Skl*Temp
      ! Multiply by SymMatrix
      ddTl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ddTl(i, j) = ddTl(i, j)+SymMatrix(i, k)*Temp(k, j)
          ENDDO
        ENDDO
      ENDDO
      ddTl = ddTl+Tkl*ddSl_over_Skl
      ! load it into dTl (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dTl(indx) = ddTl(i, j)
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
    ! and place them in R, then sum over R*ChargeMatrix to
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

      ! da is independent of i and j
      da = -TWO*Jklk

      ! diagonal 1/rij terms
      ddVk = ZERO
      DO i = 1, n
        ! Build the matrices need for this component
        ! RLk= nvAkl*Jij*invAkl*Lk
        ! RJLk=invAkl*Jij*invAklJ11*invAkl*Lk
        DO jj = 1, n
          RLk(1:n, jj) = invAkl(1:n, i)*ALk(i, jj)
          RJLk(1:n, jj) = invAkl(1:n, i)*Jklk(i, jj)
        ENDDO

        ! JRLk=invAkl*J11*invAkl*Jij*invAkl*Lk
        DO jj = 1, n
          JRLk(1:n, jj) = invAkl(1:n, 1)*RLk(1, jj)
        ENDDO

        db = -TWO*RLk
        dc = -TWO*(JRLk+RJLk)

        ! Now we have the pices to build dRk
        coabij = c(i, i)/(a*b(i, i))
        SUM0 = ZERO
        g1 = ONE-coabij
        g2 = ONE
        DO ii = 1, p
          SUM0 = SUM0+Glob_gm3(ii)*ii*g2
          g2 = g2*g1
        ENDDO
        SUM0 = SUM0/SQRT(b(i, i))
        dRk = SUM0*(coabij*(da*b(i, i)+a*db)-dc)/(a*b(i, i))
        dRk = dRk*Skl*Glob_gm2(p)-R(i, i)/(TWO*b(i, i))*db+R(i, i)*ddSk_over_Skl
        ddVk = ddVk+dRk*Glob_ChargeMatrix(i, i)
      ENDDO

      ! Now do the off diagonal terms
      DO j = 1, n-1
        DO i = j+1, n
          ! Build the matrices need for this component
          ! RLk=invAkl*Jij*invAkl*Lk
          ! RJLk=invAkl*Jij*invAklJ11*invAkl*Lk
          DO jj = 1, n
            RLk(1:n, jj) = (invAkl(1:n, i)-invAkl(1:n, j))*(ALk(i, jj)-ALk(j, jj))
            RJLk(1:n, jj) = (invAkl(1:n, i)-invAkl(1:n, j))*(Jklk(i, jj)-Jklk(j, jj))
          ENDDO
          ! JRLk=invAkl*J11*invAkl*Jij*invAkl*Lk
          DO jj = 1, n
            JRLk(1:n, jj) = invAkl(1:n, 1)*RLk(1, jj)
          ENDDO
          db = -TWO*RLk
          dc = -TWO*(JRLk+RJLk)
          ! Now we have the pices to build dRk
          coabij = c(i, j)/(a*b(i, j))
          SUM0 = ZERO
          g1 = ONE-coabij
          g2 = ONE
          DO ii = 1, p
            SUM0 = SUM0+Glob_gm3(ii)*ii*g2
            g2 = g2*g1
          ENDDO
          SUM0 = SUM0/SQRT(b(i, j))
          dRk = SUM0*(coabij *(da*b(i, j)+a*db)-dc)/(a*b(i, j))
          dRk = dRk*Skl*Glob_gm2(p)-R(i, j)/(TWO*b(i, j))*db+R(i, j)*ddSk_over_Skl
          ddVk = ddVk+dRk*Glob_ChargeMatrix(i, j)
        ENDDO
      ENDDO
      ! Load it into dVk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dVk(indx) = ddVk(i, j)
        ENDDO
      ENDDO

    ENDIF

    IF (gradflag > 1) THEN

      !********
      ! dVl
      !********

      ! da is independent of i and j
      da = -TWO*Jkll

      ! diagonal 1/rij terms
      ddVl = ZERO
      DO i = 1, n
        ! Build the matrices need for this component
        ! RL = invAkl*Jij*invAkl*Ll
        ! RJLl = invAkl*Jij*invAklJ11*invAkl*Ll
        DO jj = 1, n
          RLl(1:n, jj) = invAkl(1:n, i)*ALl(i, jj)
          RJLl(1:n, jj) = invAkl(1:n, i)*Jkll(i, jj)
        ENDDO

        ! JRLl=invAkl*J11*invAkl*Jij*invAkl*Ll
        DO jj = 1, n
          JRLl(1:n, jj) = invAkl(1:n, 1)*RLl(1, jj)
        ENDDO

        db = -TWO*RLl
        dc = -TWO*(JRLl+RJLl)

        ! now we have the pices to build dRl
        coabij = c(i, i)/(a*b(i, i))
        SUM0 = ZERO
        g1 = ONE-coabij
        g2 = ONE
        DO ii = 1, p
          SUM0 = SUM0+Glob_gm3(ii)*ii*g2
          g2 = g2*g1
        ENDDO
        SUM0 = SUM0/SQRT(b(i, i))
        dRl = SUM0*(coabij*(da*b(i, i)+a*db)-dc)/(a*b(i, i))

        dRl = dRl*Skl*Glob_gm2(p)-R(i, i)/(TWO*b(i, i))*db

        ! Multiply by SymMatrix
        Temp = dRl
        dRl = ZERO
        DO kk = 1, n
          DO jj = 1, n
            DO ii = 1, n
              dRl(ii, jj) = dRl(ii, jj)+SymMatrix(ii, kk)*Temp(kk, jj)
            ENDDO
          ENDDO
        ENDDO
        dRl = dRl+R(i, i)*ddSl_over_Skl
        ddVl = ddVl+dRl*Glob_ChargeMatrix(i, i)

      ENDDO

      ! Now do the off diagonal terms
      DO j = 1, n-1
        DO i = j+1, n
          ! Build the matrices need for this component
          ! RLl = invAkl*Jij*invAkl*Ll
          ! RJLl = invAkl*Jij*invAklJ11*invAkl*Ll
          DO jj = 1, n
            RLl(1:n, jj) = (invAkl(1:n, i)-invAkl(1:n, j))*(ALl(i, jj)-ALl(j, jj))
            RJLl(1:n, jj) = (invAkl(1:n, i)-invAkl(1:n, j))*(Jkll(i, jj)-Jkll(j, jj))
          ENDDO
          ! JRLl=invAkl*J11*invAkl*Jij*invAkl*Ll
          DO jj = 1, n
            JRLl(1:n, jj) = invAkl(1:n, 1)*RLl(1, jj)
          ENDDO
          db = -TWO*RLl
          dc = -TWO*(JRLl+RJLl)

          ! Now we have the pices to build dRl
          coabij = c(i, j)/(a*b(i, j))
          SUM0 = ZERO
          g1 = ONE-coabij
          g2 = ONE
          DO ii = 1, p
            SUM0 = SUM0+Glob_gm3(ii)*ii*g2
            g2 = g2*g1
          ENDDO
          SUM0 = SUM0/SQRT(b(i, j))
          dRl = SUM0*(coabij*(da*b(i, j)+a*db)-dc)/(a*b(i, j))
          dRl = dRl*Skl*Glob_gm2(p)-R(i, j)/(TWO*b(i, j))*db

          ! Multiply by SymMatrix
          Temp = dRl
          dRl = ZERO
          DO kk = 1, n
            DO jj = 1, n
              DO ii = 1, n
                dRl(ii, jj) = dRl(ii, jj) + SymMatrix(ii, kk) * Temp(kk, jj)
              ENDDO
            ENDDO
          ENDDO
          dRl = dRl+R(i, j)*ddSl_over_Skl
          ddVl = ddVl+dRl*Glob_ChargeMatrix(i, j)

        ENDDO
      ENDDO
      ! load it into dVk (a vector)
      indx = 0
      DO j = 1, n
        DO i = j, n
          indx = indx+1
          dVl(indx) = ddVl(i, j)
        ENDDO
      ENDDO

    ENDIF

    Hkl = Tkl+Vkl

    SELECT CASE (gradflag)
    CASE (0)
      dHk = ZERO
      dHl = ZERO
      dSk = ZERO
      dSl = ZERO
    CASE (1)
      dHk = dTk+dVk
      dHl = ZERO
      dSl = ZERO
    CASE (2)
      dHk = dTk+dVk
      dHl = dTl+dVl
    END SELECT


  END SUBROUTINE MatrixElements0





  SUBROUTINE MatElemTraMom(mk, vechLk, ml, vechLl, SymMatrix, Skl)
    ! Overlap-type kernel used for the purely rotational transition moments
    ! N = 0 -> N = 1,  <psi(N=0,v)|x_1|psi(N=1,v)>, in a basis of correlated
    ! Gaussians with r_1 premultipliers. Kept from the original code; this
    ! program does not call it.
    ! The symmetry adaption is applied to the ket:
    ! Al --> SymMatrix'*Ll*Ll'*SymMatrix.
    ! Input:
    !   mk, ml         :: premultiplier powers
    !   vechLk         :: n(n+1)/2 exponent parameters of the N = 0 function
    !   vechLl         :: n(n+1)/2 exponent parameters of the N = 1 function
    !   SymMatrix      :: the symmetry permutation matrix
    ! Output:
    !   Skl            :: overlap matrix element of the N = 1 type
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ml
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: Skl

    ! Local variables
    INTEGER :: i, j, k, m, q, t, ii, jj, kk, indx, mko2, mlo2
    REAL(wp) :: mkml, mkpml, temp1, temp2, temp3, temp4, temp5, temp6
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3, chi1
    REAL(wp) :: detLk, detLl, detAkl, a, g1, g2, g3, g4
    REAL(wp) :: Lk(n, n), Ll(n, n), Ak(n, n), Al(n, n)
    REAL(wp) :: invAkl(n, n), PLl(n, n), invL(n, n), Akl(n, n)
    REAL(wp) :: invLk(n, n), invLl(n, n), invAk(n, n), invAl(n, n), chAl(n, n)


    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Build the lower triangular matrices Lk and Ll from
    ! the elements in vechLk and vechLl
    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
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

    ! Build AK=LK*LL' and AL=P'*LL*LL'*P
    Ak = ZERO
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Ak(i, j) = Ak(i, j)+Lk(i, k)*Lk(j, k)
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
      DO i = 1, n
        Akl(i, j) = Ak(i, j)+Al(i, j)
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
    ! need to do the decomp on Al too to get invAl since
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
    ! Now do invAkl=invL'*invL
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

    !********************************************************
    ! Compute Skl
    !
    ! Skl = 2 * (3+mk+ml)/sqrt(9+6*(mk+ml)+4*mk*ml) *
    !      gm1(mk,ml) 2^(3n/2) *
    !      sqrt( (invAkl(1,1)/invAk(1,1))^(mk+1)  *
    !            (invAkl(1,1)/invAl(1,1))^(ml+2)  *
    !            (||Lk||||Ll||/|Akl|)^3            )
    !********************************************************
    mko2 = mk/2
    mlo2 = ml/2


    ! gj1(mko2,mlo2)*gj1(mko2,mlo2)

    chi1 = ZERO
    chi1 = chi1+(NINE+SIX*(mk+ml)+FOUR*mk*ml)
    chi1 = ((THREE+mk+ml)/SQRT(chi1))
    Skl = ZERO
    Skl = Skl+SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))
    Skl = Skl*((invAkl(1, 1)/invAk(1, 1))**(mko2))*((invAkl(1, 1)/invAl(1, 1))**(mlo2))
    Skl = Skl*Glob_gm1(mko2, mlo2)*Glob_gm1(mko2, mlo2)
    Skl = Skl*SQRT((invAkl(1, 1)/invAk(1, 1)))*SQRT((invAkl(1, 1)/invAl(1, 1)))
    Skl = Skl*TWO*chi1


    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! this is a precomputed constant stored in Glob_gm1 (module globvars)
    !
    ! gj1(mk,ml) = (3+mk+ml)/sqrt(9+6*(mk+ml)+4*mk*ml)

  END SUBROUTINE MatElemTraMom





  SUBROUTINE TMMatElem(mk, vechLk, ml, vechLl, SymMatrix, TM1, TM2)
    ! Kernels of the N = 0 -> N = 1 transition moments between an N = 0
    ! function r_1^mk exp[...] (bra) and an N = 1 function (ket):
    !   TM1 = <r_1^mk ...|x_1|x_1 r_1^ml ...>  and  TM2, the x_2 counterpart
    ! (unnormalized). Kept from the original code; this program does not
    ! call it.
    ! The symmetry adaption is applied to the ket:
    ! Al --> SymMatrix'*Ll*Ll'*SymMatrix.
    ! Input:
    !   mk, ml         :: premultiplier powers of the bra and the ket
    !   vechLk, vechLl :: n(n+1)/2 exponent parameters
    !   SymMatrix      :: the symmetry permutation matrix
    ! Output:
    !   TM1, TM2       :: the two transition-moment kernels
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ml
    REAL(wp), INTENT(IN) :: vechLk(Glob_npt), vechLl(Glob_npt)
    REAL(wp), INTENT(IN) :: SymMatrix(n, n)
    REAL(wp), INTENT(OUT) :: TM1, TM2

    ! Local variables
    INTEGER :: i, j, k, indx, info, jj, ii, kk, mko2, mlo2
    REAL(wp) :: detLk, detLl, detAkl, Skl
    REAL(wp) :: SUM0, SUM1, SUM2, SUM3
    REAL(wp) :: Lk(n, n), Ll(n, n), Ak(n, n), Al(n, n), invAkl(n, n), invL(n, n), PLl(n, n)
    INTEGER :: p, q, t
    REAL(wp) :: invLk(n, n), invLl(n, n), invAk(n, n), invAl(n, n), chAl(n, n)
    REAL(wp) :: g1, g2, chi1, g3
    REAL(wp) :: a, coabij, r1kl
    REAL(wp) :: Temp(n, n), R(n, n), b(n, n), c(n, n), coab(n, n), Vtmp(n, n)

    !***********************************************************
    ! Build intermediate results for Overlap matrix element
    ! Ak, Al, invAkl, detLk, detLl, detAkl, invAk, invAl
    !***********************************************************

    ! Build the lower triangular matrices Lk and Ll from
    ! the elements in vechLk and vechLl
    Lk = ZERO
    Ll = ZERO
    indx = 0
    DO j = 1, n
      DO i = j, n
        indx = indx+1
        Lk(i, j) = vechLk(indx)
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
    ! Build AK=LK*LL' and AL=P'*LL*LL'*P
    Ak = ZERO
    Al = ZERO
    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          Ak(i, j) = Ak(i, j)+Lk(i, k)*Lk(j, k)
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
    ! need to do the decomp on Al too to get invAl since
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
    ! Now do invAkl=invL'*invL
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
    ! Compute the transition-moment kernels TM1 and TM2
    !
    ! gm1(mk,ml) = gamma((mk+ml+3)/2) * 2^((mk+ml)/2)/
    !             sqrt( gamma(mk+3/2) * gamma(ml+3/2) )
    ! is a precomputed constant stored in Glob_gm1 (globvars)
    !********************************************************
    mko2 = mk/2
    mlo2 = ml/2

    TM1 = ZERO
    TM1 = SQRT((ONE*(2**(3*n)))*(ABS(detLk*detLl/detAkl)**3))
    TM1 = TM1*Glob_gm1(mko2, mlo2)*Glob_gm1(mko2, mlo2)
    TM1 = TM1*((invAkl(1, 1)/invAk(1, 1))**(mko2))*((invAkl(1, 1)/invAl(1, 1))**(mlo2))
    TM1 = TM1*((invAkl(1, 1)/SQRT(invAl(1, 1))))
    TM1 = TM1*(THREE+mk+ml)/SQRT(THREE*(THREE+TWO*ml))

    TM2 = TM1*invAkl(1, 2)/invAkl(1, 1)

  END SUBROUTINE TMMatElem


END MODULE matelem
