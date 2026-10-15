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
  !   ExpValuesMatElem - the nucleus-nucleus correlation function and the
  !                      interparticle distances <r_ij>, <r_ij^2>
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

CONTAINS


  SUBROUTINE MatrixElements(mk, ik, vechLk, ml, il, vechLl, SymMatrix, &
                            Skl, Hkl, dSk, dSl, dHk, dHl, gradflag)
    ! Symmetry adapted matrix elements in the basis of correlated Gaussians
    ! with x_i r_1^m premultipliers for N = 1 rotationally excited states,
    !   f_k = x_i * r_1^mk * exp[-r'(Lk*Lk' (x) I3) r].
    ! The symmetry adaption is applied to the ket (exponent only):
    ! Al --> SymMatrix'*Ll*Ll'*SymMatrix.
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
    REAL(wp) :: ddSk_over_Skl(n, n), ddSl_over_Skl(n, n), Templ(n, n)
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
    partj = il

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

      Templ = &
          (-TWO/(THREE*A11*A11)*(ml*A1j*Glob_MassMatrix(1, parti)+mk*A1i*Glob_MassMatrix(1, partj))*dlA11+&
          TWO/(THREE*A11)*(ml*Glob_MassMatrix(1, parti)*dlA1j+mk*Glob_MassMatrix(1, partj)*dlA1i))

      Templ = Templ-TWO*( dlAAlMji+dlAAkMij-mkpml/(THREE*A11*A11)*(A1j*AAlM1i+A1i*AAkM1j)*dlA11+&
          mkpml/(THREE*A11)*(dlA1j*AAlM1i+dlA1i*AAkM1j+A1j*dlAAlM1i+A1i*dlAAkM1j) )
      Templ = Templ+mkml*Glob_MassMatrix(1, 1)*&
          ( TWO*(dlAij - (mkpml-TWO)/(THREE*A11*A11)*A1i*A1j*dlA11 +&
          (mkpml-TWO)/(THREE*A11)*(A1i*dlA1j+dlA1i*A1j))/(A11*(mkpml+ONE)) - &
          TWO*(Aij+(mkpml-TWO)/THREE*A1i*A1j/A11)/(A11*A11*(mkpml+ONE))*dlA11 )
      Templ = Templ-TWO*( -dlA11*TWO/(THREE*A11*A11*A11)*&
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
      Templ = Templ+FOUR*( -dlA11/(THREE*A11*A11*A11)*(&
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

      Templ = Templ*Skl/(Aij+mkpml*A1i*A1j/(THREE*A11))
      Templ = Templ+(ddSl_over_Skl-(dlAij+mkpml/(THREE*A11)*(dlA1i*A1j+&
          A1i*dlA1j-A1i*A1j*dlA11/A11))/(Aij+mkpml/THREE*A1i*A1j/A11))*Tkl

      ddTl = ZERO
      DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            ddTl(i, j) = ddTl(i, j)+SymMatrix(i, k)*Templ(k, j)
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

    DO k = 2, n
      DO j = (k-1), n
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

    g4 = Skl*Glob_gm2(m)/(THREE*invAkl(parti, partj)+TWO*m*hoa)
    DO k = 2, n
      DO j = (k-1), n
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

          dkfob = fob(p, q)*((dkA(parti, p, :, :)-dkA(parti, q, :, :))/&
              (invAkl(parti, p)-invAkl(parti, q))+&
              (dkA(partj, p, :, :)-dkA(partj, q, :, :))/&
              (invAkl(partj, p)-invAkl(partj, q))-&
              dkb/b(p, q))

          dkcoab = ((invAkl(1, p)-invAkl(1, q))*(invAkl(1, p)-invAkl(1, q))/&
              invAkl(1, 1)/b(p, q))*(dkA(1, 1, :, :)/invAkl(1, 1)+&
              dkb/b(p, q)-TWO*(dkA(1, p, :, :)-dkA(1, q, :, :))/&
              (invAkl(1, p)-invAkl(1, q)))

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

          dlfob = fob(p, q)*((dlA(parti, p, :, :)-dlA(parti, q, :, :))/&
              (invAkl(parti, p)-invAkl(parti, q))+&
              (dlA(partj, p, :, :)-dlA(partj, q, :, :))/&
              (invAkl(partj, p)-invAkl(partj, q))-&
              dlb/b(p, q))

          dlcoab = ((invAkl(1, p)-invAkl(1, q))*(invAkl(1, p)-invAkl(1, q))/&
              invAkl(1, 1)/b(p, q))*(dlA(1, 1, :, :)/invAkl(1, 1)+&
              dlb/b(p, q)-TWO*(dlA(1, p, :, :)-dlA(1, q, :, :))/&
              (invAkl(1, p)-invAkl(1, q)))

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


  END SUBROUTINE MatrixElements





  SUBROUTINE ExpValuesMatElem(mk, vechLk, ml, vechLl, SymMatrix, &
                              NumCFGridPoints, CFGridX, ycf, CFkl, AreCorrFuncNeeded, &
                              ArePairDistNeeded, rkl, rkl_2)
    ! Matrix elements of the expectation-value operators between two basis
    ! functions: the nucleus-nucleus correlation function on one row of the
    ! (x,y) grid, and the interparticle distances <r_ij> and <r_ij^2>.
    ! The formulas are those of the x_1 r_1^m basis (Sharkey, Kirnosov and
    ! Adamowicz, J. Chem. Phys. 139, 164119 (2013)) with the powers mk and
    ! ml of the two functions, exactly as in the original code; they do not
    ! know about the x_i index of the prefactor. The caller passes the
    ! powers 2*(code/n) and is expected to warn when the basis holds
    ! functions with i > 1.
    ! Input:
    !   mk, ml           :: premultiplier powers of the bra and the ket
    !   vechLk, vechLl   :: n(n+1)/2 exponent parameters
    !   SymMatrix        :: the symmetry permutation matrix (ket)
    !   NumCFGridPoints  :: number of x grid points of the correlation function
    !   CFGridX          :: the x grid points
    !   ycf              :: the y grid point of this call
    !   AreCorrFuncNeeded  :: compute the correlation function at all
    !   ArePairDistNeeded  :: compute rkl and rkl_2 (they do not depend on
    !                         the grid, so once per pair is enough)
    ! Output (unnormalized):
    !   CFkl(1:NumCFGridPoints) :: correlation function at (CFGridX(i), ycf)
    !   rkl, rkl_2       :: <r_ij> and <r_ij^2>; the diagonal (i,i) is the
    !                       distance of pseudoparticle i from particle 0
    !================================================================
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: n = Glob_n  ! number of pseudoparticles

    INTEGER, INTENT(IN) :: mk, ml
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

    ! Computing nucleus-nucleus correlation function if needed

    p = (mk+ml)/2

    IF (AreCorrFuncNeeded) THEN
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

    ENDIF

  END SUBROUTINE ExpValuesMatElem





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
