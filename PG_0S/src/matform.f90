MODULE matform
  ! Module matform forms the Hamiltonian and overlap matrices (and their
  ! derivatives) from the matrix elements of module matelem.
  !
  ! one call returns all Glob_NumYHYTerms symmetry
  ! terms of a pair, so the work is dealt out per PAIR (PG_0S dealt it per
  ! pair and term). Buffers, storage routines and the layout of Glob_D are
  ! those of PG_0S. 
  !
  ! The inverse-iteration shift helpers ReportWrongStateOnce,
  ! RetargetShiftToEigenvalue, RefreshINVITShift, IsEigenpairUsable and
  ! IsRequestedEigenstate (called by the EnergyI* routines and the I
  ! drivers of workproc) close the module.
  USE globvars
  USE matelem
  USE misc
  USE linalg

  IMPLICIT NONE

CONTAINS


  SUBROUTINE CheckSelfOverlap(i, j, Sij)
    ! Reports a self-overlap that is not strictly positive. It can only
    ! come out zero or negative through catastrophic cancellation or a bug,
    ! and everything downstream divides by it. Called after the
    ! MPI_ALLREDUCE, so every rank sees the same value and the same verdict.

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: i, j
    REAL(wp), INTENT(IN) :: Sij

    IF (i == j) THEN
      ! Sij IS the self-overlap here, and it is about to become the divisor
      ! and to be stored in Glob_diagS(i) for every later element of row i.
      IF (Sij > ZERO) RETURN
      Glob_Warning = 2
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error in matform: non-positive self-overlap'
        WRITE(*, *) '  function index =', i, '  S_ii =', Sij
      ENDIF
    ELSE
      IF ((Glob_diagS(i) > ZERO) .AND. (Glob_diagS(j) > ZERO)) RETURN
      Glob_Warning = 2
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *) 'Error in matform: non-positive self-overlap'
        WRITE(*, *) '  i =', i, ' Glob_diagS(i) =', Glob_diagS(i)
        WRITE(*, *) '  j =', j, ' Glob_diagS(j) =', Glob_diagS(j)
      ENDIF
    ENDIF

  END SUBROUTINE CheckSelfOverlap


  FUNCTION NormFactor(i, j)
    ! 1/SQRT(S_ii*S_jj), the normalization factor of element (i,j).
    !
    ! The product is formed first; when it leaves the exponent range (two
    ! healthy self-overlaps beyond about 1e154 or 1e-154 at PREC=8) the
    ! result is 0 or Inf, silently zeroing a row or filling it with NaN.
    ! ONE/(SQRT(a)*SQRT(b)) halves both exponents before multiplying and
    ! stays in range whenever the answer does. It is used only as a
    ! FALLBACK, after the direct form has been seen to fail, because the
    ! two forms differ by an ulp or so and using it everywhere would shift
    ! every energy in its last digits. The test also catches NaN.

    IMPLICIT NONE

    REAL(wp)            :: NormFactor
    INTEGER, INTENT(IN) :: i, j

    NormFactor = 1/SQRT(Glob_diagS(i)*Glob_diagS(j))
    IF (.NOT. ((NormFactor > ZERO) .AND. (NormFactor <= HUGE(NormFactor)))) &
      NormFactor = ONE/(SQRT(Glob_diagS(i))*SQRT(Glob_diagS(j)))

  END FUNCTION NormFactor


  SUBROUTINE StoreHS(i, j, Hij, Sij)
    ! Routine StoreHS stores calculated matrix elements of the
    ! Hamiltonian and the overlap in proper places of global
    ! arrays. Upon doing this, the routine normalizes
    ! matrix elements.
    ! Important comment:  i must be greater or equal to j.

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: i, j
    REAL(wp), INTENT(IN) :: Hij, Sij
    REAL(wp)             :: f

    CALL CheckSelfOverlap(i, j, Sij)

    SELECT CASE (Glob_GSEPSolutionMethod)
    CASE ('I')
      ! Only the lower triangle of array Glob_H (including the diagonal) is
      ! used to store H-Glob_ApproxEnergy*S.
      ! The entire array Glob_S is used to store S.
      IF (i == j) THEN
        Glob_diagS(i) = Sij
        Glob_S(i, i) = ONE
        Glob_H(i, i) = Hij/Glob_diagS(i)-Glob_ApproxEnergy
      ELSE
        f = NormFactor(i, j)
        Glob_S(i, j) = Sij*f
        Glob_S(j, i) = Glob_S(i, j)
        Glob_H(i, j) = (Hij-Glob_ApproxEnergy*Sij)*f
      ENDIF
    CASE ('G')
      ! In the case when Glob_GSEPSolutionMethod='G' the diagonal
      ! of the Hamiltonian matrix is stored in Glob_diagH
      ! The diagonal of the overlap is stored in Glob_diagS
      ! The lower triangles of arrays Glob_H and
      ! Glob_S are used to store H and S
      IF (i == j) THEN
        Glob_diagS(i) = Sij
        Glob_diagH(i) = Hij/Glob_diagS(i)
      ELSE
        f = NormFactor(i, j)
        Glob_S(i, j) = Sij*f
        Glob_H(i, j) = Hij*f
      ENDIF
    CASE ('Q')
      ! The QR method keeps the normalized, UNSHIFTED H and S in their
      ! canonical order: the complete lower triangles, both diagonals
      ! included, are authoritative and the upper triangles are not
      ! written (qrlinalg reads only the lower triangles). Glob_diagS
      ! keeps the raw self-overlap for the normalization of the later
      ! elements of row i and of the derivatives.
      IF (i == j) THEN
        Glob_diagS(i) = Sij
        Glob_S(i, i) = ONE
        Glob_H(i, i) = Hij/Glob_diagS(i)
      ELSE
        f = NormFactor(i, j)
        Glob_S(i, j) = Sij*f
        Glob_H(i, j) = Hij*f
      ENDIF
    END SELECT

  END SUBROUTINE StoreHS


  SUBROUTINE StoreHSD(i, j, Hij, Sij, Di, Dj)
    ! Routine StoreHSD stores the normalized matrix elements of H and S
    ! and their derivatives in the global arrays. Requires i >= j. Di and
    ! Dj are the derivatives of Hij and Sij with respect to the nonlinear
    ! parameters of functions i and j:
    !   Di(1:Glob_npt)             dHij/dvechLi
    !   Di(Glob_npt+1:2*Glob_npt)  dSij/dvechLi
    !   Dj(1:Glob_npt)             dHij/dvechLj
    !   Dj(Glob_npt+1:2*Glob_npt)  dSij/dvechLj

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: i, j
    REAL(wp), INTENT(IN) :: Hij, Sij
    REAL(wp), INTENT(IN) :: Di(2*Glob_npt), Dj(2*Glob_npt)
    REAL(wp)             :: f

    CALL CheckSelfOverlap(i, j, Sij)

    SELECT CASE (Glob_GSEPSolutionMethod)
    CASE ('I')
      ! Only the lower triangle of array Glob_H (including the diagonal) is
      ! used to store H-Glob_ApproxEnergy*S.
      ! The entire array Glob_S is used to store S.
      IF (i == j) THEN
        Glob_diagS(i) = Sij
        Glob_S(i, i) = ONE
        Glob_H(i, i) = Hij/Glob_diagS(i)-Glob_ApproxEnergy
        Glob_D(1:2*Glob_npt, i-Glob_nfru, i) = TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
      ELSE
        f = NormFactor(i, j)
        Glob_S(i, j) = Sij*f
        Glob_S(j, i) = Glob_S(i, j)
        Glob_H(i, j) = (Hij-Glob_ApproxEnergy*Sij)*f
        Glob_D(1:2*Glob_npt, i-Glob_nfru, j) = Di(1:2*Glob_npt)*f
        IF (j > Glob_nfru) Glob_D(1:2*Glob_npt, j-Glob_nfru, i) = Dj(1:2*Glob_npt)*f
      ENDIF
    CASE ('G')
      ! In the case when Glob_GSEPSolutionMethod='G' the diagonal
      ! of the Hamiltonian matrix is stored in Glob_diagH
      ! The diagonal of the overlap is stored in Glob_diagS
      ! Lower triangles of arrays Glob_H and
      ! Glob_S are used to store H and S
      IF (i == j) THEN
        Glob_diagS(i) = Sij
        Glob_diagH(i) = Hij/Glob_diagS(i)
        Glob_D(1:2*Glob_npt, i-Glob_nfru, i) = TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
      ELSE
        f = NormFactor(i, j)
        Glob_S(i, j) = Sij*f
        Glob_H(i, j) = Hij*f
        Glob_D(1:2*Glob_npt, i-Glob_nfru, j) = Di(1:2*Glob_npt)*f
        IF (j > Glob_nfru) Glob_D(1:2*Glob_npt, j-Glob_nfru, i) = Dj(1:2*Glob_npt)*f
      ENDIF
    END SELECT

  END SUBROUTINE StoreHSD


  SUBROUTINE ComputeMatElem(Nmin, Nmax)
    ! Subroutine ComputeMatElem computes the matrix elements of H and S
    ! (no derivatives) for the functions Nmin..Nmax and stores them with
    ! StoreHS. The elements of the first Nmin-1 functions must already be
    ! stored; set Nmin=1 to compute everything.
    !  Input parameters :
    !   Nmin-1 :: number of functions whose matrix elements are known
    !     Nmax :: number of functions whose matrix elements are needed
    ! Each pair is evaluated by ONE rank (pair i goes to rank
    ! MOD(i,Glob_NumOfProcs)); the other ranks leave zeros in the buffer
    ! and the MPI_ALLREDUCE assembles the sum.

    IMPLICIT NONE

    ! Arguments :
    INTEGER, INTENT(IN) :: Nmin, Nmax
    ! Local variables :
    INTEGER  :: k, l, i, kk, ll, ii, j
    INTEGER  :: kstart, lstart, lstop, n, npt, nb
    INTEGER  :: mk, ml
    REAL(wp) :: Paraml(Glob_n*(Glob_n+1)/2)
    ! Bra side of the pairs (k,l): Lk, Ak = Lk*Lk' and M*Ak, built once per k
    REAL(wp) :: PreLk(Glob_n, Glob_n, 1), PreAk(Glob_n, Glob_n, 1), PreMAk(Glob_n, Glob_n, 1)
    REAL(wp) :: Ssum, Hsum
    ! Work arrays of MatrixElementsOpt, one entry per symmetry term. The
    ! derivative buffers are not used here but must be present in the call.
    REAL(wp) :: SymMatrixBuf(Glob_NumYHYTerms, Glob_n, Glob_n)
    REAL(wp) :: SklBuf(Glob_NumYHYTerms), HklBuf(Glob_NumYHYTerms)
    REAL(wp) :: dSkBuf(Glob_NumYHYTerms, Glob_npt), dSlBuf(Glob_NumYHYTerms, Glob_npt)
    REAL(wp) :: dHkBuf(Glob_NumYHYTerms, Glob_npt), dHlBuf(Glob_NumYHYTerms, Glob_npt)

    n = Glob_n
    npt = Glob_npt
    nb = Glob_HSBuffLen
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,i0,a)') 'ComputeMatElem: functions ', Nmin, '..', Nmax, &
                                          ' (', Nmax*(Nmax+1)/2-(Nmin-1)*Nmin/2, ' pairs)'
    ENDIF
    ! SymMatrixBuf is loop-invariant and MatrixElementsOpt does not modify
    ! it, so it is built once here rather than once per pair
    DO j = 1, Glob_NumYHYTerms
      SymMatrixBuf(j, :, :) = Glob_YHYMatr(1:n, 1:n, j)
    ENDDO
    Glob_HklBuff1(1:nb) = ZERO
    Glob_SklBuff1(1:nb) = ZERO
    i = 0

    DO k = Nmin, Nmax
      CALL PrecomputeMatrices_MA_L_A(npt, 1, Glob_NonlinParam(1:npt, k:k), &
                                     Glob_MassMatrix(1:n, 1:n), PreLk, PreAk, PreMAk)
      mk = Glob_PWR(k)
      DO l = k, 1, -1
        i = i+1
        IF (i == 1) THEN
          kstart = k
          lstart = l
        ENDIF
        Paraml(1:npt) = Glob_NonlinParam(1:npt, l)
        ml = Glob_PWR(l)
        Hsum = ZERO; Ssum = ZERO
        IF (MOD(i, Glob_NumOfProcs) == Glob_ProcID) THEN
          CALL MatrixElementsOpt(mk, PreLk(:, :, 1), PreAk(:, :, 1), PreMAk(:, :, 1), &
                                 ml, Paraml, SymMatrixBuf, SklBuf, HklBuf, &
                                 dSkBuf, dSlBuf, dHkBuf, dHlBuf, 0)
          DO j = 1, Glob_NumYHYTerms
            Hsum = Hsum+Glob_YHYCoeff(j)*HklBuf(j)
            Ssum = Ssum+Glob_YHYCoeff(j)*SklBuf(j)
          ENDDO
        ENDIF
        Glob_HklBuff1(i) = Hsum
        Glob_SklBuff1(i) = Ssum
        IF (i == Glob_HSBuffLen) THEN
          CALL MPI_ALLREDUCE(Glob_HklBuff1, Glob_HklBuff2, i, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_ALLREDUCE(Glob_SklBuff1, Glob_SklBuff2, i, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          ii = 0
          DO kk = kstart, k
            IF (kk == kstart) THEN
              ll = lstart
            ELSE
              ll = kk
            ENDIF
            IF (kk == k) THEN
              lstop = l
            ELSE
              lstop = 1
            ENDIF
            DO WHILE (ll >= lstop)
              ii = ii+1
              CALL StoreHS(kk, ll, Glob_HklBuff2(ii), Glob_SklBuff2(ii))
              ll = ll-1
            ENDDO
          ENDDO
          i = 0
          Glob_HklBuff1(1:nb) = ZERO
          Glob_SklBuff1(1:nb) = ZERO
        ENDIF
      ENDDO
    ENDDO
    IF (i > 0) THEN
      CALL MPI_ALLREDUCE(Glob_HklBuff1, Glob_HklBuff2, i, &
                         MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_ALLREDUCE(Glob_SklBuff1, Glob_SklBuff2, i, &
                         MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      ii = 0
      DO kk = kstart, Nmax
        IF (kk == kstart) THEN
          ll = lstart
        ELSE
          ll = kk
        ENDIF
        lstop = 1
        DO WHILE (ll >= lstop)
          ii = ii+1
          CALL StoreHS(kk, ll, Glob_HklBuff2(ii), Glob_SklBuff2(ii))
          ll = ll-1
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE ComputeMatElem


  SUBROUTINE ComputeMatElemAndDeriv(Nmin, Nmax)
    ! Subroutine ComputeMatElemAndDeriv computes the matrix elements of H
    ! and S and their derivatives for the functions Nmin..Nmax and stores
    ! them with StoreHSD. The elements of the first Nmin-1 functions must
    ! already be stored (set Nmin=1 for all), and Glob_nfo and Glob_nfru
    ! must be set.
    !  Input parameters :
    !   Nmin-1 :: number of functions whose matrix elements are known
    !     Nmax :: number of functions whose matrix elements are needed
    ! Each pair is evaluated by ONE rank and assembled by MPI_ALLREDUCE.
    ! gradflag is 1 when only the derivatives with respect to function k
    ! are wanted and 2 when those with respect to l are wanted too (l is
    ! itself optimized). Dk(1:npt) = dH/dvechLk, Dk(npt+1:2*npt) = dS/dvechLk,
    ! and likewise for Dl.

    IMPLICIT NONE

    ! Arguments :
    INTEGER, INTENT(IN) :: Nmin, Nmax
    ! Local variables :
    INTEGER  :: k, l, i, kk, ll, ii, j
    INTEGER  :: kstart, lstart, lstop, n, npt, npt2, nb
    INTEGER  :: mk, ml, gradflag
    REAL(wp) :: Paraml(Glob_n*(Glob_n+1)/2)
    ! Bra side of the pairs (k,l): Lk, Ak = Lk*Lk' and M*Ak, built once per k
    REAL(wp) :: PreLk(Glob_n, Glob_n, 1), PreAk(Glob_n, Glob_n, 1), PreMAk(Glob_n, Glob_n, 1)
    REAL(wp) :: Ssum, Hsum
    REAL(wp) :: Dksum(Glob_n*(Glob_n+1))
    REAL(wp) :: Dlsum(Glob_n*(Glob_n+1))
    LOGICAL  :: grad_l
    ! Work arrays of MatrixElementsOpt, one entry per symmetry term
    REAL(wp) :: SymMatrixBuf(Glob_NumYHYTerms, Glob_n, Glob_n)
    REAL(wp) :: SklBuf(Glob_NumYHYTerms), HklBuf(Glob_NumYHYTerms)
    REAL(wp) :: dSkBuf(Glob_NumYHYTerms, Glob_npt), dSlBuf(Glob_NumYHYTerms, Glob_npt)
    REAL(wp) :: dHkBuf(Glob_NumYHYTerms, Glob_npt), dHlBuf(Glob_NumYHYTerms, Glob_npt)

    n = Glob_n
    npt = Glob_npt
    npt2 = 2*npt  ! the derivative buffers and Glob_D are dimensioned 2*Glob_npt
    nb = Glob_HSBuffLen
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,i0,a,i0)') 'ComputeMatElemAndDeriv: functions ', Nmin, '..', Nmax, &
                                            ', derivatives for ', Glob_nfo, ' function(s) after ', Glob_nfru
    ENDIF
    ! SymMatrixBuf is loop-invariant and MatrixElementsOpt does not modify
    ! it, so it is built once here rather than once per pair
    DO j = 1, Glob_NumYHYTerms
      SymMatrixBuf(j, :, :) = Glob_YHYMatr(1:n, 1:n, j)
    ENDDO

    Glob_HklBuff1(1:nb) = ZERO
    Glob_SklBuff1(1:nb) = ZERO
    Glob_DkBuff1(1:npt2, 1:nb) = ZERO
    IF (Glob_nfo > 1) Glob_DlBuff1(1:npt2, 1:nb) = ZERO
    i = 0

    DO k = Nmin, Nmax
      CALL PrecomputeMatrices_MA_L_A(npt, 1, Glob_NonlinParam(1:npt, k:k), &
                                     Glob_MassMatrix(1:n, 1:n), PreLk, PreAk, PreMAk)
      mk = Glob_PWR(k)
      DO l = k, 1, -1
        i = i+1
        IF (i == 1) THEN
          kstart = k
          lstart = l
        ENDIF
        Paraml(1:npt) = Glob_NonlinParam(1:npt, l)
        ml = Glob_PWR(l)
        Hsum = ZERO
        Ssum = ZERO
        Dksum(1:npt2) = ZERO
        IF ((l > Glob_nfru) .AND. (l /= k)) THEN
          grad_l = .TRUE.
          gradflag = 2
          Dlsum(1:npt2) = ZERO
        ELSE
          grad_l = .FALSE.
          gradflag = 1
        ENDIF
        IF (MOD(i, Glob_NumOfProcs) == Glob_ProcID) THEN
          CALL MatrixElementsOpt(mk, PreLk(:, :, 1), PreAk(:, :, 1), PreMAk(:, :, 1), &
                                 ml, Paraml, SymMatrixBuf, SklBuf, HklBuf, &
                                 dSkBuf, dSlBuf, dHkBuf, dHlBuf, gradflag)
          DO j = 1, Glob_NumYHYTerms
            Hsum = Hsum+Glob_YHYCoeff(j)*HklBuf(j)
            Ssum = Ssum+Glob_YHYCoeff(j)*SklBuf(j)
            Dksum(1:npt) = Dksum(1:npt)+Glob_YHYCoeff(j)*dHkBuf(j, 1:npt)
            Dksum(npt+1:npt2) = Dksum(npt+1:npt2)+Glob_YHYCoeff(j)*dSkBuf(j, 1:npt)
            IF (grad_l) THEN
              Dlsum(1:npt) = Dlsum(1:npt)+Glob_YHYCoeff(j)*dHlBuf(j, 1:npt)
              Dlsum(npt+1:npt2) = Dlsum(npt+1:npt2)+Glob_YHYCoeff(j)*dSlBuf(j, 1:npt)
            ENDIF
          ENDDO
        ENDIF
        Glob_HklBuff1(i) = Hsum
        Glob_SklBuff1(i) = Ssum
        Glob_DkBuff1(1:npt2, i) = Dksum(1:npt2)
        IF (grad_l) Glob_DlBuff1(1:npt2, i) = Dlsum(1:npt2)
        IF (i == Glob_HSBuffLen) THEN
          CALL MPI_ALLREDUCE(Glob_HklBuff1, Glob_HklBuff2, i, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_ALLREDUCE(Glob_SklBuff1, Glob_SklBuff2, i, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          CALL MPI_ALLREDUCE(Glob_DkBuff1, Glob_DkBuff2, i*npt2, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          IF (Glob_nfo > 1) CALL MPI_ALLREDUCE(Glob_DlBuff1, Glob_DlBuff2, i*npt2, &
                                               MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          ii = 0
          DO kk = kstart, k
            IF (kk == kstart) THEN
              ll = lstart
            ELSE
              ll = kk
            ENDIF
            IF (kk == k) THEN
              lstop = l
            ELSE
              lstop = 1
            ENDIF
            DO WHILE (ll >= lstop)
              ii = ii+1
              CALL StoreHSD(kk, ll, Glob_HklBuff2(ii), Glob_SklBuff2(ii), &
                            Glob_DkBuff2(1:npt2, ii), Glob_DlBuff2(1:npt2, ii))
              ll = ll-1
            ENDDO
          ENDDO
          i = 0
          Glob_HklBuff1(1:nb) = ZERO
          Glob_SklBuff1(1:nb) = ZERO
          Glob_DkBuff1(1:npt2, 1:nb) = ZERO
          IF (Glob_nfo > 1) Glob_DlBuff1(1:npt2, 1:nb) = ZERO
        ENDIF
      ENDDO
    ENDDO
    IF (i > 0) THEN
      CALL MPI_ALLREDUCE(Glob_HklBuff1, Glob_HklBuff2, i, &
                         MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_ALLREDUCE(Glob_SklBuff1, Glob_SklBuff2, i, &
                         MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      CALL MPI_ALLREDUCE(Glob_DkBuff1, Glob_DkBuff2, i*npt2, &
                         MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      IF (Glob_nfo > 1) CALL MPI_ALLREDUCE(Glob_DlBuff1, Glob_DlBuff2, i*npt2, &
                                           MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      ii = 0
      DO kk = kstart, Nmax
        IF (kk == kstart) THEN
          ll = lstart
        ELSE
          ll = kk
        ENDIF
        lstop = 1
        DO WHILE (ll >= lstop)
          ii = ii+1
          CALL StoreHSD(kk, ll, Glob_HklBuff2(ii), Glob_SklBuff2(ii), &
                        Glob_DkBuff2(1:npt2, ii), Glob_DlBuff2(1:npt2, ii))
          ll = ll-1
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE ComputeMatElemAndDeriv

  FUNCTION SelfOverlapCancellation(pwr, param, Ssum, Sabs) RESULT(Cfac)
    ! Cancellation factor of the self-overlap <phi|Y+Y|phi> of ONE basis
    ! function with premultiplier power pwr and nonlinear parameters param:
    !   Ssum = sum_k c_k S_k       the self-overlap that survives Y+Y
    !   Sabs = sum_k |c_k S_k|     the size of the terms that were summed
    !   Cfac = Sabs/|Ssum| >= 1    log10(Cfac) digits are lost; HUGE if Ssum = 0
    ! One self-pair call of MatrixElementsOpt, no MPI: every rank gets the
    ! same value. Used by the acceptance test of BasisEnlG/I and by SaveHSRaw.
    ! Verbose 3: two summary lines per call on rank 0 (sums, C, digits lost,
    ! positive/negative parts, largest term); Verbose 4: also the nonlinear
    ! parameters and every term c_k*S_k.
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: pwr
    REAL(wp), INTENT(IN)  :: param(Glob_npt)
    REAL(wp), INTENT(OUT) :: Ssum, Sabs
    REAL(wp)              :: Cfac
    INTEGER  :: j
    INTEGER  :: jmax        ! index of the largest |c_k S_k|
    REAL(wp) :: Term        ! c_k S_k
    REAL(wp) :: Spos, Sneg  ! sums of the positive and of the negative terms
    REAL(wp) :: Tmax        ! largest |c_k S_k|
    REAL(wp) :: SymMatrixBuf(Glob_NumYHYTerms, Glob_n, Glob_n)
    REAL(wp) :: SklBuf(Glob_NumYHYTerms), HklBuf(Glob_NumYHYTerms)
    REAL(wp) :: dSkBuf(Glob_NumYHYTerms, Glob_npt), dSlBuf(Glob_NumYHYTerms, Glob_npt)
    REAL(wp) :: dHkBuf(Glob_NumYHYTerms, Glob_npt), dHlBuf(Glob_NumYHYTerms, Glob_npt)
    ! Bra side: Lk, Ak = Lk*Lk' and M*Ak of this function
    REAL(wp) :: ParamBuf(Glob_npt, 1)
    REAL(wp) :: PreLk(Glob_n, Glob_n, 1), PreAk(Glob_n, Glob_n, 1), PreMAk(Glob_n, Glob_n, 1)

    DO j = 1, Glob_NumYHYTerms
      SymMatrixBuf(j, :, :) = Glob_YHYMatr(1:Glob_n, 1:Glob_n, j)
    ENDDO
    ParamBuf(1:Glob_npt, 1) = param(1:Glob_npt)
    CALL PrecomputeMatrices_MA_L_A(Glob_npt, 1, ParamBuf, Glob_MassMatrix(1:Glob_n, 1:Glob_n), &
                                   PreLk, PreAk, PreMAk)
    CALL MatrixElementsOpt(pwr, PreLk(:, :, 1), PreAk(:, :, 1), PreMAk(:, :, 1), &
                           pwr, param, SymMatrixBuf, SklBuf, HklBuf, &
                           dSkBuf, dSlBuf, dHkBuf, dHlBuf, 0)
    Ssum = ZERO
    Sabs = ZERO
    Spos = ZERO
    Sneg = ZERO
    Tmax = ZERO
    jmax = 0
    DO j = 1, Glob_NumYHYTerms
      Term = Glob_YHYCoeff(j)*SklBuf(j)
      Ssum = Ssum+Term
      Sabs = Sabs+ABS(Term)
      IF (Term > ZERO) THEN
        Spos = Spos+Term
      ELSE
        Sneg = Sneg+Term
      ENDIF
      IF (ABS(Term) > Tmax) THEN
        Tmax = ABS(Term)
        jmax = j
      ENDIF
    ENDDO
    IF (Ssum /= ZERO) THEN
      Cfac = Sabs/ABS(Ssum)
    ELSE
      Cfac = HUGE(Cfac)
    ENDIF
    IF (Cfac < ONE) Cfac = ONE

    IF ((Glob_ProcID == 0) .AND. (Verbose >= 3)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,es12.5,a,es12.5,a,es10.3,a,f7.1,a)') &
        'SelfOverlapCancellation: power ', pwr, ', terms ', Glob_NumYHYTerms, ': Ssum = ', Ssum, &
        '  Sabs = ', Sabs, '  C = ', Cfac, '  (', LOG10(Cfac), ' digits lost)'
      WRITE(*, '(1x,a,es12.5,a,es12.5,a,i0,a,es12.5,a,es9.2)') &
        '  positive terms ', Spos, '  negative terms ', Sneg, '  largest |term| k=', jmax, ': ', Tmax, &
        '  limit C = ', Glob_MaxSelfOverlapCancel
      IF (Cfac > Glob_MaxSelfOverlapCancel) WRITE(*, '(1x,a)') '  C is above the limit: the function is rejected'
      IF (Verbose >= 4) THEN
        WRITE(*, '(1x,a,i0,a)') '  nonlinear parameters (', Glob_npt, '):'
        WRITE(*, '(4(1x,es16.8))') param(1:Glob_npt)
        WRITE(*, '(1x,a)') '       k               c_k               S_k           c_k*S_k'
        DO j = 1, Glob_NumYHYTerms
          WRITE(*, '(1x,i7,3(1x,es17.9))') j, Glob_YHYCoeff(j), SklBuf(j), Glob_YHYCoeff(j)*SklBuf(j)
        ENDDO
      ENDIF
    ENDIF
  END FUNCTION SelfOverlapCancellation


  SUBROUTINE ReportWrongStateOnce(Nmax, ErrorCode)
    !==================================================================
    ! Subroutine ReportWrongStateOnce
    !==================================================================
    ! Says once per run when GSEPIIS converged to a state other than the
    ! one the data file asked for (Glob_LastEigIndex from the inertia count
    ! and the sign of lambda-shift). A REPORT, not a rejection: a transient
    ! mismatch inside an optimization is normal with index targeting on; a
    ! PERSISTENT one means the wrong state is being minimized and the
    ! Hylleraas-Undheim-MacDonald bound does not apply. Called by EnergyIA,
    ! EnergyIAM and EnergyIB after every solve.
    !==================================================================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nmax       ! size of the basis just solved
    INTEGER, INTENT(IN) :: ErrorCode  ! status of that solve

    IF (Glob_WrongStateReported) RETURN
    IF (Glob_ProcID /= 0) RETURN
    IF (ErrorCode /= 0) RETURN
    IF (Glob_LastEigIndex <= 0) RETURN
    IF (Nmax < Glob_WhichEigenvalue) RETURN
    IF (Glob_LastEigIndex == Glob_WhichEigenvalue) RETURN

    Glob_WrongStateReported = .TRUE.
    WRITE(*, *)
    IF (Verbose >= 1) WRITE(*, *) '*** WARNING: inverse iteration is not on the requested state ***'
    WRITE(*, *) '  WHICH_EIGENVALUE  = ', Glob_WhichEigenvalue
    WRITE(*, *) '  eigenvalue found  = ', Glob_LastEigIndex
    WRITE(*, *) '  eigenvalues below the shift = ', Glob_NumEvalsBelowShift
    IF (Glob_EigIdxTargeting == 1) THEN
      WRITE(*, *) '  The shift is re-anchored at each BBOP step; if this'
      WRITE(*, *) '  persists, CURRENT_ENERGY is far from the wanted level.'
    ELSE
      IF (Verbose >= 2) WRITE(*, *) '  Eigenvalue-index targeting is OFF (Glob_EigIdxTargeting=0).'
      WRITE(*, *) '  Inverse iteration returns the level nearest the shift,'
      WRITE(*, *) '  so it will keep optimizing this one.'
    ENDIF
    WRITE(*, *) '  This message is printed once per run.'
    WRITE(*, *)

  END SUBROUTINE ReportWrongStateOnce


  SUBROUTINE RetargetShiftToEigenvalue(N, RoutineName)
    !==================================================================
    ! Subroutine RetargetShiftToEigenvalue
    !==================================================================
    ! Moves Glob_ApproxEnergy onto eigenvalue number Glob_WhichEigenvalue so
    ! that inverse iteration converges to the state actually asked for.
    ! GSEPIIS returns the eigenvalue NEAREST the shift; for an EXCITED state
    ! the Hylleraas-Undheim-MacDonald bound holds only for the eigenvalue
    ! selected BY INDEX, so minimizing "the eigenvalue nearest a shift"
    ! slides onto lower states (variational collapse).
    ! GSEPIIS_ShiftForIndex bisects on the inertia count until exactly
    ! Glob_WhichEigenvalue-1 eigenvalues lie below sigma; Glob_H is then
    ! re-shifted in O(N^2): H - sigma_new*S = (H - sigma_old*S) +
    ! (sigma_old - sigma_new)*S (Glob_S has a unit diagonal in 'I' mode).
    ! Cost: one LDL^T per bisection step, once per BBOP step. Inert unless
    ! Glob_EigIdxTargeting == 1.
    !==================================================================

    IMPLICIT NONE

    INTEGER, INTENT(IN)      :: N            ! current basis size
    CHARACTER(*), INTENT(IN) :: RoutineName  ! caller, for the messages

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: invDwork, wwork
    REAL(wp)                            :: NewSigma, lo, hi, width, OldSigma
    INTEGER                             :: NumFact, ErrorCode

    IF (Glob_EigIdxTargeting /= 1) RETURN

    ! A basis of N functions supports only N eigenvalues; asking for
    ! number Glob_WhichEigenvalue before the basis is that large is
    ! meaningless, so leave the shift alone until it is.
    IF (N < Glob_WhichEigenvalue) RETURN

    ! Refuse to work from a poisoned shift. When an energy evaluation is
    ! rejected the caller substitutes 1e31, and if that value reaches
    ! CURRENT_ENERGY then main.f90 turns it into Glob_ApproxEnergy on the
    ! next BBOP step. Bisecting from there is meaningless - it produced a
    ! shift of -1.19e19 in testing - and would bury the real failure.
    IF (ABS(Glob_ApproxEnergy) > 1.0E10_wp) THEN
      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 2) WRITE(*, *) 'Eigenvalue-index targeting skipped in ', RoutineName
        WRITE(*, '(a)', ADVANCE='no') ' the shift is not a physical energy: '
        CALL writerealadv(6, Glob_ApproxEnergy)
        WRITE(*, *) 'An earlier energy evaluation failed and poisoned CURRENT_ENERGY.'
        WRITE(*, *)
      ENDIF
      RETURN
    ENDIF

    OldSigma = Glob_ApproxEnergy
    ALLOCATE(invDwork(N))
    ALLOCATE(wwork(N))

    ! A deliberately modest starting bracket - GSEPIIS_ShiftForIndex widens
    ! it by itself if the inertia count says the target is outside.
    width = MAX(0.1_wp, 0.1_wp*ABS(OldSigma))
    lo = OldSigma-width
    hi = OldSigma+width

    ! Glob_H is re-shifted IN PLACE: on success its lower triangle comes
    ! back holding H - NewSigma*S, and on failure it is restored to
    ! H - OldSigma*S. Either way it stays consistent with
    ! Glob_ApproxEnergy as set below, which is what the next EnergyIA
    ! assumes. Its upper triangle is left holding an L^T that the next
    ! solve overwrites.
    CALL GSEPIIS_ShiftForIndex(N, Glob_H, Glob_HSLeadDim, OldSigma, &
                               Glob_S, Glob_HSLeadDim, &
                               Glob_WhichEigenvalue, lo, hi, NewSigma, &
                               invDwork, wwork, NumFact, ErrorCode)

    IF (ErrorCode == 0) THEN

      Glob_ApproxEnergy = NewSigma

      ! Discard the previous eigenvector as the starting guess. Essential: an
      ! exact eigenvector of the pencil is a FIXED POINT of inverse iteration
      ! at ANY shift ((H-sigma*S)^{-1} S x is parallel to x), so starting from
      ! the old state's vector would report convergence after one step and
      ! return the old eigenvalue, silently undoing the retargeting. A vector
      ! of ones has a component along every eigenvector and is what the
      ! callers use as the initial guess of a fresh step.
      Glob_LastEigvector(1:N) = ONE

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        IF (Verbose >= 2) WRITE(*, *) 'Eigenvalue-index targeting in ', RoutineName
        WRITE(*, *) '  WHICH_EIGENVALUE      = ', Glob_WhichEigenvalue
        WRITE(*, '(a26)', ADVANCE='no') '   shift before         ='
        CALL writerealadv(6, OldSigma)
        WRITE(*, '(a26)', ADVANCE='no') '   shift after          ='
        CALL writerealadv(6, NewSigma)
        WRITE(*, *) '  factorizations used   = ', NumFact
        WRITE(*, *)
      ENDIF

    ELSE

      IF (Glob_ProcID == 0) THEN
        WRITE(*, *)
        WRITE(*, *) 'Warning: eigenvalue-index targeting failed in ', RoutineName
        WRITE(*, *) 'ErrorCode = ', ErrorCode, ' - continuing with the original shift'
        WRITE(*, *)
      ENDIF

    ENDIF

    DEALLOCATE(wwork)
    DEALLOCATE(invDwork)

  END SUBROUTINE RetargetShiftToEigenvalue


  SUBROUTINE RefreshINVITShift(N)
    !==================================================================
    ! Subroutine RefreshINVITShift
    !==================================================================
    ! Re-anchors the inverse-iteration shift on the energy reached so far
    ! and re-shifts Glob_H to match. Inverse iteration converges at the rate
    ! |lambda_k - sigma| / |lambda_nearest_other - sigma|; a shift set once
    ! per step drifts away from lambda_k as the optimization pushes it down,
    ! until sigma sits between lambda_k and lambda_{k+1} at comparable
    ! distance and the iterate alternates between the two eigenvectors
    ! (|dx| rises every iteration, the solve is thrown away). Setting
    ! sigma = E*Glob_InvItParameter puts it just below lambda_k again and
    ! keeps the inertia count at k-1, so the eigenvalue INDEX is preserved
    ! without bisection. Cost: an O(N^2) update of the lower triangle,
    ! H - sigma_new*S = (H - sigma_old*S) + (sigma_old - sigma_new)*S (unit
    ! diagonal of Glob_S in 'I' mode). THE CALLER MUST REFACTORIZE FROM
    ! ROW 1 afterwards (one EnergyIA(1,N,.FALSE.,ErrCode) does it).
    !==================================================================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N  ! current basis size

    REAL(wp)                            :: SigmaNow, NewSigma, Offset
    INTEGER                             :: it, mcount, FactCode
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: invDwork, wwork

    IF (N < 1) RETURN

    ! EVERY branch below must be taken identically on every rank: the search
    ! loop calls LDLTF, which is collective, so a rank that loops a different
    ! number of times desynchronizes the MPI_ALLREDUCEs inside it and the job
    ! dies with MPI_ERR_TRUNCATE. Broadcasting the two values the branches
    ! depend on makes that impossible by construction rather than by argument.
    ! It costs one broadcast per optimized function.
    CALL MPI_BCAST(Glob_CurrEnergy, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    CALL MPI_BCAST(Glob_ApproxEnergy, 1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)

    ! Never re-anchor on a rejection sentinel - that is how a shift of
    ! -1.19e19 was produced once already.
    IF (ABS(Glob_CurrEnergy) > 1.0E10_wp) RETURN

    SigmaNow = Glob_ApproxEnergy
    Offset = ABS(Glob_CurrEnergy)*(Glob_InvItParameter-ONE)
    IF (Offset <= ZERO) RETURN

    !------------------------------------------------------------------
    ! The offset must FIT IN THE GAP below lambda_k, and INV_IT_PARAM is a
    ! fixed RELATIVE step (1e-6 of |E| at 1.000001) while the gap shrinks
    ! as the basis improves; a too large offset puts sigma below lambda_{k-1}
    ! and IsRequestedEigenstate refuses most solves. So place sigma, count
    ! the eigenvalues below it, and halve the offset until exactly
    ! WhichEigenvalue-1 lie below (one O(N^3/3) factorization per trial,
    ! once per function). With N < WhichEigenvalue keep the plain offset.
    !------------------------------------------------------------------
    IF ((Glob_EigIdxTargeting /= 1) .OR. (N < Glob_WhichEigenvalue)) THEN
      NewSigma = Glob_CurrEnergy-Offset
      CALL ReShiftOnly(N, Glob_H, Glob_HSLeadDim, SigmaNow, Glob_S, Glob_HSLeadDim, NewSigma)
      Glob_ApproxEnergy = SigmaNow
      RETURN
    ENDIF

    ALLOCATE(invDwork(N))
    ALLOCATE(wwork(N))

    DO it = 1, 60
      NewSigma = Glob_CurrEnergy-Offset
      CALL ReShiftAndFactor(N, Glob_H, Glob_HSLeadDim, SigmaNow, &
                            Glob_S, Glob_HSLeadDim, NewSigma, &
                            invDwork, wwork, mcount, FactCode)
      IF (FactCode == 0) THEN
        IF (mcount == Glob_WhichEigenvalue-1) EXIT
        IF (mcount >= Glob_WhichEigenvalue) THEN
          ! sigma came out ABOVE lambda_k - push it further down
          Offset = Offset*TWO
          CYCLE
        ENDIF
      ENDIF
      ! too few below: sigma is under lambda_{k-1}, move it up towards
      ! lambda_k. A singular factorization is handled the same way - the
      ! shift landed on an eigenvalue and any move gets off it.
      Offset = Offset/TWO
    ENDDO

    !------------------------------------------------------------------
    ! The right COUNT is necessary but not sufficient: sigma may sit in the
    ! lower half of (lambda_{k-1}, lambda_k) and inverse iteration then
    ! converges to lambda_{k-1}. Halving from above leaves the accepted
    ! offset in (gap/2, gap]; two more halvings put it in (gap/8, gap/4],
    ! so lambda_k is at most gap/4 away and lambda_{k-1} at least 3*gap/4:
    ! a convergence ratio of 1/3 or better. Shrinking the offset only moves
    ! sigma up towards lambda_k, so the inertia count cannot change.
    !------------------------------------------------------------------
    Offset = Offset/FOUR
    NewSigma = Glob_CurrEnergy-Offset
    CALL ReShiftAndFactor(N, Glob_H, Glob_HSLeadDim, SigmaNow, &
                          Glob_S, Glob_HSLeadDim, NewSigma, &
                          invDwork, wwork, mcount, FactCode)

    DEALLOCATE(wwork)
    DEALLOCATE(invDwork)

    Glob_ApproxEnergy = SigmaNow

  END SUBROUTINE RefreshINVITShift


  FUNCTION IsEigenpairUsable(ErrorCode)
    !==================================================================
    ! Function IsEigenpairUsable
    !==================================================================
    ! Verdict on a GSEPIIS result, shared by every inverse-iteration solve.
    ! ErrorCode=2 (Glob_EigvalTol not reached within Glob_MaxIterForGSEPIIS
    ! iterations) is not by itself a failure: the attainable residual is set
    ! by the conditioning of the overlap matrix and the eigenpair is usually
    ! far more accurate than the optimizer needs. The result is rejected
    ! only when the achieved residual Glob_LastEigvalTol exceeds
    ! Glob_EigvalTolUsable, or when the factorization failed (ErrorCode=1).
    !==================================================================

    IMPLICIT NONE

    LOGICAL             :: IsEigenpairUsable  ! function result
    INTEGER, INTENT(IN) :: ErrorCode          ! as returned by GSEPIIS

    SELECT CASE (ErrorCode)
    CASE (0)
      IsEigenpairUsable = .TRUE.
    CASE (2)
      IsEigenpairUsable = (Glob_LastEigvalTol < Glob_EigvalTolUsable)
    CASE DEFAULT
      IsEigenpairUsable = .FALSE.
    END SELECT

  END FUNCTION IsEigenpairUsable


  FUNCTION IsRequestedEigenstate(Nmax)
    !==================================================================
    ! Function IsRequestedEigenstate
    !==================================================================
    ! .TRUE. if the eigenpair GSEPIIS just returned is eigenvalue number
    ! Glob_WhichEigenvalue. Inverse iteration returns the level nearest the
    ! shift, which is held fixed during a step while the basis moves, so a
    ! trial point can push a LOWER level across the shift; its energy is
    ! lower, the optimizer takes the step and from then on minimizes that
    ! lower state (variational collapse - the Hylleraas-Undheim-MacDonald
    ! bound holds only while lambda_k IS lambda_k). Rejecting the point
    ! costs nothing: the optimizer backs off and the shift is re-anchored
    ! at the next step. Returns .TRUE. (does not interfere) when index
    ! targeting is off (Glob_EigIdxTargeting/=1), when the basis is smaller
    ! than the requested index, or when the index is unknown (factorization
    ! failed - already rejected by the caller).
    !==================================================================

    IMPLICIT NONE

    LOGICAL             :: IsRequestedEigenstate  ! function result
    INTEGER, INTENT(IN) :: Nmax                   ! size of the basis just solved

    IsRequestedEigenstate = .TRUE.

    IF (Glob_EigIdxTargeting /= 1) RETURN
    IF (Nmax < Glob_WhichEigenvalue) RETURN
    IF (Glob_LastEigIndex <= 0) RETURN

    IsRequestedEigenstate = (Glob_LastEigIndex == Glob_WhichEigenvalue)

  END FUNCTION IsRequestedEigenstate


END MODULE matform
