MODULE linalg  ! version 0.52
  ! Module linalg contains some linear algebra routines
  ! that are used in calculations
  USE globvars
  IMPLICIT NONE

  ! Mode switches of the routines below, set by linalg_setparam for the current problem size:
  !   *_PMode      1 = distributed over the MPI processes, 0 = every process computes the whole result
  !   *_SModeType  in the serial branch: 1 = every process computes, 0 = rank 0 computes and broadcasts
  !   *_UseBLAS    1 = call BLAS in the serial branch (only when a BLAS library is linked)

  INTEGER :: Glob_LDLTF_PMode = 1
  INTEGER :: Glob_LDLTF_SModeType = 0

  INTEGER :: Glob_LDLHF_PMode = 1
  INTEGER :: Glob_LDLHF_SModeType = 0

  INTEGER :: Glob_LDLTS_PMode = 1
  INTEGER :: Glob_LDLTS_SModeType = 0
  INTEGER :: Glob_LDLTS_UseBLAS = 0

  INTEGER :: Glob_LDLHS_PMode = 1
  INTEGER :: Glob_LDLHS_SModeType = 0
  INTEGER :: Glob_LDLHS_UseBLAS = 0

  INTEGER :: Glob_MTMVL_PMode = 1
  INTEGER :: Glob_MTMVL_UseBLAS = 0

  INTEGER :: Glob_MHMVL_PMode = 1
  INTEGER :: Glob_MHMVL_UseBLAS = 0

  INTEGER :: Glob_MTMV_PMode = 1
  INTEGER :: Glob_MTMV_UseBLAS = 0

  INTEGER :: Glob_MHMV_PMode = 1
  INTEGER :: Glob_MHMV_UseBLAS = 0

  INTEGER :: Glob_VMMTMV_PMode = 1

  INTEGER :: Glob_VMMHMV_PMode = 1

  INTEGER :: Glob_RMaxAbsEl_PMode = 0

  INTEGER :: Glob_CMaxAbsReOrIm_PMode = 0

  INTEGER :: Glob_RDotProd_PMode = 0
  INTEGER :: Glob_RDotProd_UseBLAS = 0

  INTEGER :: Glob_CDotProd_PMode = 0
  INTEGER :: Glob_CDotProd_UseBLAS = 0

  INTEGER :: Glob_RDotProdItself_PMode = 0
  INTEGER :: Glob_RDotProdItself_UseBLAS = 0

  INTEGER :: Glob_CDotProdItself_PMode = 0
  INTEGER :: Glob_CDotProdItself_UseBLAS = 0

  INTEGER :: Glob_RDotProdQuotient_PMode = 0

  INTEGER :: Glob_CDotProdQuotient_PMode = 0

  INTEGER :: Glob_RVScale_UseBLAS = 0

  INTEGER :: Glob_CVScale_UseBLAS = 0

  INTEGER :: Glob_RVDiffEucNorm_PMode = 0

  INTEGER :: Glob_CVDiffEucNorm_PMode = 0

  ! Panel width of the cache-blocked serial LDLT/LDLH factorization (measured optimum about 48 on
  ! x86-64). Updates of at least this many columns count as bulk factorizations for the routing.
  INTEGER, PARAMETER :: Glob_LDLF_BlockSize = 48

  ! Cost models routing bulk factorizations (wp=8 only; R real, C complex):
  !   t_serial = Cs*n^3,  t_parallel = Ca*n^3 + Cb*n^2.
  ! Measured once per run by linalg_ldlf_calibrate; Cs < 0 = not calibrated (parallel branch used).
  REAL(wp) :: Glob_LDLF_CsR = -1.0_wp, Glob_LDLF_CaR = 0.0_wp, Glob_LDLF_CbR = 0.0_wp
  REAL(wp) :: Glob_LDLF_CsC = -1.0_wp, Glob_LDLF_CaC = 0.0_wp, Glob_LDLF_CbC = 0.0_wp

  ! =====================================================================
  ! Eigenvalue-index targeting : the number of negative D_ii of
  ! A - sigma*B = L*D*L^T is the number of eigenvalues below sigma (Sylvester), so GSEPIIS
  ! knows the index of the eigenvalue it returns and GSEPIIS_ShiftForIndex can bisect on it.
  ! Drivers: RetargetShiftToEigenvalue, RefreshINVITShift, IsRequestedEigenstate (matform).
  ! =====================================================================

  ! Master switch. 0 = off (reference behaviour), 1 = on (default).
  ! main.f90 overrides it from the environment variable ECG_EIG_IDX_TARGETING.
  INTEGER :: Glob_EigIdxTargeting = 1

  ! Maximum number of bisection steps. Each one costs a factorization.
  ! 60 is far more than a double-precision bracket can ever need.
  INTEGER :: Glob_EigIdxMaxBisect = 60

  ! Relative width at which the bracket is considered converged.
  REAL(wp) :: Glob_EigIdxBisectTol = 1.0e-12_wp

  ! Smallest Rayleigh quotient of B, x^T B x / x^T x, that GSEPIIS accepts
  ! as coming from a genuine eigenvector - see the guard in GSEPIIS. Set
  ! it to zero to restore the reference behaviour of dividing
  ! unconditionally.
  REAL(wp) :: Glob_MinRayleighDenom = 1.0e-8_wp

CONTAINS

  SUBROUTINE linalg_setparam(n)
    ! Sets the mode switches for a problem of size n. Call it whenever n changes, and
    ! collectively at np > 1 (the calibration uses MPI). O(n^2)/O(n^3) routines run serially
    ! at np = 1 and in their parallel branches at np > 1; the O(n) helpers stay serial. For wp=8
    ! LDLTF/LDLHF route bulk factorizations between the blocked serial path and the parallel
    ! branch with the calibrated cost models; thin updates stay in the parallel branch.
    INTEGER :: n
    INTEGER :: pm, min_procs_par

    ! Smallest process count at which the parallel branches pay off, selected per
    ! working precision. Benchmarks place it at 2 for all three supported kinds;
    ! it is kept as an explicit per-wp knob so it can be retuned independently.
    SELECT CASE (wp)
    CASE (8)  ! double precision
      min_procs_par = 2
    CASE (10)  ! extended (80-bit) precision
      min_procs_par = 2
    CASE (16)  ! quadruple precision (needs Intel MPI for the parallel branches)
      min_procs_par = 2
    CASE DEFAULT
      min_procs_par = 2
    END SELECT

    IF (Glob_NumOfProcs >= min_procs_par) THEN
      pm = 1  ! parallel
    ELSE
      pm = 0  ! single process - avoid all MPI call overhead
    ENDIF

    ! Heavy routines: parallel at np > 1 (LDLTF/LDLHF decide internally between the blocked
    ! serial path and the parallel branch).
    Glob_LDLTF_PMode = pm
    Glob_LDLHF_PMode = pm
    Glob_LDLTS_PMode = pm
    Glob_LDLHS_PMode = pm
    Glob_MTMVL_PMode = pm
    Glob_MHMVL_PMode = pm
    Glob_MTMV_PMode = pm
    Glob_MHMV_PMode = pm
    Glob_VMMTMV_PMode = pm
    Glob_VMMHMV_PMode = pm

    ! SModeType 0 at np > 1: rank 0 factorizes alone and broadcasts the factor (faster than
    ! redundant computation on a shared-memory node); 1 at np = 1.
    Glob_LDLTF_SModeType = 1-pm
    Glob_LDLHF_SModeType = 1-pm

    ! O(n) vector helpers: always serial at these problem sizes.
    Glob_RMaxAbsEl_PMode = 0
    Glob_CMaxAbsReOrIm_PMode = 0
    Glob_RDotProd_PMode = 0
    Glob_CDotProd_PMode = 0
    Glob_RDotProdItself_PMode = 0
    Glob_CDotProdItself_PMode = 0
    Glob_RDotProdQuotient_PMode = 0
    Glob_CDotProdQuotient_PMode = 0
    Glob_RVDiffEucNorm_PMode = 0
    Glob_CVDiffEucNorm_PMode = 0

    ! Calibrate the routing cost models once per run (wp=8, np > 1; MPI collectives inside).
    IF ((wp == 8) .AND. (Glob_NumOfProcs > 1) .AND. (Glob_LDLF_CsR < ZERO)) THEN
      CALL linalg_ldlf_calibrate
    ENDIF

    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,i0,a,i0,a,l1)') 'linalg_setparam: n=', n, ' np=', Glob_NumOfProcs, ' PMode=', &
        pm, ' LDLTF SModeType=', Glob_LDLTF_SModeType, ' calibrated=', (Glob_LDLF_CsR >= ZERO)
    ENDIF
  END SUBROUTINE linalg_setparam

  SUBROUTINE linalg_ldlf_calibrate
    ! Times the blocked serial path and the parallel branch of LDLTF/LDLHF on synthetic
    ! matrices of size 320 and 640, fits t_serial = Cs*n^3 and t_parallel = Ca*n^3 + Cb*n^2
    ! and broadcasts the coefficients. Collective; tens of milliseconds once per run.
    INTEGER, PARAMETER       :: m1 = 320, m2 = 640, nrep = 2
    REAL(wp), ALLOCATABLE    :: AR0(:, :), AR(:, :), invDR(:), wR(:)
    COMPLEX(wp), ALLOCATABLE :: AC0(:, :), AC(:, :), invDC(:), wC(:)
    REAL(wp)                 :: dr, tp1r, tp2r, tp1c, tp2c, tsr, tsc, cf(6)
    DOUBLE PRECISION         :: tt0, tt1, tbest
    INTEGER                  :: i, j, r, ec
    INTEGER                  :: savePMr, saveSMr, savePMc, saveSMc

    savePMr = Glob_LDLTF_PMode
    saveSMr = Glob_LDLTF_SModeType
    savePMc = Glob_LDLHF_PMode
    saveSMc = Glob_LDLHF_SModeType

    ALLOCATE(AR0(m2, m2), AR(m2, m2), invDR(m2), wR(m2))
    ALLOCATE(AC0(m2, m2), AC(m2, m2), invDC(m2), wC(m2))
    ! Synthetic diagonally dominant test matrices (only the lower triangles
    ! are referenced by the factorizations)
    DO j = 1, m2
      DO i = j, m2
        dr = ONE/(ONE+ABS(i-j))
        AR0(i, j) = dr
        IF (i == j) THEN
          AC0(i, j) = CMPLX(dr+3*ONE, ZERO, wp)
        ELSE
          AC0(i, j) = CMPLX(dr, (3*dr)/(10*(1+MOD(i+j, 7))), wp)
        ENDIF
      ENDDO
      AR0(j, j) = AR0(j, j)+3*ONE
    ENDDO

    !--- Parallel branch, collectively, at sizes m1 and m2. The routing in
    ! LDLTF/LDLHF sends bulk factorizations to the parallel branch as long
    ! as the models are not calibrated yet (Cs<0), which is exactly the
    ! state here, so a plain call with PMode=1 measures the parallel branch.
    Glob_LDLTF_PMode = 1
    Glob_LDLHF_PMode = 1
    tbest = HUGE(tbest)
    DO r = 1, nrep
      AR(1:m1, 1:m1) = AR0(1:m1, 1:m1)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt0 = MPI_WTIME()
      CALL LDLTF(1, m1, AR, m2, invDR, wR, ec)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt1 = MPI_WTIME()
      IF (tt1-tt0 < tbest) tbest = tt1-tt0
    ENDDO
    tp1r = REAL(tbest, wp)
    tbest = HUGE(tbest)
    DO r = 1, nrep
      AR = AR0
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt0 = MPI_WTIME()
      CALL LDLTF(1, m2, AR, m2, invDR, wR, ec)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt1 = MPI_WTIME()
      IF (tt1-tt0 < tbest) tbest = tt1-tt0
    ENDDO
    tp2r = REAL(tbest, wp)
    tbest = HUGE(tbest)
    DO r = 1, nrep
      AC(1:m1, 1:m1) = AC0(1:m1, 1:m1)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt0 = MPI_WTIME()
      CALL LDLHF(1, m1, AC, m2, invDC, wC, ec)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt1 = MPI_WTIME()
      IF (tt1-tt0 < tbest) tbest = tt1-tt0
    ENDDO
    tp1c = REAL(tbest, wp)
    tbest = HUGE(tbest)
    DO r = 1, nrep
      AC = AC0
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt0 = MPI_WTIME()
      CALL LDLHF(1, m2, AC, m2, invDC, wC, ec)
      CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)
      tt1 = MPI_WTIME()
      IF (tt1-tt0 < tbest) tbest = tt1-tt0
    ENDDO
    tp2c = REAL(tbest, wp)

    !--- Blocked serial path at size m2, on rank 0 alone (the other ranks
    ! wait at the barrier), matching how the serial path actually runs at
    ! np>1: PMode=0 forces the serial path, SModeType=1 removes its MPI
    ! broadcasts so the call is local to each process.
    Glob_LDLTF_PMode = 0
    Glob_LDLHF_PMode = 0
    Glob_LDLTF_SModeType = 1
    Glob_LDLHF_SModeType = 1
    tsr = ZERO
    tsc = ZERO
    IF (Glob_ProcID == 0) THEN
      tbest = HUGE(tbest)
      DO r = 1, nrep
        AR = AR0
        tt0 = MPI_WTIME()
        CALL LDLTF(1, m2, AR, m2, invDR, wR, ec)
        tt1 = MPI_WTIME()
        IF (tt1-tt0 < tbest) tbest = tt1-tt0
      ENDDO
      tsr = REAL(tbest, wp)
      tbest = HUGE(tbest)
      DO r = 1, nrep
        AC = AC0
        tt0 = MPI_WTIME()
        CALL LDLHF(1, m2, AC, m2, invDC, wC, ec)
        tt1 = MPI_WTIME()
        IF (tt1-tt0 < tbest) tbest = tt1-tt0
      ENDDO
      tsc = REAL(tbest, wp)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, Glob_MPIErrCode)

    Glob_LDLTF_PMode = savePMr
    Glob_LDLTF_SModeType = saveSMr
    Glob_LDLHF_PMode = savePMc
    Glob_LDLHF_SModeType = saveSMc

    !--- Fit the models on rank 0 and broadcast the coefficients so that
    ! every process takes identical routing decisions (different decisions
    ! on different ranks would deadlock the collective branches).
    IF (Glob_ProcID == 0) THEN
      CALL fitpar(tp1r, tp2r, cf(1), cf(2))
      CALL fitpar(tp1c, tp2c, cf(4), cf(5))
      cf(3) = tsr/(REAL(m2, wp)**3)
      cf(6) = tsc/(REAL(m2, wp)**3)
    ENDIF
    CALL MPI_BCAST(cf, 6, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    Glob_LDLF_CaR = cf(1)
    Glob_LDLF_CbR = cf(2)
    Glob_LDLF_CsR = cf(3)
    Glob_LDLF_CaC = cf(4)
    Glob_LDLF_CbC = cf(5)
    Glob_LDLF_CsC = cf(6)

    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,3es11.3,a)') 'linalg_ldlf_calibrate: t_par(320), t_par(640), t_ser(640) =', tp1r, tp2r, &
        tsr, ' s'
    ENDIF
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,3es11.3)') '  real cost models Cs, Ca, Cb =', Glob_LDLF_CsR, Glob_LDLF_CaR, Glob_LDLF_CbR
    ENDIF

  CONTAINS

    SUBROUTINE fitpar(t1, t2, a, b)
      ! Fits t = a*m^3 + b*m^2 through (m1,t1) and (m2,t2), clamping both
      ! coefficients to be nonnegative (protects the extrapolation against
      ! measurement noise; when one coefficient is clamped the other is
      ! recomputed from the larger, more reliable measurement).
      REAL(wp) :: t1, t2, a, b
      REAL(wp) :: f1, f2
      f1 = REAL(m1, wp)
      f2 = REAL(m2, wp)
      a = (t1*f2**2-t2*f1**2)/(f1**3*f2**2-f2**3*f1**2)
      IF (a < ZERO) THEN
        a = ZERO
        b = t2/f2**2
      ELSE
        b = (t1-a*f1**3)/f1**2
        IF (b < ZERO) THEN
          b = ZERO
          a = t2/f2**3
        ENDIF
      ENDIF
    END SUBROUTINE fitpar

  END SUBROUTINE linalg_ldlf_calibrate

  SUBROUTINE LDLTF(m, n, A, nA, invD, w, ErrorCode)
    ! Updates the factorization A = L*D*L^T of the real symmetric n x n matrix A, given the
    ! factorization of its leading (m-1) x (m-1) block (m = 1: factorize everything).
    !   A(nA,n)    lower triangle incl. diagonal: A (unchanged); upper triangle: L^T (in: size m-1, out: n)
    !   invD(n)    1/D_ii (in: 1..m-1, out: 1..n);   w(n) work
    !   ErrorCode  0 ok, 1 zero pivot (A singular)

    ! Arguments :
    INTEGER  :: m, n, nA, ErrorCode
    REAL(wp) :: A(nA, n), invD(n), w(n)
    ! Local variables :
    INTEGER               :: i, j, jm, im, i0, ib, ibs
    INTEGER               :: q, k, ji, jf, jim, jiR, RowsPerProc, mod_im_Glob_NumOfProcs
    REAL(wp)              :: x, y, z
    REAL(wp)              :: tv(Glob_LDLF_BlockSize)
    REAL(wp), ALLOCATABLE :: Pnl(:, :)
    LOGICAL               :: serialbulk

    ErrorCode = 0
    ! Routing: PMode 0 = serial path. wp=8 with at least one panel of new columns: blocked serial
    ! path at np = 1 or when the cost models prefer it (rank 0 computes, SModeType 0), otherwise
    ! the parallel branch. wp=10/16: plain column algorithm (ibs = 1).
    ibs = Glob_LDLF_BlockSize
    IF (wp /= 8) ibs = 1
    serialbulk = (Glob_LDLTF_PMode == 0)
    IF ((.NOT. serialbulk) .AND. (wp == 8) .AND. (n-m+1 >= Glob_LDLF_BlockSize)) THEN
      IF (Glob_NumOfProcs == 1) THEN
        serialbulk = .TRUE.
      ELSEIF (Glob_LDLF_CsR >= ZERO) THEN
        ! t_ser(n)<=t_par(n) with both sides divided by n^2
        IF (Glob_LDLF_CsR*n <= Glob_LDLF_CaR*n+Glob_LDLF_CbR) serialbulk = .TRUE.
      ENDIF
    ENDIF
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,a,a,i0)') 'LDLTF: m=', m, ' n=', n, ' path=', &
        MERGE('serial  ', 'parallel', serialbulk), ' panel=', ibs
    ENDIF
    IF (serialbulk) THEN
      ! Serial path, blocked: columns m..n in panels of Glob_LDLF_BlockSize held transposed in Pnl,
      ! each completed column streamed once per panel. Bitwise identical to the unblocked algorithm.
      IF ((Glob_LDLTF_SModeType /= 0) .OR. (Glob_ProcID == 0)) THEN
        IF ((ibs > 1) .AND. (n > m)) ALLOCATE(Pnl(ibs, n))
        DO i0 = m, n, ibs
          ib = MIN(ibs, n-i0+1)
          IF (ib == 1) THEN
            ! Plain column algorithm (original unblocked code)
            i = i0
            DO j = 1, i-1
              jm = j-1
              x = A(i, j)
              DO q = 1, jm
                x = x-A(q, i)*A(q, j)
              ENDDO
              A(j, i) = x
            ENDDO
            im = i-1
            DO q = 1, im
              w(q) = A(q, i)*invD(q)
            ENDDO
            x = A(i, i)
            DO q = 1, im
              x = x-A(q, i)*w(q)
            ENDDO
            DO q = 1, im
              A(q, i) = w(q)
            ENDDO
            IF (x == ZERO) THEN
              ErrorCode = 1
              RETURN
            ENDIF
            invD(i) = ONE/x
          ELSE
            ! Seed the panel buffer with the matrix elements of rows
            ! i0..i0+ib-1 taken from the (untouched) lower triangle
            DO k = 1, ib
              i = i0+k-1
              DO j = 1, i-1
                Pnl(k, j) = A(i, j)
              ENDDO
            ENDDO
            ! Shared sweep: forward-substitute all panel columns against
            ! the completed columns 1..i0-1 (this is where blocking pays)
            DO j = 1, i0-1
              jm = j-1
              DO k = 1, ib
                tv(k) = Pnl(k, j)
              ENDDO
              DO q = 1, jm
                x = A(q, j)
                DO k = 1, ib
                  tv(k) = tv(k)-Pnl(k, q)*x
                ENDDO
              ENDDO
              DO k = 1, ib
                Pnl(k, j) = tv(k)
              ENDDO
            ENDDO
            ! Finish the panel columns one by one: substitution against the
            ! panel columns completed just before them and the diagonal
            ! step, then write the scaled column of the factor back into
            ! the upper triangle of A
            DO k = 1, ib
              i = i0+k-1
              DO j = i0, i-1
                jm = j-1
                x = Pnl(k, j)
                DO q = 1, jm
                  x = x-Pnl(k, q)*A(q, j)
                ENDDO
                Pnl(k, j) = x
              ENDDO
              im = i-1
              DO q = 1, im
                w(q) = Pnl(k, q)*invD(q)
              ENDDO
              x = A(i, i)
              DO q = 1, im
                x = x-Pnl(k, q)*w(q)
              ENDDO
              DO q = 1, im
                A(q, i) = w(q)
              ENDDO
              IF (x == ZERO) THEN
                ErrorCode = 1
                RETURN
              ENDIF
              invD(i) = ONE/x
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF (Glob_LDLTF_SModeType == 0) THEN
        DO i = m, n
          CALL MPI_BCAST(A(1:i-1, i), i-1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        ENDDO
        CALL MPI_BCAST(invD(m:n), n-m+1, MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
    ELSE
      ! Parallel version (far from being perfect...)
      DO i = m, n
        im = i-1
        A(1:im, i) = ZERO
        RowsPerProc = im/Glob_NumOfProcs
        mod_im_Glob_NumOfProcs = MOD(im, Glob_NumOfProcs)
        DO k = 1, RowsPerProc
          jf = k*Glob_NumOfProcs
          jim = jf-Glob_NumOfProcs
          ji = jim+1
          jiR = ji+Glob_ProcID
          ! A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
          z = ZERO
          DO q = 1, jim
            z = z-A(q, i)*A(q, jiR)
          ENDDO
          A(jiR, i) = z
          !
          CALL MPI_ALLREDUCE(A(ji:jf, i), w(1:Glob_NumOfProcs), Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          A(ji:jf, i) = w(1:Glob_NumOfProcs)
          DO j = ji, jf
            jm = j-1
            ! A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
            z = A(j, i)+A(i, j)
            DO q = ji, jm
              z = z-A(q, i)*A(q, j)
            ENDDO
            A(j, i) = z
            !
          ENDDO
        ENDDO
        IF (mod_im_Glob_NumOfProcs > 0) THEN
          ji = RowsPerProc*Glob_NumOfProcs+1
          jim = ji-1
          jf = im
          jiR = ji+Glob_ProcID
          IF (jiR < i) THEN
            ! A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
            z = ZERO
            DO q = 1, jim
              z = z-A(q, i)*A(q, jiR)
            ENDDO
            A(jiR, i) = z
            !
          ENDIF
          CALL MPI_ALLREDUCE(A(ji:jf, i), w(1:mod_im_Glob_NumOfProcs), mod_im_Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          A(ji:jf, i) = w(1:mod_im_Glob_NumOfProcs)
          DO j = ji, jf
            jm = j-1
            ! A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
            z = A(j, i)+A(i, j)
            DO q = ji, jm
              z = z-A(q, i)*A(q, j)
            ENDDO
            A(j, i) = z
            !
          ENDDO
        ENDIF
        ! j==1 case
        w(1:im) = A(1:im, i)*invD(1:im)
        y = ZERO
        DO k = 1+Glob_ProcID, im, Glob_NumOfProcs
          y = y+A(k, i)*w(k)
        ENDDO
        CALL MPI_ALLREDUCE(y, x, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x = A(i, i)-x
        A(1:im, i) = w(1:im)
        IF (x == ZERO) THEN
          ErrorCode = 1
          RETURN
        ENDIF
        invD(i) = ONE/x
      ENDDO
    ENDIF

  END SUBROUTINE LDLTF

  SUBROUTINE LDLHF(m, n, A, nA, invD, w, ErrorCode)
    ! Updates the factorization A = L*D*L^H of the complex hermitian n x n matrix A, given the
    ! factorization of its leading (m-1) x (m-1) block (m = 1: factorize everything).
    !   A(nA,n)    lower triangle incl. diagonal: A (unchanged); upper triangle: L^H (in: size m-1, out: n)
    !   invD(n)    1/D_ii (in: 1..m-1, out: 1..n);   w(n) work
    !   ErrorCode  0 ok, 1 zero pivot (A singular)

    ! Arguments :
    INTEGER     :: m, n, nA, ErrorCode
    COMPLEX(wp) :: A(nA, n), invD(n), w(n)
    ! Local variables :
    INTEGER               :: i, j, jm, im, i0, ib, ibs
    INTEGER               :: q, k, ji, jf, jim, jiR, RowsPerProc, mod_im_Glob_NumOfProcs
    COMPLEX(wp)           :: x, y
    REAL(wp)              :: xr, xi
    REAL(wp)              :: tvr(Glob_LDLF_BlockSize), tvi(Glob_LDLF_BlockSize)
    REAL(wp), ALLOCATABLE :: PnR(:, :), PnI(:, :)
    LOGICAL               :: serialbulk

    ErrorCode = 0
    ! Routing between the blocked serial path and the parallel branch: same
    ! rule and rationale as in LDLTF (see the comments there), including the
    ! restriction of the blocked kernel to wp=8, with the complex-case
    ! performance-model coefficients.
    ibs = Glob_LDLF_BlockSize
    IF (wp /= 8) ibs = 1
    serialbulk = (Glob_LDLHF_PMode == 0)
    IF ((.NOT. serialbulk) .AND. (wp == 8) .AND. (n-m+1 >= Glob_LDLF_BlockSize)) THEN
      IF (Glob_NumOfProcs == 1) THEN
        serialbulk = .TRUE.
      ELSEIF (Glob_LDLF_CsC >= ZERO) THEN
        ! t_ser(n)<=t_par(n) with both sides divided by n^2
        IF (Glob_LDLF_CsC*n <= Glob_LDLF_CaC*n+Glob_LDLF_CbC) serialbulk = .TRUE.
      ENDIF
    ENDIF
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,a,a,i0)') 'LDLHF: m=', m, ' n=', n, ' path=', &
        MERGE('serial  ', 'parallel', serialbulk), ' panel=', ibs
    ENDIF
    IF (serialbulk) THEN
      ! Serial path, blocked like LDLTF; the panel is kept as real and imaginary planes (PnR, PnI)
      ! so the kernels are real multiply-adds. May differ from the unblocked algorithm in the last bits.
      IF ((Glob_LDLHF_SModeType /= 0) .OR. (Glob_ProcID == 0)) THEN
        IF ((ibs > 1) .AND. (n > m)) ALLOCATE(PnR(ibs, n), PnI(ibs, n))
        DO i0 = m, n, ibs
          ib = MIN(ibs, n-i0+1)
          IF (ib == 1) THEN
            ! Plain column algorithm (original unblocked code)
            i = i0
            DO j = 1, i-1
              jm = j-1
              x = A(i, j)
              DO q = 1, jm
                x = x-CONJG(A(q, i))*A(q, j)
              ENDDO
              A(j, i) = CONJG(x)
            ENDDO
            im = i-1
            DO q = 1, im
              w(q) = A(q, i)*invD(q)
            ENDDO
            x = A(i, i)
            DO q = 1, im
              x = x-CONJG(A(q, i))*w(q)
            ENDDO
            x = CONJG(x)
            DO q = 1, im
              A(q, i) = w(q)
            ENDDO
            IF (x == ZERO) THEN
              ErrorCode = 1
              RETURN
            ENDIF
            invD(i) = ONE/x
          ELSE
            ! Seed the panel buffers with the matrix elements of rows
            ! i0..i0+ib-1 taken from the (untouched) lower triangle
            DO k = 1, ib
              i = i0+k-1
              DO j = 1, i-1
                PnR(k, j) = REAL(A(i, j), wp)
                PnI(k, j) = AIMAG(A(i, j))
              ENDDO
            ENDDO
            ! Shared sweep: forward-substitute all panel columns against
            ! the completed columns 1..i0-1 (this is where blocking pays)
            DO j = 1, i0-1
              jm = j-1
              DO k = 1, ib
                tvr(k) = PnR(k, j)
                tvi(k) = PnI(k, j)
              ENDDO
              DO q = 1, jm
                xr = REAL(A(q, j), wp)
                xi = AIMAG(A(q, j))
                DO k = 1, ib
                  tvr(k) = tvr(k)-(PnR(k, q)*xr-PnI(k, q)*xi)
                  tvi(k) = tvi(k)-(PnR(k, q)*xi+PnI(k, q)*xr)
                ENDDO
              ENDDO
              DO k = 1, ib
                PnR(k, j) = tvr(k)
                PnI(k, j) = tvi(k)
              ENDDO
            ENDDO
            ! Finish the panel columns one by one: substitution against the
            ! panel columns completed just before them and the diagonal
            ! step, then write the scaled column of the factor back into
            ! the upper triangle of A
            DO k = 1, ib
              i = i0+k-1
              DO j = i0, i-1
                jm = j-1
                xr = PnR(k, j)
                xi = PnI(k, j)
                DO q = 1, jm
                  xr = xr-(PnR(k, q)*REAL(A(q, j), wp)-PnI(k, q)*AIMAG(A(q, j)))
                  xi = xi-(PnR(k, q)*AIMAG(A(q, j))+PnI(k, q)*REAL(A(q, j), wp))
                ENDDO
                PnR(k, j) = xr
                PnI(k, j) = xi
              ENDDO
              im = i-1
              DO q = 1, im
                w(q) = CMPLX(PnR(k, q), -PnI(k, q), wp)*invD(q)
              ENDDO
              x = A(i, i)
              DO q = 1, im
                x = x-CMPLX(PnR(k, q), PnI(k, q), wp)*w(q)
              ENDDO
              x = CONJG(x)
              DO q = 1, im
                A(q, i) = w(q)
              ENDDO
              IF (x == ZERO) THEN
                ErrorCode = 1
                RETURN
              ENDIF
              invD(i) = ONE/x
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF (Glob_LDLHF_SModeType == 0) THEN
        DO i = m, n
          CALL MPI_BCAST(A(1:i-1, i), 2*(i-1), MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        ENDDO
        CALL MPI_BCAST(invD(m:n), 2*(n-m+1), MPI_WP, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
        CALL MPI_BCAST(ErrorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
      ENDIF
    ELSE
      ! Parallel version (far from being perfect...)
      DO i = m, n
        im = i-1
        A(1:im, i) = ZERO
        RowsPerProc = im/Glob_NumOfProcs
        mod_im_Glob_NumOfProcs = MOD(im, Glob_NumOfProcs)
        DO k = 1, RowsPerProc
          jf = k*Glob_NumOfProcs
          jim = jf-Glob_NumOfProcs
          ji = jim+1
          jiR = ji+Glob_ProcID
          A(jiR, i) = CONJG(-DOT_PRODUCT(A(1:jim, i), A(1:jim, jiR)))
          CALL MPI_ALLREDUCE(A(ji:jf, i), w(1:Glob_NumOfProcs), 2*Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          A(ji:jf, i) = w(1:Glob_NumOfProcs)
          DO j = ji, jf
            jm = j-1
            A(j, i) = A(j, i)+CONJG(A(i, j)-DOT_PRODUCT(A(ji:jm, i), A(ji:jm, j)))
          ENDDO
        ENDDO
        IF (mod_im_Glob_NumOfProcs > 0) THEN
          ji = RowsPerProc*Glob_NumOfProcs+1
          jim = ji-1
          jf = im
          jiR = ji+Glob_ProcID
          IF (jiR < i) THEN
            A(jiR, i) = CONJG(-DOT_PRODUCT(A(1:jim, i), A(1:jim, jiR)))
          ENDIF
          CALL MPI_ALLREDUCE(A(ji:jf, i), w(1:mod_im_Glob_NumOfProcs), 2*mod_im_Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          A(ji:jf, i) = w(1:mod_im_Glob_NumOfProcs)
          DO j = ji, jf
            jm = j-1
            A(j, i) = A(j, i)+CONJG(A(i, j)-DOT_PRODUCT(A(ji:jm, i), A(ji:jm, j)))
          ENDDO
        ENDIF
        ! j==1 case
        w(1:im) = A(1:im, i)*invD(1:im)
        y = ZERO
        DO k = 1+Glob_ProcID, im, Glob_NumOfProcs
          y = y+CONJG(A(k, i))*w(k)
        ENDDO
        y = CONJG(y)
        CALL MPI_ALLREDUCE(y, x, 2, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x = A(i, i)-x
        A(1:im, i) = w(1:im)
        IF (x == ZERO) THEN
          ErrorCode = 1
          RETURN
        ENDIF
        invD(i) = ONE/x
      ENDDO
    ENDIF

  END SUBROUTINE LDLHF

  SUBROUTINE LDLTS(n, A, nA, invD, b, x)
    ! Solves A*x = b from the factorization A = L*D*L^T of LDLTF (L*y = b, then D*L^T*x = y).
    !   A(nA,n)  L^T in the upper triangle (unit diagonal implied);   invD  1/D_ii
    !   b        right-hand side, destroyed in the parallel branch;   x  solution

    ! Arguments
    INTEGER  :: n, nA
    REAL(wp) :: A(nA, n), invD(n), b(n), x(n)
    ! Local variables
    INTEGER  :: i, j, k, ji, jf, jim, jfm, RowsPerProc, mod_n_Glob_NumOfProcs, jiR
    REAL(wp) :: t

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'LDLTS: n=', n, ' PMode=', Glob_LDLTS_PMode
    ENDIF
    IF (Glob_LDLTS_PMode == 0) THEN
      IF (Glob_LDLTS_UseBLAS == 0) THEN
        ! Solution of L*y=b
        x(1:n) = b(1:n)
        DO j = 1, n
          t = x(j)
          DO i = 1, j-1
            t = t-A(i, j)*x(i)
          ENDDO
          x(j) = t
        ENDDO
        ! Solution of D*LT*x=y
        x(1:n) = x(1:n)*invD(1:n)
        DO j = N, 1, -1
          t = x(j)
          DO i = j-1, 1, -1
            x(i) = x(i)-t*A(i, j)
          ENDDO
        ENDDO
      ELSE
        x(1:n) = b(1:n)
        ! call BLAS routine DTRSV
        CALL DTRSV('U', 'T', 'U', n, A, nA, x, 1)
        x(1:n) = x(1:n)*invD(1:n)
        ! call BLAS routine DTRSV
        CALL DTRSV('U', 'N', 'U', n, A, nA, x, 1)
      ENDIF
    ELSE
      ! Solution of L*y=b
      RowsPerProc = n/Glob_NumOfProcs
      x(1:n) = ZERO
      DO i = 1, RowsPerProc
        ji = (i-1)*Glob_NumOfProcs+1
        jim = ji-1
        jf = jim+Glob_NumOfProcs
        t = ZERO
        jiR = ji+Glob_ProcID
        DO k = 1, jim
          t = t-A(k, jiR)*x(k)
        ENDDO
        x(jiR) = t+b(jiR)
        CALL MPI_ALLREDUCE(x(ji:jf), b(ji:jf), Glob_NumOfProcs, &
                           MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x(ji:jf) = b(ji:jf)
        DO j = ji, jf
          t = x(j)
          DO k = ji, j-1
            t = t-A(k, j)*x(k)
          ENDDO
          x(j) = t
        ENDDO
      ENDDO
      mod_n_Glob_NumOfProcs = MOD(n, Glob_NumOfProcs)
      IF (mod_n_Glob_NumOfProcs > 0) THEN
        jim = RowsPerProc*Glob_NumOfProcs
        ji = jim+1
        jiR = ji+Glob_ProcID
        IF (jiR <= n) THEN
          t = ZERO
          DO k = 1, jim
            t = t-A(k, jiR)*x(k)
          ENDDO
          x(jiR) = t+b(jiR)
        ENDIF
        CALL MPI_ALLREDUCE(x(ji:n), b(ji:n), mod_n_Glob_NumOfProcs, &
                           MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x(ji:n) = b(ji:n)
        DO j = ji, n
          t = x(j)
          DO k = ji, j-1
            t = t-A(k, j)*x(k)
          ENDDO
          x(j) = t
        ENDDO
      ENDIF
      ! Solution of D*LT*x=y
      b(1:n) = x(1:n)*invD(1:n)
      x(1:n) = ZERO
      DO i = 1, RowsPerProc
        jf = n-i*Glob_NumOfProcs+1
        ji = jf-1+Glob_NumOfProcs
        jiR = ji+1+Glob_ProcID
        IF (jiR <= n) THEN
          t = x(jiR)
          DO k = ji, 1, -1
            x(k) = x(k)-t*A(k, jiR)
          ENDDO
          CALL MPI_ALLREDUCE(x(jf:ji), b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs), Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          x(jf:ji) = b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        ENDIF
        x(jf:ji) = x(jf:ji)+b(jf:ji)
        DO j = ji, jf, -1
          t = x(j)
          DO k = j-1, jf, -1
            x(k) = x(k)-t*A(k, j)
          ENDDO
        ENDDO
      ENDDO
      IF (mod_n_Glob_NumOfProcs > 0) THEN
        ji = mod_n_Glob_NumOfProcs
        jiR = ji+1+Glob_ProcID
        IF (jiR <= n) THEN
          t = x(jiR)
          DO k = ji, 1, -1
            x(k) = x(k)-t*A(k, jiR)
          ENDDO
          CALL MPI_ALLREDUCE(x(1:ji), b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs), mod_n_Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          x(1:ji) = b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        ENDIF
        x(1:ji) = x(1:ji)+b(1:ji)
        DO j = ji, 1, -1
          t = x(j)
          DO k = j-1, 1, -1
            x(k) = x(k)-t*A(k, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE LDLTS

  SUBROUTINE LDLHS(n, A, nA, invD, b, x)
    ! Solves A*x = b from the factorization A = L*D*L^H of LDLHF (L*y = b, then D*L^H*x = y).
    !   A(nA,n)  L^H in the upper triangle (unit diagonal implied);   invD  1/D_ii
    !   b        right-hand side, destroyed in the parallel branch;   x  solution

    ! Arguments
    INTEGER     :: n, nA
    COMPLEX(wp) :: A(nA, n), invD(n), b(n), x(n)
    ! Local variables
    INTEGER     :: i, j, k, ji, jf, jim, jfm, RowsPerProc, mod_n_Glob_NumOfProcs, jiR
    COMPLEX(wp) :: t

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'LDLHS: n=', n, ' PMode=', Glob_LDLHS_PMode
    ENDIF
    IF (Glob_LDLHS_PMode == 0) THEN
      IF (Glob_LDLHS_UseBLAS == 0) THEN
        ! Solution of L*y=b
        x(1:n) = b(1:n)
        DO j = 1, n
          t = x(j)
          DO i = 1, j-1
            t = t-CONJG(A(i, j))*x(i)
          ENDDO
          x(j) = t
        ENDDO
        ! Solution of D*LH*x=y
        x(1:n) = x(1:n)*invD(1:n)
        DO j = N, 1, -1
          t = x(j)
          DO i = j-1, 1, -1
            x(i) = x(i)-t*A(i, j)
          ENDDO
        ENDDO
      ELSE
        x(1:n) = b(1:n)
        ! call BLAS routine ZTRSV
        CALL ZTRSV('U', 'C', 'U', n, A, nA, x, 1)
        x(1:n) = x(1:n)*invD(1:n)
        ! call BLAS routine ZTRSV
        CALL ZTRSV('U', 'N', 'U', n, A, nA, x, 1)
      ENDIF
    ELSE
      ! Solution of L*y=b
      RowsPerProc = n/Glob_NumOfProcs
      x(1:n) = ZERO
      DO i = 1, RowsPerProc
        ji = (i-1)*Glob_NumOfProcs+1
        jim = ji-1
        jf = jim+Glob_NumOfProcs
        t = ZERO
        jiR = ji+Glob_ProcID
        DO k = 1, jim
          t = t-CONJG(A(k, jiR))*x(k)
        ENDDO
        x(jiR) = t+b(jiR)
        CALL MPI_ALLREDUCE(x(ji:jf), b(ji:jf), 2*Glob_NumOfProcs, &
                           MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x(ji:jf) = b(ji:jf)
        DO j = ji, jf
          t = x(j)
          DO k = ji, j-1
            t = t-CONJG(A(k, j))*x(k)
          ENDDO
          x(j) = t
        ENDDO
      ENDDO
      mod_n_Glob_NumOfProcs = MOD(n, Glob_NumOfProcs)
      IF (mod_n_Glob_NumOfProcs > 0) THEN
        jim = RowsPerProc*Glob_NumOfProcs
        ji = jim+1
        jiR = ji+Glob_ProcID
        IF (jiR <= n) THEN
          t = ZERO
          DO k = 1, jim
            t = t-CONJG(A(k, jiR))*x(k)
          ENDDO
          x(jiR) = t+b(jiR)
        ENDIF
        CALL MPI_ALLREDUCE(x(ji:n), b(ji:n), 2*mod_n_Glob_NumOfProcs, &
                           MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
        x(ji:n) = b(ji:n)
        DO j = ji, n
          t = x(j)
          DO k = ji, j-1
            t = t-CONJG(A(k, j))*x(k)
          ENDDO
          x(j) = t
        ENDDO
      ENDIF
      ! Solution of D*LH*x=y
      b(1:n) = x(1:n)*invD(1:n)
      x(1:n) = ZERO
      DO i = 1, RowsPerProc
        jf = n-i*Glob_NumOfProcs+1
        ji = jf-1+Glob_NumOfProcs
        jiR = ji+1+Glob_ProcID
        IF (jiR <= n) THEN
          t = x(jiR)
          DO k = ji, 1, -1
            x(k) = x(k)-t*A(k, jiR)
          ENDDO
          CALL MPI_ALLREDUCE(x(jf:ji), b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs), 2*Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          x(jf:ji) = b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        ENDIF
        x(jf:ji) = x(jf:ji)+b(jf:ji)
        DO j = ji, jf, -1
          t = x(j)
          DO k = j-1, jf, -1
            x(k) = x(k)-t*A(k, j)
          ENDDO
        ENDDO
      ENDDO
      IF (mod_n_Glob_NumOfProcs > 0) THEN
        ji = mod_n_Glob_NumOfProcs
        jiR = ji+1+Glob_ProcID
        IF (jiR <= n) THEN
          t = x(jiR)
          DO k = ji, 1, -1
            x(k) = x(k)-t*A(k, jiR)
          ENDDO
          CALL MPI_ALLREDUCE(x(1:ji), b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs), 2*mod_n_Glob_NumOfProcs, &
                             MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
          x(1:ji) = b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        ENDIF
        x(1:ji) = x(1:ji)+b(1:ji)
        DO j = ji, 1, -1
          t = x(j)
          DO k = j-1, 1, -1
            x(k) = x(k)-t*A(k, j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE LDLHS

  SUBROUTINE MTMVL(n, A, nA, x, y, w)
    ! y = A*x for a real symmetric A given by its lower triangle (leading dimension nA); w(n) work.

    INTEGER  :: n, nA, j, jm, ji, jf, p, i
    REAL(wp) :: A(nA, n), x(n), y(n), w(n), t, s

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'MTMVL: n=', n, ' PMode=', Glob_MTMVL_PMode
    ENDIF
    IF (Glob_MTMVL_PMode == 0) THEN
      IF (Glob_MTMVL_UseBLAS == 0) THEN
        IF (wp == 8) THEN
          ! wp=8: one fused pass over the lower triangle (each column serves the dot product for y(j)
          ! and the axpy for y(j+1:n)); wp=10/16: the two-sweep code below is faster.
          y(1:n) = ZERO
          DO j = 1, n
            t = x(j)
            s = y(j)+A(j, j)*t
            DO i = j+1, n
              s = s+A(i, j)*x(i)
              y(i) = y(i)+A(i, j)*t
            ENDDO
            y(j) = s
          ENDDO
        ELSE
          ! We use the property A*x=(R+L)*x = (x^T*R^T)^T + (x^T*L^T)^T
          ! where L is the lower triangle of A (without the diagonal) and
          ! R is the upper triangle (including the diagonal).
          ! Computing (x^T*R^T)^T
          DO j = 1, n
            t = ZERO
            DO i = j, n
              t = t+A(i, j)*x(i)
            ENDDO
            y(j) = t
          ENDDO
          ! Computing (x^T*L^T)^T
          DO j = n, 1, -1
            t = x(j)
            DO i = j+1, n
              y(i) = y(i)+t*A(i, j)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        ! call BLAS routine DSYMV
        CALL DSYMV('L', n, ONE, A, nA, x, 1, ZERO, y, 1)
      ENDIF
    ELSE
      w(1:n) = ZERO
      p = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) p = p+1
      ! Computing (x^T*R^T)^T
      ji = 1+p*Glob_ProcID
      jf = MIN(p*(Glob_ProcID+1), n)
      DO j = ji, jf
        t = ZERO
        DO i = j, n
          t = t+A(i, j)*x(i)
        ENDDO
        w(j) = t
      ENDDO
      ! Computing (x^T*L^T)^T
      ji = 1+p*(Glob_NumOfProcs-Glob_ProcID-1)
      jf = MIN(p*(Glob_NumOfProcs-Glob_ProcID), n)
      DO j = jf, ji, -1
        t = x(j)
        DO i = j+1, n
          w(i) = w(i)+t*A(i, j)
        ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(w, y, n, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF

  END SUBROUTINE MTMVL

  SUBROUTINE MHMVL(n, A, nA, x, y, w)
    ! y = A*x for a complex hermitian A given by its lower triangle (leading dimension nA); w(n) work.

    INTEGER     :: n, nA, j, jm, ji, jf, p, i
    COMPLEX(wp) :: A(nA, n), x(n), y(n), w(n), t, s

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'MHMVL: n=', n, ' PMode=', Glob_MHMVL_PMode
    ENDIF
    IF (Glob_MHMVL_PMode == 0) THEN
      IF (Glob_MHMVL_UseBLAS == 0) THEN
        IF (wp == 8) THEN
          ! wp=8: one fused pass over the lower triangle (conjugated dot product for y(j) and axpy
          ! for y(j+1:n)); wp=10/16: the two-sweep code below is faster.
          y(1:n) = CMPLX(ZERO, ZERO, wp)
          DO j = 1, n
            t = x(j)
            s = y(j)+CONJG(A(j, j))*t
            DO i = j+1, n
              s = s+CONJG(A(i, j))*x(i)
              y(i) = y(i)+A(i, j)*t
            ENDDO
            y(j) = s
          ENDDO
        ELSE
          ! We use the property A*x=(R+L)*x = (x^H*R^H)^H + (x^H*L^H)^H
          ! where L is the lower triangle of A (without the diagonal) and
          ! R is the upper triangle (including the diagonal).
          ! Computing (x^H*R^H)^H
          DO j = 1, n
            t = CMPLX(ZERO, ZERO, wp)
            DO i = j, n
              t = t+CONJG(A(i, j))*x(i)
            ENDDO
            y(j) = t
          ENDDO
          ! Computing (x^H*L^H)^H
          DO j = n, 1, -1
            t = x(j)
            DO i = j+1, n
              y(i) = y(i)+t*A(i, j)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        ! call BLAS routine ZHEMV
        CALL ZHEMV('L', n, CMPLX(ONE, ZERO, wp), A, nA, x, 1, CMPLX(ZERO, ZERO, wp), y, 1)
      ENDIF
    ELSE
      w(1:n) = ZERO
      p = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) p = p+1
      ! Computing (x^H*R^H)^H
      ji = 1+p*Glob_ProcID
      jf = MIN(p*(Glob_ProcID+1), n)
      DO j = ji, jf
        t = CMPLX(ZERO, ZERO, wp)
        DO i = j, n
          t = t+CONJG(A(i, j))*x(i)
        ENDDO
        w(j) = t
      ENDDO
      ! Computing (x^H*L^H)^H
      ji = 1+p*(Glob_NumOfProcs-Glob_ProcID-1)
      jf = MIN(p*(Glob_NumOfProcs-Glob_ProcID), n)
      DO j = jf, ji, -1
        t = x(j)
        DO i = j+1, n
          w(i) = w(i)+t*A(i, j)
        ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(w, y, 2*n, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF

  END SUBROUTINE MHMVL

  SUBROUTINE MTMV(n, A, nA, x, y, w)
    ! y = A^T*x for a real matrix A (leading dimension nA); w(n) work.

    INTEGER  :: n, nA, ib, ie, jb, je, j, i, k
    REAL(wp) :: A(nA, n), x(n), y(n), w(n), t
    INTEGER  :: n2, p, q

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'MTMV: n=', n, ' PMode=', Glob_MTMV_PMode
    ENDIF
    IF (Glob_MTMV_PMode == 0) THEN
      IF (Glob_MTMV_UseBLAS == 0) THEN
        DO i = 1, n
          t = ZERO
          DO j = 1, n
            t = t+A(j, i)*x(j)
          ENDDO
          y(i) = t
        ENDDO
      ELSE
        ! call BLAS routine DGEMV
        CALL DGEMV('T', n, n, ONE, A, nA, x, 1, ZERO, y, 1)
      ENDIF
    ELSE
      ! Parallel branch: the elements A(1,1)..A(n,n) in row order are split into contiguous chunks
      ! (ib,jb)..(ie,je) of k elements, one per process, as evenly as possible.
      w(1:n) = ZERO
      n2 = n*n
      p = n2/Glob_NumOfProcs
      IF (MOD(n2, Glob_NumOfProcs) /= 0) p = p+1
      q = Glob_ProcID*p
      IF (q+p > n2) THEN
        k = MAX(n2-q, 0)
      ELSE
        k = p
      ENDIF
      ib = q/n+1
      ie = MIN((q+p-1)/n+1, n)
      jb = MOD(q, n)+1
      je = MOD(q+k-1, n)+1
      IF (ie > ib) THEN
        t = ZERO
        DO j = jb, n
          t = t+A(j, ib)*x(j)
        ENDDO
        w(ib) = t
        DO i = ib+1, ie-1
          t = ZERO
          DO j = 1, n
            t = t+A(j, i)*x(j)
          ENDDO
          w(i) = t
        ENDDO
        t = ZERO
        DO j = 1, je
          t = t+A(j, i)*x(j)
        ENDDO
        w(ie) = t
      ELSE
        IF (ie == ib) THEN
          t = ZERO
          DO j = jb, je
            t = t+A(j, ib)*x(j)
          ENDDO
          w(ib) = t
        ENDIF
      ENDIF
      CALL MPI_ALLREDUCE(w, y, n, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END SUBROUTINE MTMV

  SUBROUTINE MHMV(n, A, nA, x, y, w)
    ! y = A^H*x for a complex matrix A (leading dimension nA); w(n) work.

    INTEGER     :: n, nA, ib, ie, jb, je, j, i, k
    COMPLEX(wp) :: A(nA, n), x(n), y(n), w(n), t
    INTEGER     :: n2, p, q

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'MHMV: n=', n, ' PMode=', Glob_MHMV_PMode
    ENDIF
    IF (Glob_MHMV_PMode == 0) THEN
      IF (Glob_MHMV_UseBLAS == 0) THEN
        DO i = 1, n
          t = CMPLX(ZERO, ZERO, wp)
          DO j = 1, n
            t = t+CONJG(A(j, i))*x(j)
          ENDDO
          y(i) = t
        ENDDO
      ELSE
        ! call BLAS routine ZGEMV
        CALL ZGEMV('C', n, n, CMPLX(ONE, ZERO, wp), A, nA, x, 1, CMPLX(ZERO, ZERO, wp), y, 1)
      ENDIF
    ELSE
      ! Parallel branch: the elements A(1,1)..A(n,n) in row order are split into contiguous chunks
      ! (ib,jb)..(ie,je) of k elements, one per process, as evenly as possible.
      w(1:n) = CMPLX(ZERO, ZERO, wp)
      n2 = n*n
      p = n2/Glob_NumOfProcs
      IF (MOD(n2, Glob_NumOfProcs) /= 0) p = p+1
      q = Glob_ProcID*p
      IF (q+p > n2) THEN
        k = MAX(n2-q, 0)
      ELSE
        k = p
      ENDIF
      ib = q/n+1
      ie = MIN((q+p-1)/n+1, n)
      jb = MOD(q, n)+1
      je = MOD(q+k-1, n)+1
      IF (ie > ib) THEN
        t = CMPLX(ZERO, ZERO, wp)
        DO j = jb, n
          t = t+CONJG(A(j, ib))*x(j)
        ENDDO
        w(ib) = t
        DO i = ib+1, ie-1
          t = CMPLX(ZERO, ZERO, wp)
          DO j = 1, n
            t = t+CONJG(A(j, i))*x(j)
          ENDDO
          w(i) = t
        ENDDO
        t = CMPLX(ZERO, ZERO, wp)
        DO j = 1, je
          t = t+CONJG(A(j, i))*x(j)
        ENDDO
        w(ie) = t
      ELSE
        IF (ie == ib) THEN
          t = CMPLX(ZERO, ZERO, wp)
          DO j = jb, je
            t = t+CONJG(A(j, ib))*x(j)
          ENDDO
          w(ib) = t
        ENDIF
      ENDIF
      CALL MPI_ALLREDUCE(w, y, 2*n, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END SUBROUTINE MHMV

  FUNCTION VMMTMV(n, A, nA, x)
    ! x^T*A*x for a real symmetric A given by its lower triangle (leading dimension nA).

    ! Arguments
    INTEGER  :: n, nA
    REAL(wp) :: A(nA, n), x(n)
    REAL(wp) :: VMMTMV
    ! Local variables
    INTEGER  :: j, k, q, p
    REAL(wp) :: s, sum

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'VMMTMV: n=', n, ' PMode=', Glob_VMMTMV_PMode
    ENDIF
    IF (Glob_VMMTMV_PMode == 0) THEN
      VMMTMV = ZERO
      DO k = 1, n
        VMMTMV = VMMTMV+x(k)*x(k)*A(k, k)
        s = ZERO
        DO j = k+1, n
          s = s+x(j)*A(j, k)
        ENDDO
        VMMTMV = VMMTMV+TWO*s*x(k)
      ENDDO
    ELSE
      sum = ZERO
      q = n/(2*Glob_NumOfProcs)
      DO p = 1, 2*q, 2
        k = (p-1)*Glob_NumOfProcs+Glob_ProcID+1
        sum = sum+x(k)*x(k)*A(k, k)
        s = ZERO
        DO j = k+1, n
          s = s+x(j)*A(j, k)
        ENDDO
        sum = sum+TWO*s*x(k)
        k = (p+1)*Glob_NumOfProcs-Glob_ProcID
        sum = sum+x(k)*x(k)*A(k, k)
        s = ZERO
        DO j = k+1, n
          s = s+x(j)*A(j, k)
        ENDDO
        sum = sum+TWO*s*x(k)
      ENDDO
      p = Glob_ProcID
      DO k = 2*q*Glob_NumOfProcs+1, n
        IF (p == 0) THEN
          sum = sum+x(k)*x(k)*A(k, k)
          p = p+Glob_NumOfProcs
        ENDIF
        s = ZERO
        DO j = k+p, n, Glob_NumOfProcs
          s = s+x(j)*A(j, k)
        ENDDO
        sum = sum+TWO*s*x(k)
        p = j-n-1
      ENDDO
      CALL MPI_ALLREDUCE(sum, VMMTMV, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF

  END FUNCTION VMMTMV

  FUNCTION VMMHMV(n, A, nA, x)
    ! x^H*A*x for a complex hermitian A given by its lower triangle (leading dimension nA).

    ! Arguments
    INTEGER     :: n, nA
    COMPLEX(wp) :: A(nA, n), x(n)
    REAL(wp)    :: VMMHMV
    ! Local variables
    INTEGER  :: j, k, q, p
    REAL(wp) :: sr, si, sum

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'VMMHMV: n=', n, ' PMode=', Glob_VMMHMV_PMode
    ENDIF
    IF (Glob_VMMHMV_PMode == 0) THEN
      VMMHMV = ZERO
      DO k = 1, n
        VMMHMV = VMMHMV+(REAL(x(k), wp)*REAL(x(k), wp)+imag(x(k))*imag(x(k)))*REAL(A(k, k), wp)
        sr = ZERO
        si = ZERO
        DO j = k+1, n
          sr = sr+(REAL(x(j), wp)*REAL(A(j, k), wp)+imag(x(j))*imag(A(j, k)))
          si = si+(imag(x(j))*REAL(A(j, k), wp)-REAL(x(j), wp)*imag(A(j, k)))
        ENDDO
        VMMHMV = VMMHMV+TWO*(REAL(x(k), wp)*sr+imag(x(k))*si)
      ENDDO
    ELSE
      sum = ZERO
      q = n/(2*Glob_NumOfProcs)
      DO p = 1, 2*q, 2
        k = (p-1)*Glob_NumOfProcs+Glob_ProcID+1
        sum = sum+(REAL(x(k), wp)*REAL(x(k), wp)+imag(x(k))*imag(x(k)))*REAL(A(k, k), wp)
        sr = ZERO
        si = ZERO
        DO j = k+1, n
          sr = sr+(REAL(x(j), wp)*REAL(A(j, k), wp)+imag(x(j))*imag(A(j, k)))
          si = si+(imag(x(j))*REAL(A(j, k), wp)-REAL(x(j), wp)*imag(A(j, k)))
        ENDDO
        sum = sum+TWO*(REAL(x(k), wp)*sr+imag(x(k))*si)
        k = (p+1)*Glob_NumOfProcs-Glob_ProcID
        sum = sum+(REAL(x(k), wp)*REAL(x(k), wp)+imag(x(k))*imag(x(k)))*REAL(A(k, k), wp)
        sr = ZERO
        si = ZERO
        DO j = k+1, n
          sr = sr+(REAL(x(j), wp)*REAL(A(j, k), wp)+imag(x(j))*imag(A(j, k)))
          si = si+(imag(x(j))*REAL(A(j, k), wp)-REAL(x(j), wp)*imag(A(j, k)))
        ENDDO
        sum = sum+TWO*(REAL(x(k), wp)*sr+imag(x(k))*si)
      ENDDO
      p = Glob_ProcID
      DO k = 2*q*Glob_NumOfProcs+1, n
        IF (p == 0) THEN
          sum = sum+(REAL(x(k), wp)*REAL(x(k), wp)+imag(x(k))*imag(x(k)))*REAL(A(k, k), wp)
          p = p+Glob_NumOfProcs
        ENDIF
        sr = ZERO
        si = ZERO
        DO j = k+p, n, Glob_NumOfProcs
          sr = sr+(REAL(x(j), wp)*REAL(A(j, k), wp)+imag(x(j))*imag(A(j, k)))
          si = si+(imag(x(j))*REAL(A(j, k), wp)-REAL(x(j), wp)*imag(A(j, k)))
        ENDDO
        sum = sum+TWO*(REAL(x(k), wp)*sr+imag(x(k))*si)
        p = j-n-1
      ENDDO
      CALL MPI_ALLREDUCE(sum, VMMHMV, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF

  END FUNCTION VMMHMV

  FUNCTION RMaxAbsEl(n, x)
    ! Function RMaxAbsEl returns the magnitude of the largest by magnitude element
    ! among the first n elements of real array x
    INTEGER  :: n, j, ji, jf, k
    REAL(wp) :: x(n)
    REAL(wp) :: RMaxAbsEl, MaxAE
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RMaxAbsEl: n=', n
    ENDIF
    MaxAE = ZERO
    IF (Glob_RMaxAbsEl_PMode == 0) THEN
      DO j = 1, n
        IF (ABS(x(j)) > MaxAE) MaxAE = ABS(x(j))
      ENDDO
      RMaxAbsEl = MaxAE
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      DO j = ji, jf
        IF (ABS(x(j)) > MaxAE) MaxAE = ABS(x(j))
      ENDDO
      CALL MPI_ALLREDUCE(MaxAE, RMaxAbsEl, 1, MPI_WP, MPI_MAX, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION RMaxAbsEl

  FUNCTION CMaxAbsReOrIm(n, x)
    ! Function CMaxAbsEl returns the magnitude of the largest real or imaginary part
    ! among the first n elements of complex array x
    INTEGER     :: n, j, ji, jf, k, maxj
    COMPLEX(wp) :: x(n)
    REAL(wp)    :: CMaxAbsReOrIm, MaxAE, t
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CMaxAbsReOrIm: n=', n
    ENDIF
    MaxAE = ZERO
    IF (Glob_CMaxAbsReOrIm_PMode == 0) THEN
      DO j = 1, n
        t = ABS(REAL(x(j), wp))
        IF (t > MaxAE) MaxAE = t
        t = ABS(imag(x(j)))
        IF (t > MaxAE) MaxAE = t
      ENDDO
      CMaxAbsReOrIm = MaxAE
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      DO j = ji, jf
        t = ABS(REAL(x(j), wp))
        IF (t > MaxAE) MaxAE = t
        t = ABS(imag(x(j)))
        IF (t > MaxAE) MaxAE = t
      ENDDO
      CALL MPI_ALLREDUCE(MaxAE, CMaxAbsReOrIm, 1, MPI_WP, MPI_MAX, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION CMaxAbsReOrIm

  FUNCTION RDotProd(n, x, y)
    ! Function RDotProd computes the dot product x^{T}y,
    ! for n-component real vectors x and y.
    INTEGER  :: n, k, ji, jf, j
    REAL(wp) :: RDotProd, a
    REAL(wp) :: x(n), y(n)
    REAL(wp) :: DDOT
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RDotProd: n=', n
    ENDIF
    IF (Glob_RDotProd_PMode == 0) THEN
      IF (Glob_RDotProd_UseBLAS == 0) THEN
        RDotProd = ZERO
        DO j = 1, n
          RDotProd = RDotProd+x(j)*y(j)
        ENDDO
      ELSE
        ! call BLAS function DDOT
        RDotProd = DDOT(n, x, 1, y, 1)
      ENDIF
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      a = ZERO
      DO j = ji, jf
        a = a+x(j)*y(j)
      ENDDO
      CALL MPI_ALLREDUCE(a, RDotProd, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION RDotProd

  FUNCTION CDotProd(n, x, y)
    ! Function CDotProd computes the dot product x^{H}y,
    ! for n-component complex vectors x and y.
    INTEGER     :: n, k, ji, jf, j
    COMPLEX(wp) :: CDotProd, a
    COMPLEX(wp) :: x(n), y(n)
    COMPLEX(wp) :: ZDOTC
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CDotProd: n=', n
    ENDIF
    IF (Glob_CDotProd_PMode == 0) THEN
      IF (Glob_CDotProd_UseBLAS == 0) THEN
        CDotProd = CMPLX(ZERO, ZERO, wp)
        DO j = 1, n
          CDotProd = CDotProd+CMPLX(REAL(x(j), wp)*REAL(y(j), wp)+imag(x(j))*imag(y(j)), &
                                    REAL(x(j), wp)*imag(y(j))-imag(x(j))*REAL(y(j), wp), wp)
        ENDDO
      ELSE
        ! call BLAS function ZDOTC
        CDotProd = ZDOTC(n, x, 1, y, 1)
      ENDIF
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      a = CMPLX(ZERO, ZERO, wp)
      DO j = ji, jf
        a = a+CMPLX(REAL(x(j), wp)*REAL(y(j), wp)+imag(x(j))*imag(y(j)), &
                    REAL(x(j), wp)*imag(y(j))-imag(x(j))*REAL(y(j), wp), wp)
      ENDDO
      CALL MPI_ALLREDUCE(a, CDotProd, 2, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION CDotProd

  FUNCTION RDotProdItself(n, x)
    ! Function RDotProdItself computes dot product x^{T}x,
    ! for an n-component real vector x.
    INTEGER  :: n, k, ji, jf, j
    REAL(wp) :: RDotProdItself, t
    REAL(wp) :: x(n)
    REAL(wp) :: DDOT
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RDotProdItself: n=', n
    ENDIF
    IF (Glob_RDotProdItself_PMode == 0) THEN
      IF (Glob_RDotProdItself_UseBLAS == 0) THEN
        RDotProdItself = ZERO
        DO j = 1, n
          RDotProdItself = RDotProdItself+x(j)*x(j)
        ENDDO
      ELSE
        ! call BLAS function DDOT
        RDotProdItself = DDOT(n, x, 1, x, 1)
      ENDIF
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      t = ZERO
      DO j = ji, jf
        t = t+x(j)*x(j)
      ENDDO
      CALL MPI_ALLREDUCE(t, RDotProdItself, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION RDotProdItself

  FUNCTION CDotProdItself(n, x)
    ! Function CDotProdItself computes dot product x^{H}x,
    ! for an n-component complex vector x.
    INTEGER     :: n, k, ji, jf, j
    REAL(wp)    :: CDotProdItself, t
    COMPLEX(wp) :: x(n)
    COMPLEX(wp) :: ZDOTC
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CDotProdItself: n=', n
    ENDIF
    IF (Glob_CDotProdItself_PMode == 0) THEN
      IF (Glob_CDotProdItself_UseBLAS == 0) THEN
        CDotProdItself = ZERO
        DO j = 1, n
          CDotProdItself = CDotProdItself+REAL(x(j), wp)*REAL(x(j), wp)+imag(x(j))*imag(x(j))
        ENDDO
      ELSE
        ! call BLAS function ZDOTC
        CDotProdItself = ZDOTC(n, x, 1, x, 1)
      ENDIF
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      t = ZERO
      DO j = ji, jf
        t = t+REAL(x(j), wp)*REAL(x(j), wp)+imag(x(j))*imag(x(j))
      ENDDO
      CALL MPI_ALLREDUCE(t, CDotProdItself, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
    ENDIF
  END FUNCTION CDotProdItself

  FUNCTION RDotProdQuotient(n, x, y)
    ! Function RDotProdQuotient computes the quotient (x^{T}y)/(y^{T}y),
    ! where x and y are n-component real vectors.
    INTEGER  :: n, k, ji, jf, j
    REAL(wp) :: RDotProdQuotient, a, b, t(2), tt(2)
    REAL(wp) :: x(n), y(n)
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RDotProdQuotient: n=', n
    ENDIF
    IF (Glob_RDotProdQuotient_PMode == 0) THEN
      a = ZERO
      b = ZERO
      DO j = 1, n
        a = a+x(j)*y(j)
        b = b+y(j)*y(j)
      ENDDO
      RDotProdQuotient = a/b
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      t(1) = ZERO
      t(2) = ZERO
      DO j = ji, jf
        t(1) = t(1)+x(j)*y(j)
        t(2) = t(2)+y(j)*y(j)
      ENDDO
      CALL MPI_ALLREDUCE(t, tt, 2, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      RDotProdQuotient = tt(1)/tt(2)
    ENDIF
  END FUNCTION RDotProdQuotient

  FUNCTION CDotProdQuotient(n, x, y)
    ! Function CDotProdQuotient computes the quotient (x^{H}y)/(y^{H}y),
    ! where x and y are n-component complex vectors.
    INTEGER     :: n, k, ji, jf, j
    COMPLEX(wp) :: CDotProdQuotient
    COMPLEX(wp) :: x(n), y(n)
    COMPLEX(wp) :: a
    REAL(wp)    :: r, t(3), tt(3)
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CDotProdQuotient: n=', n
    ENDIF
    IF (Glob_CDotProdQuotient_PMode == 0) THEN
      a = CMPLX(ZERO, ZERO, wp)
      r = ZERO
      DO j = 1, n
        a = a+CONJG(x(j))*y(j)
        r = r+REAL(y(j), wp)*REAL(y(j), wp)+imag(y(j))*imag(y(j))
      ENDDO
      CDotProdQuotient = CMPLX(REAL(a, wp)/r, imag(a)/r, wp)
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      a = CMPLX(ZERO, ZERO, wp)
      t(1) = ZERO
      DO j = ji, jf
        a = a+CONJG(x(j))*y(j)
        t(1) = t(1)+REAL(y(j), wp)*REAL(y(j), wp)+imag(y(j))*imag(y(j))
      ENDDO
      t(2) = REAL(a, wp)
      t(3) = imag(a)
      CALL MPI_ALLREDUCE(t, tt, 3, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      CDotProdQuotient = CMPLX(tt(2)/tt(1), tt(3)/tt(1), wp)
    ENDIF
  END FUNCTION CDotProdQuotient

  SUBROUTINE RVScale(n, x, alpha)
    ! Subroutine RVScale scales real vector x by real constant alpha
    INTEGER  :: n, j
    REAL(wp) :: x(n)
    REAL(wp) :: alpha

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RVScale: n=', n
    ENDIF
    IF (Glob_RVScale_UseBLAS == 0) THEN
      DO j = 1, n
        x(j) = alpha*x(j)
      ENDDO
    ELSE
      ! call BLAS routine DSCAL
      CALL DSCAL(n, alpha, x, 1)
    ENDIF

  END SUBROUTINE RVScale

  SUBROUTINE CVScale(n, x, alpha)
    ! Subroutine RVScale scales complex vector x by complex constant alpha
    INTEGER     :: n, j
    COMPLEX(wp) :: x(n)
    COMPLEX(wp) :: alpha
    REAL(wp)    :: ralpha

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CVScale: n=', n
    ENDIF
    IF (Glob_CVScale_UseBLAS == 0) THEN
      IF (imag(alpha) == ZERO) THEN
        ralpha = REAL(alpha, wp)
        DO j = 1, n
          x(j) = CMPLX(ralpha*REAL(x(j), wp), ralpha*imag(x(j)), wp)
        ENDDO
      ELSE
        DO j = 1, n
          x(j) = x(j)*alpha
        ENDDO
      ENDIF
    ELSE
      ! call BLAS routine ZSCAL
      CALL ZSCAL(n, alpha, x, 1)
    ENDIF

  END SUBROUTINE CVScale

  FUNCTION RVDiffEucNorm(n, x, alpha, y)
    ! Function RVDiffEucNorm returns the Euclidean norm of the difference x-alpha*y,
    ! where y and x are two real vectors and alpha is a real scalar.
    REAL(wp) :: RVDiffEucNorm
    INTEGER  :: n, j, k, ji, jf, q
    REAL(wp) :: x(n), y(n), alpha, t

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'RVDiffEucNorm: n=', n
    ENDIF
    IF (Glob_RVDiffEucNorm_PMode == 0) THEN
      RVDiffEucNorm = ZERO
      DO j = 1, n
        RVDiffEucNorm = RVDiffEucNorm+(x(j)-alpha*y(j))*(x(j)-alpha*y(j))
      ENDDO
      RVDiffEucNorm = SQRT(RVDiffEucNorm)
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      t = ZERO
      DO j = ji, jf
        t = t+(x(j)-alpha*y(j))*(x(j)-alpha*y(j))
      ENDDO
      CALL MPI_ALLREDUCE(t, RVDiffEucNorm, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      RVDiffEucNorm = SQRT(RVDiffEucNorm)
    ENDIF

  END FUNCTION RVDiffEucNorm

  FUNCTION CVDiffEucNorm(n, x, alpha, y)
    ! Function CVDiffEucNorm returns the Euclidean norm of the difference x-alpha*y,
    ! where y and x are two complex vectors and alpha is a complex scalar.
    REAL(wp)    :: CVDiffEucNorm
    INTEGER     :: n, j, k, ji, jf, q
    COMPLEX(wp) :: x(n), y(n), alpha, s
    REAL(wp)    :: t

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0)') 'CVDiffEucNorm: n=', n
    ENDIF
    IF (Glob_CVDiffEucNorm_PMode == 0) THEN
      CVDiffEucNorm = ZERO
      DO j = 1, n
        s = x(j)-alpha*y(j)
        CVDiffEucNorm = CVDiffEucNorm+REAL(s, wp)*REAL(s, wp)+imag(s)*imag(s)
      ENDDO
      CVDiffEucNorm = SQRT(CVDiffEucNorm)
    ELSE
      k = n/Glob_NumOfProcs
      IF (MOD(n, Glob_NumOfProcs) /= 0) k = k+1
      ji = 1+k*Glob_ProcID
      jf = MIN(k*(Glob_ProcID+1), n)
      t = CMPLX(ZERO, ZERO, wp)
      DO j = ji, jf
        s = x(j)-alpha*y(j)
        t = t+REAL(s, wp)*REAL(s, wp)+imag(s)*imag(s)
      ENDDO
      CALL MPI_ALLREDUCE(t, CVDiffEucNorm, 1, MPI_WP, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      CVDiffEucNorm = SQRT(CVDiffEucNorm)
    ENDIF

  END FUNCTION CVDiffEucNorm

  SUBROUTINE GSEPIIS(k, n, M, nM, invD, B, nB, apprlambda, v, w, Tol, &
                     lambda, x, RelAcc, MaxIter, SpecifNorm, NumIter, ErrorCode)
    ! Inverse iteration for one eigenpair of the real symmetric pencil A*x = lambda*B*x, extending
    ! the factorization M = A - apprlambda*B = L*D*L^T of the leading (k-1) x (k-1) block (LDLTF).
    ! Converges to the eigenvalue nearest apprlambda; eigenvalues closer than about epsilon
    ! cannot be separated (ErrorCode 2).
    !   M(nM,n)     upper triangle L^T, lower triangle A - apprlambda*B (sizes k-1 in, n out);  invD 1/D_ii
    !   B(nB,n)     overlap matrix
    !   apprlambda  shift: closer to the wanted eigenvalue than to any other and not equal to it
    !   v           start vector (destroyed);   w  work
    !   Tol         relative accuracy; negative: best attainable, at least |Tol| (one extra iteration)
    !   MaxIter     iteration limit;   SpecifNorm  0: x^T*B*x = 1, 1: x^T*x = 1, else max|x_i| = 1
    !   lambda, x   the eigenpair;   RelAcc  accuracy reached;   NumIter  iterations used
    !   ErrorCode   0 ok, 1 M singular, 2 Tol not reached in MaxIter, 3 x is a null direction of B

    ! Arguments
    INTEGER  :: k, n, nM, nB, MaxIter, SpecifNorm, NumIter, ErrorCode
    REAL(wp) :: M(nM, n), invD(n), B(nB, n), v(n), w(n), x(n)
    REAL(wp) :: apprlambda, Tol, lambda, RelAcc
    ! Local variables
    REAL(wp) :: NormOfDiff, NormOfDiffPrev, t1, t2, sqrtn
    REAL(wp) :: tc
    LOGICAL  :: notconverged
    ! Solver tracing at Verbose >= 3 (see globvars)
    REAL(wp) :: tr_maxx, tr_invDmin, tr_invDmax, tr_prev
    INTEGER  :: tr_i, tr_nrise

    ErrorCode = 0
    RelAcc = HUGE(RelAcc)
    NormOfDiffPrev = HUGE(NormOfDiffPrev)
    NumIter = 0
    ! Updating M=L*D*LH factorization up to size n using routine LDLHF
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, *)
      WRITE(*, *) 'GSEPIIS (I solver)  ------------------------------------------'
      WRITE(*, '(a20,i6,a3,i6)') '           k / n ', k, ' / ', n
      WRITE(*, '(a20,es11.4)') '           shift ', apprlambda
      WRITE(*, '(a20,es11.4,a3,i6)') '   tol / maxiter ', Tol, ' / ', MaxIter
    ENDIF
    CALL LDLTF(k, n, M, nM, invD, w, ErrorCode)

    ! Inertia of D = number of eigenvalues below the shift, taken at every solve; with the sign
    ! of lambda - apprlambda it identifies the eigenvalue returned. Meaningless after a failed
    ! factorization.
    IF (ErrorCode == 0) THEN
      Glob_NumEvalsBelowShift = InertiaCountBelow(n, invD)
    ELSE
      Glob_NumEvalsBelowShift = -1
    ENDIF
    Glob_LastEigIndex = -1
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      tr_invDmin = HUGE(tr_invDmin)
      tr_invDmax = ZERO
      DO tr_i = 1, n
        IF (ABS(invD(tr_i)) < tr_invDmin) tr_invDmin = ABS(invD(tr_i))
        IF (ABS(invD(tr_i)) > tr_invDmax) tr_invDmax = ABS(invD(tr_i))
      ENDDO
      WRITE(*, '(a20,i6)') '      error code ', ErrorCode
      WRITE(*, '(a20,es11.4,a3,es11.4)') ' invD  min / max ', tr_invDmin, ' / ', tr_invDmax
      IF (ErrorCode == 0) WRITE(*, '(a20,i6)') 'evals below shft ', Glob_NumEvalsBelowShift
    ENDIF

    IF (ErrorCode > 0) RETURN
    sqrtn = SQRT(n*ONE)
    tr_nrise = 0
    tr_prev = HUGE(tr_prev)
    ! Do inverse iterations until the process converges
    ! with relative accuracy Tol, or until the number
    ! of iterations exceeds the limit
    notconverged = .TRUE.
    DO WHILE ((notconverged) .AND. (NumIter < MaxIter))
      IF (NumIter /= 0) v(1:n) = x(1:n)
      ! wp=8: B*v with the symmetric routine MTMVL (half the memory traffic); wp=10/16: MTMV is faster.
      IF (wp == 8) THEN
        CALL MTMVL(n, B, nB, v, w, x)
      ELSE
        CALL MTMV(n, B, nB, v, w, x)
      ENDIF
      CALL LDLTS(n, M, nM, invD, w, x)
      t1 = RMaxAbsEl(n, x)
      tr_maxx = t1
      CALL RVScale(n, x, 1/t1)
      ! tc=(x^{T}v)/(v^{T}v)
      tc = RDotProdQuotient(n, x, v)
      ! NormOfDiff=sqrt(n)||x-tc*v||
      ! old version: NormOfDiff=sqrtn*RVDiffEucNorm(n,x,tc,v)
      NormOfDiff = RVDiffEucNorm(n, x, tc, v)/SQRT(RDotProdItself(n, x))
      IF (Tol > ZERO) THEN
        IF (NormOfDiff <= Tol) notconverged = .FALSE.
      ELSE
        IF ((NormOfDiff > NormOfDiffPrev) .AND. (NormOfDiff <= ABS(Tol))) notconverged = .FALSE.
        NormOfDiffPrev = NormOfDiff
      ENDIF
      IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,es11.4,a,es11.4)') '   GSEPIIS: tc = x.v/v.v =', tc, '  previous |dx| =', NormOfDiffPrev
      ENDIF
      NumIter = NumIter+1
      IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(a8,i5,a10,es11.4,a9,es11.4)') '   iter ', NumIter, '   |dx| = ', NormOfDiff, '  max|x|= ', tr_maxx
      ENDIF
      ! Consecutive INCREASES of |dx| mean two eigenvalues straddle the shift at
      ! nearly equal distance and the iterate alternates between their vectors.
      IF (NormOfDiff > tr_prev) THEN
        tr_nrise = tr_nrise+1
      ELSE
        tr_nrise = 0
      ENDIF
      tr_prev = NormOfDiff
    ENDDO
    RelAcc = NormOfDiff
    IF (notconverged) ErrorCode = 2
    ! Compute x^{T}Bx
    t1 = VMMTMV(n, B, nB, x)
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(a20,es11.4,a16)') '            x_Bx ', t1, '  (must be > 0)'
    ENDIF
    ! Rayleigh-quotient guard: for a genuine eigenvector x^T*B*x / x^T*x is bounded below by the
    ! smallest eigenvalue of B. A near-zero value means the iteration converged onto a null
    ! direction of a linearly dependent basis; the quotient would be garbage although the
    ! convergence test reports success. ErrorCode 3; every caller rejects the point.
    IF (t1 <= Glob_MinRayleighDenom*RDotProdItself(n, x)) THEN
      ErrorCode = 3
      lambda = apprlambda
      RelAcc = HUGE(RelAcc)
      IF ((Verbose >= 1) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, *) 'Warning in GSEPIIS: the eigenvector is a null direction of B'
        WRITE(*, *) '  x^T B x                =', t1
        WRITE(*, *) '  x^T x                  =', RDotProdItself(n, x)
        WRITE(*, *) '  Glob_MinRayleighDenom  =', Glob_MinRayleighDenom
        WRITE(*, *) '  the basis is linearly dependent in this direction;'
        WRITE(*, *) '  rejecting instead of returning a Rayleigh quotient'
      ENDIF
      RETURN
    ENDIF
    ! Compute x^{T}Mx
    t2 = VMMTMV(n, M, nM, x)
    ! Rayleigh quotient
    lambda = (t2/t1)+apprlambda

    ! Identify WHICH eigenvalue this is. Inverse iteration converges to the one
    ! nearest the shift, and with m eigenvalues strictly below the shift the two
    ! candidates are lambda_m (below) and lambda_{m+1} (above); the sign of
    ! lambda-apprlambda says which. Exact, and free.
    IF (Glob_NumEvalsBelowShift >= 0) THEN
      IF (lambda > apprlambda) THEN
        Glob_LastEigIndex = Glob_NumEvalsBelowShift+1
      ELSE
        Glob_LastEigIndex = Glob_NumEvalsBelowShift
      ENDIF
    ENDIF
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(a20,es11.4)') '            x_Mx ', t2
      WRITE(*, '(a20,i6,a2,i6)') '      iterations ', NumIter, ' /', MaxIter
      WRITE(*, '(a20,es11.4)') '        residual ', RelAcc
      WRITE(*, '(a20,es11.4)') '          lambda ', lambda
      WRITE(*, '(a20,es11.4)') '  lambda - shift ', lambda-apprlambda
      WRITE(*, '(a20,i6,a12,i6)') '   eigval index  ', Glob_LastEigIndex, '  wanted    ', Glob_WhichEigenvalue
      IF (ErrorCode == 0) THEN
        WRITE(*, '(a20,i6,a14)') '      error code ', ErrorCode, '  (converged)'
      ELSE
        WRITE(*, '(a20,i6,a18)') '      error code ', ErrorCode, '  (NOT converged)'
      ENDIF
      IF (tr_nrise >= 5) THEN
        WRITE(*, *) '  ** DIVERGED - |dx| rose for ', tr_nrise, ' consecutive iterations:'
        WRITE(*, *) '     two eigenvalues straddle the shift at nearly equal distance;'
        WRITE(*, *) '     only a better shift can fix that (see RefreshINVITShift)'
      ENDIF
    ENDIF
    SELECT CASE (SpecifNorm)
    CASE (0)  ! x^{H}Bx=1
      CALL RVScale(n, x, 1/SQRT(t1))
    CASE (1)  ! x^{H}x=1
      t1 = RDotProdItself(n, x)
      CALL RVScale(n, x, 1/SQRT(t1))
    ENDSELECT

  END SUBROUTINE GSEPIIS

  SUBROUTINE GHEPIIS(k, n, M, nM, invD, B, nB, apprlambda, v, w, Tol, &
                     lambda, x, RelAcc, MaxIter, SpecifNorm, NumIter, ErrorCode)
    ! Inverse iteration for one eigenpair of the complex hermitian pencil A*x = lambda*B*x,
    ! extending the factorization A - apprlambda*B = L*D*L^H of the leading (k-1) x (k-1) block
    ! (LDLHF). Converges to the eigenvalue nearest apprlambda; eigenvalues closer than about
    ! epsilon cannot be separated (ErrorCode 2).
    !   M(nM,n)     upper triangle L^H, lower triangle A - apprlambda*B (sizes k-1 in, n out);  invD 1/D_ii
    !   B(nB,n)     overlap matrix
    !   apprlambda  shift: closer to the wanted eigenvalue than to any other and not equal to it
    !   v           start vector (destroyed);   w  work
    !   Tol         relative accuracy; negative: best attainable, at least |Tol| (one extra iteration)
    !   MaxIter     iteration limit;   SpecifNorm  0: x^H*B*x = 1, 1: x^H*x = 1, else max(|Re|,|Im|) = 1
    !   lambda, x   the eigenpair;   RelAcc  accuracy reached;   NumIter  iterations used
    !   ErrorCode   0 ok, 1 singular, 2 Tol not reached in MaxIter

    ! Arguments
    INTEGER     :: k, n, nM, nB, MaxIter, SpecifNorm, NumIter, ErrorCode
    COMPLEX(wp) :: M(nM, n), invD(n), B(nB, n), v(n), w(n), x(n)
    REAL(wp)    :: apprlambda, Tol, lambda, RelAcc
    ! Local variables
    REAL(wp)    :: NormOfDiff, NormOfDiffPrev, t1, t2, sqrtn
    COMPLEX(wp) :: tc
    LOGICAL     :: notconverged

    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0,a,2es11.4,a,es9.2,a,i0)') 'GHEPIIS: k=', k, ' n=', n, ' shift=', apprlambda, &
        ' tol=', Tol, ' maxiter=', MaxIter
    ENDIF
    ErrorCode = 0
    RelAcc = HUGE(RelAcc)
    NormOfDiffPrev = HUGE(NormOfDiffPrev)
    NumIter = 0
    ! Updating M=L*D*LH factorization up to size n using routine LDLHF
    CALL LDLHF(k, n, M, nM, invD, w, ErrorCode)
    IF (ErrorCode > 0) RETURN
    sqrtn = SQRT(2*n*ONE)
    ! Do inverse iterations until the process converges
    ! with relative accuracy Tol, or until the number
    ! of iterations exceeds the limit
    notconverged = .TRUE.
    DO WHILE ((notconverged) .AND. (NumIter < MaxIter))
      IF (NumIter /= 0) v(1:n) = x(1:n)
      ! wp=8: B*v with the hermitian routine MHMVL (half the memory traffic); wp=10/16: MHMV is faster.
      IF (wp == 8) THEN
        CALL MHMVL(n, B, nB, v, w, x)
      ELSE
        CALL MHMV(n, B, nB, v, w, x)
      ENDIF
      CALL LDLHS(n, M, nM, invD, w, x)
      t1 = CMaxAbsReOrIm(n, x)
      CALL CVScale(n, x, CMPLX(1/t1, ZERO, wp))
      ! tc=(x^{H}v)/(v^{H}v)
      tc = CDotProdQuotient(n, x, v)
      ! NormOfDiff=sqrt(2*n)||x-tc*v||
      ! old version: NormOfDiff=sqrtn*CVDiffEucNorm(n,x,tc,v)
      NormOfDiff = CVDiffEucNorm(n, x, tc, v)/SQRT(CDotProdItself(n, x))
      IF (Tol > ZERO) THEN
        IF (NormOfDiff <= Tol) notconverged = .FALSE.
      ELSE
        IF ((NormOfDiff > NormOfDiffPrev) .AND. (NormOfDiff <= ABS(Tol))) notconverged = .FALSE.
        NormOfDiffPrev = NormOfDiff
      ENDIF
      NumIter = NumIter+1
    ENDDO
    RelAcc = NormOfDiff
    IF (notconverged) ErrorCode = 2
    ! Compute x^{H}Bx
    t1 = VMMHMV(n, B, nB, x)
    ! Compute x^{H}Mx
    t2 = VMMHMV(n, M, nM, x)
    ! Rayleigh quotient
    lambda = t2/t1+apprlambda
    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,2es16.8,a,i0,a,i0)') 'GHEPIIS: lambda=', lambda, '  iterations=', NumIter, &
        '  ErrorCode=', ErrorCode
    ENDIF
    SELECT CASE (SpecifNorm)
    CASE (0)  ! x^{H}Bx=1
      CALL CVScale(n, x, CMPLX(1/SQRT(t1), ZERO, wp))
    CASE (1)  ! x^{H}x=1
      t1 = CDotProdItself(n, x)
      CALL CVScale(n, x, CMPLX(1/SQRT(t1), ZERO, wp))
    ENDSELECT

  END SUBROUTINE GHEPIIS

  ! Eigenvalue-index targeting routines ; they use only LDLTF from the rest
  ! of the module.

  FUNCTION InertiaCountBelow(n, invD) RESULT(m)
    ! Number of eigenvalues of the pencil below the shift built into M = A - shift*B, from the
    ! signs of invD (= signs of D_ii) produced by LDLTF. O(n).
    IMPLICIT NONE
    INTEGER  :: n, m, i
    REAL(wp) :: invD(n)

    m = 0
    DO i = 1, n
      IF (invD(i) < ZERO) m = m + 1
    ENDDO

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,i0)') 'InertiaCountBelow: n=', n, ' below=', m
    ENDIF
  END FUNCTION InertiaCountBelow


  SUBROUTINE GSEPIIS_ShiftForIndex(n, A, nA, Sigma0, B, nB, WhichEigval, SigmaLo, SigmaHi, &
                                   Sigma, invD, w, NumFact, ErrorCode)
    ! Bisection on the inertia count until exactly WhichEigval-1 eigenvalues lie below the shift,
    ! so that GSEPIIS then converges to eigenvalue number WhichEigval.
    !   A(nA,n)   lower triangle: A - Sigma0*B in, A - Sigma*B out (restored on failure);
    !             upper triangle: L^T of the final factorization
    !   Sigma0    shift already in A;   B(nB,n) overlap;   SigmaLo, SigmaHi  initial bracket (widened)
    !   Sigma     result;   invD  its 1/D_ii;   w  work;   NumFact  factorizations used
    !   ErrorCode 0 ok, 3 no bracket within Glob_EigIdxMaxBisect steps
    ! Works in place: LDLTF writes only the strict upper triangle and invD, so the lower triangle
    ! is re-shifted between trials.
    IMPLICIT NONE
    INTEGER  :: n, nA, nB, WhichEigval, NumFact, ErrorCode
    REAL(wp) :: A(nA, n), B(nB, n), invD(n), w(n)
    REAL(wp) :: SigmaLo, SigmaHi, Sigma, Sigma0
    ! Local variables:
    INTEGER  :: it, mcount, CountLo, CountHi, FactCode
    REAL(wp) :: lo, hi, width, SigmaNow

    ErrorCode = 0
    NumFact = 0
    lo = SigmaLo
    hi = SigmaHi
    ! SigmaNow always records the shift currently built into A, so that the
    ! next ReShiftAndFactor knows what to subtract and the failure paths know
    ! how to put A back the way it was found.
    SigmaNow = Sigma0
    CountLo = WhichEigval
    CountHi = 0

    ! Establish a valid bracket first. Widen outwards until the number of
    ! eigenvalues below lo is < WhichEigval and the number below hi is >=
    ! WhichEigval. Each widening doubles the half-width.
    width = hi - lo
    IF (width <= ZERO) width = ONE

    DO it = 1, Glob_EigIdxMaxBisect
      CALL ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, lo, invD, w, CountLo, FactCode)
      NumFact = NumFact + 1
      IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8,a,es16.8)') '   ShiftForIndex: factorization ', NumFact, &
          ' shift=', SigmaNow, ' lo=', lo, ' hi=', hi
      ENDIF
      IF (FactCode == 0) THEN
        IF (CountLo < WhichEigval) EXIT
      ENDIF
      lo = lo - width
      width = width + width
    ENDDO
    IF (CountLo >= WhichEigval) THEN
      ErrorCode = 3
      IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,i0,a,i0)') 'GSEPIIS_ShiftForIndex: no bracket for eigenvalue ', WhichEigval, &
          ' after factorizations: ', NumFact
      ENDIF
      CALL ReShiftOnly(n, A, nA, SigmaNow, B, nB, Sigma0)
      RETURN
    ENDIF

    width = hi - lo
    DO it = 1, Glob_EigIdxMaxBisect
      CALL ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, hi, invD, w, CountHi, FactCode)
      NumFact = NumFact + 1
      IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8,a,es16.8)') '   ShiftForIndex: factorization ', NumFact, &
          ' shift=', SigmaNow, ' lo=', lo, ' hi=', hi
      ENDIF
      IF (FactCode == 0) THEN
        IF (CountHi >= WhichEigval) EXIT
      ENDIF
      hi = hi + width
      width = width + width
    ENDDO
    IF (CountHi < WhichEigval) THEN
      ErrorCode = 3
      IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,i0,a,i0)') 'GSEPIIS_ShiftForIndex: no bracket for eigenvalue ', WhichEigval, &
          ' after factorizations: ', NumFact
      ENDIF
      CALL ReShiftOnly(n, A, nA, SigmaNow, B, nB, Sigma0)
      RETURN
    ENDIF

    ! Bisect until the interval is negligible. Any point in the final
    ! interval has exactly WhichEigval-1 eigenvalues below it.
    DO it = 1, Glob_EigIdxMaxBisect
      IF (ABS(hi-lo) <= Glob_EigIdxBisectTol*(ABS(lo)+ABS(hi)+ONE)) EXIT
      Sigma = (lo+hi)/TWO
      CALL ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, Sigma, invD, w, mcount, FactCode)
      NumFact = NumFact + 1
      IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
        WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8,a,es16.8)') '   ShiftForIndex: factorization ', NumFact, &
          ' shift=', SigmaNow, ' lo=', lo, ' hi=', hi
      ENDIF
      IF (FactCode /= 0) THEN
        ! A singular factorization means the trial shift hit an eigenvalue: nudge it and retry once,
        ! then fall through to the bracket update (a CYCLE here would repeat the same midpoint).
        Sigma = Sigma + Glob_EigIdxBisectTol*(ABS(Sigma)+ONE)
        CALL ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, Sigma, invD, w, mcount, FactCode)
        NumFact = NumFact + 1
        IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
          WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8,a,es16.8)') '   ShiftForIndex: factorization ', NumFact, &
            ' shift=', SigmaNow, ' lo=', lo, ' hi=', hi
        ENDIF
        IF (FactCode /= 0) THEN
          ! Still singular: give up on this midpoint, pull the top of the bracket
          ! down to it and carry on. The invariant CountLo < k <= CountHi is
          ! preserved because a singular shift cannot lie outside the bracket.
          hi = Sigma
          CYCLE
        ENDIF
      ENDIF
      IF (mcount < WhichEigval) THEN
        lo = Sigma
      ELSE
        hi = Sigma
      ENDIF
    ENDDO

    ! Return the low end of the final bracket: it has exactly WhichEigval-1
    ! eigenvalues below it, so lambda_WhichEigval is the nearest eigenvalue
    ! above, which is what inverse iteration will lock onto. Leave A holding
    ! the matching shift and factorization.
    Sigma = lo
    CALL ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, Sigma, invD, w, mcount, FactCode)
    NumFact = NumFact + 1
    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8,a,es16.8)') '   ShiftForIndex: factorization ', NumFact, ' shift=', &
        SigmaNow, ' lo=', lo, ' hi=', hi
    ENDIF

    IF ((Verbose >= 3) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,es16.8,a,i0,a,i0)') 'GSEPIIS_ShiftForIndex: eigenvalue ', WhichEigval, ' shift=', &
        Sigma, ' below=', mcount, ' factorizations=', NumFact
    ENDIF
  END SUBROUTINE GSEPIIS_ShiftForIndex


  SUBROUTINE ReShiftAndFactor(n, A, nA, SigmaNow, B, nB, Sigma, invD, w, mcount, ErrorCode)
    ! Re-shifts the lower triangle of A from SigmaNow to Sigma,
    !   A - Sigma*B = (A - SigmaNow*B) + (SigmaNow - Sigma)*B,
    ! factorizes it and returns the inertia count; SigmaNow is updated. Used by GSEPIIS_ShiftForIndex.
    IMPLICIT NONE
    INTEGER  :: n, nA, nB, mcount, ErrorCode
    REAL(wp) :: A(nA, n), B(nB, n), invD(n), w(n), Sigma, SigmaNow

    CALL ReShiftOnly(n, A, nA, SigmaNow, B, nB, Sigma)

    CALL LDLTF(1, n, A, nA, invD, w, ErrorCode)
    IF (ErrorCode > 0) THEN
      mcount = 0
      RETURN
    ENDIF
    mcount = InertiaCountBelow(n, invD)

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,es16.8,a,i0,a,i0)') 'ReShiftAndFactor: n=', n, ' shift=', SigmaNow, ' below=', &
        mcount, ' ErrorCode=', ErrorCode
    ENDIF
  END SUBROUTINE ReShiftAndFactor


  SUBROUTINE ReShiftOnly(n, A, nA, SigmaNow, B, nB, Sigma)
    ! Replaces the shift built into the lower triangle of A, from SigmaNow to
    ! Sigma, and updates SigmaNow to match. No factorization.
    IMPLICIT NONE
    INTEGER  :: n, nA, nB
    REAL(wp) :: A(nA, n), B(nB, n), Sigma, SigmaNow
    INTEGER  :: i, j
    REAL(wp) :: d

    IF ((Verbose >= 4) .AND. (Glob_ProcID == 0)) THEN
      WRITE(*, '(1x,a,i0,a,es16.8,a,es16.8)') 'ReShiftOnly: n=', n, ' from=', SigmaNow, ' to=', Sigma
    ENDIF
    d = SigmaNow - Sigma
    IF (d /= ZERO) THEN
      DO j = 1, n
        DO i = j, n
          A(i, j) = A(i, j) + d*B(i, j)
        ENDDO
      ENDDO
    ENDIF
    SigmaNow = Sigma

  END SUBROUTINE ReShiftOnly

END MODULE linalg
