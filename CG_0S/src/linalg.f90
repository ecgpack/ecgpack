module linalg !version 0.52
!Module linalg contains some linear algebra routines
!that are used in calculations
  use globvars
  implicit none

!These global variables define in what mode (parallel or using
!a single process only) certain routines should work. They
!also define whether BLAS or generic routines should be
!used (it concerns only the single process mode). The
!values of these variables may change depending on the size
!of the eigenvalue problem, architecture, etc.
!Basically, "0" stands for "no" and "1" stands for "yes".

  integer :: Glob_LDLTF_PMode=1
  integer :: Glob_LDLTF_SModeType=0

  integer :: Glob_LDLHF_PMode=1
  integer :: Glob_LDLHF_SModeType=0

  integer :: Glob_LDLTS_PMode=1
  integer :: Glob_LDLTS_SModeType=0
  integer :: Glob_LDLTS_UseBLAS=0

  integer :: Glob_LDLHS_PMode=1
  integer :: Glob_LDLHS_SModeType=0
  integer :: Glob_LDLHS_UseBLAS=0

  integer :: Glob_MTMVL_PMode=1
  integer :: Glob_MTMVL_UseBLAS=0

  integer :: Glob_MHMVL_PMode=1
  integer :: Glob_MHMVL_UseBLAS=0

  integer :: Glob_MTMV_PMode=1
  integer :: Glob_MTMV_UseBLAS=0

  integer :: Glob_MHMV_PMode=1
  integer :: Glob_MHMV_UseBLAS=0

  integer :: Glob_VMMTMV_PMode=1

  integer :: Glob_VMMHMV_PMode=1

  integer :: Glob_RMaxAbsEl_PMode=0

  integer :: Glob_CMaxAbsReOrIm_PMode=0

  integer :: Glob_RDotProd_PMode=0
  integer :: Glob_RDotProd_UseBLAS=0

  integer :: Glob_CDotProd_PMode=0
  integer :: Glob_CDotProd_UseBLAS=0

  integer :: Glob_RDotProdItself_PMode=0
  integer :: Glob_RDotProdItself_UseBLAS=0

  integer :: Glob_CDotProdItself_PMode=0
  integer :: Glob_CDotProdItself_UseBLAS=0

  integer :: Glob_RDotProdQuotient_PMode=0

  integer :: Glob_CDotProdQuotient_PMode=0

  integer :: Glob_RVScale_UseBLAS=0

  integer :: Glob_CVScale_UseBLAS=0

  integer :: Glob_RVDiffEucNorm_PMode=0

  integer :: Glob_CVDiffEucNorm_PMode=0

!Column-panel width used by the blocked single-process branches of the
!LDLT/LDLH factorization routines (LDLTF and LDLHF). Processing a panel of
!Glob_LDLF_BlockSize columns per sweep lets every already-completed column
!of the factor be streamed from memory once per panel instead of once per
!new column, cutting the memory traffic of the O(n^3) factorization by
!roughly this factor (the algorithm is memory-bandwidth-bound for matrices
!that exceed the cache, i.e. for n of a few hundred and beyond). The value
!is a trade-off: it must be large enough to amortize the streaming but small
!enough that the panel buffer stays cache-resident. Measured on x86-64
!(gfortran -O3 -march=native, n=2000): 8 and 16 are clearly suboptimal,
!the optimum is around 48, and the curve is flat between 32 and 64.
!The routines also use this value as the routing threshold between the
!blocked serial algorithm and the parallel branch (see LDLTF/LDLHF).
  integer,parameter :: Glob_LDLF_BlockSize=48

!Performance-model coefficients used to route bulk LDLT/LDLH factorizations
!between the blocked serial path and the parallel branch when several MPI
!processes are available (wp=8 only; suffix R = real/LDLTF, C =
!complex/LDLHF). They are measured once per run by linalg_setparam, which
!must then be called COLLECTIVELY on all processes: the serial blocked path
!is modeled as t_ser(n) = Cs*n^3 and the parallel branch as
!t_par(n) = Ca*n^3 + Cb*n^2 (compute + communication-latency terms), both
!fitted to timings of small synthetic factorizations executed at the actual
!process count on the actual machine. A negative Cs means "not calibrated":
!bulk factorizations then stay in the parallel branch whenever more than one
!process runs (the behavior of the original code), so skipping the
!calibration can cost performance but can never lose to the original code.
  real(wp) :: Glob_LDLF_CsR=-1.0_wp,Glob_LDLF_CaR=0.0_wp,Glob_LDLF_CbR=0.0_wp
  real(wp) :: Glob_LDLF_CsC=-1.0_wp,Glob_LDLF_CaC=0.0_wp,Glob_LDLF_CbC=0.0_wp

contains

  subroutine linalg_setparam(n)
!Subroutine linalg_setparam sets the values of the global parameters
!declared in module linalg that provide maximum performance (roughly).
!The user needs to pass the dimension of the linear algebra problems, n.
!The optimal values of the global parameters are specific for each
!particular n. Therefore it is important that the user calls linalg_setparam
!each time the dimensionality of the linear algebra problem is changed
!in his/her calculations.
!The optimal values of the global parameters also quite depend on the
!computer architecture and the processor interconnect speed. Below one can
!find an implementation that is based on rough tests made on a SGI Altix 4700
!supercomputer and on an a AMD Athlon MP PC cluster with Gigabit Ethernet.
!Each time the present code is used on a different machine, it is highly
!recommended to run performance tests and determine the sets of optimal
!parameters as a function of n.
!
!Present implementation. The heavy O(n^2)/O(n^3) routines (the LDL
!factorization/solve and the matrix-vector / quadratic-form routines) carry
!the whole cost of solving or updating an eigenvalue problem of size
!n = 1000..10000 on 1..128 cores. Benchmarks on such problems show that the
!parallel branches (PMode=1) become faster than running the routines serially
!(redundantly on every process) once there is more than one MPI process: the
!O(n^2/P) work distribution then outweighs the MPI_ALLREDUCE/MPI_BCAST traffic,
!and the advantage grows with n and with P. On a single process the parallel
!branches only add communication-call overhead (the parallel factorization
!alone issues ~n^2/2 tiny reductions), so running serially is at least as fast.
!Hence the crossover sits between 1 and 2 processes for every supported
!precision, and the rule is: serial when Glob_NumOfProcs==1, parallel otherwise.
!
!Exception - the LDL factorizations in double precision (wp=8). Their
!single-process branches are cache blocked (see Glob_LDLF_BlockSize), which
!makes the bulk factorization several times faster than the fine-grained
!parallel branch up to a machine-, compiler- and size-dependent process
!count (e.g. up to ~48 processes at n=3000 with gfortran on an EPYC node,
!but only up to ~4 with ifort, whose code generation for the blocked kernel
!is several times slower). LDLTF/LDLHF therefore route internally between
!the blocked serial path (executed by rank 0 alone and broadcast,
!SModeType=0) and the parallel branch, using small performance models
!calibrated at first call here (see Glob_LDLF_Cs*/Ca*/Cb* and
!linalg_ldlf_calibrate below). Because the calibration involves collective
!MPI operations, linalg_setparam MUST be called collectively by all
!processes when Glob_NumOfProcs>1. Thin incremental updates (fewer than one
!panel of new columns, e.g. the m=n single-column update) always keep the
!parallel branch at np>1. For wp=10/16 the blocked kernel does not pay
!(scalar x87 / software-emulated arithmetic), so those builds keep the
!plain algorithm and the original serial-vs-parallel rule for the
!factorizations as well.
!
!How strongly this matters, however, depends on the working-precision kind wp,
!because the cost of the arithmetic between communication calls scales with it
!while the communication cost does not. This routine therefore selects the
!behaviour per wp:
!  * wp=8  (double, hardware fp64): serial at one process is a clear win
!          (~12% on the full factorization) - fp64 arithmetic is fast, so the
!          fine-grained reductions of the parallel branches weigh heavily.
!  * wp=10 (extended, 80-bit): arithmetic is slower, so the single-process
!          serial advantage shrinks to a few percent; crossover still at 2.
!          MPI_REAL16 reductions work under both OpenMPI and Intel MPI here.
!  * wp=16 (quadruple, software-emulated, ~1-2 orders of magnitude slower):
!          communication is negligible next to the compute, so serial and
!          parallel are within noise at one process while parallel scales
!          almost ideally beyond it (~1.9x on 2 cores). Serial at one process
!          is kept because it is no slower AND it avoids the MPI_REAL16
!          MPI_ALLREDUCE calls, which are BROKEN in OpenMPI - quad builds must
!          use an Intel-MPI toolchain (e.g. intel-2025b/intel-2023b) to run the
!          parallel branches at >1 process.
!The O(n) vector helpers (dot products, scalings, norms, max-element) are never
!worth parallelizing at these sizes for any precision, so they stay serial.
!UseBLAS is left at 0 because this standalone harness links no BLAS library (the
!BLAS paths are satisfied by the stub routines in main.f90); set it only in
!builds that actually link BLAS.
    integer n
    integer :: pm, min_procs_par

    !Smallest process count at which the parallel branches pay off, selected per
    !working precision. Benchmarks place it at 2 for all three supported kinds;
    !it is kept as an explicit per-wp knob so it can be retuned independently.
    select case (wp)
    case (8)             !double precision
      min_procs_par=2
    case (10)            !extended (80-bit) precision
      min_procs_par=2
    case (16)            !quadruple precision (needs Intel MPI for the parallel branches)
      min_procs_par=2
    case default
      min_procs_par=2
    end select

    if (Glob_NumOfProcs>=min_procs_par) then
      pm=1   !parallel
    else
      pm=0   !single process - avoid all MPI call overhead
    endif

    !Heavy routines: parallel on >1 process, serial on a single process.
    !For the factorizations, PMode=1 only means "the parallel branch may be
    !used": LDLTF/LDLHF route internally between the cache-blocked serial
    !algorithm (bulk factorizations, where it beats the fine-grained
    !parallel branch on at least 1-8 processes) and the parallel branch
    !(thin incremental updates of fewer than Glob_LDLF_BlockSize columns).
    Glob_LDLTF_PMode=pm
    Glob_LDLHF_PMode=pm
    Glob_LDLTS_PMode=pm
    Glob_LDLHS_PMode=pm
    Glob_MTMVL_PMode=pm
    Glob_MHMVL_PMode=pm
    Glob_MTMV_PMode=pm
    Glob_MHMV_PMode=pm
    Glob_VMMTMV_PMode=pm
    Glob_VMMHMV_PMode=pm

    !When the blocked serial factorization runs while several MPI processes
    !are active (the bulk-factorization routing above), rank 0 computes it
    !alone and broadcasts the factor (SModeType=0): measured faster than
    !all processes computing it redundantly, because a single process keeps
    !the whole memory bandwidth to itself and the broadcast of the factor
    !is comparatively cheap on a shared-memory node. On a single process
    !the broadcast step is pointless, so SModeType=1 there.
    Glob_LDLTF_SModeType=1-pm
    Glob_LDLHF_SModeType=1-pm

    !O(n) vector helpers: always serial at these problem sizes.
    Glob_RMaxAbsEl_PMode=0
    Glob_CMaxAbsReOrIm_PMode=0
    Glob_RDotProd_PMode=0
    Glob_CDotProd_PMode=0
    Glob_RDotProdItself_PMode=0
    Glob_CDotProdItself_PMode=0
    Glob_RDotProdQuotient_PMode=0
    Glob_CDotProdQuotient_PMode=0
    Glob_RVDiffEucNorm_PMode=0
    Glob_CVDiffEucNorm_PMode=0

    !Calibrate the factorization-routing performance models once per run
    !(wp=8 with several processes only; see the comments at
    !Glob_LDLF_Cs*/Ca*/Cb* and in linalg_ldlf_calibrate). This involves
    !collective MPI operations - hence the collective-call requirement
    !stated above.
    if ((wp==8).and.(Glob_NumOfProcs>1).and.(Glob_LDLF_CsR<ZERO)) then
      call linalg_ldlf_calibrate
    endif

  end subroutine linalg_setparam

  subroutine linalg_ldlf_calibrate
!Subroutine linalg_ldlf_calibrate measures, on the actual machine and at the
!actual MPI process count, how fast the blocked serial path and the parallel
!branch of the LDLT/LDLH factorizations are, and stores the coefficients of
!two small performance models:
!    t_serial(n)   = Cs*n^3
!    t_parallel(n) = Ca*n^3 + Cb*n^2
!(the n^3 term is arithmetic, the n^2 term is the per-column communication
!latency of the parallel branch). LDLTF/LDLHF compare the two models at the
!actual problem size to decide where to run each bulk factorization, which
!makes the routing robust across machines (allreduce latency, memory
!bandwidth) and compilers (the blocked kernel is 5-7x faster than the plain
!one under gfortran/ifx, but ifort generates much slower code for both).
!The measurements factorize small synthetic well-conditioned matrices with
!the real code paths of LDLTF/LDLHF: the parallel branch collectively at
!sizes m1 and m2 (fitting Ca and Cb), the blocked serial path on rank 0
!alone at size m2 (fitting Cs). The mode flags are saved and restored, and
!the resulting coefficients are broadcast from rank 0, so all processes are
!guaranteed to take identical routing decisions afterwards.
!This routine is COLLECTIVE over MPI_COMM_WORLD; it costs a few tens of
!milliseconds (a few hundred with slow compilers) once per run.
    integer,parameter :: m1=320,m2=640,nrep=2
    real(wp),allocatable    :: AR0(:,:),AR(:,:),invDR(:),wR(:)
    complex(wp),allocatable :: AC0(:,:),AC(:,:),invDC(:),wC(:)
    real(wp)         :: dr,tp1r,tp2r,tp1c,tp2c,tsr,tsc,cf(6)
    double precision :: tt0,tt1,tbest
    integer          :: i,j,r,ec
    integer          :: savePMr,saveSMr,savePMc,saveSMc

    savePMr=Glob_LDLTF_PMode
    saveSMr=Glob_LDLTF_SModeType
    savePMc=Glob_LDLHF_PMode
    saveSMc=Glob_LDLHF_SModeType

    allocate(AR0(m2,m2),AR(m2,m2),invDR(m2),wR(m2))
    allocate(AC0(m2,m2),AC(m2,m2),invDC(m2),wC(m2))
    !Synthetic diagonally dominant test matrices (only the lower triangles
    !are referenced by the factorizations)
    do j=1,m2
      do i=j,m2
        dr=ONE/(ONE+abs(i-j))
        AR0(i,j)=dr
        if (i==j) then
          AC0(i,j)=cmplx(dr+3*ONE,ZERO,wp)
        else
          AC0(i,j)=cmplx(dr,(3*dr)/(10*(1+mod(i+j,7))),wp)
        endif
      enddo
      AR0(j,j)=AR0(j,j)+3*ONE
    enddo

    !--- Parallel branch, collectively, at sizes m1 and m2. The routing in
    !LDLTF/LDLHF sends bulk factorizations to the parallel branch as long
    !as the models are not calibrated yet (Cs<0), which is exactly the
    !state here, so a plain call with PMode=1 measures the parallel branch.
    Glob_LDLTF_PMode=1
    Glob_LDLHF_PMode=1
    tbest=huge(tbest)
    do r=1,nrep
      AR(1:m1,1:m1)=AR0(1:m1,1:m1)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt0=MPI_WTIME()
      call LDLTF(1,m1,AR,m2,invDR,wR,ec)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt1=MPI_WTIME()
      if (tt1-tt0<tbest) tbest=tt1-tt0
    enddo
    tp1r=real(tbest,wp)
    tbest=huge(tbest)
    do r=1,nrep
      AR=AR0
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt0=MPI_WTIME()
      call LDLTF(1,m2,AR,m2,invDR,wR,ec)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt1=MPI_WTIME()
      if (tt1-tt0<tbest) tbest=tt1-tt0
    enddo
    tp2r=real(tbest,wp)
    tbest=huge(tbest)
    do r=1,nrep
      AC(1:m1,1:m1)=AC0(1:m1,1:m1)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt0=MPI_WTIME()
      call LDLHF(1,m1,AC,m2,invDC,wC,ec)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt1=MPI_WTIME()
      if (tt1-tt0<tbest) tbest=tt1-tt0
    enddo
    tp1c=real(tbest,wp)
    tbest=huge(tbest)
    do r=1,nrep
      AC=AC0
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt0=MPI_WTIME()
      call LDLHF(1,m2,AC,m2,invDC,wC,ec)
      call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)
      tt1=MPI_WTIME()
      if (tt1-tt0<tbest) tbest=tt1-tt0
    enddo
    tp2c=real(tbest,wp)

    !--- Blocked serial path at size m2, on rank 0 alone (the other ranks
    !wait at the barrier), matching how the serial path actually runs at
    !np>1: PMode=0 forces the serial path, SModeType=1 removes its MPI
    !broadcasts so the call is local to each process.
    Glob_LDLTF_PMode=0
    Glob_LDLHF_PMode=0
    Glob_LDLTF_SModeType=1
    Glob_LDLHF_SModeType=1
    tsr=ZERO
    tsc=ZERO
    if (Glob_ProcID==0) then
      tbest=huge(tbest)
      do r=1,nrep
        AR=AR0
        tt0=MPI_WTIME()
        call LDLTF(1,m2,AR,m2,invDR,wR,ec)
        tt1=MPI_WTIME()
        if (tt1-tt0<tbest) tbest=tt1-tt0
      enddo
      tsr=real(tbest,wp)
      tbest=huge(tbest)
      do r=1,nrep
        AC=AC0
        tt0=MPI_WTIME()
        call LDLHF(1,m2,AC,m2,invDC,wC,ec)
        tt1=MPI_WTIME()
        if (tt1-tt0<tbest) tbest=tt1-tt0
      enddo
      tsc=real(tbest,wp)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIErrCode)

    Glob_LDLTF_PMode=savePMr
    Glob_LDLTF_SModeType=saveSMr
    Glob_LDLHF_PMode=savePMc
    Glob_LDLHF_SModeType=saveSMc

    !--- Fit the models on rank 0 and broadcast the coefficients so that
    !every process takes identical routing decisions (different decisions
    !on different ranks would deadlock the collective branches).
    if (Glob_ProcID==0) then
      call fitpar(tp1r,tp2r,cf(1),cf(2))
      call fitpar(tp1c,tp2c,cf(4),cf(5))
      cf(3)=tsr/(real(m2,wp)**3)
      cf(6)=tsc/(real(m2,wp)**3)
    endif
    call MPI_BCAST(cf,6,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
    Glob_LDLF_CaR=cf(1)
    Glob_LDLF_CbR=cf(2)
    Glob_LDLF_CsR=cf(3)
    Glob_LDLF_CaC=cf(4)
    Glob_LDLF_CbC=cf(5)
    Glob_LDLF_CsC=cf(6)

  contains

    subroutine fitpar(t1,t2,a,b)
!Fits t = a*m^3 + b*m^2 through (m1,t1) and (m2,t2), clamping both
!coefficients to be nonnegative (protects the extrapolation against
!measurement noise; when one coefficient is clamped the other is
!recomputed from the larger, more reliable measurement).
      real(wp) t1,t2,a,b
      real(wp) f1,f2
      f1=real(m1,wp)
      f2=real(m2,wp)
      a=(t1*f2**2-t2*f1**2)/(f1**3*f2**2-f2**3*f1**2)
      if (a<ZERO) then
        a=ZERO
        b=t2/f2**2
      else
        b=(t1-a*f1**3)/f1**2
        if (b<ZERO) then
          b=ZERO
          a=t2/f2**3
        endif
      endif
    end subroutine fitpar

  end subroutine linalg_ldlf_calibrate

  subroutine LDLTF(m,n,A,nA,invD,w,ErrorCode)
!Subroutine LDLTF updates the factorization of a real symmetric
!matrix A of size n in the form A=L*D*LT using a known factorization
!of its submatrix of size m-1. Here L is a lower triangular matrix
!with diagonal elements equal to 1, D is a diagonal matrix, LT is
!the transpose of L. If no factorization of a submatrix is
!known, one should set m=1. The lower triangle of A (including the
!diagonal) remains unchanged on exit.
!  Input parameters :
!     m-1 - The size of the submatrix of A, whose factorization is
!           known (m>=1);
!       n - The size of matrix A;
!      nA - The leading dimension of A;
!       A - A two-dimensional array. Its lower triangle (including
!           the diagonal) contains matrix A and the upper one
!           contains the elements of matrix LT (of size m-1);
!    invD - An array containing inverse values of
!           the diagonal elements of D (up to size m-1);
!       w - A work vector (at least of size n)
!  Output parameters:
!       A - A two-dimensional array. The elements of matrix A
!           (including the diagonal) are stored in the lower
!           triangle. Matrix LT is returned in the upper
!           triangle (updated to size n);
!    invD - An array containing inverse values of
!           the diagonal elements of D (updated to size n);
!ErrorCode- The error flag. If ErrorCode=0 then the procedure
!           finished successfully. If ErrorCode=1 then the
!           procedure cannot be completed because matrix A is
!           singular

!Arguments :
    integer               m,n,nA,ErrorCode
    real(wp)    A(nA,n),invD(n),w(n)
!Local variables :
    integer        i,j,jm,im,i0,ib,ibs
    integer        q,k,ji,jf,jim,jiR,RowsPerProc,mod_im_Glob_NumOfProcs
    real(wp)    x,y,z
    real(wp)    tv(Glob_LDLF_BlockSize)
    real(wp),allocatable :: Pnl(:,:)
    logical        serialbulk

    ErrorCode=0
    !Routing between the serial path and the parallel branch:
    ! - PMode==0 always selects the serial path (as in the original code);
    ! - for hardware-SIMD double precision (wp=8), a bulk factorization (at
    !   least one full panel of new columns) additionally runs on the
    !   blocked serial path when only one process is available, or when the
    !   performance models calibrated by linalg_setparam predict the
    !   blocked serial path (rank 0 computes alone and broadcasts,
    !   SModeType=0) to be faster than the parallel branch at this n and
    !   process count. Without calibration, bulk factorizations stay in
    !   the parallel branch at np>1 - i.e. the original behavior.
    !Thin incremental updates (fewer than one panel of new columns, e.g.
    !the m=n single-column update) keep the parallel branch at np>1, where
    !their O((n-m+1)*n/P) work split still pays. For extended (wp=10,
    !scalar x87) and quadruple (wp=16, software-emulated) precision the
    !blocked kernel does not pay - the plain column algorithm is used
    !serially and the parallel branch keeps its full role at more than one
    !process (ibs=1 below selects the plain path).
    ibs=Glob_LDLF_BlockSize
    if (wp/=8) ibs=1
    serialbulk=(Glob_LDLTF_PMode==0)
    if ((.not.serialbulk).and.(wp==8).and.(n-m+1>=Glob_LDLF_BlockSize)) then
      if (Glob_NumOfProcs==1) then
        serialbulk=.true.
      elseif (Glob_LDLF_CsR>=ZERO) then
        !t_ser(n)<=t_par(n) with both sides divided by n^2
        if (Glob_LDLF_CsR*n<=Glob_LDLF_CaR*n+Glob_LDLF_CbR) serialbulk=.true.
      endif
    endif
    if (serialbulk) then
      !Single process version, blocked for cache efficiency. The columns
      !i=m..n are processed in panels of Glob_LDLF_BlockSize columns held
      !in the transposed buffer Pnl (Pnl(k,q) is element q of the k-th
      !panel column). All forward-substitution work of the panel columns
      !against the already completed part of the factor (columns 1..i0-1)
      !is done inside one sweep over that part, so each completed column is
      !streamed from memory once per panel instead of once per new column,
      !and the innermost loop of the sweep runs over the independent panel
      !entries Pnl(1:ib,..), which are contiguous and therefore vectorize
      !well. The arithmetic performed for every column (and the order of
      !the operations within each column) is identical to the unblocked
      !algorithm, so the results are bitwise identical; only the memory
      !access pattern changes. Single-column panels (in particular the
      !incremental update case m=n) take the original unblocked path.
      if ((Glob_LDLTF_SModeType/=0).or.(Glob_ProcID==0)) then
        if ((ibs>1).and.(n>m)) allocate(Pnl(ibs,n))
        do i0=m,n,ibs
          ib=min(ibs,n-i0+1)
          if (ib==1) then
            !Plain column algorithm (original unblocked code)
            i=i0
            do j=1,i-1
              jm=j-1
              x=A(i,j)
              do q=1,jm
                x=x-A(q,i)*A(q,j)
              enddo
              A(j,i)=x
            enddo
            im=i-1
            do q=1,im
              w(q)=A(q,i)*invD(q)
            enddo
            x=A(i,i)
            do q=1,im
              x=x-A(q,i)*w(q)
            enddo
            do q=1,im
              A(q,i)=w(q)
            enddo
            if (x==ZERO) then
              ErrorCode=1
              return
            endif
            invD(i)=ONE/x
          else
            !Seed the panel buffer with the matrix elements of rows
            !i0..i0+ib-1 taken from the (untouched) lower triangle
            do k=1,ib
              i=i0+k-1
              do j=1,i-1
                Pnl(k,j)=A(i,j)
              enddo
            enddo
            !Shared sweep: forward-substitute all panel columns against
            !the completed columns 1..i0-1 (this is where blocking pays)
            do j=1,i0-1
              jm=j-1
              do k=1,ib
                tv(k)=Pnl(k,j)
              enddo
              do q=1,jm
                x=A(q,j)
                do k=1,ib
                  tv(k)=tv(k)-Pnl(k,q)*x
                enddo
              enddo
              do k=1,ib
                Pnl(k,j)=tv(k)
              enddo
            enddo
            !Finish the panel columns one by one: substitution against the
            !panel columns completed just before them and the diagonal
            !step, then write the scaled column of the factor back into
            !the upper triangle of A
            do k=1,ib
              i=i0+k-1
              do j=i0,i-1
                jm=j-1
                x=Pnl(k,j)
                do q=1,jm
                  x=x-Pnl(k,q)*A(q,j)
                enddo
                Pnl(k,j)=x
              enddo
              im=i-1
              do q=1,im
                w(q)=Pnl(k,q)*invD(q)
              enddo
              x=A(i,i)
              do q=1,im
                x=x-Pnl(k,q)*w(q)
              enddo
              do q=1,im
                A(q,i)=w(q)
              enddo
              if (x==ZERO) then
                ErrorCode=1
                return
              endif
              invD(i)=ONE/x
            enddo
          endif
        enddo
      endif
      if (Glob_LDLTF_SModeType==0) then
        do i=m,n
          call MPI_BCAST(A(1:i-1,i),i-1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        enddo
        call MPI_BCAST(invD(m:n),n-m+1,MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      endif
    else
      !Parallel version (far from being perfect...)
      do i=m,n
        im=i-1
        A(1:im,i)=ZERO
        RowsPerProc=im/Glob_NumOfProcs
        mod_im_Glob_NumOfProcs=mod(im,Glob_NumOfProcs)
        do k=1,RowsPerProc
          jf=k*Glob_NumOfProcs
          jim=jf-Glob_NumOfProcs
          ji=jim+1
          jiR=ji+Glob_ProcID
          !A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
          z=ZERO
          do q=1,jim
            z=z-A(q,i)*A(q,jiR)
          enddo
          A(jiR,i)=z
          !
          call MPI_ALLREDUCE(A(ji:jf,i),w(1:Glob_NumOfProcs),Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          A(ji:jf,i)=w(1:Glob_NumOfProcs)
          do j=ji,jf
            jm=j-1
            !A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
            z=A(j,i)+A(i,j)
            do q=ji,jm
              z=z-A(q,i)*A(q,j)
            enddo
            A(j,i)=z
            !
          enddo
        enddo
        if (mod_im_Glob_NumOfProcs>0) then
          ji=RowsPerProc*Glob_NumOfProcs+1
          jim=ji-1
          jf=im
          jiR=ji+Glob_ProcID
          if (jiR<i) then
            !A(jiR,i)=-dot_product(A(1:jim,i),A(1:jim,jiR))
            z=ZERO
            do q=1,jim
              z=z-A(q,i)*A(q,jiR)
            enddo
            A(jiR,i)=z
            !
          endif
          call MPI_ALLREDUCE(A(ji:jf,i),w(1:mod_im_Glob_NumOfProcs),mod_im_Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          A(ji:jf,i)=w(1:mod_im_Glob_NumOfProcs)
          do j=ji,jf
            jm=j-1
            !A(j,i)=A(j,i)+A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j))
            z=A(j,i)+A(i,j)
            do q=ji,jm
              z=z-A(q,i)*A(q,j)
            enddo
            A(j,i)=z
            !
          enddo
        endif
        !j==1 case
        w(1:im)=A(1:im,i)*invD(1:im)
        y=ZERO
        do k=1+Glob_ProcID,im,Glob_NumOfProcs
          y=y+A(k,i)*w(k)
        enddo
        call MPI_ALLREDUCE(y,x,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x=A(i,i)-x
        A(1:im,i)=w(1:im)
        if (x==ZERO) then
          ErrorCode=1
          return
        endif
        invD(i)=ONE/x
      enddo
    endif

  end subroutine LDLTF

  subroutine LDLHF(m,n,A,nA,invD,w,ErrorCode)
!Subroutine LDLHF updates the factorization of a complex hermitian
!matrix A of size n in the form A=L*D*LH using a known factorization
!of its submatrix of size m-1. Here L is a lower triangular matrix
!with diagonal elements equal to 1, D is a diagonal matrix, LH is
!the hermitian conjugate of L. If no factorization of a submatrix is
!known, one should set m=1. The lower triangle of A (including the
!diagonal) remains unchanged on exit.
!  Input parameters :
!     m-1 - The size of the submatrix of A, whose factorization is
!           known (m>=1);
!       n - The size of matrix A;
!      nA - The leading dimension of A;
!       A - A two-dimensional array. Its lower triangle (including
!           the diagonal) contains matrix A and the upper one
!           contains the elements of matrix LH (of size m-1);
!    invD - An array containing inverse values of
!           the diagonal elements of D (up to size m-1);
!       w - A work vector (at least of size n)
!  Output parameters:
!       A - A two-dimensional array. The elements of matrix A
!           (including the diagonal) are stored in the lower
!           triangle. Matrix LH is returned in the upper
!           triangle (updated to size n);
!    invD - An array containing inverse values of
!           the diagonal elements of D (updated to size n);
!ErrorCode- The error flag. If ErrorCode=0 then the procedure
!           finished successfully. If ErrorCode=1 then the
!           procedure cannot be completed because matrix A is
!           singular

!Arguments :
    integer               m,n,nA,ErrorCode
    complex(wp) A(nA,n),invD(n),w(n)
!Local variables :
    integer        i,j,jm,im,i0,ib,ibs
    integer        q,k,ji,jf,jim,jiR,RowsPerProc,mod_im_Glob_NumOfProcs
    complex(wp) x,y
    real(wp)    xr,xi
    real(wp)    tvr(Glob_LDLF_BlockSize),tvi(Glob_LDLF_BlockSize)
    real(wp),allocatable :: PnR(:,:),PnI(:,:)
    logical        serialbulk

    ErrorCode=0
    !Routing between the blocked serial path and the parallel branch: same
    !rule and rationale as in LDLTF (see the comments there), including the
    !restriction of the blocked kernel to wp=8, with the complex-case
    !performance-model coefficients.
    ibs=Glob_LDLF_BlockSize
    if (wp/=8) ibs=1
    serialbulk=(Glob_LDLHF_PMode==0)
    if ((.not.serialbulk).and.(wp==8).and.(n-m+1>=Glob_LDLF_BlockSize)) then
      if (Glob_NumOfProcs==1) then
        serialbulk=.true.
      elseif (Glob_LDLF_CsC>=ZERO) then
        !t_ser(n)<=t_par(n) with both sides divided by n^2
        if (Glob_LDLF_CsC*n<=Glob_LDLF_CaC*n+Glob_LDLF_CbC) serialbulk=.true.
      endif
    endif
    if (serialbulk) then
      !Single process version, blocked for cache efficiency exactly like
      !the single process branch of LDLTF (see the comments there). The
      !panel buffers PnR/PnI hold the real and imaginary planes of the
      !unconjugated substitution intermediates (the original algorithm
      !stores their conjugates in the upper triangle of A) as two separate
      !real arrays: the sweep kernels then consist of pure real
      !multiply-add loops, which every supported compiler vectorizes well,
      !whereas complex-typed inner loops are pessimized by some (measured:
      !the plane-split kernel is 1.1-1.8x faster than the complex-typed one
      !depending on the compiler). The arithmetic follows the same formulas
      !as the unblocked algorithm; only the compiler's fused-multiply-add
      !contraction pattern of the complex products may differ, so results
      !can deviate in the last bits (as they already do between compilers).
      if ((Glob_LDLHF_SModeType/=0).or.(Glob_ProcID==0)) then
        if ((ibs>1).and.(n>m)) allocate(PnR(ibs,n),PnI(ibs,n))
        do i0=m,n,ibs
          ib=min(ibs,n-i0+1)
          if (ib==1) then
            !Plain column algorithm (original unblocked code)
            i=i0
            do j=1,i-1
              jm=j-1
              x=A(i,j)
              do q=1,jm
                x=x-conjg(A(q,i))*A(q,j)
              enddo
              A(j,i)=conjg(x)
            enddo
            im=i-1
            do q=1,im
              w(q)=A(q,i)*invD(q)
            enddo
            x=A(i,i)
            do q=1,im
              x=x-conjg(A(q,i))*w(q)
            enddo
            x=conjg(x)
            do q=1,im
              A(q,i)=w(q)
            enddo
            if (x==ZERO) then
              ErrorCode=1
              return
            endif
            invD(i)=ONE/x
          else
            !Seed the panel buffers with the matrix elements of rows
            !i0..i0+ib-1 taken from the (untouched) lower triangle
            do k=1,ib
              i=i0+k-1
              do j=1,i-1
                PnR(k,j)=real(A(i,j),wp)
                PnI(k,j)=aimag(A(i,j))
              enddo
            enddo
            !Shared sweep: forward-substitute all panel columns against
            !the completed columns 1..i0-1 (this is where blocking pays)
            do j=1,i0-1
              jm=j-1
              do k=1,ib
                tvr(k)=PnR(k,j)
                tvi(k)=PnI(k,j)
              enddo
              do q=1,jm
                xr=real(A(q,j),wp)
                xi=aimag(A(q,j))
                do k=1,ib
                  tvr(k)=tvr(k)-(PnR(k,q)*xr-PnI(k,q)*xi)
                  tvi(k)=tvi(k)-(PnR(k,q)*xi+PnI(k,q)*xr)
                enddo
              enddo
              do k=1,ib
                PnR(k,j)=tvr(k)
                PnI(k,j)=tvi(k)
              enddo
            enddo
            !Finish the panel columns one by one: substitution against the
            !panel columns completed just before them and the diagonal
            !step, then write the scaled column of the factor back into
            !the upper triangle of A
            do k=1,ib
              i=i0+k-1
              do j=i0,i-1
                jm=j-1
                xr=PnR(k,j)
                xi=PnI(k,j)
                do q=1,jm
                  xr=xr-(PnR(k,q)*real(A(q,j),wp)-PnI(k,q)*aimag(A(q,j)))
                  xi=xi-(PnR(k,q)*aimag(A(q,j))+PnI(k,q)*real(A(q,j),wp))
                enddo
                PnR(k,j)=xr
                PnI(k,j)=xi
              enddo
              im=i-1
              do q=1,im
                w(q)=cmplx(PnR(k,q),-PnI(k,q),wp)*invD(q)
              enddo
              x=A(i,i)
              do q=1,im
                x=x-cmplx(PnR(k,q),PnI(k,q),wp)*w(q)
              enddo
              x=conjg(x)
              do q=1,im
                A(q,i)=w(q)
              enddo
              if (x==ZERO) then
                ErrorCode=1
                return
              endif
              invD(i)=ONE/x
            enddo
          endif
        enddo
      endif
      if (Glob_LDLHF_SModeType==0) then
        do i=m,n
          call MPI_BCAST(A(1:i-1,i),2*(i-1),MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        enddo
        call MPI_BCAST(invD(m:n),2*(n-m+1),MPI_WP,0,MPI_COMM_WORLD,Glob_MPIErrCode)
        call MPI_BCAST(ErrorCode,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
      endif
    else
      !Parallel version (far from being perfect...)
      do i=m,n
        im=i-1
        A(1:im,i)=ZERO
        RowsPerProc=im/Glob_NumOfProcs
        mod_im_Glob_NumOfProcs=mod(im,Glob_NumOfProcs)
        do k=1,RowsPerProc
          jf=k*Glob_NumOfProcs
          jim=jf-Glob_NumOfProcs
          ji=jim+1
          jiR=ji+Glob_ProcID
          A(jiR,i)=conjg(-dot_product(A(1:jim,i),A(1:jim,jiR)))
          call MPI_ALLREDUCE(A(ji:jf,i),w(1:Glob_NumOfProcs),2*Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          A(ji:jf,i)=w(1:Glob_NumOfProcs)
          do j=ji,jf
            jm=j-1
            A(j,i)=A(j,i)+conjg(A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j)))
          enddo
        enddo
        if (mod_im_Glob_NumOfProcs>0) then
          ji=RowsPerProc*Glob_NumOfProcs+1
          jim=ji-1
          jf=im
          jiR=ji+Glob_ProcID
          if (jiR<i) then
            A(jiR,i)=conjg(-dot_product(A(1:jim,i),A(1:jim,jiR)))
          endif
          call MPI_ALLREDUCE(A(ji:jf,i),w(1:mod_im_Glob_NumOfProcs),2*mod_im_Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          A(ji:jf,i)=w(1:mod_im_Glob_NumOfProcs)
          do j=ji,jf
            jm=j-1
            A(j,i)=A(j,i)+conjg(A(i,j)-dot_product(A(ji:jm,i),A(ji:jm,j)))
          enddo
        endif
        !j==1 case
        w(1:im)=A(1:im,i)*invD(1:im)
        y=ZERO
        do k=1+Glob_ProcID,im,Glob_NumOfProcs
          y=y+conjg(A(k,i))*w(k)
        enddo
        y=conjg(y)
        call MPI_ALLREDUCE(y,x,2,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x=A(i,i)-x
        A(1:im,i)=w(1:im)
        if (x==ZERO) then
          ErrorCode=1
          return
        endif
        invD(i)=ONE/x
      enddo
    endif

  end subroutine LDLHF

  subroutine LDLTS(n,A,nA,invD,b,x)
!Subroutine LDLTS finds the solution of a linear system
!A*x=b, where A is a real symmetric matrix, b is a right-hand
!side. Prior to calling LDLTS one must factorize A in the
!form A=L*D*LT using routine LDLTF. The solution of the
!system is found as a result of consecutive solutions of
!systems L*y=b and D*LT*x=y. Vector b is destroyed on exit (only
!if the routine runs in parallel mode).
!  Input parameters :
!       n - The size of the linear system;
!      nA - the leading dimension of A;
!       A - A two-dimensional array that contains the elements
!           of matrix LT (the result of routine LDLTF) in the upper
!           triangle (not including the diagonal elements, which
!           are assumed to be ones). The lower part of A is not used.
!    invD - An array containing inverse values of
!           the diagonal elements of D (the result of LDLTF);
!       b - An array containing the right-hand side. It
!           is also used as workspace so it is destroyed
!           on exit.
!  Output parameters :
!       x - An array containing the solution;

!Arguments
    integer     n,nA
    real(wp) A(nA,n),invD(n),b(n),x(n)
!Local variables
    integer     i,j,k,ji,jf,jim,jfm,RowsPerProc,mod_n_Glob_NumOfProcs,jiR
    real(wp) t

    if (Glob_LDLTS_PMode==0) then
      if (Glob_LDLTS_UseBLAS==0) then
        !Solution of L*y=b
        x(1:n)=b(1:n)
        do j=1,n
          t=x(j)
          do i=1,j-1
            t=t-A(i,j)*x(i)
          enddo
          x(j)=t
        enddo
        !Solution of D*LT*x=y
        x(1:n)=x(1:n)*invD(1:n)
        do j=N,1,-1
          t=x(j)
          do i=j-1,1,-1
            x(i)=x(i)-t*A(i,j)
          enddo
        enddo
      else
        x(1:n)=b(1:n)
        !call BLAS routine DTRSV
        call DTRSV('U','T','U',n,A,nA,x,1)
        x(1:n)=x(1:n)*invD(1:n)
        !call BLAS routine DTRSV
        call DTRSV('U','N','U',n,A,nA,x,1)
      endif
    else
      !Solution of L*y=b
      RowsPerProc=n/Glob_NumOfProcs
      x(1:n)=ZERO
      do i=1,RowsPerProc
        ji=(i-1)*Glob_NumOfProcs+1
        jim=ji-1
        jf=jim+Glob_NumOfProcs
        t=ZERO
        jiR=ji+Glob_ProcID
        do k=1,jim
          t=t-A(k,jiR)*x(k)
        enddo
        x(jiR)=t+b(jiR)
        call MPI_ALLREDUCE(x(ji:jf),b(ji:jf),Glob_NumOfProcs, &
                           MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x(ji:jf)=b(ji:jf)
        do j=ji,jf
          t=x(j)
          do k=ji,j-1
            t=t-A(k,j)*x(k)
          enddo
          x(j)=t
        enddo
      enddo
      mod_n_Glob_NumOfProcs=mod(n,Glob_NumOfProcs)
      if (mod_n_Glob_NumOfProcs>0) then
        jim=RowsPerProc*Glob_NumOfProcs
        ji=jim+1
        jiR=ji+Glob_ProcID
        if (jiR<=n) then
          t=ZERO
          do k=1,jim
            t=t-A(k,jiR)*x(k)
          enddo
          x(jiR)=t+b(jiR)
        endif
        call MPI_ALLREDUCE(x(ji:n),b(ji:n),mod_n_Glob_NumOfProcs, &
                           MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x(ji:n)=b(ji:n)
        do j=ji,n
          t=x(j)
          do k=ji,j-1
            t=t-A(k,j)*x(k)
          enddo
          x(j)=t
        enddo
      endif
      !Solution of D*LT*x=y
      b(1:n)=x(1:n)*invD(1:n)
      x(1:n)=ZERO
      do i=1,RowsPerProc
        jf=n-i*Glob_NumOfProcs+1
        ji=jf-1+Glob_NumOfProcs
        jiR=ji+1+Glob_ProcID
        if (jiR<=n) then
          t=x(jiR)
          do k=ji,1,-1
            x(k)=x(k)-t*A(k,jiR)
          enddo
          call MPI_ALLREDUCE(x(jf:ji),b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs),Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          x(jf:ji)=b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        endif
        x(jf:ji)=x(jf:ji)+b(jf:ji)
        do j=ji,jf,-1
          t=x(j)
          do k=j-1,jf,-1
            x(k)=x(k)-t*A(k,j)
          enddo
        enddo
      enddo
      if (mod_n_Glob_NumOfProcs>0) then
        ji=mod_n_Glob_NumOfProcs
        jiR=ji+1+Glob_ProcID
        if (jiR<=n) then
          t=x(jiR)
          do k=ji,1,-1
            x(k)=x(k)-t*A(k,jiR)
          enddo
          call MPI_ALLREDUCE(x(1:ji),b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs),mod_n_Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          x(1:ji)=b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        endif
        x(1:ji)=x(1:ji)+b(1:ji)
        do j=ji,1,-1
          t=x(j)
          do k=j-1,1,-1
            x(k)=x(k)-t*A(k,j)
          enddo
        enddo
      endif
    endif

  end subroutine LDLTS

  subroutine LDLHS(n,A,nA,invD,b,x)
!Subroutine LDLHS finds the solution of a linear system
!A*x=b, where A is a complex hermitian matrix, b is a right-hand
!side. Prior to calling LDLHS one must factorize A in the
!form A=L*D*LH using routine LDLHF. The solution of the
!system is found as a result of consecutive solutions of
!systems L*y=b and D*LH*x=y. Vector b is destroyed on exit (only
!if the routine runs in parallel mode).
!  Input parameters :
!       n - The size of the linear system;
!      nA - the leading dimension of A;
!       A - A two-dimensional array that contains the elements
!           of matrix LH (the result of routine LDLHF) in the upper
!           triangle (not including the diagonal elements, which
!           are assumed to be ones). The lower part of A is not used.
!    invD - An array containing inverse values of
!           the diagonal elements of D (the result of LDLHF);
!       b - An array containing the right-hand side. It
!           is also used as workspace so it is destroyed
!           on exit.
!  Output parameters :
!       x - An array containing the solution;

!Arguments
    integer        n,nA
    complex(wp) A(nA,n),invD(n),b(n),x(n)
!Local variables
    integer        i,j,k,ji,jf,jim,jfm,RowsPerProc,mod_n_Glob_NumOfProcs,jiR
    complex(wp) t

    if (Glob_LDLHS_PMode==0) then
      if (Glob_LDLHS_UseBLAS==0) then
        !Solution of L*y=b
        x(1:n)=b(1:n)
        do j=1,n
          t=x(j)
          do i=1,j-1
            t=t-conjg(A(i,j))*x(i)
          enddo
          x(j)=t
        enddo
        !Solution of D*LH*x=y
        x(1:n)=x(1:n)*invD(1:n)
        do j=N,1,-1
          t=x(j)
          do i=j-1,1,-1
            x(i)=x(i)-t*A(i,j)
          enddo
        enddo
      else
        x(1:n)=b(1:n)
        !call BLAS routine ZTRSV
        call ZTRSV('U','C','U',n,A,nA,x,1)
        x(1:n)=x(1:n)*invD(1:n)
        !call BLAS routine ZTRSV
        call ZTRSV('U','N','U',n,A,nA,x,1)
      endif
    else
      !Solution of L*y=b
      RowsPerProc=n/Glob_NumOfProcs
      x(1:n)=ZERO
      do i=1,RowsPerProc
        ji=(i-1)*Glob_NumOfProcs+1
        jim=ji-1
        jf=jim+Glob_NumOfProcs
        t=ZERO
        jiR=ji+Glob_ProcID
        do k=1,jim
          t=t-conjg(A(k,jiR))*x(k)
        enddo
        x(jiR)=t+b(jiR)
        call MPI_ALLREDUCE(x(ji:jf),b(ji:jf),2*Glob_NumOfProcs, &
                           MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x(ji:jf)=b(ji:jf)
        do j=ji,jf
          t=x(j)
          do k=ji,j-1
            t=t-conjg(A(k,j))*x(k)
          enddo
          x(j)=t
        enddo
      enddo
      mod_n_Glob_NumOfProcs=mod(n,Glob_NumOfProcs)
      if (mod_n_Glob_NumOfProcs>0) then
        jim=RowsPerProc*Glob_NumOfProcs
        ji=jim+1
        jiR=ji+Glob_ProcID
        if (jiR<=n) then
          t=ZERO
          do k=1,jim
            t=t-conjg(A(k,jiR))*x(k)
          enddo
          x(jiR)=t+b(jiR)
        endif
        call MPI_ALLREDUCE(x(ji:n),b(ji:n),2*mod_n_Glob_NumOfProcs, &
                           MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
        x(ji:n)=b(ji:n)
        do j=ji,n
          t=x(j)
          do k=ji,j-1
            t=t-conjg(A(k,j))*x(k)
          enddo
          x(j)=t
        enddo
      endif
      !Solution of D*LH*x=y
      b(1:n)=x(1:n)*invD(1:n)
      x(1:n)=ZERO
      do i=1,RowsPerProc
        jf=n-i*Glob_NumOfProcs+1
        ji=jf-1+Glob_NumOfProcs
        jiR=ji+1+Glob_ProcID
        if (jiR<=n) then
          t=x(jiR)
          do k=ji,1,-1
            x(k)=x(k)-t*A(k,jiR)
          enddo
          call MPI_ALLREDUCE(x(jf:ji),b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs),2*Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          x(jf:ji)=b(jf+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        endif
        x(jf:ji)=x(jf:ji)+b(jf:ji)
        do j=ji,jf,-1
          t=x(j)
          do k=j-1,jf,-1
            x(k)=x(k)-t*A(k,j)
          enddo
        enddo
      enddo
      if (mod_n_Glob_NumOfProcs>0) then
        ji=mod_n_Glob_NumOfProcs
        jiR=ji+1+Glob_ProcID
        if (jiR<=n) then
          t=x(jiR)
          do k=ji,1,-1
            x(k)=x(k)-t*A(k,jiR)
          enddo
          call MPI_ALLREDUCE(x(1:ji),b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs),2*mod_n_Glob_NumOfProcs, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          x(1:ji)=b(1+Glob_NumOfProcs:ji+Glob_NumOfProcs)
        endif
        x(1:ji)=x(1:ji)+b(1:ji)
        do j=ji,1,-1
          t=x(j)
          do k=j-1,1,-1
            x(k)=x(k)-t*A(k,j)
          enddo
        enddo
      endif
    endif

  end subroutine LDLHS

  subroutine MTMVL(n,A,nA,x,y,w)
!Subroutine MTMVL computes the product of real symmetric matrix
!A and vector x using only the lower triangle of A. The upper
!triangle is not referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A
!           (its lower triangle, including the diagonal);
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of length n)
!  Output parameters :
!       y - An array (vector) containing the result;

    integer      n,nA,j,jm,ji,jf,p,i
    real(wp)  A(nA,n),x(n),y(n),w(n),t,s

    if (Glob_MTMVL_PMode==0) then
      if (Glob_MTMVL_UseBLAS==0) then
        if (wp==8) then
          !Single fused pass over the lower triangle: each column
          !A(j:n,j) is read from memory once and used for both the
          !dot-product contribution to y(j) and the axpy contribution
          !to y(j+1:n), which halves the memory traffic compared with
          !making two separate sweeps over the triangle. Worthwhile
          !only for hardware double precision, where the routine is
          !memory-bandwidth-bound; for wp=10/16 the scalar arithmetic
          !dominates and the two-sweep code below is faster.
          y(1:n)=ZERO
          do j=1,n
            t=x(j)
            s=y(j)+A(j,j)*t
            do i=j+1,n
              s=s+A(i,j)*x(i)
              y(i)=y(i)+A(i,j)*t
            enddo
            y(j)=s
          enddo
        else
          !We use the property A*x=(R+L)*x = (x^T*R^T)^T + (x^T*L^T)^T
          !where L is the lower triangle of A (without the diagonal) and
          !R is the upper triangle (including the diagonal).
          !Computing (x^T*R^T)^T
          do j=1,n
            t=ZERO
            do i=j,n
              t=t+A(i,j)*x(i)
            enddo
            y(j)=t
          enddo
          !Computing (x^T*L^T)^T
          do j=n,1,-1
            t=x(j)
            do i=j+1,n
              y(i)=y(i)+t*A(i,j)
            enddo
          enddo
        endif
      else
        !call BLAS routine DSYMV
        call DSYMV('L',n,ONE,A,nA,x,1,ZERO,y,1)
      endif
    else
      w(1:n)=ZERO
      p=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) p=p+1
      !Computing (x^T*R^T)^T
      ji=1+p*Glob_ProcID
      jf=min(p*(Glob_ProcID+1),n)
      do j=ji,jf
        t=ZERO
        do i=j,n
          t=t+A(i,j)*x(i)
        enddo
        w(j)=t
      enddo
      !Computing (x^T*L^T)^T
      ji=1+p*(Glob_NumOfProcs-Glob_ProcID-1)
      jf=min(p*(Glob_NumOfProcs-Glob_ProcID),n)
      do j=jf,ji,-1
        t=x(j)
        do i=j+1,n
          w(i)=w(i)+t*A(i,j)
        enddo
      enddo
      call MPI_ALLREDUCE(w,y,n,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

  end subroutine MTMVL

  subroutine MHMVL(n,A,nA,x,y,w)
!Subroutine MHMVL computes the product of complex hermitian matrix
!A and vector x using only the lower triangle of A. The upper
!triangle is not referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A
!           (its lower triangle, including the diagonal);
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of length n)
!  Output parameters :
!       y - An array (vector) containing the result;

    integer         n,nA,j,jm,ji,jf,p,i
    complex(wp)  A(nA,n),x(n),y(n),w(n),t,s

    if (Glob_MHMVL_PMode==0) then
      if (Glob_MHMVL_UseBLAS==0) then
        if (wp==8) then
          !Single fused pass over the lower triangle: each column
          !A(j:n,j) is read from memory once and used for both the
          !(conjugated) dot-product contribution to y(j) and the axpy
          !contribution to y(j+1:n), which halves the memory traffic
          !compared with making two separate sweeps over the triangle.
          !Worthwhile only for hardware double precision, where the
          !routine is memory-bandwidth-bound; for wp=10/16 the scalar
          !arithmetic dominates and the two-sweep code below is faster.
          y(1:n)=cmplx(ZERO,ZERO,wp)
          do j=1,n
            t=x(j)
            s=y(j)+conjg(A(j,j))*t
            do i=j+1,n
              s=s+conjg(A(i,j))*x(i)
              y(i)=y(i)+A(i,j)*t
            enddo
            y(j)=s
          enddo
        else
          !We use the property A*x=(R+L)*x = (x^H*R^H)^H + (x^H*L^H)^H
          !where L is the lower triangle of A (without the diagonal) and
          !R is the upper triangle (including the diagonal).
          !Computing (x^H*R^H)^H
          do j=1,n
            t=cmplx(ZERO,ZERO,wp)
            do i=j,n
              t=t+conjg(A(i,j))*x(i)
            enddo
            y(j)=t
          enddo
          !Computing (x^H*L^H)^H
          do j=n,1,-1
            t=x(j)
            do i=j+1,n
              y(i)=y(i)+t*A(i,j)
            enddo
          enddo
        endif
      else
        !call BLAS routine ZHEMV
        call ZHEMV('L',n,cmplx(ONE,ZERO,wp),A,nA,x,1,cmplx(ZERO,ZERO,wp),y,1)
      endif
    else
      w(1:n)=ZERO
      p=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) p=p+1
      !Computing (x^H*R^H)^H
      ji=1+p*Glob_ProcID
      jf=min(p*(Glob_ProcID+1),n)
      do j=ji,jf
        t=cmplx(ZERO,ZERO,wp)
        do i=j,n
          t=t+conjg(A(i,j))*x(i)
        enddo
        w(j)=t
      enddo
      !Computing (x^H*L^H)^H
      ji=1+p*(Glob_NumOfProcs-Glob_ProcID-1)
      jf=min(p*(Glob_NumOfProcs-Glob_ProcID),n)
      do j=jf,ji,-1
        t=x(j)
        do i=j+1,n
          w(i)=w(i)+t*A(i,j)
        enddo
      enddo
      call MPI_ALLREDUCE(w,y,2*n,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

  end subroutine MHMVL

  subroutine MTMV(n,A,nA,x,y,w)
!Subroutine MTMV computes the product of transposed
!real matrix A and vector x.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A;
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of length n)
!  Output parameters :
!       y - An array (vector) containing the result, y=A^T*x

    integer      n,nA,ib,ie,jb,je,j,i,k
    real(wp)  A(nA,n),x(n),y(n),w(n),t
    integer      n2,p,q

    if (Glob_MTMV_PMode==0) then
      if (Glob_MTMV_UseBLAS==0) then
        do i=1,n
          t=ZERO
          do j=1,n
            t=t+A(j,i)*x(j)
          enddo
          y(i)=t
        enddo
      else
        !call BLAS routine DGEMV
        call DGEMV('T',n,n,ONE,A,nA,x,1,ZERO,y,1)
      endif
    else
      !Each process is assigned to deal with some part of
      !matrix A. This part includes all matrix elements between
      !ib,jb and  ie,je (i and j stand for row and index
      !respectively). The word "between" means in natural order,
      !A11, A12, ..., A1n, A21, A22, ... Since we actually need the
      !transposed matrix, the access to matrix elements will therefore
      !be by rows, which suitable for Fortran type of storage. The
      !actual number of elements assigned to a particular process
      !is k. The algorithm divides the labor between processes as
      !uniformly as possible.
      w(1:n)=ZERO
      n2=n*n
      p=n2/Glob_NumOfProcs
      if (mod(n2,Glob_NumOfProcs)/=0) p=p+1
      q=Glob_ProcID*p
      if (q+p>n2) then
        k=max(n2-q,0)
      else
        k=p
      endif
      ib=q/n+1
      ie=min((q+p-1)/n+1,n)
      jb=mod(q,n)+1
      je=mod(q+k-1,n)+1
      if (ie>ib) then
        t=ZERO
        do j=jb,n
          t=t+A(j,ib)*x(j)
        enddo
        w(ib)=t
        do i=ib+1,ie-1
          t=ZERO
          do j=1,n
            t=t+A(j,i)*x(j)
          enddo
          w(i)=t
        enddo
        t=ZERO
        do j=1,je
          t=t+A(j,i)*x(j)
        enddo
        w(ie)=t
      else
        if (ie==ib) then
          t=ZERO
          do j=jb,je
            t=t+A(j,ib)*x(j)
          enddo
          w(ib)=t
        endif
      endif
      call MPI_ALLREDUCE(w,y,n,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end subroutine MTMV

  subroutine MHMV(n,A,nA,x,y,w)
!Subroutine MHMV computes the product of conjugated and transposed
!complex matrix A and vector x.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A;
!      nA - The leading dimension of A;
!       x - An array containing vector x;
!       w - A work array (at least of length n)
!  Output parameters :
!       y - An array (vector) containing the result, y=A^H*x

    integer         n,nA,ib,ie,jb,je,j,i,k
    complex(wp)  A(nA,n),x(n),y(n),w(n),t
    integer         n2,p,q

    if (Glob_MHMV_PMode==0) then
      if (Glob_MHMV_UseBLAS==0) then
        do i=1,n
          t=cmplx(ZERO,ZERO,wp)
          do j=1,n
            t=t+conjg(A(j,i))*x(j)
          enddo
          y(i)=t
        enddo
      else
        !call BLAS routine ZGEMV
        call ZGEMV('C',n,n,cmplx(ONE,ZERO,wp),A,nA,x,1,cmplx(ZERO,ZERO,wp),y,1)
      endif
    else
      !Each process is assigned to deal with some part of
      !matrix A. This part includes all matrix elements between
      !ib,jb and  ie,je (i and j stand for row and index
      !respectively). The word "between" means in natural order,
      !A11, A12, ..., A1n, A21, A22, ...  Since we actually need
      !the transposed (and conjugated) matrix, the access to
      !matrix elements will therefore be by rows, which suitable
      !for Fortran type of storage. The actual number of elements
      !assigned to a particular process is k. The algorithm
      !divides the labor between processes as uniformly as possible.
      w(1:n)=cmplx(ZERO,ZERO,wp)
      n2=n*n
      p=n2/Glob_NumOfProcs
      if (mod(n2,Glob_NumOfProcs)/=0) p=p+1
      q=Glob_ProcID*p
      if (q+p>n2) then
        k=max(n2-q,0)
      else
        k=p
      endif
      ib=q/n+1
      ie=min((q+p-1)/n+1,n)
      jb=mod(q,n)+1
      je=mod(q+k-1,n)+1
      if (ie>ib) then
        t=cmplx(ZERO,ZERO,wp)
        do j=jb,n
          t=t+conjg(A(j,ib))*x(j)
        enddo
        w(ib)=t
        do i=ib+1,ie-1
          t=cmplx(ZERO,ZERO,wp)
          do j=1,n
            t=t+conjg(A(j,i))*x(j)
          enddo
          w(i)=t
        enddo
        t=cmplx(ZERO,ZERO,wp)
        do j=1,je
          t=t+conjg(A(j,i))*x(j)
        enddo
        w(ie)=t
      else
        if (ie==ib) then
          t=cmplx(ZERO,ZERO,wp)
          do j=jb,je
            t=t+conjg(A(j,ib))*x(j)
          enddo
          w(ib)=t
        endif
      endif
      call MPI_ALLREDUCE(w,y,2*n,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end subroutine MHMV

  function VMMTMV(n,A,nA,x)
!Function VMMTMV computes the product x^T*A*x, where x is a real vector,
!x^T is transposed x, and A is a real symmetric matrix.
!Only the lower triangle of A (including the diagonal) is referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A (only the lower
!           triangle is needed);
!      nA - The leading dimension of A;
!       x - An array containing vector x;

!Arguments
    integer        n,nA
    real(wp)    A(nA,n),x(n)
    real(wp)    VMMTMV
!Local variables
    integer      j,k,q,p
    real(wp)  s,sum

    if (Glob_VMMTMV_PMode==0) then
      VMMTMV=ZERO
      do k=1,n
        VMMTMV=VMMTMV+x(k)*x(k)*A(k,k)
        s=ZERO
        do j=k+1,n
          s=s+x(j)*A(j,k)
        enddo
        VMMTMV=VMMTMV+TWO*s*x(k)
      enddo
    else
      sum=ZERO
      q=n/(2*Glob_NumOfProcs)
      do p=1,2*q,2
        k=(p-1)*Glob_NumOfProcs+Glob_ProcID+1
        sum=sum+x(k)*x(k)*A(k,k)
        s=ZERO
        do j=k+1,n
          s=s+x(j)*A(j,k)
        enddo
        sum=sum+TWO*s*x(k)
        k=(p+1)*Glob_NumOfProcs-Glob_ProcID
        sum=sum+x(k)*x(k)*A(k,k)
        s=ZERO
        do j=k+1,n
          s=s+x(j)*A(j,k)
        enddo
        sum=sum+TWO*s*x(k)
      enddo
      p=Glob_ProcID
      do k=2*q*Glob_NumOfProcs+1,n
        if (p==0) then
          sum=sum+x(k)*x(k)*A(k,k)
          p=p+Glob_NumOfProcs
        endif
        s=ZERO
        do j=k+p,n,Glob_NumOfProcs
          s=s+x(j)*A(j,k)
        enddo
        sum=sum+TWO*s*x(k)
        p=j-n-1
      enddo
      call MPI_ALLREDUCE(sum,VMMTMV,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

  end function VMMTMV

  function VMMHMV(n,A,nA,x)
!Function VMMHMV computes the product x^H*A*x, where x is a complex vector,
!x^H is conjugated and transposed x, and A is a complex hermitian matrix.
!Only the lower triangle of A (including the diagonal) is referenced.
!  Input parameters :
!       n - The size of matrix A;
!       A - A two-dimensional array containing matrix A (only the lower
!           triangle is needed);
!      nA - The leading dimension of A;
!       x - An array containing vector x;

!Arguments
    integer        n,nA
    complex(wp) A(nA,n),x(n)
    real(wp)    VMMHMV
!Local variables
    integer      j,k,q,p
    real(wp)  sr,si,sum

    if (Glob_VMMHMV_PMode==0) then
      VMMHMV=ZERO
      do k=1,n
        VMMHMV=VMMHMV+(real(x(k),wp)*real(x(k),wp)+imag(x(k))*imag(x(k)))*real(A(k,k),wp)
        sr=ZERO
        si=ZERO
        do j=k+1,n
          sr=sr+(real(x(j),wp)*real(A(j,k),wp)+imag(x(j))*imag(A(j,k)))
          si=si+(imag(x(j))*real(A(j,k),wp)-real(x(j),wp)*imag(A(j,k)))
        enddo
        VMMHMV=VMMHMV+TWO*(real(x(k),wp)*sr+imag(x(k))*si)
      enddo
    else
      sum=ZERO
      q=n/(2*Glob_NumOfProcs)
      do p=1,2*q,2
        k=(p-1)*Glob_NumOfProcs+Glob_ProcID+1
        sum=sum+(real(x(k),wp)*real(x(k),wp)+imag(x(k))*imag(x(k)))*real(A(k,k),wp)
        sr=ZERO
        si=ZERO
        do j=k+1,n
          sr=sr+(real(x(j),wp)*real(A(j,k),wp)+imag(x(j))*imag(A(j,k)))
          si=si+(imag(x(j))*real(A(j,k),wp)-real(x(j),wp)*imag(A(j,k)))
        enddo
        sum=sum+TWO*(real(x(k),wp)*sr+imag(x(k))*si)
        k=(p+1)*Glob_NumOfProcs-Glob_ProcID
        sum=sum+(real(x(k),wp)*real(x(k),wp)+imag(x(k))*imag(x(k)))*real(A(k,k),wp)
        sr=ZERO
        si=ZERO
        do j=k+1,n
          sr=sr+(real(x(j),wp)*real(A(j,k),wp)+imag(x(j))*imag(A(j,k)))
          si=si+(imag(x(j))*real(A(j,k),wp)-real(x(j),wp)*imag(A(j,k)))
        enddo
        sum=sum+TWO*(real(x(k),wp)*sr+imag(x(k))*si)
      enddo
      p=Glob_ProcID
      do k=2*q*Glob_NumOfProcs+1,n
        if (p==0) then
          sum=sum+(real(x(k),wp)*real(x(k),wp)+imag(x(k))*imag(x(k)))*real(A(k,k),wp)
          p=p+Glob_NumOfProcs
        endif
        sr=ZERO
        si=ZERO
        do j=k+p,n,Glob_NumOfProcs
          sr=sr+(real(x(j),wp)*real(A(j,k),wp)+imag(x(j))*imag(A(j,k)))
          si=si+(imag(x(j))*real(A(j,k),wp)-real(x(j),wp)*imag(A(j,k)))
        enddo
        sum=sum+TWO*(real(x(k),wp)*sr+imag(x(k))*si)
        p=j-n-1
      enddo
      call MPI_ALLREDUCE(sum,VMMHMV,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif

  end function VMMHMV

  function RMaxAbsEl(n,x)
!Function RMaxAbsEl returns the magnitude of the largest by magnitude element
!among the first n elements of real array x
    integer        n,j,ji,jf,k
    real(wp)    x(n)
    real(wp)    RMaxAbsEl,MaxAE
    MaxAE=ZERO
    if (Glob_RMaxAbsEl_PMode==0) then
      do j=1,n
        if (abs(x(j))>MaxAE) MaxAE=abs(x(j))
      enddo
      RMaxAbsEl=MaxAE
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      do j=ji,jf
        if (abs(x(j))>MaxAE) MaxAE=abs(x(j))
      enddo
      call MPI_ALLREDUCE(MaxAE,RMaxAbsEl,1,MPI_WP,MPI_MAX,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function RMaxAbsEl

  function CMaxAbsReOrIm(n,x)
!Function CMaxAbsEl returns the magnitude of the largest real or imaginary part
!among the first n elements of complex array x
    integer  n,j,ji,jf,k,maxj
    complex(wp)  x(n)
    real(wp)     CMaxAbsReOrIm,MaxAE,t
    MaxAE=ZERO
    if (Glob_CMaxAbsReOrIm_PMode==0) then
      do j=1,n
        t=abs(real(x(j),wp))
        if (t>MaxAE) MaxAE=t
        t=abs(imag(x(j)))
        if (t>MaxAE) MaxAE=t
      enddo
      CMaxAbsReOrIm=MaxAE
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      do j=ji,jf
        t=abs(real(x(j),wp))
        if (t>MaxAE) MaxAE=t
        t=abs(imag(x(j)))
        if (t>MaxAE) MaxAE=t
      enddo
      call MPI_ALLREDUCE(MaxAE,CMaxAbsReOrIm,1,MPI_WP,MPI_MAX,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function CMaxAbsReOrIm

  function RDotProd(n,x,y)
!Function RDotProd computes the dot product x^{T}y,
!for n-component real vectors x and y.
    integer           n,k,ji,jf,j
    real(wp)    RDotProd,a
    real(wp)    x(n),y(n)
    real(wp)    DDOT
    if (Glob_RDotProd_PMode==0) then
      if (Glob_RDotProd_UseBLAS==0) then
        RDotProd=ZERO
        do j=1,n
          RDotProd=RDotProd+x(j)*y(j)
        enddo
      else
        !call BLAS function DDOT
        RDotProd=DDOT(n,x,1,y,1)
      endif
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      a=ZERO
      do j=ji,jf
        a=a+x(j)*y(j)
      enddo
      call MPI_ALLREDUCE(a,RDotProd,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function RDotProd

  function CDotProd(n,x,y)
!Function CDotProd computes the dot product x^{H}y,
!for n-component complex vectors x and y.
    integer           n,k,ji,jf,j
    complex(wp)    CDotProd,a
    complex(wp)    x(n),y(n)
    complex(wp)    ZDOTC
    if (Glob_CDotProd_PMode==0) then
      if (Glob_CDotProd_UseBLAS==0) then
        CDotProd=cmplx(ZERO,ZERO,wp)
        do j=1,n
          CDotProd=CDotProd+cmplx(real(x(j),wp)*real(y(j),wp)+imag(x(j))*imag(y(j)), &
                                  real(x(j),wp)*imag(y(j))-imag(x(j))*real(y(j),wp),wp)
        enddo
      else
        !call BLAS function ZDOTC
        CDotProd=ZDOTC(n,x,1,y,1)
      endif
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      a=cmplx(ZERO,ZERO,wp)
      do j=ji,jf
        a=a+cmplx(real(x(j),wp)*real(y(j),wp)+imag(x(j))*imag(y(j)), &
                  real(x(j),wp)*imag(y(j))-imag(x(j))*real(y(j),wp),wp)
      enddo
      call MPI_ALLREDUCE(a,CDotProd,2,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function CDotProd

  function RDotProdItself(n,x)
!Function RDotProdItself computes dot product x^{T}x,
!for an n-component real vector x.
    integer           n,k,ji,jf,j
    real(wp)       RDotProdItself,t
    real(wp)       x(n)
    real(wp)       DDOT
    if (Glob_RDotProdItself_PMode==0) then
      if (Glob_RDotProdItself_UseBLAS==0) then
        RDotProdItself=ZERO
        do j=1,n
          RDotProdItself=RDotProdItself+x(j)*x(j)
        enddo
      else
        !call BLAS function DDOT
        RDotProdItself=DDOT(n,x,1,x,1)
      endif
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      t=ZERO
      do j=ji,jf
        t=t+x(j)*x(j)
      enddo
      call MPI_ALLREDUCE(t,RDotProdItself,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function RDotProdItself

  function CDotProdItself(n,x)
!Function CDotProdItself computes dot product x^{H}x,
!for an n-component complex vector x.
    integer           n,k,ji,jf,j
    real(wp)       CDotProdItself,t
    complex(wp)    x(n)
    complex(wp)    ZDOTC
    if (Glob_CDotProdItself_PMode==0) then
      if (Glob_CDotProdItself_UseBLAS==0) then
        CDotProdItself=ZERO
        do j=1,n
          CDotProdItself=CDotProdItself+real(x(j),wp)*real(x(j),wp)+imag(x(j))*imag(x(j))
        enddo
      else
        !call BLAS function ZDOTC
        CDotProdItself=ZDOTC(n,x,1,x,1)
      endif
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      t=ZERO
      do j=ji,jf
        t=t+real(x(j),wp)*real(x(j),wp)+imag(x(j))*imag(x(j))
      enddo
      call MPI_ALLREDUCE(t,CDotProdItself,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
    endif
  end function CDotProdItself

  function RDotProdQuotient(n,x,y)
!Function RDotProdQuotient computes the quotient (x^{T}y)/(y^{T}y),
!where x and y are n-component real vectors.
    integer        n,k,ji,jf,j
    real(wp)    RDotProdQuotient,a,b,t(2),tt(2)
    real(wp)    x(n),y(n)
    if (Glob_RDotProdQuotient_PMode==0) then
      a=ZERO
      b=ZERO
      do j=1,n
        a=a+x(j)*y(j)
        b=b+y(j)*y(j)
      enddo
      RDotProdQuotient=a/b
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      t(1)=ZERO
      t(2)=ZERO
      do j=ji,jf
        t(1)=t(1)+x(j)*y(j)
        t(2)=t(2)+y(j)*y(j)
      enddo
      call MPI_ALLREDUCE(t,tt,2,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      RDotProdQuotient=tt(1)/tt(2)
    endif
  end function RDotProdQuotient

  function CDotProdQuotient(n,x,y)
!Function CDotProdQuotient computes the quotient (x^{H}y)/(y^{H}y),
!where x and y are n-component complex vectors.
    integer        n,k,ji,jf,j
    complex(wp) CDotProdQuotient
    complex(wp) x(n),y(n)
    complex(wp) a
    real(wp)    r,t(3),tt(3)
    if (Glob_CDotProdQuotient_PMode==0) then
      a=cmplx(ZERO,ZERO,wp)
      r=ZERO
      do j=1,n
        a=a+conjg(x(j))*y(j)
        r=r+real(y(j),wp)*real(y(j),wp)+imag(y(j))*imag(y(j))
      enddo
      CDotProdQuotient=cmplx(real(a,wp)/r,imag(a)/r,wp)
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      a=cmplx(ZERO,ZERO,wp)
      t(1)=ZERO
      do j=ji,jf
        a=a+conjg(x(j))*y(j)
        t(1)=t(1)+real(y(j),wp)*real(y(j),wp)+imag(y(j))*imag(y(j))
      enddo
      t(2)=real(a,wp)
      t(3)=imag(a)
      call MPI_ALLREDUCE(t,tt,3,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      CDotProdQuotient=cmplx(tt(2)/tt(1),tt(3)/tt(1),wp)
    endif
  end function CDotProdQuotient

  subroutine RVScale(n,x,alpha)
!Subroutine RVScale scales real vector x by real constant alpha
    integer      n,j
    real(wp)  x(n)
    real(wp)  alpha

    if (Glob_RVScale_UseBLAS==0) then
      do j=1,n
        x(j)=alpha*x(j)
      enddo
    else
      !call BLAS routine DSCAL
      call DSCAL(n,alpha,x,1)
    endif

  end subroutine RVScale

  subroutine CVScale(n,x,alpha)
!Subroutine RVScale scales complex vector x by complex constant alpha
    integer         n,j
    complex(wp)  x(n)
    complex(wp)  alpha
    real(wp)     ralpha

    if (Glob_CVScale_UseBLAS==0) then
      if (imag(alpha)==ZERO) then
        ralpha=real(alpha,wp)
        do j=1,n
          x(j)=cmplx(ralpha*real(x(j),wp),ralpha*imag(x(j)),wp)
        enddo
      else
        do j=1,n
          x(j)=x(j)*alpha
        enddo
      endif
    else
      !call BLAS routine ZSCAL
      call ZSCAL(n,alpha,x,1)
    endif

  end subroutine CVScale

  function RVDiffEucNorm(n,x,alpha,y)
!Function RVDiffEucNorm returns the Euclidean norm of the difference x-alpha*y,
!where y and x are two real vectors and alpha is a real scalar.
    real(wp)  RVDiffEucNorm
    integer      n,j,k,ji,jf,q
    real(wp)  x(n),y(n),alpha,t

    if (Glob_RVDiffEucNorm_PMode==0) then
      RVDiffEucNorm=ZERO
      do j=1,n
        RVDiffEucNorm=RVDiffEucNorm+(x(j)-alpha*y(j))*(x(j)-alpha*y(j))
      enddo
      RVDiffEucNorm=sqrt(RVDiffEucNorm)
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      t=ZERO
      do j=ji,jf
        t=t+(x(j)-alpha*y(j))*(x(j)-alpha*y(j))
      enddo
      call MPI_ALLREDUCE(t,RVDiffEucNorm,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      RVDiffEucNorm=sqrt(RVDiffEucNorm)
    endif

  end function RVDiffEucNorm

  function CVDiffEucNorm(n,x,alpha,y)
!Function CVDiffEucNorm returns the Euclidean norm of the difference x-alpha*y,
!where y and x are two complex vectors and alpha is a complex scalar.
    real(wp)     CVDiffEucNorm
    integer         n,j,k,ji,jf,q
    complex(wp)  x(n),y(n),alpha,s
    real(wp)     t

    if (Glob_CVDiffEucNorm_PMode==0) then
      CVDiffEucNorm=ZERO
      do j=1,n
        s=x(j)-alpha*y(j)
        CVDiffEucNorm=CVDiffEucNorm+real(s,wp)*real(s,wp)+imag(s)*imag(s)
      enddo
      CVDiffEucNorm=sqrt(CVDiffEucNorm)
    else
      k=n/Glob_NumOfProcs
      if (mod(n,Glob_NumOfProcs)/=0) k=k+1
      ji=1+k*Glob_ProcID
      jf=min(k*(Glob_ProcID+1),n)
      t=cmplx(ZERO,ZERO,wp)
      do j=ji,jf
        s=x(j)-alpha*y(j)
        t=t+real(s,wp)*real(s,wp)+imag(s)*imag(s)
      enddo
      call MPI_ALLREDUCE(t,CVDiffEucNorm,1,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      CVDiffEucNorm=sqrt(CVDiffEucNorm)
    endif

  end function CVDiffEucNorm

  subroutine GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
                     lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
!Subroutine GSEPIIS finds a nondegenerate eigenvalue and the
!corresponding eigenvector of the generalized real symmetric eigenvalue
!problem A*x=lambda*B*x with symmetric matrices A and B using the
!inverse iteration method. The subroutine makes use of a
!known factorization of A-apprlambda*B submatrix of size k-1
!in the form M=A-apprlambda*B=L*D*LT which is done by
!routine LDLTF. The desired eigenvalue must not have any other
!pathologically close eigenvalues (i.e. with relative separation
!of the order MachineEpsilon and less), otherwise the solution process
!cannot separate them from each other and they will be interpreted
!as one degenerate eigenvalue. In the last case the eigenvector will
!lie in the subspace spanned by the corresponding eigenvectors
!and the subroutine can generate ErrorCode=2 upon exit. In
!the case of degenerate eigenvalues the accuracy is not guaranteed,
!however if a group of eigenvalues is much closer
!to apprlambda than other eigenvalues then after sufficient
!number of iterations desirable eigenvector will be found,
!even though there may be generated ErrorCode=2.
!  Input parameters :
!     k-1 - The size of M=A-apprlambda*B submatrix, whose
!           factorization is known (k>=1);
!       n - The size of matrices A and B;
!       M - A two-dimensional array containing known submatrix LH
!           (of size k-1) in the upper triangle (excluding the
!           diagonal, the latter is assumed to consist of ones)
!           and submatrix A-apprlambda*B (of the size k-1) in
!           the lower triangle (including the diagonal);
!      nM - The leading dimension of M;
!    invD - An array containing the inverse values of the
!           diagonal matrix D (of size k-1);
!       B - A two-dimensional array containing matrix B;
!      nB - The leading dimension of B;
!apprlambda - An approximation to the desired eigenvalue. To
!           get a converging process, it is necessary to have
!           the distance between apprlambda and the desired
!           eigenvalue smaller than the distance between apprlambda
!           and any other eigenvalue. The closer apprlambda to
!           the desired eigenvalue the faster the iterative
!           process will converge. However, apprlambda cannot be
!           exactly equal to the real eigenvalue as in this case the
!           matrices the routine deals with are going to be singular;
!       v - An array (at least of size n) where an approximation
!           to the desired eigenvector should be placed. Having
!           a good approximation helps to speed up the iterative process.
!           If such an approximation is not known, it is quite
!           acceptable to use a random vector (but not a zero vector!).
!           This array is also used as a work array. Thus, it is
!           destroyed upon exit;
!       w - A work array (at least of size n);
!     Tol - The desired relative accuracy of the calculations (i.e. the
!           relative accuracy of the eigenvector and the eigenvalue).
!           If this value is set to be small and negative then the
!           routine finds as accurate solution as possible, but not less
!           accurate than abs(Tol). In the case of negative Tol
!           the calculations require at least one more iteration (but potentially
!           can result in many more iterations, up to the iteration limit, MaxIter)
!           and the total computation time may increase substantially.
! MaxIter - The maximum number of iterations allowed;
!SpecifNorm - The parameter that determines how the eigenvector is
!           normalized upon exit. If SpecifNorm=0 then it is normalized
!           in such a way that xT*B*x=1. If SpecifNorm=1 then xT*x=1.
!           Otherwise the eigenvector is normalized in such a way that
!           the largest by magnitude component is one;
!  Output parameters:
!       M - A two-dimensional array whose upper triangle (excluding
!           the diagonal) contains the updated matrix LT (of
!           size n) and the lower triangle (including the
!           diagonal) contains the lower triangle of matrix
!           A-apprlambda*B (of size n);
!    invD - An array containing the inverse values of
!           the diagonal matrix D (of size n);
!  lambda - The eigenvalue found;
!       x - The eigenvector found;
!  RelAcc - A rough estimate of the relative accuracy reached.
!           It may be somewhat inaccurate. A more accurate estimate
!           would require additional calculations and for the
!           sake of efficiency was not implemented.
! NumIter - The number of iterations performed;
!ErrorCode - An error flag. If ErrorCode=0 then the routine
!           finished successfully. There may be the following
!           errors on exit:
!           ErrorCode=1 - matrix A-apprlambda*B is singular
!           or almost singular;
!           ErrorCode=2 - the process did not converge to the
!           accuracy Tol after the maximum number of
!           iterations was performed;

!Arguments
    integer         k,n,nM,nB,MaxIter,SpecifNorm,NumIter,ErrorCode
    real(wp)     M(nM,n),invD(n),B(nB,n),v(n),w(n),x(n)
    real(wp)     apprlambda,Tol,lambda,RelAcc
!Local variables
    real(wp)     NormOfDiff,NormOfDiffPrev,t1,t2,sqrtn
    real(wp)     tc
    logical         notconverged

    ErrorCode=0
    RelAcc=huge(RelAcc)
    NormOfDiffPrev=huge(NormOfDiffPrev)
    NumIter=0
!Updating M=L*D*LH factorization up to size n using routine LDLHF
    call LDLTF(k,n,M,nM,invD,w,ErrorCode)
    if (ErrorCode>0) return
    sqrtn=sqrt(n*ONE)
!Do inverse iterations until the process converges
!with relative accuracy Tol, or until the number
!of iterations exceeds the limit
    notconverged=.true.
    do while ((notconverged).AND.(NumIter<MaxIter))
      if (NumIter/=0) v(1:n)=x(1:n)
      !B is symmetric, so for hardware double precision B*v is computed
      !with the symmetric routine MTMVL, which reads only the lower
      !triangle of B (half the memory traffic of the general
      !transposed-matrix routine MTMV; this routine is bandwidth-bound at
      !wp=8). For wp=10/16 the scalar arithmetic dominates instead and the
      !pure streaming dot products of MTMV are measurably faster, so those
      !builds keep the original call (wp is a compile-time constant, hence
      !the branch costs nothing).
      if (wp==8) then
        call MTMVL(n,B,nB,v,w,x)
      else
        call MTMV(n,B,nB,v,w,x)
      endif
      call LDLTS(n,M,nM,invD,w,x)
      t1=RMaxAbsEl(n,x)
      call RVScale(n,x,1/t1)
      !tc=(x^{T}v)/(v^{T}v)
      tc=RDotProdQuotient(n,x,v)
      !NormOfDiff=sqrt(n)||x-tc*v||
      !old version: NormOfDiff=sqrtn*RVDiffEucNorm(n,x,tc,v)
      NormOfDiff=RVDiffEucNorm(n,x,tc,v)/sqrt(RDotProdItself(n,x))
      if (Tol>ZERO) then
        if (NormOfDiff<=Tol) notconverged=.false.
      else
        if ((NormOfDiff>NormOfDiffPrev).AND.(NormOfDiff<=abs(Tol))) notconverged=.false.
        NormOfDiffPrev=NormOfDiff
      endif
      NumIter=NumIter+1
    enddo
    RelAcc=NormOfDiff
    if (notconverged) ErrorCode=2
!Compute x^{T}Bx
    t1=VMMTMV(n,B,nB,x)
!Compute x^{T}Mx
    t2=VMMTMV(n,M,nM,x)
!Rayleigh quotient
    lambda=(t2/t1)+apprlambda
    select case        (SpecifNorm)
    case (0) !x^{H}Bx=1
      call RVScale(n,x,1/sqrt(t1))
    case (1) !x^{H}x=1
      t1=RDotProdItself(n,x)
      call RVScale(n,x,1/sqrt(t1))
    endselect

  end subroutine GSEPIIS

  subroutine GHEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
                     lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
!Subroutine GHEPIIS finds a nondegenerate eigenvalue and the
!corresponding eigenvector of the generalized complex hermitian eigenvalue
!problem A*x=lambda*B*x with hermitian matrices A and B using the
!inverse iteration method. The subroutine makes use of a
!known factorization of A-apprlambda*B submatrix of size k-1
!in the form A-apprlambda*B=L*D*LH which is done by
!routine LDLHF. The desired eigenvalue must not have any other
!pathologically close eigenvalues (i.e. with relative separation
!of the order MachineEpsilon and less), otherwise the solution process
!cannot separate them from each other and they will be interpreted
!as one degenerate eigenvalue. In the last case the eigenvector will
!lie in the subspace spanned by the corresponding eigenvectors
!and the subroutine can generate ErrorCode=2 upon exit. In
!the case of degenerate eigenvalues the accuracy is not guaranteed,
!however if a group of eigenvalues is much closer
!to apprlambda than other eigenvalues then after sufficient
!number of iterations desirable eigenvector will be found,
!even though there may be generated ErrorCode=2.
!  Input parameters :
!     k-1 - The size of A-apprlambda*B submatrix, whose
!           factorization is known (k>=1);
!       n - The size of matrices A and B;
!       M - A two-dimensional array containing known submatrix LH
!           (of size k-1) in the upper triangle (excluding the
!           diagonal, the latter is assumed to consist of ones)
!           and submatrix A-apprlambda*B (of the size k-1) in
!           the lower triangle (including the diagonal);
!      nM - The leading dimension of M;
!    invD - An array containing the inverse values of the
!           diagonal matrix D (of size k-1);
!       B - A two-dimensional array containing matrix B;
!      nB - The leading dimension of B;
!apprlambda - An approximation to the desired eigenvalue. To
!           get a converging process, it is necessary to have
!           the distance between apprlambda and the desired
!           eigenvalue smaller than the distance between apprlambda
!           and any other eigenvalue. The closer apprlambda to
!           the desired eigenvalue the faster the iterative
!           process will converge. However, apprlambda cannot be
!           exactly equal to the real eigenvalue as in this case the
!           matrices the routine deals with are going to be singular;
!       v - An array (at least of size n) where an approximation
!           to the desired eigenvector should be placed. Having
!           a good approximation helps to speed up the iterative process.
!           If such an approximation is not known, it is quite
!           acceptable to use a random vector (but not a zero vector!).
!           This array is also used as a work array. Thus, it is
!           destroyed upon exit;
!       w - A work array (at least of size n);
!     Tol - The desired relative accuracy of the calculations (i.e. the
!           relative accuracy of the eigenvector and the eigenvalue).
!           If this value is set to be small and negative then the
!           routine finds as accurate solution as possible, but not less
!           accurate than abs(Tol). In the case of negative Tol
!           the calculations require at least one more iteration (but potentially
!           can result in many more iterations, up to the iteration limit, MaxIter)
!           and the total computation time may increase substantially.
! MaxIter - The maximum number of iterations allowed;
!SpecifNorm - The parameter that determines how the eigenvector is
!           normalized upon exit. If SpecifNorm=0 then it is normalized
!           in such a way that xH*B*x=1. If SpecifNorm=1 then xH*x=1.
!           Otherwise the eigenvector is normalized in such a way that
!           the magnitude of the largest either real or imaginary part
!           is one;
!  Output parameters:
!       M - A two-dimensional array whose upper triangle (excluding
!           the diagonal) contains the updated matrix LH (of
!           size n) and the lower triangle (including the
!           diagonal) contains the lower triangle of matrix
!           A-apprlambda*B (of size n);
!    invD - An array containing the inverse values of
!           the diagonal matrix D (of size n);
!  lambda - The eigenvalue found;
!       x - The eigenvector found;
!  RelAcc - A rough estimate of the relative accuracy reached.
!           It may be somewhat inaccurate. A more accurate estimate
!           would require additional calculations and for the
!           sake of efficiency was not implemented.
! NumIter - The number of iterations performed;
!ErrorCode - An error flag. If ErrorCode=0 then the routine
!           finished successfully. There may be the following
!           errors on exit:
!           ErrorCode=1 - matrix A-apprlambda*B is singular
!           or almost singular;
!           ErrorCode=2 - the process did not converge to the
!           accuracy Tol after the maximum number of
!           iterations was performed;

!Arguments
    integer         k,n,nM,nB,MaxIter,SpecifNorm,NumIter,ErrorCode
    complex(wp)  M(nM,n),invD(n),B(nB,n),v(n),w(n),x(n)
    real(wp)     apprlambda,Tol,lambda,RelAcc
!Local variables
    real(wp)     NormOfDiff,NormOfDiffPrev,t1,t2,sqrtn
    complex(wp)  tc
    logical         notconverged

    ErrorCode=0
    RelAcc=huge(RelAcc)
    NormOfDiffPrev=huge(NormOfDiffPrev)
    NumIter=0
!Updating M=L*D*LH factorization up to size n using routine LDLHF
    call LDLHF(k,n,M,nM,invD,w,ErrorCode)
    if (ErrorCode>0) return
    sqrtn=sqrt(2*n*ONE)
!Do inverse iterations until the process converges
!with relative accuracy Tol, or until the number
!of iterations exceeds the limit
    notconverged=.true.
    do while ((notconverged).AND.(NumIter<MaxIter))
      if (NumIter/=0) v(1:n)=x(1:n)
      !B is hermitian, so for hardware double precision B*v is computed
      !with the hermitian routine MHMVL, which reads only the lower
      !triangle of B (half the memory traffic of the general
      !conjugate-transposed-matrix routine MHMV; this routine is
      !bandwidth-bound at wp=8). For wp=10/16 the scalar arithmetic
      !dominates instead and the pure streaming dot products of MHMV are
      !measurably faster, so those builds keep the original call (wp is a
      !compile-time constant, hence the branch costs nothing).
      if (wp==8) then
        call MHMVL(n,B,nB,v,w,x)
      else
        call MHMV(n,B,nB,v,w,x)
      endif
      call LDLHS(n,M,nM,invD,w,x)
      t1=CMaxAbsReOrIm(n,x)
      call CVScale(n,x,cmplx(1/t1,ZERO,wp))
      !tc=(x^{H}v)/(v^{H}v)
      tc=CDotProdQuotient(n,x,v)
      !NormOfDiff=sqrt(2*n)||x-tc*v||
      !old version: NormOfDiff=sqrtn*CVDiffEucNorm(n,x,tc,v)
      NormOfDiff=CVDiffEucNorm(n,x,tc,v)/sqrt(CDotProdItself(n,x))
      if (Tol>ZERO) then
        if (NormOfDiff<=Tol) notconverged=.false.
      else
        if ((NormOfDiff>NormOfDiffPrev).AND.(NormOfDiff<=abs(Tol))) notconverged=.false.
        NormOfDiffPrev=NormOfDiff
      endif
      NumIter=NumIter+1
    enddo
    RelAcc=NormOfDiff
    if (notconverged) ErrorCode=2
!Compute x^{H}Bx
    t1=VMMHMV(n,B,nB,x)
!Compute x^{H}Mx
    t2=VMMHMV(n,M,nM,x)
!Rayleigh quotient
    lambda=t2/t1+apprlambda
    select case        (SpecifNorm)
    case (0) !x^{H}Bx=1
      call CVScale(n,x,cmplx(1/sqrt(t1),ZERO,wp))
    case (1) !x^{H}x=1
      t1=CDotProdItself(n,x)
      call CVScale(n,x,cmplx(1/sqrt(t1),ZERO,wp))
    endselect

  end subroutine GHEPIIS

end module linalg
