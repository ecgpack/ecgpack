module gpu_backend
!CUDA backend for RG_0S matrix-element assembly. Compiled only with -DUSE_CUDA.
!
!The kernels call the original host/device matrix-element routine in matelem.f90.
!The original matrix-element routine reads managed mass and charge module arrays.
!Basis and permutation arrays are staged explicitly on the device.
  use globvars        !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem,   only: MatrixElementsHS_RG_0S, PrecomputeMatrices_L_A
  use cudafor
  use cublas
  use cusolverDn
  implicit none
  private

  !Public surface used by the shared sources (every call site is #ifdef USE_CUDA):
  public :: gpu_backend_init      !collective one-time startup            (main)
  public :: gpu_finalize          !release device/handle at shutdown      (main)
  public :: gpu_eig_active
  public :: gpu_dsygvx
  public :: gpu_active            !.true. => use the CUDA ME backend      (matform)
  public :: gpu_build_HS          !energy matrix-element build            (matform)
  public :: gpu_build_HS_deriv    !energy+gradient matrix-element build   (matform)
  public :: gpu_la_active         !.true. => method-I solve on the GPU    (workproc)
  public :: gpu_la_check          !.true. => shadow-check it on the host  (workproc)
  public :: gpu_la_invalidate     !force a full re-upload of the mirror   (main, workproc)
  public :: gpu_gsepiis           !GPU drop-in for linalg's GSEPIIS       (workproc)

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1
  logical, save :: use_la    = .false.   !ECG_GPU_LA=1: method-I eigensolve (GSEPIIS) on the GPU
  logical, save :: la_check  = .false.   !ECG_GPU_LA_CHECK=1: verify the mirror and shadow the host solve
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call
  logical, save :: ctx_ready = .false.   !device context successfully initialised
  logical, save :: use_determ = .false.  !ECG_DETERM=1: deterministic term reduction
                                         !(ordered, CPU-matching) instead of atomicAdd

  type(cusolverDnHandle), save :: eig_hdl
  logical, save :: eig_created = .false.
  integer, save :: eig_cap_n = 0, eig_cap_work = 0
  real(wp), device, allocatable, save :: dA(:,:), dB(:,:), dW(:), dwork(:)
  integer, device, allocatable, save :: dinfo(:)
  real(wp), allocatable, save :: hA(:,:), hB(:,:)

  !PERSISTENT MATRIX-ELEMENT WORKSPACE.
  !Everything below used to be allocated, filled and freed on EVERY call to the
  !chunk launchers. That is ruinous for the regime production actually runs in:
  !during an optimization sweep the changed functions are permuted to the end of
  !the basis, so one objective evaluation rebuilds only ~K pairs -- yet the old
  !code paid ~11 cudaMalloc + 11 cudaFree (cudaFree SYNCHRONIZES the device) and
  !re-uploaded the whole basis, all NumYHYTerms permutation matrices and the
  !physics constants for that small amount of work. Measured effect: 64 CPU ranks
  !beat one V100 by ~7x on small-K growth (see data/results/rg2p_h2h).
  !Now: allocated once, grown on demand, released in gpu_finalize. The constants
  !are uploaded once; the basis once per build instead of once per chunk.
  real(wp), device, allocatable, save :: d_Lm(:,:,:), d_Am(:,:,:)
  real(wp), device, allocatable, save :: d_YHY(:), d_coeff(:)
  integer,  device, allocatable, save :: d_k(:), d_l(:), d_gl(:)
  real(wp), device, allocatable, save :: d_H(:), d_S(:), d_Dk(:), d_Dl(:)
  real(wp), allocatable, save         :: Lh(:,:,:), Ah(:,:,:)  !host workspace
  integer, save :: cap_basis = 0    !3rd extent of the staged basis arrays
  integer, save :: cap_pairs = 0    !length of the per-pair buffers
  integer, save :: cap_grad  = 0    !length of the gradient slabs
  integer, save :: cap_terms = 0    !NumYHYTerms the constants were staged for
  logical, save :: consts_ready = .false.

  !METHOD-I EIGENSOLVE ON THE GPU (gpu_gsepiis).
  !The device keeps a MIRROR of the host's Glob_H (lower triangle: H-lambda*S,
  !upper triangle: the unit factor U=L' of the LDL' factorization), Glob_S and
  !Glob_invD, with the host leading dimension. The mirror is kept equal to the
  !host arrays by one invariant, borrowed from the host algorithm itself:
  !GSEPIIS(k,n) trusts the factorization of the leading k-1 block, so nothing
  !with both indices < k can have changed since that block was factorized.
  !Each call therefore uploads only the entries with max(row,col) >= d, where
  !d = min(k, la_valid+1) and la_valid is how many leading indices the device
  !is known to hold. Anything that may break the invariant (a host solve, a new
  !BBOP step - the shift changes - or a reallocation) resets la_valid to 0.
  !ECG_GPU_LA_CHECK=1 verifies the untouched leading block bitwise on every call.
  !la_L holds the same factor transposed, L=U' in its strict lower triangle (its
  !upper triangle and diagonal are unused), so that both triangular solves of an
  !inverse iteration run as cuBLAS OP_N: U'y=b as Ly=b, then Ux=y. cuBLAS dtrsv
  !OP_T streams at about half the OP_N rate on a V100 (221 vs 442 GB/s, K=10000).
  !It is valid wherever la_M's upper triangle is (same la_valid).
  type(cublasHandle), save :: la_hdl
  logical, save :: la_created = .false.
  integer, save :: la_ld    = 0          !leading dimension the mirror was allocated for
  integer, save :: la_valid = 0          !leading indices known equal on host and device
  integer, save :: la_bufcap = 0
  real(wp), device, allocatable, save :: la_M(:,:), la_B(:,:), la_L(:,:)
  real(wp), device, allocatable, save :: la_invD(:), la_v(:), la_w(:), la_x(:), la_t(:), la_buf(:)
  integer,  device, allocatable, save :: la_err(:)
  real(wp), allocatable, save         :: la_hbuf(:)
  integer, parameter :: LA_RED = 256     !threads of the single-block reduction kernel (power of 2)

  integer, parameter :: NN  = Glob_AllowedNumOfPseudoParticles
  integer, parameter :: NNP = NN*(NN+1)/2   !max packed vech length; Dk/Dl are 2*np long (np=n(n+1)/2)
  integer, parameter :: CUF_BLK = 128       !threads/block. The term loop STRIDES (see the kernels), so
                                            !block size is independent of NumYHYTerms. The shared
                                            !MatrixElements costs ~254 regs/thread on device; 254*128 <
                                            !65536 (V100 per-block reg limit), so heavy atoms (large
                                            !NumYHYTerms, e.g. Carbon=720) launch instead of failing.
                                            !One-thread-per-term busted this at ~258 terms.

  !Signatures of matform's StoreHS / StoreHSD, matched exactly so those routines
  !can be handed in as callbacks. The chunk loops call them one pair at a time.
  abstract interface
    subroutine store_hs_cb(i,j,Hij,Sij)
      import :: wp
      integer,  intent(in) :: i,j
      real(wp), intent(in) :: Hij,Sij
    end subroutine store_hs_cb
    subroutine store_hsd_cb(i,j,Hij,Sij,Di,Dj)
      import :: wp, Glob_npt
      integer,  intent(in) :: i,j
      real(wp), intent(in) :: Hij,Sij
      real(wp), intent(in) :: Di(2*Glob_npt),Dj(2*Glob_npt)
    end subroutine store_hsd_cb
  end interface

contains

  ! ==========================================================================
  !  (0) FAIL-CLOSED STATUS CHECKING
  !  Every CUDA / cuSOLVER call returns a status that must be inspected: on the
  !  GPU a failure does not raise -- it silently leaves garbage in a buffer that
  !  then flows into the MPI reduction and corrupts EVERY rank. So on any error
  !  we print rich context (op, rank, node, device, error text) and MPI_ABORT --
  !  a single rank limping on with bad data poisons the whole run.
  ! ==========================================================================

  subroutine cuda_check(istat, ctx)
  !Abort the whole job if a CUDA runtime call failed (istat /= cudaSuccess).
    integer,      intent(in) :: istat
    character(*), intent(in) :: ctx
    integer :: dev, ierr, slen
    character(MPI_MAX_PROCESSOR_NAME) :: node
    if (istat == cudaSuccess) return
    node=''; slen=0; call MPI_Get_processor_name(node, slen, ierr)
    dev=-1;  ierr = cudaGetDevice(dev)
    write(*,'(1x,a)')            '================= GPU FATAL ================='
    write(*,'(1x,a,i0,a,a,a,i0)')'rank ',Glob_ProcID,'  node ',trim(node),'  device ',dev
    write(*,'(1x,a,a)')          'context: ',trim(ctx)
    write(*,'(1x,a,i0,a,a)')     'cuda status ',istat,': ',trim(cudaGetErrorString(istat))
    write(*,'(1x,a)')            '============================================='
    flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine cuda_check

  subroutine solver_check(istat, info, ctx)
  !Abort on API failures. The device-side eigenproblem info is returned to
  !the G-method caller so it can apply its existing recovery policy.
    integer,      intent(in) :: istat, info
    character(*), intent(in) :: ctx
    integer :: dev, ierr, slen
    character(MPI_MAX_PROCESSOR_NAME) :: node
    if ((istat == CUSOLVER_STATUS_SUCCESS) .and. (info == 0)) return
    node=''; slen=0; call MPI_Get_processor_name(node, slen, ierr)
    dev=-1;  ierr = cudaGetDevice(dev)
    write(*,'(1x,a)')             '=============== cuSOLVER FATAL =============='
    write(*,'(1x,a,i0,a,a,a,i0)') 'rank ',Glob_ProcID,'  node ',trim(node),'  device ',dev
    write(*,'(1x,a,a)')           'context: ',trim(ctx)
    write(*,'(1x,a,i0,a,i0)')     'cusolver status ',istat,'   info ',info
    write(*,'(1x,a)')             '============================================'
    flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine solver_check

  ! ==========================================================================
  !  (1) ORCHESTRATION -- env selection, lifecycle, chunked MPI builds
  ! ==========================================================================

  subroutine gpu_backend_init()
  !Read the env knobs on rank 0 and broadcast them so every rank agrees (the
  !on/off flag and the chunk size must match, or the per-chunk collective
  !ALLREDUCEs would deadlock/mismatch). Bring up the device context if enabled.
  !Called once, collectively, at startup.
    character(64) :: val
    integer       :: ios, tmp, local_rank, allok
    logical       :: ok
    if (Glob_ProcID==0) then
      call get_environment_variable('ECG_GPU', val, status=ios)
      use_me  = (ios==0) .and. (val(1:1)=='1')
      call get_environment_variable('ECG_GPU_EIG', val, status=ios)
      use_eig = (ios==0) .and. (val(1:1)=='1')
      if (use_eig) use_me = .true.
      call get_environment_variable('ECG_GPU_BATCH', val, status=ios)
      if (ios==0) then
        read(val,*,iostat=ios) tmp
        if ((ios==0).and.(tmp>0)) batch_cap = tmp
      endif
      call get_environment_variable('ECG_DETERM', val, status=ios)
      use_determ = (ios==0).and.(val(1:1)=='1')
      call get_environment_variable('ECG_GPU_LA', val, status=ios)
      use_la = (ios==0).and.(val(1:1)=='1')
      call get_environment_variable('ECG_GPU_LA_CHECK', val, status=ios)
      la_check = (ios==0).and.(val(1:1)=='1')
      if (use_la.and.(Glob_NumOfProcs>1)) then
        !The host method-I routines take MPI-parallel branches at np>1; the GPU
        !solve is written for the single-rank case only.
        write(*,'(1x,a)') '  ECG_GPU_LA ignored: it needs exactly 1 MPI process'
        use_la = .false.
      endif
      if (.not.use_la) la_check = .false.
      write(*,'(1x,a)')             '---------------- GPU backend ----------------'
      write(*,'(1x,a,l1,a,i0,a,l1)')'  ME on GPU=',use_me, &
                                    '   batch=',batch_cap,'   determ=',use_determ
      if (.not.use_me) write(*,'(1x,a)') '  (GPU disabled: matrix elements run on CPU)'
      write(*,'(1x,a,l1)')          '  eig on GPU=',use_eig
      write(*,'(1x,a,l1,a,l1)')     '  method-I solve on GPU=',use_la,'   check=',la_check
      write(*,'(1x,a)')             '---------------------------------------------'
    endif
    call MPI_BCAST(use_eig,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_me,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(batch_cap,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_determ, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_la,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(la_check,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    if (.not.(use_me.or.use_la)) return

    !Bring up the device context on every rank, then agree COLLECTIVELY: if ANY
    !rank failed to acquire a GPU, abort everywhere. Never let some ranks run the
    !CUDA path while others silently fall back -- that diverges the per-chunk
    !ALLREDUCE and hangs the job.
    local_rank = node_local_rank()
    ok = gpu_init(local_rank)
    allok = 0; if (ok) allok = 1
    call MPI_ALLREDUCE(MPI_IN_PLACE, allok, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, Glob_MPIErrCode)
    if (allok /= 1) then
      if (Glob_ProcID==0) write(*,'(1x,a)') &
        'FATAL gpu_backend_init: GPU requested but >=1 rank has no usable device'
      flush(6); call MPI_ABORT(MPI_COMM_WORLD, 1, Glob_MPIErrCode)
    endif
    ctx_ready = .true.
  end subroutine gpu_backend_init

  integer function node_local_rank()
  !Rank index among the ranks that share this physical node -- the correct basis
  !for GPU selection (world rank % nGPUs is wrong across multiple nodes).
    integer :: comm, lrank, ierr
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                             MPI_INFO_NULL, comm, ierr)
    call MPI_Comm_rank(comm, lrank, ierr)
    call MPI_Comm_free(comm, ierr)
    node_local_rank = lrank
  end function node_local_rank

  logical function gpu_active()
    gpu_active = use_me
  end function gpu_active

  logical function gpu_eig_active()
    gpu_eig_active = use_eig .and. ctx_ready
  end function gpu_eig_active

  logical function gpu_la_active()
    gpu_la_active = use_la .and. ctx_ready
  end function gpu_la_active

  logical function gpu_la_check()
    gpu_la_check = la_check
  end function gpu_la_check

  subroutine gpu_la_invalidate()
  !The next gpu_gsepiis call re-uploads the whole mirror (see the comment at
  !la_valid). Cheap: it only resets a counter.
    la_valid = 0
  end subroutine gpu_la_invalidate

  subroutine stage_constants()
  !Upload the per-run invariants: the Y+Y permutation matrices and their
  !coefficients. Mass and scaled charge use shared managed module arrays.
  !The staged invariants do not change during a basis build. The term count is
  !remembered because the symmetry block can be rebuilt mid-run, which would
  !make the staged copy the wrong size.
    integer :: n, nterms
    n = Glob_n; nterms = Glob_NumYHYTerms
    if (consts_ready .and. (cap_terms == nterms)) return
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff)
    allocate(d_YHY(n*n*nterms), d_coeff(nterms))
    call cuda_check(cudaMemcpy(d_YHY,   Glob_YHYMatr,  n*n*nterms), 'stage_constants: YHY matrices')
    call cuda_check(cudaMemcpy(d_coeff, Glob_YHYCoeff, nterms),     'stage_constants: YHY coefficients')
    cap_terms = nterms; consts_ready = .true.
  end subroutine stage_constants

  subroutine stage_basis(n, np, Kmax)
  !Recompute the hoisted per-function matrices and upload them ONCE per build.
  !(They used to be re-uploaded for every chunk, so a 12-chunk build shipped the
  !whole basis 12 times.) The buffers only ever grow.
    integer, intent(in) :: n, np, Kmax
    if (Kmax > cap_basis) then
      if (allocated(Lh))   deallocate(Lh, Ah)
      if (allocated(d_Lm)) deallocate(d_Lm, d_Am)
      allocate(Lh(n,n,Kmax), Ah(n,n,Kmax))
      allocate(d_Lm(n,n,Kmax), d_Am(n,n,Kmax))
      cap_basis = Kmax
    endif
    call PrecomputeMatrices_L_A(np,Kmax,Glob_NonlinParam(1:np,1:Kmax),Lh,Ah)
    call cuda_check(cudaMemcpy(d_Lm,  Lh,  n*n*Kmax), 'stage_basis: L')
    call cuda_check(cudaMemcpy(d_Am,  Ah,  n*n*Kmax), 'stage_basis: A')
  end subroutine stage_basis

  subroutine ensure_pair_capacity(npairs, npt2, need_grad)
  !Grow the per-pair buffers if this build needs more room than the last one.
    integer, intent(in) :: npairs, npt2
    logical, intent(in) :: need_grad
    if (npairs > cap_pairs) then
      if (allocated(d_k)) deallocate(d_k, d_l, d_gl, d_H, d_S)
      allocate(d_k(npairs), d_l(npairs), d_gl(npairs), d_H(npairs), d_S(npairs))
      cap_pairs = npairs
    endif
    if (need_grad) then
      if (npairs*npt2 > cap_grad) then
        if (allocated(d_Dk)) deallocate(d_Dk, d_Dl)
        allocate(d_Dk(npairs*npt2), d_Dl(npairs*npt2))
        cap_grad = npairs*npt2
      endif
    endif
  end subroutine ensure_pair_capacity

  subroutine gpu_build_HS(Nmin,Nmax,store)
  !Energy-only matrix-element build over the basis pair triangle [Nmin,Nmax],
  !processed in fixed-size chunks so neither host nor device memory scales with
  !the whole K(K+1)/2 pair list. Each pair's summed H/S is handed back via the
  !store callback (matform's StoreHS). Every rank uses the same batch_cap, so the
  !per-chunk ALLREDUCE stays collective.
    integer, intent(in)    :: Nmin,Nmax
    procedure(store_hs_cb) :: store
    integer               :: n,np,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable  :: k_list(:),l_list(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    n=Glob_n; np=Glob_np
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (energy): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_constants()                      !once per run
    call stage_basis(n, np, Nmax)               !once per build, not per chunk
    call ensure_pair_capacity(batch, 0, .false.)
    allocate(k_list(batch),l_list(batch),Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    ipair=0; nf=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((nf==batch).or.(ipair==npairs)) then
          call cuf_compute_matelem_batch(Nmax, k_list, l_list, nf, &
              Glob_NumYHYTerms, n, np, Glob_NumOfProcs, Glob_ProcID, &
              Hout, Sout)
          call MPI_ALLREDUCE(Hout,Hout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout,Sout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip))
          enddo
          nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,Hout,Sout,Hout_r,Sout_r)
  end subroutine gpu_build_HS

  subroutine gpu_build_HS_deriv(Nmin,Nmax,store)
  !Energy+gradient build. The (2*np x pairs) gradient buffers are the memory
  !wall at large K, so the same chunking applies. Each pair's H/S and gradient
  !slabs are handed back via the store callback (matform's StoreHSD).
    integer, intent(in)     :: Nmin,Nmax
    procedure(store_hsd_cb) :: store
    integer               :: n,np,npt2,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable  :: k_list(:),l_list(:),grad_l(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    real(wp),allocatable  :: Dkout(:,:),Dlout(:,:),Dkout_r(:,:),Dlout_r(:,:)
    n=Glob_n; np=Glob_np; npt2=np*2
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (grad): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_constants()                      !once per run
    call stage_basis(n, np, Nmax)               !once per build, not per chunk
    call ensure_pair_capacity(batch, npt2, .true.)
    allocate(k_list(batch),l_list(batch),grad_l(batch))
    allocate(Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    allocate(Dkout(npt2,batch),Dlout(npt2,batch),Dkout_r(npt2,batch),Dlout_r(npt2,batch))
    ipair=0; nf=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((l>Glob_nfru).and.(l/=k)) then
          grad_l(nf)=1
        else
          grad_l(nf)=0
        endif
        if ((nf==batch).or.(ipair==npairs)) then
          Dkout(1:npt2,1:nf)=ZERO
          Dlout(1:npt2,1:nf)=ZERO
          call cuf_compute_matelem_deriv_batch(Nmax, k_list, l_list, grad_l, nf, &
              Glob_NumYHYTerms, n, np, Glob_NumOfProcs, Glob_ProcID, &
              Hout, Sout, Dkout, Dlout)
          call MPI_ALLREDUCE(Hout, Hout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout, Sout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Dkout,Dkout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) &
          call MPI_ALLREDUCE(Dlout,Dlout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip), &
                       Dkout_r(1:npt2,ip),Dlout_r(1:npt2,ip))
          enddo
          nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,grad_l)
    deallocate(Hout,Sout,Hout_r,Sout_r)
    deallocate(Dkout,Dlout,Dkout_r,Dlout_r)
  end subroutine gpu_build_HS_deriv

  ! ==========================================================================
  !  (2) DEVICE CODE -- kernels and host launchers
  ! ==========================================================================

  !Energy kernel: one block per (k,l) pair; threads STRIDE over the permutation
  !terms. Each thread calls the shared MatrixElements (grad=false) for its terms.
  attributes(global) subroutine me_energy_kernel(Lm, Am, np, Kmax, &
      k_list, l_list, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, &
      determ, Hout, Sout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid, determ
    real(wp)         :: Lm(n,n,Kmax), Am(n,n,Kmax)
    integer          :: k_list(*), l_list(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    real(wp)         :: Hout(*), Sout(*)
    real(wp), shared :: sh_H, sh_S
    real(wp), shared :: shm(*)                !dynamic: 2*nterms, deterministic term buffer
    integer  :: pair_idx, j, qbase, k0, l0, istat
    real(wp) :: Hkl, Skl, coeff, hsum, ssum
    real(wp) :: dDk(2*NNP), dDl(2*NNP)   !unused energy-path gradient outputs

    pair_idx = blockIdx%x
    qbase = (pair_idx-1)*nterms - 1          !matches Fortran q=(i-1)*nterms-1
    k0    = k_list(pair_idx)                 !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)

    if (determ /= 0) then
      !DETERMINISTIC: write each term's contribution to shm(j) (non-owned -> 0.0, exact),
      !then thread 1 sums shm(1..nterms) IN ORDER -- bit-identical to the CPU term loop,
      !so the result is reproducible run-to-run and matches the CPU sum order.
      do j = threadIdx%x, nterms, blockDim%x
        if (mod(qbase+j, nprocs) == procid) then
          call MatrixElementsHS_RG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                              Am(1,1,k0), Am(1,1,l0), &
                              YHYMatr((j-1)*n*n+1), &
                              Hkl, Skl, dDk, dDl, .false., .false.)
          coeff = YHYCoeff(j)
          shm(j)        = coeff*Hkl
          shm(nterms+j) = coeff*Skl
        else
          shm(j)        = 0.0_wp
          shm(nterms+j) = 0.0_wp
        endif
      enddo
      call syncthreads()
      if (threadIdx%x==1) then
        hsum=0.0_wp; ssum=0.0_wp
        do j=1,nterms
          hsum = hsum + shm(j)
          ssum = ssum + shm(nterms+j)
        enddo
        Hout(pair_idx)=hsum; Sout(pair_idx)=ssum
      endif
    else
      !FAST (default): atomicAdd reduction -- non-deterministic term order.
      if (threadIdx%x==1) then
        sh_H=0.0_wp; sh_S=0.0_wp
      endif
      call syncthreads()
      do j = threadIdx%x, nterms, blockDim%x   !STRIDED: block size independent of nterms
        if (mod(qbase+j, nprocs) == procid) then
          call MatrixElementsHS_RG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                              Am(1,1,k0), Am(1,1,l0), &
                              YHYMatr((j-1)*n*n+1), &
                              Hkl, Skl, dDk, dDl, .false., .false.)
          coeff = YHYCoeff(j)
          istat = atomicadd(sh_H, coeff*Hkl)
          istat = atomicadd(sh_S, coeff*Skl)
        endif
      enddo
      call syncthreads()
      if (threadIdx%x==1) then
        Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
      endif
    endif
  end subroutine me_energy_kernel

  !Host launcher for the energy build. The basis and the constants are already
  !resident (stage_basis / stage_constants), so this only ships the chunk's pair
  !list, launches, and collects H/S. No device allocation happens here.
  subroutine cuf_compute_matelem_batch(Kmax, k_list, l_list, npairs, &
      nterms, n, np, nprocs, procid, Hout, Sout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    integer,  intent(in) :: k_list(npairs), l_list(npairs)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    integer  :: blk, determ, shmem

    call cuda_check(cudaMemcpy(d_k, k_list, npairs), 'cuf_compute_matelem_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l, l_list, npairs), 'cuf_compute_matelem_batch: H2D l list')

    blk = min(nterms, CUF_BLK)
    determ = 0; shmem = 0
    if (use_determ) then; determ = 1; shmem = 2*nterms*8; endif   !2*nterms real*8 dynamic shared
    call me_energy_kernel<<<npairs, blk, shmem>>>(d_Lm, d_Am, np, Kmax, d_k, d_l, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        determ, d_H, d_S)
    call cuda_check(cudaGetLastError(),      'cuf_compute_matelem_batch: energy kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'cuf_compute_matelem_batch: energy kernel execution')

    call cuda_check(cudaMemcpy(Hout, d_H, npairs), 'cuf_compute_matelem_batch: D2H H')
    call cuda_check(cudaMemcpy(Sout, d_S, npairs), 'cuf_compute_matelem_batch: D2H S')
  end subroutine cuf_compute_matelem_batch

  !Energy+gradient kernel: one block per (k,l) pair; threads STRIDE over the terms.
  !Each thread calls the shared MatrixElements (grad_k=true, grad_l=flag) and
  !atomic-accumulates H,S and the gradient slabs Dk,Dl into shared memory.
  attributes(global) subroutine me_grad_kernel(Lm, Am, np, Kmax, &
      k_list, l_list, grad_l_flag, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, &
      Hout, Sout, Dkout, Dlout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid
    real(wp)         :: Lm(n,n,Kmax), Am(n,n,Kmax)
    integer          :: k_list(*), l_list(*), grad_l_flag(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    real(wp)         :: Hout(*), Sout(*), Dkout(*), Dlout(*)   !Dk/Dl slabs (2*np,npairs) flattened
    real(wp), shared :: sh_H, sh_S
    real(wp), shared :: sh_Dk(2*NNP), sh_Dl(2*NNP)
    integer  :: pair_idx, j, qbase, k0, l0, istat, q, npt2, nthr
    real(wp) :: Hkl, Skl, coeff, Dk(2*NNP), Dl(2*NNP)
    logical  :: gl

    pair_idx = blockIdx%x
    nthr     = blockDim%x
    npt2     = 2*np
    if (threadIdx%x==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    do q=threadIdx%x,npt2,nthr    !stride-init the shared gradient slabs
      sh_Dk(q)=0.0_wp; sh_Dl(q)=0.0_wp
    enddo
    call syncthreads()

    qbase = (pair_idx-1)*nterms - 1
    k0    = k_list(pair_idx)                !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)
    gl    = (grad_l_flag(pair_idx)==1)
    do j = threadIdx%x, nterms, nthr        !STRIDED: block size independent of nterms
      if (mod(qbase+j, nprocs) == procid) then
        call MatrixElementsHS_RG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                            Am(1,1,k0), Am(1,1,l0), &
                            YHYMatr((j-1)*n*n+1), &
                            Hkl, Skl, Dk, Dl, .true., gl)
        coeff = YHYCoeff(j)
        istat = atomicadd(sh_H, coeff*Hkl)
        istat = atomicadd(sh_S, coeff*Skl)
        do q=1,npt2
          istat = atomicadd(sh_Dk(q), coeff*Dk(q))
          if (gl) istat = atomicadd(sh_Dl(q), coeff*Dl(q))
        enddo
      endif
    enddo
    call syncthreads()

    if (threadIdx%x==1) then
      Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
    endif
    do q=threadIdx%x,npt2,nthr    !stride-write the slabs back to global
      Dkout((pair_idx-1)*npt2+q)=sh_Dk(q)
      Dlout((pair_idx-1)*npt2+q)=sh_Dl(q)
    enddo
  end subroutine me_grad_kernel

  !Host launcher for the gradient build. Same as above: resident basis and
  !constants, so only the chunk's pair list moves in and the results move out.
  subroutine cuf_compute_matelem_deriv_batch(Kmax, k_list, l_list, grad_l_flag, &
      npairs, nterms, n, np, nprocs, procid, Hout, Sout, Dkout, Dlout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    integer,  intent(in) :: k_list(npairs), l_list(npairs), grad_l_flag(npairs)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp), intent(out):: Dkout(2*np*npairs), Dlout(2*np*npairs)
    integer  :: npt2, blk

    npt2   = 2*np

    call cuda_check(cudaMemcpy(d_k,  k_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l,  l_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D l list')
    call cuda_check(cudaMemcpy(d_gl, grad_l_flag, npairs), 'cuf_compute_matelem_deriv_batch: H2D grad flags')

    blk = min(nterms, CUF_BLK)
    call me_grad_kernel<<<npairs, blk>>>(d_Lm, d_Am, np, Kmax, d_k, d_l, d_gl, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        d_H, d_S, d_Dk, d_Dl)
    call cuda_check(cudaGetLastError(),      'cuf_compute_matelem_deriv_batch: grad kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'cuf_compute_matelem_deriv_batch: grad kernel execution')

    call cuda_check(cudaMemcpy(Hout,  d_H,  npairs),      'cuf_compute_matelem_deriv_batch: D2H H')
    call cuda_check(cudaMemcpy(Sout,  d_S,  npairs),      'cuf_compute_matelem_deriv_batch: D2H S')
    call cuda_check(cudaMemcpy(Dkout, d_Dk, npt2*npairs), 'cuf_compute_matelem_deriv_batch: D2H Dk')
    call cuda_check(cudaMemcpy(Dlout, d_Dl, npt2*npairs), 'cuf_compute_matelem_deriv_batch: D2H Dl')
  end subroutine cuf_compute_matelem_deriv_batch

  !Device lifecycle. gpu_init selects device = node-local-rank % nGPUs and returns
  !.false. (does NOT abort) on any failure, so the caller decides collectively.
  logical function gpu_init(local_rank)
    integer, intent(in)  :: local_rank
    integer :: ndev, dev, istat
    type(cudaDeviceProp) :: prop
    gpu_init = .false.
    istat = cudaGetDeviceCount(ndev)
    if ((istat /= cudaSuccess) .or. (ndev <= 0)) then
      write(*,'(1x,a,i0,a)') 'gpu_init: rank ',Glob_ProcID,' found no CUDA device'
      return
    endif
    dev   = mod(local_rank, ndev)
    istat = cudaSetDevice(dev)                  ; if (istat /= cudaSuccess) return
    istat = cudaGetDeviceProperties(prop, dev)  ; if (istat /= cudaSuccess) return
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') 'gpu_init: rank ',Glob_ProcID, &
        ' (node-local ',local_rank,') -> device ',dev,' ',trim(prop%name)
    gpu_init = .true.
  end function gpu_init

  subroutine gpu_finalize()
    integer :: istat
  !Release application-owned device and host buffers at shutdown.
  !Called once at shutdown (main, before MPI_FINALIZE).
    if (.not. ctx_ready) return
    if (eig_created) then
      istat = cusolverDnDestroy(eig_hdl)
      eig_created = .false.
    endif
    if (allocated(dA)) deallocate(dA)
    if (allocated(dB)) deallocate(dB)
    if (allocated(dW)) deallocate(dW)
    if (allocated(dwork)) deallocate(dwork)
    if (allocated(dinfo)) deallocate(dinfo)
    if (allocated(hA)) deallocate(hA)
    if (allocated(hB)) deallocate(hB)
    eig_cap_n=0; eig_cap_work=0
    !persistent matrix-element workspace
    if (allocated(d_Lm))  deallocate(d_Lm, d_Am)
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff)
    if (allocated(d_k))   deallocate(d_k, d_l, d_gl, d_H, d_S)
    if (allocated(d_Dk))  deallocate(d_Dk, d_Dl)
    if (allocated(Lh))    deallocate(Lh, Ah)
    cap_basis=0; cap_pairs=0; cap_grad=0; cap_terms=0; consts_ready=.false.
    !method-I mirror
    if (la_created) then
      istat = cublasDestroy(la_hdl)
      la_created = .false.
    endif
    if (allocated(la_M))    deallocate(la_M, la_B, la_L, la_invD, la_v, la_w, la_x, la_t, la_err)
    if (allocated(la_buf))  deallocate(la_buf)
    if (allocated(la_hbuf)) deallocate(la_hbuf)
    la_ld=0; la_valid=0; la_bufcap=0
  end subroutine gpu_finalize

  subroutine gpu_dsygvx(jobz_vec, n, H, ldH, S, ldS, iwhich, eval, Z, info_out)
    integer  :: jobz_vec, n, ldH, ldS, iwhich, info_out
    real(wp) :: H(ldH,n), S(ldS,n), eval, Z(*)
    !Handle, device buffers and host staging are MODULE-scope + save (see top):
    !host staging is heap so it is NOT a K*K automatic array on the stack, and the
    !handle can be destroyed in gpu_finalize.
    integer  :: lwork, meig, jobz, uplo
    real(wp) :: hZ(n), hW1(1)    !hZ is a length-n vector (small), safe on the stack
    integer  :: hinfo(1)

    if (.not. eig_created) then
      call solver_check(cusolverDnCreate(eig_hdl), 0, 'gpu_dsygvx: cusolverDnCreate')
      allocate(dinfo(1)); eig_created = .true.
    endif
    if (n /= eig_cap_n) then
      if (allocated(dA)) deallocate(dA, dB, dW)
      if (allocated(hA)) deallocate(hA, hB)
      allocate(dA(n,n), dB(n,n), dW(n))
      allocate(hA(n,n), hB(n,n))       !exact leading dimension for cudaMemcpy
      eig_cap_n = n
    endif
    hA(1:n,1:n) = H(1:n,1:n)           !host de-stride (ld -> n), then explicit H2D
    hB(1:n,1:n) = S(1:n,1:n)
    call cuda_check(cudaMemcpy(dA, hA, n*n), 'gpu_dsygvx: H2D A')   !array-section device
    call cuda_check(cudaMemcpy(dB, hB, n*n), 'gpu_dsygvx: H2D B')   !assignment ICEs nvfortran

    jobz = CUSOLVER_EIG_MODE_NOVECTOR
    if (jobz_vec /= 0) jobz = CUSOLVER_EIG_MODE_VECTOR
    uplo = CUBLAS_FILL_MODE_UPPER

    call solver_check(cusolverDnDsygvdx_bufferSize(eig_hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, lwork), 0, 'gpu_dsygvx: Dsygvdx_bufferSize')
    if (lwork > eig_cap_work) then
      if (allocated(dwork)) deallocate(dwork)
      allocate(dwork(lwork)); eig_cap_work = lwork
    endif
    call solver_check(cusolverDnDsygvdx(eig_hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, dwork, lwork, dinfo(1)), 0, 'gpu_dsygvx: Dsygvdx')
    call cuda_check(cudaDeviceSynchronize(), 'gpu_dsygvx: solve sync')

    call cuda_check(cudaMemcpy(hinfo, dinfo, 1), 'gpu_dsygvx: D2H info')
    info_out = hinfo(1)
    !Return numerical failures to the existing G-method recovery path.
    !API/status failures still abort in solver_check.
    if (info_out /= 0) return
    if (meig < 1) then
      info_out = -999
      return
    endif
    call cuda_check(cudaMemcpy(hW1, dW, 1), 'gpu_dsygvx: D2H eigenvalue')
    eval = hW1(1)               !selected eigenvalue (first of meig)
    if (jobz_vec /= 0) then
      call cuda_check(cudaMemcpy(hZ, dA, n), 'gpu_dsygvx: D2H eigenvector')
      Z(1:n) = hZ(1:n)          !first column of A (the eigenvector)
    endif
  end subroutine gpu_dsygvx


  ! ==========================================================================
  !  (3) METHOD-I EIGENSOLVE ON THE GPU  (see the comment at la_valid)
  ! ==========================================================================

  subroutine blas_check(istat, ctx)
  !Abort the whole job if a cuBLAS call failed.
    integer,      intent(in) :: istat
    character(*), intent(in) :: ctx
    integer :: ierr
    if (istat == CUBLAS_STATUS_SUCCESS) return
    write(*,'(1x,a)')        '================ cuBLAS FATAL ==============='
    write(*,'(1x,a,i0,a,a)') 'rank ',Glob_ProcID,'  context: ',trim(ctx)
    write(*,'(1x,a,i0)')     'cublas status ',istat
    write(*,'(1x,a)')        '============================================='
    flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine blas_check

  subroutine la_prepare(ld)
  !Create the cuBLAS handle once; (re)allocate the mirror when the host leading
  !dimension changes (a reallocated basis), which also invalidates it.
    integer, intent(in) :: ld
    if (.not.la_created) then
      call blas_check(cublasCreate(la_hdl), 'la_prepare: cublasCreate')
      la_created = .true.
    endif
    if (ld /= la_ld) then
      if (allocated(la_M)) deallocate(la_M, la_B, la_L, la_invD, la_v, la_w, la_x, la_t, la_err)
      allocate(la_M(ld,ld), la_B(ld,ld), la_L(ld,ld))
      allocate(la_invD(ld), la_v(ld), la_w(ld), la_x(ld), la_t(ld), la_err(1))
      la_ld = ld
      la_valid = 0
    endif
  end subroutine la_prepare

  subroutine la_upload_rows(A, dA, ld, d, n)
  !Copy A(d:n, 1:d-1) (the strided part of the trailing region: rows d..n of
  !the leading columns) to the device: pack on the host, one copy, scatter.
    integer,  intent(in) :: ld, d, n
    real(wp), intent(in) :: A(ld,*)
    real(wp), device     :: dA(ld,*)
    integer :: r, c, j, cnt
    r = n-d+1; c = d-1; cnt = r*c
    if (cnt <= 0) return
    if (cnt > la_bufcap) then
      if (allocated(la_buf))  deallocate(la_buf)
      if (allocated(la_hbuf)) deallocate(la_hbuf)
      allocate(la_buf(cnt), la_hbuf(cnt))
      la_bufcap = cnt
    endif
    do j = 1, c
      la_hbuf((j-1)*r+1:j*r) = A(d:n, j)
    enddo
    call cuda_check(cudaMemcpy(la_buf, la_hbuf, cnt), 'la_upload_rows: H2D')
    call la_scatter_kernel<<<(cnt+255)/256, 256>>>(dA, ld, la_buf, r, c, d)
    call cuda_check(cudaGetLastError(), 'la_upload_rows: scatter launch')
  end subroutine la_upload_rows

  attributes(global) subroutine la_scatter_kernel(A, ld, buf, r, c, d)
  !A(d+i-1, j) = buf(i,j), i=1..r, j=1..c
    integer, value :: ld, r, c, d
    real(wp) :: A(ld,*), buf(*)
    integer  :: idx, i, j
    idx = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (idx > r*c) return
    j = (idx-1)/r + 1
    i = idx - (j-1)*r
    A(d+i-1, j) = buf(idx)
  end subroutine la_scatter_kernel

  attributes(global) subroutine ldl_finish_kernel(A, L, ld, invD, i, err)
  !Last step of one row of LDLTF (linalg.f90, unblocked path). On entry
  !A(1:i-1,i) holds y with U'y = A(i,1:i-1). On exit A(1:i-1,i) = y*invD,
  !the same values in row i of L (L=U'), invD(i) = 1/(A(i,i) - sum y*y*invD);
  !err(1)=i marks a zero pivot, as the host's ErrorCode=1. One block; the sum
  !is a fixed-order tree, so the result is reproducible run to run.
    integer, value :: ld, i
    real(wp) :: A(ld,*), L(ld,*), invD(*)
    integer  :: err(*)
    real(wp), shared :: red(LA_RED)
    integer  :: q, t, h
    real(wp) :: s, y, wq
    if (err(1) /= 0) return
    t = threadIdx%x
    s = 0.0_wp
    do q = t, i-1, LA_RED
      y  = A(q,i)
      wq = y*invD(q)
      s  = s + y*wq
      A(q,i) = wq
      L(i,q) = wq
    enddo
    red(t) = s
    call syncthreads()
    h = LA_RED/2
    do while (h >= 1)
      if (t <= h) red(t) = red(t) + red(t+h)
      call syncthreads()
      h = h/2
    enddo
    if (t == 1) then
      s = A(i,i) - red(1)
      if (s == 0.0_wp) then
        err(1) = i
      else
        invD(i) = 1.0_wp/s
      endif
    endif
  end subroutine ldl_finish_kernel

  attributes(global) subroutine vmul_kernel(n, x, d)
  !x = x*d elementwise (the D^-1 step of LDLTS)
    integer, value :: n
    real(wp) :: x(*), d(*)
    integer  :: i
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i <= n) x(i) = x(i)*d(i)
  end subroutine vmul_kernel

  subroutine la_verify_leading(c, M, ld, invD, B)
  !ECG_GPU_LA_CHECK: the leading c x c blocks of the mirror must equal the host
  !arrays bit for bit (that is the invariant the partial upload relies on).
    integer,  intent(in) :: c, ld
    real(wp), intent(in) :: M(ld,*), invD(*), B(ld,*)
    real(wp), allocatable :: h(:,:), hv(:)
    integer :: nbad, ierr, j
    allocate(h(ld,c), hv(c))
    call cuda_check(cudaMemcpy(h, la_M, ld*c), 'la_verify: D2H M')
    nbad = count(h(1:c,1:c) /= M(1:c,1:c))
    call cuda_check(cudaMemcpy(h, la_L, ld*c), 'la_verify: D2H L')
    do j = 1, c-1                                  !L(i,j) = U(j,i), i > j
      nbad = nbad + count(h(j+1:c,j) /= M(j,j+1:c))
    enddo
    call cuda_check(cudaMemcpy(h, la_B, ld*c), 'la_verify: D2H B')
    nbad = nbad + count(h(1:c,1:c) /= B(1:c,1:c))
    call cuda_check(cudaMemcpy(hv, la_invD, c), 'la_verify: D2H invD')
    nbad = nbad + count(hv(1:c) /= invD(1:c))
    deallocate(h, hv)
    if (nbad /= 0) then
      write(*,'(1x,a,i0,a,i0)') 'FATAL ECG_GPU_LA_CHECK: device mirror differs from host in ', &
        nbad,' entries of the leading block of order ',c
      flush(6); call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
  end subroutine la_verify_leading

  subroutine gpu_gsepiis(k,n,M,nM,invD,B,nB,apprlambda,v,w,Tol, &
                         lambda,x,RelAcc,MaxIter,SpecifNorm,NumIter,ErrorCode)
  !Drop-in for GSEPIIS (linalg.f90): same arguments, same algorithm, same
  !outputs, including what GSEPIIS leaves in M (factor columns k..n), invD,
  !v, w and x. The O(n^2) work runs on the GPU mirror; only the trailing part
  !of M and B goes up and only columns k..n of M come back.
  !  LDLTF   rows k..n: dcopy (row -> column), dtrsv(L=U'), ldl_finish_kernel
  !  MTMVL   dsymv on the lower triangle of B
  !  LDLTS   dtrsv(U'), x*invD, dtrsv(U)
  !  VMMTMV  dsymv + ddot (lower triangles of B and M)
    integer  :: k,n,nM,nB,MaxIter,SpecifNorm,NumIter,ErrorCode
    real(wp) :: M(nM,*),invD(*),B(nB,*),v(*),w(*),x(*)
    real(wp) :: apprlambda,Tol,lambda,RelAcc
    integer  :: d, i, im, imax, herr(1), ierr, ncols, k1   !not "nb": that IS the dummy nB
    real(wp) :: NormOfDiff,NormOfDiffPrev,t1,t2,tc,a,b2,xx,hv(1),tlast
    logical  :: notconverged

    if (nB /= nM) then
      write(*,'(1x,a)') 'FATAL gpu_gsepiis: H and S must share the leading dimension'
      flush(6); call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    k1 = max(k,1)
    if (la_check) then
      !surface any earlier asynchronous error here rather than in a later copy
      call cuda_check(cudaDeviceSynchronize(), 'gpu_gsepiis: entry sync')
      call cuda_check(cudaGetLastError(), 'gpu_gsepiis: entry last error')
    endif
    call la_prepare(nM)
    d = min(k1, la_valid+1)
    if (la_check .and. (d > 1)) call la_verify_leading(d-1, M, nM, invD, B)

    !Upload everything with max(row,col) >= d: columns d..n (contiguous) and
    !rows d..n of the leading columns (strided). invD(d..k-1): factor rows the
    !host holds but the device may not.
    if (d <= n) then
      call cuda_check(cudaMemcpy(la_M(1,d), M(1,d), nM*(n-d+1)), 'gpu_gsepiis: H2D M columns')
      call cuda_check(cudaMemcpy(la_B(1,d), B(1,d), nM*(n-d+1)), 'gpu_gsepiis: H2D B columns')
      call la_upload_rows(M, la_M, nM, d, n)
      call la_upload_rows(B, la_B, nM, d, n)
    endif
    if (k1 > d) then
      call cuda_check(cudaMemcpy(la_invD(d), invD(d), min(k1,n+1)-d), 'gpu_gsepiis: H2D invD')
      !factor columns d..k1-1 came from the host: L(d:k1-1, 1:k1-1) = U(1:k1-1, d:k1-1)'
      !(this also fills L's unused upper part of those rows, harmlessly)
      call blas_check(cublasDgeam(la_hdl, CUBLAS_OP_T, CUBLAS_OP_N, k1-d, k1-1, 1.0_wp, la_M(1,d), nM, &
                                  0.0_wp, la_L(d,1), nM, la_L(d,1), nM), 'gpu_gsepiis: geam L')
    endif
    call cuda_check(cudaMemcpy(la_v, v, n), 'gpu_gsepiis: H2D v')

    ErrorCode=0
    RelAcc=huge(RelAcc)
    NormOfDiffPrev=huge(NormOfDiffPrev)
    NormOfDiff=huge(NormOfDiff)
    NumIter=0

    !LDLTF(k,n): rows k..n of the factor, one row at a time, as the host's
    !unblocked path (the one the m=n single-row update always takes)
    herr(1)=0
    call cuda_check(cudaMemcpy(la_err, herr, 1), 'gpu_gsepiis: H2D err')
    do i = k1, n
      im = i-1
      if (im > 0) then
        call blas_check(cublasDcopy_v2(la_hdl, im, la_M(i,1), nM, la_M(1,i), 1), 'LDLTF dcopy')
        call blas_check(cublasDtrsv_v2(la_hdl, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, &
                                       im, la_L, nM, la_M(1,i), 1), 'LDLTF dtrsv')
      endif
      call ldl_finish_kernel<<<1, LA_RED>>>(la_M, la_L, nM, la_invD, i, la_err)
    enddo
    call cuda_check(cudaGetLastError(), 'gpu_gsepiis: ldl_finish launch')
    call cuda_check(cudaMemcpy(herr, la_err, 1), 'gpu_gsepiis: D2H err')
    if (herr(1) /= 0) then
      !Zero pivot at row herr(1): the host returns at once with columns k..i
      !written and invD(k..i-1) set. Reproduce exactly that, resync next time.
      i = herr(1)
      call cuda_check(cudaMemcpy(M(1,k1), la_M(1,k1), nM*(i-k1+1)), 'gpu_gsepiis: D2H M (error)')
      if (i > k1) call cuda_check(cudaMemcpy(invD(k1), la_invD(k1), i-k1), 'gpu_gsepiis: D2H invD (error)')
      ErrorCode=1
      la_valid=0
      return
    endif

    !inverse iteration
    notconverged=.true.
    do while ((notconverged).and.(NumIter<MaxIter))
      if (NumIter/=0) call blas_check(cublasDcopy_v2(la_hdl, n, la_x, 1, la_v, 1), 'II copy x->v')
      call blas_check(cublasDsymv_v2(la_hdl, CUBLAS_FILL_MODE_LOWER, n, 1.0_wp, la_B, nM, &
                                     la_v, 1, 0.0_wp, la_w, 1), 'II dsymv B')
      call blas_check(cublasDcopy_v2(la_hdl, n, la_w, 1, la_x, 1), 'II copy w->x')
      call blas_check(cublasDtrsv_v2(la_hdl, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, &
                                     n, la_L, nM, la_x, 1), 'II dtrsv L')
      call vmul_kernel<<<(n+255)/256, 256>>>(n, la_x, la_invD)
      call blas_check(cublasDtrsv_v2(la_hdl, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, &
                                     n, la_M, nM, la_x, 1), 'II dtrsv U')
      call blas_check(cublasIdamax_v2(la_hdl, n, la_x, 1, imax), 'II idamax')
      call cuda_check(cudaMemcpy(hv, la_x(imax), 1), 'gpu_gsepiis: D2H max element')
      t1 = abs(hv(1))
      tlast = t1
      call blas_check(cublasDscal_v2(la_hdl, n, 1.0_wp/t1, la_x, 1), 'II dscal')
      call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_v, 1, a),  'II dot x.v')
      call blas_check(cublasDdot_v2(la_hdl, n, la_v, 1, la_v, 1, b2), 'II dot v.v')
      tc = a/b2
      call blas_check(cublasDcopy_v2(la_hdl, n, la_x, 1, la_t, 1), 'II copy x->t')
      call blas_check(cublasDaxpy_v2(la_hdl, n, -tc, la_v, 1, la_t, 1), 'II axpy')
      call blas_check(cublasDnrm2_v2(la_hdl, n, la_t, 1, t2), 'II nrm2')
      call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_x, 1, xx), 'II dot x.x')
      NormOfDiff = t2/sqrt(xx)
      if (Tol>0.0_wp) then
        if (NormOfDiff<=Tol) notconverged=.false.
      else
        if ((NormOfDiff>NormOfDiffPrev).and.(NormOfDiff<=abs(Tol))) notconverged=.false.
        NormOfDiffPrev=NormOfDiff
      endif
      NumIter=NumIter+1
    enddo
    RelAcc=NormOfDiff
    if (notconverged) ErrorCode=2
    !Rayleigh quotient: t1=x'Bx, t2=x'(H-lambda S)x. The last iteration solved
    !M*x = w/tlast (w = B*v), so x'Mx = x'w/tlast without another sweep over
    !M. Both forms carry rounding of order eps*|M|*|x|^2 in t2.
    call blas_check(cublasDsymv_v2(la_hdl, CUBLAS_FILL_MODE_LOWER, n, 1.0_wp, la_B, nM, &
                                   la_x, 1, 0.0_wp, la_t, 1), 'RQ dsymv B')
    call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_t, 1, t1), 'RQ dot B')
    if (NumIter > 0) then
      call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_w, 1, t2), 'RQ dot x.w')
      t2 = t2/tlast
    else
      call blas_check(cublasDsymv_v2(la_hdl, CUBLAS_FILL_MODE_LOWER, n, 1.0_wp, la_M, nM, &
                                     la_x, 1, 0.0_wp, la_t, 1), 'RQ dsymv M')
      call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_t, 1, t2), 'RQ dot M')
    endif
    lambda=(t2/t1)+apprlambda
    select case (SpecifNorm)
    case (0) !x'Bx=1
      call blas_check(cublasDscal_v2(la_hdl, n, 1.0_wp/sqrt(t1), la_x, 1), 'norm B')
    case (1) !x'x=1
      call blas_check(cublasDdot_v2(la_hdl, n, la_x, 1, la_x, 1, xx), 'norm dot')
      call blas_check(cublasDscal_v2(la_hdl, n, 1.0_wp/sqrt(xx), la_x, 1), 'norm 2')
    end select

    !Back to the host: exactly what GSEPIIS itself would have changed.
    ncols = n-k1+1
    if (ncols > 0) then
      call cuda_check(cudaMemcpy(M(1,k1), la_M(1,k1), nM*ncols), 'gpu_gsepiis: D2H M columns')
      call cuda_check(cudaMemcpy(invD(k1), la_invD(k1), ncols),   'gpu_gsepiis: D2H invD')
    endif
    call cuda_check(cudaMemcpy(x, la_x, n), 'gpu_gsepiis: D2H x')
    call cuda_check(cudaMemcpy(v, la_v, n), 'gpu_gsepiis: D2H v')
    call cuda_check(cudaMemcpy(w, la_w, n), 'gpu_gsepiis: D2H w')
    la_valid = n
  end subroutine gpu_gsepiis

end module gpu_backend
