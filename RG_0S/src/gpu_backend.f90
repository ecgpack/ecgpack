module gpu_backend
!The complete GPU backend for RG_0S -- both the high-level orchestration AND the
!native CUDA Fortran device code, in one module. Compiled ONLY with -DUSE_CUDA,
!so the shared Fortran sources (globvars, matform, workproc, main) carry only a
!thin, #ifdef-guarded dispatch and are otherwise identical to a CPU-only build.
!
!Two layers live here, separated by the section banners below:
!  (1) ORCHESTRATION -- run-time selection (env vars), the device-context
!      lifecycle, and the MPI-parallel chunked matrix-element builds. The chunk
!      loops call back into matform's StoreHS/StoreHSD through procedure arguments,
!      so this module does not depend on matform (which would be circular) yet
!      reuses its storage routines.
!  (2) DEVICE CODE -- the attributes(global) matrix-element kernels (which call the
!      SHARED MatrixElements device routine in matelem.f90, so the physics lives
!      once, not duplicated), the device lifecycle, and the cuSOLVER eigensolver.
!      Physics constants are passed as kernel arguments -- device code cannot read
!      host module globals.
  use globvars        !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem,   only: MatrixElementsHS_RG_0S, Precompute_LAMA, Precompute_PermutationMaps
  use cudafor
  use cublas
  use cusolverDn
  implicit none
  private

  !Public surface used by the shared sources (every call site is #ifdef USE_CUDA):
  public :: gpu_backend_init      !collective one-time startup            (main)
  public :: gpu_finalize          !release device/handle at shutdown      (main)
  public :: gpu_active            !.true. => use the CUDA ME backend      (matform)
  public :: gpu_eig_active        !.true. => use the cuSOLVER eigensolver (workproc)
  public :: gpu_build_HS          !energy matrix-element build            (matform)
  public :: gpu_build_HS_deriv    !energy+gradient matrix-element build   (matform)
  public :: gpu_dsygvx            !cuSOLVER eigensolve                    (workproc)

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1 (implies use_me)
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call
  logical, save :: ctx_ready = .false.   !device context successfully initialised
  logical, save :: use_determ = .false.  !ECG_DETERM=1: deterministic term reduction
                                         !(ordered, CPU-matching) instead of atomicAdd

  !cuSOLVER eigensolver state. Module-scope + save so gpu_finalize can release
  !the handle/buffers and so host staging is heap (NOT K*K automatic stack arrays,
  !which overflow at large K). All capacity-cached across calls.
  type(cusolverDnHandle), save   :: eig_hdl
  logical, save                  :: eig_created = .false.
  integer, save                  :: eig_cap_n = 0, eig_cap_work = 0
  real(wp), device, allocatable, save :: dA(:,:), dB(:,:), dW(:), dwork(:)
  integer,  device, allocatable, save :: dinfo(:)
  real(wp), allocatable, save    :: hA(:,:), hB(:,:)   !host staging (heap)

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
  real(wp), device, allocatable, save :: d_Lm(:,:,:), d_Am(:,:,:), d_MAm(:,:,:)
  real(wp), device, allocatable, save :: d_YHY(:), d_coeff(:)
  integer,  device, allocatable, save :: d_perm(:), d_iperm(:)
  logical,  device, allocatable, save :: d_isperm(:)
  real(wp), device, allocatable, save :: d_mass(:,:), d_chargeM(:,:)
  integer,  device, allocatable, save :: d_k(:), d_l(:), d_gl(:)
  real(wp), device, allocatable, save :: d_H(:), d_S(:), d_Dk(:), d_Dl(:)
  real(wp), allocatable, save         :: Lh(:,:,:), Ah(:,:,:), MAh(:,:,:)  !host workspace
  integer, save :: cap_basis = 0    !3rd extent of the staged basis arrays
  integer, save :: cap_pairs = 0    !length of the per-pair buffers
  integer, save :: cap_grad  = 0    !length of the gradient slabs
  integer, save :: cap_terms = 0    !NumYHYTerms the constants were staged for
  logical, save :: consts_ready = .false.

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
  !Abort if a cuSOLVER call failed. cuSOLVER has TWO error channels: the routine
  !status (istat) AND a device-side info flag (like LAPACK INFO: info>0 => the
  !info-th leading minor of S is not positive definite -- a near-linear-dependent
  !ECG basis; info<0 => bad argument).
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
      if (use_eig) use_me = .true.       !GPU eigensolver implies GPU matrix elements
      call get_environment_variable('ECG_GPU_BATCH', val, status=ios)
      if (ios==0) then
        read(val,*,iostat=ios) tmp
        if ((ios==0).and.(tmp>0)) batch_cap = tmp
      endif
      call get_environment_variable('ECG_DETERM', val, status=ios)
      use_determ = (ios==0).and.(val(1:1)=='1')
      write(*,'(1x,a)')             '---------------- GPU backend ----------------'
      write(*,'(1x,a,l1,a,l1,a,i0,a,l1)')'  ME on GPU=',use_me,'   eig on GPU=',use_eig, &
                                    '   batch=',batch_cap,'   determ=',use_determ
      if (.not.use_me) write(*,'(1x,a)') '  (GPU disabled: ME and eigensolve run on CPU)'
      write(*,'(1x,a)')             '---------------------------------------------'
    endif
    call MPI_BCAST(use_me,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_eig,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(batch_cap,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_determ, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    if (.not.use_me) return

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
    gpu_eig_active = use_eig
  end function gpu_eig_active

  subroutine stage_constants()
  !Upload the per-run invariants: the Y+Y permutation matrices and their
  !coefficients, the mass matrix, and the scaled charge matrix. None of these
  !change while a basis is being built, so this runs once. The term count is
  !remembered because the symmetry block can be rebuilt mid-run, which would
  !make the staged copy the wrong size.
    integer :: n, nterms
    integer, allocatable :: perm(:,:), iperm(:,:)
    logical, allocatable :: isperm(:)
    n = Glob_n; nterms = Glob_NumYHYTerms
    if (consts_ready .and. (cap_terms == nterms)) return
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff, d_perm, d_iperm, d_isperm)
    allocate(d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_perm(n*nterms), d_iperm(n*nterms), d_isperm(nterms))
    allocate(perm(n,nterms), iperm(n,nterms), isperm(nterms))
    call Precompute_PermutationMaps(n,nterms,Glob_YHYMatr(1:n,1:n,1:nterms), &
                                    perm,iperm,isperm)
    if (.not.allocated(d_mass)) allocate(d_mass(n,n), d_chargeM(0:n,0:n))
    call cuda_check(cudaMemcpy(d_YHY,   Glob_YHYMatr,  n*n*nterms), 'stage_constants: YHY matrices')
    call cuda_check(cudaMemcpy(d_coeff, Glob_YHYCoeff, nterms),     'stage_constants: YHY coefficients')
    call cuda_check(cudaMemcpy(d_perm, perm, n*nterms),             'stage_constants: permutation maps')
    call cuda_check(cudaMemcpy(d_iperm, iperm, n*nterms),           'stage_constants: inverse permutation maps')
    call cuda_check(cudaMemcpy(d_isperm, isperm, nterms),           'stage_constants: permutation flags')
    call cuda_check(cudaMemcpy(d_mass,  Glob_MassMatrix, n*n),      'stage_constants: mass matrix')
    call cuda_check(cudaMemcpy(d_chargeM, Glob_ScaledPseudoChargeMatrix, (n+1)*(n+1)), &
                    'stage_constants: charge matrix')
    cap_terms = nterms; consts_ready = .true.
  end subroutine stage_constants

  subroutine stage_basis(n, np, Kmax)
  !Recompute the hoisted per-function matrices and upload them ONCE per build.
  !(They used to be re-uploaded for every chunk, so a 12-chunk build shipped the
  !whole basis 12 times.) The buffers only ever grow.
    integer, intent(in) :: n, np, Kmax
    if (Kmax > cap_basis) then
      if (allocated(Lh))   deallocate(Lh, Ah, MAh)
      if (allocated(d_Lm)) deallocate(d_Lm, d_Am, d_MAm)
      allocate(Lh(n,n,Kmax), Ah(n,n,Kmax), MAh(n,n,Kmax))
      allocate(d_Lm(n,n,Kmax), d_Am(n,n,Kmax), d_MAm(n,n,Kmax))
      cap_basis = Kmax
    endif
    call Precompute_LAMA(n,np,Kmax,Glob_NonlinParam(1:np,1:Kmax), &
                         Glob_MassMatrix(1:n,1:n),Lh,Ah,MAh)
    call cuda_check(cudaMemcpy(d_Lm,  Lh,  n*n*Kmax), 'stage_basis: L')
    call cuda_check(cudaMemcpy(d_Am,  Ah,  n*n*Kmax), 'stage_basis: A')
    call cuda_check(cudaMemcpy(d_MAm, MAh, n*n*Kmax), 'stage_basis: MA')
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
              Glob_PiRaised3n2, Hout, Sout)
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
              Glob_PiRaised3n2, Hout, Sout, Dkout, Dlout)
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
  !  (2) DEVICE CODE -- kernels, host launchers, lifecycle, eigensolver
  !  (all private except gpu_dsygvx)
  ! ==========================================================================

  !Energy kernel: one block per (k,l) pair; threads STRIDE over the permutation
  !terms. Each thread calls the shared MatrixElements (grad=false) for its terms.
  attributes(global) subroutine me_energy_kernel(Lm, Am, MAm, np, Kmax, &
      k_list, l_list, YHYMatr, YHYCoeff, perm, iperm, isperm, nterms, n, nprocs, procid, &
      mass, chargeM, sqrtpi, pir3n2, determ, Hout, Sout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid, determ
    real(wp)         :: Lm(n,n,Kmax), Am(n,n,Kmax), MAm(n,n,Kmax)
    integer          :: k_list(*), l_list(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    integer          :: perm(*), iperm(*)
    logical          :: isperm(*)
    real(wp)         :: mass(n,n), chargeM(0:n,0:n)
    real(wp), value  :: sqrtpi, pir3n2
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
          call MatrixElementsHS_RG_0S(n, np, Lm(1,1,k0), Lm(1,1,l0), &
                              Am(1,1,k0), Am(1,1,l0), MAm(1,1,k0), &
                              YHYMatr((j-1)*n*n+1), perm((j-1)*n+1), iperm((j-1)*n+1), isperm(j), &
                              mass, chargeM, &
                              sqrtpi, pir3n2, &
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
          call MatrixElementsHS_RG_0S(n, np, Lm(1,1,k0), Lm(1,1,l0), &
                              Am(1,1,k0), Am(1,1,l0), MAm(1,1,k0), &
                              YHYMatr((j-1)*n*n+1), perm((j-1)*n+1), iperm((j-1)*n+1), isperm(j), &
                              mass, chargeM, &
                              sqrtpi, pir3n2, &
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
      nterms, n, np, nprocs, procid, pir3n2, Hout, Sout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    integer,  intent(in) :: k_list(npairs), l_list(npairs)
    real(wp), intent(in) :: pir3n2
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp) :: sqrtpi
    integer  :: blk, determ, shmem

    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    call cuda_check(cudaMemcpy(d_k, k_list, npairs), 'cuf_compute_matelem_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l, l_list, npairs), 'cuf_compute_matelem_batch: H2D l list')

    blk = min(nterms, CUF_BLK)
    determ = 0; shmem = 0
    if (use_determ) then; determ = 1; shmem = 2*nterms*8; endif   !2*nterms real*8 dynamic shared
    call me_energy_kernel<<<npairs, blk, shmem>>>(d_Lm, d_Am, d_MAm, np, Kmax, d_k, d_l, &
        d_YHY, d_coeff, d_perm, d_iperm, d_isperm, nterms, n, nprocs, procid, &
        d_mass, d_chargeM, sqrtpi, pir3n2, determ, d_H, d_S)
    call cuda_check(cudaGetLastError(),      'cuf_compute_matelem_batch: energy kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'cuf_compute_matelem_batch: energy kernel execution')

    call cuda_check(cudaMemcpy(Hout, d_H, npairs), 'cuf_compute_matelem_batch: D2H H')
    call cuda_check(cudaMemcpy(Sout, d_S, npairs), 'cuf_compute_matelem_batch: D2H S')
  end subroutine cuf_compute_matelem_batch

  !Energy+gradient kernel: one block per (k,l) pair; threads STRIDE over the terms.
  !Each thread calls the shared MatrixElements (grad_k=true, grad_l=flag) and
  !atomic-accumulates H,S and the gradient slabs Dk,Dl into shared memory.
  attributes(global) subroutine me_grad_kernel(Lm, Am, MAm, np, Kmax, &
      k_list, l_list, grad_l_flag, YHYMatr, YHYCoeff, perm, iperm, isperm, nterms, n, nprocs, procid, &
      mass, chargeM, sqrtpi, pir3n2, &
      Hout, Sout, Dkout, Dlout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid
    real(wp)         :: Lm(n,n,Kmax), Am(n,n,Kmax), MAm(n,n,Kmax)
    integer          :: k_list(*), l_list(*), grad_l_flag(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    integer          :: perm(*), iperm(*)
    logical          :: isperm(*)
    real(wp)         :: mass(n,n), chargeM(0:n,0:n)
    real(wp), value  :: sqrtpi, pir3n2
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
        call MatrixElementsHS_RG_0S(n, np, Lm(1,1,k0), Lm(1,1,l0), &
                            Am(1,1,k0), Am(1,1,l0), MAm(1,1,k0), &
                            YHYMatr((j-1)*n*n+1), perm((j-1)*n+1), iperm((j-1)*n+1), isperm(j), &
                            mass, chargeM, &
                            sqrtpi, pir3n2, &
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
      npairs, nterms, n, np, nprocs, procid, pir3n2, Hout, Sout, Dkout, Dlout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    integer,  intent(in) :: k_list(npairs), l_list(npairs), grad_l_flag(npairs)
    real(wp), intent(in) :: pir3n2
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp), intent(out):: Dkout(2*np*npairs), Dlout(2*np*npairs)
    real(wp) :: sqrtpi
    integer  :: npt2, blk

    npt2   = 2*np
    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    call cuda_check(cudaMemcpy(d_k,  k_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l,  l_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D l list')
    call cuda_check(cudaMemcpy(d_gl, grad_l_flag, npairs), 'cuf_compute_matelem_deriv_batch: H2D grad flags')

    blk = min(nterms, CUF_BLK)
    call me_grad_kernel<<<npairs, blk>>>(d_Lm, d_Am, d_MAm, np, Kmax, d_k, d_l, d_gl, &
        d_YHY, d_coeff, d_perm, d_iperm, d_isperm, nterms, n, nprocs, procid, &
        d_mass, d_chargeM, sqrtpi, pir3n2, &
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
  !Release the cuSOLVER handle and application-owned device/host buffers.
  !Called once at shutdown (main, before MPI_FINALIZE).
  !Do not call cudaDeviceReset here. NVHPC 26.5/CUDA 13.2 performs additional
  !CUDA Fortran runtime cleanup after this routine; destroying the context first
  !makes that cleanup fail with CUDA error 709 (CONTEXT_IS_DESTROYED). Process
  !teardown releases the context after the compiler runtime has finished.
    integer :: istat
    if (.not. ctx_ready) return
    if (eig_created) then
      istat = cusolverDnDestroy(eig_hdl); eig_created = .false.
    endif
    if (allocated(dA))    deallocate(dA)
    if (allocated(dB))    deallocate(dB)
    if (allocated(dW))    deallocate(dW)
    if (allocated(dwork)) deallocate(dwork)
    if (allocated(dinfo)) deallocate(dinfo)
    if (allocated(hA))    deallocate(hA)
    if (allocated(hB))    deallocate(hB)
    !persistent matrix-element workspace
    if (allocated(d_Lm))  deallocate(d_Lm, d_Am, d_MAm)
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff, d_perm, d_iperm, d_isperm)
    if (allocated(d_mass))deallocate(d_mass, d_chargeM)
    if (allocated(d_k))   deallocate(d_k, d_l, d_gl, d_H, d_S)
    if (allocated(d_Dk))  deallocate(d_Dk, d_Dl)
    if (allocated(Lh))    deallocate(Lh, Ah, MAh)
    cap_basis=0; cap_pairs=0; cap_grad=0; cap_terms=0; consts_ready=.false.
  end subroutine gpu_finalize

  !cuSOLVER generalized symmetric eigensolver. Solves H x = lambda S x (itype 1),
  !returns the iwhich-th eigenvalue and, if jobz_vec/=0, its eigenvector. Handle +
  !device buffers cached across calls.
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
    if (n > eig_cap_n) then
      if (allocated(dA)) deallocate(dA, dB, dW)
      if (allocated(hA)) deallocate(hA, hB)
      allocate(dA(n,n), dB(n,n), dW(n))
      allocate(hA(n,n), hB(n,n))       !heap host staging (capacity-cached)
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
    !info/=0 => cuSOLVER failed (e.g. S not positive definite); meig<1 => no
    !eigenpair in the requested range. Both are fail-closed.
    call solver_check(CUSOLVER_STATUS_SUCCESS, info_out, 'gpu_dsygvx: eigensolver info flag')
    if (meig < 1) call solver_check(CUSOLVER_STATUS_SUCCESS, -999, 'gpu_dsygvx: meig<1 (no eigenpair)')
    call cuda_check(cudaMemcpy(hW1, dW, 1), 'gpu_dsygvx: D2H eigenvalue')
    eval = hW1(1)               !selected eigenvalue (first of meig)
    if (jobz_vec /= 0) then
      call cuda_check(cudaMemcpy(hZ, dA, n), 'gpu_dsygvx: D2H eigenvector')
      Z(1:n) = hZ(1:n)          !first column of A (the eigenvector)
    endif
  end subroutine gpu_dsygvx

end module gpu_backend
