module gpu_backend
!CUDA backend for PG_0S: matrix-element assembly (energy and gradients) and the
!cuSOLVER eigensolver. Compiled only with -DUSE_CUDA.
!
!It follows the RG backend (RG_1P/src/gpu_backend.f90) with one structural
!difference. MatrixElementsOpt evaluates ALL symmetry terms of a pair in one
!call, and matform deals the work out per PAIR (not per pair and term). So here
!one thread handles one pair: it calls the same host/device MatrixElementsOpt
!the CPU path uses, then sums the terms in their natural order. That keeps the
!CPU's MPI split, and the term sum is deterministic (no atomics).
!
!Performance note. MatrixElementsOpt declares ~110 arrays whose extent is the
!run-time Glob_NumYHYTerms. In device code nvfortran turns each of them into a
!device malloc/free per call, which makes the matrix elements far slower on the
!GPU than on one CPU core (measured on a V100 for HD: 40-75x slower; see the
!matelem repository, branch_pg). The backend is correct, and the eigensolver
!path is unaffected, but ECG_GPU=1 alone is not a speed-up for PG_0S today.
!Those mallocs also need a device heap far larger than CUDA's 8 MB default:
!gpu_init sets it (ECG_GPU_HEAP, default 1 GiB). An exhausted heap makes the
!kernel fault, which cuda_check turns into an abort.
!
!PrecomputeMatrices_MA_L_A stays host-only: the bra-side matrices of every
!function are built on the host once per build (stage_basis) and uploaded.
  use globvars        !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem,   only: MatrixElementsOpt, PrecomputeMatrices_MA_L_A
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

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call
  integer(kind=cuda_count_kind), save :: heap_bytes = 1073741824_cuda_count_kind  !ECG_GPU_HEAP
  logical, save :: ctx_ready = .false.   !device context successfully initialised

  type(cusolverDnHandle), save :: eig_hdl
  logical, save :: eig_created = .false.
  integer, save :: eig_cap_n = 0, eig_cap_work = 0
  real(wp), device, allocatable, save :: dA(:,:), dB(:,:), dW(:), dwork(:)
  integer, device, allocatable, save :: dinfo(:)
  real(wp), allocatable, save :: hA(:,:), hB(:,:)

  !PERSISTENT MATRIX-ELEMENT WORKSPACE (as in the RG backend): allocated once,
  !grown on demand, released in gpu_finalize. The symmetry terms are uploaded
  !once per run; the basis once per build, not once per chunk.
  !  bra side  : Lm, Am, MAm (n,n,K) from PrecomputeMatrices_MA_L_A
  !  ket side  : par (npt,K) nonlinear parameters, pwr (K) r1 powers
  !  terms     : sym (nterms,n,n) in MatrixElementsOpt's term-first layout
  !  per pair  : MatrixElementsOpt's term-indexed outputs live in global memory
  !              (bS, bH: nterms x pairs; bdSk..bdHl: nterms x npt x pairs) and
  !              are passed as contiguous columns
  real(wp), device, allocatable, save :: d_Lm(:,:,:), d_Am(:,:,:), d_MAm(:,:,:)
  real(wp), device, allocatable, save :: d_par(:,:)
  integer,  device, allocatable, save :: d_pwr(:)
  real(wp), device, allocatable, save :: d_sym(:,:,:), d_coeff(:)
  integer,  device, allocatable, save :: d_k(:), d_l(:), d_gf(:)
  real(wp), device, allocatable, save :: d_H(:), d_S(:), d_Dk(:,:), d_Dl(:,:)
  real(wp), device, allocatable, save :: d_bS(:,:), d_bH(:,:)
  real(wp), device, allocatable, save :: d_bdSk(:,:,:), d_bdSl(:,:,:), d_bdHk(:,:,:), d_bdHl(:,:,:)
  real(wp), allocatable, save         :: Lh(:,:,:), Ah(:,:,:), MAh(:,:,:)  !host workspace
  integer, save :: cap_basis = 0    !3rd extent of the staged basis arrays
  integer, save :: cap_pairs = 0    !length of the per-pair buffers
  integer, save :: cap_terms = 0    !NumYHYTerms the constants and scratch were staged for
  logical, save :: consts_ready = .false.

  integer, parameter :: NN  = Glob_AllowedNumOfPseudoParticles
  integer, parameter :: NPT = NN*(NN+1)/2   !Glob_npt (real ECGs)
  integer, parameter :: CUF_BLK = 64        !threads/block: one thread per pair, 255 regs/thread
  !Upper bound on the per-pair scratch (term-indexed outputs) held on the device;
  !the chunk size is reduced below ECG_GPU_BATCH if needed to respect it.
  integer(8), parameter :: SCRATCH_BYTES = 536870912_8

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
  !  A failed CUDA call does not raise: it leaves garbage in a buffer that then
  !  flows into the MPI reduction and corrupts every rank. So on any error print
  !  the context and MPI_ABORT.
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
  !Called once, collectively, at startup, before ProgramDataInit allocates the
  !managed globals.
    character(64) :: val
    integer       :: ios, tmp, local_rank, allok
    integer(8)    :: tmp8
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
      call get_environment_variable('ECG_GPU_HEAP', val, status=ios)
      if (ios==0) then
        read(val,*,iostat=ios) tmp8
        if ((ios==0).and.(tmp8>0)) heap_bytes = tmp8
      endif
      write(*,'(1x,a)')             '---------------- GPU backend ----------------'
      write(*,'(1x,a,l1,a,i0,a,i0,a)')'  ME on GPU=',use_me, &
                                    '   batch=',batch_cap,'   heap=',heap_bytes/1048576,' MiB'
      if (.not.use_me) write(*,'(1x,a)') '  (GPU disabled: matrix elements run on CPU)'
      write(*,'(1x,a,l1)')          '  eig on GPU=',use_eig
      write(*,'(1x,a)')             '---------------------------------------------'
    endif
    call MPI_BCAST(use_eig,    1, MPI_LOGICAL,  0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_me,     1, MPI_LOGICAL,  0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(batch_cap,  1, MPI_INTEGER,  0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(heap_bytes, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
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
    gpu_eig_active = use_eig .and. ctx_ready
  end function gpu_eig_active

  subroutine stage_constants()
  !Upload the per-run invariants: the Y+Y matrices (re-laid out term-first, as
  !MatrixElementsOpt takes them) and their coefficients. Mass and charge are
  !managed module arrays. The term count is remembered because the symmetry
  !block can be rebuilt mid-run, which would make the staged copy (and the
  !term-indexed scratch) the wrong size.
    integer :: n, nterms, j
    real(wp), allocatable :: symT(:,:,:)
    n = Glob_n; nterms = Glob_NumYHYTerms
    if (consts_ready .and. (cap_terms == nterms)) return
    if (allocated(d_sym)) deallocate(d_sym, d_coeff)
    allocate(d_sym(nterms,n,n), d_coeff(nterms), symT(nterms,n,n))
    do j = 1, nterms
      symT(j,1:n,1:n) = Glob_YHYMatr(1:n,1:n,j)
    enddo
    call cuda_check(cudaMemcpy(d_sym,   symT,          nterms*n*n), 'stage_constants: YHY matrices')
    call cuda_check(cudaMemcpy(d_coeff, Glob_YHYCoeff, nterms),     'stage_constants: YHY coefficients')
    deallocate(symT)
    !term-indexed scratch is sized by nterms: force a regrow
    if (cap_terms /= nterms) cap_pairs = 0
    cap_terms = nterms; consts_ready = .true.
  end subroutine stage_constants

  subroutine stage_basis(n, npt_, Kmax)
  !Recompute the hoisted bra-side matrices of every function and upload them,
  !with the ket-side parameters and powers, ONCE per build. Buffers only grow.
    integer, intent(in) :: n, npt_, Kmax
    if (Kmax > cap_basis) then
      if (allocated(Lh))   deallocate(Lh, Ah, MAh)
      if (allocated(d_Lm)) deallocate(d_Lm, d_Am, d_MAm, d_par, d_pwr)
      allocate(Lh(n,n,Kmax), Ah(n,n,Kmax), MAh(n,n,Kmax))
      allocate(d_Lm(n,n,Kmax), d_Am(n,n,Kmax), d_MAm(n,n,Kmax), d_par(npt_,Kmax), d_pwr(Kmax))
      cap_basis = Kmax
    endif
    call PrecomputeMatrices_MA_L_A(npt_, Kmax, Glob_NonlinParam(1:npt_,1:Kmax), &
                                   Glob_MassMatrix(1:n,1:n), Lh, Ah, MAh)
    call cuda_check(cudaMemcpy(d_Lm,  Lh,  n*n*Kmax), 'stage_basis: L')
    call cuda_check(cudaMemcpy(d_Am,  Ah,  n*n*Kmax), 'stage_basis: A')
    call cuda_check(cudaMemcpy(d_MAm, MAh, n*n*Kmax), 'stage_basis: MA')
    call cuda_check(cudaMemcpy(d_par, Glob_NonlinParam(1:npt_,1:Kmax), npt_*Kmax), 'stage_basis: parameters')
    call cuda_check(cudaMemcpy(d_pwr, Glob_PWR(1:Kmax), Kmax), 'stage_basis: powers')
  end subroutine stage_basis

  integer function chunk_size(npairs)
  !ECG_GPU_BATCH, reduced if the term-indexed scratch of that many pairs would
  !exceed SCRATCH_BYTES. Identical on every rank (same inputs), so the
  !per-chunk ALLREDUCEs stay matched.
    integer, intent(in) :: npairs
    integer(8) :: per_pair
    per_pair = 8_8*Glob_NumYHYTerms*(2_8 + 4_8*Glob_npt)
    chunk_size = int(min(int(min(npairs, batch_cap),8), max(1_8, SCRATCH_BYTES/per_pair)))
  end function chunk_size

  subroutine ensure_pair_capacity(npairs)
  !Grow the per-pair buffers if this build needs more room than the last one.
    integer, intent(in) :: npairs
    integer :: nterms
    if (npairs <= cap_pairs) return
    nterms = Glob_NumYHYTerms
    if (allocated(d_k)) deallocate(d_k, d_l, d_gf, d_H, d_S, d_Dk, d_Dl, &
                                   d_bS, d_bH, d_bdSk, d_bdSl, d_bdHk, d_bdHl)
    allocate(d_k(npairs), d_l(npairs), d_gf(npairs), d_H(npairs), d_S(npairs))
    allocate(d_Dk(2*Glob_npt,npairs), d_Dl(2*Glob_npt,npairs))
    allocate(d_bS(nterms,npairs), d_bH(nterms,npairs))
    allocate(d_bdSk(nterms,Glob_npt,npairs), d_bdSl(nterms,Glob_npt,npairs), &
             d_bdHk(nterms,Glob_npt,npairs), d_bdHl(nterms,Glob_npt,npairs))
    call cuda_check(cudaGetLastError(), 'ensure_pair_capacity: device allocation')
    cap_pairs = npairs
  end subroutine ensure_pair_capacity

  subroutine gpu_build_HS(Nmin,Nmax,store)
  !Energy-only build over the pair triangle [Nmin,Nmax], in fixed-size chunks
  !so neither host nor device memory scales with the whole K(K+1)/2 pair list.
  !Pair ipair is evaluated by rank MOD(ipair,Glob_NumOfProcs), the others
  !return zeros, and the per-chunk ALLREDUCE assembles the sum -- the same
  !per-pair split ComputeMatElem uses. Each pair's H/S goes back through the
  !store callback (matform's StoreHS).
    integer, intent(in)    :: Nmin,Nmax
    procedure(store_hs_cb) :: store
    integer               :: k,l,npairs,ipair,nf,ip,batch,gbase
    integer, allocatable  :: k_list(:),l_list(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    call stage_constants()                      !once per run
    batch  = chunk_size(npairs)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (energy): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_basis(Glob_n, Glob_npt, Nmax)    !once per build, not per chunk
    call ensure_pair_capacity(batch)
    allocate(k_list(batch),l_list(batch),Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    ipair=0; nf=0; gbase=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((nf==batch).or.(ipair==npairs)) then
          call launch_energy(nf, gbase, Nmax, k_list, l_list, Hout, Sout)
          call MPI_ALLREDUCE(Hout,Hout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout,Sout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip))
          enddo
          gbase=ipair; nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,Hout,Sout,Hout_r,Sout_r)
  end subroutine gpu_build_HS

  subroutine gpu_build_HS_deriv(Nmin,Nmax,store)
  !Energy+gradient build, same chunking. gradflag per pair follows
  !ComputeMatElemAndDeriv: 2 (derivatives for k AND l) when l is itself being
  !optimized and differs from k, else 1. Each pair's H/S and gradient slabs go
  !back through the store callback (matform's StoreHSD), in its layout:
  !D(1:npt) = dH/dvechL, D(npt+1:2*npt) = dS/dvechL.
    integer, intent(in)     :: Nmin,Nmax
    procedure(store_hsd_cb) :: store
    integer               :: npt2,k,l,npairs,ipair,nf,ip,batch,gbase
    integer, allocatable  :: k_list(:),l_list(:),gf_list(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    real(wp),allocatable  :: Dkout(:,:),Dlout(:,:),Dkout_r(:,:),Dlout_r(:,:)
    npt2 = 2*Glob_npt
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    call stage_constants()                      !once per run
    batch  = chunk_size(npairs)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (grad): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_basis(Glob_n, Glob_npt, Nmax)    !once per build, not per chunk
    call ensure_pair_capacity(batch)
    allocate(k_list(batch),l_list(batch),gf_list(batch))
    allocate(Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    allocate(Dkout(npt2,batch),Dlout(npt2,batch),Dkout_r(npt2,batch),Dlout_r(npt2,batch))
    Dlout_r = ZERO
    ipair=0; nf=0; gbase=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((l>Glob_nfru).and.(l/=k)) then
          gf_list(nf)=2
        else
          gf_list(nf)=1
        endif
        if ((nf==batch).or.(ipair==npairs)) then
          call launch_grad(nf, gbase, Nmax, k_list, l_list, gf_list, Hout, Sout, Dkout, Dlout)
          call MPI_ALLREDUCE(Hout, Hout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout, Sout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Dkout,Dkout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) &
          call MPI_ALLREDUCE(Dlout,Dlout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip), &
                       Dkout_r(1:npt2,ip),Dlout_r(1:npt2,ip))
          enddo
          gbase=ipair; nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,gf_list)
    deallocate(Hout,Sout,Hout_r,Sout_r)
    deallocate(Dkout,Dlout,Dkout_r,Dlout_r)
  end subroutine gpu_build_HS_deriv

  ! ==========================================================================
  !  (2) DEVICE CODE -- kernels and host launchers
  ! ==========================================================================

  !One thread per pair. Thread i owns pair i of the chunk, whose global index
  !is gbase+i; only the rank MOD(gbase+i,nprocs) evaluates it, the others write
  !zeros (summed away by the ALLREDUCE). The pair's terms are summed in order,
  !so the result does not depend on scheduling.
  attributes(global) subroutine me_energy_kernel(npairs, Kmax, nterms, gbase, nprocs, procid, &
      k_list, l_list, Lm, Am, MAm, par, pwr, sym, coeff, bS, bH, bdSk, bdSl, bdHk, bdHl, &
      Hout, Sout)
    integer, value  :: npairs, Kmax, nterms, gbase, nprocs, procid
    integer         :: k_list(npairs), l_list(npairs), pwr(Kmax)
    real(wp)        :: Lm(NN,NN,Kmax), Am(NN,NN,Kmax), MAm(NN,NN,Kmax), par(NPT,Kmax)
    real(wp)        :: sym(nterms,NN,NN), coeff(nterms)
    real(wp)        :: bS(nterms,npairs), bH(nterms,npairs)
    real(wp)        :: bdSk(nterms,NPT,npairs), bdSl(nterms,NPT,npairs)
    real(wp)        :: bdHk(nterms,NPT,npairs), bdHl(nterms,NPT,npairs)
    real(wp)        :: Hout(npairs), Sout(npairs)
    integer  :: i, k, l, t
    real(wp) :: hsum, ssum

    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > npairs) return
    hsum = 0.0_wp; ssum = 0.0_wp
    if (mod(gbase+i, nprocs) == procid) then
      k = k_list(i); l = l_list(i)
      call MatrixElementsOpt(pwr(k), Lm(:,:,k), Am(:,:,k), MAm(:,:,k), pwr(l), par(:,l), sym, &
                             bS(:,i), bH(:,i), bdSk(:,:,i), bdSl(:,:,i), bdHk(:,:,i), bdHl(:,:,i), 0)
      do t = 1, nterms
        hsum = hsum + coeff(t)*bH(t,i)
        ssum = ssum + coeff(t)*bS(t,i)
      enddo
    endif
    Hout(i) = hsum
    Sout(i) = ssum
  end subroutine me_energy_kernel

  !Host launcher for the energy build. The basis and the constants are already
  !resident, so this only ships the chunk's pair list, launches, and collects
  !H/S. No device allocation happens here.
  subroutine launch_energy(npairs, gbase, Kmax, k_list, l_list, Hout, Sout)
    integer,  intent(in) :: npairs, gbase, Kmax
    integer,  intent(in) :: k_list(npairs), l_list(npairs)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    integer :: nblk
    call cuda_check(cudaMemcpy(d_k, k_list, npairs), 'launch_energy: H2D k list')
    call cuda_check(cudaMemcpy(d_l, l_list, npairs), 'launch_energy: H2D l list')
    nblk = (npairs + CUF_BLK - 1)/CUF_BLK
    call me_energy_kernel<<<nblk, CUF_BLK>>>(npairs, Kmax, cap_terms, gbase, Glob_NumOfProcs, Glob_ProcID, &
        d_k, d_l, d_Lm, d_Am, d_MAm, d_par, d_pwr, d_sym, d_coeff, &
        d_bS, d_bH, d_bdSk, d_bdSl, d_bdHk, d_bdHl, d_H, d_S)
    call cuda_check(cudaGetLastError(),      'launch_energy: kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'launch_energy: kernel execution')
    call cuda_check(cudaMemcpy(Hout, d_H, npairs), 'launch_energy: D2H H')
    call cuda_check(cudaMemcpy(Sout, d_S, npairs), 'launch_energy: D2H S')
  end subroutine launch_energy

  !Energy+gradient kernel: as the energy kernel, plus the coefficient-weighted
  !term sums of the derivative buffers into this pair's Dk/Dl slabs.
  attributes(global) subroutine me_grad_kernel(npairs, Kmax, nterms, gbase, nprocs, procid, &
      k_list, l_list, gf_list, Lm, Am, MAm, par, pwr, sym, coeff, bS, bH, bdSk, bdSl, bdHk, bdHl, &
      Hout, Sout, Dkout, Dlout)
    integer, value  :: npairs, Kmax, nterms, gbase, nprocs, procid
    integer         :: k_list(npairs), l_list(npairs), gf_list(npairs), pwr(Kmax)
    real(wp)        :: Lm(NN,NN,Kmax), Am(NN,NN,Kmax), MAm(NN,NN,Kmax), par(NPT,Kmax)
    real(wp)        :: sym(nterms,NN,NN), coeff(nterms)
    real(wp)        :: bS(nterms,npairs), bH(nterms,npairs)
    real(wp)        :: bdSk(nterms,NPT,npairs), bdSl(nterms,NPT,npairs)
    real(wp)        :: bdHk(nterms,NPT,npairs), bdHl(nterms,NPT,npairs)
    real(wp)        :: Hout(npairs), Sout(npairs), Dkout(2*NPT,npairs), Dlout(2*NPT,npairs)
    integer  :: i, k, l, t, p, gf
    real(wp) :: hsum, ssum, dhk, dsk, dhl, dsl

    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > npairs) return
    hsum = 0.0_wp; ssum = 0.0_wp
    do p = 1, 2*NPT
      Dkout(p,i) = 0.0_wp
      Dlout(p,i) = 0.0_wp
    enddo
    if (mod(gbase+i, nprocs) == procid) then
      k = k_list(i); l = l_list(i); gf = gf_list(i)
      call MatrixElementsOpt(pwr(k), Lm(:,:,k), Am(:,:,k), MAm(:,:,k), pwr(l), par(:,l), sym, &
                             bS(:,i), bH(:,i), bdSk(:,:,i), bdSl(:,:,i), bdHk(:,:,i), bdHl(:,:,i), gf)
      do t = 1, nterms
        hsum = hsum + coeff(t)*bH(t,i)
        ssum = ssum + coeff(t)*bS(t,i)
      enddo
      do p = 1, NPT
        dhk = 0.0_wp; dsk = 0.0_wp; dhl = 0.0_wp; dsl = 0.0_wp
        do t = 1, nterms
          dhk = dhk + coeff(t)*bdHk(t,p,i)
          dsk = dsk + coeff(t)*bdSk(t,p,i)
          if (gf == 2) then
            dhl = dhl + coeff(t)*bdHl(t,p,i)
            dsl = dsl + coeff(t)*bdSl(t,p,i)
          endif
        enddo
        Dkout(p,i)     = dhk
        Dkout(NPT+p,i) = dsk
        Dlout(p,i)     = dhl
        Dlout(NPT+p,i) = dsl
      enddo
    endif
    Hout(i) = hsum
    Sout(i) = ssum
  end subroutine me_grad_kernel

  !Host launcher for the gradient build.
  subroutine launch_grad(npairs, gbase, Kmax, k_list, l_list, gf_list, Hout, Sout, Dkout, Dlout)
    integer,  intent(in) :: npairs, gbase, Kmax
    integer,  intent(in) :: k_list(npairs), l_list(npairs), gf_list(npairs)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp), intent(out):: Dkout(2*Glob_npt,npairs), Dlout(2*Glob_npt,npairs)
    integer :: nblk
    call cuda_check(cudaMemcpy(d_k,  k_list,  npairs), 'launch_grad: H2D k list')
    call cuda_check(cudaMemcpy(d_l,  l_list,  npairs), 'launch_grad: H2D l list')
    call cuda_check(cudaMemcpy(d_gf, gf_list, npairs), 'launch_grad: H2D gradflags')
    nblk = (npairs + CUF_BLK - 1)/CUF_BLK
    call me_grad_kernel<<<nblk, CUF_BLK>>>(npairs, Kmax, cap_terms, gbase, Glob_NumOfProcs, Glob_ProcID, &
        d_k, d_l, d_gf, d_Lm, d_Am, d_MAm, d_par, d_pwr, d_sym, d_coeff, &
        d_bS, d_bH, d_bdSk, d_bdSl, d_bdHk, d_bdHl, d_H, d_S, d_Dk, d_Dl)
    call cuda_check(cudaGetLastError(),      'launch_grad: kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'launch_grad: kernel execution')
    call cuda_check(cudaMemcpy(Hout,  d_H,  npairs),            'launch_grad: D2H H')
    call cuda_check(cudaMemcpy(Sout,  d_S,  npairs),            'launch_grad: D2H S')
    call cuda_check(cudaMemcpy(Dkout, d_Dk, 2*Glob_npt*npairs), 'launch_grad: D2H Dk')
    call cuda_check(cudaMemcpy(Dlout, d_Dl, 2*Glob_npt*npairs), 'launch_grad: D2H Dl')
  end subroutine launch_grad

  !Device lifecycle. gpu_init selects device = node-local-rank % nGPUs, sets the
  !device heap (needed by MatrixElementsOpt's run-time-sized arrays, see the
  !header), and returns .false. (does NOT abort) on any failure, so the caller
  !decides collectively.
  logical function gpu_init(local_rank)
    integer, intent(in)  :: local_rank
    integer :: ndev, dev, istat
    integer(kind=cuda_count_kind) :: heap
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
    !The heap limit only takes effect before the first kernel launch
    istat = cudaDeviceSetLimit(cudaLimitMallocHeapSize, heap_bytes)
    if (istat /= cudaSuccess) then
      write(*,'(1x,a,i0,a,a)') 'gpu_init: rank ',Glob_ProcID,' cannot set the device heap: ', &
        trim(cudaGetErrorString(istat))
      return
    endif
    istat = cudaDeviceGetLimit(heap, cudaLimitMallocHeapSize) ; if (istat /= cudaSuccess) return
    write(*,'(1x,a,i0,a,i0,a,i0,a,a,a,i0,a)') 'gpu_init: rank ',Glob_ProcID, &
        ' (node-local ',local_rank,') -> device ',dev,' ',trim(prop%name),', heap ',heap/1048576,' MiB'
    gpu_init = .true.
  end function gpu_init

  subroutine gpu_finalize()
    integer :: istat
  !Release application-owned device and host buffers at shutdown.
  !Called once at shutdown (main, before MPI_FINALIZE). No cudaDeviceReset:
  !NVHPC performs runtime cleanup later and fails with CUDA error 709 otherwise.
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
    if (allocated(d_Lm))  deallocate(d_Lm, d_Am, d_MAm, d_par, d_pwr)
    if (allocated(d_sym)) deallocate(d_sym, d_coeff)
    if (allocated(d_k))   deallocate(d_k, d_l, d_gf, d_H, d_S, d_Dk, d_Dl, &
                                     d_bS, d_bH, d_bdSk, d_bdSl, d_bdHk, d_bdHl)
    if (allocated(Lh))    deallocate(Lh, Ah, MAh)
    cap_basis=0; cap_pairs=0; cap_terms=0; consts_ready=.false.
  end subroutine gpu_finalize

  !cuSOLVER replacement for the DSYGVX calls of the G method: computes the
  !iwhich-th eigenvalue (and, if jobz_vec/=0, its eigenvector) of H c = E S c.
  !Unchanged from the RG backend.
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

end module gpu_backend
