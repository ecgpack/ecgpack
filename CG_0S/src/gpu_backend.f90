module gpu_backend
!CUDA backend for CG_0S matrix-element assembly. Compiled only with -DUSE_CUDA.
!
!This is the complex-ECG counterpart of RG_0S/src/gpu_backend.f90 and keeps its
!structure and public API. The differences all follow from CG_0S being complex:
!  * the per-function hoisted inputs are L (real) and C=L*L'+i*B (complex), built
!    by PrecomputeMatrices_L_C instead of PrecomputeMatrices_L_A;
!  * H, S and the gradient slabs Dk, Dl are complex. CUDA Fortran has no complex
!    atomicAdd, so every complex accumulator is held as separate real and
!    imaginary shared slots. Complex addition is componentwise, so the
!    deterministic path (ECG_DETERM=1) still sums in exactly the CPU order;
!  * the gradient slabs are 2*Glob_npt = 4*Glob_np complex numbers per function;
!  * the generalized eigenproblem is complex Hermitian: cusolverDnZhegvdx replaces
!    cusolverDnDsygvdx (gpu_zhegvx replaces gpu_dsygvx).
!
!The kernels call the original host/device matrix-element routine in matelem.f90.
!The original matrix-element routine reads managed mass and charge module arrays.
!Basis and permutation arrays are staged explicitly on the device.
  use globvars        !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem,   only: MatrixElementsHS_CG_0S, PrecomputeMatrices_L_C
  use cudafor
  use cublas
  use cusolverDn
  implicit none
  private

  !Public surface used by the shared sources (every call site is #ifdef USE_CUDA):
  public :: gpu_backend_init      !collective one-time startup            (main)
  public :: gpu_finalize          !release device/handle at shutdown      (main)
  public :: gpu_eig_active
  public :: gpu_zhegvx
  public :: gpu_active            !.true. => use the CUDA ME backend      (matform)
  public :: gpu_build_HS          !energy matrix-element build            (matform)
  public :: gpu_build_HS_deriv    !energy+gradient matrix-element build   (matform)

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call
  logical, save :: ctx_ready = .false.   !device context successfully initialised
  logical, save :: use_determ = .false.  !ECG_DETERM=1: deterministic term reduction
                                         !(ordered, CPU-matching) instead of atomicAdd

  type(cusolverDnHandle), save :: eig_hdl
  logical, save :: eig_created = .false.
  integer, save :: eig_cap_n = 0, eig_cap_work = 0
  complex(wp), device, allocatable, save :: dA(:,:), dB(:,:), dwork(:)
  real(wp),    device, allocatable, save :: dW(:)
  integer,     device, allocatable, save :: dinfo(:)
  complex(wp), allocatable, save :: hA(:,:), hB(:,:)

  !PERSISTENT MATRIX-ELEMENT WORKSPACE (see RG_0S/src/gpu_backend.f90 for the
  !measurements behind it): allocated once, grown on demand, released in
  !gpu_finalize. The constants are uploaded once; the basis once per build.
  real(wp),    device, allocatable, save :: d_Lm(:,:,:)
  complex(wp), device, allocatable, save :: d_Cm(:,:,:)
  real(wp),    device, allocatable, save :: d_YHY(:), d_coeff(:)
  integer,     device, allocatable, save :: d_k(:), d_l(:), d_gl(:)
  complex(wp), device, allocatable, save :: d_H(:), d_S(:), d_Dk(:), d_Dl(:)
  real(wp),    allocatable, save         :: Lh(:,:,:)   !host workspace
  complex(wp), allocatable, save         :: Ch(:,:,:)   !host workspace
  integer, save :: cap_basis = 0    !3rd extent of the staged basis arrays
  integer, save :: cap_pairs = 0    !length of the per-pair buffers
  integer, save :: cap_grad  = 0    !length of the gradient slabs
  integer, save :: cap_terms = 0    !NumYHYTerms the constants were staged for
  logical, save :: consts_ready = .false.

  integer, parameter :: NN  = Glob_AllowedNumOfPseudoParticles
  integer, parameter :: NNP = NN*(NN+1)/2   !np=n(n+1)/2; Dk/Dl hold 4*np complex numbers
  integer, parameter :: NG_MAX = 4*NNP      !complex gradient entries per function (2*Glob_npt);
                                            !distinct from the run-time ng (case-insensitive!)
  integer, parameter :: CUF_BLK = 128       !threads/block. The term loop STRIDES (see the kernels), so
                                            !block size is independent of NumYHYTerms.

  !Signatures of matform's StoreHS / StoreHSD, matched exactly so those routines
  !can be handed in as callbacks. The chunk loops call them one pair at a time.
  abstract interface
    subroutine store_hs_cb(i,j,Hij,Sij)
      import :: wp
      integer,     intent(in) :: i,j
      complex(wp), intent(in) :: Hij,Sij
    end subroutine store_hs_cb
    subroutine store_hsd_cb(i,j,Hij,Sij,Di,Dj)
      import :: wp, Glob_npt
      integer,     intent(in) :: i,j
      complex(wp), intent(in) :: Hij,Sij
      complex(wp), intent(in) :: Di(2*Glob_npt),Dj(2*Glob_npt)
    end subroutine store_hsd_cb
  end interface

contains

  ! ==========================================================================
  !  (0) FAIL-CLOSED STATUS CHECKING
  !  On the GPU a failure does not raise -- it silently leaves garbage in a
  !  buffer that then flows into the MPI reduction and corrupts every rank. So
  !  on any error we print context and MPI_ABORT.
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
      write(*,'(1x,a)')             '---------------- GPU backend ----------------'
      write(*,'(1x,a,l1,a,i0,a,l1)')'  ME on GPU=',use_me, &
                                    '   batch=',batch_cap,'   determ=',use_determ
      if (.not.use_me) write(*,'(1x,a)') '  (GPU disabled: matrix elements run on CPU)'
      write(*,'(1x,a,l1)')          '  eig on GPU=',use_eig
      write(*,'(1x,a)')             '---------------------------------------------'
    endif
    call MPI_BCAST(use_eig,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_me,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(batch_cap,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_determ, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    if (.not.use_me) return

    !Bring up the device context on every rank, then agree COLLECTIVELY: if ANY
    !rank failed to acquire a GPU, abort everywhere.
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
  !Upload the per-run invariants: the Y+Y permutation matrices and their
  !coefficients (both real in CG_0S as well). Mass and scaled charge use shared
  !managed module arrays. The term count is remembered because the symmetry
  !block can be rebuilt mid-run, which would make the staged copy the wrong size.
    integer :: n, nterms
    n = Glob_n; nterms = Glob_NumYHYTerms
    if (consts_ready .and. (cap_terms == nterms)) return
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff)
    allocate(d_YHY(n*n*nterms), d_coeff(nterms))
    call cuda_check(cudaMemcpy(d_YHY,   Glob_YHYMatr,  n*n*nterms), 'stage_constants: YHY matrices')
    call cuda_check(cudaMemcpy(d_coeff, Glob_YHYCoeff, nterms),     'stage_constants: YHY coefficients')
    cap_terms = nterms; consts_ready = .true.
  end subroutine stage_constants

  subroutine stage_basis(n, npt, Kmax)
  !Recompute the hoisted per-function matrices L and C=L*L'+i*B and upload them
  !ONCE per build. The buffers only ever grow.
    integer, intent(in) :: n, npt, Kmax
    if (Kmax > cap_basis) then
      if (allocated(Lh))   deallocate(Lh, Ch)
      if (allocated(d_Lm)) deallocate(d_Lm, d_Cm)
      allocate(Lh(n,n,Kmax), Ch(n,n,Kmax))
      allocate(d_Lm(n,n,Kmax), d_Cm(n,n,Kmax))
      cap_basis = Kmax
    endif
    call PrecomputeMatrices_L_C(npt,Kmax,Glob_NonlinParam(1:npt,1:Kmax),Lh,Ch)
    call cuda_check(cudaMemcpy(d_Lm, Lh, n*n*Kmax), 'stage_basis: L')
    call cuda_check(cudaMemcpy(d_Cm, Ch, n*n*Kmax), 'stage_basis: C')
  end subroutine stage_basis

  subroutine ensure_pair_capacity(npairs, ng, need_grad)
  !Grow the per-pair buffers if this build needs more room than the last one.
  !ng is the number of complex gradient entries per pair (2*Glob_npt).
    integer, intent(in) :: npairs, ng
    logical, intent(in) :: need_grad
    if (npairs > cap_pairs) then
      if (allocated(d_k)) deallocate(d_k, d_l, d_gl, d_H, d_S)
      allocate(d_k(npairs), d_l(npairs), d_gl(npairs), d_H(npairs), d_S(npairs))
      cap_pairs = npairs
    endif
    if (need_grad) then
      if (npairs*ng > cap_grad) then
        if (allocated(d_Dk)) deallocate(d_Dk, d_Dl)
        allocate(d_Dk(npairs*ng), d_Dl(npairs*ng))
        cap_grad = npairs*ng
      endif
    endif
  end subroutine ensure_pair_capacity

  subroutine gpu_build_HS(Nmin,Nmax,store)
  !Energy-only matrix-element build over the basis pair triangle [Nmin,Nmax],
  !processed in fixed-size chunks so neither host nor device memory scales with
  !the whole K(K+1)/2 pair list. Each pair's summed H/S is handed back via the
  !store callback (matform's StoreHS). Every rank uses the same batch_cap, so the
  !per-chunk ALLREDUCE stays collective. As in matform, the complex buffers are
  !reduced as real arrays of twice the length.
    integer, intent(in)    :: Nmin,Nmax
    procedure(store_hs_cb) :: store
    integer                  :: n,npt,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable     :: k_list(:),l_list(:)
    complex(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    n=Glob_n; npt=Glob_npt
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (energy): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_constants()                      !once per run
    call stage_basis(n, npt, Nmax)              !once per build, not per chunk
    call ensure_pair_capacity(batch, 0, .false.)
    allocate(k_list(batch),l_list(batch),Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    ipair=0; nf=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((nf==batch).or.(ipair==npairs)) then
          call cuf_compute_matelem_batch(Nmax, k_list, l_list, nf, &
              Glob_NumYHYTerms, n, Glob_NumOfProcs, Glob_ProcID, &
              Hout, Sout)
          call MPI_ALLREDUCE(Hout,Hout_r,2*nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout,Sout_r,2*nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
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
  !Energy+gradient build. The (2*npt x pairs) complex gradient buffers are the
  !memory wall at large K, so the same chunking applies. Each pair's H/S and
  !gradient slabs are handed back via the store callback (matform's StoreHSD).
    integer, intent(in)     :: Nmin,Nmax
    procedure(store_hsd_cb) :: store
    integer                  :: n,npt,ng,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable     :: k_list(:),l_list(:),grad_l(:)
    complex(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    complex(wp),allocatable  :: Dkout(:,:),Dlout(:,:),Dkout_r(:,:),Dlout_r(:,:)
    n=Glob_n; npt=Glob_npt; ng=2*npt
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (grad): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    call stage_constants()                      !once per run
    call stage_basis(n, npt, Nmax)              !once per build, not per chunk
    call ensure_pair_capacity(batch, ng, .true.)
    allocate(k_list(batch),l_list(batch),grad_l(batch))
    allocate(Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    allocate(Dkout(ng,batch),Dlout(ng,batch),Dkout_r(ng,batch),Dlout_r(ng,batch))
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
          Dkout(1:ng,1:nf)=ZERO
          Dlout(1:ng,1:nf)=ZERO
          call cuf_compute_matelem_deriv_batch(Nmax, k_list, l_list, grad_l, nf, &
              Glob_NumYHYTerms, n, ng, Glob_NumOfProcs, Glob_ProcID, &
              Hout, Sout, Dkout, Dlout)
          call MPI_ALLREDUCE(Hout, Hout_r, 2*nf,    MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout, Sout_r, 2*nf,    MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Dkout,Dkout_r,2*nf*ng, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) &
          call MPI_ALLREDUCE(Dlout,Dlout_r,2*nf*ng, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip), &
                       Dkout_r(1:ng,ip),Dlout_r(1:ng,ip))
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
  !terms. Each thread calls the shared MatrixElementsHS_CG_0S (no gradients).
  attributes(global) subroutine me_energy_kernel(Lm, Cm, Kmax, &
      k_list, l_list, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, &
      determ, Hout, Sout)
    integer, value   :: Kmax, nterms, n, nprocs, procid, determ
    real(wp)         :: Lm(n,n,Kmax)
    complex(wp)      :: Cm(n,n,Kmax)
    integer          :: k_list(*), l_list(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    complex(wp)      :: Hout(*), Sout(*)
    real(wp), shared :: sh_Hr, sh_Hi, sh_Sr, sh_Si
    real(wp), shared :: shm(*)                !dynamic: 4*nterms, deterministic term buffer
    integer     :: pair_idx, j, qbase, k0, l0, istat
    real(wp)    :: coeff, hr, hi, sr, si
    complex(wp) :: Hkl, Skl
    complex(wp) :: dDk(NG_MAX), dDl(NG_MAX)           !unused energy-path gradient outputs

    pair_idx = blockIdx%x
    qbase = (pair_idx-1)*nterms - 1          !matches Fortran q=(i-1)*nterms-1
    k0    = k_list(pair_idx)                 !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)

    if (determ /= 0) then
      !DETERMINISTIC: each term's contribution goes to shm (non-owned -> 0.0,
      !exact), then thread 1 sums them IN ORDER, matching the CPU term loop.
      !shm layout: [Re H | Im H | Re S | Im S], nterms each.
      do j = threadIdx%x, nterms, blockDim%x
        if (mod(qbase+j, nprocs) == procid) then
          call MatrixElementsHS_CG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                              Cm(1,1,k0), Cm(1,1,l0), &
                              YHYMatr((j-1)*n*n+1), &
                              Hkl, Skl, dDk, dDl, .false., .false.)
          coeff = YHYCoeff(j)
          shm(j)          = coeff*real(Hkl,wp)
          shm(nterms+j)   = coeff*aimag(Hkl)
          shm(2*nterms+j) = coeff*real(Skl,wp)
          shm(3*nterms+j) = coeff*aimag(Skl)
        else
          shm(j)          = 0.0_wp
          shm(nterms+j)   = 0.0_wp
          shm(2*nterms+j) = 0.0_wp
          shm(3*nterms+j) = 0.0_wp
        endif
      enddo
      call syncthreads()
      if (threadIdx%x==1) then
        hr=0.0_wp; hi=0.0_wp; sr=0.0_wp; si=0.0_wp
        do j=1,nterms
          hr = hr + shm(j)
          hi = hi + shm(nterms+j)
          sr = sr + shm(2*nterms+j)
          si = si + shm(3*nterms+j)
        enddo
        Hout(pair_idx)=cmplx(hr,hi,wp); Sout(pair_idx)=cmplx(sr,si,wp)
      endif
    else
      !FAST (default): atomicAdd reduction -- non-deterministic term order.
      if (threadIdx%x==1) then
        sh_Hr=0.0_wp; sh_Hi=0.0_wp; sh_Sr=0.0_wp; sh_Si=0.0_wp
      endif
      call syncthreads()
      do j = threadIdx%x, nterms, blockDim%x   !STRIDED: block size independent of nterms
        if (mod(qbase+j, nprocs) == procid) then
          call MatrixElementsHS_CG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                              Cm(1,1,k0), Cm(1,1,l0), &
                              YHYMatr((j-1)*n*n+1), &
                              Hkl, Skl, dDk, dDl, .false., .false.)
          coeff = YHYCoeff(j)
          istat = atomicadd(sh_Hr, coeff*real(Hkl,wp))
          istat = atomicadd(sh_Hi, coeff*aimag(Hkl))
          istat = atomicadd(sh_Sr, coeff*real(Skl,wp))
          istat = atomicadd(sh_Si, coeff*aimag(Skl))
        endif
      enddo
      call syncthreads()
      if (threadIdx%x==1) then
        Hout(pair_idx)=cmplx(sh_Hr,sh_Hi,wp); Sout(pair_idx)=cmplx(sh_Sr,sh_Si,wp)
      endif
    endif
  end subroutine me_energy_kernel

  !Host launcher for the energy build. The basis and the constants are already
  !resident (stage_basis / stage_constants), so this only ships the chunk's pair
  !list, launches, and collects H/S. No device allocation happens here.
  subroutine cuf_compute_matelem_batch(Kmax, k_list, l_list, npairs, &
      nterms, n, nprocs, procid, Hout, Sout)
    integer,     intent(in) :: Kmax, npairs, nterms, n, nprocs, procid
    integer,     intent(in) :: k_list(npairs), l_list(npairs)
    complex(wp), intent(out):: Hout(npairs), Sout(npairs)
    integer  :: blk, determ, shmem

    call cuda_check(cudaMemcpy(d_k, k_list, npairs), 'cuf_compute_matelem_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l, l_list, npairs), 'cuf_compute_matelem_batch: H2D l list')

    blk = min(nterms, CUF_BLK)
    determ = 0; shmem = 0
    if (use_determ) then; determ = 1; shmem = 4*nterms*8; endif   !4*nterms real*8 dynamic shared
    call me_energy_kernel<<<npairs, blk, shmem>>>(d_Lm, d_Cm, Kmax, d_k, d_l, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        determ, d_H, d_S)
    call cuda_check(cudaGetLastError(),      'cuf_compute_matelem_batch: energy kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'cuf_compute_matelem_batch: energy kernel execution')

    call cuda_check(cudaMemcpy(Hout, d_H, npairs), 'cuf_compute_matelem_batch: D2H H')
    call cuda_check(cudaMemcpy(Sout, d_S, npairs), 'cuf_compute_matelem_batch: D2H S')
  end subroutine cuf_compute_matelem_batch

  !Energy+gradient kernel: one block per (k,l) pair; threads STRIDE over the terms.
  !Each thread calls the shared MatrixElementsHS_CG_0S (grad_k=true, grad_l=flag)
  !and atomic-accumulates H, S and the gradient slabs Dk, Dl into shared memory,
  !real and imaginary parts in separate arrays (sh_Dkr=Re Dk, sh_Dki=Im Dk).
  !All of them are sized by the compile-time NG_MAX. They must not be sized by the
  !run-time argument ng: that makes them automatic arrays, which nvfortran puts in
  !the DYNAMIC shared segment that this launch does not size (illegal address on
  !V100/H100, silently tolerated on sm_86). Fortran is case-insensitive, so the
  !parameter must not be named NG either, or the dummy ng shadows it.
  attributes(global) subroutine me_grad_kernel(Lm, Cm, Kmax, &
      k_list, l_list, grad_l_flag, YHYMatr, YHYCoeff, nterms, n, ng, nprocs, procid, &
      Hout, Sout, Dkout, Dlout)
    integer, value   :: Kmax, nterms, n, ng, nprocs, procid
    real(wp)         :: Lm(n,n,Kmax)
    complex(wp)      :: Cm(n,n,Kmax)
    integer          :: k_list(*), l_list(*), grad_l_flag(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    complex(wp)      :: Hout(*), Sout(*), Dkout(*), Dlout(*)   !Dk/Dl slabs (ng,npairs) flattened
    real(wp), shared :: sh_Hr, sh_Hi, sh_Sr, sh_Si
    real(wp), shared :: sh_Dkr(NG_MAX), sh_Dki(NG_MAX), sh_Dlr(NG_MAX), sh_Dli(NG_MAX)
    integer     :: pair_idx, j, qbase, k0, l0, istat, q, nthr
    real(wp)    :: coeff
    complex(wp) :: Hkl, Skl, Dk(NG_MAX), Dl(NG_MAX)
    logical     :: gl

    pair_idx = blockIdx%x
    nthr     = blockDim%x
    if (threadIdx%x==1) then
      sh_Hr=0.0_wp; sh_Hi=0.0_wp; sh_Sr=0.0_wp; sh_Si=0.0_wp
    endif
    do q=threadIdx%x,ng,nthr      !stride-init the shared gradient slabs
      sh_Dkr(q)=0.0_wp; sh_Dki(q)=0.0_wp
      sh_Dlr(q)=0.0_wp; sh_Dli(q)=0.0_wp
    enddo
    call syncthreads()

    qbase = (pair_idx-1)*nterms - 1
    k0    = k_list(pair_idx)                !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)
    gl    = (grad_l_flag(pair_idx)==1)
    do j = threadIdx%x, nterms, nthr        !STRIDED: block size independent of nterms
      if (mod(qbase+j, nprocs) == procid) then
        call MatrixElementsHS_CG_0S(Lm(1,1,k0), Lm(1,1,l0), &
                            Cm(1,1,k0), Cm(1,1,l0), &
                            YHYMatr((j-1)*n*n+1), &
                            Hkl, Skl, Dk, Dl, .true., gl)
        coeff = YHYCoeff(j)
        istat = atomicadd(sh_Hr, coeff*real(Hkl,wp))
        istat = atomicadd(sh_Hi, coeff*aimag(Hkl))
        istat = atomicadd(sh_Sr, coeff*real(Skl,wp))
        istat = atomicadd(sh_Si, coeff*aimag(Skl))
        do q=1,ng
          istat = atomicadd(sh_Dkr(q), coeff*real(Dk(q),wp))
          istat = atomicadd(sh_Dki(q), coeff*aimag(Dk(q)))
          if (gl) then
            istat = atomicadd(sh_Dlr(q), coeff*real(Dl(q),wp))
            istat = atomicadd(sh_Dli(q), coeff*aimag(Dl(q)))
          endif
        enddo
      endif
    enddo
    call syncthreads()

    if (threadIdx%x==1) then
      Hout(pair_idx)=cmplx(sh_Hr,sh_Hi,wp); Sout(pair_idx)=cmplx(sh_Sr,sh_Si,wp)
    endif
    do q=threadIdx%x,ng,nthr      !stride-write the slabs back to global
      Dkout((pair_idx-1)*ng+q)=cmplx(sh_Dkr(q),sh_Dki(q),wp)
      Dlout((pair_idx-1)*ng+q)=cmplx(sh_Dlr(q),sh_Dli(q),wp)
    enddo
  end subroutine me_grad_kernel

  !Host launcher for the gradient build. Same as above: resident basis and
  !constants, so only the chunk's pair list moves in and the results move out.
  subroutine cuf_compute_matelem_deriv_batch(Kmax, k_list, l_list, grad_l_flag, &
      npairs, nterms, n, ng, nprocs, procid, Hout, Sout, Dkout, Dlout)
    integer,     intent(in) :: Kmax, npairs, nterms, n, ng, nprocs, procid
    integer,     intent(in) :: k_list(npairs), l_list(npairs), grad_l_flag(npairs)
    complex(wp), intent(out):: Hout(npairs), Sout(npairs)
    complex(wp), intent(out):: Dkout(ng*npairs), Dlout(ng*npairs)
    integer  :: blk

    call cuda_check(cudaMemcpy(d_k,  k_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D k list')
    call cuda_check(cudaMemcpy(d_l,  l_list,      npairs), 'cuf_compute_matelem_deriv_batch: H2D l list')
    call cuda_check(cudaMemcpy(d_gl, grad_l_flag, npairs), 'cuf_compute_matelem_deriv_batch: H2D grad flags')

    blk = min(nterms, CUF_BLK)
    call me_grad_kernel<<<npairs, blk>>>(d_Lm, d_Cm, Kmax, d_k, d_l, d_gl, &
        d_YHY, d_coeff, nterms, n, ng, nprocs, procid, &
        d_H, d_S, d_Dk, d_Dl)
    call cuda_check(cudaGetLastError(),      'cuf_compute_matelem_deriv_batch: grad kernel launch')
    call cuda_check(cudaDeviceSynchronize(), 'cuf_compute_matelem_deriv_batch: grad kernel execution')

    call cuda_check(cudaMemcpy(Hout,  d_H,  npairs),    'cuf_compute_matelem_deriv_batch: D2H H')
    call cuda_check(cudaMemcpy(Sout,  d_S,  npairs),    'cuf_compute_matelem_deriv_batch: D2H S')
    call cuda_check(cudaMemcpy(Dkout, d_Dk, ng*npairs), 'cuf_compute_matelem_deriv_batch: D2H Dk')
    call cuda_check(cudaMemcpy(Dlout, d_Dl, ng*npairs), 'cuf_compute_matelem_deriv_batch: D2H Dl')
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
    if (allocated(d_Lm))  deallocate(d_Lm, d_Cm)
    if (allocated(d_YHY)) deallocate(d_YHY, d_coeff)
    if (allocated(d_k))   deallocate(d_k, d_l, d_gl, d_H, d_S)
    if (allocated(d_Dk))  deallocate(d_Dk, d_Dl)
    if (allocated(Lh))    deallocate(Lh, Ch)
    cap_basis=0; cap_pairs=0; cap_grad=0; cap_terms=0; consts_ready=.false.
  end subroutine gpu_finalize

  subroutine gpu_zhegvx(jobz_vec, n, H, ldH, S, ldS, iwhich, eval, Z, info_out)
  !Complex Hermitian counterpart of RG_0S's gpu_dsygvx: computes the iwhich-th
  !eigenvalue (and, if jobz_vec/=0, its eigenvector) of H*x=E*S*x with
  !cusolverDnZhegvdx. Like ZHEGVX(...,'U',...) in workproc, it reads the upper
  !triangles of H and S. The eigenvalues of a Hermitian problem are real.
    integer     :: jobz_vec, n, ldH, ldS, iwhich, info_out
    complex(wp) :: H(ldH,n), S(ldS,n), Z(*)
    real(wp)    :: eval
    !Handle, device buffers and host staging are MODULE-scope + save (see top):
    !host staging is heap so it is NOT a K*K automatic array on the stack, and the
    !handle can be destroyed in gpu_finalize.
    integer     :: lwork, meig, jobz, uplo
    complex(wp) :: hZ(n)    !length-n vector (small), safe on the stack
    real(wp)    :: hW1(1)
    integer     :: hinfo(1)

    if (.not. eig_created) then
      call solver_check(cusolverDnCreate(eig_hdl), 0, 'gpu_zhegvx: cusolverDnCreate')
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
    call cuda_check(cudaMemcpy(dA, hA, n*n), 'gpu_zhegvx: H2D A')   !array-section device
    call cuda_check(cudaMemcpy(dB, hB, n*n), 'gpu_zhegvx: H2D B')   !assignment ICEs nvfortran

    jobz = CUSOLVER_EIG_MODE_NOVECTOR
    if (jobz_vec /= 0) jobz = CUSOLVER_EIG_MODE_VECTOR
    uplo = CUBLAS_FILL_MODE_UPPER

    call solver_check(cusolverDnZhegvdx_bufferSize(eig_hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, lwork), 0, 'gpu_zhegvx: Zhegvdx_bufferSize')
    if (lwork > eig_cap_work) then
      if (allocated(dwork)) deallocate(dwork)
      allocate(dwork(lwork)); eig_cap_work = lwork
    endif
    call solver_check(cusolverDnZhegvdx(eig_hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, dwork, lwork, dinfo(1)), 0, 'gpu_zhegvx: Zhegvdx')
    call cuda_check(cudaDeviceSynchronize(), 'gpu_zhegvx: solve sync')

    call cuda_check(cudaMemcpy(hinfo, dinfo, 1), 'gpu_zhegvx: D2H info')
    info_out = hinfo(1)
    !Return numerical failures to the existing G-method recovery path.
    !API/status failures still abort in solver_check.
    if (info_out /= 0) return
    if (meig < 1) then
      info_out = -999
      return
    endif
    call cuda_check(cudaMemcpy(hW1, dW, 1), 'gpu_zhegvx: D2H eigenvalue')
    eval = hW1(1)               !selected eigenvalue (first of meig)
    if (jobz_vec /= 0) then
      call cuda_check(cudaMemcpy(hZ, dA, n), 'gpu_zhegvx: D2H eigenvector')
      Z(1:n) = hZ(1:n)          !first column of A (the eigenvector)
    endif
  end subroutine gpu_zhegvx

end module gpu_backend
