!Module qrlinalg provides a QR-based shifted inverse-iteration eigensolver for
!the generalized symmetric eigenvalue problem
!
!                 H*x = lambda*S*x .
!
!For a real shift sigma the module represents the shifted matrix
!
!                 A = H - sigma*S
!
!of the current order n by two ingredients:
!
!  1. a BASE FACTORIZATION  A_base = Q*R  of order n0, computed at a snapshot
!     and never modified afterwards. The rows of the explicit orthogonal factor
!     Q and of the upper-triangular factor R are distributed cyclically over the
!     MPI ranks of MPI_COMM_WORLD (row i lives on rank mod(i-1,P));
!
!  2. a BORDER of at most max_pending "pending" functions whose current row and
!     column differ from the snapshot, or that were appended after it. Their
!     current columns are stored and are never applied to Q and R. Linear
!     systems with the current A are solved by block elimination on the
!     unmodified base matrix.
!
!Row and column updates are therefore deferred (batched): replacing or
!appending one function costs one distributed solve with the base factors and
!O(max_pending*n) local work, instead of O(n**2) Givens rotations applied to
!explicit Q and R. The only O(n**3) operation is a fresh factorization, which
!the caller requests when the border is full ("flush").
!
!Mathematical scheme. Let D0 be the set of base slots (indices not greater than
!n0) whose current row and column differ from the base: every pending function
!with an index not greater than n0 and every trailing base slot deleted by
!truncate. Let I = {1..n0} \ D0 be the clean set. The functions appended
!beyond the base (index greater than n0) are pending as well. On the base
!indices the current matrix is a symmetric low-rank modification of A_base,
!and the appended functions border it. Both are expressed by one bordered
!("augmented") system on the full base matrix:
!
!   [ A_base   E     D     B  ] [ x0 ]   [ r0 ]
!   [ E^T      0    -I     0  ] [ yE ] = [ 0  ]
!   [ D^T     -I     DD    K  ] [ yD ]   [ 0  ]
!   [ B^T      0     K^T   CA ] [ xA ]   [ rA ]
!
!   E   the unit vectors e_k of the members k of D0               (n0 by m0)
!   D   current minus base column of every member of D0 on the rows of I,
!       zero on the rows of D0                                    (n0 by m0)
!   B   current column of every appended function on the rows of I, zero on
!       the rows of D0                                            (n0 by mA)
!   DD  current minus base values of the D0-D0 block
!   K   current values of the D0-appended block
!   CA  current values of the appended-appended block.
!
!Eliminating yE = D^T*x0 + DD*yD + K*xA and yD = E^T*x0 shows that x0 and xA
!solve the current system A*x = r: x0 is the solution on the base indices
!(r0 is r on the base indices, zero on deleted slots) and xA the solution on
!the appended functions. A deleted trailing base slot is a member of D0 whose
!current row and column are zero except for the diagonal, which decouples it.
!
!With the border columns Bt = [E | D | B] the module stores, replicated on
!every rank,
!
!        Zb(:,q)  = A_base^{-1}*e_{D0(q)}                 (n0, per member of D0)
!        Ab(:,q)  = A_base*e_{D0(q)}, the base column     (n0, per member of D0)
!        Bt(:,s)  column of D or of B                     (n0, per border entry)
!        Gs(:,s)  = A_base^{-1}*Bt(:,s)                   (n0, per border entry)
!        C        current border-border values
!        BG = Bt^T*Gs,  BZ = Bt^T*Zb
!        Sc = Caug - [E D B]^T * A_base^{-1} * [E D B]    (Schur complement),
!
!where Caug is the lower-right 3 by 3 block of the augmented matrix. A*x = r
!is solved by
!
!        g  = A_base^{-1}*r0                           (one base solve)
!        t  = [0; 0; rA] - [E D B]^T*g
!        y  = Sc^{-1}*t = [yE; yD; xA]
!        x0 = g - Zb*yE - Gs*[yD; xA] .
!
!No explicit inverse of a large matrix is formed; all corrections are applied
!to the result of the base solve. The last Ma border entries form the ACTIVE
!TAIL (the functions being optimized now); the other border entries and all
!unit-vector unknowns yE are COMMITTED. Sc is partitioned accordingly; LU(Scc)
!is computed when the committed set changes, and X = Scc^{-1}*Sca,
!Sigma = Saa - Sac*X and LU(Sigma) are recomputed lazily at the first solve
!after a tail column changed.
!
!Why the base matrix is modified and not restricted. An alternative is to
!exclude the rows and columns of D0 from the base and to border the clean
!block A_II with the full current columns. The inverse of A_II must then be
!expressed through A_base^{-1} as a difference of terms that are larger than
!the result by the ratio of the near-singularity of A_base to that of A_II,
!whatever the size of the pending changes; this loses many digits when the
!pending functions carry a substantial part of the wanted eigenvector. In the
!scheme above the unknowns yE = D^T*x0 + ... multiply the large columns Zb,
!so the cancelling terms are proportional to the size of the pending CHANGES,
!which is small for the trial steps of an optimizer.
!
!Distribution. Every procedure of the state type is COLLECTIVE: all ranks call
!it with identical (replicated) arguments and obtain identical results and
!status values. H, S, right-hand sides, eigenvectors, the border data and all
!small matrices are replicated; only the rows of Q and R are distributed. Each
!rank stores its rows transposed, so that a row is a contiguous column of the
!local arrays. Reductions over the ranks use MPI_ALLREDUCE on MPI_WP data. When
!the program runs with one process, no MPI routine is called with wp data
!(gfortran quadruple precision cannot be reduced by Open MPI), so the serial
!path works in every precision.
!
!Numerical regime. The shifted matrix is deliberately nearly singular (the
!shift lies very close to the wanted eigenvalue). Solutions of linear systems
!then carry large errors along the wanted eigenvector, which inverse iteration
!tolerates, and small errors in all other directions. A small diagonal element
!of R or a small pivot of the active Schur complement Sigma is therefore
!accepted; only exact zeros, subnormal and non-finite values are rejected
!there. The committed block Scc is tested with a relative pivot threshold
!that detects pivots lost to cancellation, and a failure is reported as
!QR_ERR_SINGULAR so the caller can refactorize. Inverse iteration applies the
!bordered solve in residual-correction form after the first iteration (see
!real_solve), which makes the attainable accuracy of the eigenpair
!independent of the error growth of block elimination.
!
!The caller owns H, S, and all vectors passed through the public interface.
!Only the canonical LOWER triangles of H and S are ever referenced. A state
!does not retain pointers to caller arrays. All storage is allocated by
!initialize, except for the temporaries that rank zero needs during a fresh
!factorization; replace_symmetric, append_symmetric, truncate, begin_active
!and solve do not allocate memory.
!
!Fresh factorization. Two algorithms fill the same distributed layout. The
!serial one runs the bundled LAPACK routines on rank zero and ships the rows
!of the factors to their owners; its cost does not decrease with the number
!of ranks. The distributed one is a blocked Householder QR factorization
!that works on the distributed rows: every rank forms its rows of the shifted
!matrix from the replicated lower triangles, the matrix is reduced to R in
!place, and Q is accumulated from the reflectors in the same layout, so that
!no rank holds an n by n temporary (a rank allocates work arrays of about
!2*(1+1/P)*QR_FACTOR_BLOCK*n elements instead). qr_set_factorization, or the
!environment variable ECG_QR_FACTOR, selects one of the algorithms or leaves
!the choice to a timing model that is calibrated once per run (the default;
!see real_factor_select). With one process the default is the serial
!algorithm; the distributed algorithm, when it is requested, then makes no
!MPI call either.
!
!The working kind wp is selected when the library is compiled. The bundled
!LAPACK routines DGEQRF and DORGQR used by the fresh factorization are generic
!in wp. The real-specific kernels are kept in procedures named real_*, so that
!a complex Hermitian twin can be added mechanically.
!
!Complex Hermitian twin. Type qr_complex_state implements the same scheme with
!the same interface for the generalized Hermitian eigenvalue problem with
!complex Hermitian H and S (complex basis functions, code CG_0S). The shift,
!the eigenvalue, the tolerances and the residuals remain real. Its procedures
!are named complex_* and correspond one to one to the real_* procedures; the
!border bookkeeping is identical. The differences are the following.
!
!  - Storage convention. Only the canonical lower triangles of H and S are
!    referenced. The stored element (i,j), i>=j, is <phi_i|O|phi_j>, with the
!    row function as the bra; the conceptual element (i,j) with i<j is the
!    complex conjugate of the stored element (j,i). The columns passed to
!    replace_symmetric and append_symmetric are conceptual columns: element i
!    of the column of function k is <phi_i|O|phi_k>. The row of the function
!    is the conjugate of its column and is never supplied.
!  - Diagonal elements are real. factorize_fresh, replace_symmetric and
!    append_symmetric return QR_ERR_INVALID_ARGUMENT, leaving the state
!    unchanged, when the imaginary part of a diagonal element of H or of S
!    exceeds QR_DIAGONAL_TOLERANCE times the scale of the supplied data;
!    smaller imaginary parts are discarded. The products with H and S ignore
!    the imaginary parts of the diagonals.
!  - Every transpose of the real scheme becomes a conjugate transpose: the
!    base solve is R^{-1}*(Q^H*b), the augmented matrix has the blocks E^H,
!    D^H, B^H and K^H below the diagonal, BG = Bt^H*Gs, BZ = Bt^H*Zb, the
!    right-hand side of the small system is [0; 0; rA] - [E D B]^H*g, and the
!    border-border block is filled as a Hermitian matrix. The Schur
!    complement and its two LU factorizations are general complex matrices.
!  - Each rank stores the complex conjugates of its rows of Q, that is, the
!    corresponding columns of Q^H, so that the base solve needs no
!    conjugation in its inner loop.
!  - Inner products conjugate their first argument. The Rayleigh quotients
!    and the quadratic form x^H*S*x are the real parts of the computed
!    values. Iterates are scaled with the largest absolute real or imaginary
!    component, the direction change is measured with the complex projection
!    coefficient (v^H*x)/(v^H*v), and the normalizations are x^H*S*x = 1,
!    x^H*x = 1, or unit largest real or imaginary component.
!  - Complex data are sent and reduced as real data of twice the length with
!    the datatype MPI_WP. The serial fresh factorization uses the bundled
!    ZGEQRF and ZUNGQR; the distributed one generates its reflectors with the
!    conventions of ZLARFG (real diagonal of R), applies them to Q itself,
!    and conjugates the stored rows at the end. The two algorithms are
!    selected for both types together (qr_set_factorization, ECG_QR_FACTOR),
!    but each type has its own calibrated timing model.
module qrlinalg
  use iso_fortran_env, only: int64
  use wp_def
  implicit none
  private

  integer, parameter, public :: QR_SUCCESS = 0
    !The requested operation completed successfully.
  integer, parameter, public :: QR_ERR_INVALID_ARGUMENT = 1
    !A scalar control, index, or mathematical input property is invalid.
  integer, parameter, public :: QR_ERR_ALLOCATION = 2
    !Initialization or a fresh factorization could not allocate storage.
  integer, parameter, public :: QR_ERR_NOT_IMPLEMENTED = 3
    !The requested numerical operation is not implemented in this version.
  integer, parameter, public :: QR_ERR_FACTORIZATION = 4
    !A bundled LAPACK factorization or workspace query reported an error.
  integer, parameter, public :: QR_ERR_SINGULAR = 5
    !The base factor, a committed small matrix, or an iterate is unusable.
  integer, parameter, public :: QR_ERR_NO_CONVERGENCE = 6
    !Inverse iteration reached max_iter before satisfying the tolerance.
  integer, parameter, public :: QR_ERR_INVALID_STATE = 7
    !The state is uninitialized, invalid, or missing required owned storage.
  integer, parameter, public :: QR_ERR_DIMENSION_MISMATCH = 8
    !A caller-owned array has an incompatible or empty extent.
  integer, parameter, public :: QR_ERR_CAPACITY_EXCEEDED = 9
    !The requested order is larger than the initialized capacity.
  integer, parameter, public :: QR_ERR_ZERO_INITIAL_VECTOR = 10
    !The supplied inverse-iteration starting vector is numerically zero.
  integer, parameter, public :: QR_ERR_NONPOSITIVE_OVERLAP = 11
    !The final vector has a non-positive or numerically zero overlap norm.
  integer, parameter, public :: QR_ERR_BATCH_FULL = 12
    !The border (batch of deferred updates) cannot hold the request; the
    !caller must refactorize with factorize_fresh.
  integer, parameter :: QR_STATUS_MESSAGE_LENGTH = 96
    !Fixed result length used by qr_status_message; returned text is padded.

  integer, parameter :: QR_TAG_R = 7411
    !Message tag of the rows of R shipped by rank zero after DGEQRF.
  integer, parameter :: QR_TAG_Q = 7412
    !Message tag of the rows of Q shipped by rank zero after DORGQR.
  integer, parameter :: QR_NUM_WORK_VECTORS = 6
    !Number of capacity-length work vectors owned by a state.
  integer, parameter :: QR_MAX_POLISH = 4
    !Largest number of residual-correction steps of one solve.
  integer, parameter :: QR_NUM_SMALL_VECTORS = 3
    !Number of small work vectors (length 2*max_pending) owned by a state.
  logical, parameter :: QR_PRODUCT_PAIRS = (digits(1.0_wp) == 64)
    !True for 80-bit extended precision: the symmetric product then processes
    !two columns per pass instead of four.
  real(wp), parameter :: QR_RELATIVE_PIVOT = 100.0_wp * epsilon(1.0_wp)
    !Relative pivot threshold of the LU factorization of Scc.
  real(wp), parameter :: QR_DIAGONAL_TOLERANCE = 100.0_wp * epsilon(1.0_wp)
    !Complex state only: largest accepted imaginary part of a diagonal element
    !of H or S, relative to the scale of the supplied matrix or column.

  integer, parameter, public :: QR_FACTOR_AUTO = 0
    !Fresh factorization: the module chooses the algorithm (the default).
  integer, parameter, public :: QR_FACTOR_SERIAL = 1
    !Fresh factorization by LAPACK on rank zero; rows shipped afterwards.
  integer, parameter, public :: QR_FACTOR_DISTRIBUTED = 2
    !Fresh factorization by the distributed blocked Householder algorithm.
  integer, parameter :: QR_FACTOR_ENV_UNREAD = -1
    !The environment variable ECG_QR_FACTOR has not been read yet.
  integer, parameter :: QR_FACTOR_ENV_UNSET = -2
    !ECG_QR_FACTOR is absent or does not name an algorithm.
  integer, parameter :: QR_FACTOR_BLOCK = 32
    !Number of columns of a panel of the distributed factorization.
  integer, parameter :: QR_FACTOR_CACHE_BITS = 8 * 65536
    !Size in bits of the part of the replicated matrix W that the kernels of
    !the distributed factorization keep in work at a time.
  logical, parameter :: QR_FACTOR_DOT_KERNELS = (digits(1.0_wp) == 64)
    !True for 80-bit extended precision: the two kernels of the distributed
    !factorization then accumulate inner products in registers
    !(real_dot_multiply), because storing an 80-bit value is much more
    !expensive than loading one; otherwise they update whole columns, which
    !the compiler vectorizes.
  integer, parameter :: QR_FACTOR_DOT_ROWS = 128
    !Number of local rows that the inner-product kernels process at a time.
  integer, parameter :: QR_CALIBRATION_SMALL = &
    merge(128, 64, digits(1.0_wp) <= 64)
  integer, parameter :: QR_CALIBRATION_LARGE = 2 * QR_CALIBRATION_SMALL
    !The two matrix orders of the calibration of the timing models (smaller
    !for software-emulated quadruple precision).
  real(wp), parameter :: QR_SAFE_MINIMUM = tiny(1.0_wp) / epsilon(1.0_wp)
    !Smallest magnitude of a diagonal element of R that a Householder
    !reflector generates without rescaling (as in the LAPACK routine DLARFG).
  real(wp), parameter :: QR_SMALL_NORM2 = QR_SAFE_MINIMUM / epsilon(1.0_wp)
    !A squared column norm below this value is recomputed with scaling.
  integer, parameter :: QR_SCALE_EXPONENT_HIGH = maxexponent(1.0_wp) / 2 - 32
  integer, parameter :: QR_SCALE_EXPONENT_LOW2 = &
    minexponent(1.0_wp) + 2 * digits(1.0_wp) + 32
    !The distributed factorization scales the matrix by a power of two when
    !the binary exponent e of its largest element is above the first value
    !or 2*e is below the second one.

  integer, save :: qr_factor_request = QR_FACTOR_AUTO
    !Algorithm requested by qr_set_factorization.
  integer, save :: qr_factor_environment = QR_FACTOR_ENV_UNREAD
    !Algorithm requested by ECG_QR_FACTOR, or one of the two codes above.
  real(kind(1.0d0)), save :: real_factor_model(4) = -1.0d0
    !Timing models of the two factorization algorithms of the real state,
    !see real_factor_select; a negative first element means not calibrated.
  real(kind(1.0d0)), save :: complex_factor_model(4) = -1.0d0
    !The same timing models for the complex state (complex_factor_select).

  public :: wp
  public :: qr_real_state
  public :: qr_complex_state
  public :: qr_status_message
  public :: qr_set_factorization
#ifdef QRLINALG_TESTING
  public :: real_factor_model
  public :: complex_factor_model
#endif

  !Interfaces to the external LAPACK routines supplied in src/LAPACK.f. The
  !routine names retain their traditional D prefixes, but their scalar kind is
  !the compile-time wp.
  interface
    ! Compute a compact real QR factorization in A and TAU.
    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, lda, lwork
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: tau(*)
      real(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dgeqrf

    ! Generate explicit real Q from the reflectors returned by DGEQRF.
    subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, k, lda, lwork
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(in) :: tau(*)
      real(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dorgqr
  end interface

  !Interfaces to the complex LAPACK routines supplied in src/LAPACK.f. Their
  !scalar kind is the compile-time wp as well.
  interface
    ! Compute a compact complex QR factorization in A and TAU.
    subroutine zgeqrf(m, n, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, lda, lwork
      complex(wp), intent(inout) :: a(lda, *)
      complex(wp), intent(out) :: tau(*)
      complex(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zgeqrf

    ! Generate explicit complex Q from the reflectors returned by ZGEQRF.
    subroutine zungqr(m, n, k, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, k, lda, lwork
      complex(wp), intent(inout) :: a(lda, *)
      complex(wp), intent(in) :: tau(*)
      complex(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zungqr
  end interface

  type :: qr_real_state
#ifndef QRLINALG_TESTING
    private
#endif
    !The components are private in a normal library build. QRLINALG_TESTING is
    !defined only for a white-box test build.
    !
    !--- Orders, shift, validity -------------------------------------------
    ! Current order of the represented matrix. Zero means no factorization.
    integer :: n = 0
    ! Order of the base factorization A_base = Q*R.
    integer :: n0 = 0
    ! Largest order that fits in the allocated storage.
    integer :: capacity = 0
    ! Largest number of border entries.
    integer :: mmax = 0
    ! Shift represented by the base factorization and the border data.
    real(wp) :: shift = 0.0_wp
    ! True only after a successful fresh factorization.
    logical :: valid = .false.
    ! True when a diagonal element of R is zero, subnormal, or not finite.
    logical :: r_singular = .false.
    ! Algorithm that computed the base factorization: QR_FACTOR_SERIAL or
    ! QR_FACTOR_DISTRIBUTED; zero before the first fresh factorization.
    integer :: factor_kind = 0
    !--- Distribution -------------------------------------------------------
    ! Communicator (MPI_COMM_WORLD), number of ranks P, and own rank.
    integer :: comm = 0
    integer :: nprocs = 1
    integer :: rank = 0
    ! Number of local rows for the capacity and for the base order.
    integer :: nloc_cap = 0
    integer :: nloc = 0
    ! Number of P by P diagonal blocks of the base R.
    integer :: nblk = 0
    ! qt(1:n0,l) is row i of Q and rt(1:n0,l) is row i of R, where
    ! i = rank+1+(l-1)*P. Leading dimension capacity.
    real(wp), allocatable :: qt(:,:)
    real(wp), allocatable :: rt(:,:)
    ! rd(:,:,b) is the diagonal block of R on rows and columns
    ! (b-1)*P+1..min(b*P,n0), replicated on all ranks.
    real(wp), allocatable :: rd(:,:,:)
    !--- Border bookkeeping -------------------------------------------------
    ! m: number of border entries. An entry belongs to a pending function or
    ! to a deleted trailing base slot (index above n, not above n0).
    ! ma: length of the active tail (the last ma entries).
    ! nres: number of announced appended indices beyond the base that have
    ! not been supplied yet.
    ! m0: number of members of D0 (the entries with index not above n0).
    integer :: m = 0
    integer :: ma = 0
    integer :: nres = 0
    integer :: m0 = 0
    ! bidx(p) is the canonical index of the border entry at position p and
    ! bslot(p) its storage slot (column of bt and gs, row/column of cbb and
    ! bg). bslot(m+1:mmax) lists the free slots. bpos(i) is the border
    ! position of canonical index i, or zero.
    integer, allocatable :: bidx(:)
    integer, allocatable :: bslot(:)
    integer, allocatable :: bpos(:)
    ! d0idx(q) is the base slot of member q of D0; d0pos(i) is the position
    ! of base slot i in D0, or zero. Member q owns column q of zb and ab.
    integer, allocatable :: d0idx(:)
    integer, allocatable :: d0pos(:)
    ! stale(s) is nonzero when the entry in slot s is a deleted base slot
    ! whose decoupling ("dummy") column has not been installed yet.
    integer, allocatable :: stale(:)
    ! cu(1:m0+mc) lists the committed unknowns of the small system in the
    ! order of LU(Scc); valid while scc_ok is true.
    integer, allocatable :: cu(:)
    !--- Replicated border data ---------------------------------------------
    ! bt(:,s): border column of entry slot s (current minus base column on
    ! the clean rows for a member of D0, current column on the clean rows
    ! for an appended function; zero on the rows of D0).
    ! gs(:,s) = A_base^{-1}*bt(:,s).
    ! zb(:,q) = A_base^{-1}*e and ab(:,q) = A_base*e for member q of D0.
    real(wp), allocatable :: bt(:,:)
    real(wp), allocatable :: gs(:,:)
    real(wp), allocatable :: zb(:,:)
    real(wp), allocatable :: ab(:,:)
    ! cbb: current border-border values (slot by slot); bg = Bt^T*Gs (slot by
    ! slot); bz = Bt^T*Zb (slot by D0 position).
    real(wp), allocatable :: cbb(:,:)
    real(wp), allocatable :: bg(:,:)
    real(wp), allocatable :: bz(:,:)
    ! sc: Schur complement of the augmented system. Unknown q (1..m0) is the
    ! unit-vector unknown yE of member q of D0; unknown mmax+s belongs to
    ! the border entry in slot s. Dimension 2*mmax by 2*mmax.
    real(wp), allocatable :: sc(:,:)
    ! LU factors with partial pivoting of Scc (order of cu) and of Sigma
    ! (tail order); scc_mag is the magnitude bookkeeping of the pivot test;
    ! xca = Scc^{-1}*Sca.
    real(wp), allocatable :: scc_lu(:,:)
    real(wp), allocatable :: scc_mag(:,:)
    real(wp), allocatable :: sig_lu(:,:)
    real(wp), allocatable :: xca(:,:)
    integer, allocatable :: scc_piv(:)
    integer, allocatable :: sig_piv(:)
    ! scc_ok: LU(Scc) matches the committed set. tail_ok: xca and LU(Sigma)
    ! match the tail columns. scc_singular records a failed relative pivot
    ! test of Scc; sig_singular an unusable pivot of Sigma.
    logical :: scc_ok = .true.
    logical :: tail_ok = .true.
    logical :: scc_singular = .false.
    logical :: sig_singular = .false.
    !--- Workspace ------------------------------------------------------------
    ! wn: capacity-length vectors. Columns 1 and 2 belong to the bordered
    ! solve and to the bookkeeping; columns 3 to 6 to the public procedures.
    real(wp), allocatable :: wn(:,:)
    ! wsm: three vectors of length 2*mmax of the bordered solve.
    real(wp), allocatable :: wsm(:,:)
    ! wloc: one value per local row. wblk: one value per rank.
    real(wp), allocatable :: wloc(:)
    real(wp), allocatable :: wblk(:)
    ! wred: reduction buffer of the BG/BZ rebuild (2*mmax**2).
    real(wp), allocatable :: wred(:)
    ! imark: capacity-length marks, zero between calls. iw: two mmax-length
    ! integer work vectors.
    integer, allocatable :: imark(:)
    integer, allocatable :: iw(:,:)
    !--- Counters ---------------------------------------------------------------
    ! Count of successful replace, append, and truncate/delete operations.
    integer(int64) :: structural_updates = 0_int64
    ! Successful structural updates since the last fresh factorization.
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => real_initialize
    procedure :: clear => clear_real_state
    procedure :: factorize_fresh => real_factorize_fresh
    procedure :: factorization_kind => real_factorization_kind
    procedure :: begin_active => real_begin_active
    procedure :: replace_symmetric => real_replace_symmetric
    procedure :: append_symmetric => real_append_symmetric
    procedure :: truncate => real_truncate
    procedure :: delete_symmetric => real_delete_symmetric
    procedure :: solve => real_solve
    procedure :: solve_shifted => real_solve_shifted
    procedure :: symmetric_product => real_symmetric_product
    procedure :: shifted_quotient => real_shifted_quotient
    procedure :: factorization_residual => real_factorization_residual
    procedure :: is_valid => real_is_valid
    procedure :: order => real_order
    procedure :: base_order => real_base_order
    procedure :: get_capacity => real_get_capacity
    procedure :: get_shift => real_get_shift
    procedure :: pending_count => real_pending_count
    procedure :: active_count => real_active_count
    procedure :: excluded_count => real_excluded_count
    procedure :: max_pending => real_max_pending
    procedure :: get_update_count => real_get_update_count
    procedure :: get_updates_since_fresh => real_get_updates_since_fresh
  end type qr_real_state

  type :: qr_complex_state
#ifndef QRLINALG_TESTING
    private
#endif
    !The components are private in a normal library build. QRLINALG_TESTING is
    !defined only for a white-box test build.
    !
    !--- Orders, shift, validity -------------------------------------------
    ! Current order of the represented matrix. Zero means no factorization.
    integer :: n = 0
    ! Order of the base factorization A_base = Q*R.
    integer :: n0 = 0
    ! Largest order that fits in the allocated storage.
    integer :: capacity = 0
    ! Largest number of border entries.
    integer :: mmax = 0
    ! Shift represented by the base factorization and the border data.
    real(wp) :: shift = 0.0_wp
    ! True only after a successful fresh factorization.
    logical :: valid = .false.
    ! True when a diagonal element of R is zero, subnormal, or not finite.
    logical :: r_singular = .false.
    ! Algorithm that computed the base factorization: QR_FACTOR_SERIAL or
    ! QR_FACTOR_DISTRIBUTED; zero before the first fresh factorization.
    integer :: factor_kind = 0
    !--- Distribution -------------------------------------------------------
    ! Communicator (MPI_COMM_WORLD), number of ranks P, and own rank.
    integer :: comm = 0
    integer :: nprocs = 1
    integer :: rank = 0
    ! Number of local rows for the capacity and for the base order.
    integer :: nloc_cap = 0
    integer :: nloc = 0
    ! Number of P by P diagonal blocks of the base R.
    integer :: nblk = 0
    ! qt(1:n0,l) is the complex conjugate of row i of Q (column i of Q^H)
    ! and rt(1:n0,l) is row i of R, where i = rank+1+(l-1)*P. Leading
    ! dimension capacity.
    complex(wp), allocatable :: qt(:,:)
    complex(wp), allocatable :: rt(:,:)
    ! rd(:,:,b) is the diagonal block of R on rows and columns
    ! (b-1)*P+1..min(b*P,n0), replicated on all ranks.
    complex(wp), allocatable :: rd(:,:,:)
    !--- Border bookkeeping -------------------------------------------------
    ! m: number of border entries. An entry belongs to a pending function or
    ! to a deleted trailing base slot (index above n, not above n0).
    ! ma: length of the active tail (the last ma entries).
    ! nres: number of announced appended indices beyond the base that have
    ! not been supplied yet.
    ! m0: number of members of D0 (the entries with index not above n0).
    integer :: m = 0
    integer :: ma = 0
    integer :: nres = 0
    integer :: m0 = 0
    ! bidx(p) is the canonical index of the border entry at position p and
    ! bslot(p) its storage slot (column of bt and gs, row/column of cbb and
    ! bg). bslot(m+1:mmax) lists the free slots. bpos(i) is the border
    ! position of canonical index i, or zero.
    integer, allocatable :: bidx(:)
    integer, allocatable :: bslot(:)
    integer, allocatable :: bpos(:)
    ! d0idx(q) is the base slot of member q of D0; d0pos(i) is the position
    ! of base slot i in D0, or zero. Member q owns column q of zb and ab.
    integer, allocatable :: d0idx(:)
    integer, allocatable :: d0pos(:)
    ! stale(s) is nonzero when the entry in slot s is a deleted base slot
    ! whose decoupling ("dummy") column has not been installed yet.
    integer, allocatable :: stale(:)
    ! cu(1:m0+mc) lists the committed unknowns of the small system in the
    ! order of LU(Scc); valid while scc_ok is true.
    integer, allocatable :: cu(:)
    !--- Replicated border data ---------------------------------------------
    ! bt(:,s): border column of entry slot s (current minus base column on
    ! the clean rows for a member of D0, current column on the clean rows
    ! for an appended function; zero on the rows of D0).
    ! gs(:,s) = A_base^{-1}*bt(:,s).
    ! zb(:,q) = A_base^{-1}*e and ab(:,q) = A_base*e for member q of D0.
    complex(wp), allocatable :: bt(:,:)
    complex(wp), allocatable :: gs(:,:)
    complex(wp), allocatable :: zb(:,:)
    complex(wp), allocatable :: ab(:,:)
    ! cbb: current border-border values (slot by slot), a Hermitian matrix
    ! whose element (s1,s2) is the matrix element (bidx of s1, bidx of s2);
    ! bg = Bt^H*Gs (slot by slot); bz = Bt^H*Zb (slot by D0 position).
    complex(wp), allocatable :: cbb(:,:)
    complex(wp), allocatable :: bg(:,:)
    complex(wp), allocatable :: bz(:,:)
    ! sc: Schur complement of the augmented system. Unknown q (1..m0) is the
    ! unit-vector unknown yE of member q of D0; unknown mmax+s belongs to
    ! the border entry in slot s. Dimension 2*mmax by 2*mmax.
    complex(wp), allocatable :: sc(:,:)
    ! LU factors with partial pivoting of Scc (order of cu) and of Sigma
    ! (tail order); scc_mag is the real magnitude bookkeeping of the pivot
    ! test; xca = Scc^{-1}*Sca.
    complex(wp), allocatable :: scc_lu(:,:)
    real(wp), allocatable :: scc_mag(:,:)
    complex(wp), allocatable :: sig_lu(:,:)
    complex(wp), allocatable :: xca(:,:)
    integer, allocatable :: scc_piv(:)
    integer, allocatable :: sig_piv(:)
    ! scc_ok: LU(Scc) matches the committed set. tail_ok: xca and LU(Sigma)
    ! match the tail columns. scc_singular records a failed relative pivot
    ! test of Scc; sig_singular an unusable pivot of Sigma.
    logical :: scc_ok = .true.
    logical :: tail_ok = .true.
    logical :: scc_singular = .false.
    logical :: sig_singular = .false.
    !--- Workspace ------------------------------------------------------------
    ! wn: capacity-length vectors. Columns 1 and 2 belong to the bordered
    ! solve and to the bookkeeping; columns 3 to 6 to the public procedures.
    complex(wp), allocatable :: wn(:,:)
    ! wsm: three vectors of length 2*mmax of the bordered solve.
    complex(wp), allocatable :: wsm(:,:)
    ! wloc: one value per local row. wblk: one value per rank.
    complex(wp), allocatable :: wloc(:)
    complex(wp), allocatable :: wblk(:)
    ! wred: reduction buffer of the BG/BZ rebuild (2*mmax**2).
    complex(wp), allocatable :: wred(:)
    ! imark: capacity-length marks, zero between calls. iw: two mmax-length
    ! integer work vectors.
    integer, allocatable :: imark(:)
    integer, allocatable :: iw(:,:)
    !--- Counters
    ! --------------------------------------------------------------- Count of
    ! successful replace, append, and truncate/delete operations.
    integer(int64) :: structural_updates = 0_int64
    ! Successful structural updates since the last fresh factorization.
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => complex_initialize
    procedure :: clear => clear_complex_state
    procedure :: factorize_fresh => complex_factorize_fresh
    procedure :: factorization_kind => complex_factorization_kind
    procedure :: begin_active => complex_begin_active
    procedure :: replace_symmetric => complex_replace_symmetric
    procedure :: append_symmetric => complex_append_symmetric
    procedure :: truncate => complex_truncate
    procedure :: delete_symmetric => complex_delete_symmetric
    procedure :: solve => complex_solve
    procedure :: solve_shifted => complex_solve_shifted
    procedure :: symmetric_product => complex_symmetric_product
    procedure :: shifted_quotient => complex_shifted_quotient
    procedure :: factorization_residual => complex_factorization_residual
    procedure :: is_valid => complex_is_valid
    procedure :: order => complex_order
    procedure :: base_order => complex_base_order
    procedure :: get_capacity => complex_get_capacity
    procedure :: get_shift => complex_get_shift
    procedure :: pending_count => complex_pending_count
    procedure :: active_count => complex_active_count
    procedure :: excluded_count => complex_excluded_count
    procedure :: max_pending => complex_max_pending
    procedure :: get_update_count => complex_get_update_count
    procedure :: get_updates_since_fresh => complex_get_updates_since_fresh
  end type qr_complex_state

contains

  pure function qr_status_message(info) result(message)
  !Function qr_status_message returns a stable, human-readable description of
  !a public qrlinalg status value. The result describes the library condition
  !only; it does not classify the condition as fatal, retryable, or otherwise
  !impose caller error-handling policy.
  !
  !  Input parameter:
  !    info - Integer status returned by a qrlinalg operation or another
  !           integer that is to be interpreted as a qrlinalg status.
  !
  !  Result:
  !    message - Blank-padded status description with fixed length
  !              QR_STATUS_MESSAGE_LENGTH. Callers normally use TRIM(message)
  !              when writing or composing diagnostics. An unrecognized value
  !              returns a stable fallback description rather than failing.
  !
  !The function performs no allocation or input/output, changes no state, and
  !is safe to call from pure procedures. Programs must continue to branch on
  !the integer status symbols rather than parse the descriptive text.
    integer, intent(in) :: info
    character(len=QR_STATUS_MESSAGE_LENGTH) :: message

    select case (info)
    case (QR_SUCCESS)
      message = 'operation completed successfully'
    case (QR_ERR_INVALID_ARGUMENT)
      message = 'invalid scalar control, index, or mathematical input property'
    case (QR_ERR_ALLOCATION)
      message = 'state storage allocation failed'
    case (QR_ERR_NOT_IMPLEMENTED)
      message = 'requested numerical operation is not implemented'
    case (QR_ERR_FACTORIZATION)
      message = 'LAPACK workspace query or QR factorization failed'
    case (QR_ERR_SINGULAR)
      message = 'shifted factorization or generated iterate is numerically unusable'
    case (QR_ERR_NO_CONVERGENCE)
      message = 'inverse iteration did not converge within max_iter; approximation is available'
    case (QR_ERR_INVALID_STATE)
      message = 'QR state is uninitialized, invalid, or missing owned storage'
    case (QR_ERR_DIMENSION_MISMATCH)
      message = 'caller array extents are empty or incompatible'
    case (QR_ERR_CAPACITY_EXCEEDED)
      message = 'requested active order exceeds initialized capacity'
    case (QR_ERR_ZERO_INITIAL_VECTOR)
      message = 'inverse-iteration starting vector is numerically zero'
    case (QR_ERR_NONPOSITIVE_OVERLAP)
      message = 'overlap quadratic form is non-positive or numerically zero'
    case (QR_ERR_BATCH_FULL)
      message = 'batch of deferred updates is full; a fresh factorization is required'
    case default
      message = 'unrecognized qrlinalg status code'
    end select
  end function qr_status_message

  subroutine clear_real_state(self)
  !Subroutine clear_real_state releases every allocation owned by a real QR
  !state and restores the state to its default, uninitialized condition. It is
  !the implementation of the public clear() binding and is also used before
  !reinitialization and after an allocation failure. Calling clear() for an
  !already empty state is valid.
  !
  !  Input/output parameter:
  !    self - The real QR state. On exit all allocatable components are
  !           unallocated; all orders, counts and counters are zero; valid is
  !           false.
  !
  !The routine makes no MPI call and has no failure status. Like every other
  !procedure it must be called on all ranks so that the replicated states stay
  !identical. After clear(), initialize must be called before another
  !factorization can be constructed.
    class(qr_real_state), intent(inout) :: self

    if (allocated(self%qt)) deallocate(self%qt)
    if (allocated(self%rt)) deallocate(self%rt)
    if (allocated(self%rd)) deallocate(self%rd)
    if (allocated(self%bidx)) deallocate(self%bidx)
    if (allocated(self%bslot)) deallocate(self%bslot)
    if (allocated(self%bpos)) deallocate(self%bpos)
    if (allocated(self%d0idx)) deallocate(self%d0idx)
    if (allocated(self%d0pos)) deallocate(self%d0pos)
    if (allocated(self%stale)) deallocate(self%stale)
    if (allocated(self%cu)) deallocate(self%cu)
    if (allocated(self%bt)) deallocate(self%bt)
    if (allocated(self%gs)) deallocate(self%gs)
    if (allocated(self%zb)) deallocate(self%zb)
    if (allocated(self%ab)) deallocate(self%ab)
    if (allocated(self%cbb)) deallocate(self%cbb)
    if (allocated(self%bg)) deallocate(self%bg)
    if (allocated(self%bz)) deallocate(self%bz)
    if (allocated(self%sc)) deallocate(self%sc)
    if (allocated(self%scc_lu)) deallocate(self%scc_lu)
    if (allocated(self%scc_mag)) deallocate(self%scc_mag)
    if (allocated(self%sig_lu)) deallocate(self%sig_lu)
    if (allocated(self%xca)) deallocate(self%xca)
    if (allocated(self%scc_piv)) deallocate(self%scc_piv)
    if (allocated(self%sig_piv)) deallocate(self%sig_piv)
    if (allocated(self%wn)) deallocate(self%wn)
    if (allocated(self%wsm)) deallocate(self%wsm)
    if (allocated(self%wloc)) deallocate(self%wloc)
    if (allocated(self%wblk)) deallocate(self%wblk)
    if (allocated(self%wred)) deallocate(self%wred)
    if (allocated(self%imark)) deallocate(self%imark)
    if (allocated(self%iw)) deallocate(self%iw)
    self%n = 0
    self%n0 = 0
    self%capacity = 0
    self%mmax = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%r_singular = .false.
    self%comm = 0
    self%nprocs = 1
    self%rank = 0
    self%nloc_cap = 0
    self%nloc = 0
    self%nblk = 0
    self%m = 0
    self%ma = 0
    self%nres = 0
    self%m0 = 0
    self%scc_ok = .true.
    self%tail_ok = .true.
    self%scc_singular = .false.
    self%sig_singular = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_real_state

  subroutine real_initialize(self, capacity, max_pending, info)
  !Subroutine real_initialize prepares a real QR state for symmetric matrices
  !of order not greater than capacity and for at most max_pending deferred
  !(border) functions. Any factorization and counters previously held by the
  !state are discarded. The routine is collective.
  !
  !The rank and the number of processes of MPI_COMM_WORLD are obtained here and
  !fix the row distribution for the lifetime of the state. With P ranks, the
  !following storage is allocated on every rank:
  !  - two capacity by ceiling(capacity/P) arrays for the local rows of Q and
  !    of R (transposed layout);
  !  - the replicated P by P diagonal blocks of R (about P*capacity elements);
  !  - four capacity by mmax arrays for the border columns Bt, their base
  !    solves Gs, and the base solves Zb and base columns Ab of the members
  !    of D0;
  !  - small matrices of order mmax and 2*mmax and their pivot vectors;
  !  - six capacity-length work vectors and small work arrays,
  !where mmax = min(max_pending,capacity). Initialization does not form a QR
  !factorization; therefore the order is zero and valid is false on exit.
  !
  !  Input parameters:
  !    capacity    - Maximum matrix order supported by the state; positive.
  !    max_pending - Maximum number of border entries (functions whose updates
  !                  are deferred, and deleted trailing base slots); positive.
  !                  Values above capacity are reduced to capacity;
  !                  max_pending() returns the value in effect.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Existing allocations and factorization
  !            metadata are destroyed before the arguments are validated.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when all storage is ready on all ranks;
  !            QR_ERR_INVALID_ARGUMENT when capacity or max_pending is not
  !            positive;
  !            QR_ERR_INVALID_STATE when MPI has not been initialized;
  !            QR_ERR_ALLOCATION when an allocation fails on any rank.
  !
  !The allocation status is reduced with MPI_MAX over the ranks before any
  !rank branches on it, so that all ranks return the same info. If
  !initialization fails, self is returned in the empty state described for
  !clear_real_state on every rank.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: capacity, max_pending
    integer, intent(out) :: info
    integer :: allocation_status, block_capacity, i, mpi_err, pending
    logical :: mpi_is_ready

    !Discard all previous storage before constructing the new state.
    call clear_real_state(self)
    if (capacity <= 0 .or. max_pending <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    call MPI_Initialized(mpi_is_ready, mpi_err)
    if (.not. mpi_is_ready) then
      info = QR_ERR_INVALID_STATE
      return
    end if

    self%comm = MPI_COMM_WORLD
    call MPI_Comm_size(self%comm, self%nprocs, mpi_err)
    call MPI_Comm_rank(self%comm, self%rank, mpi_err)
    pending = min(max_pending, capacity)
    self%nloc_cap = qr_local_count(capacity, self%rank, self%nprocs)
    block_capacity = (capacity + self%nprocs - 1) / self%nprocs

    allocate(self%qt(capacity, self%nloc_cap), &
             self%rt(capacity, self%nloc_cap), &
             self%rd(self%nprocs, self%nprocs, block_capacity), &
             self%bidx(pending), self%bslot(pending), self%bpos(capacity), &
             self%d0idx(pending), self%d0pos(capacity), &
             self%stale(pending), self%cu(2 * pending), &
             self%bt(capacity, pending), self%gs(capacity, pending), &
             self%zb(capacity, pending), self%ab(capacity, pending), &
             self%cbb(pending, pending), self%bg(pending, pending), &
             self%bz(pending, pending), &
             self%sc(2 * pending, 2 * pending), &
             self%scc_lu(2 * pending, 2 * pending), &
             self%scc_mag(2 * pending, 2 * pending), &
             self%sig_lu(pending, pending), self%xca(2 * pending, pending), &
             self%scc_piv(2 * pending), self%sig_piv(pending), &
             self%wn(capacity, QR_NUM_WORK_VECTORS), &
             self%wsm(2 * pending, QR_NUM_SMALL_VECTORS), &
             self%wloc(max(1, self%nloc_cap)), self%wblk(self%nprocs), &
             self%wred(2 * pending * pending), &
             self%imark(capacity), self%iw(pending, 2), &
             stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    !A rank-local failure must become a global one before anyone branches.
    if (self%nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%qt = 0.0_wp
    self%rt = 0.0_wp
    self%rd = 0.0_wp
    self%bt = 0.0_wp
    self%gs = 0.0_wp
    self%zb = 0.0_wp
    self%ab = 0.0_wp
    self%cbb = 0.0_wp
    self%bg = 0.0_wp
    self%bz = 0.0_wp
    self%sc = 0.0_wp
    self%scc_lu = 0.0_wp
    self%scc_mag = 0.0_wp
    self%sig_lu = 0.0_wp
    self%xca = 0.0_wp
    self%scc_piv = 0
    self%sig_piv = 0
    self%wn = 0.0_wp
    self%wsm = 0.0_wp
    self%wloc = 0.0_wp
    self%wblk = 0.0_wp
    self%wred = 0.0_wp
    self%imark = 0
    self%iw = 0
    self%bidx = 0
    self%bpos = 0
    self%d0idx = 0
    self%d0pos = 0
    self%stale = 0
    self%cu = 0
    do i = 1, pending
      self%bslot(i) = i
    end do

    !Only the capacities are committed here. The order, the represented shift,
    !and the validity flag are committed by a successful factorize_fresh.
    self%capacity = capacity
    self%mmax = pending
    info = QR_SUCCESS
  end subroutine real_initialize

  subroutine real_factorize_fresh(self, h, s, shift, info, active_order)
  !Subroutine real_factorize_fresh constructs the base QR factorization of the
  !shifted real symmetric matrix
  !
  !                     A_base = H - shift*S = Q*R
  !
  !of order n and empties the border: after a successful call the base order
  !equals the current order, nothing is pending, and no base slot is excluded.
  !This is the only O(n**3) operation of the module ("flush"). The routine is
  !collective.
  !
  !Two algorithms compute the factors and fill the same distributed layout:
  !the serial one (real_factor_serial: LAPACK on rank zero, rows shipped to
  !their owners) and the distributed blocked Householder factorization
  !(real_factor_distributed). real_factor_select decides which of them runs;
  !its decision is identical on all ranks. See qr_set_factorization.
  !
  !Neither H nor S is modified or retained. When active_order is absent, H and
  !S must be square and have equal order. When active_order is present, H and
  !S may have larger physical storage and need only have at least n rows and n
  !columns. In both forms, only the lower triangles of the leading n by n
  !principal blocks are referenced.
  !
  !  Input parameters:
  !    h     - Real array whose leading n by n lower triangle contains the
  !            symmetric Hamiltonian matrix H (replicated).
  !    s     - Real array whose leading n by n lower triangle contains the
  !            symmetric overlap matrix S (replicated).
  !    shift - The real shift represented by the factorization.
  !    active_order - Optional order n. It must be positive, must not exceed
  !                   the capacity, and must fit within both extents of H and
  !                   S. If absent, n is the common square order of H and S.
  !
  !  Input/output parameter:
  !    self  - An initialized real QR state. On successful exit the order and
  !            the base order are n, the shift is stored, valid is true, the
  !            border and D0 are empty, and updates_since_fresh is zero.
  !            structural_updates is a lifetime counter and is not reset.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when the factorization succeeds;
  !            QR_ERR_INVALID_STATE when self has not been initialized;
  !            QR_ERR_INVALID_ARGUMENT when active_order is present but is
  !            not positive;
  !            QR_ERR_DIMENSION_MISMATCH when the original no-active-order
  !            matrix rules are violated, or when an explicit active block
  !            does not fit within H or S;
  !            QR_ERR_CAPACITY_EXCEEDED when n is greater than the capacity;
  !            QR_ERR_ALLOCATION when a rank cannot allocate the temporaries
  !            of the factorization (reported on all ranks);
  !            QR_ERR_FACTORIZATION when a LAPACK routine of the serial
  !            algorithm reports an error.
  !
  !All argument checks and the allocation of the temporaries are completed
  !before existing factors are overwritten; such an error preserves the
  !previous state. A LAPACK error leaves the state invalid (order zero).
  !An exactly singular or non-finite R does not fail here: as in the previous
  !implementation it is reported as QR_ERR_SINGULAR by the procedures that
  !need to solve with R.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer, intent(in), optional :: active_order
    integer :: blk, ii, matrix_n, np
    logical :: distributed

    !Check initialization, matrix shapes, active order, and capacity before
    !overwriting any component of a previously valid factorization.
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%qt)) return
    if (present(active_order)) then
      matrix_n = active_order
      if (matrix_n <= 0) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      if (matrix_n > self%capacity) then
        info = QR_ERR_CAPACITY_EXCEEDED
        return
      end if
      info = QR_ERR_DIMENSION_MISMATCH
      if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n) return
      if (size(s, 1) < matrix_n .or. size(s, 2) < matrix_n) return
    else
      info = QR_ERR_DIMENSION_MISMATCH
      if (size(h, 1) /= size(h, 2)) return
      if (size(s, 1) /= size(s, 2)) return
      if (size(h, 1) /= size(s, 1)) return
      matrix_n = size(h, 1)
      if (matrix_n <= 0) return
    end if
    if (matrix_n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if

    np = self%nprocs

    !Choose the algorithm and compute the factors. Both algorithms fill
    !qt, rt, rd, nloc and nblk, and both leave the previous state untouched
    !when they fail before the factors are overwritten.
    call real_factor_select(self%comm, np, self%rank, matrix_n, distributed)
    if (distributed) then
      call real_factor_distributed(self, matrix_n, h, s, shift, info)
    else
      call real_factor_serial(self, matrix_n, h, s, shift, info)
    end if
    if (info /= QR_SUCCESS) return

    !A QR factorization can complete for an exactly rank-deficient matrix.
    !Record whether the diagonal of R can be used by an unguarded triangular
    !solve. Every rank evaluates the replicated diagonal blocks and obtains
    !the same flag.
    self%r_singular = .false.
    do blk = 1, self%nblk
      do ii = 1, min(np, matrix_n - (blk - 1) * np)
        !The negated form also rejects NaN.
        if (.not. (abs(self%rd(ii, ii, blk)) > tiny(1.0_wp) .and. &
                   abs(self%rd(ii, ii, blk)) <= huge(1.0_wp))) then
          self%r_singular = .true.
        end if
      end do
    end do

    !Commit the factorization only after every numerical stage succeeds.
    self%n = matrix_n
    self%n0 = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    if (distributed) then
      self%factor_kind = QR_FACTOR_DISTRIBUTED
    else
      self%factor_kind = QR_FACTOR_SERIAL
    end if
    info = QR_SUCCESS
  end subroutine real_factorize_fresh

  !===========================================================================
  ! Fresh factorization: choice of the algorithm, the serial algorithm, and
  ! the distributed blocked Householder algorithm with its kernels
  !===========================================================================

  subroutine qr_set_factorization(choice, info)
  !Subroutine qr_set_factorization selects the algorithm that factorize_fresh
  !uses from now on, for all states of the module (real and complex):
  !
  !   QR_FACTOR_SERIAL      - LAPACK on rank zero; the rows of R and Q are
  !                           shipped to their owners (stage A);
  !   QR_FACTOR_DISTRIBUTED - distributed blocked Householder factorization
  !                           (stage B);
  !   QR_FACTOR_AUTO        - the default: with one process the serial
  !                           algorithm; with several processes the faster
  !                           one according to a timing model calibrated at
  !                           the first factorization of the run (see
  !                           real_factor_select).
  !
  !The environment variable ECG_QR_FACTOR (serial, distributed, or auto),
  !when it is set to one of these words, takes precedence over the value set
  !here, so that the algorithm can be changed without recompiling.
  !
  !The routine makes no MPI call. Like every procedure of the module it must
  !be called with the same argument on all ranks; ranks with different
  !choices would deadlock in the next fresh factorization.
  !
  !  Input parameter:
  !    choice - One of the three QR_FACTOR_* values.
  !
  !  Output parameter:
  !    info   - Optional. QR_SUCCESS, or QR_ERR_INVALID_ARGUMENT when choice is
  !             none of the three values (the setting is then unchanged).
    integer, intent(in) :: choice
    integer, intent(out), optional :: info

    if (choice == QR_FACTOR_AUTO .or. choice == QR_FACTOR_SERIAL .or. &
        choice == QR_FACTOR_DISTRIBUTED) then
      qr_factor_request = choice
      if (present(info)) info = QR_SUCCESS
    else
      if (present(info)) info = QR_ERR_INVALID_ARGUMENT
    end if
  end subroutine qr_set_factorization

  subroutine qr_factor_resolve(comm, nprocs, rank, choice)
  !Subroutine qr_factor_resolve returns the requested factorization algorithm
  !in effect: the value of the environment variable ECG_QR_FACTOR if it is
  !set to serial, distributed, or auto (letter case is ignored), otherwise
  !the value set by qr_set_factorization. The environment is read once per
  !run, on rank zero, at the first call, and its integer code is broadcast,
  !so that all ranks use the same value whatever their own environments
  !contain. The routine is collective.
  !
  !  Input parameters:
  !    comm, nprocs, rank - Communicator, number of ranks, and own rank.
  !
  !  Output parameter:
  !    choice - QR_FACTOR_AUTO, QR_FACTOR_SERIAL, or QR_FACTOR_DISTRIBUTED.
    integer, intent(in) :: comm, nprocs, rank
    integer, intent(out) :: choice
    character(len=32) :: text
    integer :: code, i, length, mpi_err, status

    if (qr_factor_environment == QR_FACTOR_ENV_UNREAD) then
      code = QR_FACTOR_ENV_UNSET
      if (rank == 0) then
        call get_environment_variable('ECG_QR_FACTOR', text, length, status)
        !A nonzero status means that the variable does not exist or does not
        !fit in text; neither can be one of the three words.
        if (status == 0 .and. length > 0) then
          do i = 1, len(text)
            if (text(i:i) >= 'A' .and. text(i:i) <= 'Z') then
              text(i:i) = achar(iachar(text(i:i)) + 32)
            end if
          end do
          select case (trim(adjustl(text)))
          case ('serial')
            code = QR_FACTOR_SERIAL
          case ('distributed')
            code = QR_FACTOR_DISTRIBUTED
          case ('auto')
            code = QR_FACTOR_AUTO
          end select
        end if
      end if
      if (nprocs > 1) then
        call MPI_Bcast(code, 1, MPI_INTEGER, 0, comm, mpi_err)
      end if
      qr_factor_environment = code
    end if
    choice = qr_factor_request
    if (qr_factor_environment /= QR_FACTOR_ENV_UNSET) then
      choice = qr_factor_environment
    end if
  end subroutine qr_factor_resolve

  pure function qr_first_local(k, rank, nprocs) result(l)
  !Function qr_first_local returns the position l of the first local row of a
  !rank whose global index rank+1+(l-1)*nprocs is not smaller than k. The
  !result exceeds the number of local rows when the rank owns no such row.
  !
  !  Input parameters:
  !    k      - Global row index, k >= 1.
  !    rank   - Rank, 0 <= rank < nprocs.
  !    nprocs - Number of ranks.
  !
  !  Result:
  !    l - Local position, l >= 1.
    integer, intent(in) :: k, rank, nprocs
    integer :: l

    if (k <= rank + 1) then
      l = 1
    else
      l = (k - rank - 2) / nprocs + 2
    end if
  end function qr_first_local

  pure function real_factorization_kind(self) result(kind_used)
  !Function real_factorization_kind reports which algorithm computed the base
  !factorization that the state holds.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    kind_used - QR_FACTOR_SERIAL or QR_FACTOR_DISTRIBUTED; zero when the
  !                state holds no valid factorization.
    class(qr_real_state), intent(in) :: self
    integer :: kind_used

    kind_used = 0
    if (self%valid) kind_used = self%factor_kind
  end function real_factorization_kind

  subroutine real_factor_select(comm, nprocs, rank, n, distributed)
  !Subroutine real_factor_select decides whether a fresh factorization of
  !order n runs the distributed algorithm. The routine is collective and its
  !result is identical on all ranks: it depends only on the requested
  !algorithm (qr_factor_resolve), on the number of ranks, and, for the
  !automatic choice with several ranks, on the timing model
  !
  !        t_serial(n)      = model(1)*n**3 + model(2)*n**2
  !        t_distributed(n) = model(3)*n**3 + model(4)*n
  !
  !whose coefficients are measured once per run by real_factor_calibrate, on
  !rank zero, and broadcast. With one process the automatic choice is the
  !serial algorithm (LAPACK without any message).
  !
  !  Input parameters:
  !    comm, nprocs, rank - Communicator, number of ranks, and own rank.
  !    n - Order of the matrix to factorize.
  !
  !  Output parameter:
  !    distributed - True when the distributed algorithm is to be used.
    integer, intent(in) :: comm, nprocs, rank, n
    logical, intent(out) :: distributed
    integer :: choice
    real(kind(1.0d0)) :: order

    call qr_factor_resolve(comm, nprocs, rank, choice)
    if (choice == QR_FACTOR_SERIAL) then
      distributed = .false.
    else if (choice == QR_FACTOR_DISTRIBUTED) then
      distributed = .true.
    else if (nprocs == 1) then
      distributed = .false.
    else
      if (real_factor_model(1) < 0.0d0) then
        call real_factor_calibrate(comm, nprocs, rank)
      end if
      order = real(n, kind(1.0d0))
      distributed = (real_factor_model(3) * order**3 + &
                     real_factor_model(4) * order < &
                     real_factor_model(1) * order**3 + &
                     real_factor_model(2) * order**2)
    end if
  end subroutine real_factor_select

  subroutine real_factor_calibrate(comm, nprocs, rank)
  !Subroutine real_factor_calibrate measures, on the actual machine and with
  !the actual number of ranks, how fast the two factorization algorithms are,
  !and stores the coefficients of the two timing models of real_factor_select
  !in real_factor_model. It follows the calibration of the LDL factorizations
  !in module linalg (linalg_ldlf_calibrate).
  !
  !Both algorithms factorize a synthetic well-conditioned matrix of two small
  !orders in a scratch state; the shortest wall time of rank zero, taken
  !between two barriers, is used. The serial algorithm is timed NUM_REPEAT
  !times and the much cheaper distributed one twice as often, in alternating
  !passes, so that a short disturbance of the node cannot spoil all timings
  !of one kind; an overestimated distributed time is the error that could
  !cost most, because it would keep large matrices on one rank.
  !
  !The cubic terms describe the arithmetic. The quadratic term of the serial
  !model is the shipping of the rows, and the linear term of the distributed
  !model is the latency of the reductions (a fixed number per column).
  !Negative coefficients, which can result from timing noise, are set to zero
  !and the other coefficient is then refitted to the larger order. Rank zero
  !fits the models and broadcasts them, so that every rank takes identical
  !decisions afterwards.
  !
  !If the scratch state cannot be allocated, the model is set such that the
  !distributed algorithm is always chosen.
  !
  !The routine is collective and must be called with more than one rank. It
  !costs a few factorizations of order QR_CALIBRATION_LARGE.
  !
  !  Input parameters:
  !    comm, nprocs, rank - Communicator, number of ranks, and own rank.
    integer, intent(in) :: comm, nprocs, rank
    integer, parameter :: NUM_REPEAT = 2
    type(qr_real_state) :: scratch
    real(wp), allocatable :: a(:,:)
    real(kind(1.0d0)) :: elapsed(2, 2), model(4), start
    integer :: algorithm, allocation_status, i, info, isize, j, mpi_err
    integer :: order(2), rep
    logical :: failed

    order(1) = QR_CALIBRATION_SMALL
    order(2) = QR_CALIBRATION_LARGE
    model(1) = 1.0d0
    model(2) = 0.0d0
    model(3) = 0.0d0
    model(4) = 0.0d0

    allocate(a(order(2), order(2)), stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, comm, mpi_err)
    end if
    info = QR_ERR_ALLOCATION
    if (allocation_status == 0) then
      call real_initialize(scratch, order(2), 1, info)
    end if
    if (info == QR_SUCCESS) then
      !Synthetic diagonally dominant matrix; every leading block is well
      !conditioned. It serves as H, and the shift is zero.
      do j = 1, order(2)
        do i = 1, order(2)
          a(i, j) = 1.0_wp / real(1 + abs(i - j), wp)
        end do
        a(j, j) = a(j, j) + 3.0_wp
      end do
      failed = .false.
      elapsed = huge(1.0d0)
      do rep = 1, 2 * NUM_REPEAT
        do isize = 1, 2
          do algorithm = 1, 2
            if (algorithm == 1 .and. rep > NUM_REPEAT) cycle
            if (nprocs > 1) call MPI_Barrier(comm, mpi_err)
            start = MPI_Wtime()
            if (algorithm == 1) then
              call real_factor_serial(scratch, order(isize), a, a, 0.0_wp, &
                                      info)
            else
              call real_factor_distributed(scratch, order(isize), a, a, &
                                           0.0_wp, info)
            end if
            if (nprocs > 1) call MPI_Barrier(comm, mpi_err)
            elapsed(isize, algorithm) = min(elapsed(isize, algorithm), &
                                            MPI_Wtime() - start)
            if (info /= QR_SUCCESS) failed = .true.
          end do
        end do
      end do
      if (rank == 0 .and. .not. failed) then
        call fit(elapsed(1, 1), elapsed(2, 1), 2, model(1), model(2))
        call fit(elapsed(1, 2), elapsed(2, 2), 1, model(3), model(4))
      end if
      if (nprocs > 1) then
        call MPI_Bcast(model, 4, MPI_DOUBLE_PRECISION, 0, comm, mpi_err)
      end if
    end if
    call clear_real_state(scratch)
    if (allocated(a)) deallocate(a)
    real_factor_model = model

  contains

    subroutine fit(t1, t2, low_power, cubic, low)
    !Fits t = cubic*n**3 + low*n**low_power to the two measurements.
      real(kind(1.0d0)), intent(in) :: t1, t2
      integer, intent(in) :: low_power
      real(kind(1.0d0)), intent(out) :: cubic, low
      real(kind(1.0d0)) :: n1, n2

      n1 = real(order(1), kind(1.0d0))
      n2 = real(order(2), kind(1.0d0))
      cubic = (t2 / n2**low_power - t1 / n1**low_power) / &
              (n2**(3 - low_power) - n1**(3 - low_power))
      if (cubic < 0.0d0) then
        cubic = 0.0d0
        low = t2 / n2**low_power
      else
        low = t1 / n1**low_power - cubic * n1**(3 - low_power)
        if (low < 0.0d0) then
          low = 0.0d0
          cubic = t2 / n2**3
        end if
      end if
    end subroutine fit

  end subroutine real_factor_calibrate

  subroutine real_factor_serial(self, matrix_n, h, s, shift, info)
  !Subroutine real_factor_serial computes the base factors of the shifted
  !matrix H - shift*S of order matrix_n with LAPACK on rank zero and
  !distributes them (stage A). It is called by real_factorize_fresh, which
  !has validated the arguments and which commits the factorization, and by
  !real_factor_calibrate.
  !
  !Stage A implementation: rank zero forms the full matrix from the two lower
  !triangles in a temporary n by n array, calls the bundled DGEQRF, ships to
  !every rank its rows of R (row i goes to rank mod(i-1,P)), broadcasts the
  !P by P diagonal blocks of R, calls DORGQR to overwrite the temporary by the
  !explicit Q, ships the rows of Q, and releases its temporaries (one n by n
  !array, tau, one row buffer, and the LAPACK workspace). The other ranks only
  !receive. With one process no message is sent.
  !
  !  Input parameters:
  !    matrix_n - Order n; positive and not above the capacity.
  !    h, s     - Arrays whose leading n by n lower triangles hold H and S.
  !    shift    - The shift.
  !
  !  Input/output parameter:
  !    self - Initialized state. On success qt, rt, rd, nloc and nblk describe
  !           the new factors; the order, the validity flag and r_singular are
  !           left for the caller to set (the state is marked invalid). The
  !           border is emptied.
  !
  !  Output parameter:
  !    info - QR_SUCCESS; QR_ERR_ALLOCATION when rank zero cannot allocate its
  !           temporaries (state unchanged); QR_ERR_FACTORIZATION when a
  !           workspace query (state unchanged) or DGEQRF or DORGQR (state
  !           invalid) reports an error. Identical on all ranks.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: matrix_n
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    real(wp), allocatable :: a(:,:), tau(:), work(:), row(:)
    real(wp) :: work_query(1)
    integer :: allocation_status, blk, dest, i, ii, j, l, lapack_info, lwork
    integer :: mpi_err, np
    logical :: is_root

    np = self%nprocs
    is_root = (self%rank == 0)

    !Rank zero allocates the temporaries of the serial factorization. The
    !status is made global before any rank branches on it.
    allocation_status = 0
    lapack_info = 0
    lwork = 1
    if (is_root) then
      allocate(a(matrix_n, matrix_n), tau(matrix_n), row(matrix_n), &
               stat=allocation_status)
      if (allocation_status == 0) then
        !LWORK=-1 performs no factorization and returns the recommended
        !workspace length in WORK(1).
        call dgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work_query, -1, &
                    lapack_info)
        if (lapack_info == 0) then
          lwork = max(matrix_n, ceiling(work_query(1)))
          call dorgqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, &
                      work_query, -1, lapack_info)
        end if
        if (lapack_info == 0) then
          lwork = max(lwork, ceiling(work_query(1)))
          allocate(work(lwork), stat=allocation_status)
        end if
      end if
      if (allocation_status /= 0) allocation_status = 1
      !Encode a failed workspace query so that one reduction carries both.
      if (lapack_info /= 0) allocation_status = 2
    end if
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status == 1) then
      info = QR_ERR_ALLOCATION
      return
    else if (allocation_status /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !The factors are about to be overwritten. Mark the state invalid until
    !every stage has completed.
    self%valid = .false.
    self%n = 0
    self%n0 = 0
    self%nloc = 0
    self%nblk = 0
    self%r_singular = .false.
    call real_reset_border(self)

    !Read each stored lower-triangular column contiguously and reflect every
    !off-diagonal value into the opposite triangle. LAPACK requires the full
    !shifted matrix, but callers need define only the canonical lower triangle.
    lapack_info = 0
    if (is_root) then
      do j = 1, matrix_n
        do i = j, matrix_n
          a(i, j) = h(i, j) - shift * s(i, j)
          if (i /= j) a(j, i) = a(i, j)
        end do
      end do
      !On exit from DGEQRF the upper triangle contains R; the strict lower
      !triangle and tau describe Q.
      call dgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work, lwork, &
                  lapack_info)
    end if
    if (np > 1) then
      call MPI_Bcast(lapack_info, 1, MPI_INTEGER, 0, self%comm, mpi_err)
    end if
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    self%nloc = qr_local_count(matrix_n, self%rank, np)
    self%nblk = (matrix_n + np - 1) / np

    !Distribute the rows of R before DORGQR destroys the compact reflector
    !representation. Row i is padded with zeros on the left, so that every
    !stored row has a defined triangular meaning.
    if (is_root) then
      self%rd(:, :, 1:self%nblk) = 0.0_wp
      do i = 1, matrix_n
        do j = 1, i - 1
          row(j) = 0.0_wp
        end do
        do j = i, matrix_n
          row(j) = a(i, j)
        end do
        blk = (i - 1) / np + 1
        ii = i - (blk - 1) * np
        do j = i, min(blk * np, matrix_n)
          self%rd(ii, j - (blk - 1) * np, blk) = a(i, j)
        end do
        dest = mod(i - 1, np)
        if (dest == 0) then
          self%rt(1:matrix_n, blk) = row(1:matrix_n)
        else
          call MPI_Send(row, matrix_n, MPI_WP, dest, QR_TAG_R, self%comm, &
                        mpi_err)
        end if
      end do
    else
      do l = 1, self%nloc
        call MPI_Recv(self%rt(1, l), matrix_n, MPI_WP, 0, QR_TAG_R, &
                      self%comm, MPI_STATUS_IGNORE, mpi_err)
      end do
    end if
    if (np > 1) then
      call MPI_Bcast(self%rd, np * np * self%nblk, MPI_WP, 0, self%comm, &
                     mpi_err)
    end if

    !Expand the Householder vectors and tau into the explicit orthogonal
    !matrix Q and distribute its rows.
    if (is_root) then
      call dorgqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, work, &
                  lwork, lapack_info)
    end if
    if (np > 1) then
      call MPI_Bcast(lapack_info, 1, MPI_INTEGER, 0, self%comm, mpi_err)
    end if
    if (lapack_info /= 0) then
      self%nloc = 0
      self%nblk = 0
      info = QR_ERR_FACTORIZATION
      return
    end if
    if (is_root) then
      do i = 1, matrix_n
        dest = mod(i - 1, np)
        if (dest == 0) then
          blk = (i - 1) / np + 1
          do j = 1, matrix_n
            self%qt(j, blk) = a(i, j)
          end do
        else
          do j = 1, matrix_n
            row(j) = a(i, j)
          end do
          call MPI_Send(row, matrix_n, MPI_WP, dest, QR_TAG_Q, self%comm, &
                        mpi_err)
        end if
      end do
    else
      do l = 1, self%nloc
        call MPI_Recv(self%qt(1, l), matrix_n, MPI_WP, 0, QR_TAG_Q, &
                      self%comm, MPI_STATUS_IGNORE, mpi_err)
      end do
    end if

    info = QR_SUCCESS
  end subroutine real_factor_serial

  subroutine real_factor_distributed(self, matrix_n, h, s, shift, info)
  !Subroutine real_factor_distributed computes the base factors of the shifted
  !matrix A = H - shift*S of order n = matrix_n by a blocked Householder QR
  !factorization that works on the distributed rows (stage B). It is called
  !by real_factorize_fresh, which has validated the arguments and which
  !commits the factorization, and by real_factor_calibrate. It fills exactly
  !what real_factor_serial fills.
  !
  !Layout. Row i of the matrix lives on rank mod(i-1,P) and is stored as the
  !contiguous column l = (i-1)/P+1 of the local array rt, during and after
  !the factorization. Every rank forms its rows directly from the replicated
  !lower triangles of H and S; there is no n by n temporary on any rank.
  !
  !Reduction to R. The columns are processed in panels of QR_FACTOR_BLOCK
  !columns. A panel is copied to a compact buffer and factorized by
  !real_panel_factor (one reduction per column). The reflectors
  !H(j) = I - tau(j)*v(j)*v(j)^T of a panel, with the columns v(j) collected
  !in V (unit diagonal), form the block reflector
  !H(j1)*...*H(j2) = I - V*T*V^T with an upper triangular T (the compact WY
  !form of LAPACK). The columns to the right of the panel are updated as
  !
  !        A := (I - V*T^T*V^T)*A = A - (V*T^T)*W ,     W = V^T*A :
  !
  !every rank multiplies its rows of V and of A (real_block_product), one
  !reduction makes W global, and every rank subtracts its rows of V*T^T
  !times W (real_block_update). The same reduction carries V^T*V, from which
  !T is computed. A panel therefore costs one small reduction per column and
  !one reduction of W. On exit R occupies the upper triangle; the vectors v
  !are kept below the diagonal until Q has been formed, and T is saved.
  !
  !Explicit Q. Q = H(1)*H(2)*...*H(n) is accumulated in qt by applying the
  !block reflectors in reverse order to the identity,
  !
  !        Q := (I - V*T*V^T)*Q = Q - (V*T)*W ,          W = V^T*Q ,
  !
  !with the same two kernels and one reduction per panel. Finally the part of
  !rt below the diagonal is cleared and the diagonal blocks of R are
  !collected in rd by one reduction.
  !
  !Scaling. The reflectors are generated with the conventions and the
  !safeguards of the LAPACK routine DLARFG (real_reflector). Sums of squares
  !are formed without scaling, which is safe because the matrix is first
  !multiplied by a power of two when its largest element is so large or so
  !small that the squares could overflow or underflow (R is scaled back at
  !the end; both operations are exact), and because a column whose part below
  !the diagonal is tiny relative to the matrix is treated separately by
  !real_panel_factor.
  !
  !All quantities on which a rank branches (statuses, scale codes, tau, the
  !diagonal of R) come from reductions or from replicated data, so all ranks
  !take identical branches. With one process no MPI routine is called.
  !
  !  Input parameters:
  !    matrix_n - Order n; positive and not above the capacity.
  !    h, s     - Arrays whose leading n by n lower triangles hold H and S.
  !    shift    - The shift.
  !
  !  Input/output parameter:
  !    self - Initialized state. On success qt, rt, rd, nloc and nblk describe
  !           the new factors; the order, the validity flag and r_singular are
  !           left for the caller to set (the state is marked invalid). The
  !           border is emptied.
  !
  !  Output parameter:
  !    info - QR_SUCCESS, or QR_ERR_ALLOCATION when a rank cannot allocate its
  !           work arrays (about 2*(1+1/P)*QR_FACTOR_BLOCK*n elements; the
  !           state is then unchanged). Identical on all ranks.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: matrix_n
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    real(wp), allocatable :: pan(:,:), red(:), tau(:), tsave(:,:), wbuf(:)
    real(wp), allocatable :: work(:), ymat(:,:)
    real(wp) :: gscale
    integer :: allocation_status, code, i, j, j1, j2, jj, l, l1, la, lwork, m
    integer :: me, mpi_err, n, nb, nblk, nla, nloc, np, pb, power

    n = matrix_n
    np = self%nprocs
    me = self%rank
    nloc = qr_local_count(n, me, np)
    nb = min(QR_FACTOR_BLOCK, n)

    !Work arrays: the panel (later the explicit V of the panel) and V*T, one
    !column per local row; the T factors of all panels; the reduction
    !buffers; the transposed copies of the inner-product kernels. The status
    !is made global before any rank branches on it.
    lwork = 1
    if (QR_FACTOR_DOT_KERNELS) then
      lwork = QR_FACTOR_DOT_ROWS * nb + &
              max(16 * max(QR_FACTOR_DOT_ROWS, nb), &
                  QR_FACTOR_CACHE_BITS / storage_size(work))
    end if
    allocate(pan(nb, max(1, nloc)), ymat(nb, max(1, nloc)), tsave(nb, n), &
             tau(nb), red(2 * nb), wbuf((n + nb) * nb), work(lwork), &
             stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status /= 0) then
      info = QR_ERR_ALLOCATION
      return
    end if

    !The factors are about to be overwritten. Mark the state invalid until
    !every stage has completed.
    self%valid = .false.
    self%n = 0
    self%n0 = 0
    self%nloc = 0
    self%nblk = 0
    self%r_singular = .false.
    call real_reset_border(self)

    !Form the local rows of the shifted matrix. The part of row i to the left
    !of the diagonal is row i of the stored lower triangles; the part to the
    !right is column i of the lower triangles.
    code = -huge(0)
    do l = 1, nloc
      i = me + 1 + (l - 1) * np
      do j = 1, i - 1
        self%rt(j, l) = h(i, j) - shift * s(i, j)
      end do
      self%rt(i, l) = h(i, i) - shift * s(i, i)
      do j = i + 1, n
        self%rt(j, l) = h(j, i) - shift * s(j, i)
      end do
      code = max(code, real_scale_code(n, self%rt(1, l), 1))
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, code, 1, MPI_INTEGER, MPI_MAX, &
                         self%comm, mpi_err)
    end if

    !code is the binary exponent of the largest element. Scale the matrix by
    !a power of two when unscaled sums of squares of its elements could
    !overflow or underflow. A zero or non-finite matrix is not scaled.
    gscale = 1.0_wp
    if (abs(code) /= huge(0)) then
      if (code > QR_SCALE_EXPONENT_HIGH .or. &
          2 * code < QR_SCALE_EXPONENT_LOW2) then
        power = max(min(-code, maxexponent(1.0_wp) - 4), &
                    minexponent(1.0_wp) + 4)
        gscale = scale(1.0_wp, power)
        do l = 1, nloc
          do j = 1, n
            self%rt(j, l) = gscale * self%rt(j, l)
          end do
        end do
      end if
    end if

    !Reduction to upper triangular form, panel by panel. l1 is the first
    !local row that the panel touches and nla the number of such rows (not
    !positive when the rank owns none).
    do j1 = 1, n, nb
      pb = min(nb, n - j1 + 1)
      j2 = j1 + pb - 1
      l1 = qr_first_local(j1, me, np)
      nla = nloc - l1 + 1
      do la = 1, nla
        do jj = 1, pb
          pan(jj, la) = self%rt(j1 + jj - 1, l1 + la - 1)
        end do
      end do
      call real_panel_factor(self%comm, np, me, j1, pb, l1, nla, pan, nb, &
                             tau, red)
      !Columns j1..j2 are final: R on and above the diagonal, the reflector
      !vectors below it.
      do la = 1, nla
        do jj = 1, pb
          self%rt(j1 + jj - 1, l1 + la - 1) = pan(jj, la)
        end do
      end do
      !Turn the panel buffer into the explicit V: unit diagonal and zeros
      !above it in the rows j1..j2.
      do jj = 1, pb
        j = j1 + jj - 1
        if (mod(j - 1, np) == me) then
          la = (j - 1 - me) / np + 2 - l1
          pan(jj, la) = 1.0_wp
          do i = jj + 1, pb
            pan(i, la) = 0.0_wp
          end do
        end if
      end do
      !W = V^T*A for the m columns to the right of the panel, followed in the
      !same buffer (leading dimension m+pb) by the pb by pb matrix V^T*V.
      m = n - j2
      do i = 1, (m + pb) * pb
        wbuf(i) = 0.0_wp
      end do
      if (nla > 0) then
        if (m > 0) then
          call real_block_product(m, self%rt(j2 + 1, l1), self%capacity, &
                                  nla, pan, nb, pb, wbuf(1), m + pb, work)
        end if
        call real_block_gram(pb, nla, pan, nb, wbuf(m + 1), m + pb)
      end if
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, wbuf, (m + pb) * pb, MPI_WP, &
                           MPI_SUM, self%comm, mpi_err)
      end if
      call real_block_tfactor(pb, tau, wbuf(m + 1), m + pb, tsave(1, j1), nb)
      if (nla > 0 .and. m > 0) then
        call real_block_vt(pb, nla, pan, nb, tsave(1, j1), nb, .true., ymat, &
                           nb)
        call real_block_update(m, self%rt(j2 + 1, l1), self%capacity, nla, &
                               ymat, nb, pb, wbuf(1), m + pb, work)
      end if
    end do

    !Accumulate the explicit Q: start from the identity and apply the block
    !reflectors of the panels in reverse order. A block reflector acts on the
    !rows and columns j1..n only.
    do l = 1, nloc
      do j = 1, n
        self%qt(j, l) = 0.0_wp
      end do
      self%qt(me + 1 + (l - 1) * np, l) = 1.0_wp
    end do
    do j1 = ((n - 1) / nb) * nb + 1, 1, -nb
      pb = min(nb, n - j1 + 1)
      l1 = qr_first_local(j1, me, np)
      nla = nloc - l1 + 1
      !Explicit V of the panel from the part of rt below the diagonal.
      do la = 1, nla
        l = l1 + la - 1
        i = me + 1 + (l - 1) * np
        do jj = 1, pb
          j = j1 + jj - 1
          if (i > j) then
            pan(jj, la) = self%rt(j, l)
          else if (i == j) then
            pan(jj, la) = 1.0_wp
          else
            pan(jj, la) = 0.0_wp
          end if
        end do
      end do
      m = n - j1 + 1
      do i = 1, m * pb
        wbuf(i) = 0.0_wp
      end do
      if (nla > 0) then
        call real_block_product(m, self%qt(j1, l1), self%capacity, nla, pan, &
                                nb, pb, wbuf(1), m, work)
      end if
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, wbuf, m * pb, MPI_WP, MPI_SUM, &
                           self%comm, mpi_err)
      end if
      if (nla > 0) then
        call real_block_vt(pb, nla, pan, nb, tsave(1, j1), nb, .false., &
                           ymat, nb)
        call real_block_update(m, self%qt(j1, l1), self%capacity, nla, ymat, &
                               nb, pb, wbuf(1), m, work)
      end if
    end do

    !Clear the reflector vectors stored below the diagonal of R, so that
    !every stored row has a defined triangular meaning, and undo the scaling.
    do l = 1, nloc
      i = me + 1 + (l - 1) * np
      do j = 1, i - 1
        self%rt(j, l) = 0.0_wp
      end do
      if (gscale /= 1.0_wp) then
        do j = i, n
          self%rt(j, l) = self%rt(j, l) / gscale
        end do
      end if
    end do

    !Replicated diagonal blocks of R: every rank owns one row of each block
    !(its local row l belongs to block l); one reduction, in which each
    !element is the sum of one value and zeros, completes them.
    nblk = (n + np - 1) / np
    self%rd(:, :, 1:nblk) = 0.0_wp
    do l = 1, nloc
      do j = me + 1 + (l - 1) * np, min(l * np, n)
        self%rd(me + 1, j - (l - 1) * np, l) = self%rt(j, l)
      end do
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, self%rd, np * np * nblk, MPI_WP, &
                         MPI_SUM, self%comm, mpi_err)
    end if

    self%nloc = nloc
    self%nblk = nblk
    info = QR_SUCCESS
  end subroutine real_factor_distributed

  subroutine real_panel_factor(comm, np, me, j1, pb, l1, nla, pan, ldp, tau, &
                               red)
  !Subroutine real_panel_factor reduces one panel of the distributed matrix
  !to upper triangular form by Householder reflections. The panel consists of
  !the global columns j1..j1+pb-1 and of the rows j1..n; a rank holds its rows
  !of the panel as the columns of the compact buffer pan.
  !
  !Column j is processed as follows.
  !  1. Every rank accumulates, over its rows below the diagonal,
  !     z(c) = sum_i a(i,j)*a(i,c) for the panel columns c >= j (z(j) is the
  !     squared norm of the part of column j below the diagonal). The owner
  !     of row j adds the pivot row a(j,c). One reduction makes both global.
  !  2. Every rank generates the reflector from the pivot element and the
  !     norm (real_reflector): the new diagonal element beta, tau, and the
  !     scale factor of the vector v = (1, x/(alpha-beta)). The products
  !     v^T*a(:,c) follow from z and the pivot row without further
  !     communication.
  !  3. Every rank scales its part of v in place and applies the reflector to
  !     its rows of the panel columns to the right of j.
  !
  !A column whose squared norm below the diagonal is smaller than
  !QR_SMALL_NORM2 (in particular zero) or not finite is treated separately,
  !because its unscaled sum of squares is inaccurate: one integer reduction
  !gives the binary exponent of its largest element, and step 1 is repeated
  !with the column scaled by a power of two. This reproduces the safeguards
  !of DLARFG at the cost of two more reductions for such a column only.
  !
  !  Input parameters:
  !    comm, np, me - Communicator, number of ranks, and own rank.
  !    j1   - Global index of the first column of the panel.
  !    pb   - Number of columns of the panel.
  !    l1   - Local position of the first local row of the panel (global
  !           index not below j1).
  !    nla  - Number of local rows of the panel; may be zero or negative.
  !    ldp  - Leading dimension of pan, at least pb.
  !
  !  Input/output parameter:
  !    pan  - pan(jj,la) is the element in column j1+jj-1 of the local row
  !           l1+la-1. On exit the upper triangle of the panel holds R and
  !           the part below the diagonal the reflector vectors (without
  !           their unit diagonal elements).
  !
  !  Output parameter:
  !    tau  - The pb scalar factors of the reflectors; replicated.
  !
  !  Workspace:
  !    red  - Reduction buffer of length 2*pb.
    integer, intent(in) :: comm, np, me, j1, pb, l1, nla, ldp
    real(wp), intent(inout) :: pan(ldp, *)
    real(wp), intent(out) :: tau(pb)
    real(wp), intent(inout) :: red(2 * pb)
    real(wp) :: alpha, factor2, vl
    real(wp) :: beta, factor1, factor3, sc, xnorm2
    integer :: c, code, j, jj, la, lgt, lown, mpi_err, nc, nrows, power
    logical :: owner

    do jj = 1, pb
      j = j1 + jj - 1
      nc = pb - jj + 1
      !The local rows below the diagonal are the positions lgt..nla; the
      !owner of row j holds it at position lown.
      lgt = qr_first_local(j + 1, me, np) - l1 + 1
      nrows = max(nla - lgt + 1, 0)
      owner = (mod(j - 1, np) == me)
      lown = (j - 1 - me) / np + 2 - l1

      sc = 1.0_wp
      call accumulate()
      xnorm2 = red(1)
      !The negated form also takes NaN to the separate treatment.
      if (.not. (xnorm2 >= QR_SMALL_NORM2 .and. &
                 xnorm2 <= huge(1.0_wp))) then
        code = -huge(0)
        if (nrows > 0) code = real_scale_code(nrows, pan(jj, lgt), ldp)
        if (np > 1) then
          call MPI_Allreduce(MPI_IN_PLACE, code, 1, MPI_INTEGER, MPI_MAX, &
                             comm, mpi_err)
        end if
        if (code == -huge(0)) then
          !The column is exactly zero below the diagonal.
          xnorm2 = 0.0_wp
        else if (code /= huge(0)) then
          power = max(min(-code, maxexponent(1.0_wp) - 4), &
                      minexponent(1.0_wp) + 4)
          sc = scale(1.0_wp, power)
          call accumulate()
          xnorm2 = red(1)
        end if
      end if

      alpha = red(nc + 1)
      call real_reflector(alpha, xnorm2, sc, beta, tau(jj), factor1, &
                          factor2, factor3)
      !red(c) becomes tau*(v^T*a(:,c)), the row that the reflector subtracts
      !from the pivot row, for the panel columns to the right of j.
      do c = 2, nc
        red(c) = tau(jj) * (red(nc + c) + (factor2 * red(c)) * factor3)
      end do
      do la = lgt, nla
        vl = (pan(jj, la) * factor1) * factor2
        pan(jj, la) = vl
        do c = 2, nc
          pan(jj + c - 1, la) = pan(jj + c - 1, la) - vl * red(c)
        end do
      end do
      if (owner) then
        pan(jj, lown) = beta
        do c = 2, nc
          pan(jj + c - 1, lown) = pan(jj + c - 1, lown) - red(c)
        end do
      end if
    end do

  contains

    subroutine accumulate()
    !Step 1 for the current column with the scale factor sc; on exit red
    !holds z (nc values) followed by the pivot row (nc values), reduced.
      integer :: k

      if (nrows > 0) then
        call real_panel_column(nc, nrows, pan(jj, lgt), ldp, sc, red)
      else
        do k = 1, nc
          red(k) = 0.0_wp
        end do
      end if
      do k = 1, nc
        if (owner) then
          red(nc + k) = pan(jj + k - 1, lown)
        else
          red(nc + k) = 0.0_wp
        end if
      end do
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, red, 2 * nc, MPI_WP, MPI_SUM, comm, &
                           mpi_err)
      end if
    end subroutine accumulate

  end subroutine real_panel_factor

  pure subroutine real_panel_column(nc, nrows, x, ldx, sc, z)
  !Subroutine real_panel_column computes, for one column of a panel and the
  !local rows below its diagonal element, the products of the scaled column
  !with itself and with the panel columns to its right:
  !
  !        z(1) = sum_i (sc*x(1,i))**2 ,
  !        z(c) = sum_i sc*x(1,i)*x(c,i) ,   c = 2..nc ,
  !
  !where x(c,i) is the element of local row i in the c-th column counted from
  !the current one.
  !
  !  Input parameters:
  !    nc    - Number of panel columns from the current one to the last.
  !    nrows - Number of local rows below the diagonal; positive.
  !    x     - Rows of the panel, stored as columns, starting at the current
  !            panel column; leading dimension ldx.
  !    ldx   - Leading dimension of x.
  !    sc    - Scale factor of the current column (a power of two).
  !
  !  Output parameter:
  !    z     - The nc local sums.
    integer, intent(in) :: nc, nrows, ldx
    real(wp), intent(in) :: x(ldx, *)
    real(wp), intent(in) :: sc
    real(wp), intent(out) :: z(nc)
    integer :: c, i
    real(wp) :: xs

    do c = 1, nc
      z(c) = 0.0_wp
    end do
    do i = 1, nrows
      xs = x(1, i) * sc
      z(1) = z(1) + xs * xs
      do c = 2, nc
        z(c) = z(c) + xs * x(c, i)
      end do
    end do
  end subroutine real_panel_column

  pure function real_scale_code(n, x, incx) result(code)
  !Function real_scale_code describes the magnitude of the largest element of
  !a vector by an integer, which can be reduced over the ranks with MPI_MAX
  !without sending floating-point data:
  !
  !        huge(0)     when an element is not finite;
  !        -huge(0)    when all elements are zero (or n is zero);
  !        otherwise the binary exponent of the largest absolute value.
  !
  !  Input parameters:
  !    n    - Number of elements.
  !    x    - Vector with increment incx.
  !    incx - Increment, positive.
  !
  !  Result:
  !    code - The integer described above.
    integer, intent(in) :: n, incx
    real(wp), intent(in) :: x(*)
    integer :: code
    integer :: i, k
    real(wp) :: largest, magnitude
    logical :: finite

    largest = 0.0_wp
    finite = .true.
    k = 1
    do i = 1, n
      magnitude = abs(x(k))
      if (magnitude > largest) largest = magnitude
      !The negated form also detects NaN.
      if (.not. (magnitude <= huge(1.0_wp))) finite = .false.
      k = k + incx
    end do
    if (.not. finite) then
      code = huge(0)
    else if (largest == 0.0_wp) then
      code = -huge(0)
    else
      code = exponent(largest)
    end if
  end function real_scale_code

  pure subroutine real_reflector(alpha, xnorm2, sc, beta, tau, factor1, &
                                 factor2, factor3)
  !Subroutine real_reflector generates an elementary Householder reflector
  !H = I - tau*v*v^T with v = (1, v2) such that H^T*(alpha, x) = (beta, 0),
  !with the conventions of the LAPACK routine DLARFG: beta has the sign
  !opposite to alpha, tau = (beta-alpha)/beta, v2 = x/(alpha-beta), and
  !H = I (tau = 0) when x is zero.
  !
  !The vector x is distributed and is not an argument. The caller supplies
  !the squared norm of sc*x, where sc is a power of two (one unless the
  !column needed scaling), and receives the factors with which it completes
  !the operation:
  !
  !        v2 = (x*factor1)*factor2 ,
  !        v^T*a(:,c) = a(j,c) + (factor2*z(c))*factor3 ,
  !
  !where z(c) = sum_i sc*x(i)*a(i,c). When alpha and x are both so small
  !that beta is below QR_SAFE_MINIMUM, the reflector is generated from the
  !scaled quantities, as DLARFG does.
  !
  !  Input parameters:
  !    alpha  - The pivot element.
  !    xnorm2 - Squared Euclidean norm of sc*x.
  !    sc     - Scale factor of x.
  !
  !  Output parameters:
  !    beta   - The new diagonal element.
  !    tau    - The scalar factor of the reflector.
  !    factor1, factor2, factor3 - See above.
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: xnorm2, sc
    real(wp), intent(out) :: tau, factor2
    real(wp), intent(out) :: beta, factor1, factor3
    real(wp) :: alpha_scaled
    real(wp) :: beta_scaled, norm, xnorm, xnorm_scaled

    xnorm_scaled = sqrt(xnorm2)
    xnorm = xnorm_scaled / sc
    if (xnorm == 0.0_wp) then
      beta = alpha
      tau = 0.0_wp
      factor1 = 1.0_wp
      factor2 = 0.0_wp
      factor3 = 0.0_wp
      return
    end if
    norm = hypot(abs(alpha), xnorm)
    if (norm >= QR_SAFE_MINIMUM) then
      beta = -sign(norm, alpha)
      tau = (beta - alpha) / beta
      factor1 = 1.0_wp
      factor2 = 1.0_wp / (alpha - beta)
      factor3 = 1.0_wp / sc
    else
      alpha_scaled = alpha * sc
      norm = hypot(abs(alpha_scaled), xnorm_scaled)
      beta_scaled = -sign(norm, alpha_scaled)
      tau = (beta_scaled - alpha_scaled) / beta_scaled
      factor1 = sc
      factor2 = 1.0_wp / (alpha_scaled - beta_scaled)
      factor3 = 1.0_wp
      beta = beta_scaled / sc
    end if
  end subroutine real_reflector

  pure subroutine real_block_product(m, b, ldb, ncols, v, ldv, k, w, ldw, &
                                     work)
  !Subroutine real_block_product adds to w the product of the local rows of a
  !distributed matrix with the local rows of the reflector block V,
  !
  !        w(i,jj) = w(i,jj) + sum_l b(i,l)*v(jj,l) ,
  !
  !i = 1..m, jj = 1..k, l = 1..ncols. Column l of b is a local row of the
  !matrix (transposed layout) and column l of v the same local row of V, so
  !that w, summed over the ranks, is the transpose of V^T*A.
  !
  !The rows i are processed in chunks, such that the part of w in work fits
  !in the cache, and four columns of b at a time, so that a column of w is
  !read and written once per four columns of b.
  !
  !With QR_FACTOR_DOT_KERNELS the product is formed from inner products
  !instead: blocks of QR_FACTOR_DOT_ROWS local rows of b and v are copied
  !transposed to work, chunk by chunk, and real_dot_multiply accumulates the
  !sums over the local rows in registers.
  !
  !  Input parameters:
  !    m     - Number of rows of b and w (columns of the distributed matrix).
  !    b     - Local rows of the matrix; leading dimension ldb.
  !    ldb   - Leading dimension of b.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V: v(jj,l) is the element of reflector jj in
  !            local row l; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    k     - Number of reflectors of the block.
  !    ldw   - Leading dimension of w.
  !
  !  Input/output parameter:
  !    w     - Accumulated product, m by k.
  !
  !  Workspace:
  !    work  - Referenced only with QR_FACTOR_DOT_KERNELS; length as allocated
  !            by real_factor_distributed.
    integer, intent(in) :: m, ldb, ncols, ldv, k, ldw
    real(wp), intent(in) :: b(ldb, *), v(ldv, *)
    real(wp), intent(inout) :: w(ldw, *)
    real(wp), intent(inout) :: work(*)
    integer :: chunk, i, i0, i1, jj, l, l0, l4, nl
    real(wp) :: c1, c2, c3, c4

    if (QR_FACTOR_DOT_KERNELS) then
      !work(1:nl*k) holds the block of v transposed (nl by k), and the rest
      !a chunk of the block of b transposed (nl by at most chunk).
      chunk = max(16, QR_FACTOR_CACHE_BITS / &
                      (QR_FACTOR_DOT_ROWS * storage_size(c1)))
      do l0 = 1, ncols, QR_FACTOR_DOT_ROWS
        nl = min(QR_FACTOR_DOT_ROWS, ncols - l0 + 1)
        do jj = 1, k
          do l = 1, nl
            work(l + (jj - 1) * nl) = v(jj, l0 + l - 1)
          end do
        end do
        do i0 = 1, m, chunk
          i1 = min(i0 + chunk - 1, m)
          do l = 1, nl
            do i = i0, i1
              work(nl * k + l + (i - i0) * nl) = b(i, l0 + l - 1)
            end do
          end do
          call real_dot_multiply(nl, i1 - i0 + 1, k, work(nl * k + 1), nl, &
                                 work(1), nl, 1.0_wp, w(i0, 1), ldw)
        end do
      end do
      return
    end if

    chunk = max(16, QR_FACTOR_CACHE_BITS / (k * storage_size(c1)))
    l4 = ncols - mod(ncols, 4)
    do i0 = 1, m, chunk
      i1 = min(i0 + chunk - 1, m)
      do l = 1, l4, 4
        do jj = 1, k
          c1 = v(jj, l)
          c2 = v(jj, l + 1)
          c3 = v(jj, l + 2)
          c4 = v(jj, l + 3)
          do i = i0, i1
            w(i, jj) = w(i, jj) + c1 * b(i, l) + c2 * b(i, l + 1) + &
                       c3 * b(i, l + 2) + c4 * b(i, l + 3)
          end do
        end do
      end do
      do l = l4 + 1, ncols
        do jj = 1, k
          c1 = v(jj, l)
          do i = i0, i1
            w(i, jj) = w(i, jj) + c1 * b(i, l)
          end do
        end do
      end do
    end do
  end subroutine real_block_product

  pure subroutine real_block_update(m, b, ldb, ncols, y, ldy, k, w, ldw, &
                                    work)
  !Subroutine real_block_update subtracts from the local rows of a distributed
  !matrix the product of the local rows of Y = V*T^T (or V*T) with the
  !replicated matrix W,
  !
  !        b(i,l) = b(i,l) - sum_jj w(i,jj)*y(jj,l) ,
  !
  !i = 1..m, l = 1..ncols, jj = 1..k. The rows i are processed in chunks, as
  !in real_block_product, and four columns of w at a time (eight when
  !EIGHT_COLUMNS is set; four were measured to be slightly faster for real
  !data in double precision).
  !
  !With QR_FACTOR_DOT_KERNELS every chunk of w is copied transposed to work,
  !and real_dot_multiply forms the sums over jj as inner products.
  !
  !  Input parameters:
  !    m     - Number of rows of b and w (columns of the distributed matrix).
  !    ldb   - Leading dimension of b.
  !    ncols - Number of local rows; positive.
  !    y     - Local rows of Y: y(jj,l) is the element of column jj of Y in
  !            local row l; leading dimension ldy.
  !    ldy   - Leading dimension of y.
  !    k     - Number of reflectors of the block.
  !    w     - Replicated m by k matrix; leading dimension ldw.
  !    ldw   - Leading dimension of w.
  !
  !  Input/output parameter:
  !    b     - Local rows of the matrix; leading dimension ldb.
  !
  !  Workspace:
  !    work  - Referenced only with QR_FACTOR_DOT_KERNELS; length as allocated
  !            by real_factor_distributed.
    integer, intent(in) :: m, ldb, ncols, ldy, k, ldw
    real(wp), intent(in) :: y(ldy, *), w(ldw, *)
    real(wp), intent(inout) :: b(ldb, *)
    real(wp), intent(inout) :: work(*)
    logical, parameter :: EIGHT_COLUMNS = .false.
    integer :: chunk, i, i0, i1, jj, k4, k8, l
    real(wp) :: c1, c2, c3, c4, c5, c6, c7, c8

    chunk = max(16, QR_FACTOR_CACHE_BITS / (k * storage_size(c1)))
    if (QR_FACTOR_DOT_KERNELS) then
      do i0 = 1, m, chunk
        i1 = min(i0 + chunk - 1, m)
        do jj = 1, k
          do i = i0, i1
            work(jj + (i - i0) * k) = w(i, jj)
          end do
        end do
        call real_dot_multiply(k, i1 - i0 + 1, ncols, work(1), k, y, ldy, &
                               -1.0_wp, b(i0, 1), ldb)
      end do
      return
    end if

    k8 = 0
    if (EIGHT_COLUMNS) k8 = k - mod(k, 8)
    k4 = k - mod(k - k8, 4)
    do i0 = 1, m, chunk
      i1 = min(i0 + chunk - 1, m)
      do l = 1, ncols
        do jj = 1, k8, 8
          c1 = y(jj, l)
          c2 = y(jj + 1, l)
          c3 = y(jj + 2, l)
          c4 = y(jj + 3, l)
          c5 = y(jj + 4, l)
          c6 = y(jj + 5, l)
          c7 = y(jj + 6, l)
          c8 = y(jj + 7, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj) - c2 * w(i, jj + 1) - &
                      c3 * w(i, jj + 2) - c4 * w(i, jj + 3) - &
                      c5 * w(i, jj + 4) - c6 * w(i, jj + 5) - &
                      c7 * w(i, jj + 6) - c8 * w(i, jj + 7)
          end do
        end do
        do jj = k8 + 1, k4, 4
          c1 = y(jj, l)
          c2 = y(jj + 1, l)
          c3 = y(jj + 2, l)
          c4 = y(jj + 3, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj) - c2 * w(i, jj + 1) - &
                      c3 * w(i, jj + 2) - c4 * w(i, jj + 3)
          end do
        end do
        do jj = k4 + 1, k
          c1 = y(jj, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj)
          end do
        end do
      end do
    end do
  end subroutine real_block_update

  pure subroutine real_dot_multiply(n, ma, mb, a, lda, b, ldb, alpha, c, ldc)
  !Subroutine real_dot_multiply adds a multiple of the product of the
  !transpose of one matrix with another matrix to a third,
  !
  !        c(i,j) = c(i,j) + alpha * sum_l a(l,i)*b(l,j) ,
  !
  !i = 1..ma, j = 1..mb, l = 1..n (no element is conjugated). Every element of
  !the product is an inner product of two contiguous columns; four inner
  !products with the same column of a are accumulated in one pass, so that
  !each element of c is read and written only once.
  !
  !  Input parameters:
  !    n     - Number of rows of a and b.
  !    ma    - Number of columns of a (rows of c).
  !    mb    - Number of columns of b (columns of c).
  !    a     - n by ma matrix; leading dimension lda.
  !    lda   - Leading dimension of a.
  !    b     - n by mb matrix; leading dimension ldb.
  !    ldb   - Leading dimension of b.
  !    alpha - Real multiplier of the product.
  !    ldc   - Leading dimension of c.
  !
  !  Input/output parameter:
  !    c     - ma by mb matrix; leading dimension ldc.
    integer, intent(in) :: n, ma, mb, lda, ldb, ldc
    real(wp), intent(in) :: a(lda, *), b(ldb, *)
    real(wp), intent(in) :: alpha
    real(wp), intent(inout) :: c(ldc, *)
    integer :: i, j, l, mb4
    real(wp) :: t1, t2, t3, t4, x

    mb4 = mb - mod(mb, 4)
    do j = 1, mb4, 4
      do i = 1, ma
        t1 = 0.0_wp
        t2 = 0.0_wp
        t3 = 0.0_wp
        t4 = 0.0_wp
        do l = 1, n
          x = a(l, i)
          t1 = t1 + x * b(l, j)
          t2 = t2 + x * b(l, j + 1)
          t3 = t3 + x * b(l, j + 2)
          t4 = t4 + x * b(l, j + 3)
        end do
        c(i, j) = c(i, j) + alpha * t1
        c(i, j + 1) = c(i, j + 1) + alpha * t2
        c(i, j + 2) = c(i, j + 2) + alpha * t3
        c(i, j + 3) = c(i, j + 3) + alpha * t4
      end do
    end do
    do j = mb4 + 1, mb
      do i = 1, ma
        t1 = 0.0_wp
        do l = 1, n
          t1 = t1 + a(l, i) * b(l, j)
        end do
        c(i, j) = c(i, j) + alpha * t1
      end do
    end do
  end subroutine real_dot_multiply

  pure subroutine real_block_gram(k, ncols, v, ldv, g, ldg)
  !Subroutine real_block_gram adds to g the local part of the strict upper
  !triangle of V^T*V,
  !
  !        g(a,b) = g(a,b) + sum_l v(a,l)*v(b,l) ,   a < b ,
  !
  !from which real_block_tfactor computes the triangular factor T.
  !
  !  Input parameters:
  !    k     - Number of reflectors of the block.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    ldg   - Leading dimension of g.
  !
  !  Input/output parameter:
  !    g     - Accumulated strict upper triangle, k by k.
    integer, intent(in) :: k, ncols, ldv, ldg
    real(wp), intent(in) :: v(ldv, *)
    real(wp), intent(inout) :: g(ldg, *)
    integer :: a, b, l
    real(wp) :: vb

    do l = 1, ncols
      do b = 2, k
        vb = v(b, l)
        do a = 1, b - 1
          g(a, b) = g(a, b) + v(a, l) * vb
        end do
      end do
    end do
  end subroutine real_block_gram

  pure subroutine real_block_tfactor(k, tau, g, ldg, t, ldt)
  !Subroutine real_block_tfactor forms the upper triangular factor T of the
  !block reflector H(1)*...*H(k) = I - V*T*V^T from the scalar factors tau
  !and the strict upper triangle of G = V^T*V, by the recurrence of the
  !LAPACK routine DLARFT (forward direction, columnwise storage):
  !
  !        T(1:b-1,b) = -tau(b)*T(1:b-1,1:b-1)*G(1:b-1,b) ,   T(b,b) = tau(b) .
  !
  !  Input parameters:
  !    k   - Number of reflectors.
  !    tau - Scalar factors.
  !    g   - Strict upper triangle of V^T*V; leading dimension ldg.
  !    ldg - Leading dimension of g.
  !    ldt - Leading dimension of t.
  !
  !  Output parameter:
  !    t   - The k by k factor; the part below the diagonal is set to zero.
    integer, intent(in) :: k, ldg, ldt
    real(wp), intent(in) :: tau(k), g(ldg, *)
    real(wp), intent(out) :: t(ldt, *)
    integer :: a, a2, b
    real(wp) :: tsum

    do b = 1, k
      do a = 1, b - 1
        t(a, b) = -tau(b) * g(a, b)
      end do
      !Multiply by the leading triangle in place; row a needs only the
      !elements a..b-1 of the column, which are still unchanged.
      do a = 1, b - 1
        tsum = 0.0_wp
        do a2 = a, b - 1
          tsum = tsum + t(a, a2) * t(a2, b)
        end do
        t(a, b) = tsum
      end do
      t(b, b) = tau(b)
      do a = b + 1, k
        t(a, b) = 0.0_wp
      end do
    end do
  end subroutine real_block_tfactor

  pure subroutine real_block_vt(k, ncols, v, ldv, t, ldt, adjoint, y, ldy)
  !Subroutine real_block_vt multiplies the local rows of the reflector block
  !V by the triangular factor of the block reflector,
  !
  !        Y = V*T^T   (adjoint true:  y(jj,l) = sum_{a>=jj} t(jj,a)*v(a,l))
  !        Y = V*T     (adjoint false: y(jj,l) = sum_{a<=jj} t(a,jj)*v(a,l)).
  !
  !The first form belongs to the reduction of the matrix (application of the
  !adjoint of the block reflector), the second to the accumulation of Q.
  !
  !  Input parameters:
  !    k     - Number of reflectors of the block.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    t     - Upper triangular factor; leading dimension ldt.
  !    ldt   - Leading dimension of t.
  !    adjoint - Selects the form.
  !    ldy   - Leading dimension of y.
  !
  !  Output parameter:
  !    y     - Local rows of Y, k by ncols.
    integer, intent(in) :: k, ncols, ldv, ldt, ldy
    real(wp), intent(in) :: v(ldv, *), t(ldt, *)
    logical, intent(in) :: adjoint
    real(wp), intent(out) :: y(ldy, *)
    integer :: a, jj, l
    real(wp) :: tsum

    do l = 1, ncols
      if (adjoint) then
        do jj = 1, k
          tsum = 0.0_wp
          do a = jj, k
            tsum = tsum + t(jj, a) * v(a, l)
          end do
          y(jj, l) = tsum
        end do
      else
        do jj = 1, k
          tsum = 0.0_wp
          do a = 1, jj
            tsum = tsum + t(a, jj) * v(a, l)
          end do
          y(jj, l) = tsum
        end do
      end if
    end do
  end subroutine real_block_vt

  subroutine real_begin_active(self, indices, need_flush, info)
  !Subroutine real_begin_active declares the set of functions that the
  !following trials will modify (the ACTIVE set). The module never guesses the
  !active set: every index passed later to replace_symmetric must have been
  !announced here. The routine is collective.
  !
  !The announced indices are either existing functions (index not greater than
  !order()) or the consecutive range order()+1, order()+2, ... of functions
  !that append_symmetric will supply; both kinds may be mixed, in any order.
  !
  !  1. If the announced set equals the current active set (the entries of
  !     the active tail plus the announced appended indices beyond the base
  !     that are still unsupplied), nothing happens.
  !  2. Otherwise the border must be able to hold the request. When it cannot
  !     (pending_count() plus the new entries would exceed max_pending()),
  !     need_flush is returned true and the state is unchanged: the caller
  !     calls factorize_fresh and then begin_active again.
  !  3. The previous tail is committed and the announced indices that already
  !     have a border entry are moved to the end of the border list. An
  !     announced existing function without an entry is necessarily a clean
  !     base slot; it joins D0 with a zero modification: one base solve for
  !     its column of Zb, its base column recovered from the factors as
  !     Q*(R*e), the entries of the other border columns on its row moved
  !     into the border-border block, and one distributed rebuild of the
  !     small matrices.
  !  4. If the committed set changed, LU(Scc) is recomputed. This is the only
  !     place, besides the lazy refresh after truncate changed a committed
  !     entry, where O(m**3) work on the committed block happens.
  !
  !The represented matrix does not change. The tail columns stay as they are
  !until replace_symmetric or append_symmetric supplies new ones. Announcing
  !an appended range reserves border entries; it does not change order().
  !
  !Cost: nothing in case 1; O(m) bookkeeping plus LU(Scc) when only the
  !partition changes; for each newly pending base function one base solve and
  !one column recovery, plus one distributed O(m**2*n0/P) rebuild.
  !
  !  Input parameter:
  !    indices - Distinct canonical indices of the active set. The size may be
  !              zero, which commits the current tail.
  !
  !  Input/output parameter:
  !    self    - A valid real QR state.
  !
  !  Output parameters:
  !    need_flush - True when the request does not fit in the border; the
  !                 state is then unchanged and info is QR_SUCCESS.
  !    info    - QR_SUCCESS;
  !              QR_ERR_INVALID_STATE when self holds no valid factorization;
  !              QR_ERR_INVALID_ARGUMENT when an index is not positive, when
  !              indices are repeated, or when the indices above order() do
  !              not form the range order()+1..order()+k;
  !              QR_ERR_CAPACITY_EXCEEDED when order()+k exceeds the capacity;
  !              QR_ERR_BATCH_FULL when the set itself is larger than
  !              max_pending() and can therefore never be held (the caller
  !              applies such a trial by a fresh factorization);
  !              QR_ERR_SINGULAR when R or the committed block Scc is
  !              numerically singular. The border bookkeeping is then
  !              consistent, but the state cannot solve with this partition;
  !              the caller refactorizes.
  !
  !Argument errors leave the state unchanged.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: indices(:)
    logical, intent(out) :: need_flush
    integer, intent(out) :: info
    integer :: count_active, count_append, count_beyond, count_entry
    integer :: count_new, i, idx, k, mc, p
    logical :: committed_changed, same_set

    need_flush = .false.
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return

    !Validate the indices and classify them.
    count_active = size(indices)
    count_append = 0
    do i = 1, count_active
      if (indices(i) < 1) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      if (indices(i) > self%n) count_append = count_append + 1
    end do
    do i = 1, count_active
      if (indices(i) > self%n + count_append) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
    end do
    if (self%n + count_append > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    do i = 1, count_active
      idx = indices(i)
      if (self%imark(idx) /= 0) then
        !Repeated index. Remove the marks set so far.
        do k = 1, i - 1
          self%imark(indices(k)) = 0
        end do
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      self%imark(idx) = 1
    end do
    !From here on every return path must clear the marks.

    if (count_active > self%mmax) then
      call clear_marks()
      info = QR_ERR_BATCH_FULL
      return
    end if

    !Compare with the current active set. An announced index is beyond the
    !base (pure reservation), has a border entry (a pending function or a
    !deleted base slot that an append will reuse), or is a clean base slot.
    mc = self%m - self%ma
    count_new = 0
    count_beyond = 0
    count_entry = 0
    same_set = .true.
    do i = 1, count_active
      idx = indices(i)
      if (idx > self%n .and. idx > self%n0) then
        count_beyond = count_beyond + 1
      else if (self%bpos(idx) == 0) then
        count_new = count_new + 1
        same_set = .false.
      else
        count_entry = count_entry + 1
        if (self%bpos(idx) <= mc) same_set = .false.
      end if
    end do
    if (count_entry /= self%ma .or. count_beyond /= self%nres) same_set = .false.
    if (same_set) then
      call clear_marks()
      info = QR_SUCCESS
      return
    end if

    if (self%m + count_new + count_beyond > self%mmax) then
      call clear_marks()
      need_flush = .true.
      info = QR_SUCCESS
      return
    end if
    if (count_new > 0 .and. self%r_singular) then
      call clear_marks()
      info = QR_ERR_SINGULAR
      return
    end if

    !Commit the previous tail and move the announced entries to the end of
    !the border list, keeping the relative order of both groups.
    committed_changed = .false.
    k = 0
    do p = 1, self%m
      if (self%imark(self%bidx(p)) == 0) then
        k = k + 1
        self%iw(k, 1) = self%bidx(p)
        self%iw(k, 2) = self%bslot(p)
        if (p > mc) committed_changed = .true.
      end if
    end do
    self%ma = self%m - k
    do p = 1, self%m
      if (self%imark(self%bidx(p)) /= 0) then
        k = k + 1
        self%iw(k, 1) = self%bidx(p)
        self%iw(k, 2) = self%bslot(p)
        if (p <= mc) committed_changed = .true.
      end if
    end do
    do p = 1, self%m
      self%bidx(p) = self%iw(p, 1)
      self%bslot(p) = self%iw(p, 2)
      self%bpos(self%bidx(p)) = p
    end do
    self%nres = count_beyond
    self%tail_ok = .false.
    call clear_marks()

    !Clean base slots join D0 with a zero modification and enter the tail.
    if (count_new > 0) then
      do i = 1, count_active
        idx = indices(i)
        if (idx > self%n) cycle
        if (self%bpos(idx) /= 0) cycle
        call real_join_base_slot(self, idx, .false.)
        self%ma = self%ma + 1
      end do
      call real_rebuild_small(self)
      committed_changed = .true.
    end if

    if (committed_changed) self%scc_ok = .false.
    call real_refresh_committed(self)
    if (self%r_singular .or. self%scc_singular) then
      info = QR_ERR_SINGULAR
    else
      info = QR_SUCCESS
    end if

  contains

    subroutine clear_marks()
      integer :: mark
      do mark = 1, count_active
        self%imark(indices(mark)) = 0
      end do
    end subroutine clear_marks

  end subroutine real_begin_active

  subroutine real_replace_symmetric(self, idx, h_col, s_col, info)
  !Subroutine real_replace_symmetric supplies the NEW full conceptual column
  !of the active function idx. In contrast to the previous implementation the
  !arguments are the new columns themselves, not differences. The represented
  !matrix becomes the symmetric matrix whose row and column idx equal
  !
  !                 a = h_col - shift*s_col ,
  !
  !all other entries being unchanged. Nothing in Q or R changes. The routine
  !is collective.
  !
  !The column is split into its part on the clean rows (stored in Bt, as the
  !difference from the base column when idx is a base slot, and zero on the
  !rows of D0) and its entries on the border indices (stored in row and
  !column of the border-border block). One distributed base solve gives
  !Gs(:,k). The row and column of BG, the row of BZ, and the row and column of
  !the Schur complement that belong to the entry are updated with
  !O((2*m+m0)*n0) local work. The tail factors X and LU(Sigma) are recomputed
  !lazily by the next solve.
  !
  !When several active functions change in one trial, the caller supplies the
  !final column of each of them in any order; the entry shared by two active
  !functions is taken from the column supplied last, so the columns must be
  !mutually consistent (symmetric) at the end of the trial.
  !
  !  Input parameters:
  !    idx   - Canonical index of the function, 1 <= idx <= order(). It must
  !            belong to the active set announced by begin_active.
  !    h_col - New column idx of H, length order().
  !    s_col - New column idx of S, length order().
  !
  !  Input/output parameter:
  !    self  - A valid real QR state. On success both structural-update
  !            counters are incremented.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS;
  !            QR_ERR_INVALID_STATE when self holds no valid factorization;
  !            QR_ERR_INVALID_ARGUMENT when idx is out of range or is not an
  !            active border entry;
  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from
  !            order();
  !            QR_ERR_SINGULAR when R is numerically singular; the caller
  !            must refactorize.
  !
  !No allocation is performed. An error leaves the state unchanged.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(in) :: h_col(:), s_col(:)
    integer, intent(out) :: info
    integer :: i, p

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    if (idx < 1 .or. idx > self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(h_col) /= self%n .or. size(s_col) /= self%n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    p = self%bpos(idx)
    if (p <= self%m - self%ma) then
      !Not pending, or committed: the function was not announced.
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    do i = 1, self%n
      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)
    end do
    call real_install_column(self, p, self%wn(:, 3), .false.)

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_replace_symmetric

  subroutine real_append_symmetric(self, h_col, s_col, info)
  !Subroutine real_append_symmetric adds the function order()+1 with its full
  !conceptual column. The order grows by one and the new function becomes an
  !entry of the active tail. Nothing in Q or R changes. The routine is
  !collective.
  !
  !The new index normally was announced by begin_active as part of the
  !appended range. An unannounced append is also accepted when the border can
  !hold it and simply extends the active tail. Two cases exist:
  !
  !  - the index lies beyond the base (the usual case): a new border entry is
  !    created (an announced reservation is consumed);
  !  - the index is a trailing base slot that truncate deleted earlier: the
  !    border entry of the deleted slot receives the new column. If that
  !    entry was committed it moves to the tail, which makes the next solve
  !    refactorize the committed block.
  !
  !Both cases cost one distributed base solve and O((2*m+m0)*n0) local work.
  !
  !  Input parameters:
  !    h_col - Column order()+1 of H, length order()+1 (the last element is
  !            the new diagonal element).
  !    s_col - Column order()+1 of S, length order()+1.
  !
  !  Input/output parameter:
  !    self  - A valid real QR state. On success order() grows by one and both
  !            structural-update counters are incremented.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS;
  !            QR_ERR_INVALID_STATE when self holds no valid factorization;
  !            QR_ERR_CAPACITY_EXCEEDED when order() equals the capacity;
  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from
  !            order()+1;
  !            QR_ERR_BATCH_FULL when the index lies beyond the base, was not
  !            announced, and the border is full; the caller must
  !            refactorize;
  !            QR_ERR_SINGULAR when R is numerically singular.
  !
  !No allocation is performed. An error leaves the state unchanged.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h_col(:), s_col(:)
    integer, intent(out) :: info
    integer :: i, new_n, p, p2, slot, slot2

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    new_n = self%n + 1
    if (new_n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    if (size(h_col) /= new_n .or. size(s_col) /= new_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (new_n <= self%n0) then
      !A deleted base slot always owns a border entry.
      if (self%bpos(new_n) == 0) return
    else if (self%nres == 0 .and. self%m >= self%mmax) then
      info = QR_ERR_BATCH_FULL
      return
    end if
    if (self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    if (new_n <= self%n0) then
      !Reuse the entry of the deleted slot; move it to the end of the tail
      !if it is committed.
      p = self%bpos(new_n)
      if (p <= self%m - self%ma) then
        slot = self%bslot(p)
        do p2 = p, self%m - 1
          self%bidx(p2) = self%bidx(p2 + 1)
          self%bslot(p2) = self%bslot(p2 + 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        p = self%m
        self%bidx(p) = new_n
        self%bslot(p) = slot
        self%bpos(new_n) = p
        self%ma = self%ma + 1
        self%scc_ok = .false.
      end if
    else
      !Create the border entry at the end of the tail.
      if (self%nres > 0) self%nres = self%nres - 1
      p = self%m + 1
      slot = self%bslot(p)
      self%bidx(p) = new_n
      self%bpos(new_n) = p
      self%m = p
      self%ma = self%ma + 1
      self%stale(slot) = 0
      do p2 = 1, self%m
        slot2 = self%bslot(p2)
        self%cbb(slot, slot2) = 0.0_wp
        self%cbb(slot2, slot) = 0.0_wp
      end do
    end if
    self%n = new_n

    do i = 1, new_n
      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)
    end do
    call real_install_column(self, p, self%wn(:, 3), .false.)

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_append_symmetric

  subroutine real_truncate(self, target_order, info)
  !Subroutine real_truncate removes the trailing functions
  !target_order+1..order(). The routine is collective.
  !
  !  - A removed function beyond the base loses its border entry (O(m) work;
  !    if the entry was committed, LU(Scc) is recomputed lazily).
  !  - A removed pending base function keeps its border entry, which is
  !    marked as a deleted slot. Its decoupling column is installed lazily
  !    (one base solve) by the next solve, unless append_symmetric reuses the
  !    index first, in which case no work at all is spent on the deletion.
  !  - A removed clean base function joins D0 as a deleted slot: one base
  !    solve, one column recovery, an O(m*n0) correction of the stored border
  !    columns, and one distributed rebuild of the small matrices after all
  !    removals. Its entry is committed.
  !
  !A base slot stays in D0 until the next fresh factorization. Announced but
  !unsupplied appended indices beyond the base are forgotten when the order
  !actually decreases.
  !
  !Typical use (basis enlargement): truncate(base) drops the previous
  !candidate, whose border entry is the active tail, without touching the
  !committed set or LU(Scc).
  !
  !  Input parameter:
  !    target_order - New order, 1 <= target_order <= order(). A value equal
  !                   to order() is a no-op.
  !
  !  Input/output parameter:
  !    self - A valid real QR state. On success the structural-update counters
  !           grow by the number of removed functions.
  !
  !  Output parameter:
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self holds no valid factorization;
  !           QR_ERR_INVALID_ARGUMENT when target_order is out of range;
  !           QR_ERR_BATCH_FULL when the border cannot take the clean base
  !           slots that would be deleted; the state is unchanged and the
  !           caller refactorizes at the target order;
  !           QR_ERR_SINGULAR when R is numerically singular and a base solve
  !           would be required (state unchanged).
  !
  !No allocation is performed.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: target_order
    integer, intent(out) :: info
    integer :: count_new, i, mc, old_n, p, p2, slot
    logical :: grew

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    if (target_order < 1 .or. target_order > self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (target_order == self%n) then
      info = QR_SUCCESS
      return
    end if

    !Count the clean base slots that need a border entry before changing
    !anything.
    count_new = 0
    do i = target_order + 1, self%n
      if (self%bpos(i) == 0) count_new = count_new + 1
    end do
    if (self%m + count_new > self%mmax) then
      info = QR_ERR_BATCH_FULL
      return
    end if
    if (count_new > 0 .and. self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    old_n = self%n
    grew = .false.
    do i = old_n, target_order + 1, -1
      p = self%bpos(i)
      mc = self%m - self%ma
      if (i > self%n0) then
        !Appended function: drop the border entry and free its slot.
        if (p <= mc) then
          self%scc_ok = .false.
        else
          self%ma = self%ma - 1
          self%tail_ok = .false.
        end if
        slot = self%bslot(p)
        do p2 = p, self%m - 1
          self%bidx(p2) = self%bidx(p2 + 1)
          self%bslot(p2) = self%bslot(p2 + 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        self%bslot(self%m) = slot
        self%bidx(self%m) = 0
        self%bpos(i) = 0
        self%m = self%m - 1
      else if (p > 0) then
        !Pending base function: its entry becomes a deleted slot.
        self%stale(self%bslot(p)) = 1
        if (p <= mc) then
          self%scc_ok = .false.
        else
          self%tail_ok = .false.
        end if
      else
        !Clean base function: it joins D0 as a deleted slot. The new entry
        !is created at the end of the list and moved in front of the tail.
        call real_join_base_slot(self, i, .true.)
        p = self%m
        slot = self%bslot(p)
        do p2 = p, mc + 2, -1
          self%bidx(p2) = self%bidx(p2 - 1)
          self%bslot(p2) = self%bslot(p2 - 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        self%bidx(mc + 1) = i
        self%bslot(mc + 1) = slot
        self%bpos(i) = mc + 1
        self%scc_ok = .false.
        grew = .true.
      end if
    end do
    self%n = target_order
    self%nres = 0
    self%structural_updates = self%structural_updates + &
                              int(old_n - target_order, int64)
    self%updates_since_fresh = self%updates_since_fresh + &
                               int(old_n - target_order, int64)
    if (grew) call real_rebuild_small(self)
    info = QR_SUCCESS
  end subroutine real_truncate

  subroutine real_delete_symmetric(self, idx, info)
  !Subroutine real_delete_symmetric removes the LAST function. It is kept for
  !interface compatibility and is equivalent to truncate(order()-1). Interior
  !deletion is not supported by the deferred scheme: the caller compacts the
  !physical matrices and refactorizes.
  !
  !  Input parameter:
  !    idx  - Index of the function to delete; must equal order(), and the
  !           order must be at least two.
  !
  !  Input/output parameter:
  !    self - A valid real QR state.
  !
  !  Output parameter:
  !    info - QR_ERR_INVALID_ARGUMENT when idx differs from order() or the
  !           order is one; otherwise the status of truncate.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (self%n < 2 .or. idx /= self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    call real_truncate(self, self%n - 1, info)
  end subroutine real_delete_symmetric

  subroutine real_solve(self, h, s, v_initial, x, lambda, tol, max_iter, &
                        norm_mode, rel_acc, num_iter, eig_residual, info, &
                        residual_goal, num_polish)
  !Subroutine real_solve finds one eigenvalue and its eigenvector for the real
  !generalized symmetric eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration with the represented matrix H-shift*S (base
  !factors plus border). H and S are the CURRENT physical matrices; the caller
  !keeps them consistent with the columns supplied to replace_symmetric and
  !append_symmetric. The routine is collective.
  !
  !The desired eigenvalue should be closer to shift than every other
  !eigenvalue. The convergence rate is governed principally by the ratio of
  !the distance from shift to the desired eigenvalue and the distance from
  !shift to the next closest eigenvalue. A nondegenerate eigenvalue is
  !assumed.
  !
  !For a current vector v, one iteration computes the solution of
  !
  !             (H-shift*S)*x = S*v ,
  !
  !up to a scalar factor. The first iteration solves this system directly:
  !the right-hand side is a distributed symmetric product with the lower
  !triangle of S and the solve is the bordered solve described in the module
  !header (one distributed base solve). Every later iteration obtains the
  !same direction in residual-correction form,
  !
  !             mu = (v^T*H*v)/(v^T*S*v),
  !             x  = v - (H-shift*S)^{-1} * (H*v - mu*S*v)
  !                = (mu-shift) * (H-shift*S)^{-1} * S*v ,
  !
  !at the cost of one additional symmetric product. The reason is numerical.
  !The bordered solve eliminates on the clean block of the base matrix, whose
  !inverse is represented as a difference of large terms when the base matrix
  !is nearly singular; its error is proportional to the size of the solution
  !and can exceed that of a solve with fresh factors by many orders of
  !magnitude when the pending functions carry much of the eigenvector. In the
  !correction form the solve is applied to the small residual of the current
  !iterate, so this error only slows the convergence slightly and does not
  !limit the attainable accuracy; the fixed point is the eigenvector of the
  !physical pair (H,S). The new vector is divided by max(abs(x)). Its change
  !of direction relative to v is measured by
  !
  !             alpha   = (x^T*v)/(v^T*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0, iteration stops when rel_acc<=tol. For tol<=0, iteration
  !continues until rel_acc no longer decreases (it begins to increase, or the
  !iterate has reached a fixed point or a rounding cycle of the correction
  !form, where rel_acc stays zero or constant) and the current value is not
  !larger than abs(tol). The latter rule requests the most
  !accurate attainable direction subject to the floor abs(tol), and
  !necessarily performs at least one additional comparison iteration.
  !
  !After iteration, two distributed products hs = H*x and ss = S*x give the
  !physical Rayleigh quotient and the physical eigenpair residual
  !
  !             lambda       = (x^T*hs)/(x^T*ss),
  !             eig_residual = ||hs-lambda*ss||_2 /
  !                            (||hs||_2+abs(lambda)*||ss||_2+tiny(1.0_wp)) .
  !
  !The residual is independent of the factorization and of the border data,
  !so it also detects a damaged representation; the caller decides whether to
  !refactorize. Only the lower triangles of the leading n by n blocks of H and
  !S are referenced.
  !
  !Optional additional correction steps ("polish"). When residual_goal is
  !present and eig_residual exceeds it after the iteration has stopped (for
  !example because the stopping rule was satisfied after the first, directly
  !solved iteration, or because max_iter was reached), the routine applies up
  !to QR_MAX_POLISH further steps
  !
  !             x <- x - (H-shift*S)^{-1} * (H*x - lambda*S*x) ,
  !
  !reusing the two products of the residual evaluation. Each step costs one
  !bordered solve and two symmetric products and is accepted only if it
  !lowers the physical residual; the steps stop when the goal is met, when a
  !step fails to lower the residual, or when it lowers it by less than a
  !factor of two. This gives the caller a cheap alternative to a fresh
  !factorization when the residual check fails. Without residual_goal no such
  !step is made. The stopping rule, rel_acc, num_iter and info are never
  !affected by these steps.
  !
  !  Input parameters:
  !    h         - Real array with at least n rows and n columns holding the
  !                lower triangle of the current H.
  !    s         - Real array with at least n rows and n columns holding the
  !                lower triangle of the current S, expected to be positive
  !                definite.
  !    v_initial - A nonzero starting approximation of length n. The array is
  !                not modified.
  !    tol       - Directional convergence tolerance. Its sign selects the
  !                stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^T*S*x = 1;
  !                  1: x^T*x = 1;
  !                other: max(abs(x)) = 1.
  !
  !  Input/output parameter:
  !    self      - A valid real QR state of order n. The represented matrix
  !                and the counters are unchanged; the lazily maintained tail
  !                factors and the workspace are updated.
  !
  !  Output parameters:
  !    x         - The final eigenvector approximation in the normalization
  !                selected by norm_mode. It has length n.
  !    lambda    - The Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - The directional difference from the final iteration. This
  !                is a convergence estimate, not a rigorously bounded error.
  !    num_iter  - Number of inverse iterations performed.
  !    eig_residual - Relative physical eigenpair residual defined above.
  !    info      - QR_SUCCESS when the stopping criterion is satisfied;
  !                QR_ERR_INVALID_STATE when self has no valid factorization;
  !                QR_ERR_DIMENSION_MISMATCH when H, S, v_initial, or x has an
  !                extent inconsistent with the order;
  !                QR_ERR_INVALID_ARGUMENT when max_iter is not positive;
  !                QR_ERR_ZERO_INITIAL_VECTOR when v_initial is numerically
  !                zero in the working precision;
  !                QR_ERR_NONPOSITIVE_OVERLAP when the final x^T*S*x is not
  !                greater than the precision-scaled zero threshold;
  !                QR_ERR_SINGULAR when a diagonal element of R is zero,
  !                subnormal, or not finite, when a small matrix of the border
  !                (Scc or Sigma) cannot be factorized, or when an
  !                iteration produces a numerically zero or non-finite vector;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. In this
  !                case x, lambda, rel_acc, num_iter and eig_residual describe
  !                the best approximation reached.
  !
  !  Optional parameters:
  !    residual_goal - Input. Requested bound of eig_residual; enables the
  !                    residual-correction steps described above. A value of
  !                    zero requests the best attainable residual.
  !    num_polish    - Output. Number of accepted correction steps (zero when
  !                    residual_goal is absent).
  !
  !No allocation is performed. On an error detected before the first
  !iteration, x and lambda are zero, num_iter is zero, and rel_acc and
  !eig_residual are huge.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: v_initial(:)
    real(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter
    real(wp), intent(out) :: eig_residual
    integer, intent(out) :: info
    real(wp), intent(in), optional :: residual_goal
    integer, intent(out), optional :: num_polish
    integer :: i, matrix_n, polish_count
    real(wp) :: goal

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    eig_residual = huge(1.0_wp)
    if (present(num_polish)) num_polish = 0
    info = QR_ERR_INVALID_STATE

    !Verify the state, array dimensions, iteration limit, and starting-vector
    !norm before any numerical work.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v_initial) /= matrix_n .or. size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (max_iter <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    do i = 1, matrix_n
      self%wn(i, 3) = v_initial(i)
    end do
    if (real_norm_squared(matrix_n, self%wn(:, 3)) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    !Bring the lazily maintained small factors up to date and refuse to
    !divide by unusable pivots.
    call real_refresh_small(self)
    if (.not. real_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !A negative goal disables the correction steps.
    goal = -1.0_wp
    if (present(residual_goal)) goal = max(residual_goal, 0.0_wp)
    call real_inverse_iteration(self, matrix_n, h, size(h, 1), s, size(s, 1), &
                                self%wn(:, 3), self%wn(:, 4), self%wn(:, 5), &
                                self%wn(:, 6), lambda, tol, max_iter, &
                                norm_mode, goal, rel_acc, num_iter, &
                                eig_residual, polish_count, info)
    if (present(num_polish)) num_polish = polish_count
    do i = 1, matrix_n
      x(i) = self%wn(i, 5)
    end do
  end subroutine real_solve

  subroutine real_inverse_iteration(self, n, h, ldh, s, lds, v, w, x, y, &
                                    lambda, tol, max_iter, norm_mode, &
                                    residual_goal, rel_acc, num_iter, &
                                    eig_residual, num_polish, info)
  !Subroutine real_inverse_iteration is the numerical body of real_solve. The
  !work vectors are passed as explicit-shape arguments so that the loops
  !operate on plain arrays.
  !
  !  Input parameters:
  !    n, h, ldh, s, lds - Order and the two physical matrices with their
  !                        leading dimensions (lower triangles referenced).
  !    tol, max_iter, norm_mode - As for real_solve.
  !    residual_goal - Goal of the residual-correction steps; a negative value
  !                    disables them.
  !
  !  Input/output parameters:
  !    self - Valid state whose small factors are up to date.
  !    v    - On entry the starting vector; overwritten (finally by H*x).
  !    w    - Work vector (S*v during iteration, finally S*x).
  !    y    - Work vector (H*v and the residual right-hand side).
  !
  !  Output parameters:
  !    x, lambda, rel_acc, num_iter, eig_residual, info - As for real_solve;
  !    x is zero when an error prevents a usable iterate.
  !    num_polish - Number of accepted correction steps.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: n, ldh, lds
    real(wp), intent(in) :: h(ldh, *), s(lds, *)
    real(wp), intent(inout) :: v(n), w(n), y(n)
    real(wp), intent(out) :: x(n)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(in) :: residual_goal
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter
    real(wp), intent(out) :: eig_residual
    integer, intent(out) :: num_polish, info
    integer :: i
    real(wp) :: coefficient, current_norm_squared, eigenvector_norm_squared
    real(wp) :: max_component, new_lambda, new_overlap, new_residual
    real(wp) :: norm_of_diff, norm_of_diff_previous, overlap_norm_squared
    real(wp) :: previous_residual
    logical :: corrected, not_converged

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    num_polish = 0
    eig_residual = huge(1.0_wp)
    norm_of_diff = huge(1.0_wp)
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform inverse iterations until the direction criterion is satisfied or
    !the caller-supplied iteration limit is exhausted.
    do while (not_converged .and. num_iter < max_iter)
      !Form w=S*v from the lower triangle of S.
      call real_sym_product(self%rank, self%nprocs, self%comm, n, s, lds, &
                            v, w)
      corrected = .false.
      if (num_iter > 0) then
        !Residual-correction form of the same iterate. v is the previous
        !iterate, an approximate eigenvector; with its Rayleigh quotient mu,
        !    v - A^{-1}*(H*v - mu*S*v) = (mu-shift)*A^{-1}*S*v .
        !The solve is applied to a small right-hand side, so the errors of
        !the bordered solve are proportional to the correction and not to v.
        overlap_norm_squared = real_dot(n, v, w)
        if (overlap_norm_squared > tiny(1.0_wp)) then
          call real_sym_product(self%rank, self%nprocs, self%comm, n, h, ldh, &
                                v, y)
          coefficient = real_dot(n, v, y) / overlap_norm_squared
          do i = 1, n
            y(i) = y(i) - coefficient * w(i)
          end do
          call real_bordered_solve(self, y, x)
          do i = 1, n
            x(i) = v(i) - x(i)
          end do
          corrected = .true.
        end if
      end if
      !First iteration (or unusable quotient): solve (H-shift*S)*x=S*v.
      if (.not. corrected) call real_bordered_solve(self, w, x)

      !Scale the solution so its largest absolute component is one. Scaling
      !does not change the eigendirection and limits growth in repeated solves.
      max_component = 0.0_wp
      do i = 1, n
        max_component = max(max_component, abs(x(i)))
      end do
      if (.not. (max_component > tiny(1.0_wp) .and. &
                 max_component <= huge(1.0_wp))) then
        x = 0.0_wp
        info = QR_ERR_SINGULAR
        return
      end if
      do i = 1, n
        x(i) = x(i) / max_component
      end do

      !Compute the least-direction-change coefficient and the relative norm of
      !the component of x not parallel to the previous iterate v.
      current_norm_squared = real_norm_squared(n, v)
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = real_dot(n, x, v) / current_norm_squared
      eigenvector_norm_squared = real_norm_squared(n, x)
      norm_of_diff = real_direction_difference(n, x, coefficient, v, &
                                               eigenvector_norm_squared)
      !A partially non-finite iterate is detected through its direction error.
      if (.not. (norm_of_diff <= huge(1.0_wp))) then
        x = 0.0_wp
        info = QR_ERR_SINGULAR
        return
      end if

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        !For a non-positive tolerance, accept only after the direction error
        !has passed its minimum and remains within abs(tol). In the
        !residual-correction form the iterate can reach a fixed point or a
        !rounding cycle, where the direction error is zero or constant and
        !never increases; a value that no longer decreases is therefore the
        !minimum as well.
        if (norm_of_diff >= norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        do i = 1, n
          v(i) = x(i)
        end do
      end if
    end do

    rel_acc = norm_of_diff
    !Failure to meet the stopping rule is nonfatal: the final iterate is still
    !used to calculate and normalize an eigenpair approximation.
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    !Physical products for the Rayleigh quotient, the residual, and, for
    !norm_mode=0, the final normalization. A non-positive x^T*S*x violates the
    !positive-definite S precondition or indicates unusable numerical data.
    call real_physical_quotient(self%rank, self%nprocs, self%comm, n, h, ldh, &
                                s, lds, x, v, w, lambda, &
                                overlap_norm_squared, eig_residual)
    if (.not. (overlap_norm_squared > tiny(1.0_wp))) then
      lambda = 0.0_wp
      eig_residual = huge(1.0_wp)
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if

    !Optional residual-correction steps: y = x - A^{-1}*(H*x - lambda*S*x).
    !On entry to every step v = H*x and w = S*x for the current x.
    if (residual_goal >= 0.0_wp) then
      do while (eig_residual > residual_goal .and. num_polish < QR_MAX_POLISH)
        do i = 1, n
          v(i) = v(i) - lambda * w(i)
        end do
        call real_bordered_solve(self, v, y)
        max_component = 0.0_wp
        do i = 1, n
          y(i) = x(i) - y(i)
          max_component = max(max_component, abs(y(i)))
        end do
        if (.not. (max_component > tiny(1.0_wp) .and. &
                   max_component <= huge(1.0_wp))) exit
        do i = 1, n
          y(i) = y(i) / max_component
        end do
        call real_physical_quotient(self%rank, self%nprocs, self%comm, n, h, &
                                    ldh, s, lds, y, v, w, new_lambda, &
                                    new_overlap, new_residual)
        !Reject a step that does not lower the physical residual (the negated
        !comparison also rejects NaN) and keep the current x.
        if (.not. (new_overlap > tiny(1.0_wp))) exit
        if (.not. (new_residual < eig_residual)) exit
        do i = 1, n
          x(i) = y(i)
        end do
        previous_residual = eig_residual
        lambda = new_lambda
        overlap_norm_squared = new_overlap
        eig_residual = new_residual
        num_polish = num_polish + 1
        !Stop when the steps no longer make substantial progress.
        if (eig_residual > 0.5_wp * previous_residual) exit
      end do
    end if

    !Apply the normalization requested by the caller. No scaling is required
    !for other norm_mode values because each iterate already has unit
    !largest-component magnitude.
    select case (norm_mode)
    case (0) !Normalize so that x^T*S*x=1.
      coefficient = sqrt(overlap_norm_squared)
      do i = 1, n
        x(i) = x(i) / coefficient
      end do
    case (1) !Normalize so that x^T*x=1.
      eigenvector_norm_squared = real_norm_squared(n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = sqrt(eigenvector_norm_squared)
      do i = 1, n
        x(i) = x(i) / coefficient
      end do
    end select
  end subroutine real_inverse_iteration

  subroutine real_physical_quotient(rank, nprocs, comm, n, h, ldh, s, lds, &
                                    x, hx, sx, lambda, overlap, residual)
  !Subroutine real_physical_quotient evaluates, with two distributed symmetric
  !products, the physical Rayleigh quotient and the relative eigenpair
  !residual of a vector:
  !
  !        hx = H*x,  sx = S*x,  overlap = x^T*sx,  lambda = (x^T*hx)/overlap,
  !        residual = ||hx-lambda*sx||_2 /
  !                   (||hx||_2+abs(lambda)*||sx||_2+tiny(1.0_wp)) .
  !
  !  Input parameters:
  !    rank, nprocs, comm - Own rank, number of ranks, and communicator.
  !    n, h, ldh, s, lds - Order and the physical matrices (lower triangles).
  !    x    - Replicated vector of length n.
  !
  !  Output parameters:
  !    hx, sx   - The two products.
  !    lambda   - Rayleigh quotient; zero when overlap is not positive.
  !    overlap  - x^T*S*x.
  !    residual - Relative residual; huge when overlap is not positive.
    integer, intent(in) :: rank, nprocs, comm, n, ldh, lds
    real(wp), intent(in) :: h(ldh, *), s(lds, *), x(n)
    real(wp), intent(out) :: hx(n), sx(n)
    real(wp), intent(out) :: lambda, overlap, residual
    integer :: i
    real(wp) :: difference, residual_norm_squared

    call real_sym_product(rank, nprocs, comm, n, h, ldh, x, hx)
    call real_sym_product(rank, nprocs, comm, n, s, lds, x, sx)
    overlap = real_dot(n, x, sx)
    if (.not. (overlap > tiny(1.0_wp))) then
      lambda = 0.0_wp
      residual = huge(1.0_wp)
      return
    end if
    lambda = real_dot(n, x, hx) / overlap
    residual_norm_squared = 0.0_wp
    do i = 1, n
      difference = hx(i) - lambda * sx(i)
      residual_norm_squared = residual_norm_squared + difference * difference
    end do
    residual = sqrt(residual_norm_squared) / &
               (sqrt(real_norm_squared(n, hx)) + &
                abs(lambda) * sqrt(real_norm_squared(n, sx)) + tiny(1.0_wp))
  end subroutine real_physical_quotient

  subroutine real_solve_shifted(self, b, x, info)
  !Subroutine real_solve_shifted solves one linear system
  !
  !                     (H-shift*S)*x = b
  !
  !with the represented matrix of the current order (base factors plus
  !border). It exposes the bordered solve used by inverse iteration for tests
  !and diagnostics. The routine is collective; b must be replicated.
  !
  !  Input parameter:
  !    b    - Right-hand side of length order().
  !
  !  Input/output parameter:
  !    self - A valid real QR state. The represented matrix is unchanged; the
  !           lazily maintained tail factors and the workspace are updated.
  !
  !  Output parameters:
  !    x    - Solution of length order(); zero on error.
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self holds no valid factorization;
  !           QR_ERR_DIMENSION_MISMATCH when b or x does not have length
  !           order();
  !           QR_ERR_SINGULAR when R or a small matrix of the border cannot
  !           be used.
  !
  !No allocation is performed.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: b(:)
    real(wp), intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: i

    x = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    if (size(b) /= self%n .or. size(x) /= self%n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    call real_refresh_small(self)
    if (.not. real_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if
    do i = 1, self%n
      self%wn(i, 3) = b(i)
    end do
    call real_bordered_solve(self, self%wn(:, 3), self%wn(:, 4))
    do i = 1, self%n
      x(i) = self%wn(i, 4)
    end do
    info = QR_SUCCESS
  end subroutine real_solve_shifted

  subroutine real_symmetric_product(self, a, x, y, info)
  !Subroutine real_symmetric_product computes y = A*x for a replicated real
  !symmetric matrix of which only the lower triangle is referenced, with the
  !distributed kernel used by inverse iteration (cyclic column ownership, one
  !reduction; no MPI call with one process). It is offered to callers that
  !need products with the physical matrices, for example for their own
  !residual checks, and to the timing driver. The routine is collective.
  !
  !The order of the product is size(x); it is independent of order() and a
  !factorization need not exist, but the state must be initialized and the
  !order must not exceed its capacity.
  !
  !  Input parameters:
  !    a    - Real array with at least size(x) rows and columns holding the
  !           lower triangle of the matrix.
  !    x    - Replicated vector.
  !
  !  Input/output parameter:
  !    self - An initialized real QR state; two work vectors are overwritten.
  !
  !  Output parameters:
  !    y    - Replicated product, same length as x; zero on error.
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self is not initialized;
  !           QR_ERR_DIMENSION_MISMATCH when x is empty, y has a different
  !           length, or a is too small;
  !           QR_ERR_CAPACITY_EXCEEDED when size(x) exceeds the capacity.
  !
  !No allocation is performed.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: a(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: y(:)
    integer, intent(out) :: info
    integer :: i, n

    y = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%wn)) return
    n = size(x)
    if (n <= 0 .or. size(y) /= n .or. size(a, 1) < n .or. size(a, 2) < n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    do i = 1, n
      self%wn(i, 3) = x(i)
    end do
    call real_sym_product(self%rank, self%nprocs, self%comm, n, a, size(a, 1), &
                          self%wn(:, 3), self%wn(:, 4))
    do i = 1, n
      y(i) = self%wn(i, 4)
    end do
    info = QR_SUCCESS
  end subroutine real_symmetric_product

  subroutine real_shifted_quotient(self, h, s, x, lambda, info)
  !Subroutine real_shifted_quotient evaluates the Rayleigh quotient of a vector
  !for the physical pair (H,S) in the shifted form
  !
  !        lambda = shift + (x^T*(H-shift*S)*x)/(x^T*S*x) ,
  !
  !where every element of H-shift*S is formed before it is multiplied by the
  !components of x. For an eigenvector whose eigenvalue is close to the shift
  !this is far more accurate than the plain quotient (x^T*H*x)/(x^T*S*x) that
  !solve returns: the rounding errors of the sums act on the small difference
  !lambda-shift only, and H and shift*S cancel element by element before the
  !products with the large, sign-alternating coefficients of a nearly
  !linearly dependent basis are accumulated. Method I obtains its energy in
  !this way from the stored shifted matrix. For a 2850-function basis in
  !extended precision the plain quotient was measured to be wrong by up to 15
  !units of the last place and the shifted quotient by at most one half,
  !which matters to an optimizer that resolves energy differences of that
  !size.
  !
  !The caller applies this routine to the vector returned by solve when the
  !eigenvalue is needed to the last digit. The routine is collective: the
  !columns of the lower triangles are assigned to the ranks cyclically, each
  !rank accumulates the two quadratic forms over its columns, and one
  !reduction of two values completes them. With one process no MPI call is
  !made.
  !
  !  Input parameters:
  !    h - Real array with at least n rows and n columns holding the lower
  !        triangle of the current H.
  !    s - Real array with at least n rows and n columns holding the lower
  !        triangle of the current S.
  !    x - Replicated vector of length n = order(), normally an approximate
  !        eigenvector. It need not be normalized.
  !
  !  Input parameter (state):
  !    self - A valid real QR state; only its order, shift, and distribution
  !           data are used.
  !
  !  Output parameters:
  !    lambda - The shifted Rayleigh quotient; zero on error.
  !    info   - QR_SUCCESS;
  !             QR_ERR_INVALID_STATE when self holds no valid factorization;
  !             QR_ERR_DIMENSION_MISMATCH when H, S, or x has an extent
  !             inconsistent with the order;
  !             QR_ERR_NONPOSITIVE_OVERLAP when x^T*S*x is not positive.
  !
  !No allocation is performed and the state is not modified.
    class(qr_real_state), intent(in) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: lambda
    integer, intent(out) :: info
    integer :: j, matrix_n, mpi_err
    real(wp) :: column_form, column_overlap, forms(2)

    lambda = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if

    !forms(1) accumulates x^T*(H-shift*S)*x and forms(2) x^T*S*x over the
    !columns owned by this rank. A column contributes its diagonal term and
    !twice its part below the diagonal.
    forms = 0.0_wp
    do j = self%rank + 1, matrix_n, self%nprocs
      call real_shifted_column(matrix_n - j, h(j:matrix_n, j), &
                               s(j:matrix_n, j), self%shift, &
                               x(j:matrix_n), column_form, column_overlap)
      forms(1) = forms(1) + column_form
      forms(2) = forms(2) + column_overlap
    end do
    if (self%nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, forms, 2, MPI_WP, MPI_SUM, self%comm, &
                         mpi_err)
    end if
    !The negated comparison also rejects NaN.
    if (.not. (forms(2) > tiny(1.0_wp))) then
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if
    lambda = self%shift + forms(1) / forms(2)
    info = QR_SUCCESS
  end subroutine real_shifted_quotient

  pure subroutine real_shifted_column(m, h_col, s_col, shift, x, form, overlap)
  !Subroutine real_shifted_column returns the contributions of one column of
  !the lower triangles to the quadratic forms x^T*(H-shift*S)*x and x^T*S*x:
  !the diagonal term plus twice the part below the diagonal.
  !
  !  Input parameters:
  !    m     - Number of elements below the diagonal.
  !    h_col - Column of H from its diagonal element downwards (m+1 elements).
  !    s_col - The same part of the column of S.
  !    shift - The shift.
  !    x     - The matching part of the vector (m+1 elements).
  !
  !  Output parameters:
  !    form    - Contribution to x^T*(H-shift*S)*x.
  !    overlap - Contribution to x^T*S*x.
    integer, intent(in) :: m
    real(wp), intent(in) :: h_col(0:m), s_col(0:m), shift, x(0:m)
    real(wp), intent(out) :: form, overlap
    integer :: i
    real(wp) :: form_below, overlap_below

    form_below = 0.0_wp
    overlap_below = 0.0_wp
    do i = 1, m
      form_below = form_below + (h_col(i) - shift * s_col(i)) * x(i)
      overlap_below = overlap_below + s_col(i) * x(i)
    end do
    form = x(0) * ((h_col(0) - shift * s_col(0)) * x(0) + 2.0_wp * form_below)
    overlap = x(0) * (s_col(0) * x(0) + 2.0_wp * overlap_below)
  end subroutine real_shifted_column

  subroutine real_factorization_residual(self, h, s, v, absolute_residual, &
                                         relative_residual, info, &
                                         forward_error)
  !Subroutine real_factorization_residual measures how well the represented
  !matrix (base factors plus border) reproduces the action of the caller-owned
  !physical shifted matrix A = H-shift*S on a probe vector. The comparison of
  !Q*(R*v) with A*v used by the previous implementation is not possible once
  !updates are deferred, because Q*R is only the base matrix. The check is
  !therefore made through the bordered solve:
  !
  !       b = A*v                      (physical action, two products)
  !       x = solve(represented A, b)  (bordered solve)
  !       difference = A*x - b         (physical action again)
  !
  !       absolute_residual = ||difference||_2,
  !       relative_residual = ||difference||_2 /
  !                           (||A*x||_2+||b||_2+tiny(1.0_wp)),
  !       forward_error     = ||x-v||_2/||v||_2      (optional).
  !
  !The relative residual has the same scaling as the action residual of the
  !previous implementation and is of the order of the working precision times
  !a modest factor when the representation is healthy, also for the nearly
  !singular shifted matrices of inverse iteration, because the large error of
  !x along the near-null vector is annihilated by A. The optional forward
  !error is the plain relative difference between x and v; it is proportional
  !to the condition number of A and is therefore large (harmlessly) when the
  !shift is close to an eigenvalue. Refresh policy belongs to the caller, who
  !should judge relative_residual. The routine is collective and costs one
  !bordered solve and four distributed symmetric products.
  !
  !  Input parameters:
  !    h - Real array with at least n rows and n columns holding the lower
  !        triangle of the current H.
  !    s - Real array with at least n rows and n columns holding the lower
  !        triangle of the current S.
  !    v - Nonzero real probe vector with exactly n elements. A fixed dense
  !        vector is appropriate; the eigenvector itself is not, since A
  !        nearly annihilates it.
  !
  !  Input/output parameter:
  !    self - Valid real QR state. The represented matrix and the counters are
  !           unchanged; the workspace is overwritten.
  !
  !  Output parameters:
  !    absolute_residual - Euclidean norm of difference; zero on error.
  !    relative_residual - Scaled action error defined above; zero on error.
  !    info              - QR_SUCCESS;
  !                        QR_ERR_INVALID_STATE when no valid factorization is
  !                        available;
  !                        QR_ERR_DIMENSION_MISMATCH when H or S cannot contain
  !                        the active block or V does not have length n;
  !                        QR_ERR_ZERO_INITIAL_VECTOR when V is numerically
  !                        zero;
  !                        QR_ERR_SINGULAR when R or a small matrix of the
  !                        border cannot be used.
  !    forward_error     - Optional relative forward error defined above; zero
  !                        on error.
  !
  !The routine performs no allocation and no input/output.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: v(:)
    real(wp), intent(out) :: absolute_residual, relative_residual
    integer, intent(out) :: info
    real(wp), intent(out), optional :: forward_error
    integer :: i, matrix_n
    real(wp) :: action_norm, difference, difference_norm_squared
    real(wp) :: probe_norm_squared, rhs_norm

    absolute_residual = 0.0_wp
    relative_residual = 0.0_wp
    if (present(forward_error)) forward_error = 0.0_wp
    info = QR_ERR_INVALID_STATE

    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    do i = 1, matrix_n
      self%wn(i, 5) = v(i)
    end do
    probe_norm_squared = real_norm_squared(matrix_n, self%wn(:, 5))
    if (probe_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if
    call real_refresh_small(self)
    if (.not. real_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !b = (H-shift*S)*v in work vector 3; x in work vector 4.
    call real_sym_product(self%rank, self%nprocs, self%comm, matrix_n, h, &
                          size(h, 1), self%wn(:, 5), self%wn(:, 3))
    call real_sym_product(self%rank, self%nprocs, self%comm, matrix_n, s, &
                          size(s, 1), self%wn(:, 5), self%wn(:, 1))
    do i = 1, matrix_n
      self%wn(i, 3) = self%wn(i, 3) - self%shift * self%wn(i, 1)
    end do
    call real_bordered_solve(self, self%wn(:, 3), self%wn(:, 4))

    if (present(forward_error)) then
      difference_norm_squared = 0.0_wp
      do i = 1, matrix_n
        difference = self%wn(i, 4) - self%wn(i, 5)
        difference_norm_squared = difference_norm_squared + &
                                  difference * difference
      end do
      forward_error = sqrt(difference_norm_squared / probe_norm_squared)
    end if

    !A*x in work vector 1, using work vector 2 for S*x.
    call real_sym_product(self%rank, self%nprocs, self%comm, matrix_n, h, &
                          size(h, 1), self%wn(:, 4), self%wn(:, 1))
    call real_sym_product(self%rank, self%nprocs, self%comm, matrix_n, s, &
                          size(s, 1), self%wn(:, 4), self%wn(:, 2))
    do i = 1, matrix_n
      self%wn(i, 1) = self%wn(i, 1) - self%shift * self%wn(i, 2)
    end do
    rhs_norm = sqrt(real_norm_squared(matrix_n, self%wn(:, 3)))
    action_norm = sqrt(real_norm_squared(matrix_n, self%wn(:, 1)))
    difference_norm_squared = 0.0_wp
    do i = 1, matrix_n
      difference = self%wn(i, 1) - self%wn(i, 3)
      difference_norm_squared = difference_norm_squared + &
                                difference * difference
    end do
    absolute_residual = sqrt(difference_norm_squared)
    relative_residual = absolute_residual / &
                        (action_norm + rhs_norm + tiny(1.0_wp))
    info = QR_SUCCESS
  end subroutine real_factorization_residual

  pure function real_is_valid(self) result(factors_are_valid)
  !Function real_is_valid reports whether a real QR state currently contains
  !a completed base factorization. It returns true after a successful
  !factorize_fresh and false for a default, cleared, or merely initialized
  !state, and after a failed factorization. It does not test numerical
  !singularity or the health of the border.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    factors_are_valid - Stored validity flag.
    class(qr_real_state), intent(in) :: self
    logical :: factors_are_valid

    factors_are_valid = self%valid
  end function real_is_valid

  pure function real_order(self) result(active_order)
  !Function real_order returns the current order n of the represented matrix.
  !It grows with append_symmetric, shrinks with truncate and
  !delete_symmetric, and is zero when no factorization is held. It may be
  !smaller or larger than base_order().
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    active_order - Current order n.
    class(qr_real_state), intent(in) :: self
    integer :: active_order

    active_order = self%n
  end function real_order

  pure function real_base_order(self) result(order_of_base)
  !Function real_base_order returns the order n0 of the base factorization,
  !that is, order() at the time of the last successful factorize_fresh.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    order_of_base - Base order n0; zero when no factorization is held.
    class(qr_real_state), intent(in) :: self
    integer :: order_of_base

    order_of_base = self%n0
  end function real_base_order

  pure function real_get_capacity(self) result(allocated_capacity)
  !Function real_get_capacity returns the maximum matrix order for which the
  !state has allocated storage. A cleared or failed-initialization state
  !returns zero.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    allocated_capacity - Maximum supported order.
    class(qr_real_state), intent(in) :: self
    integer :: allocated_capacity

    allocated_capacity = self%capacity
  end function real_get_capacity

  pure function real_get_shift(self) result(represented_shift)
  !Function real_get_shift returns the shift of the represented matrix
  !H-shift*S. The value is updated only by a successful factorize_fresh and
  !is zero in an initialized state that has not been factorized.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    represented_shift - Stored shift.
    class(qr_real_state), intent(in) :: self
    real(wp) :: represented_shift

    represented_shift = self%shift
  end function real_get_shift

  pure function real_pending_count(self) result(count)
  !Function real_pending_count returns the number m of border entries: the
  !functions whose current row and column are held outside the base
  !factorization, plus the trailing base slots deleted by truncate (each of
  !them occupies an entry until the next fresh factorization). Announced
  !appended indices that have not been supplied yet are not counted. The
  !border is full when the value reaches max_pending().
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    count - Number of border entries, 0 <= count <= max_pending().
    class(qr_real_state), intent(in) :: self
    integer :: count

    count = self%m
  end function real_pending_count

  pure function real_active_count(self) result(count)
  !Function real_active_count returns the size of the current active set: the
  !entries of the active tail (including deleted base slots that are waiting
  !to be reused by an append) plus the announced appended indices beyond the
  !base that have not been supplied yet.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    count - Size of the active set.
    class(qr_real_state), intent(in) :: self
    integer :: count

    count = self%ma + self%nres
  end function real_active_count

  pure function real_excluded_count(self) result(count)
  !Function real_excluded_count returns the number m0 of members of D0, the
  !base slots whose row and column are modified: pending functions with an
  !index not greater than base_order() and trailing base slots deleted by
  !truncate. Every member owns a border entry, so the value never exceeds
  !pending_count(); only a fresh factorization reduces it.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    count - Number of modified base slots.
    class(qr_real_state), intent(in) :: self
    integer :: count

    count = self%m0
  end function real_excluded_count

  pure function real_max_pending(self) result(count)
  !Function real_max_pending returns the largest number of border entries the
  !state can hold, min(max_pending,capacity) of the initialize call.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    count - Border capacity; zero for an uninitialized state.
    class(qr_real_state), intent(in) :: self
    integer :: count

    count = self%mmax
  end function real_max_pending

  pure function real_get_update_count(self) result(update_count)
  !Function real_get_update_count returns the lifetime number of successful
  !structural updates (replace_symmetric and append_symmetric calls, and
  !functions removed by truncate or delete_symmetric). Fresh factorization
  !does not reset the value; clear and initialize do.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    update_count - Lifetime count as a 64-bit integer.
    class(qr_real_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%structural_updates
  end function real_get_update_count

  pure function real_get_updates_since_fresh(self) result(update_count)
  !Function real_get_updates_since_fresh returns the number of successful
  !structural updates since the most recent successful factorize_fresh.
  !
  !  Input parameter:
  !    self - The real QR state to query.
  !
  !  Result:
  !    update_count - Count since the last fresh factorization.
    class(qr_real_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%updates_since_fresh
  end function real_get_updates_since_fresh

  !===========================================================================
  ! Internal procedures: border bookkeeping and small matrices
  !===========================================================================

  subroutine real_reset_border(self)
  !Subroutine real_reset_border empties the border list and D0 and restores
  !the free-slot list. It is called when the base factorization is replaced.
  !
  !  Input/output parameter:
  !    self - State with allocated bookkeeping arrays.
    class(qr_real_state), intent(inout) :: self
    integer :: i

    self%m = 0
    self%ma = 0
    self%nres = 0
    self%m0 = 0
    self%bidx = 0
    self%bpos = 0
    self%d0idx = 0
    self%d0pos = 0
    self%stale = 0
    self%cu = 0
    self%imark = 0
    do i = 1, self%mmax
      self%bslot(i) = i
    end do
    self%scc_ok = .true.
    self%tail_ok = .true.
    self%scc_singular = .false.
    self%sig_singular = .false.
  end subroutine real_reset_border

  pure function real_can_solve(self) result(usable)
  !Function real_can_solve reports whether the bordered solve may divide by
  !the stored pivots: R, Scc and Sigma are all usable and up to date.
  !
  !  Input parameter:
  !    self - State after real_refresh_small.
  !
  !  Result:
  !    usable - True when a bordered solve can be performed.
    class(qr_real_state), intent(in) :: self
    logical :: usable

    usable = .not. (self%r_singular .or. self%scc_singular .or. &
                    self%sig_singular) .and. self%scc_ok .and. self%tail_ok
  end function real_can_solve

  subroutine real_join_base_slot(self, idx, deleted)
  !Subroutine real_join_base_slot makes the clean base slot idx a member of D0
  !and creates its border entry at the end of the border list. It computes
  !the column Zb(:,q) = A_base^{-1}*e_idx by one distributed base solve and
  !recovers the base column Ab(:,q) = Q*(R*e_idx) from the factors. The new
  !entry starts with a zero modification (the function is unchanged), so its
  !border column and its base solve are zero. The entry of every older border
  !column on row idx moves into the border-border block, and the stored base
  !solve of that column is corrected by a multiple of Zb(:,q).
  !
  !The caller is responsible for the position of the entry in the border list
  !(tail or committed), for the counter of the tail, and for calling
  !real_rebuild_small afterwards.
  !
  !  Input parameters:
  !    idx     - Base slot, 1 <= idx <= n0, not yet a member of D0. The caller
  !              has verified that the border has room and that R is usable.
  !    deleted - True when the slot is being deleted by truncate. The entry
  !              is then marked stale; its decoupling column is installed
  !              later by real_install_column.
  !
  !  Input/output parameter:
  !    self - State; m and m0 grow by one. Work vector 1 is overwritten.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    logical, intent(in) :: deleted
    integer :: i, idx2, p, p2, q, slot, slot2
    real(wp) :: current, entry

    q = self%m0 + 1
    self%m0 = q
    self%d0idx(q) = idx
    self%d0pos(idx) = q
    do i = 1, self%n0
      self%wn(i, 1) = 0.0_wp
    end do
    self%wn(idx, 1) = 1.0_wp
    call real_base_solve(self, self%wn(:, 1), self%zb(:, q))
    call real_base_column(self, idx, self%ab(:, q), self%wn(:, 1))

    p = self%m + 1
    slot = self%bslot(p)
    self%bidx(p) = idx
    self%bpos(idx) = p
    self%m = p
    do i = 1, self%n0
      self%bt(i, slot) = 0.0_wp
      self%gs(i, slot) = 0.0_wp
    end do
    if (deleted) then
      self%stale(slot) = 1
    else
      self%stale(slot) = 0
    end if
    do p2 = 1, p
      slot2 = self%bslot(p2)
      self%cbb(slot, slot2) = 0.0_wp
      self%cbb(slot2, slot) = 0.0_wp
    end do
    self%cbb(slot, slot) = self%ab(idx, q)

    do p2 = 1, p - 1
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      entry = self%bt(idx, slot2)
      !Current value of the matrix element (idx,idx2): the stored column of a
      !member of D0 is a difference from its base column.
      current = entry
      if (idx2 <= self%n0) current = current + self%ab(idx, self%d0pos(idx2))
      if (deleted .or. self%stale(slot2) /= 0) current = 0.0_wp
      self%cbb(slot, slot2) = current
      self%cbb(slot2, slot) = current
      if (entry /= 0.0_wp) then
        self%bt(idx, slot2) = 0.0_wp
        call real_axpy(self%n0, -entry, self%zb(:, q), self%gs(:, slot2))
      end if
    end do
  end subroutine real_join_base_slot

  subroutine real_install_column(self, p, a, dummy)
  !Subroutine real_install_column installs the current shifted column of the
  !border entry at position p. Its part on the clean rows goes to Bt (as the
  !difference from the base column when the entry is a member of D0), its
  !values on the border indices go to row and column of the border-border
  !block, Gs is obtained by one distributed base solve, and the row and column
  !of BG, the row of BZ and the row and column of the Schur complement that
  !belong to the entry are updated.
  !
  !With dummy true the entry is a deleted base slot and receives its
  !decoupling column: zero off the diagonal and the base diagonal element on
  !the diagonal (one if that element is zero).
  !
  !  Input parameters:
  !    p     - Border position of the entry.
  !    a     - Shifted column of length n (current order); not referenced
  !            when dummy is true.
  !    dummy - Selects the decoupling column of a deleted base slot.
  !
  !  Input/output parameter:
  !    self - State. The entry is no longer stale; tail_ok (entry in the tail)
  !           or scc_ok (committed entry) becomes false.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: p
    real(wp), intent(in) :: a(self%n)
    logical, intent(in) :: dummy
    integer :: clean_rows, i, idx, idx2, p2, q, slot, slot2
    real(wp) :: current, value

    slot = self%bslot(p)
    idx = self%bidx(p)
    q = 0
    if (idx <= self%n0) q = self%d0pos(idx)
    clean_rows = min(self%n, self%n0)
    do i = 1, clean_rows
      if (self%d0pos(i) /= 0) then
        self%bt(i, slot) = 0.0_wp
      else
        value = 0.0_wp
        if (.not. dummy) value = a(i)
        if (q > 0) value = value - self%ab(i, q)
        self%bt(i, slot) = value
      end if
    end do
    !Base rows beyond the current order are deleted slots, members of D0.
    do i = clean_rows + 1, self%n0
      self%bt(i, slot) = 0.0_wp
    end do
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      if (p2 == p) then
        if (dummy) then
          current = self%ab(idx, q)
          if (current == 0.0_wp) current = 1.0_wp
        else
          current = a(idx)
        end if
      else if (dummy .or. idx2 > self%n) then
        !A deleted slot is decoupled from everything else.
        current = 0.0_wp
      else
        current = a(idx2)
      end if
      self%cbb(slot, slot2) = current
      self%cbb(slot2, slot) = current
    end do

    call real_base_solve(self, self%bt(:, slot), self%gs(:, slot))

    !Row and column of BG and row of BZ.
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      self%bg(slot, slot2) = real_dot(self%n0, self%bt(:, slot), &
                                      self%gs(:, slot2))
      if (slot2 /= slot) then
        self%bg(slot2, slot) = real_dot(self%n0, self%bt(:, slot2), &
                                        self%gs(:, slot))
      end if
    end do
    do q = 1, self%m0
      self%bz(slot, q) = real_dot(self%n0, self%bt(:, slot), self%zb(:, q))
    end do
    call real_fill_schur(self, p)

    self%stale(slot) = 0
    if (p > self%m - self%ma) then
      self%tail_ok = .false.
    else
      self%scc_ok = .false.
    end if
  end subroutine real_install_column

  subroutine real_fill_schur(self, p)
  !Subroutine real_fill_schur recomputes the row and the column of the Schur
  !complement Sc that belong to the border entry at position p from the
  !stored small matrices. With u = mmax+slot the unknown of the entry, own(q)
  !equal to one when member q of D0 is the entry itself and zero otherwise,
  !and base(i,j) the base matrix element of two members of D0 (zero when
  !either entry is an appended function),
  !
  !        Sc(q,u)   = -own(q) - Gs(D0(q),slot)
  !        Sc(u,q)   = -own(q) - BZ(slot,q)
  !        Sc(u,u2)  = C(slot,slot2) - base(idx,idx2)  - BG(slot,slot2)
  !        Sc(u2,u)  = C(slot2,slot) - base(idx2,idx)  - BG(slot2,slot) .
  !
  !  Input parameter:
  !    p    - Border position of the entry.
  !
  !  Input/output parameter:
  !    self - State with valid cbb, bg, bz, gs and ab for the entry.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: p
    integer :: idx, idx2, p2, q, q2, qs, slot, slot2, u, u2
    real(wp) :: base12, base21, own

    slot = self%bslot(p)
    idx = self%bidx(p)
    u = self%mmax + slot
    qs = 0
    if (idx <= self%n0) qs = self%d0pos(idx)
    do q = 1, self%m0
      own = 0.0_wp
      if (q == qs) own = 1.0_wp
      self%sc(q, u) = -own - self%gs(self%d0idx(q), slot)
      self%sc(u, q) = -own - self%bz(slot, q)
    end do
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      u2 = self%mmax + slot2
      base12 = 0.0_wp
      base21 = 0.0_wp
      if (qs > 0 .and. idx2 <= self%n0) then
        q2 = self%d0pos(idx2)
        base12 = self%ab(idx, q2)
        base21 = self%ab(idx2, qs)
      end if
      self%sc(u, u2) = self%cbb(slot, slot2) - base12 - self%bg(slot, slot2)
      self%sc(u2, u) = self%cbb(slot2, slot) - base21 - self%bg(slot2, slot)
    end do
  end subroutine real_fill_schur

  subroutine real_rebuild_small(self)
  !Subroutine real_rebuild_small recomputes all small matrices from the stored
  !border columns: BG = Bt^T*Gs and BZ = Bt^T*Zb by a distributed pass in which
  !every rank takes a contiguous slice of rows (one reduction), and the whole
  !Schur complement Sc. It is called whenever D0 grows, which changes the
  !stored columns of all border entries; it is never called per trial.
  !
  !  Input/output parameter:
  !    self - State. The factors of Scc and of the tail are marked out of
  !           date.
    class(qr_real_state), intent(inout) :: self
    integer :: chunk, count, hi, k, lo, m, m0, mpi_err, p, p2, q, q2
    integer :: slot, slot2

    m = self%m
    m0 = self%m0
    self%scc_ok = .false.
    self%tail_ok = .false.

    if (self%nprocs > 1) then
      chunk = (self%n0 + self%nprocs - 1) / self%nprocs
      lo = self%rank * chunk + 1
      hi = min(self%n0, (self%rank + 1) * chunk)
    else
      lo = 1
      hi = self%n0
    end if
    count = m * m + m * m0
    if (count > 0) then
      k = 0
      do p2 = 1, m
        slot2 = self%bslot(p2)
        do p = 1, m
          slot = self%bslot(p)
          k = k + 1
          if (hi >= lo) then
            self%wred(k) = real_dot(hi - lo + 1, self%bt(lo:hi, slot), &
                                    self%gs(lo:hi, slot2))
          else
            self%wred(k) = 0.0_wp
          end if
        end do
      end do
      do q = 1, m0
        do p = 1, m
          slot = self%bslot(p)
          k = k + 1
          if (hi >= lo) then
            self%wred(k) = real_dot(hi - lo + 1, self%bt(lo:hi, slot), &
                                    self%zb(lo:hi, q))
          else
            self%wred(k) = 0.0_wp
          end if
        end do
      end do
      if (self%nprocs > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, self%wred, count, MPI_WP, MPI_SUM, &
                           self%comm, mpi_err)
      end if
      k = 0
      do p2 = 1, m
        slot2 = self%bslot(p2)
        do p = 1, m
          k = k + 1
          self%bg(self%bslot(p), slot2) = self%wred(k)
        end do
      end do
      do q = 1, m0
        do p = 1, m
          k = k + 1
          self%bz(self%bslot(p), q) = self%wred(k)
        end do
      end do
    end if

    !Unit-vector block of the Schur complement, -E^T*Zb, and the rows and
    !columns of all border entries.
    do q2 = 1, m0
      do q = 1, m0
        self%sc(q, q2) = -self%zb(self%d0idx(q), q2)
      end do
    end do
    do p = 1, m
      call real_fill_schur(self, p)
    end do
  end subroutine real_rebuild_small

  subroutine real_refresh_committed(self)
  !Subroutine real_refresh_committed brings the committed part of the small
  !system up to date: it installs the decoupling columns of committed deleted
  !slots that are still stale, and recomputes the LU factorization of the
  !committed block Scc when it is out of date. The committed unknowns are the
  !unit-vector unknowns of all members of D0 followed by the committed border
  !entries in border order.
  !
  !  Input/output parameter:
  !    self - State. scc_singular records the result of the relative pivot
  !           test; scc_ok becomes true when the factors are usable.
    class(qr_real_state), intent(inout) :: self
    integer :: k, k2, mc, ncu, p

    if (self%r_singular) return
    mc = self%m - self%ma
    do p = 1, mc
      if (self%stale(self%bslot(p)) /= 0) then
        call real_install_column(self, p, self%wn(:, 3), .true.)
      end if
    end do
    if (self%scc_ok) return
    ncu = self%m0 + mc
    do k = 1, self%m0
      self%cu(k) = k
    end do
    do p = 1, mc
      self%cu(self%m0 + p) = self%mmax + self%bslot(p)
    end do
    do k2 = 1, ncu
      do k = 1, ncu
        self%scc_lu(k, k2) = self%sc(self%cu(k), self%cu(k2))
      end do
    end do
    call real_lu_factor(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                        QR_RELATIVE_PIVOT, self%scc_mag, self%scc_singular)
    self%scc_ok = .not. self%scc_singular
    self%tail_ok = .false.
  end subroutine real_refresh_committed

  subroutine real_refresh_small(self)
  !Subroutine real_refresh_small brings the lazily maintained factors up to
  !date before a bordered solve: the committed block (see
  !real_refresh_committed), the decoupling columns of stale deleted slots in
  !the tail, and for the active tail X = Scc^{-1}*Sca, Sigma = Saa - Sac*X
  !and LU(Sigma) if a tail column changed. All work is O(m**2*Ma) or less,
  !except a refactorization of Scc and one base solve per stale slot.
  !
  !Sigma is the Schur complement of everything else in the nearly singular
  !shifted matrix and is therefore legitimately tiny; its pivots are rejected
  !only when they are zero, subnormal, or not finite.
  !
  !  Input/output parameter:
  !    self - State; see real_can_solve for the resulting flags.
    class(qr_real_state), intent(inout) :: self
    integer :: a, a2, k, mc, ncu, p, u
    real(wp) :: tsum

    if (self%r_singular) return
    mc = self%m - self%ma
    do p = mc + 1, self%m
      if (self%stale(self%bslot(p)) /= 0) then
        call real_install_column(self, p, self%wn(:, 3), .true.)
      end if
    end do
    call real_refresh_committed(self)
    if (.not. self%scc_ok) return
    if (self%tail_ok) return
    ncu = self%m0 + mc
    do a = 1, self%ma
      u = self%mmax + self%bslot(mc + a)
      do k = 1, ncu
        self%xca(k, a) = self%sc(self%cu(k), u)
      end do
      if (ncu > 0) then
        call real_lu_solve(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                           self%xca(:, a))
      end if
    end do
    do a2 = 1, self%ma
      do a = 1, self%ma
        u = self%mmax + self%bslot(mc + a)
        tsum = 0.0_wp
        do k = 1, ncu
          tsum = tsum + self%sc(u, self%cu(k)) * self%xca(k, a2)
        end do
        self%sig_lu(a, a2) = self%sc(u, self%mmax + self%bslot(mc + a2)) - tsum
      end do
    end do
    call real_lu_factor(self%ma, self%sig_lu, self%mmax, self%sig_piv, &
                        0.0_wp, self%scc_mag, self%sig_singular)
    self%tail_ok = .not. self%sig_singular
  end subroutine real_refresh_small

  !===========================================================================
  ! Internal procedures: distributed kernels
  !===========================================================================

  subroutine real_bordered_solve(self, r, x)
  !Subroutine real_bordered_solve solves A*x = r for the represented matrix of
  !the current order n by block elimination of the augmented system described
  !in the module header:
  !
  !        r0 = r on the base indices, zero on deleted slots
  !        g  = A_base^{-1}*r0
  !        t  = [0; 0; rA] - [E D B]^T*g
  !        y  = Sc^{-1}*t      (block elimination on the committed block)
  !        x0 = g - Zb*yE - Gs*[yD; xA] .
  !
  !One distributed base solve is performed; the remaining work,
  !O((2*m+m0)*n0+(m+m0)**2), is executed identically on every rank.
  !
  !  Input parameter:
  !    r - Replicated right-hand side of length n. It must not be associated
  !        with work vector 1 or 2 of the state.
  !
  !  Input/output parameter:
  !    self - State for which real_can_solve is true. Work vectors 1 and 2 and
  !           the small work vectors are overwritten.
  !
  !  Output parameter:
  !    x - Replicated solution of length n; not associated with work vector 1
  !        or 2 of the state.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: r(self%n)
    real(wp), intent(out) :: x(self%n)
    integer :: a, clean_rows, i, k, m, m0, mc, ncu, p, q, slot, u
    real(wp) :: tsum

    m = self%m
    m0 = self%m0
    mc = m - self%ma
    ncu = m0 + mc
    clean_rows = min(self%n, self%n0)

    !Right-hand side on the base indices and base solve. Work vector 1 is
    !r0, work vector 2 is g and finally x0.
    do i = 1, clean_rows
      self%wn(i, 1) = r(i)
    end do
    do i = clean_rows + 1, self%n0
      self%wn(i, 1) = 0.0_wp
    end do
    call real_base_solve(self, self%wn(:, 1), self%wn(:, 2))

    if (m > 0) then
      !Right-hand side of the small system: committed unknowns in small
      !vector 1 (order of cu), tail unknowns in small vector 2.
      do q = 1, m0
        self%wsm(q, 1) = -self%wn(self%d0idx(q), 2)
      end do
      do p = 1, m
        slot = self%bslot(p)
        tsum = -real_dot(self%n0, self%bt(:, slot), self%wn(:, 2))
        if (self%bidx(p) > self%n0) tsum = tsum + r(self%bidx(p))
        if (p <= mc) then
          self%wsm(m0 + p, 1) = tsum
        else
          self%wsm(p - mc, 2) = tsum
        end if
      end do

      !yc = Scc^{-1}*tc, ya = Sigma^{-1}*(ta - Sac*yc), yc = yc - X*ya.
      if (ncu > 0) then
        call real_lu_solve(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                           self%wsm(:, 1))
      end if
      do a = 1, self%ma
        u = self%mmax + self%bslot(mc + a)
        tsum = self%wsm(a, 2)
        do k = 1, ncu
          tsum = tsum - self%sc(u, self%cu(k)) * self%wsm(k, 1)
        end do
        self%wsm(a, 2) = tsum
      end do
      if (self%ma > 0) then
        call real_lu_solve(self%ma, self%sig_lu, self%mmax, self%sig_piv, &
                           self%wsm(:, 2))
        do a = 1, self%ma
          tsum = self%wsm(a, 2)
          do k = 1, ncu
            self%wsm(k, 1) = self%wsm(k, 1) - self%xca(k, a) * tsum
          end do
        end do
      end if

      !x0 = g - Zb*yE - Gs*y, with the negated coefficients in small vector 3.
      if (m0 > 0) then
        do q = 1, m0
          self%wsm(q, 3) = -self%wsm(q, 1)
        end do
        call real_add_leading_columns(self%n0, self%zb, self%capacity, m0, &
                                      self%wsm(:, 3), self%wn(:, 2))
      end if
      if (mc > 0) then
        do p = 1, mc
          self%wsm(p, 3) = -self%wsm(m0 + p, 1)
        end do
        call real_add_columns(self%n0, self%gs, self%capacity, mc, &
                              self%bslot(1:mc), self%wsm(:, 3), self%wn(:, 2))
      end if
      if (self%ma > 0) then
        do a = 1, self%ma
          self%wsm(a, 3) = -self%wsm(a, 2)
        end do
        call real_add_columns(self%n0, self%gs, self%capacity, self%ma, &
                              self%bslot(mc + 1:m), self%wsm(:, 3), &
                              self%wn(:, 2))
      end if
    end if

    do i = 1, clean_rows
      x(i) = self%wn(i, 2)
    end do
    do p = 1, m
      if (self%bidx(p) > self%n0) then
        if (p <= mc) then
          x(self%bidx(p)) = self%wsm(m0 + p, 1)
        else
          x(self%bidx(p)) = self%wsm(p - mc, 2)
        end if
      end if
    end do
  end subroutine real_bordered_solve

  subroutine real_base_solve(self, b, x)
  !Subroutine real_base_solve computes x = A_base^{-1}*b = R^{-1}*(Q^T*b) with
  !the distributed base factors.
  !
  !  1. y = Q^T*b: every rank accumulates b(i) times its rows of Q (contiguous
  !     columns of qt), followed by one reduction of n0 values.
  !  2. Back substitution by blocks of P rows from the bottom. Every rank owns
  !     exactly one row of a block and computes the part of its row of R to
  !     the right of the block times the already known solution (one dot
  !     product on a contiguous column of rt). One reduction of at most P
  !     values makes these products global; then every rank solves the P by P
  !     triangular block from the replicated diagonal blocks. The bottom block
  !     needs no reduction.
  !
  !With one process the routine performs an ordinary back substitution and
  !makes no MPI call.
  !
  !  Input parameter:
  !    b - Replicated right-hand side of length n0.
  !
  !  Input/output parameter:
  !    self - State with usable R (r_singular is false). The local work
  !           arrays wloc and wblk are overwritten.
  !
  !  Output parameter:
  !    x - Replicated solution of length n0. It must not overlap b.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: b(self%n0)
    real(wp), intent(out) :: x(self%n0)
    integer :: blk, i, ie, ii, is, jj, l, mpi_err, n0, nb, np
    real(wp) :: tsum

    n0 = self%n0
    np = self%nprocs

    !y = Q^T*b, accumulated in x.
    do l = 1, self%nloc
      self%wloc(l) = b(self%rank + 1 + (l - 1) * np)
    end do
    do i = 1, n0
      x(i) = 0.0_wp
    end do
    call real_add_leading_columns(n0, self%qt, self%capacity, self%nloc, &
                                  self%wloc, x)
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, x, n0, MPI_WP, MPI_SUM, self%comm, &
                         mpi_err)
    end if

    if (np == 1) then
      do i = n0, 1, -1
        tsum = x(i)
        if (i < n0) then
          tsum = tsum - real_dot(n0 - i, self%rt(i + 1:n0, i), x(i + 1:n0))
        end if
        x(i) = tsum / self%rt(i, i)
      end do
    else
      do blk = self%nblk, 1, -1
        is = (blk - 1) * np + 1
        ie = min(blk * np, n0)
        nb = ie - is + 1
        if (ie < n0) then
          do ii = 1, nb
            self%wblk(ii) = 0.0_wp
          end do
          if (self%rank < nb) then
            self%wblk(self%rank + 1) = real_dot(n0 - ie, &
                                                self%rt(ie + 1:n0, blk), &
                                                x(ie + 1:n0))
          end if
          call MPI_Allreduce(MPI_IN_PLACE, self%wblk, nb, MPI_WP, MPI_SUM, &
                             self%comm, mpi_err)
          do ii = 1, nb
            x(is + ii - 1) = x(is + ii - 1) - self%wblk(ii)
          end do
        end if
        do ii = nb, 1, -1
          tsum = x(is + ii - 1)
          do jj = ii + 1, nb
            tsum = tsum - self%rd(ii, jj, blk) * x(is + jj - 1)
          end do
          x(is + ii - 1) = tsum / self%rd(ii, ii, blk)
        end do
      end do
    end if
  end subroutine real_base_solve

  subroutine real_base_column(self, idx, column, work)
  !Subroutine real_base_column recovers column idx of the base matrix from the
  !distributed factors, column = Q*(R*e_idx). Every rank contributes its
  !entries of column idx of R (one reduction) and then the products of its
  !rows of Q with that column (one reduction).
  !
  !  Input parameter:
  !    idx - Column index, 1 <= idx <= n0.
  !
  !  Input/output parameter:
  !    self - State with valid base factors (not modified; declared
  !           intent(inout) because the output arguments are normally work
  !           vectors owned by the state).
  !
  !  Output parameters:
  !    column - Replicated column of length n0.
  !    work   - Work vector of length n0 (column idx of R on exit).
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(out) :: column(self%n0), work(self%n0)
    integer :: i, l, mpi_err, n0, np

    n0 = self%n0
    np = self%nprocs
    do i = 1, n0
      work(i) = 0.0_wp
      column(i) = 0.0_wp
    end do
    do l = 1, self%nloc
      i = self%rank + 1 + (l - 1) * np
      if (i <= idx) work(i) = self%rt(idx, l)
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, work, n0, MPI_WP, MPI_SUM, self%comm, &
                         mpi_err)
    end if
    do l = 1, self%nloc
      i = self%rank + 1 + (l - 1) * np
      column(i) = real_dot(idx, self%qt(1:idx, l), work(1:idx))
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, column, n0, MPI_WP, MPI_SUM, &
                         self%comm, mpi_err)
    end if
  end subroutine real_base_column

  subroutine real_sym_product(rank, nprocs, comm, n, a, lda, x, y)
  !Subroutine real_sym_product computes y = A*x for a real symmetric matrix of
  !which only the lower triangle is stored (and replicated on all ranks). The
  !columns are assigned to the ranks cyclically; a rank accumulates both the
  !column contribution (below the diagonal) and the mirrored row contribution
  !of its columns in one pass over each column, and one reduction of n values
  !completes the product. With one process no MPI call is made.
  !
  !The owned columns are processed several at a time over their common rows,
  !so that y is read and written once per group; the few rows above the
  !common range are handled column by column. Groups of four are used except
  !for 80-bit extended precision (QR_PRODUCT_PAIRS), where pairs keep all
  !running values in the eight registers of the x87 unit and were measured
  !to be faster.
  !
  !  Input parameters:
  !    rank, nprocs, comm - Own rank, number of ranks, and communicator.
  !    n   - Order of the matrix.
  !    a   - Array with leading dimension lda whose leading n by n lower
  !          triangle holds the matrix. The upper triangle is not referenced.
  !    lda - Leading dimension of a.
  !    x   - Replicated vector of length n.
  !
  !  Output parameter:
  !    y   - Replicated product of length n. It must not overlap x.
    integer, intent(in) :: rank, nprocs, comm, n, lda
    real(wp), intent(in) :: a(lda, *), x(n)
    real(wp), intent(out) :: y(n)
    integer :: i, j, j2, j3, j4, mpi_err
    real(wp) :: t1, t2, t3, t4, u1, u2, u3

    do i = 1, n
      y(i) = 0.0_wp
    end do
    j = rank + 1
    if (QR_PRODUCT_PAIRS) then
      !Pairs of owned columns j < j2 with common rows j2+1..n.
      do while (j + nprocs < n)
        j2 = j + nprocs
        call real_symv_column(j2 - j, a(j + 1, j), x(j), x(j + 1), y(j + 1), &
                              u1)
        call real_symv_columns2(n - j2, a(j2 + 1, j), a(j2 + 1, j2), x(j), &
                                x(j2), x(j2 + 1), y(j2 + 1), t1, t2)
        y(j) = y(j) + a(j, j) * x(j) + (u1 + t1)
        y(j2) = y(j2) + a(j2, j2) * x(j2) + t2
        j = j2 + nprocs
      end do
    else
      !Groups of four owned columns j < j2 < j3 < j4 with common rows
      !j4+1..n.
      do while (j + 3 * nprocs < n)
        j2 = j + nprocs
        j3 = j2 + nprocs
        j4 = j3 + nprocs
        call real_symv_column(j4 - j, a(j + 1, j), x(j), x(j + 1), y(j + 1), &
                              u1)
        call real_symv_column(j4 - j2, a(j2 + 1, j2), x(j2), x(j2 + 1), &
                              y(j2 + 1), u2)
        call real_symv_column(j4 - j3, a(j3 + 1, j3), x(j3), x(j3 + 1), &
                              y(j3 + 1), u3)
        call real_symv_columns4(n - j4, a(j4 + 1, j), a(j4 + 1, j2), &
                                a(j4 + 1, j3), a(j4 + 1, j4), x(j), x(j2), &
                                x(j3), x(j4), x(j4 + 1), y(j4 + 1), t1, t2, &
                                t3, t4)
        y(j) = y(j) + a(j, j) * x(j) + (u1 + t1)
        y(j2) = y(j2) + a(j2, j2) * x(j2) + (u2 + t2)
        y(j3) = y(j3) + a(j3, j3) * x(j3) + (u3 + t3)
        y(j4) = y(j4) + a(j4, j4) * x(j4) + t4
        j = j4 + nprocs
      end do
    end if
    !Remaining owned columns one at a time.
    do while (j <= n)
      if (j < n) then
        call real_symv_column(n - j, a(j + 1, j), x(j), x(j + 1), y(j + 1), t1)
        y(j) = y(j) + a(j, j) * x(j) + t1
      else
        y(j) = y(j) + a(j, j) * x(j)
      end if
      j = j + nprocs
    end do
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, y, n, MPI_WP, MPI_SUM, comm, mpi_err)
    end if
  end subroutine real_sym_product

  !===========================================================================
  ! Internal procedures: local kernels
  !===========================================================================

  pure function qr_local_count(order, rank, nprocs) result(count)
  !Function qr_local_count returns the number of rows i in 1..order that the
  !cyclic distribution assigns to a rank, that is, with mod(i-1,nprocs)=rank.
  !
  !  Input parameters:
  !    order  - Number of rows.
  !    rank   - Rank, 0 <= rank < nprocs.
  !    nprocs - Number of ranks.
  !
  !  Result:
  !    count - Number of local rows.
    integer, intent(in) :: order, rank, nprocs
    integer :: count

    if (order > rank) then
      count = (order - rank - 1) / nprocs + 1
    else
      count = 0
    end if
  end function qr_local_count

  pure subroutine real_symv_column(n, a, xj, x, y, tsum)
  !Subroutine real_symv_column processes the part of one column of a symmetric
  !matrix below the diagonal: it adds xj*a to y (column contribution) and
  !returns tsum = a^T*x (mirrored row contribution) in a single pass. Four
  !independent partial sums shorten the dependency chain of the dot product.
  !
  !  Input parameters:
  !    n  - Number of elements.
  !    a  - Column part of length n.
  !    xj - Multiplier of the column.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameter:
  !    tsum - a^T*x.
    integer, intent(in) :: n
    real(wp), intent(in) :: a(n), xj, x(n)
    real(wp), intent(inout) :: y(n)
    real(wp), intent(out) :: tsum
    integer :: i, n4
    real(wp) :: t1, t2, t3, t4

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    n4 = n - mod(n, 4)
    do i = 1, n4, 4
      y(i) = y(i) + xj * a(i)
      y(i + 1) = y(i + 1) + xj * a(i + 1)
      y(i + 2) = y(i + 2) + xj * a(i + 2)
      y(i + 3) = y(i + 3) + xj * a(i + 3)
      t1 = t1 + a(i) * x(i)
      t2 = t2 + a(i + 1) * x(i + 1)
      t3 = t3 + a(i + 2) * x(i + 2)
      t4 = t4 + a(i + 3) * x(i + 3)
    end do
    do i = n4 + 1, n
      y(i) = y(i) + xj * a(i)
      t1 = t1 + a(i) * x(i)
    end do
    tsum = (t1 + t2) + (t3 + t4)
  end subroutine real_symv_column

  pure subroutine real_symv_columns4(n, a1, a2, a3, a4, x1, x2, x3, x4, x, y, &
                                     t1, t2, t3, t4)
  !Subroutine real_symv_columns4 processes the common rows of four columns of
  !a symmetric matrix below their diagonals in a single pass: it adds
  !x1*a1+x2*a2+x3*a3+x4*a4 to y (column contributions) and returns the four
  !inner products tk = ak^T*x (mirrored row contributions).
  !
  !  Input parameters:
  !    n  - Number of rows.
  !    a1, a2, a3, a4 - Column parts of length n.
  !    x1, x2, x3, x4 - Multipliers of the columns.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameters:
  !    t1, t2, t3, t4 - The four inner products.
    integer, intent(in) :: n
    real(wp), intent(in) :: a1(n), a2(n), a3(n), a4(n), x1, x2, x3, x4, x(n)
    real(wp), intent(inout) :: y(n)
    real(wp), intent(out) :: t1, t2, t3, t4
    integer :: i
    real(wp) :: xi

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    do i = 1, n
      xi = x(i)
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i) + x3 * a3(i) + x4 * a4(i)
      t1 = t1 + a1(i) * xi
      t2 = t2 + a2(i) * xi
      t3 = t3 + a3(i) * xi
      t4 = t4 + a4(i) * xi
    end do
  end subroutine real_symv_columns4

  pure subroutine real_symv_columns2(n, a1, a2, x1, x2, x, y, t1, t2)
  !Subroutine real_symv_columns2 processes the common rows of two columns of a
  !symmetric matrix below their diagonals in a single pass: it adds
  !x1*a1+x2*a2 to y (column contributions) and returns the two inner products
  !tk = ak^T*x (mirrored row contributions). Two partial sums per inner
  !product shorten the dependency chains.
  !
  !  Input parameters:
  !    n  - Number of rows.
  !    a1, a2 - Column parts of length n.
  !    x1, x2 - Multipliers of the columns.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameters:
  !    t1, t2 - The two inner products.
    integer, intent(in) :: n
    real(wp), intent(in) :: a1(n), a2(n), x1, x2, x(n)
    real(wp), intent(inout) :: y(n)
    real(wp), intent(out) :: t1, t2
    integer :: i, n2
    real(wp) :: s1, s2

    t1 = 0.0_wp
    t2 = 0.0_wp
    s1 = 0.0_wp
    s2 = 0.0_wp
    n2 = n - mod(n, 2)
    do i = 1, n2, 2
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i)
      y(i + 1) = y(i + 1) + x1 * a1(i + 1) + x2 * a2(i + 1)
      t1 = t1 + a1(i) * x(i)
      t2 = t2 + a2(i) * x(i)
      s1 = s1 + a1(i + 1) * x(i + 1)
      s2 = s2 + a2(i + 1) * x(i + 1)
    end do
    do i = n2 + 1, n
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i)
      t1 = t1 + a1(i) * x(i)
      t2 = t2 + a2(i) * x(i)
    end do
    t1 = t1 + s1
    t2 = t2 + s2
  end subroutine real_symv_columns2

  pure function real_dot(n, x, y) result(dot)
  !Function real_dot computes the real inner product x^T*y of the first n
  !elements with four independent partial sums.
  !
  !  Input parameters:
  !    n - Number of elements.
  !    x - Real vector of length n.
  !    y - Real vector of length n.
  !
  !  Result:
  !    dot - x^T*y.
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n), y(n)
    real(wp) :: dot
    integer :: i, n4
    real(wp) :: t1, t2, t3, t4

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    n4 = n - mod(n, 4)
    do i = 1, n4, 4
      t1 = t1 + x(i) * y(i)
      t2 = t2 + x(i + 1) * y(i + 1)
      t3 = t3 + x(i + 2) * y(i + 2)
      t4 = t4 + x(i + 3) * y(i + 3)
    end do
    do i = n4 + 1, n
      t1 = t1 + x(i) * y(i)
    end do
    dot = (t1 + t2) + (t3 + t4)
  end function real_dot

  pure subroutine real_axpy(n, alpha, x, y)
  !Subroutine real_axpy adds alpha*x to y for the first n elements.
  !
  !  Input parameters:
  !    n     - Number of elements.
  !    alpha - Multiplier.
  !    x     - Real vector of length n.
  !
  !  Input/output parameter:
  !    y     - Real vector of length n.
    integer, intent(in) :: n
    real(wp), intent(in) :: alpha, x(n)
    real(wp), intent(inout) :: y(n)
    integer :: i

    do i = 1, n
      y(i) = y(i) + alpha * x(i)
    end do
  end subroutine real_axpy

  pure subroutine real_add_leading_columns(n, a, lda, ncols, coef, y)
  !Subroutine real_add_leading_columns adds a linear combination of the first
  !ncols columns of a to y,
  !
  !        y(1:n) = y(1:n) + sum_k coef(k)*a(1:n,k) ,
  !
  !four columns at a time so that y is traversed once per four columns.
  !
  !  Input parameters:
  !    n     - Number of rows.
  !    a     - Matrix with leading dimension lda.
  !    lda   - Leading dimension of a.
  !    ncols - Number of columns.
  !    coef  - Coefficients, length ncols.
  !
  !  Input/output parameter:
  !    y     - Accumulated vector of length n.
    integer, intent(in) :: n, lda, ncols
    real(wp), intent(in) :: a(lda, *), coef(*)
    real(wp), intent(inout) :: y(n)
    integer :: i, k, k4
    real(wp) :: c1, c2, c3, c4

    k4 = ncols - mod(ncols, 4)
    do k = 1, k4, 4
      c1 = coef(k)
      c2 = coef(k + 1)
      c3 = coef(k + 2)
      c4 = coef(k + 3)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, k) + c2 * a(i, k + 1) + c3 * a(i, k + 2) + &
               c4 * a(i, k + 3)
      end do
    end do
    do k = k4 + 1, ncols
      c1 = coef(k)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, k)
      end do
    end do
  end subroutine real_add_leading_columns

  pure subroutine real_add_columns(n, a, lda, ncols, cols, coef, y)
  !Subroutine real_add_columns adds a linear combination of selected columns
  !of a to y,
  !
  !        y(1:n) = y(1:n) + sum_k coef(k)*a(1:n,cols(k)) ,
  !
  !four columns at a time.
  !
  !  Input parameters:
  !    n     - Number of rows.
  !    a     - Matrix with leading dimension lda.
  !    lda   - Leading dimension of a.
  !    ncols - Number of selected columns.
  !    cols  - Column indices, length ncols.
  !    coef  - Coefficients, length ncols.
  !
  !  Input/output parameter:
  !    y     - Accumulated vector of length n.
    integer, intent(in) :: n, lda, ncols, cols(*)
    real(wp), intent(in) :: a(lda, *), coef(*)
    real(wp), intent(inout) :: y(n)
    integer :: i, j1, j2, j3, j4, k, k4
    real(wp) :: c1, c2, c3, c4

    k4 = ncols - mod(ncols, 4)
    do k = 1, k4, 4
      j1 = cols(k)
      j2 = cols(k + 1)
      j3 = cols(k + 2)
      j4 = cols(k + 3)
      c1 = coef(k)
      c2 = coef(k + 1)
      c3 = coef(k + 2)
      c4 = coef(k + 3)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, j1) + c2 * a(i, j2) + c3 * a(i, j3) + &
               c4 * a(i, j4)
      end do
    end do
    do k = k4 + 1, ncols
      j1 = cols(k)
      c1 = coef(k)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, j1)
      end do
    end do
  end subroutine real_add_columns

  pure subroutine real_lu_factor(n, a, lda, piv, relative_tolerance, mag, &
                                 singular)
  !Subroutine real_lu_factor computes the LU factorization with partial
  !pivoting of a small dense real matrix in place, P*A = L*U, with unit lower
  !triangular L stored below the diagonal and the row interchanges applied to
  !whole rows (LAPACK convention).
  !
  !A pivot is rejected, and the matrix reported singular, when its magnitude
  !is zero, subnormal, or not finite. With a positive relative_tolerance a
  !pivot is also rejected when it has been lost to cancellation: the routine
  !accumulates for every element the magnitude sum
  !
  !        mag(i,j) = |a(i,j)| + sum_k |l(i,k)|*|u(k,j)|
  !
  !of the terms that were combined to produce it, and rejects a pivot whose
  !magnitude does not exceed relative_tolerance*mag. Such a pivot retains
  !fewer than about log10(1/relative_tolerance) of its digits relative to the
  !working precision. A pivot that is small without cancellation, as occurs
  !for badly scaled but well determined matrices, is accepted. Passing zero as
  !relative_tolerance disables the test and mag is not referenced.
  !
  !  Input parameters:
  !    n   - Order of the matrix; zero is allowed.
  !    lda - Leading dimension of a and mag.
  !    relative_tolerance - Relative pivot threshold (non-negative).
  !
  !  Input/output parameter:
  !    a   - On entry the matrix, on exit its factors. Undefined content is
  !          left behind when singular is true.
  !
  !  Output parameters:
  !    piv - Pivot row of every elimination step, length n.
  !    mag - Work array (lda by n) of the magnitude sums.
  !    singular - True when a pivot was rejected.
    integer, intent(in) :: n, lda
    real(wp), intent(inout) :: a(lda, *)
    integer, intent(out) :: piv(*)
    real(wp), intent(in) :: relative_tolerance
    real(wp), intent(inout) :: mag(lda, *)
    logical, intent(out) :: singular
    integer :: i, ip, j, k
    real(wp) :: factor, pivot_magnitude, swap
    logical :: track

    singular = .false.
    track = (relative_tolerance > 0.0_wp)
    if (track) then
      do j = 1, n
        do i = 1, n
          mag(i, j) = abs(a(i, j))
        end do
      end do
    end if
    do k = 1, n
      ip = k
      pivot_magnitude = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > pivot_magnitude) then
          ip = i
          pivot_magnitude = abs(a(i, k))
        end if
      end do
      piv(k) = ip
      !The negated form also rejects NaN.
      if (.not. (pivot_magnitude > tiny(1.0_wp) .and. &
                 pivot_magnitude <= huge(1.0_wp))) then
        singular = .true.
        return
      end if
      if (track) then
        if (.not. (pivot_magnitude > relative_tolerance * mag(ip, k))) then
          singular = .true.
          return
        end if
      end if
      if (ip /= k) then
        do j = 1, n
          swap = a(k, j)
          a(k, j) = a(ip, j)
          a(ip, j) = swap
        end do
        if (track) then
          do j = 1, n
            swap = mag(k, j)
            mag(k, j) = mag(ip, j)
            mag(ip, j) = swap
          end do
        end if
      end if
      do i = k + 1, n
        a(i, k) = a(i, k) / a(k, k)
      end do
      do j = k + 1, n
        factor = a(k, j)
        do i = k + 1, n
          a(i, j) = a(i, j) - a(i, k) * factor
        end do
        if (track) then
          do i = k + 1, n
            mag(i, j) = mag(i, j) + abs(a(i, k) * factor)
          end do
        end if
      end do
    end do
  end subroutine real_lu_factor

  pure subroutine real_lu_solve(n, a, lda, piv, b)
  !Subroutine real_lu_solve solves A*x = b in place with the factors computed
  !by real_lu_factor.
  !
  !  Input parameters:
  !    n   - Order of the matrix; zero is allowed.
  !    a   - LU factors with leading dimension lda.
  !    lda - Leading dimension of a.
  !    piv - Pivot rows.
  !
  !  Input/output parameter:
  !    b   - Right-hand side on entry, solution on exit (length n).
    integer, intent(in) :: n, lda
    real(wp), intent(in) :: a(lda, *)
    integer, intent(in) :: piv(*)
    real(wp), intent(inout) :: b(*)
    integer :: i, k
    real(wp) :: swap, value

    do k = 1, n
      if (piv(k) /= k) then
        swap = b(k)
        b(k) = b(piv(k))
        b(piv(k)) = swap
      end if
    end do
    do k = 1, n
      value = b(k)
      do i = k + 1, n
        b(i) = b(i) - a(i, k) * value
      end do
    end do
    do k = n, 1, -1
      b(k) = b(k) / a(k, k)
      value = b(k)
      do i = 1, k - 1
        b(i) = b(i) - a(i, k) * value
      end do
    end do
  end subroutine real_lu_solve

  pure function real_norm_squared(n, x) result(norm_squared)
  !Function real_norm_squared computes the real Euclidean inner product
  !
  !                         x^T*x = sum(x(i)^2)
  !
  !over the first n elements of x. The result is the squared norm; the square
  !root is not taken.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Real vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - x^T*x for x(1:n).
    integer, intent(in) :: n
    real(wp), intent(in) :: x(*)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + x(i) * x(i)
    end do
  end function real_norm_squared

  pure function real_direction_difference(n, x, alpha, y, x_norm_squared) &
      result(relative_difference)
  !Function real_direction_difference computes the scale-independent change of
  !direction between two real vectors,
  !
  !                 ||x-alpha*y||_2 / ||x||_2 .
  !
  !In inverse iteration alpha=(x^T*y)/(y^T*y), so alpha*y is the component of
  !x parallel to the previous iterate y. The result measures only the remaining
  !directional change and is unaffected by real rescaling of an eigenvector.
  !
  !  Input parameters:
  !    n              - Number of elements included in the calculation.
  !    x              - New real iterate, containing at least n elements.
  !    alpha          - Scalar projection coefficient multiplying y.
  !    y              - Previous real iterate, containing at least n elements.
  !    x_norm_squared - Precomputed positive value x^T*x.
  !
  !  Result:
  !    relative_difference - Relative Euclidean norm shown above.
    integer, intent(in) :: n
    real(wp), intent(in) :: x(*), alpha, y(*), x_norm_squared
    real(wp) :: relative_difference
    real(wp) :: difference, difference_norm_squared
    integer :: i

    difference_norm_squared = 0.0_wp
    do i = 1, n
      difference = x(i) - alpha * y(i)
      difference_norm_squared = difference_norm_squared + &
                                difference * difference
    end do
    relative_difference = sqrt(difference_norm_squared / x_norm_squared)
  end function real_direction_difference

  !===========================================================================
  !===========================================================================
  ! Complex Hermitian twin: type qr_complex_state
  !
  ! Everything below mirrors the real procedures above one to one (same
  ! order, same structure, names complex_* instead of real_*). The conventions
  ! that differ from the real case are listed in the module header.
  !===========================================================================
  !===========================================================================

  subroutine clear_complex_state(self)
  !Subroutine clear_complex_state releases every allocation owned by a complex
  !QR state and restores the state to its default, uninitialized condition. It
  !is the implementation of the public clear() binding and is also used before
  !reinitialization and after an allocation failure. Calling clear() for an
  !already empty state is valid.
  !
  !  Input/output parameter:
  !    self - The complex QR state. On exit all allocatable components are
  !           unallocated; all orders, counts and counters are zero; valid is
  !           false.
  !
  !The routine makes no MPI call and has no failure status. Like every other
  !procedure it must be called on all ranks so that the replicated states stay
  !identical. After clear(), initialize must be called before another
  !factorization can be constructed.
    class(qr_complex_state), intent(inout) :: self

    if (allocated(self%qt)) deallocate(self%qt)
    if (allocated(self%rt)) deallocate(self%rt)
    if (allocated(self%rd)) deallocate(self%rd)
    if (allocated(self%bidx)) deallocate(self%bidx)
    if (allocated(self%bslot)) deallocate(self%bslot)
    if (allocated(self%bpos)) deallocate(self%bpos)
    if (allocated(self%d0idx)) deallocate(self%d0idx)
    if (allocated(self%d0pos)) deallocate(self%d0pos)
    if (allocated(self%stale)) deallocate(self%stale)
    if (allocated(self%cu)) deallocate(self%cu)
    if (allocated(self%bt)) deallocate(self%bt)
    if (allocated(self%gs)) deallocate(self%gs)
    if (allocated(self%zb)) deallocate(self%zb)
    if (allocated(self%ab)) deallocate(self%ab)
    if (allocated(self%cbb)) deallocate(self%cbb)
    if (allocated(self%bg)) deallocate(self%bg)
    if (allocated(self%bz)) deallocate(self%bz)
    if (allocated(self%sc)) deallocate(self%sc)
    if (allocated(self%scc_lu)) deallocate(self%scc_lu)
    if (allocated(self%scc_mag)) deallocate(self%scc_mag)
    if (allocated(self%sig_lu)) deallocate(self%sig_lu)
    if (allocated(self%xca)) deallocate(self%xca)
    if (allocated(self%scc_piv)) deallocate(self%scc_piv)
    if (allocated(self%sig_piv)) deallocate(self%sig_piv)
    if (allocated(self%wn)) deallocate(self%wn)
    if (allocated(self%wsm)) deallocate(self%wsm)
    if (allocated(self%wloc)) deallocate(self%wloc)
    if (allocated(self%wblk)) deallocate(self%wblk)
    if (allocated(self%wred)) deallocate(self%wred)
    if (allocated(self%imark)) deallocate(self%imark)
    if (allocated(self%iw)) deallocate(self%iw)
    self%n = 0
    self%n0 = 0
    self%capacity = 0
    self%mmax = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%r_singular = .false.
    self%comm = 0
    self%nprocs = 1
    self%rank = 0
    self%nloc_cap = 0
    self%nloc = 0
    self%nblk = 0
    self%m = 0
    self%ma = 0
    self%nres = 0
    self%m0 = 0
    self%scc_ok = .true.
    self%tail_ok = .true.
    self%scc_singular = .false.
    self%sig_singular = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_complex_state

  subroutine complex_initialize(self, capacity, max_pending, info)
  !Subroutine complex_initialize prepares a complex QR state for Hermitian
  !matrices of order not greater than capacity and for at most max_pending
  !deferred (border) functions. Any factorization and counters previously held
  !by the state are discarded. The routine is collective.
  !
  !The rank and the number of processes of MPI_COMM_WORLD are obtained here and
  !fix the row distribution for the lifetime of the state. With P ranks, the
  !following storage is allocated on every rank:
  !  - two capacity by ceiling(capacity/P) arrays for the local rows of Q and
  !    of R (transposed layout);
  !  - the replicated P by P diagonal blocks of R (about P*capacity elements);
  !  - four capacity by mmax arrays for the border columns Bt, their base
  !    solves Gs, and the base solves Zb and base columns Ab of the members
  !    of D0;
  !  - small matrices of order mmax and 2*mmax and their pivot vectors;
  !  - six capacity-length work vectors and small work arrays,
  !where mmax = min(max_pending,capacity). Initialization does not form a QR
  !factorization; therefore the order is zero and valid is false on exit.
  !
  !  Input parameters:
  !    capacity    - Maximum matrix order supported by the state; positive.
  !    max_pending - Maximum number of border entries (functions whose updates
  !                  are deferred, and deleted trailing base slots); positive.
  !                  Values above capacity are reduced to capacity;
  !                  max_pending() returns the value in effect.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Existing allocations and factorization
  !            metadata are destroyed before the arguments are validated.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when all storage is ready on all ranks;
  !            QR_ERR_INVALID_ARGUMENT when capacity or max_pending is not
  !            positive;
  !            QR_ERR_INVALID_STATE when MPI has not been initialized;
  !            QR_ERR_ALLOCATION when an allocation fails on any rank.
  !
  !The allocation status is reduced with MPI_MAX over the ranks before any
  !rank branches on it, so that all ranks return the same info. If
  !initialization fails, self is returned in the empty state described for
  !clear_complex_state on every rank.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: capacity, max_pending
    integer, intent(out) :: info
    integer :: allocation_status, block_capacity, i, mpi_err, pending
    logical :: mpi_is_ready

    !Discard all previous storage before constructing the new state.
    call clear_complex_state(self)
    if (capacity <= 0 .or. max_pending <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    call MPI_Initialized(mpi_is_ready, mpi_err)
    if (.not. mpi_is_ready) then
      info = QR_ERR_INVALID_STATE
      return
    end if

    self%comm = MPI_COMM_WORLD
    call MPI_Comm_size(self%comm, self%nprocs, mpi_err)
    call MPI_Comm_rank(self%comm, self%rank, mpi_err)
    pending = min(max_pending, capacity)
    self%nloc_cap = qr_local_count(capacity, self%rank, self%nprocs)
    block_capacity = (capacity + self%nprocs - 1) / self%nprocs

    allocate(self%qt(capacity, self%nloc_cap), &
             self%rt(capacity, self%nloc_cap), &
             self%rd(self%nprocs, self%nprocs, block_capacity), &
             self%bidx(pending), self%bslot(pending), self%bpos(capacity), &
             self%d0idx(pending), self%d0pos(capacity), &
             self%stale(pending), self%cu(2 * pending), &
             self%bt(capacity, pending), self%gs(capacity, pending), &
             self%zb(capacity, pending), self%ab(capacity, pending), &
             self%cbb(pending, pending), self%bg(pending, pending), &
             self%bz(pending, pending), &
             self%sc(2 * pending, 2 * pending), &
             self%scc_lu(2 * pending, 2 * pending), &
             self%scc_mag(2 * pending, 2 * pending), &
             self%sig_lu(pending, pending), self%xca(2 * pending, pending), &
             self%scc_piv(2 * pending), self%sig_piv(pending), &
             self%wn(capacity, QR_NUM_WORK_VECTORS), &
             self%wsm(2 * pending, QR_NUM_SMALL_VECTORS), &
             self%wloc(max(1, self%nloc_cap)), self%wblk(self%nprocs), &
             self%wred(2 * pending * pending), &
             self%imark(capacity), self%iw(pending, 2), &
             stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    !A rank-local failure must become a global one before anyone branches.
    if (self%nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%qt = 0.0_wp
    self%rt = 0.0_wp
    self%rd = 0.0_wp
    self%bt = 0.0_wp
    self%gs = 0.0_wp
    self%zb = 0.0_wp
    self%ab = 0.0_wp
    self%cbb = 0.0_wp
    self%bg = 0.0_wp
    self%bz = 0.0_wp
    self%sc = 0.0_wp
    self%scc_lu = 0.0_wp
    self%scc_mag = 0.0_wp
    self%sig_lu = 0.0_wp
    self%xca = 0.0_wp
    self%scc_piv = 0
    self%sig_piv = 0
    self%wn = 0.0_wp
    self%wsm = 0.0_wp
    self%wloc = 0.0_wp
    self%wblk = 0.0_wp
    self%wred = 0.0_wp
    self%imark = 0
    self%iw = 0
    self%bidx = 0
    self%bpos = 0
    self%d0idx = 0
    self%d0pos = 0
    self%stale = 0
    self%cu = 0
    do i = 1, pending
      self%bslot(i) = i
    end do

    !Only the capacities are committed here. The order, the represented shift,
    !and the validity flag are committed by a successful factorize_fresh.
    self%capacity = capacity
    self%mmax = pending
    info = QR_SUCCESS
  end subroutine complex_initialize

  subroutine complex_factorize_fresh(self, h, s, shift, info, active_order)
  !Subroutine complex_factorize_fresh constructs the base QR factorization of
  !the shifted complex Hermitian matrix
  !
  !                     A_base = H - shift*S = Q*R
  !
  !of order n and empties the border: after a successful call the base order
  !equals the current order, nothing is pending, and no base slot is excluded.
  !This is the only O(n**3) operation of the module ("flush"). The routine is
  !collective.
  !
  !Two algorithms compute the factors and fill the same distributed layout:
  !the serial one (complex_factor_serial: LAPACK on rank zero, rows shipped to
  !their owners) and the distributed blocked Householder factorization
  !(complex_factor_distributed). complex_factor_select decides which of them
  !runs; its decision is identical on all ranks. See qr_set_factorization.
  !
  !Neither H nor S is modified or retained. When active_order is absent, H and
  !S must be square and have equal order. When active_order is present, H and
  !S may have larger physical storage and need only have at least n rows and n
  !columns. In both forms, only the lower triangles of the leading n by n
  !principal blocks are referenced; the upper triangle of the shifted matrix
  !is generated by conjugating its lower triangle. Each diagonal element of H
  !and of S must be real to within
  !
  !   QR_DIAGONAL_TOLERANCE*max(1,max(abs(Re)+abs(Im)) over the triangle) .
  !
  !H and S are tested separately, so that cancellation cannot make two
  !non-Hermitian inputs acceptable. Roundoff-sized imaginary parts of the
  !diagonals are accepted and discarded.
  !
  !  Input parameters:
  !    h     - Complex array whose leading n by n lower triangle contains the
  !            Hermitian Hamiltonian matrix H (replicated).
  !    s     - Complex array whose leading n by n lower triangle contains the
  !            Hermitian overlap matrix S (replicated).
  !    shift - The real shift represented by the factorization.
  !    active_order - Optional order n. It must be positive, must not exceed
  !                   the capacity, and must fit within both extents of H and
  !                   S. If absent, n is the common square order of H and S.
  !
  !  Input/output parameter:
  !    self  - An initialized complex QR state. On successful exit the order
  !            and the base order are n, the shift is stored, valid is true,
  !            the border and D0 are empty, and updates_since_fresh is zero.
  !            structural_updates is a lifetime counter and is not reset.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when the factorization succeeds;
  !            QR_ERR_INVALID_STATE when self has not been initialized;
  !            QR_ERR_INVALID_ARGUMENT when active_order is present but is
  !            not positive, or when a diagonal element of H or of S is
  !            not real within the tolerance above;
  !            QR_ERR_DIMENSION_MISMATCH when the original no-active-order
  !            matrix rules are violated, or when an explicit active block
  !            does not fit within H or S;
  !            QR_ERR_CAPACITY_EXCEEDED when n is greater than the capacity;
  !            QR_ERR_ALLOCATION when a rank cannot allocate the temporaries
  !            of the factorization (reported on all ranks);
  !            QR_ERR_FACTORIZATION when a LAPACK routine of the serial
  !            algorithm reports an error.
  !
  !All argument checks and the allocation of the temporaries are completed
  !before existing factors are overwritten; such an error preserves the
  !previous state. A LAPACK error leaves the state invalid (order zero).
  !An exactly singular or non-finite R does not fail here: as in the previous
  !implementation it is reported as QR_ERR_SINGULAR by the procedures that
  !need to solve with R.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer, intent(in), optional :: active_order
    integer :: blk, i, ii, j, matrix_n, np
    real(wp) :: h_scale, s_scale
    logical :: distributed

    !Check initialization, matrix shapes, active order, and capacity before
    !overwriting any component of a previously valid factorization.
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%qt)) return
    if (present(active_order)) then
      matrix_n = active_order
      if (matrix_n <= 0) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      if (matrix_n > self%capacity) then
        info = QR_ERR_CAPACITY_EXCEEDED
        return
      end if
      info = QR_ERR_DIMENSION_MISMATCH
      if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n) return
      if (size(s, 1) < matrix_n .or. size(s, 2) < matrix_n) return
    else
      info = QR_ERR_DIMENSION_MISMATCH
      if (size(h, 1) /= size(h, 2)) return
      if (size(s, 1) /= size(s, 2)) return
      if (size(h, 1) /= size(s, 1)) return
      matrix_n = size(h, 1)
      if (matrix_n <= 0) return
    end if
    if (matrix_n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if

    !Each physical diagonal must be real to within roundoff. H and S are
    !tested separately, so that cancellation cannot make two non-Hermitian
    !inputs acceptable, and only the lower triangles are read. Every rank
    !evaluates the same replicated data and takes the same branch.
    h_scale = 1.0_wp
    s_scale = 1.0_wp
    do j = 1, matrix_n
      do i = j, matrix_n
        h_scale = max(h_scale, abs(real(h(i, j), wp)) + abs(aimag(h(i, j))))
        s_scale = max(s_scale, abs(real(s(i, j), wp)) + abs(aimag(s(i, j))))
      end do
    end do
    do i = 1, matrix_n
      !The negated form also rejects NaN.
      if (.not. (abs(aimag(h(i, i))) <= QR_DIAGONAL_TOLERANCE * h_scale .and. &
                 abs(aimag(s(i, i))) <= QR_DIAGONAL_TOLERANCE * s_scale)) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
    end do

    np = self%nprocs

    !Choose the algorithm and compute the factors. Both algorithms fill
    !qt, rt, rd, nloc and nblk, and both leave the previous state untouched
    !when they fail before the factors are overwritten.
    call complex_factor_select(self%comm, np, self%rank, matrix_n, distributed)
    if (distributed) then
      call complex_factor_distributed(self, matrix_n, h, s, shift, info)
    else
      call complex_factor_serial(self, matrix_n, h, s, shift, info)
    end if
    if (info /= QR_SUCCESS) return

    !A QR factorization can complete for an exactly rank-deficient matrix.
    !Record whether the diagonal of R can be used by an unguarded triangular
    !solve. Every rank evaluates the replicated diagonal blocks and obtains
    !the same flag.
    self%r_singular = .false.
    do blk = 1, self%nblk
      do ii = 1, min(np, matrix_n - (blk - 1) * np)
        !The negated form also rejects NaN.
        if (.not. (abs(self%rd(ii, ii, blk)) > tiny(1.0_wp) .and. &
                   abs(self%rd(ii, ii, blk)) <= huge(1.0_wp))) then
          self%r_singular = .true.
        end if
      end do
    end do

    !Commit the factorization only after every numerical stage succeeds.
    self%n = matrix_n
    self%n0 = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    if (distributed) then
      self%factor_kind = QR_FACTOR_DISTRIBUTED
    else
      self%factor_kind = QR_FACTOR_SERIAL
    end if
    info = QR_SUCCESS
  end subroutine complex_factorize_fresh

  !===========================================================================
  ! Fresh factorization: choice of the algorithm, the serial algorithm, and
  ! the distributed blocked Householder algorithm with its kernels
  !===========================================================================

  pure function complex_factorization_kind(self) result(kind_used)
  !Function complex_factorization_kind reports which algorithm computed the
  !base factorization that the state holds.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    kind_used - QR_FACTOR_SERIAL or QR_FACTOR_DISTRIBUTED; zero when the
  !                state holds no valid factorization.
    class(qr_complex_state), intent(in) :: self
    integer :: kind_used

    kind_used = 0
    if (self%valid) kind_used = self%factor_kind
  end function complex_factorization_kind

  subroutine complex_factor_select(comm, nprocs, rank, n, distributed)
  !Subroutine complex_factor_select decides whether a fresh factorization of
  !order n runs the distributed algorithm. The routine is collective and its
  !result is identical on all ranks: it depends only on the requested
  !algorithm (qr_factor_resolve), on the number of ranks, and, for the
  !automatic choice with several ranks, on the timing model
  !
  !        t_serial(n)      = model(1)*n**3 + model(2)*n**2
  !        t_distributed(n) = model(3)*n**3 + model(4)*n
  !
  !whose coefficients are measured once per run by complex_factor_calibrate, on
  !rank zero, and broadcast. With one process the automatic choice is the
  !serial algorithm (LAPACK without any message).
  !
  !  Input parameters:
  !    comm, nprocs, rank - Communicator, number of ranks, and own rank.
  !    n - Order of the matrix to factorize.
  !
  !  Output parameter:
  !    distributed - True when the distributed algorithm is to be used.
    integer, intent(in) :: comm, nprocs, rank, n
    logical, intent(out) :: distributed
    integer :: choice
    real(kind(1.0d0)) :: order

    call qr_factor_resolve(comm, nprocs, rank, choice)
    if (choice == QR_FACTOR_SERIAL) then
      distributed = .false.
    else if (choice == QR_FACTOR_DISTRIBUTED) then
      distributed = .true.
    else if (nprocs == 1) then
      distributed = .false.
    else
      if (complex_factor_model(1) < 0.0d0) then
        call complex_factor_calibrate(comm, nprocs, rank)
      end if
      order = real(n, kind(1.0d0))
      distributed = (complex_factor_model(3) * order**3 + &
                     complex_factor_model(4) * order < &
                     complex_factor_model(1) * order**3 + &
                     complex_factor_model(2) * order**2)
    end if
  end subroutine complex_factor_select

  subroutine complex_factor_calibrate(comm, nprocs, rank)
  !Subroutine complex_factor_calibrate measures, on the actual machine and with
  !the actual number of ranks, how fast the two factorization algorithms are,
  !and stores the coefficients of the two timing models of
  !complex_factor_select in complex_factor_model. It follows the calibration of
  !the LDL factorizations in module linalg (linalg_ldlf_calibrate).
  !
  !Both algorithms factorize a synthetic well-conditioned matrix of two small
  !orders in a scratch state; the shortest wall time of rank zero, taken
  !between two barriers, is used. The serial algorithm is timed NUM_REPEAT
  !times and the much cheaper distributed one twice as often, in alternating
  !passes, so that a short disturbance of the node cannot spoil all timings
  !of one kind; an overestimated distributed time is the error that could
  !cost most, because it would keep large matrices on one rank.
  !
  !The cubic terms describe the arithmetic. The quadratic term of the serial
  !model is the shipping of the rows, and the linear term of the distributed
  !model is the latency of the reductions (a fixed number per column).
  !Negative coefficients, which can result from timing noise, are set to zero
  !and the other coefficient is then refitted to the larger order. Rank zero
  !fits the models and broadcasts them, so that every rank takes identical
  !decisions afterwards.
  !
  !If the scratch state cannot be allocated, the model is set such that the
  !distributed algorithm is always chosen.
  !
  !The routine is collective and must be called with more than one rank. It
  !costs a few factorizations of order QR_CALIBRATION_LARGE.
  !
  !  Input parameters:
  !    comm, nprocs, rank - Communicator, number of ranks, and own rank.
    integer, intent(in) :: comm, nprocs, rank
    integer, parameter :: NUM_REPEAT = 2
    type(qr_complex_state) :: scratch
    complex(wp), allocatable :: a(:,:)
    real(kind(1.0d0)) :: elapsed(2, 2), model(4), start
    integer :: algorithm, allocation_status, i, info, isize, j, mpi_err
    integer :: order(2), rep
    logical :: failed

    order(1) = QR_CALIBRATION_SMALL
    order(2) = QR_CALIBRATION_LARGE
    model(1) = 1.0d0
    model(2) = 0.0d0
    model(3) = 0.0d0
    model(4) = 0.0d0

    allocate(a(order(2), order(2)), stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, comm, mpi_err)
    end if
    info = QR_ERR_ALLOCATION
    if (allocation_status == 0) then
      call complex_initialize(scratch, order(2), 1, info)
    end if
    if (info == QR_SUCCESS) then
      !Synthetic Hermitian matrix with a dominant diagonal; every leading
      !block is well conditioned. It serves as H, and the shift is zero.
      do j = 1, order(2)
        do i = 1, order(2)
          a(i, j) = cmplx(1.0_wp / real(1 + abs(i - j), wp), &
                          0.25_wp * real(i - j, wp) / &
                          real((1 + abs(i - j))**2, wp), kind=wp)
        end do
        a(j, j) = a(j, j) + 3.0_wp
      end do
      failed = .false.
      elapsed = huge(1.0d0)
      do rep = 1, 2 * NUM_REPEAT
        do isize = 1, 2
          do algorithm = 1, 2
            if (algorithm == 1 .and. rep > NUM_REPEAT) cycle
            if (nprocs > 1) call MPI_Barrier(comm, mpi_err)
            start = MPI_Wtime()
            if (algorithm == 1) then
              call complex_factor_serial(scratch, order(isize), a, a, 0.0_wp, &
                                         info)
            else
              call complex_factor_distributed(scratch, order(isize), a, a, &
                                              0.0_wp, info)
            end if
            if (nprocs > 1) call MPI_Barrier(comm, mpi_err)
            elapsed(isize, algorithm) = min(elapsed(isize, algorithm), &
                                            MPI_Wtime() - start)
            if (info /= QR_SUCCESS) failed = .true.
          end do
        end do
      end do
      if (rank == 0 .and. .not. failed) then
        call fit(elapsed(1, 1), elapsed(2, 1), 2, model(1), model(2))
        call fit(elapsed(1, 2), elapsed(2, 2), 1, model(3), model(4))
      end if
      if (nprocs > 1) then
        call MPI_Bcast(model, 4, MPI_DOUBLE_PRECISION, 0, comm, mpi_err)
      end if
    end if
    call clear_complex_state(scratch)
    if (allocated(a)) deallocate(a)
    complex_factor_model = model

  contains

    subroutine fit(t1, t2, low_power, cubic, low)
    !Fits t = cubic*n**3 + low*n**low_power to the two measurements.
      real(kind(1.0d0)), intent(in) :: t1, t2
      integer, intent(in) :: low_power
      real(kind(1.0d0)), intent(out) :: cubic, low
      real(kind(1.0d0)) :: n1, n2

      n1 = real(order(1), kind(1.0d0))
      n2 = real(order(2), kind(1.0d0))
      cubic = (t2 / n2**low_power - t1 / n1**low_power) / &
              (n2**(3 - low_power) - n1**(3 - low_power))
      if (cubic < 0.0d0) then
        cubic = 0.0d0
        low = t2 / n2**low_power
      else
        low = t1 / n1**low_power - cubic * n1**(3 - low_power)
        if (low < 0.0d0) then
          low = 0.0d0
          cubic = t2 / n2**3
        end if
      end if
    end subroutine fit

  end subroutine complex_factor_calibrate

  subroutine complex_factor_serial(self, matrix_n, h, s, shift, info)
  !Subroutine complex_factor_serial computes the base factors of the shifted
  !matrix H - shift*S of order matrix_n with LAPACK on rank zero and
  !distributes them (stage A). It is called by complex_factorize_fresh, which
  !has validated the arguments and which commits the factorization, and by
  !complex_factor_calibrate.
  !
  !Stage A implementation: rank zero forms the full matrix from the two lower
  !triangles in a temporary n by n array, calls the bundled ZGEQRF, ships to
  !every rank its rows of R (row i goes to rank mod(i-1,P)), broadcasts the
  !P by P diagonal blocks of R, calls ZUNGQR to overwrite the temporary by the
  !explicit Q, ships the complex conjugates of the rows of Q, and releases its
  !temporaries (one n by n array, tau, one row buffer, and the LAPACK
  !workspace). The other ranks only receive. With one process no message is
  !sent. The diagonal of the shifted matrix is formed from the real parts of
  !the diagonals of H and S.
  !
  !  Input parameters:
  !    matrix_n - Order n; positive and not above the capacity.
  !    h, s     - Arrays whose leading n by n lower triangles hold H and S.
  !    shift    - The shift.
  !
  !  Input/output parameter:
  !    self - Initialized state. On success qt, rt, rd, nloc and nblk describe
  !           the new factors; the order, the validity flag and r_singular are
  !           left for the caller to set (the state is marked invalid). The
  !           border is emptied.
  !
  !  Output parameter:
  !    info - QR_SUCCESS; QR_ERR_ALLOCATION when rank zero cannot allocate its
  !           temporaries (state unchanged); QR_ERR_FACTORIZATION when a
  !           workspace query (state unchanged) or ZGEQRF or ZUNGQR (state
  !           invalid) reports an error. Identical on all ranks.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: matrix_n
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    complex(wp), allocatable :: a(:,:), tau(:), work(:), row(:)
    complex(wp) :: shifted_value, work_query(1)
    integer :: allocation_status, blk, dest, i, ii, j, l, lapack_info, lwork
    integer :: mpi_err, np
    logical :: is_root

    np = self%nprocs
    is_root = (self%rank == 0)

    !Rank zero allocates the temporaries of the serial factorization. The
    !status is made global before any rank branches on it.
    allocation_status = 0
    lapack_info = 0
    lwork = 1
    if (is_root) then
      allocate(a(matrix_n, matrix_n), tau(matrix_n), row(matrix_n), &
               stat=allocation_status)
      if (allocation_status == 0) then
        !LWORK=-1 performs no factorization and returns the recommended
        !workspace length in WORK(1).
        call zgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work_query, -1, &
                    lapack_info)
        if (lapack_info == 0) then
          lwork = max(matrix_n, ceiling(real(work_query(1), wp)))
          call zungqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, &
                      work_query, -1, lapack_info)
        end if
        if (lapack_info == 0) then
          lwork = max(lwork, ceiling(real(work_query(1), wp)))
          allocate(work(lwork), stat=allocation_status)
        end if
      end if
      if (allocation_status /= 0) allocation_status = 1
      !Encode a failed workspace query so that one reduction carries both.
      if (lapack_info /= 0) allocation_status = 2
    end if
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status == 1) then
      info = QR_ERR_ALLOCATION
      return
    else if (allocation_status /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !The factors are about to be overwritten. Mark the state invalid until
    !every stage has completed.
    self%valid = .false.
    self%n = 0
    self%n0 = 0
    self%nloc = 0
    self%nblk = 0
    self%r_singular = .false.
    call complex_reset_border(self)

    !Read each stored lower-triangular column contiguously and reflect the
    !complex conjugate of every off-diagonal value into the opposite triangle.
    !LAPACK requires the full shifted matrix, but callers need define only the
    !canonical lower triangle. The diagonal is formed from the real parts, so
    !that the factorized matrix is exactly Hermitian.
    lapack_info = 0
    if (is_root) then
      do j = 1, matrix_n
        a(j, j) = cmplx(real(h(j, j), wp) - shift * real(s(j, j), wp), &
                        0.0_wp, kind=wp)
        do i = j + 1, matrix_n
          shifted_value = h(i, j) - shift * s(i, j)
          a(i, j) = shifted_value
          a(j, i) = conjg(shifted_value)
        end do
      end do
      !On exit from ZGEQRF the upper triangle contains R; the strict lower
      !triangle and tau describe Q.
      call zgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work, lwork, &
                  lapack_info)
    end if
    if (np > 1) then
      call MPI_Bcast(lapack_info, 1, MPI_INTEGER, 0, self%comm, mpi_err)
    end if
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    self%nloc = qr_local_count(matrix_n, self%rank, np)
    self%nblk = (matrix_n + np - 1) / np

    !Distribute the rows of R before ZUNGQR destroys the compact reflector
    !representation. Row i is padded with zeros on the left, so that every
    !stored row has a defined triangular meaning.
    if (is_root) then
      self%rd(:, :, 1:self%nblk) = 0.0_wp
      do i = 1, matrix_n
        do j = 1, i - 1
          row(j) = 0.0_wp
        end do
        do j = i, matrix_n
          row(j) = a(i, j)
        end do
        blk = (i - 1) / np + 1
        ii = i - (blk - 1) * np
        do j = i, min(blk * np, matrix_n)
          self%rd(ii, j - (blk - 1) * np, blk) = a(i, j)
        end do
        dest = mod(i - 1, np)
        if (dest == 0) then
          self%rt(1:matrix_n, blk) = row(1:matrix_n)
        else
          call MPI_Send(row, 2 * matrix_n, MPI_WP, dest, QR_TAG_R, &
                        self%comm, &
                        mpi_err)
        end if
      end do
    else
      do l = 1, self%nloc
        call MPI_Recv(self%rt(1, l), 2 * matrix_n, MPI_WP, 0, QR_TAG_R, &
                      self%comm, MPI_STATUS_IGNORE, mpi_err)
      end do
    end if
    if (np > 1) then
      call MPI_Bcast(self%rd, 2 * np * np * self%nblk, MPI_WP, 0, &
                     self%comm, mpi_err)
    end if

    !Expand the Householder vectors and tau into the explicit unitary matrix
    !Q and distribute its rows. A rank stores the complex conjugate of each of
    !its rows, that is, a column of Q^H.
    if (is_root) then
      call zungqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, work, &
                  lwork, lapack_info)
    end if
    if (np > 1) then
      call MPI_Bcast(lapack_info, 1, MPI_INTEGER, 0, self%comm, mpi_err)
    end if
    if (lapack_info /= 0) then
      self%nloc = 0
      self%nblk = 0
      info = QR_ERR_FACTORIZATION
      return
    end if
    if (is_root) then
      do i = 1, matrix_n
        dest = mod(i - 1, np)
        if (dest == 0) then
          blk = (i - 1) / np + 1
          do j = 1, matrix_n
            self%qt(j, blk) = conjg(a(i, j))
          end do
        else
          do j = 1, matrix_n
            row(j) = conjg(a(i, j))
          end do
          call MPI_Send(row, 2 * matrix_n, MPI_WP, dest, QR_TAG_Q, &
                        self%comm, &
                        mpi_err)
        end if
      end do
    else
      do l = 1, self%nloc
        call MPI_Recv(self%qt(1, l), 2 * matrix_n, MPI_WP, 0, QR_TAG_Q, &
                      self%comm, MPI_STATUS_IGNORE, mpi_err)
      end do
    end if

    info = QR_SUCCESS
  end subroutine complex_factor_serial

  subroutine complex_factor_distributed(self, matrix_n, h, s, shift, info)
  !Subroutine complex_factor_distributed computes the base factors of the
  !shifted matrix A = H - shift*S of order n = matrix_n by a blocked
  !Householder QR factorization that works on the distributed rows (stage B).
  !It is called by complex_factorize_fresh, which has validated the arguments
  !and which commits the factorization, and by complex_factor_calibrate. It
  !fills exactly what complex_factor_serial fills.
  !
  !Layout. Row i of the matrix lives on rank mod(i-1,P) and is stored as the
  !contiguous column l = (i-1)/P+1 of the local array rt, during and after
  !the factorization. Every rank forms its rows directly from the replicated
  !lower triangles of H and S (the part of a row to the right of the diagonal
  !is the complex conjugate of a stored column, and the diagonal is formed
  !from the real parts); there is no n by n temporary on any rank.
  !
  !Reduction to R. The columns are processed in panels of QR_FACTOR_BLOCK
  !columns. A panel is copied to a compact buffer and factorized by
  !complex_panel_factor (one reduction per column). The reflectors
  !H(j) = I - tau(j)*v(j)*v(j)^H of a panel, with the columns v(j) collected
  !in V (unit diagonal), form the block reflector
  !H(j1)*...*H(j2) = I - V*T*V^H with an upper triangular T (the compact WY
  !form of LAPACK). The columns to the right of the panel are updated as
  !
  !        A := (I - V*T^H*V^H)*A = A - (V*T^H)*W ,     W = V^H*A :
  !
  !every rank multiplies its rows of V and of A (complex_block_product), one
  !reduction makes W global, and every rank subtracts its rows of V*T^H
  !times W (complex_block_update). The same reduction carries V^H*V, from which
  !T is computed. A panel therefore costs one small reduction per column and
  !one reduction of W. On exit R occupies the upper triangle; the vectors v
  !are kept below the diagonal until Q has been formed, and T is saved.
  !
  !Explicit Q. Q = H(1)*H(2)*...*H(n) is accumulated in qt by applying the
  !block reflectors in reverse order to the identity,
  !
  !        Q := (I - V*T*V^H)*Q = Q - (V*T)*W ,          W = V^H*Q ,
  !
  !with the same two kernels and one reduction per panel. Finally qt is
  !replaced by its complex conjugate (a rank stores columns of Q^H), the part
  !of rt below the diagonal is cleared, and the diagonal blocks of R are
  !collected in rd by one reduction. The diagonal of R is real, as in ZGEQRF.
  !
  !Scaling. The reflectors are generated with the conventions and the
  !safeguards of the LAPACK routine ZLARFG (complex_reflector). Sums of squares
  !are formed without scaling, which is safe because the matrix is first
  !multiplied by a power of two when its largest element is so large or so
  !small that the squares could overflow or underflow (R is scaled back at
  !the end; both operations are exact), and because a column whose part below
  !the diagonal is tiny relative to the matrix is treated separately by
  !complex_panel_factor.
  !
  !All quantities on which a rank branches (statuses, scale codes, tau, the
  !diagonal of R) come from reductions or from replicated data, so all ranks
  !take identical branches. With one process no MPI routine is called.
  !
  !  Input parameters:
  !    matrix_n - Order n; positive and not above the capacity.
  !    h, s     - Arrays whose leading n by n lower triangles hold H and S.
  !    shift    - The shift.
  !
  !  Input/output parameter:
  !    self - Initialized state. On success qt, rt, rd, nloc and nblk describe
  !           the new factors; the order, the validity flag and r_singular are
  !           left for the caller to set (the state is marked invalid). The
  !           border is emptied.
  !
  !  Output parameter:
  !    info - QR_SUCCESS, or QR_ERR_ALLOCATION when a rank cannot allocate its
  !           work arrays (about 2*(1+1/P)*QR_FACTOR_BLOCK*n elements; the
  !           state is then unchanged). Identical on all ranks.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: matrix_n
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    complex(wp), allocatable :: pan(:,:), red(:), tau(:), tsave(:,:), wbuf(:)
    complex(wp), allocatable :: work(:), ymat(:,:)
    real(wp) :: gscale
    integer :: allocation_status, code, i, j, j1, j2, jj, l, l1, la, lwork, m
    integer :: me, mpi_err, n, nb, nblk, nla, nloc, np, pb, power

    n = matrix_n
    np = self%nprocs
    me = self%rank
    nloc = qr_local_count(n, me, np)
    nb = min(QR_FACTOR_BLOCK, n)

    !Work arrays: the panel (later the explicit V of the panel) and V*T, one
    !column per local row; the T factors of all panels; the reduction
    !buffers; the transposed copies of the inner-product kernels. The status
    !is made global before any rank branches on it.
    lwork = 1
    if (QR_FACTOR_DOT_KERNELS) then
      lwork = QR_FACTOR_DOT_ROWS * nb + &
              max(16 * max(QR_FACTOR_DOT_ROWS, nb), &
                  QR_FACTOR_CACHE_BITS / storage_size(work))
    end if
    allocate(pan(nb, max(1, nloc)), ymat(nb, max(1, nloc)), tsave(nb, n), &
             tau(nb), red(2 * nb), wbuf((n + nb) * nb), work(lwork), &
             stat=allocation_status)
    if (allocation_status /= 0) allocation_status = 1
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, allocation_status, 1, MPI_INTEGER, &
                         MPI_MAX, self%comm, mpi_err)
    end if
    if (allocation_status /= 0) then
      info = QR_ERR_ALLOCATION
      return
    end if

    !The factors are about to be overwritten. Mark the state invalid until
    !every stage has completed.
    self%valid = .false.
    self%n = 0
    self%n0 = 0
    self%nloc = 0
    self%nblk = 0
    self%r_singular = .false.
    call complex_reset_border(self)

    !Form the local rows of the shifted matrix. The part of row i to the left
    !of the diagonal is row i of the stored lower triangles; the part to the
    !right is the complex conjugate of column i of the lower triangles.
    code = -huge(0)
    do l = 1, nloc
      i = me + 1 + (l - 1) * np
      do j = 1, i - 1
        self%rt(j, l) = h(i, j) - shift * s(i, j)
      end do
      self%rt(i, l) = cmplx(real(h(i, i), wp) - shift * real(s(i, i), wp), &
                            0.0_wp, kind=wp)
      do j = i + 1, n
        self%rt(j, l) = conjg(h(j, i) - shift * s(j, i))
      end do
      code = max(code, complex_scale_code(n, self%rt(1, l), 1))
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, code, 1, MPI_INTEGER, MPI_MAX, &
                         self%comm, mpi_err)
    end if

    !code is the binary exponent of the largest element. Scale the matrix by
    !a power of two when unscaled sums of squares of its elements could
    !overflow or underflow. A zero or non-finite matrix is not scaled.
    gscale = 1.0_wp
    if (abs(code) /= huge(0)) then
      if (code > QR_SCALE_EXPONENT_HIGH .or. &
          2 * code < QR_SCALE_EXPONENT_LOW2) then
        power = max(min(-code, maxexponent(1.0_wp) - 4), &
                    minexponent(1.0_wp) + 4)
        gscale = scale(1.0_wp, power)
        do l = 1, nloc
          do j = 1, n
            self%rt(j, l) = gscale * self%rt(j, l)
          end do
        end do
      end if
    end if

    !Reduction to upper triangular form, panel by panel. l1 is the first
    !local row that the panel touches and nla the number of such rows (not
    !positive when the rank owns none).
    do j1 = 1, n, nb
      pb = min(nb, n - j1 + 1)
      j2 = j1 + pb - 1
      l1 = qr_first_local(j1, me, np)
      nla = nloc - l1 + 1
      do la = 1, nla
        do jj = 1, pb
          pan(jj, la) = self%rt(j1 + jj - 1, l1 + la - 1)
        end do
      end do
      call complex_panel_factor(self%comm, np, me, j1, pb, l1, nla, pan, nb, &
                                tau, red)
      !Columns j1..j2 are final: R on and above the diagonal, the reflector
      !vectors below it.
      do la = 1, nla
        do jj = 1, pb
          self%rt(j1 + jj - 1, l1 + la - 1) = pan(jj, la)
        end do
      end do
      !Turn the panel buffer into the explicit V: unit diagonal and zeros
      !above it in the rows j1..j2.
      do jj = 1, pb
        j = j1 + jj - 1
        if (mod(j - 1, np) == me) then
          la = (j - 1 - me) / np + 2 - l1
          pan(jj, la) = 1.0_wp
          do i = jj + 1, pb
            pan(i, la) = 0.0_wp
          end do
        end if
      end do
      !W = V^H*A for the m columns to the right of the panel, followed in the
      !same buffer (leading dimension m+pb) by the pb by pb matrix V^H*V.
      m = n - j2
      do i = 1, (m + pb) * pb
        wbuf(i) = 0.0_wp
      end do
      if (nla > 0) then
        if (m > 0) then
          call complex_block_product(m, self%rt(j2 + 1, l1), self%capacity, &
                                     nla, pan, nb, pb, wbuf(1), m + pb, work)
        end if
        call complex_block_gram(pb, nla, pan, nb, wbuf(m + 1), m + pb)
      end if
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, wbuf, 2 * (m + pb) * pb, MPI_WP, &
                           MPI_SUM, self%comm, mpi_err)
      end if
      call complex_block_tfactor(pb, tau, wbuf(m + 1), m + pb, tsave(1, j1), nb)
      if (nla > 0 .and. m > 0) then
        call complex_block_vt(pb, nla, pan, nb, tsave(1, j1), nb, .true., &
                              ymat, nb)
        call complex_block_update(m, self%rt(j2 + 1, l1), self%capacity, nla, &
                                  ymat, nb, pb, wbuf(1), m + pb, work)
      end if
    end do

    !Accumulate the explicit Q: start from the identity and apply the block
    !reflectors of the panels in reverse order. A block reflector acts on the
    !rows and columns j1..n only.
    do l = 1, nloc
      do j = 1, n
        self%qt(j, l) = 0.0_wp
      end do
      self%qt(me + 1 + (l - 1) * np, l) = 1.0_wp
    end do
    do j1 = ((n - 1) / nb) * nb + 1, 1, -nb
      pb = min(nb, n - j1 + 1)
      l1 = qr_first_local(j1, me, np)
      nla = nloc - l1 + 1
      !Explicit V of the panel from the part of rt below the diagonal.
      do la = 1, nla
        l = l1 + la - 1
        i = me + 1 + (l - 1) * np
        do jj = 1, pb
          j = j1 + jj - 1
          if (i > j) then
            pan(jj, la) = self%rt(j, l)
          else if (i == j) then
            pan(jj, la) = 1.0_wp
          else
            pan(jj, la) = 0.0_wp
          end if
        end do
      end do
      m = n - j1 + 1
      do i = 1, m * pb
        wbuf(i) = 0.0_wp
      end do
      if (nla > 0) then
        call complex_block_product(m, self%qt(j1, l1), self%capacity, nla, &
                                   pan, nb, pb, wbuf(1), m, work)
      end if
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, wbuf, 2 * m * pb, MPI_WP, MPI_SUM, &
                           self%comm, mpi_err)
      end if
      if (nla > 0) then
        call complex_block_vt(pb, nla, pan, nb, tsave(1, j1), nb, .false., &
                              ymat, nb)
        call complex_block_update(m, self%qt(j1, l1), self%capacity, nla, &
                                  ymat, nb, pb, wbuf(1), m, work)
      end if
    end do

    !The reflectors were applied to Q itself. A rank stores the complex
    !conjugate of each of its rows, that is, a column of Q^H.
    do l = 1, nloc
      do j = 1, n
        self%qt(j, l) = conjg(self%qt(j, l))
      end do
    end do

    !Clear the reflector vectors stored below the diagonal of R, so that
    !every stored row has a defined triangular meaning, and undo the scaling.
    do l = 1, nloc
      i = me + 1 + (l - 1) * np
      do j = 1, i - 1
        self%rt(j, l) = 0.0_wp
      end do
      if (gscale /= 1.0_wp) then
        do j = i, n
          self%rt(j, l) = self%rt(j, l) / gscale
        end do
      end if
    end do

    !Replicated diagonal blocks of R: every rank owns one row of each block
    !(its local row l belongs to block l); one reduction, in which each
    !element is the sum of one value and zeros, completes them.
    nblk = (n + np - 1) / np
    self%rd(:, :, 1:nblk) = 0.0_wp
    do l = 1, nloc
      do j = me + 1 + (l - 1) * np, min(l * np, n)
        self%rd(me + 1, j - (l - 1) * np, l) = self%rt(j, l)
      end do
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, self%rd, 2 * np * np * nblk, MPI_WP, &
                         MPI_SUM, self%comm, mpi_err)
    end if

    self%nloc = nloc
    self%nblk = nblk
    info = QR_SUCCESS
  end subroutine complex_factor_distributed

  subroutine complex_panel_factor(comm, np, me, j1, pb, l1, nla, pan, ldp, &
                                  tau, red)
  !Subroutine complex_panel_factor reduces one panel of the distributed matrix
  !to upper triangular form by Householder reflections. The panel consists of
  !the global columns j1..j1+pb-1 and of the rows j1..n; a rank holds its rows
  !of the panel as the columns of the compact buffer pan.
  !
  !Column j is processed as follows.
  !  1. Every rank accumulates, over its rows below the diagonal,
  !     z(c) = sum_i conjg(a(i,j))*a(i,c) for the panel columns c >= j (z(j) is
  !     the squared norm of the part of column j below the diagonal). The owner
  !     of row j adds the pivot row a(j,c). One reduction makes both global.
  !  2. Every rank generates the reflector from the pivot element and the
  !     norm (complex_reflector): the new real diagonal element beta, tau, and
  !     the scale factor of the vector v = (1, x/(alpha-beta)). The products
  !     v^H*a(:,c) follow from z and the pivot row without further
  !     communication.
  !  3. Every rank scales its part of v in place and applies the reflector to
  !     its rows of the panel columns to the right of j.
  !
  !A column whose squared norm below the diagonal is smaller than
  !QR_SMALL_NORM2 (in particular zero) or not finite is treated separately,
  !because its unscaled sum of squares is inaccurate: one integer reduction
  !gives the binary exponent of its largest element, and step 1 is repeated
  !with the column scaled by a power of two. This reproduces the safeguards
  !of ZLARFG at the cost of two more reductions for such a column only.
  !
  !  Input parameters:
  !    comm, np, me - Communicator, number of ranks, and own rank.
  !    j1   - Global index of the first column of the panel.
  !    pb   - Number of columns of the panel.
  !    l1   - Local position of the first local row of the panel (global
  !           index not below j1).
  !    nla  - Number of local rows of the panel; may be zero or negative.
  !    ldp  - Leading dimension of pan, at least pb.
  !
  !  Input/output parameter:
  !    pan  - pan(jj,la) is the element in column j1+jj-1 of the local row
  !           l1+la-1. On exit the upper triangle of the panel holds R and
  !           the part below the diagonal the reflector vectors (without
  !           their unit diagonal elements).
  !
  !  Output parameter:
  !    tau  - The pb scalar factors of the reflectors; replicated.
  !
  !  Workspace:
  !    red  - Reduction buffer of length 2*pb.
    integer, intent(in) :: comm, np, me, j1, pb, l1, nla, ldp
    complex(wp), intent(inout) :: pan(ldp, *)
    complex(wp), intent(out) :: tau(pb)
    complex(wp), intent(inout) :: red(2 * pb)
    complex(wp) :: alpha, factor2, vl
    real(wp) :: beta, factor1, factor3, sc, xnorm2
    integer :: c, code, j, jj, la, lgt, lown, mpi_err, nc, nrows, power
    logical :: owner

    do jj = 1, pb
      j = j1 + jj - 1
      nc = pb - jj + 1
      !The local rows below the diagonal are the positions lgt..nla; the
      !owner of row j holds it at position lown.
      lgt = qr_first_local(j + 1, me, np) - l1 + 1
      nrows = max(nla - lgt + 1, 0)
      owner = (mod(j - 1, np) == me)
      lown = (j - 1 - me) / np + 2 - l1

      sc = 1.0_wp
      call accumulate()
      xnorm2 = real(red(1), wp)
      !The negated form also takes NaN to the separate treatment.
      if (.not. (xnorm2 >= QR_SMALL_NORM2 .and. &
                 xnorm2 <= huge(1.0_wp))) then
        code = -huge(0)
        if (nrows > 0) code = complex_scale_code(nrows, pan(jj, lgt), ldp)
        if (np > 1) then
          call MPI_Allreduce(MPI_IN_PLACE, code, 1, MPI_INTEGER, MPI_MAX, &
                             comm, mpi_err)
        end if
        if (code == -huge(0)) then
          !The column is exactly zero below the diagonal.
          xnorm2 = 0.0_wp
        else if (code /= huge(0)) then
          power = max(min(-code, maxexponent(1.0_wp) - 4), &
                      minexponent(1.0_wp) + 4)
          sc = scale(1.0_wp, power)
          call accumulate()
          xnorm2 = real(red(1), wp)
        end if
      end if

      alpha = red(nc + 1)
      call complex_reflector(alpha, xnorm2, sc, beta, tau(jj), factor1, &
                             factor2, factor3)
      !red(c) becomes conjg(tau)*(v^H*a(:,c)), the row that the adjoint of the
      !reflector subtracts
      !from the pivot row, for the panel columns to the right of j.
      do c = 2, nc
        red(c) = conjg(tau(jj)) * (red(nc + c) + &
                                   (conjg(factor2) * red(c)) * factor3)
      end do
      do la = lgt, nla
        vl = (pan(jj, la) * factor1) * factor2
        pan(jj, la) = vl
        do c = 2, nc
          pan(jj + c - 1, la) = pan(jj + c - 1, la) - vl * red(c)
        end do
      end do
      if (owner) then
        pan(jj, lown) = beta
        do c = 2, nc
          pan(jj + c - 1, lown) = pan(jj + c - 1, lown) - red(c)
        end do
      end if
    end do

  contains

    subroutine accumulate()
    !Step 1 for the current column with the scale factor sc; on exit red
    !holds z (nc values) followed by the pivot row (nc values), reduced.
      integer :: k

      if (nrows > 0) then
        call complex_panel_column(nc, nrows, pan(jj, lgt), ldp, sc, red)
      else
        do k = 1, nc
          red(k) = 0.0_wp
        end do
      end if
      do k = 1, nc
        if (owner) then
          red(nc + k) = pan(jj + k - 1, lown)
        else
          red(nc + k) = 0.0_wp
        end if
      end do
      if (np > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, red, 4 * nc, MPI_WP, MPI_SUM, comm, &
                           mpi_err)
      end if
    end subroutine accumulate

  end subroutine complex_panel_factor

  pure subroutine complex_panel_column(nc, nrows, x, ldx, sc, z)
  !Subroutine complex_panel_column computes, for one column of a panel and the
  !local rows below its diagonal element, the products of the scaled column
  !with itself and with the panel columns to its right:
  !
  !        z(1) = sum_i abs(sc*x(1,i))**2 ,
  !        z(c) = sum_i sc*conjg(x(1,i))*x(c,i) ,   c = 2..nc ,
  !
  !where x(c,i) is the element of local row i in the c-th column counted from
  !the current one.
  !
  !  Input parameters:
  !    nc    - Number of panel columns from the current one to the last.
  !    nrows - Number of local rows below the diagonal; positive.
  !    x     - Rows of the panel, stored as columns, starting at the current
  !            panel column; leading dimension ldx.
  !    ldx   - Leading dimension of x.
  !    sc    - Scale factor of the current column (a power of two).
  !
  !  Output parameter:
  !    z     - The nc local sums.
    integer, intent(in) :: nc, nrows, ldx
    complex(wp), intent(in) :: x(ldx, *)
    real(wp), intent(in) :: sc
    complex(wp), intent(out) :: z(nc)
    integer :: c, i
    complex(wp) :: xs

    do c = 1, nc
      z(c) = 0.0_wp
    end do
    do i = 1, nrows
      xs = conjg(x(1, i)) * sc
      z(1) = z(1) + (real(xs, wp)**2 + aimag(xs)**2)
      do c = 2, nc
        z(c) = z(c) + xs * x(c, i)
      end do
    end do
  end subroutine complex_panel_column

  pure function complex_scale_code(n, x, incx) result(code)
  !Function complex_scale_code describes the magnitude of the largest element
  !of a vector by an integer, which can be reduced over the ranks with MPI_MAX
  !without sending floating-point data:
  !
  !        huge(0)     when an element is not finite;
  !        -huge(0)    when all elements are zero (or n is zero);
  !        otherwise the binary exponent of the largest absolute value of a
  !        real or imaginary part.
  !
  !  Input parameters:
  !    n    - Number of elements.
  !    x    - Vector with increment incx.
  !    incx - Increment, positive.
  !
  !  Result:
  !    code - The integer described above.
    integer, intent(in) :: n, incx
    complex(wp), intent(in) :: x(*)
    integer :: code
    integer :: i, k
    real(wp) :: largest, magnitude
    logical :: finite

    largest = 0.0_wp
    finite = .true.
    k = 1
    do i = 1, n
      magnitude = max(abs(real(x(k), wp)), abs(aimag(x(k))))
      if (magnitude > largest) largest = magnitude
      !The negated form also detects NaN.
      if (.not. (magnitude <= huge(1.0_wp))) finite = .false.
      k = k + incx
    end do
    if (.not. finite) then
      code = huge(0)
    else if (largest == 0.0_wp) then
      code = -huge(0)
    else
      code = exponent(largest)
    end if
  end function complex_scale_code

  pure subroutine complex_reflector(alpha, xnorm2, sc, beta, tau, factor1, &
                                    factor2, factor3)
  !Subroutine complex_reflector generates an elementary Householder reflector
  !H = I - tau*v*v^H with v = (1, v2) such that H^H*(alpha, x) = (beta, 0),
  !with the conventions of the LAPACK routine ZLARFG: beta is real and has
  !the sign opposite to the real part of alpha, tau = (beta-alpha)/beta is
  !complex, v2 = x/(alpha-beta), and H = I (tau = 0) when x is zero and alpha
  !is real. The diagonal of R is therefore real.
  !
  !The vector x is distributed and is not an argument. The caller supplies
  !the squared norm of sc*x, where sc is a power of two (one unless the
  !column needed scaling), and receives the factors with which it completes
  !the operation:
  !
  !        v2 = (x*factor1)*factor2 ,
  !        v^H*a(:,c) = a(j,c) + (conjg(factor2)*z(c))*factor3 ,
  !
  !where z(c) = sum_i sc*conjg(x(i))*a(i,c). When alpha and x are both so small
  !that beta is below QR_SAFE_MINIMUM, the reflector is generated from the
  !scaled quantities, as ZLARFG does. Moduli are formed with HYPOT and the
  !reciprocal of alpha-beta with complex_reciprocal, so that no real or
  !imaginary part is squared.
  !
  !  Input parameters:
  !    alpha  - The pivot element.
  !    xnorm2 - Squared Euclidean norm of sc*x.
  !    sc     - Scale factor of x.
  !
  !  Output parameters:
  !    beta   - The new diagonal element (real).
  !    tau    - The scalar factor of the reflector.
  !    factor1, factor2, factor3 - See above.
    complex(wp), intent(in) :: alpha
    real(wp), intent(in) :: xnorm2, sc
    complex(wp), intent(out) :: tau, factor2
    real(wp), intent(out) :: beta, factor1, factor3
    complex(wp) :: alpha_scaled
    real(wp) :: beta_scaled, norm, xnorm, xnorm_scaled

    xnorm_scaled = sqrt(xnorm2)
    xnorm = xnorm_scaled / sc
    if (xnorm == 0.0_wp .and. aimag(alpha) == 0.0_wp) then
      beta = real(alpha, wp)
      tau = 0.0_wp
      factor1 = 1.0_wp
      factor2 = 0.0_wp
      factor3 = 0.0_wp
      return
    end if
    norm = hypot(hypot(real(alpha, wp), aimag(alpha)), xnorm)
    if (norm >= QR_SAFE_MINIMUM) then
      beta = -sign(norm, real(alpha, wp))
      tau = cmplx((beta - real(alpha, wp)) / beta, -aimag(alpha) / beta, &
                  kind=wp)
      factor1 = 1.0_wp
      factor2 = complex_reciprocal(alpha - beta)
      factor3 = 1.0_wp / sc
    else
      alpha_scaled = alpha * sc
      norm = hypot(hypot(real(alpha_scaled, wp), aimag(alpha_scaled)), &
                   xnorm_scaled)
      beta_scaled = -sign(norm, real(alpha_scaled, wp))
      tau = cmplx((beta_scaled - real(alpha_scaled, wp)) / beta_scaled, &
                  -aimag(alpha_scaled) / beta_scaled, kind=wp)
      factor1 = sc
      factor2 = complex_reciprocal(alpha_scaled - beta_scaled)
      factor3 = 1.0_wp
      beta = beta_scaled / sc
    end if
  end subroutine complex_reflector

  pure subroutine complex_block_product(m, b, ldb, ncols, v, ldv, k, w, ldw, &
                                        work)
  !Subroutine complex_block_product adds to w the product of the local rows of
  !a distributed matrix with the local rows of the reflector block V,
  !
  !        w(i,jj) = w(i,jj) + sum_l b(i,l)*conjg(v(jj,l)) ,
  !
  !i = 1..m, jj = 1..k, l = 1..ncols. Column l of b is a local row of the
  !matrix (transposed layout) and column l of v the same local row of V, so
  !that w, summed over the ranks, is the transpose of V^H*A.
  !
  !The rows i are processed in chunks, such that the part of w in work fits
  !in the cache, and four columns of b at a time, so that a column of w is
  !read and written once per four columns of b.
  !
  !With QR_FACTOR_DOT_KERNELS the product is formed from inner products
  !instead: blocks of QR_FACTOR_DOT_ROWS local rows of b and v are copied
  !transposed to work, chunk by chunk, and complex_dot_multiply accumulates the
  !sums over the local rows in registers.
  !
  !  Input parameters:
  !    m     - Number of rows of b and w (columns of the distributed matrix).
  !    b     - Local rows of the matrix; leading dimension ldb.
  !    ldb   - Leading dimension of b.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V: v(jj,l) is the element of reflector jj in
  !            local row l; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    k     - Number of reflectors of the block.
  !    ldw   - Leading dimension of w.
  !
  !  Input/output parameter:
  !    w     - Accumulated product, m by k.
  !
  !  Workspace:
  !    work  - Referenced only with QR_FACTOR_DOT_KERNELS; length as allocated
  !            by complex_factor_distributed.
    integer, intent(in) :: m, ldb, ncols, ldv, k, ldw
    complex(wp), intent(in) :: b(ldb, *), v(ldv, *)
    complex(wp), intent(inout) :: w(ldw, *)
    complex(wp), intent(inout) :: work(*)
    integer :: chunk, i, i0, i1, jj, l, l0, l4, nl
    complex(wp) :: c1, c2, c3, c4

    if (QR_FACTOR_DOT_KERNELS) then
      !work(1:nl*k) holds the block of v transposed (nl by k), and the rest
      !a chunk of the block of b transposed (nl by at most chunk).
      chunk = max(16, QR_FACTOR_CACHE_BITS / &
                      (QR_FACTOR_DOT_ROWS * storage_size(c1)))
      do l0 = 1, ncols, QR_FACTOR_DOT_ROWS
        nl = min(QR_FACTOR_DOT_ROWS, ncols - l0 + 1)
        do jj = 1, k
          do l = 1, nl
            work(l + (jj - 1) * nl) = conjg(v(jj, l0 + l - 1))
          end do
        end do
        do i0 = 1, m, chunk
          i1 = min(i0 + chunk - 1, m)
          do l = 1, nl
            do i = i0, i1
              work(nl * k + l + (i - i0) * nl) = b(i, l0 + l - 1)
            end do
          end do
          call complex_dot_multiply(nl, i1 - i0 + 1, k, work(nl * k + 1), nl, &
                                    work(1), nl, 1.0_wp, w(i0, 1), ldw)
        end do
      end do
      return
    end if

    chunk = max(16, QR_FACTOR_CACHE_BITS / (k * storage_size(c1)))
    l4 = ncols - mod(ncols, 4)
    do i0 = 1, m, chunk
      i1 = min(i0 + chunk - 1, m)
      do l = 1, l4, 4
        do jj = 1, k
          c1 = conjg(v(jj, l))
          c2 = conjg(v(jj, l + 1))
          c3 = conjg(v(jj, l + 2))
          c4 = conjg(v(jj, l + 3))
          do i = i0, i1
            w(i, jj) = w(i, jj) + c1 * b(i, l) + c2 * b(i, l + 1) + &
                       c3 * b(i, l + 2) + c4 * b(i, l + 3)
          end do
        end do
      end do
      do l = l4 + 1, ncols
        do jj = 1, k
          c1 = conjg(v(jj, l))
          do i = i0, i1
            w(i, jj) = w(i, jj) + c1 * b(i, l)
          end do
        end do
      end do
    end do
  end subroutine complex_block_product

  pure subroutine complex_block_update(m, b, ldb, ncols, y, ldy, k, w, ldw, &
                                       work)
  !Subroutine complex_block_update subtracts from the local rows of a
  !distributed matrix the product of the local rows of Y = V*T^H (or V*T) with
  !the replicated matrix W,
  !
  !        b(i,l) = b(i,l) - sum_jj w(i,jj)*y(jj,l) ,
  !
  !i = 1..m, l = 1..ncols, jj = 1..k. The rows i are processed in chunks, as
  !in complex_block_product, and eight columns of w at a time (EIGHT_COLUMNS;
  !for complex data in double precision eight were measured to be clearly
  !faster than four).
  !
  !With QR_FACTOR_DOT_KERNELS every chunk of w is copied transposed to work,
  !and complex_dot_multiply forms the sums over jj as inner products.
  !
  !  Input parameters:
  !    m     - Number of rows of b and w (columns of the distributed matrix).
  !    ldb   - Leading dimension of b.
  !    ncols - Number of local rows; positive.
  !    y     - Local rows of Y: y(jj,l) is the element of column jj of Y in
  !            local row l; leading dimension ldy.
  !    ldy   - Leading dimension of y.
  !    k     - Number of reflectors of the block.
  !    w     - Replicated m by k matrix; leading dimension ldw.
  !    ldw   - Leading dimension of w.
  !
  !  Input/output parameter:
  !    b     - Local rows of the matrix; leading dimension ldb.
  !
  !  Workspace:
  !    work  - Referenced only with QR_FACTOR_DOT_KERNELS; length as allocated
  !            by complex_factor_distributed.
    integer, intent(in) :: m, ldb, ncols, ldy, k, ldw
    complex(wp), intent(in) :: y(ldy, *), w(ldw, *)
    complex(wp), intent(inout) :: b(ldb, *)
    complex(wp), intent(inout) :: work(*)
    logical, parameter :: EIGHT_COLUMNS = .true.
    integer :: chunk, i, i0, i1, jj, k4, k8, l
    complex(wp) :: c1, c2, c3, c4, c5, c6, c7, c8

    chunk = max(16, QR_FACTOR_CACHE_BITS / (k * storage_size(c1)))
    if (QR_FACTOR_DOT_KERNELS) then
      do i0 = 1, m, chunk
        i1 = min(i0 + chunk - 1, m)
        do jj = 1, k
          do i = i0, i1
            work(jj + (i - i0) * k) = w(i, jj)
          end do
        end do
        call complex_dot_multiply(k, i1 - i0 + 1, ncols, work(1), k, y, ldy, &
                                  -1.0_wp, b(i0, 1), ldb)
      end do
      return
    end if

    k8 = 0
    if (EIGHT_COLUMNS) k8 = k - mod(k, 8)
    k4 = k - mod(k - k8, 4)
    do i0 = 1, m, chunk
      i1 = min(i0 + chunk - 1, m)
      do l = 1, ncols
        do jj = 1, k8, 8
          c1 = y(jj, l)
          c2 = y(jj + 1, l)
          c3 = y(jj + 2, l)
          c4 = y(jj + 3, l)
          c5 = y(jj + 4, l)
          c6 = y(jj + 5, l)
          c7 = y(jj + 6, l)
          c8 = y(jj + 7, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj) - c2 * w(i, jj + 1) - &
                      c3 * w(i, jj + 2) - c4 * w(i, jj + 3) - &
                      c5 * w(i, jj + 4) - c6 * w(i, jj + 5) - &
                      c7 * w(i, jj + 6) - c8 * w(i, jj + 7)
          end do
        end do
        do jj = k8 + 1, k4, 4
          c1 = y(jj, l)
          c2 = y(jj + 1, l)
          c3 = y(jj + 2, l)
          c4 = y(jj + 3, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj) - c2 * w(i, jj + 1) - &
                      c3 * w(i, jj + 2) - c4 * w(i, jj + 3)
          end do
        end do
        do jj = k4 + 1, k
          c1 = y(jj, l)
          do i = i0, i1
            b(i, l) = b(i, l) - c1 * w(i, jj)
          end do
        end do
      end do
    end do
  end subroutine complex_block_update

  pure subroutine complex_dot_multiply(n, ma, mb, a, lda, b, ldb, alpha, c, ldc)
  !Subroutine complex_dot_multiply adds a multiple of the product of the plain
  !(unconjugated) transpose of one matrix with another matrix to a third,
  !
  !        c(i,j) = c(i,j) + alpha * sum_l a(l,i)*b(l,j) ,
  !
  !i = 1..ma, j = 1..mb, l = 1..n (no element is conjugated). Every element of
  !the product is an inner product of two contiguous columns; four inner
  !products with the same column of a are accumulated in one pass, so that
  !each element of c is read and written only once.
  !
  !  Input parameters:
  !    n     - Number of rows of a and b.
  !    ma    - Number of columns of a (rows of c).
  !    mb    - Number of columns of b (columns of c).
  !    a     - n by ma matrix; leading dimension lda.
  !    lda   - Leading dimension of a.
  !    b     - n by mb matrix; leading dimension ldb.
  !    ldb   - Leading dimension of b.
  !    alpha - Real multiplier of the product.
  !    ldc   - Leading dimension of c.
  !
  !  Input/output parameter:
  !    c     - ma by mb matrix; leading dimension ldc.
    integer, intent(in) :: n, ma, mb, lda, ldb, ldc
    complex(wp), intent(in) :: a(lda, *), b(ldb, *)
    real(wp), intent(in) :: alpha
    complex(wp), intent(inout) :: c(ldc, *)
    integer :: i, j, l, mb4
    complex(wp) :: t1, t2, t3, t4, x

    mb4 = mb - mod(mb, 4)
    do j = 1, mb4, 4
      do i = 1, ma
        t1 = 0.0_wp
        t2 = 0.0_wp
        t3 = 0.0_wp
        t4 = 0.0_wp
        do l = 1, n
          x = a(l, i)
          t1 = t1 + x * b(l, j)
          t2 = t2 + x * b(l, j + 1)
          t3 = t3 + x * b(l, j + 2)
          t4 = t4 + x * b(l, j + 3)
        end do
        c(i, j) = c(i, j) + alpha * t1
        c(i, j + 1) = c(i, j + 1) + alpha * t2
        c(i, j + 2) = c(i, j + 2) + alpha * t3
        c(i, j + 3) = c(i, j + 3) + alpha * t4
      end do
    end do
    do j = mb4 + 1, mb
      do i = 1, ma
        t1 = 0.0_wp
        do l = 1, n
          t1 = t1 + a(l, i) * b(l, j)
        end do
        c(i, j) = c(i, j) + alpha * t1
      end do
    end do
  end subroutine complex_dot_multiply

  pure subroutine complex_block_gram(k, ncols, v, ldv, g, ldg)
  !Subroutine complex_block_gram adds to g the local part of the strict upper
  !triangle of V^H*V,
  !
  !        g(a,b) = g(a,b) + sum_l conjg(v(a,l))*v(b,l) ,   a < b ,
  !
  !from which complex_block_tfactor computes the triangular factor T.
  !
  !  Input parameters:
  !    k     - Number of reflectors of the block.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    ldg   - Leading dimension of g.
  !
  !  Input/output parameter:
  !    g     - Accumulated strict upper triangle, k by k.
    integer, intent(in) :: k, ncols, ldv, ldg
    complex(wp), intent(in) :: v(ldv, *)
    complex(wp), intent(inout) :: g(ldg, *)
    integer :: a, b, l
    complex(wp) :: vb

    do l = 1, ncols
      do b = 2, k
        vb = v(b, l)
        do a = 1, b - 1
          g(a, b) = g(a, b) + conjg(v(a, l)) * vb
        end do
      end do
    end do
  end subroutine complex_block_gram

  pure subroutine complex_block_tfactor(k, tau, g, ldg, t, ldt)
  !Subroutine complex_block_tfactor forms the upper triangular factor T of the
  !block reflector H(1)*...*H(k) = I - V*T*V^H from the scalar factors tau
  !and the strict upper triangle of G = V^H*V, by the recurrence of the
  !LAPACK routine ZLARFT (forward direction, columnwise storage):
  !
  !        T(1:b-1,b) = -tau(b)*T(1:b-1,1:b-1)*G(1:b-1,b) ,   T(b,b) = tau(b) .
  !
  !  Input parameters:
  !    k   - Number of reflectors.
  !    tau - Scalar factors.
  !    g   - Strict upper triangle of V^H*V; leading dimension ldg.
  !    ldg - Leading dimension of g.
  !    ldt - Leading dimension of t.
  !
  !  Output parameter:
  !    t   - The k by k factor; the part below the diagonal is set to zero.
    integer, intent(in) :: k, ldg, ldt
    complex(wp), intent(in) :: tau(k), g(ldg, *)
    complex(wp), intent(out) :: t(ldt, *)
    integer :: a, a2, b
    complex(wp) :: tsum

    do b = 1, k
      do a = 1, b - 1
        t(a, b) = -tau(b) * g(a, b)
      end do
      !Multiply by the leading triangle in place; row a needs only the
      !elements a..b-1 of the column, which are still unchanged.
      do a = 1, b - 1
        tsum = 0.0_wp
        do a2 = a, b - 1
          tsum = tsum + t(a, a2) * t(a2, b)
        end do
        t(a, b) = tsum
      end do
      t(b, b) = tau(b)
      do a = b + 1, k
        t(a, b) = 0.0_wp
      end do
    end do
  end subroutine complex_block_tfactor

  pure subroutine complex_block_vt(k, ncols, v, ldv, t, ldt, adjoint, y, ldy)
  !Subroutine complex_block_vt multiplies the local rows of the reflector block
  !V by the triangular factor of the block reflector,
  !
  !        Y = V*T^H   (adjoint true:  y(jj,l) = sum_{a>=jj}
  !        conjg(t(jj,a))*v(a,l)) Y = V*T     (adjoint false: y(jj,l) =
  !        sum_{a<=jj} t(a,jj)*v(a,l)).
  !
  !The first form belongs to the reduction of the matrix (application of the
  !adjoint of the block reflector), the second to the accumulation of Q.
  !
  !  Input parameters:
  !    k     - Number of reflectors of the block.
  !    ncols - Number of local rows; positive.
  !    v     - Local rows of V; leading dimension ldv.
  !    ldv   - Leading dimension of v.
  !    t     - Upper triangular factor; leading dimension ldt.
  !    ldt   - Leading dimension of t.
  !    adjoint - Selects the form.
  !    ldy   - Leading dimension of y.
  !
  !  Output parameter:
  !    y     - Local rows of Y, k by ncols.
    integer, intent(in) :: k, ncols, ldv, ldt, ldy
    complex(wp), intent(in) :: v(ldv, *), t(ldt, *)
    logical, intent(in) :: adjoint
    complex(wp), intent(out) :: y(ldy, *)
    integer :: a, jj, l
    complex(wp) :: tsum

    do l = 1, ncols
      if (adjoint) then
        do jj = 1, k
          tsum = 0.0_wp
          do a = jj, k
            tsum = tsum + conjg(t(jj, a)) * v(a, l)
          end do
          y(jj, l) = tsum
        end do
      else
        do jj = 1, k
          tsum = 0.0_wp
          do a = 1, jj
            tsum = tsum + t(a, jj) * v(a, l)
          end do
          y(jj, l) = tsum
        end do
      end if
    end do
  end subroutine complex_block_vt

  subroutine complex_begin_active(self, indices, need_flush, info)
  !Subroutine complex_begin_active declares the set of functions that the
  !following trials will modify (the ACTIVE set). The module never guesses the
  !active set: every index passed later to replace_symmetric must have been
  !announced here. The routine is collective.
  !
  !The announced indices are either existing functions (index not greater than
  !order()) or the consecutive range order()+1, order()+2, ... of functions
  !that append_symmetric will supply; both kinds may be mixed, in any order.
  !
  !  1. If the announced set equals the current active set (the entries of
  !     the active tail plus the announced appended indices beyond the base
  !     that are still unsupplied), nothing happens.
  !  2. Otherwise the border must be able to hold the request. When it cannot
  !     (pending_count() plus the new entries would exceed max_pending()),
  !     need_flush is returned true and the state is unchanged: the caller
  !     calls factorize_fresh and then begin_active again.
  !  3. The previous tail is committed and the announced indices that already
  !     have a border entry are moved to the end of the border list. An
  !     announced existing function without an entry is necessarily a clean
  !     base slot; it joins D0 with a zero modification: one base solve for
  !     its column of Zb, its base column recovered from the factors as
  !     Q*(R*e), the entries of the other border columns on its row moved
  !     into the border-border block, and one distributed rebuild of the
  !     small matrices.
  !  4. If the committed set changed, LU(Scc) is recomputed. This is the only
  !     place, besides the lazy refresh after truncate changed a committed
  !     entry, where O(m**3) work on the committed block happens.
  !
  !The represented matrix does not change. The tail columns stay as they are
  !until replace_symmetric or append_symmetric supplies new ones. Announcing
  !an appended range reserves border entries; it does not change order().
  !
  !Cost: nothing in case 1; O(m) bookkeeping plus LU(Scc) when only the
  !partition changes; for each newly pending base function one base solve and
  !one column recovery, plus one distributed O(m**2*n0/P) rebuild.
  !
  !  Input parameter:
  !    indices - Distinct canonical indices of the active set. The size may be
  !              zero, which commits the current tail.
  !
  !  Input/output parameter:
  !    self    - A valid complex QR state.
  !
  !  Output parameters:
  !    need_flush - True when the request does not fit in the border; the
  !                 state is then unchanged and info is QR_SUCCESS.
  !    info    - QR_SUCCESS;
  !              QR_ERR_INVALID_STATE when self holds no valid factorization;
  !              QR_ERR_INVALID_ARGUMENT when an index is not positive, when
  !              indices are repeated, or when the indices above order() do
  !              not form the range order()+1..order()+k;
  !              QR_ERR_CAPACITY_EXCEEDED when order()+k exceeds the capacity;
  !              QR_ERR_BATCH_FULL when the set itself is larger than
  !              max_pending() and can therefore never be held (the caller
  !              applies such a trial by a fresh factorization);
  !              QR_ERR_SINGULAR when R or the committed block Scc is
  !              numerically singular. The border bookkeeping is then
  !              consistent, but the state cannot solve with this partition;
  !              the caller refactorizes.
  !
  !Argument errors leave the state unchanged.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: indices(:)
    logical, intent(out) :: need_flush
    integer, intent(out) :: info
    integer :: count_active, count_append, count_beyond, count_entry
    integer :: count_new, i, idx, k, mc, p
    logical :: committed_changed, same_set

    need_flush = .false.
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return

    !Validate the indices and classify them.
    count_active = size(indices)
    count_append = 0
    do i = 1, count_active
      if (indices(i) < 1) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      if (indices(i) > self%n) count_append = count_append + 1
    end do
    do i = 1, count_active
      if (indices(i) > self%n + count_append) then
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
    end do
    if (self%n + count_append > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    do i = 1, count_active
      idx = indices(i)
      if (self%imark(idx) /= 0) then
        !Repeated index. Remove the marks set so far.
        do k = 1, i - 1
          self%imark(indices(k)) = 0
        end do
        info = QR_ERR_INVALID_ARGUMENT
        return
      end if
      self%imark(idx) = 1
    end do
    !From here on every return path must clear the marks.

    if (count_active > self%mmax) then
      call clear_marks()
      info = QR_ERR_BATCH_FULL
      return
    end if

    !Compare with the current active set. An announced index is beyond the
    !base (pure reservation), has a border entry (a pending function or a
    !deleted base slot that an append will reuse), or is a clean base slot.
    mc = self%m - self%ma
    count_new = 0
    count_beyond = 0
    count_entry = 0
    same_set = .true.
    do i = 1, count_active
      idx = indices(i)
      if (idx > self%n .and. idx > self%n0) then
        count_beyond = count_beyond + 1
      else if (self%bpos(idx) == 0) then
        count_new = count_new + 1
        same_set = .false.
      else
        count_entry = count_entry + 1
        if (self%bpos(idx) <= mc) same_set = .false.
      end if
    end do
    if (count_entry /= self%ma .or. count_beyond /= self%nres) same_set = .false.
    if (same_set) then
      call clear_marks()
      info = QR_SUCCESS
      return
    end if

    if (self%m + count_new + count_beyond > self%mmax) then
      call clear_marks()
      need_flush = .true.
      info = QR_SUCCESS
      return
    end if
    if (count_new > 0 .and. self%r_singular) then
      call clear_marks()
      info = QR_ERR_SINGULAR
      return
    end if

    !Commit the previous tail and move the announced entries to the end of
    !the border list, keeping the relative order of both groups.
    committed_changed = .false.
    k = 0
    do p = 1, self%m
      if (self%imark(self%bidx(p)) == 0) then
        k = k + 1
        self%iw(k, 1) = self%bidx(p)
        self%iw(k, 2) = self%bslot(p)
        if (p > mc) committed_changed = .true.
      end if
    end do
    self%ma = self%m - k
    do p = 1, self%m
      if (self%imark(self%bidx(p)) /= 0) then
        k = k + 1
        self%iw(k, 1) = self%bidx(p)
        self%iw(k, 2) = self%bslot(p)
        if (p <= mc) committed_changed = .true.
      end if
    end do
    do p = 1, self%m
      self%bidx(p) = self%iw(p, 1)
      self%bslot(p) = self%iw(p, 2)
      self%bpos(self%bidx(p)) = p
    end do
    self%nres = count_beyond
    self%tail_ok = .false.
    call clear_marks()

    !Clean base slots join D0 with a zero modification and enter the tail.
    if (count_new > 0) then
      do i = 1, count_active
        idx = indices(i)
        if (idx > self%n) cycle
        if (self%bpos(idx) /= 0) cycle
        call complex_join_base_slot(self, idx, .false.)
        self%ma = self%ma + 1
      end do
      call complex_rebuild_small(self)
      committed_changed = .true.
    end if

    if (committed_changed) self%scc_ok = .false.
    call complex_refresh_committed(self)
    if (self%r_singular .or. self%scc_singular) then
      info = QR_ERR_SINGULAR
    else
      info = QR_SUCCESS
    end if

  contains

    subroutine clear_marks()
      integer :: mark
      do mark = 1, count_active
        self%imark(indices(mark)) = 0
      end do
    end subroutine clear_marks

  end subroutine complex_begin_active

  subroutine complex_replace_symmetric(self, idx, h_col, s_col, info)
  !Subroutine complex_replace_symmetric supplies the NEW full conceptual column
  !of the active function idx. In contrast to the previous implementation the
  !arguments are the new columns themselves, not differences. The represented
  !matrix becomes the Hermitian matrix whose column idx equals
  !
  !                 a = h_col - shift*s_col ,
  !
  !and whose row idx is the complex conjugate of a, all other entries being
  !unchanged. The arguments are conceptual columns: h_col(i) is the matrix
  !element <phi_i|H|phi_idx> with the row function as the bra, whether i is
  !smaller or larger than idx. The diagonal elements h_col(idx) and
  !s_col(idx) must be real to within QR_DIAGONAL_TOLERANCE times
  !max(1,max(abs(Re)+abs(Im))) of their column; their imaginary parts are
  !discarded. Nothing in Q or R changes. The routine is collective.
  !
  !The column is split into its part on the clean rows (stored in Bt, as the
  !difference from the base column when idx is a base slot, and zero on the
  !rows of D0) and its entries on the border indices (stored in row and
  !column of the border-border block). One distributed base solve gives
  !Gs(:,k). The row and column of BG, the row of BZ, and the row and column of
  !the Schur complement that belong to the entry are updated with
  !O((2*m+m0)*n0) local work. The tail factors X and LU(Sigma) are recomputed
  !lazily by the next solve.
  !
  !When several active functions change in one trial, the caller supplies the
  !final column of each of them in any order; the entry shared by two active
  !functions is taken from the column supplied last, so the columns must be
  !mutually consistent (Hermitian) at the end of the trial.
  !
  !  Input parameters:
  !    idx   - Canonical index of the function, 1 <= idx <= order(). It must
  !            belong to the active set announced by begin_active.
  !    h_col - New conceptual column idx of H, length order().
  !    s_col - New conceptual column idx of S, length order().
  !
  !  Input/output parameter:
  !    self  - A valid complex QR state. On success both structural-update
  !            counters are incremented.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS;
  !            QR_ERR_INVALID_STATE when self holds no valid factorization;
  !            QR_ERR_INVALID_ARGUMENT when idx is out of range, or when the
  !            diagonal element of a column is not real within the
  !            tolerance above, or when idx is not an active border entry;
  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from
  !            order();
  !            QR_ERR_SINGULAR when R is numerically singular; the caller
  !            must refactorize.
  !
  !No allocation is performed. An error leaves the state unchanged.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(in) :: h_col(:), s_col(:)
    integer, intent(out) :: info
    integer :: i, p

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    if (idx < 1 .or. idx > self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(h_col) /= self%n .or. size(s_col) /= self%n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    !The diagonal elements of H and of S must be real to within roundoff.
    if (.not. (complex_diagonal_is_real(self%n, h_col, idx) .and. &
               complex_diagonal_is_real(self%n, s_col, idx))) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    p = self%bpos(idx)
    if (p <= self%m - self%ma) then
      !Not pending, or committed: the function was not announced.
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    do i = 1, self%n
      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)
    end do
    !The diagonal element is real; a roundoff-sized imaginary part of the
    !supplied values is discarded.
    self%wn(idx, 3) = cmplx(real(h_col(idx), wp) - &
                           self%shift * real(s_col(idx), wp), 0.0_wp, kind=wp)
    call complex_install_column(self, p, self%wn(:, 3), .false.)

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_replace_symmetric

  subroutine complex_append_symmetric(self, h_col, s_col, info)
  !Subroutine complex_append_symmetric adds the function order()+1 with its
  !full conceptual column. The order grows by one and the new function becomes
  !an entry of the active tail. Nothing in Q or R changes. The routine is
  !collective.
  !
  !The new index normally was announced by begin_active as part of the
  !appended range. An unannounced append is also accepted when the border can
  !hold it and simply extends the active tail. Two cases exist:
  !
  !  - the index lies beyond the base (the usual case): a new border entry is
  !    created (an announced reservation is consumed);
  !  - the index is a trailing base slot that truncate deleted earlier: the
  !    border entry of the deleted slot receives the new column. If that
  !    entry was committed it moves to the tail, which makes the next solve
  !    refactorize the committed block.
  !
  !Both cases cost one distributed base solve and O((2*m+m0)*n0) local work.
  !
  !  Input parameters:
  !    h_col - Conceptual column order()+1 of H, length order()+1: element i
  !            is <phi_i|H|phi_new> with the row function as the bra. The
  !            last element is the new diagonal element; it must be real
  !            to within QR_DIAGONAL_TOLERANCE times
  !            max(1,max(abs(Re)+abs(Im))) of the column, and its imaginary
  !            part is discarded.
  !    s_col - Conceptual column order()+1 of S, length order()+1, with the
  !            same convention.
  !
  !  Input/output parameter:
  !    self  - A valid complex QR state. On success order() grows by one and
  !            both structural-update counters are incremented.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS;
  !            QR_ERR_INVALID_STATE when self holds no valid factorization;
  !            QR_ERR_CAPACITY_EXCEEDED when order() equals the capacity;
  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from
  !            order()+1;
  !            QR_ERR_INVALID_ARGUMENT when the new diagonal element of H or
  !            of S is not real within the tolerance above;
  !            QR_ERR_BATCH_FULL when the index lies beyond the base, was not
  !            announced, and the border is full; the caller must
  !            refactorize;
  !            QR_ERR_SINGULAR when R is numerically singular.
  !
  !No allocation is performed. An error leaves the state unchanged.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h_col(:), s_col(:)
    integer, intent(out) :: info
    integer :: i, new_n, p, p2, slot, slot2

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    new_n = self%n + 1
    if (new_n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    if (size(h_col) /= new_n .or. size(s_col) /= new_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    !The diagonal elements of H and of S must be real to within roundoff.
    if (.not. (complex_diagonal_is_real(new_n, h_col, new_n) .and. &
               complex_diagonal_is_real(new_n, s_col, new_n))) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (new_n <= self%n0) then
      !A deleted base slot always owns a border entry.
      if (self%bpos(new_n) == 0) return
    else if (self%nres == 0 .and. self%m >= self%mmax) then
      info = QR_ERR_BATCH_FULL
      return
    end if
    if (self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    if (new_n <= self%n0) then
      !Reuse the entry of the deleted slot; move it to the end of the tail
      !if it is committed.
      p = self%bpos(new_n)
      if (p <= self%m - self%ma) then
        slot = self%bslot(p)
        do p2 = p, self%m - 1
          self%bidx(p2) = self%bidx(p2 + 1)
          self%bslot(p2) = self%bslot(p2 + 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        p = self%m
        self%bidx(p) = new_n
        self%bslot(p) = slot
        self%bpos(new_n) = p
        self%ma = self%ma + 1
        self%scc_ok = .false.
      end if
    else
      !Create the border entry at the end of the tail.
      if (self%nres > 0) self%nres = self%nres - 1
      p = self%m + 1
      slot = self%bslot(p)
      self%bidx(p) = new_n
      self%bpos(new_n) = p
      self%m = p
      self%ma = self%ma + 1
      self%stale(slot) = 0
      do p2 = 1, self%m
        slot2 = self%bslot(p2)
        self%cbb(slot, slot2) = 0.0_wp
        self%cbb(slot2, slot) = 0.0_wp
      end do
    end if
    self%n = new_n

    do i = 1, new_n
      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)
    end do
    !The diagonal element is real; a roundoff-sized imaginary part of the
    !supplied values is discarded.
    self%wn(new_n, 3) = cmplx(real(h_col(new_n), wp) - &
                           self%shift * real(s_col(new_n), wp), 0.0_wp, kind=wp)
    call complex_install_column(self, p, self%wn(:, 3), .false.)

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_append_symmetric

  subroutine complex_truncate(self, target_order, info)
  !Subroutine complex_truncate removes the trailing functions
  !target_order+1..order(). The routine is collective.
  !
  !  - A removed function beyond the base loses its border entry (O(m) work;
  !    if the entry was committed, LU(Scc) is recomputed lazily).
  !  - A removed pending base function keeps its border entry, which is
  !    marked as a deleted slot. Its decoupling column is installed lazily
  !    (one base solve) by the next solve, unless append_symmetric reuses the
  !    index first, in which case no work at all is spent on the deletion.
  !  - A removed clean base function joins D0 as a deleted slot: one base
  !    solve, one column recovery, an O(m*n0) correction of the stored border
  !    columns, and one distributed rebuild of the small matrices after all
  !    removals. Its entry is committed.
  !
  !A base slot stays in D0 until the next fresh factorization. Announced but
  !unsupplied appended indices beyond the base are forgotten when the order
  !actually decreases.
  !
  !Typical use (basis enlargement): truncate(base) drops the previous
  !candidate, whose border entry is the active tail, without touching the
  !committed set or LU(Scc).
  !
  !  Input parameter:
  !    target_order - New order, 1 <= target_order <= order(). A value equal
  !                   to order() is a no-op.
  !
  !  Input/output parameter:
  !    self - A valid complex QR state. On success the structural-update
  !           counters grow by the number of removed functions.
  !
  !  Output parameter:
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self holds no valid factorization;
  !           QR_ERR_INVALID_ARGUMENT when target_order is out of range;
  !           QR_ERR_BATCH_FULL when the border cannot take the clean base
  !           slots that would be deleted; the state is unchanged and the
  !           caller refactorizes at the target order;
  !           QR_ERR_SINGULAR when R is numerically singular and a base solve
  !           would be required (state unchanged).
  !
  !No allocation is performed.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: target_order
    integer, intent(out) :: info
    integer :: count_new, i, mc, old_n, p, p2, slot
    logical :: grew

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%bt)) return
    if (target_order < 1 .or. target_order > self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (target_order == self%n) then
      info = QR_SUCCESS
      return
    end if

    !Count the clean base slots that need a border entry before changing
    !anything.
    count_new = 0
    do i = target_order + 1, self%n
      if (self%bpos(i) == 0) count_new = count_new + 1
    end do
    if (self%m + count_new > self%mmax) then
      info = QR_ERR_BATCH_FULL
      return
    end if
    if (count_new > 0 .and. self%r_singular) then
      info = QR_ERR_SINGULAR
      return
    end if

    old_n = self%n
    grew = .false.
    do i = old_n, target_order + 1, -1
      p = self%bpos(i)
      mc = self%m - self%ma
      if (i > self%n0) then
        !Appended function: drop the border entry and free its slot.
        if (p <= mc) then
          self%scc_ok = .false.
        else
          self%ma = self%ma - 1
          self%tail_ok = .false.
        end if
        slot = self%bslot(p)
        do p2 = p, self%m - 1
          self%bidx(p2) = self%bidx(p2 + 1)
          self%bslot(p2) = self%bslot(p2 + 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        self%bslot(self%m) = slot
        self%bidx(self%m) = 0
        self%bpos(i) = 0
        self%m = self%m - 1
      else if (p > 0) then
        !Pending base function: its entry becomes a deleted slot.
        self%stale(self%bslot(p)) = 1
        if (p <= mc) then
          self%scc_ok = .false.
        else
          self%tail_ok = .false.
        end if
      else
        !Clean base function: it joins D0 as a deleted slot. The new entry
        !is created at the end of the list and moved in front of the tail.
        call complex_join_base_slot(self, i, .true.)
        p = self%m
        slot = self%bslot(p)
        do p2 = p, mc + 2, -1
          self%bidx(p2) = self%bidx(p2 - 1)
          self%bslot(p2) = self%bslot(p2 - 1)
          self%bpos(self%bidx(p2)) = p2
        end do
        self%bidx(mc + 1) = i
        self%bslot(mc + 1) = slot
        self%bpos(i) = mc + 1
        self%scc_ok = .false.
        grew = .true.
      end if
    end do
    self%n = target_order
    self%nres = 0
    self%structural_updates = self%structural_updates + &
                              int(old_n - target_order, int64)
    self%updates_since_fresh = self%updates_since_fresh + &
                               int(old_n - target_order, int64)
    if (grew) call complex_rebuild_small(self)
    info = QR_SUCCESS
  end subroutine complex_truncate

  subroutine complex_delete_symmetric(self, idx, info)
  !Subroutine complex_delete_symmetric removes the LAST function. It is kept
  !for interface compatibility and is equivalent to truncate(order()-1).
  !Interior deletion is not supported by the deferred scheme: the caller
  !compacts the physical matrices and refactorizes.
  !
  !  Input parameter:
  !    idx  - Index of the function to delete; must equal order(), and the
  !           order must be at least two.
  !
  !  Input/output parameter:
  !    self - A valid complex QR state.
  !
  !  Output parameter:
  !    info - QR_ERR_INVALID_ARGUMENT when idx differs from order() or the
  !           order is one; otherwise the status of truncate.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (self%n < 2 .or. idx /= self%n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    call complex_truncate(self, self%n - 1, info)
  end subroutine complex_delete_symmetric

  subroutine complex_solve(self, h, s, v_initial, x, lambda, tol, max_iter, &
                           norm_mode, rel_acc, num_iter, eig_residual, info, &
                           residual_goal, num_polish)
  !Subroutine complex_solve finds one eigenvalue and its eigenvector for the
  !complex generalized Hermitian eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration with the represented matrix H-shift*S (base
  !factors plus border). H and S are the CURRENT physical matrices; the caller
  !keeps them consistent with the columns supplied to replace_symmetric and
  !append_symmetric. The routine is collective.
  !
  !The desired eigenvalue should be closer to shift than every other
  !eigenvalue. The convergence rate is governed principally by the ratio of
  !the distance from shift to the desired eigenvalue and the distance from
  !shift to the next closest eigenvalue. A nondegenerate eigenvalue is
  !assumed.
  !
  !For a current vector v, one iteration computes the solution of
  !
  !             (H-shift*S)*x = S*v ,
  !
  !up to a scalar factor. The first iteration solves this system directly:
  !the right-hand side is a distributed Hermitian product with the lower
  !triangle of S and the solve is the bordered solve described in the module
  !header (one distributed base solve). Every later iteration obtains the
  !same direction in residual-correction form,
  !
  !             mu = (v^H*H*v)/(v^H*S*v),
  !             x  = v - (H-shift*S)^{-1} * (H*v - mu*S*v)
  !                = (mu-shift) * (H-shift*S)^{-1} * S*v ,
  !
  !at the cost of one additional Hermitian product; mu is real (the real
  !parts of the two computed Hermitian forms are used). The reason is
  !numerical. The bordered solve eliminates on the clean block of the base
  !matrix, whose inverse is represented as a difference of large terms when the
  !base matrix is nearly singular; its error is proportional to the size of the
  !solution and can exceed that of a solve with fresh factors by many orders of
  !magnitude when the pending functions carry much of the eigenvector. In the
  !correction form the solve is applied to the small residual of the current
  !iterate, so this error only slows the convergence slightly and does not
  !limit the attainable accuracy; the fixed point is the eigenvector of the
  !physical pair (H,S). The new vector is divided by its largest absolute
  !real or imaginary component. Its change of direction relative to v is
  !measured, independently of the complex phase of the iterates, by
  !
  !             alpha   = (v^H*x)/(v^H*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0, iteration stops when rel_acc<=tol. For tol<=0, iteration
  !continues until rel_acc no longer decreases (it begins to increase, or the
  !iterate has reached a fixed point or a rounding cycle of the correction
  !form, where rel_acc stays zero or constant) and the current value is not
  !larger than abs(tol). The latter rule requests the most
  !accurate attainable direction subject to the floor abs(tol), and
  !necessarily performs at least one additional comparison iteration.
  !
  !After iteration, two distributed products hs = H*x and ss = S*x give the
  !physical Rayleigh quotient and the physical eigenpair residual
  !
  !             lambda       = Re(x^H*hs)/Re(x^H*ss),
  !             eig_residual = ||hs-lambda*ss||_2 /
  !                            (||hs||_2+abs(lambda)*||ss||_2+tiny(1.0_wp)) .
  !
  !The residual is independent of the factorization and of the border data,
  !so it also detects a damaged representation; the caller decides whether to
  !refactorize. Only the lower triangles of the leading n by n blocks of H and
  !S are referenced.
  !
  !Optional additional correction steps ("polish"). When residual_goal is
  !present and eig_residual exceeds it after the iteration has stopped (for
  !example because the stopping rule was satisfied after the first, directly
  !solved iteration, or because max_iter was reached), the routine applies up
  !to QR_MAX_POLISH further steps
  !
  !             x <- x - (H-shift*S)^{-1} * (H*x - lambda*S*x) ,
  !
  !reusing the two products of the residual evaluation. Each step costs one
  !bordered solve and two Hermitian products and is accepted only if it
  !lowers the physical residual; the steps stop when the goal is met, when a
  !step fails to lower the residual, or when it lowers it by less than a
  !factor of two. This gives the caller a cheap alternative to a fresh
  !factorization when the residual check fails. Without residual_goal no such
  !step is made. The stopping rule, rel_acc, num_iter and info are never
  !affected by these steps.
  !
  !  Input parameters:
  !    h         - Complex array with at least n rows and n columns holding the
  !                lower triangle of the current H.
  !    s         - Complex array with at least n rows and n columns holding the
  !                lower triangle of the current S, expected to be positive
  !                definite.
  !    v_initial - A nonzero starting approximation of length n. The array is
  !                not modified.
  !    tol       - Directional convergence tolerance. Its sign selects the
  !                stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^H*S*x = 1;
  !                  1: x^H*x = 1;
  !                other: largest absolute real or imaginary component of
  !                       x equal to one.
  !                The complex phase of x is not normalized.
  !
  !  Input/output parameter:
  !    self      - A valid complex QR state of order n. The represented matrix
  !                and the counters are unchanged; the lazily maintained tail
  !                factors and the workspace are updated.
  !
  !  Output parameters:
  !    x         - The final eigenvector approximation in the normalization
  !                selected by norm_mode. It has length n.
  !    lambda    - The real Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - The directional difference from the final iteration. This
  !                is a convergence estimate, not a rigorously bounded error.
  !    num_iter  - Number of inverse iterations performed.
  !    eig_residual - Relative physical eigenpair residual defined above.
  !    info      - QR_SUCCESS when the stopping criterion is satisfied;
  !                QR_ERR_INVALID_STATE when self has no valid factorization;
  !                QR_ERR_DIMENSION_MISMATCH when H, S, v_initial, or x has an
  !                extent inconsistent with the order;
  !                QR_ERR_INVALID_ARGUMENT when max_iter is not positive;
  !                QR_ERR_ZERO_INITIAL_VECTOR when v_initial is numerically
  !                zero in the working precision;
  !                QR_ERR_NONPOSITIVE_OVERLAP when the final x^H*S*x is not
  !                greater than the precision-scaled zero threshold;
  !                QR_ERR_SINGULAR when a diagonal element of R is zero,
  !                subnormal, or not finite, when a small matrix of the border
  !                (Scc or Sigma) cannot be factorized, or when an
  !                iteration produces a numerically zero or non-finite vector;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. In this
  !                case x, lambda, rel_acc, num_iter and eig_residual describe
  !                the best approximation reached.
  !
  !  Optional parameters:
  !    residual_goal - Input. Requested bound of eig_residual; enables the
  !                    residual-correction steps described above. A value of
  !                    zero requests the best attainable residual.
  !    num_polish    - Output. Number of accepted correction steps (zero when
  !                    residual_goal is absent).
  !
  !No allocation is performed. On an error detected before the first
  !iteration, x and lambda are zero, num_iter is zero, and rel_acc and
  !eig_residual are huge.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    complex(wp), intent(in) :: v_initial(:)
    complex(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter
    real(wp), intent(out) :: eig_residual
    integer, intent(out) :: info
    real(wp), intent(in), optional :: residual_goal
    integer, intent(out), optional :: num_polish
    integer :: i, matrix_n, polish_count
    real(wp) :: goal

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    eig_residual = huge(1.0_wp)
    if (present(num_polish)) num_polish = 0
    info = QR_ERR_INVALID_STATE

    !Verify the state, array dimensions, iteration limit, and starting-vector
    !norm before any numerical work.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v_initial) /= matrix_n .or. size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (max_iter <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    do i = 1, matrix_n
      self%wn(i, 3) = v_initial(i)
    end do
    if (complex_norm_squared(matrix_n, self%wn(:, 3)) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    !Bring the lazily maintained small factors up to date and refuse to
    !divide by unusable pivots.
    call complex_refresh_small(self)
    if (.not. complex_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !A negative goal disables the correction steps.
    goal = -1.0_wp
    if (present(residual_goal)) goal = max(residual_goal, 0.0_wp)
    call complex_inverse_iteration(self, matrix_n, h, size(h, 1), s, &
                                   size(s, 1), self%wn(:, 3), self%wn(:, 4), &
                                   self%wn(:, 5), self%wn(:, 6), lambda, tol, &
                                   max_iter, norm_mode, goal, rel_acc, &
                                   num_iter, eig_residual, polish_count, info)
    if (present(num_polish)) num_polish = polish_count
    do i = 1, matrix_n
      x(i) = self%wn(i, 5)
    end do
  end subroutine complex_solve

  subroutine complex_inverse_iteration(self, n, h, ldh, s, lds, v, w, x, y, &
                                       lambda, tol, max_iter, norm_mode, &
                                       residual_goal, rel_acc, num_iter, &
                                       eig_residual, num_polish, info)
  !Subroutine complex_inverse_iteration is the numerical body of complex_solve.
  !The work vectors are passed as explicit-shape arguments so that the loops
  !operate on plain arrays.
  !
  !  Input parameters:
  !    n, h, ldh, s, lds - Order and the two physical matrices with their
  !                        leading dimensions (lower triangles referenced).
  !    tol, max_iter, norm_mode - As for complex_solve.
  !    residual_goal - Goal of the residual-correction steps; a negative value
  !                    disables them.
  !
  !  Input/output parameters:
  !    self - Valid state whose small factors are up to date.
  !    v    - On entry the starting vector; overwritten (finally by H*x).
  !    w    - Work vector (S*v during iteration, finally S*x).
  !    y    - Work vector (H*v and the residual right-hand side).
  !
  !  Output parameters:
  !    x, lambda, rel_acc, num_iter, eig_residual, info - As for complex_solve;
  !    x is zero when an error prevents a usable iterate.
  !    num_polish - Number of accepted correction steps.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: n, ldh, lds
    complex(wp), intent(in) :: h(ldh, *), s(lds, *)
    complex(wp), intent(inout) :: v(n), w(n), y(n)
    complex(wp), intent(out) :: x(n)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(in) :: residual_goal
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter
    real(wp), intent(out) :: eig_residual
    integer, intent(out) :: num_polish, info
    integer :: i
    complex(wp) :: coefficient
    real(wp) :: current_norm_squared, eigenvector_norm_squared, quotient
    real(wp) :: max_component, new_lambda, new_overlap, new_residual
    real(wp) :: norm_of_diff, norm_of_diff_previous, overlap_norm_squared
    real(wp) :: previous_residual
    logical :: corrected, not_converged

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    num_polish = 0
    eig_residual = huge(1.0_wp)
    norm_of_diff = huge(1.0_wp)
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform inverse iterations until the direction criterion is satisfied or
    !the caller-supplied iteration limit is exhausted.
    do while (not_converged .and. num_iter < max_iter)
      !Form w=S*v from the lower triangle of S.
      call complex_sym_product(self%rank, self%nprocs, self%comm, n, s, lds, &
                               v, w)
      corrected = .false.
      if (num_iter > 0) then
        !Residual-correction form of the same iterate. v is the previous
        !iterate, an approximate eigenvector; with its Rayleigh quotient mu,
        !    v - A^{-1}*(H*v - mu*S*v) = (mu-shift)*A^{-1}*S*v .
        !The solve is applied to a small right-hand side, so the errors of
        !the bordered solve are proportional to the correction and not to v.
        overlap_norm_squared = real(complex_dotc(n, v, w), wp)
        if (overlap_norm_squared > tiny(1.0_wp)) then
          call complex_sym_product(self%rank, self%nprocs, self%comm, n, h, &
                                   ldh, v, y)
          quotient = real(complex_dotc(n, v, y), wp) / overlap_norm_squared
          do i = 1, n
            y(i) = y(i) - quotient * w(i)
          end do
          call complex_bordered_solve(self, y, x)
          do i = 1, n
            x(i) = v(i) - x(i)
          end do
          corrected = .true.
        end if
      end if
      !First iteration (or unusable quotient): solve (H-shift*S)*x=S*v.
      if (.not. corrected) call complex_bordered_solve(self, w, x)

      !Scale the solution with its largest real or imaginary component. This
      !bounds both parts of every element without changing the complex
      !eigendirection and limits growth in repeated solves.
      max_component = complex_max_abs_real_or_imag(n, x)
      if (.not. (max_component > tiny(1.0_wp) .and. &
                 max_component <= huge(1.0_wp))) then
        x = 0.0_wp
        info = QR_ERR_SINGULAR
        return
      end if
      do i = 1, n
        x(i) = x(i) / max_component
      end do

      !Compute the least-direction-change coefficient (v^H*x)/(v^H*v) and the
      !relative norm of the component of x not parallel to the previous
      !iterate v. Conjugating the previous iterate makes the comparison
      !invariant under the arbitrary complex phase of an eigenvector.
      current_norm_squared = complex_norm_squared(n, v)
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = complex_dotc(n, v, x) / current_norm_squared
      eigenvector_norm_squared = complex_norm_squared(n, x)
      norm_of_diff = complex_direction_difference(n, x, coefficient, v, &
                                                  eigenvector_norm_squared)
      !A partially non-finite iterate is detected through its direction error.
      if (.not. (norm_of_diff <= huge(1.0_wp))) then
        x = 0.0_wp
        info = QR_ERR_SINGULAR
        return
      end if

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        !For a non-positive tolerance, accept only after the direction error
        !has passed its minimum and remains within abs(tol). In the
        !residual-correction form the iterate can reach a fixed point or a
        !rounding cycle, where the direction error is zero or constant and
        !never increases; a value that no longer decreases is therefore the
        !minimum as well.
        if (norm_of_diff >= norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        do i = 1, n
          v(i) = x(i)
        end do
      end if
    end do

    rel_acc = norm_of_diff
    !Failure to meet the stopping rule is nonfatal: the final iterate is still
    !used to calculate and normalize an eigenpair approximation.
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    !Physical products for the Rayleigh quotient, the residual, and, for
    !norm_mode=0, the final normalization. A non-positive x^H*S*x violates the
    !positive-definite S precondition or indicates unusable numerical data.
    call complex_physical_quotient(self%rank, self%nprocs, self%comm, n, h, &
                                   ldh, s, lds, x, v, w, lambda, &
                                   overlap_norm_squared, eig_residual)
    if (.not. (overlap_norm_squared > tiny(1.0_wp))) then
      lambda = 0.0_wp
      eig_residual = huge(1.0_wp)
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if

    !Optional residual-correction steps: y = x - A^{-1}*(H*x - lambda*S*x).
    !On entry to every step v = H*x and w = S*x for the current x.
    if (residual_goal >= 0.0_wp) then
      do while (eig_residual > residual_goal .and. num_polish < QR_MAX_POLISH)
        do i = 1, n
          v(i) = v(i) - lambda * w(i)
        end do
        call complex_bordered_solve(self, v, y)
        do i = 1, n
          y(i) = x(i) - y(i)
        end do
        max_component = complex_max_abs_real_or_imag(n, y)
        if (.not. (max_component > tiny(1.0_wp) .and. &
                   max_component <= huge(1.0_wp))) exit
        do i = 1, n
          y(i) = y(i) / max_component
        end do
        call complex_physical_quotient(self%rank, self%nprocs, self%comm, n, &
                                       h, ldh, s, lds, y, v, w, new_lambda, &
                                       new_overlap, new_residual)
        !Reject a step that does not lower the physical residual (the negated
        !comparison also rejects NaN) and keep the current x.
        if (.not. (new_overlap > tiny(1.0_wp))) exit
        if (.not. (new_residual < eig_residual)) exit
        do i = 1, n
          x(i) = y(i)
        end do
        previous_residual = eig_residual
        lambda = new_lambda
        overlap_norm_squared = new_overlap
        eig_residual = new_residual
        num_polish = num_polish + 1
        !Stop when the steps no longer make substantial progress.
        if (eig_residual > 0.5_wp * previous_residual) exit
      end do
    end if

    !Apply the normalization requested by the caller. No scaling is required
    !for other norm_mode values because the largest real or imaginary
    !component of each iterate already has unit magnitude.
    select case (norm_mode)
    case (0) !Normalize so that x^T*S*x=1.
      quotient = sqrt(overlap_norm_squared)
      do i = 1, n
        x(i) = x(i) / quotient
      end do
    case (1) !Normalize so that x^T*x=1.
      eigenvector_norm_squared = complex_norm_squared(n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      quotient = sqrt(eigenvector_norm_squared)
      do i = 1, n
        x(i) = x(i) / quotient
      end do
    end select
  end subroutine complex_inverse_iteration

  subroutine complex_physical_quotient(rank, nprocs, comm, n, h, ldh, s, lds, &
                                       x, hx, sx, lambda, overlap, residual)
  !Subroutine complex_physical_quotient evaluates, with two distributed
  !Hermitian products, the physical Rayleigh quotient and the relative
  !eigenpair residual of a vector:
  !
  !        hx = H*x,  sx = S*x,  overlap = Re(x^H*sx),
  !        lambda = Re(x^H*hx)/overlap,
  !        residual = ||hx-lambda*sx||_2 /
  !                   (||hx||_2+abs(lambda)*||sx||_2+tiny(1.0_wp)) .
  !
  !  Input parameters:
  !    rank, nprocs, comm - Own rank, number of ranks, and communicator.
  !    n, h, ldh, s, lds - Order and the physical matrices (lower triangles).
  !    x    - Replicated vector of length n.
  !
  !  Output parameters:
  !    hx, sx   - The two products.
  !    lambda   - Rayleigh quotient; zero when overlap is not positive.
  !    overlap  - Real part of x^H*S*x.
  !    residual - Relative residual; huge when overlap is not positive.
    integer, intent(in) :: rank, nprocs, comm, n, ldh, lds
    complex(wp), intent(in) :: h(ldh, *), s(lds, *), x(n)
    complex(wp), intent(out) :: hx(n), sx(n)
    real(wp), intent(out) :: lambda, overlap, residual
    integer :: i
    complex(wp) :: difference
    real(wp) :: residual_norm_squared

    call complex_sym_product(rank, nprocs, comm, n, h, ldh, x, hx)
    call complex_sym_product(rank, nprocs, comm, n, s, lds, x, sx)
    overlap = real(complex_dotc(n, x, sx), wp)
    if (.not. (overlap > tiny(1.0_wp))) then
      lambda = 0.0_wp
      residual = huge(1.0_wp)
      return
    end if
    lambda = real(complex_dotc(n, x, hx), wp) / overlap
    residual_norm_squared = 0.0_wp
    do i = 1, n
      difference = hx(i) - lambda * sx(i)
      residual_norm_squared = residual_norm_squared + &
                              real(difference, wp)**2 + aimag(difference)**2
    end do
    residual = sqrt(residual_norm_squared) / &
               (sqrt(complex_norm_squared(n, hx)) + &
                abs(lambda) * sqrt(complex_norm_squared(n, sx)) + tiny(1.0_wp))
  end subroutine complex_physical_quotient

  subroutine complex_solve_shifted(self, b, x, info)
  !Subroutine complex_solve_shifted solves one linear system
  !
  !                     (H-shift*S)*x = b
  !
  !with the represented matrix of the current order (base factors plus
  !border). It exposes the bordered solve used by inverse iteration for tests
  !and diagnostics. The routine is collective; b must be replicated.
  !
  !  Input parameter:
  !    b    - Right-hand side of length order().
  !
  !  Input/output parameter:
  !    self - A valid complex QR state. The represented matrix is unchanged;
  !           the lazily maintained tail factors and the workspace are updated.
  !
  !  Output parameters:
  !    x    - Solution of length order(); zero on error.
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self holds no valid factorization;
  !           QR_ERR_DIMENSION_MISMATCH when b or x does not have length
  !           order();
  !           QR_ERR_SINGULAR when R or a small matrix of the border cannot
  !           be used.
  !
  !No allocation is performed.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: b(:)
    complex(wp), intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: i

    x = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    if (size(b) /= self%n .or. size(x) /= self%n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    call complex_refresh_small(self)
    if (.not. complex_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if
    do i = 1, self%n
      self%wn(i, 3) = b(i)
    end do
    call complex_bordered_solve(self, self%wn(:, 3), self%wn(:, 4))
    do i = 1, self%n
      x(i) = self%wn(i, 4)
    end do
    info = QR_SUCCESS
  end subroutine complex_solve_shifted

  subroutine complex_symmetric_product(self, a, x, y, info)
  !Subroutine complex_symmetric_product computes y = A*x for a replicated
  !complex Hermitian matrix of which only the lower triangle is referenced (the
  !strictly upper triangle is taken to be the conjugate transpose of the
  !strictly lower one, and the imaginary parts of the diagonal elements are
  !ignored), with the
  !distributed kernel used by inverse iteration (cyclic column ownership, one
  !reduction; no MPI call with one process). It is offered to callers that
  !need products with the physical matrices, for example for their own
  !residual checks, and to the timing driver. The routine is collective. The
  !binding keeps the name symmetric_product of the real state.
  !
  !The order of the product is size(x); it is independent of order() and a
  !factorization need not exist, but the state must be initialized and the
  !order must not exceed its capacity.
  !
  !  Input parameters:
  !    a    - Complex array with at least size(x) rows and columns holding the
  !           lower triangle of the matrix.
  !    x    - Replicated vector.
  !
  !  Input/output parameter:
  !    self - An initialized complex QR state; two work vectors are
  !           overwritten.
  !
  !  Output parameters:
  !    y    - Replicated product, same length as x; zero on error.
  !    info - QR_SUCCESS;
  !           QR_ERR_INVALID_STATE when self is not initialized;
  !           QR_ERR_DIMENSION_MISMATCH when x is empty, y has a different
  !           length, or a is too small;
  !           QR_ERR_CAPACITY_EXCEEDED when size(x) exceeds the capacity.
  !
  !No allocation is performed.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: a(:,:)
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(out) :: y(:)
    integer, intent(out) :: info
    integer :: i, n

    y = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%wn)) return
    n = size(x)
    if (n <= 0 .or. size(y) /= n .or. size(a, 1) < n .or. size(a, 2) < n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (n > self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    do i = 1, n
      self%wn(i, 3) = x(i)
    end do
    call complex_sym_product(self%rank, self%nprocs, self%comm, n, a, &
                             size(a, 1), self%wn(:, 3), self%wn(:, 4))
    do i = 1, n
      y(i) = self%wn(i, 4)
    end do
    info = QR_SUCCESS
  end subroutine complex_symmetric_product

  subroutine complex_shifted_quotient(self, h, s, x, lambda, info)
  !Subroutine complex_shifted_quotient evaluates the Rayleigh quotient of a
  !vector for the physical pair (H,S) in the shifted form
  !
  !        lambda = shift + (x^H*(H-shift*S)*x)/(x^H*S*x) ,
  !
  !where every element of H-shift*S is formed before it is multiplied by the
  !components of x. For an eigenvector whose eigenvalue is close to the shift
  !this is far more accurate than the plain quotient (x^H*H*x)/(x^H*S*x) that
  !solve returns: the rounding errors of the sums act on the small difference
  !lambda-shift only, and H and shift*S cancel element by element before the
  !products with the large, sign-alternating coefficients of a nearly
  !linearly dependent basis are accumulated. Method I obtains its energy in
  !this way from the stored shifted matrix. For a 2850-function real basis in
  !extended precision the plain quotient of the real state was measured to be
  !wrong by up to 15 units of the last place and the shifted quotient by at
  !most one half, which matters to an optimizer that resolves energy
  !differences of that size. Both quadratic forms are real: the strictly
  !upper triangles are taken to be the conjugate transposes of the strictly
  !lower ones and the imaginary parts of the diagonal elements are ignored,
  !as in symmetric_product.
  !
  !The caller applies this routine to the vector returned by solve when the
  !eigenvalue is needed to the last digit. The routine is collective: the
  !columns of the lower triangles are assigned to the ranks cyclically, each
  !rank accumulates the two quadratic forms over its columns, and one
  !reduction of two values completes them. With one process no MPI call is
  !made.
  !
  !  Input parameters:
  !    h - Complex array with at least n rows and n columns holding the lower
  !        triangle of the current H.
  !    s - Complex array with at least n rows and n columns holding the lower
  !        triangle of the current S.
  !    x - Replicated vector of length n = order(), normally an approximate
  !        eigenvector. It need not be normalized.
  !
  !  Input parameter (state):
  !    self - A valid complex QR state; only its order, shift, and distribution
  !           data are used.
  !
  !  Output parameters:
  !    lambda - The shifted Rayleigh quotient; zero on error.
  !    info   - QR_SUCCESS;
  !             QR_ERR_INVALID_STATE when self holds no valid factorization;
  !             QR_ERR_DIMENSION_MISMATCH when H, S, or x has an extent
  !             inconsistent with the order;
  !             QR_ERR_NONPOSITIVE_OVERLAP when x^H*S*x is not positive.
  !
  !No allocation is performed and the state is not modified.
    class(qr_complex_state), intent(in) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    complex(wp), intent(in) :: x(:)
    real(wp), intent(out) :: lambda
    integer, intent(out) :: info
    integer :: j, matrix_n, mpi_err
    real(wp) :: column_form, column_overlap, forms(2)

    lambda = 0.0_wp
    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if

    !forms(1) accumulates x^H*(H-shift*S)*x and forms(2) x^H*S*x over the
    !columns owned by this rank. A column contributes its real diagonal term
    !and twice the real part of its part below the diagonal.
    forms = 0.0_wp
    do j = self%rank + 1, matrix_n, self%nprocs
      call complex_shifted_column(matrix_n - j, h(j:matrix_n, j), &
                                  s(j:matrix_n, j), self%shift, &
                                  x(j:matrix_n), column_form, column_overlap)
      forms(1) = forms(1) + column_form
      forms(2) = forms(2) + column_overlap
    end do
    if (self%nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, forms, 2, MPI_WP, MPI_SUM, self%comm, &
                         mpi_err)
    end if
    !The negated comparison also rejects NaN.
    if (.not. (forms(2) > tiny(1.0_wp))) then
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if
    lambda = self%shift + forms(1) / forms(2)
    info = QR_SUCCESS
  end subroutine complex_shifted_quotient

  pure subroutine complex_shifted_column(m, h_col, s_col, shift, x, form, &
                                         overlap)
  !Subroutine complex_shifted_column returns the contributions of one column of
  !the lower triangles to the real quadratic forms x^H*(H-shift*S)*x and
  !x^H*S*x. With the stored elements a(i,j), i>=j, of a Hermitian matrix,
  !column j contributes
  !
  !        Re(a(j,j))*abs(x(j))**2 + 2*Re(x(j)*sum_{i>j} a(i,j)*conj(x(i))) ,
  !
  !the second term being the element (i,j) and its conjugate (j,i) taken
  !together.
  !
  !  Input parameters:
  !    m     - Number of elements below the diagonal.
  !    h_col - Column of H from its diagonal element downwards (m+1 elements).
  !    s_col - The same part of the column of S.
  !    shift - The real shift.
  !    x     - The matching part of the vector (m+1 elements).
  !
  !  Output parameters:
  !    form    - Contribution to x^H*(H-shift*S)*x.
  !    overlap - Contribution to x^H*S*x.
    integer, intent(in) :: m
    complex(wp), intent(in) :: h_col(0:m), s_col(0:m), x(0:m)
    real(wp), intent(in) :: shift
    real(wp), intent(out) :: form, overlap
    integer :: i
    complex(wp) :: form_below, overlap_below
    real(wp) :: x_squared

    form_below = 0.0_wp
    overlap_below = 0.0_wp
    do i = 1, m
      form_below = form_below + (h_col(i) - shift * s_col(i)) * conjg(x(i))
      overlap_below = overlap_below + s_col(i) * conjg(x(i))
    end do
    !The diagonal elements are real; the imaginary parts of the supplied
    !values are ignored, as in the products with H and S.
    x_squared = real(x(0), wp)**2 + aimag(x(0))**2
    form = (real(h_col(0), wp) - shift * real(s_col(0), wp)) * x_squared + &
           2.0_wp * real(x(0) * form_below, wp)
    overlap = real(s_col(0), wp) * x_squared + &
              2.0_wp * real(x(0) * overlap_below, wp)
  end subroutine complex_shifted_column

  subroutine complex_factorization_residual(self, h, s, v, absolute_residual, &
                                            relative_residual, info, &
                                            forward_error)
  !Subroutine complex_factorization_residual measures how well the represented
  !matrix (base factors plus border) reproduces the action of the caller-owned
  !physical shifted matrix A = H-shift*S on a probe vector. The comparison of
  !Q*(R*v) with A*v used by the previous implementation is not possible once
  !updates are deferred, because Q*R is only the base matrix. The check is
  !therefore made through the bordered solve:
  !
  !       b = A*v                      (physical action, two products)
  !       x = solve(represented A, b)  (bordered solve)
  !       difference = A*x - b         (physical action again)
  !
  !       absolute_residual = ||difference||_2,
  !       relative_residual = ||difference||_2 /
  !                           (||A*x||_2+||b||_2+tiny(1.0_wp)),
  !       forward_error     = ||x-v||_2/||v||_2      (optional).
  !
  !The relative residual has the same scaling as the action residual of the
  !previous implementation and is of the order of the working precision times
  !a modest factor when the representation is healthy, also for the nearly
  !singular shifted matrices of inverse iteration, because the large error of
  !x along the near-null vector is annihilated by A. The optional forward
  !error is the plain relative difference between x and v; it is proportional
  !to the condition number of A and is therefore large (harmlessly) when the
  !shift is close to an eigenvalue. Refresh policy belongs to the caller, who
  !should judge relative_residual. The routine is collective and costs one
  !bordered solve and four distributed Hermitian products.
  !
  !  Input parameters:
  !    h - Complex array with at least n rows and n columns holding the lower
  !        triangle of the current H.
  !    s - Complex array with at least n rows and n columns holding the lower
  !        triangle of the current S.
  !    v - Nonzero complex probe vector with exactly n elements. A fixed dense
  !        vector is appropriate; the eigenvector itself is not, since A
  !        nearly annihilates it.
  !
  !  Input/output parameter:
  !    self - Valid complex QR state. The represented matrix and the counters
  !           are unchanged; the workspace is overwritten.
  !
  !  Output parameters:
  !    absolute_residual - Euclidean norm of difference; zero on error.
  !    relative_residual - Scaled action error defined above; zero on error.
  !    info              - QR_SUCCESS;
  !                        QR_ERR_INVALID_STATE when no valid factorization is
  !                        available;
  !                        QR_ERR_DIMENSION_MISMATCH when H or S cannot contain
  !                        the active block or V does not have length n;
  !                        QR_ERR_ZERO_INITIAL_VECTOR when V is numerically
  !                        zero;
  !                        QR_ERR_SINGULAR when R or a small matrix of the
  !                        border cannot be used.
  !    forward_error     - Optional relative forward error defined above; zero
  !                        on error.
  !
  !The routine performs no allocation and no input/output.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    complex(wp), intent(in) :: v(:)
    real(wp), intent(out) :: absolute_residual, relative_residual
    integer, intent(out) :: info
    real(wp), intent(out), optional :: forward_error
    integer :: i, matrix_n
    complex(wp) :: difference
    real(wp) :: action_norm, difference_norm_squared
    real(wp) :: probe_norm_squared, rhs_norm

    absolute_residual = 0.0_wp
    relative_residual = 0.0_wp
    if (present(forward_error)) forward_error = 0.0_wp
    info = QR_ERR_INVALID_STATE

    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%qt) .or. .not. allocated(self%wn)) return
    matrix_n = self%n
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    do i = 1, matrix_n
      self%wn(i, 5) = v(i)
    end do
    probe_norm_squared = complex_norm_squared(matrix_n, self%wn(:, 5))
    if (probe_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if
    call complex_refresh_small(self)
    if (.not. complex_can_solve(self)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !b = (H-shift*S)*v in work vector 3; x in work vector 4.
    call complex_sym_product(self%rank, self%nprocs, self%comm, matrix_n, h, &
                             size(h, 1), self%wn(:, 5), self%wn(:, 3))
    call complex_sym_product(self%rank, self%nprocs, self%comm, matrix_n, s, &
                             size(s, 1), self%wn(:, 5), self%wn(:, 1))
    do i = 1, matrix_n
      self%wn(i, 3) = self%wn(i, 3) - self%shift * self%wn(i, 1)
    end do
    call complex_bordered_solve(self, self%wn(:, 3), self%wn(:, 4))

    if (present(forward_error)) then
      difference_norm_squared = 0.0_wp
      do i = 1, matrix_n
        difference = self%wn(i, 4) - self%wn(i, 5)
        difference_norm_squared = difference_norm_squared + &
                                  real(difference, wp)**2 + &
                                  aimag(difference)**2
      end do
      forward_error = sqrt(difference_norm_squared / probe_norm_squared)
    end if

    !A*x in work vector 1, using work vector 2 for S*x.
    call complex_sym_product(self%rank, self%nprocs, self%comm, matrix_n, h, &
                             size(h, 1), self%wn(:, 4), self%wn(:, 1))
    call complex_sym_product(self%rank, self%nprocs, self%comm, matrix_n, s, &
                             size(s, 1), self%wn(:, 4), self%wn(:, 2))
    do i = 1, matrix_n
      self%wn(i, 1) = self%wn(i, 1) - self%shift * self%wn(i, 2)
    end do
    rhs_norm = sqrt(complex_norm_squared(matrix_n, self%wn(:, 3)))
    action_norm = sqrt(complex_norm_squared(matrix_n, self%wn(:, 1)))
    difference_norm_squared = 0.0_wp
    do i = 1, matrix_n
      difference = self%wn(i, 1) - self%wn(i, 3)
      difference_norm_squared = difference_norm_squared + &
                                real(difference, wp)**2 + &
                                aimag(difference)**2
    end do
    absolute_residual = sqrt(difference_norm_squared)
    relative_residual = absolute_residual / &
                        (action_norm + rhs_norm + tiny(1.0_wp))
    info = QR_SUCCESS
  end subroutine complex_factorization_residual

  pure function complex_is_valid(self) result(factors_are_valid)
  !Function complex_is_valid reports whether a complex QR state currently
  !contains a completed base factorization. It returns true after a successful
  !factorize_fresh and false for a default, cleared, or merely initialized
  !state, and after a failed factorization. It does not test numerical
  !singularity or the health of the border.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    factors_are_valid - Stored validity flag.
    class(qr_complex_state), intent(in) :: self
    logical :: factors_are_valid

    factors_are_valid = self%valid
  end function complex_is_valid

  pure function complex_order(self) result(active_order)
  !Function complex_order returns the current order n of the represented
  !matrix. It grows with append_symmetric, shrinks with truncate and
  !delete_symmetric, and is zero when no factorization is held. It may be
  !smaller or larger than base_order().
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    active_order - Current order n.
    class(qr_complex_state), intent(in) :: self
    integer :: active_order

    active_order = self%n
  end function complex_order

  pure function complex_base_order(self) result(order_of_base)
  !Function complex_base_order returns the order n0 of the base factorization,
  !that is, order() at the time of the last successful factorize_fresh.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    order_of_base - Base order n0; zero when no factorization is held.
    class(qr_complex_state), intent(in) :: self
    integer :: order_of_base

    order_of_base = self%n0
  end function complex_base_order

  pure function complex_get_capacity(self) result(allocated_capacity)
  !Function complex_get_capacity returns the maximum matrix order for which the
  !state has allocated storage. A cleared or failed-initialization state
  !returns zero.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    allocated_capacity - Maximum supported order.
    class(qr_complex_state), intent(in) :: self
    integer :: allocated_capacity

    allocated_capacity = self%capacity
  end function complex_get_capacity

  pure function complex_get_shift(self) result(represented_shift)
  !Function complex_get_shift returns the shift of the represented matrix
  !H-shift*S. The value is updated only by a successful factorize_fresh and
  !is zero in an initialized state that has not been factorized.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    represented_shift - Stored shift.
    class(qr_complex_state), intent(in) :: self
    real(wp) :: represented_shift

    represented_shift = self%shift
  end function complex_get_shift

  pure function complex_pending_count(self) result(count)
  !Function complex_pending_count returns the number m of border entries: the
  !functions whose current row and column are held outside the base
  !factorization, plus the trailing base slots deleted by truncate (each of
  !them occupies an entry until the next fresh factorization). Announced
  !appended indices that have not been supplied yet are not counted. The
  !border is full when the value reaches max_pending().
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    count - Number of border entries, 0 <= count <= max_pending().
    class(qr_complex_state), intent(in) :: self
    integer :: count

    count = self%m
  end function complex_pending_count

  pure function complex_active_count(self) result(count)
  !Function complex_active_count returns the size of the current active set:
  !the entries of the active tail (including deleted base slots that are
  !waiting to be reused by an append) plus the announced appended indices
  !beyond the base that have not been supplied yet.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    count - Size of the active set.
    class(qr_complex_state), intent(in) :: self
    integer :: count

    count = self%ma + self%nres
  end function complex_active_count

  pure function complex_excluded_count(self) result(count)
  !Function complex_excluded_count returns the number m0 of members of D0, the
  !base slots whose row and column are modified: pending functions with an
  !index not greater than base_order() and trailing base slots deleted by
  !truncate. Every member owns a border entry, so the value never exceeds
  !pending_count(); only a fresh factorization reduces it.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    count - Number of modified base slots.
    class(qr_complex_state), intent(in) :: self
    integer :: count

    count = self%m0
  end function complex_excluded_count

  pure function complex_max_pending(self) result(count)
  !Function complex_max_pending returns the largest number of border entries
  !the state can hold, min(max_pending,capacity) of the initialize call.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    count - Border capacity; zero for an uninitialized state.
    class(qr_complex_state), intent(in) :: self
    integer :: count

    count = self%mmax
  end function complex_max_pending

  pure function complex_get_update_count(self) result(update_count)
  !Function complex_get_update_count returns the lifetime number of successful
  !structural updates (replace_symmetric and append_symmetric calls, and
  !functions removed by truncate or delete_symmetric). Fresh factorization
  !does not reset the value; clear and initialize do.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    update_count - Lifetime count as a 64-bit integer.
    class(qr_complex_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%structural_updates
  end function complex_get_update_count

  pure function complex_get_updates_since_fresh(self) result(update_count)
  !Function complex_get_updates_since_fresh returns the number of successful
  !structural updates since the most recent successful factorize_fresh.
  !
  !  Input parameter:
  !    self - The complex QR state to query.
  !
  !  Result:
  !    update_count - Count since the last fresh factorization.
    class(qr_complex_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%updates_since_fresh
  end function complex_get_updates_since_fresh

  !===========================================================================
  ! Internal procedures: border bookkeeping and small matrices
  !===========================================================================

  subroutine complex_reset_border(self)
  !Subroutine complex_reset_border empties the border list and D0 and restores
  !the free-slot list. It is called when the base factorization is replaced.
  !
  !  Input/output parameter:
  !    self - State with allocated bookkeeping arrays.
    class(qr_complex_state), intent(inout) :: self
    integer :: i

    self%m = 0
    self%ma = 0
    self%nres = 0
    self%m0 = 0
    self%bidx = 0
    self%bpos = 0
    self%d0idx = 0
    self%d0pos = 0
    self%stale = 0
    self%cu = 0
    self%imark = 0
    do i = 1, self%mmax
      self%bslot(i) = i
    end do
    self%scc_ok = .true.
    self%tail_ok = .true.
    self%scc_singular = .false.
    self%sig_singular = .false.
  end subroutine complex_reset_border

  pure function complex_can_solve(self) result(usable)
  !Function complex_can_solve reports whether the bordered solve may divide by
  !the stored pivots: R, Scc and Sigma are all usable and up to date.
  !
  !  Input parameter:
  !    self - State after complex_refresh_small.
  !
  !  Result:
  !    usable - True when a bordered solve can be performed.
    class(qr_complex_state), intent(in) :: self
    logical :: usable

    usable = .not. (self%r_singular .or. self%scc_singular .or. &
                    self%sig_singular) .and. self%scc_ok .and. self%tail_ok
  end function complex_can_solve

  subroutine complex_join_base_slot(self, idx, deleted)
  !Subroutine complex_join_base_slot makes the clean base slot idx a member of
  !D0 and creates its border entry at the end of the border list. It computes
  !the column Zb(:,q) = A_base^{-1}*e_idx by one distributed base solve and
  !recovers the base column Ab(:,q) = Q*(R*e_idx) from the factors. The new
  !entry starts with a zero modification (the function is unchanged), so its
  !border column and its base solve are zero. The entry of every older border
  !column on row idx moves into the border-border block, and the stored base
  !solve of that column is corrected by a multiple of Zb(:,q).
  !
  !The caller is responsible for the position of the entry in the border list
  !(tail or committed), for the counter of the tail, and for calling
  !complex_rebuild_small afterwards.
  !
  !  Input parameters:
  !    idx     - Base slot, 1 <= idx <= n0, not yet a member of D0. The caller
  !              has verified that the border has room and that R is usable.
  !    deleted - True when the slot is being deleted by truncate. The entry
  !              is then marked stale; its decoupling column is installed
  !              later by complex_install_column.
  !
  !  Input/output parameter:
  !    self - State; m and m0 grow by one. Work vector 1 is overwritten.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    logical, intent(in) :: deleted
    integer :: i, idx2, p, p2, q, slot, slot2
    complex(wp) :: current, entry

    q = self%m0 + 1
    self%m0 = q
    self%d0idx(q) = idx
    self%d0pos(idx) = q
    do i = 1, self%n0
      self%wn(i, 1) = 0.0_wp
    end do
    self%wn(idx, 1) = 1.0_wp
    call complex_base_solve(self, self%wn(:, 1), self%zb(:, q))
    call complex_base_column(self, idx, self%ab(:, q), self%wn(:, 1))

    p = self%m + 1
    slot = self%bslot(p)
    self%bidx(p) = idx
    self%bpos(idx) = p
    self%m = p
    do i = 1, self%n0
      self%bt(i, slot) = 0.0_wp
      self%gs(i, slot) = 0.0_wp
    end do
    if (deleted) then
      self%stale(slot) = 1
    else
      self%stale(slot) = 0
    end if
    do p2 = 1, p
      slot2 = self%bslot(p2)
      self%cbb(slot, slot2) = 0.0_wp
      self%cbb(slot2, slot) = 0.0_wp
    end do
    self%cbb(slot, slot) = self%ab(idx, q)

    do p2 = 1, p - 1
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      entry = self%bt(idx, slot2)
      !Current value of the matrix element (idx,idx2): the stored column of a
      !member of D0 is a difference from its base column. The element
      !(idx2,idx) of the Hermitian border-border block is its conjugate.
      current = entry
      if (idx2 <= self%n0) current = current + self%ab(idx, self%d0pos(idx2))
      if (deleted .or. self%stale(slot2) /= 0) current = 0.0_wp
      self%cbb(slot, slot2) = current
      self%cbb(slot2, slot) = conjg(current)
      if (entry /= 0.0_wp) then
        self%bt(idx, slot2) = 0.0_wp
        call complex_axpy(self%n0, -entry, self%zb(:, q), self%gs(:, slot2))
      end if
    end do
  end subroutine complex_join_base_slot

  subroutine complex_install_column(self, p, a, dummy)
  !Subroutine complex_install_column installs the current shifted column of the
  !border entry at position p. Its part on the clean rows goes to Bt (as the
  !difference from the base column when the entry is a member of D0), its
  !values on the border indices go to the column, and their complex
  !conjugates to the row, of the Hermitian border-border block, Gs is
  !obtained by one distributed base solve, and the row and column
  !of BG, the row of BZ and the row and column of the Schur complement that
  !belong to the entry are updated.
  !
  !With dummy true the entry is a deleted base slot and receives its
  !decoupling column: zero off the diagonal and the real part of the base
  !diagonal element on the diagonal (one if that value is zero).
  !
  !  Input parameters:
  !    p     - Border position of the entry.
  !    a     - Shifted conceptual column of length n (current order): a(i)
  !            is the element (i,idx). Not referenced when dummy is true.
  !    dummy - Selects the decoupling column of a deleted base slot.
  !
  !  Input/output parameter:
  !    self - State. The entry is no longer stale; tail_ok (entry in the tail)
  !           or scc_ok (committed entry) becomes false.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: p
    complex(wp), intent(in) :: a(self%n)
    logical, intent(in) :: dummy
    integer :: clean_rows, i, idx, idx2, p2, q, slot, slot2
    complex(wp) :: current, value

    slot = self%bslot(p)
    idx = self%bidx(p)
    q = 0
    if (idx <= self%n0) q = self%d0pos(idx)
    clean_rows = min(self%n, self%n0)
    do i = 1, clean_rows
      if (self%d0pos(i) /= 0) then
        self%bt(i, slot) = 0.0_wp
      else
        value = 0.0_wp
        if (.not. dummy) value = a(i)
        if (q > 0) value = value - self%ab(i, q)
        self%bt(i, slot) = value
      end if
    end do
    !Base rows beyond the current order are deleted slots, members of D0.
    do i = clean_rows + 1, self%n0
      self%bt(i, slot) = 0.0_wp
    end do
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      if (p2 == p) then
        if (dummy) then
          current = cmplx(real(self%ab(idx, q), wp), 0.0_wp, kind=wp)
          if (current == 0.0_wp) current = 1.0_wp
        else
          current = cmplx(real(a(idx), wp), 0.0_wp, kind=wp)
        end if
      else if (dummy .or. idx2 > self%n) then
        !A deleted slot is decoupled from everything else.
        current = 0.0_wp
      else
        current = a(idx2)
      end if
      !a(idx2) is the element (idx2,idx) of the conceptual column idx; the
      !element (idx,idx2) of the Hermitian matrix is its conjugate.
      self%cbb(slot2, slot) = current
      self%cbb(slot, slot2) = conjg(current)
    end do

    call complex_base_solve(self, self%bt(:, slot), self%gs(:, slot))

    !Row and column of BG and row of BZ.
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      self%bg(slot, slot2) = complex_dotc(self%n0, self%bt(:, slot), &
                                          self%gs(:, slot2))
      if (slot2 /= slot) then
        self%bg(slot2, slot) = complex_dotc(self%n0, self%bt(:, slot2), &
                                            self%gs(:, slot))
      end if
    end do
    do q = 1, self%m0
      self%bz(slot, q) = complex_dotc(self%n0, self%bt(:, slot), self%zb(:, q))
    end do
    call complex_fill_schur(self, p)

    self%stale(slot) = 0
    if (p > self%m - self%ma) then
      self%tail_ok = .false.
    else
      self%scc_ok = .false.
    end if
  end subroutine complex_install_column

  subroutine complex_fill_schur(self, p)
  !Subroutine complex_fill_schur recomputes the row and the column of the Schur
  !complement Sc that belong to the border entry at position p from the
  !stored small matrices. With u = mmax+slot the unknown of the entry, own(q)
  !equal to one when member q of D0 is the entry itself and zero otherwise,
  !and base(i,j) the base matrix element of two members of D0 (zero when
  !either entry is an appended function),
  !
  !        Sc(q,u)   = -own(q) - Gs(D0(q),slot)
  !        Sc(u,q)   = -own(q) - BZ(slot,q)
  !        Sc(u,u2)  = C(slot,slot2) - base(idx,idx2)  - BG(slot,slot2)
  !        Sc(u2,u)  = C(slot2,slot) - base(idx2,idx)  - BG(slot2,slot) .
  !
  !BG and BZ hold inner products that conjugate the border columns Bt, and C
  !is Hermitian. The Schur complement is treated as a general complex matrix:
  !the base elements recovered from the factors are Hermitian only to
  !roundoff.
  !
  !  Input parameter:
  !    p    - Border position of the entry.
  !
  !  Input/output parameter:
  !    self - State with valid cbb, bg, bz, gs and ab for the entry.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: p
    integer :: idx, idx2, p2, q, q2, qs, slot, slot2, u, u2
    complex(wp) :: base12, base21
    real(wp) :: own

    slot = self%bslot(p)
    idx = self%bidx(p)
    u = self%mmax + slot
    qs = 0
    if (idx <= self%n0) qs = self%d0pos(idx)
    do q = 1, self%m0
      own = 0.0_wp
      if (q == qs) own = 1.0_wp
      self%sc(q, u) = -own - self%gs(self%d0idx(q), slot)
      self%sc(u, q) = -own - self%bz(slot, q)
    end do
    do p2 = 1, self%m
      slot2 = self%bslot(p2)
      idx2 = self%bidx(p2)
      u2 = self%mmax + slot2
      base12 = 0.0_wp
      base21 = 0.0_wp
      if (qs > 0 .and. idx2 <= self%n0) then
        q2 = self%d0pos(idx2)
        base12 = self%ab(idx, q2)
        base21 = self%ab(idx2, qs)
      end if
      self%sc(u, u2) = self%cbb(slot, slot2) - base12 - self%bg(slot, slot2)
      self%sc(u2, u) = self%cbb(slot2, slot) - base21 - self%bg(slot2, slot)
    end do
  end subroutine complex_fill_schur

  subroutine complex_rebuild_small(self)
  !Subroutine complex_rebuild_small recomputes all small matrices from the
  !stored border columns: BG = Bt^H*Gs and BZ = Bt^H*Zb by a distributed pass
  !in which every rank takes a contiguous slice of rows (one reduction), and
  !the whole Schur complement Sc. It is called whenever D0 grows, which changes
  !the stored columns of all border entries; it is never called per trial.
  !
  !  Input/output parameter:
  !    self - State. The factors of Scc and of the tail are marked out of
  !           date.
    class(qr_complex_state), intent(inout) :: self
    integer :: chunk, count, hi, k, lo, m, m0, mpi_err, p, p2, q, q2
    integer :: slot, slot2

    m = self%m
    m0 = self%m0
    self%scc_ok = .false.
    self%tail_ok = .false.

    if (self%nprocs > 1) then
      chunk = (self%n0 + self%nprocs - 1) / self%nprocs
      lo = self%rank * chunk + 1
      hi = min(self%n0, (self%rank + 1) * chunk)
    else
      lo = 1
      hi = self%n0
    end if
    count = m * m + m * m0
    if (count > 0) then
      k = 0
      do p2 = 1, m
        slot2 = self%bslot(p2)
        do p = 1, m
          slot = self%bslot(p)
          k = k + 1
          if (hi >= lo) then
            self%wred(k) = complex_dotc(hi - lo + 1, self%bt(lo:hi, slot), &
                                        self%gs(lo:hi, slot2))
          else
            self%wred(k) = 0.0_wp
          end if
        end do
      end do
      do q = 1, m0
        do p = 1, m
          slot = self%bslot(p)
          k = k + 1
          if (hi >= lo) then
            self%wred(k) = complex_dotc(hi - lo + 1, self%bt(lo:hi, slot), &
                                        self%zb(lo:hi, q))
          else
            self%wred(k) = 0.0_wp
          end if
        end do
      end do
      if (self%nprocs > 1) then
        call MPI_Allreduce(MPI_IN_PLACE, self%wred, 2 * count, MPI_WP, &
                           MPI_SUM, self%comm, mpi_err)
      end if
      k = 0
      do p2 = 1, m
        slot2 = self%bslot(p2)
        do p = 1, m
          k = k + 1
          self%bg(self%bslot(p), slot2) = self%wred(k)
        end do
      end do
      do q = 1, m0
        do p = 1, m
          k = k + 1
          self%bz(self%bslot(p), q) = self%wred(k)
        end do
      end do
    end if

    !Unit-vector block of the Schur complement, -E^H*Zb, and the rows and
    !columns of all border entries.
    do q2 = 1, m0
      do q = 1, m0
        self%sc(q, q2) = -self%zb(self%d0idx(q), q2)
      end do
    end do
    do p = 1, m
      call complex_fill_schur(self, p)
    end do
  end subroutine complex_rebuild_small

  subroutine complex_refresh_committed(self)
  !Subroutine complex_refresh_committed brings the committed part of the small
  !system up to date: it installs the decoupling columns of committed deleted
  !slots that are still stale, and recomputes the LU factorization of the
  !committed block Scc when it is out of date. The committed unknowns are the
  !unit-vector unknowns of all members of D0 followed by the committed border
  !entries in border order.
  !
  !  Input/output parameter:
  !    self - State. scc_singular records the result of the relative pivot
  !           test; scc_ok becomes true when the factors are usable.
    class(qr_complex_state), intent(inout) :: self
    integer :: k, k2, mc, ncu, p

    if (self%r_singular) return
    mc = self%m - self%ma
    do p = 1, mc
      if (self%stale(self%bslot(p)) /= 0) then
        call complex_install_column(self, p, self%wn(:, 3), .true.)
      end if
    end do
    if (self%scc_ok) return
    ncu = self%m0 + mc
    do k = 1, self%m0
      self%cu(k) = k
    end do
    do p = 1, mc
      self%cu(self%m0 + p) = self%mmax + self%bslot(p)
    end do
    do k2 = 1, ncu
      do k = 1, ncu
        self%scc_lu(k, k2) = self%sc(self%cu(k), self%cu(k2))
      end do
    end do
    call complex_lu_factor(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                           QR_RELATIVE_PIVOT, self%scc_mag, self%scc_singular)
    self%scc_ok = .not. self%scc_singular
    self%tail_ok = .false.
  end subroutine complex_refresh_committed

  subroutine complex_refresh_small(self)
  !Subroutine complex_refresh_small brings the lazily maintained factors up to
  !date before a bordered solve: the committed block (see
  !complex_refresh_committed), the decoupling columns of stale deleted slots in
  !the tail, and for the active tail X = Scc^{-1}*Sca, Sigma = Saa - Sac*X
  !and LU(Sigma) if a tail column changed. All work is O(m**2*Ma) or less,
  !except a refactorization of Scc and one base solve per stale slot.
  !
  !Sigma is the Schur complement of everything else in the nearly singular
  !shifted matrix and is therefore legitimately tiny; its pivots are rejected
  !only when they are zero, subnormal, or not finite.
  !
  !  Input/output parameter:
  !    self - State; see complex_can_solve for the resulting flags.
    class(qr_complex_state), intent(inout) :: self
    integer :: a, a2, k, mc, ncu, p, u
    complex(wp) :: tsum

    if (self%r_singular) return
    mc = self%m - self%ma
    do p = mc + 1, self%m
      if (self%stale(self%bslot(p)) /= 0) then
        call complex_install_column(self, p, self%wn(:, 3), .true.)
      end if
    end do
    call complex_refresh_committed(self)
    if (.not. self%scc_ok) return
    if (self%tail_ok) return
    ncu = self%m0 + mc
    do a = 1, self%ma
      u = self%mmax + self%bslot(mc + a)
      do k = 1, ncu
        self%xca(k, a) = self%sc(self%cu(k), u)
      end do
      if (ncu > 0) then
        call complex_lu_solve(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                              self%xca(:, a))
      end if
    end do
    do a2 = 1, self%ma
      do a = 1, self%ma
        u = self%mmax + self%bslot(mc + a)
        tsum = 0.0_wp
        do k = 1, ncu
          tsum = tsum + self%sc(u, self%cu(k)) * self%xca(k, a2)
        end do
        self%sig_lu(a, a2) = self%sc(u, self%mmax + self%bslot(mc + a2)) - tsum
      end do
    end do
    call complex_lu_factor(self%ma, self%sig_lu, self%mmax, self%sig_piv, &
                           0.0_wp, self%scc_mag, self%sig_singular)
    self%tail_ok = .not. self%sig_singular
  end subroutine complex_refresh_small

  !===========================================================================
  ! Internal procedures: distributed kernels
  !===========================================================================

  subroutine complex_bordered_solve(self, r, x)
  !Subroutine complex_bordered_solve solves A*x = r for the represented matrix
  !of the current order n by block elimination of the augmented system
  !described in the module header:
  !
  !        r0 = r on the base indices, zero on deleted slots
  !        g  = A_base^{-1}*r0
  !        t  = [0; 0; rA] - [E D B]^H*g
  !        y  = Sc^{-1}*t      (block elimination on the committed block)
  !        x0 = g - Zb*yE - Gs*[yD; xA] .
  !
  !One distributed base solve is performed; the remaining work,
  !O((2*m+m0)*n0+(m+m0)**2), is executed identically on every rank.
  !
  !  Input parameter:
  !    r - Replicated right-hand side of length n. It must not be associated
  !        with work vector 1 or 2 of the state.
  !
  !  Input/output parameter:
  !    self - State for which complex_can_solve is true. Work vectors 1 and 2
  !           and the small work vectors are overwritten.
  !
  !  Output parameter:
  !    x - Replicated solution of length n; not associated with work vector 1
  !        or 2 of the state.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: r(self%n)
    complex(wp), intent(out) :: x(self%n)
    integer :: a, clean_rows, i, k, m, m0, mc, ncu, p, q, slot, u
    complex(wp) :: tsum

    m = self%m
    m0 = self%m0
    mc = m - self%ma
    ncu = m0 + mc
    clean_rows = min(self%n, self%n0)

    !Right-hand side on the base indices and base solve. Work vector 1 is
    !r0, work vector 2 is g and finally x0.
    do i = 1, clean_rows
      self%wn(i, 1) = r(i)
    end do
    do i = clean_rows + 1, self%n0
      self%wn(i, 1) = 0.0_wp
    end do
    call complex_base_solve(self, self%wn(:, 1), self%wn(:, 2))

    if (m > 0) then
      !Right-hand side of the small system: committed unknowns in small
      !vector 1 (order of cu), tail unknowns in small vector 2.
      do q = 1, m0
        self%wsm(q, 1) = -self%wn(self%d0idx(q), 2)
      end do
      do p = 1, m
        slot = self%bslot(p)
        tsum = -complex_dotc(self%n0, self%bt(:, slot), self%wn(:, 2))
        if (self%bidx(p) > self%n0) tsum = tsum + r(self%bidx(p))
        if (p <= mc) then
          self%wsm(m0 + p, 1) = tsum
        else
          self%wsm(p - mc, 2) = tsum
        end if
      end do

      !yc = Scc^{-1}*tc, ya = Sigma^{-1}*(ta - Sac*yc), yc = yc - X*ya.
      if (ncu > 0) then
        call complex_lu_solve(ncu, self%scc_lu, 2 * self%mmax, self%scc_piv, &
                              self%wsm(:, 1))
      end if
      do a = 1, self%ma
        u = self%mmax + self%bslot(mc + a)
        tsum = self%wsm(a, 2)
        do k = 1, ncu
          tsum = tsum - self%sc(u, self%cu(k)) * self%wsm(k, 1)
        end do
        self%wsm(a, 2) = tsum
      end do
      if (self%ma > 0) then
        call complex_lu_solve(self%ma, self%sig_lu, self%mmax, self%sig_piv, &
                              self%wsm(:, 2))
        do a = 1, self%ma
          tsum = self%wsm(a, 2)
          do k = 1, ncu
            self%wsm(k, 1) = self%wsm(k, 1) - self%xca(k, a) * tsum
          end do
        end do
      end if

      !x0 = g - Zb*yE - Gs*y, with the negated coefficients in small vector 3.
      if (m0 > 0) then
        do q = 1, m0
          self%wsm(q, 3) = -self%wsm(q, 1)
        end do
        call complex_add_leading_columns(self%n0, self%zb, self%capacity, m0, &
                                         self%wsm(:, 3), self%wn(:, 2))
      end if
      if (mc > 0) then
        do p = 1, mc
          self%wsm(p, 3) = -self%wsm(m0 + p, 1)
        end do
        call complex_add_columns(self%n0, self%gs, self%capacity, mc, &
                                 self%bslot(1:mc), self%wsm(:, 3), &
                                 self%wn(:, 2))
      end if
      if (self%ma > 0) then
        do a = 1, self%ma
          self%wsm(a, 3) = -self%wsm(a, 2)
        end do
        call complex_add_columns(self%n0, self%gs, self%capacity, self%ma, &
                                 self%bslot(mc + 1:m), self%wsm(:, 3), &
                                 self%wn(:, 2))
      end if
    end if

    do i = 1, clean_rows
      x(i) = self%wn(i, 2)
    end do
    do p = 1, m
      if (self%bidx(p) > self%n0) then
        if (p <= mc) then
          x(self%bidx(p)) = self%wsm(m0 + p, 1)
        else
          x(self%bidx(p)) = self%wsm(p - mc, 2)
        end if
      end if
    end do
  end subroutine complex_bordered_solve

  subroutine complex_base_solve(self, b, x)
  !Subroutine complex_base_solve computes x = A_base^{-1}*b = R^{-1}*(Q^H*b)
  !with the distributed base factors.
  !
  !  1. y = Q^H*b: every rank accumulates b(i) times the complex conjugates
  !     of its rows of Q (the contiguous columns of qt, which stores columns
  !     of Q^H), followed by one reduction of n0 complex values.
  !  2. Back substitution by blocks of P rows from the bottom. Every rank owns
  !     exactly one row of a block and computes the part of its row of R to
  !     the right of the block times the already known solution (one dot
  !     product on a contiguous column of rt). One reduction of at most P
  !     values makes these products global; then every rank solves the P by P
  !     triangular block from the replicated diagonal blocks. The bottom block
  !     needs no reduction.
  !
  !With one process the routine performs an ordinary back substitution and
  !makes no MPI call.
  !
  !  Input parameter:
  !    b - Replicated right-hand side of length n0.
  !
  !  Input/output parameter:
  !    self - State with usable R (r_singular is false). The local work
  !           arrays wloc and wblk are overwritten.
  !
  !  Output parameter:
  !    x - Replicated solution of length n0. It must not overlap b.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: b(self%n0)
    complex(wp), intent(out) :: x(self%n0)
    integer :: blk, i, ie, ii, is, jj, l, mpi_err, n0, nb, np
    complex(wp) :: tsum

    n0 = self%n0
    np = self%nprocs

    !y = Q^H*b, accumulated in x. The local columns of qt are columns of Q^H,
    !so no conjugation is needed here.
    do l = 1, self%nloc
      self%wloc(l) = b(self%rank + 1 + (l - 1) * np)
    end do
    do i = 1, n0
      x(i) = 0.0_wp
    end do
    call complex_add_leading_columns(n0, self%qt, self%capacity, self%nloc, &
                                     self%wloc, x)
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, x, 2 * n0, MPI_WP, MPI_SUM, &
                         self%comm, mpi_err)
    end if

    if (np == 1) then
      do i = n0, 1, -1
        tsum = x(i)
        if (i < n0) then
          tsum = tsum - complex_dotu(n0 - i, self%rt(i + 1:n0, i), &
                                     x(i + 1:n0))
        end if
        x(i) = tsum / self%rt(i, i)
      end do
    else
      do blk = self%nblk, 1, -1
        is = (blk - 1) * np + 1
        ie = min(blk * np, n0)
        nb = ie - is + 1
        if (ie < n0) then
          do ii = 1, nb
            self%wblk(ii) = 0.0_wp
          end do
          if (self%rank < nb) then
            self%wblk(self%rank + 1) = complex_dotu(n0 - ie, &
                                                    self%rt(ie + 1:n0, blk), &
                                                    x(ie + 1:n0))
          end if
          call MPI_Allreduce(MPI_IN_PLACE, self%wblk, 2 * nb, MPI_WP, &
                             MPI_SUM, self%comm, mpi_err)
          do ii = 1, nb
            x(is + ii - 1) = x(is + ii - 1) - self%wblk(ii)
          end do
        end if
        do ii = nb, 1, -1
          tsum = x(is + ii - 1)
          do jj = ii + 1, nb
            tsum = tsum - self%rd(ii, jj, blk) * x(is + jj - 1)
          end do
          x(is + ii - 1) = tsum / self%rd(ii, ii, blk)
        end do
      end do
    end if
  end subroutine complex_base_solve

  subroutine complex_base_column(self, idx, column, work)
  !Subroutine complex_base_column recovers column idx of the base matrix from
  !the distributed factors, column = Q*(R*e_idx). Every rank contributes its
  !entries of column idx of R (one reduction) and then the products of its
  !rows of Q, the conjugates of the stored columns of qt, with that column
  !(one reduction).
  !
  !  Input parameter:
  !    idx - Column index, 1 <= idx <= n0.
  !
  !  Input/output parameter:
  !    self - State with valid base factors (not modified; declared
  !           intent(inout) because the output arguments are normally work
  !           vectors owned by the state).
  !
  !  Output parameters:
  !    column - Replicated column of length n0.
  !    work   - Work vector of length n0 (column idx of R on exit).
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(out) :: column(self%n0), work(self%n0)
    integer :: i, l, mpi_err, n0, np

    n0 = self%n0
    np = self%nprocs
    do i = 1, n0
      work(i) = 0.0_wp
      column(i) = 0.0_wp
    end do
    do l = 1, self%nloc
      i = self%rank + 1 + (l - 1) * np
      if (i <= idx) work(i) = self%rt(idx, l)
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, work, 2 * n0, MPI_WP, MPI_SUM, &
                         self%comm, mpi_err)
    end if
    do l = 1, self%nloc
      i = self%rank + 1 + (l - 1) * np
      !Row i of Q is the conjugate of the stored column of Q^H.
      column(i) = complex_dotc(idx, self%qt(1:idx, l), work(1:idx))
    end do
    if (np > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, column, 2 * n0, MPI_WP, MPI_SUM, &
                         self%comm, mpi_err)
    end if
  end subroutine complex_base_column

  subroutine complex_sym_product(rank, nprocs, comm, n, a, lda, x, y)
  !Subroutine complex_sym_product computes y = A*x for a complex Hermitian
  !matrix of which only the lower triangle is stored (and replicated on all
  !ranks). The columns are assigned to the ranks cyclically; a rank accumulates
  !both the column contribution (below the diagonal) and the conjugated row
  !contribution of its columns in one pass over each column, and one reduction
  !of n values completes the product. With one process no MPI call is made.
  !
  !The owned columns are processed several at a time over their common rows,
  !so that y is read and written once per group; the few rows above the
  !common range are handled column by column. Groups of four are used except
  !for 80-bit extended precision (QR_PRODUCT_PAIRS), where the pairs of the
  !real kernel are kept: for complex data no significant difference between
  !pairs and groups of four was measured, and single columns were slower.
  !
  !  Input parameters:
  !    rank, nprocs, comm - Own rank, number of ranks, and communicator.
  !    n   - Order of the matrix.
  !    a   - Array with leading dimension lda whose leading n by n lower
  !          triangle holds the matrix. The upper triangle is not referenced,
  !          and the imaginary parts of the diagonal elements are ignored.
  !    lda - Leading dimension of a.
  !    x   - Replicated vector of length n.
  !
  !  Output parameter:
  !    y   - Replicated product of length n. It must not overlap x.
    integer, intent(in) :: rank, nprocs, comm, n, lda
    complex(wp), intent(in) :: a(lda, *), x(n)
    complex(wp), intent(out) :: y(n)
    integer :: i, j, j2, j3, j4, mpi_err
    complex(wp) :: t1, t2, t3, t4, u1, u2, u3

    do i = 1, n
      y(i) = 0.0_wp
    end do
    j = rank + 1
    if (QR_PRODUCT_PAIRS) then
      !Pairs of owned columns j < j2 with common rows j2+1..n.
      do while (j + nprocs < n)
        j2 = j + nprocs
        call complex_symv_column(j2 - j, a(j + 1, j), x(j), x(j + 1), &
                                 y(j + 1), u1)
        call complex_symv_columns2(n - j2, a(j2 + 1, j), a(j2 + 1, j2), x(j), &
                                   x(j2), x(j2 + 1), y(j2 + 1), t1, t2)
        y(j) = y(j) + real(a(j, j), wp) * x(j) + (u1 + t1)
        y(j2) = y(j2) + real(a(j2, j2), wp) * x(j2) + t2
        j = j2 + nprocs
      end do
    else
      !Groups of four owned columns j < j2 < j3 < j4 with common rows
      !j4+1..n.
      do while (j + 3 * nprocs < n)
        j2 = j + nprocs
        j3 = j2 + nprocs
        j4 = j3 + nprocs
        call complex_symv_column(j4 - j, a(j + 1, j), x(j), x(j + 1), &
                                 y(j + 1), u1)
        call complex_symv_column(j4 - j2, a(j2 + 1, j2), x(j2), x(j2 + 1), &
                                 y(j2 + 1), u2)
        call complex_symv_column(j4 - j3, a(j3 + 1, j3), x(j3), x(j3 + 1), &
                                 y(j3 + 1), u3)
        call complex_symv_columns4(n - j4, a(j4 + 1, j), a(j4 + 1, j2), &
                                   a(j4 + 1, j3), a(j4 + 1, j4), x(j), x(j2), &
                                   x(j3), x(j4), x(j4 + 1), y(j4 + 1), t1, t2, &
                                   t3, t4)
        y(j) = y(j) + real(a(j, j), wp) * x(j) + (u1 + t1)
        y(j2) = y(j2) + real(a(j2, j2), wp) * x(j2) + (u2 + t2)
        y(j3) = y(j3) + real(a(j3, j3), wp) * x(j3) + (u3 + t3)
        y(j4) = y(j4) + real(a(j4, j4), wp) * x(j4) + t4
        j = j4 + nprocs
      end do
    end if
    !Remaining owned columns one at a time.
    do while (j <= n)
      if (j < n) then
        call complex_symv_column(n - j, a(j + 1, j), x(j), x(j + 1), y(j + 1), &
                                 t1)
        y(j) = y(j) + real(a(j, j), wp) * x(j) + t1
      else
        y(j) = y(j) + real(a(j, j), wp) * x(j)
      end if
      j = j + nprocs
    end do
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, y, 2 * n, MPI_WP, MPI_SUM, comm, &
                         mpi_err)
    end if
  end subroutine complex_sym_product

  !===========================================================================
  ! Internal procedures: local kernels
  !===========================================================================

  pure subroutine complex_symv_column(n, a, xj, x, y, tsum)
  !Subroutine complex_symv_column processes the part of one column of a
  !Hermitian matrix below the diagonal: it adds xj*a to y (column contribution)
  !and returns tsum = a^H*x (conjugated row contribution) in a single pass.
  !Four independent partial sums shorten the dependency chain of the dot
  !product.
  !
  !  Input parameters:
  !    n  - Number of elements.
  !    a  - Column part of length n.
  !    xj - Multiplier of the column.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameter:
  !    tsum - a^H*x.
    integer, intent(in) :: n
    complex(wp), intent(in) :: a(n), xj, x(n)
    complex(wp), intent(inout) :: y(n)
    complex(wp), intent(out) :: tsum
    integer :: i, n4
    complex(wp) :: t1, t2, t3, t4

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    n4 = n - mod(n, 4)
    do i = 1, n4, 4
      y(i) = y(i) + xj * a(i)
      y(i + 1) = y(i + 1) + xj * a(i + 1)
      y(i + 2) = y(i + 2) + xj * a(i + 2)
      y(i + 3) = y(i + 3) + xj * a(i + 3)
      t1 = t1 + conjg(a(i)) * x(i)
      t2 = t2 + conjg(a(i + 1)) * x(i + 1)
      t3 = t3 + conjg(a(i + 2)) * x(i + 2)
      t4 = t4 + conjg(a(i + 3)) * x(i + 3)
    end do
    do i = n4 + 1, n
      y(i) = y(i) + xj * a(i)
      t1 = t1 + conjg(a(i)) * x(i)
    end do
    tsum = (t1 + t2) + (t3 + t4)
  end subroutine complex_symv_column

  pure subroutine complex_symv_columns4(n, a1, a2, a3, a4, x1, x2, x3, x4, x, &
                                        y, t1, t2, t3, t4)
  !Subroutine complex_symv_columns4 processes the common rows of four columns
  !of a Hermitian matrix below their diagonals in a single pass: it adds
  !x1*a1+x2*a2+x3*a3+x4*a4 to y (column contributions) and returns the four
  !inner products tk = ak^H*x (conjugated row contributions).
  !
  !  Input parameters:
  !    n  - Number of rows.
  !    a1, a2, a3, a4 - Column parts of length n.
  !    x1, x2, x3, x4 - Multipliers of the columns.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameters:
  !    t1, t2, t3, t4 - The four inner products.
    integer, intent(in) :: n
    complex(wp), intent(in) :: a1(n), a2(n), a3(n), a4(n), x1, x2, x3, x4
    complex(wp), intent(in) :: x(n)
    complex(wp), intent(inout) :: y(n)
    complex(wp), intent(out) :: t1, t2, t3, t4
    integer :: i
    complex(wp) :: xi

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    do i = 1, n
      xi = x(i)
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i) + x3 * a3(i) + x4 * a4(i)
      t1 = t1 + conjg(a1(i)) * xi
      t2 = t2 + conjg(a2(i)) * xi
      t3 = t3 + conjg(a3(i)) * xi
      t4 = t4 + conjg(a4(i)) * xi
    end do
  end subroutine complex_symv_columns4

  pure subroutine complex_symv_columns2(n, a1, a2, x1, x2, x, y, t1, t2)
  !Subroutine complex_symv_columns2 processes the common rows of two columns of
  !a Hermitian matrix below their diagonals in a single pass: it adds
  !x1*a1+x2*a2 to y (column contributions) and returns the two inner products
  !tk = ak^H*x (conjugated row contributions). Two partial sums per inner
  !product shorten the dependency chains.
  !
  !  Input parameters:
  !    n  - Number of rows.
  !    a1, a2 - Column parts of length n.
  !    x1, x2 - Multipliers of the columns.
  !    x  - Vector part of length n.
  !
  !  Input/output parameter:
  !    y  - Accumulated product part of length n.
  !
  !  Output parameters:
  !    t1, t2 - The two inner products.
    integer, intent(in) :: n
    complex(wp), intent(in) :: a1(n), a2(n), x1, x2, x(n)
    complex(wp), intent(inout) :: y(n)
    complex(wp), intent(out) :: t1, t2
    integer :: i, n2
    complex(wp) :: s1, s2

    t1 = 0.0_wp
    t2 = 0.0_wp
    s1 = 0.0_wp
    s2 = 0.0_wp
    n2 = n - mod(n, 2)
    do i = 1, n2, 2
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i)
      y(i + 1) = y(i + 1) + x1 * a1(i + 1) + x2 * a2(i + 1)
      t1 = t1 + conjg(a1(i)) * x(i)
      t2 = t2 + conjg(a2(i)) * x(i)
      s1 = s1 + conjg(a1(i + 1)) * x(i + 1)
      s2 = s2 + conjg(a2(i + 1)) * x(i + 1)
    end do
    do i = n2 + 1, n
      y(i) = y(i) + x1 * a1(i) + x2 * a2(i)
      t1 = t1 + conjg(a1(i)) * x(i)
      t2 = t2 + conjg(a2(i)) * x(i)
    end do
    t1 = t1 + s1
    t2 = t2 + s2
  end subroutine complex_symv_columns2

  pure function complex_dotc(n, x, y) result(dot)
  !Function complex_dotc computes the Hermitian inner product x^H*y (the first
  !argument is conjugated) of the first n elements with four independent
  !partial sums.
  !
  !  Input parameters:
  !    n - Number of elements.
  !    x - Complex vector of length n; conjugated.
  !    y - Complex vector of length n.
  !
  !  Result:
  !    dot - x^H*y.
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(n), y(n)
    complex(wp) :: dot
    integer :: i, n4
    complex(wp) :: t1, t2, t3, t4

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    n4 = n - mod(n, 4)
    do i = 1, n4, 4
      t1 = t1 + conjg(x(i)) * y(i)
      t2 = t2 + conjg(x(i + 1)) * y(i + 1)
      t3 = t3 + conjg(x(i + 2)) * y(i + 2)
      t4 = t4 + conjg(x(i + 3)) * y(i + 3)
    end do
    do i = n4 + 1, n
      t1 = t1 + conjg(x(i)) * y(i)
    end do
    dot = (t1 + t2) + (t3 + t4)
  end function complex_dotc

  pure function complex_dotu(n, x, y) result(dot)
  !Function complex_dotu computes the unconjugated product x^T*y of the
  !first n elements of two complex vectors with four independent partial
  !sums. It is the kernel of the back substitution with the rows of R.
  !
  !  Input parameters:
  !    n - Number of elements.
  !    x - Complex vector of length n.
  !    y - Complex vector of length n.
  !
  !  Result:
  !    dot - x^T*y (no conjugation).
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(n), y(n)
    complex(wp) :: dot
    integer :: i, n4
    complex(wp) :: t1, t2, t3, t4

    t1 = 0.0_wp
    t2 = 0.0_wp
    t3 = 0.0_wp
    t4 = 0.0_wp
    n4 = n - mod(n, 4)
    do i = 1, n4, 4
      t1 = t1 + x(i) * y(i)
      t2 = t2 + x(i + 1) * y(i + 1)
      t3 = t3 + x(i + 2) * y(i + 2)
      t4 = t4 + x(i + 3) * y(i + 3)
    end do
    do i = n4 + 1, n
      t1 = t1 + x(i) * y(i)
    end do
    dot = (t1 + t2) + (t3 + t4)
  end function complex_dotu

  pure subroutine complex_axpy(n, alpha, x, y)
  !Subroutine complex_axpy adds alpha*x to y for the first n elements.
  !
  !  Input parameters:
  !    n     - Number of elements.
  !    alpha - Multiplier.
  !    x     - Complex vector of length n.
  !
  !  Input/output parameter:
  !    y     - Complex vector of length n.
    integer, intent(in) :: n
    complex(wp), intent(in) :: alpha, x(n)
    complex(wp), intent(inout) :: y(n)
    integer :: i

    do i = 1, n
      y(i) = y(i) + alpha * x(i)
    end do
  end subroutine complex_axpy

  pure subroutine complex_add_leading_columns(n, a, lda, ncols, coef, y)
  !Subroutine complex_add_leading_columns adds a linear combination of the
  !first ncols columns of a to y,
  !
  !        y(1:n) = y(1:n) + sum_k coef(k)*a(1:n,k) ,
  !
  !four columns at a time so that y is traversed once per four columns.
  !
  !  Input parameters:
  !    n     - Number of rows.
  !    a     - Matrix with leading dimension lda.
  !    lda   - Leading dimension of a.
  !    ncols - Number of columns.
  !    coef  - Coefficients, length ncols.
  !
  !  Input/output parameter:
  !    y     - Accumulated vector of length n.
    integer, intent(in) :: n, lda, ncols
    complex(wp), intent(in) :: a(lda, *), coef(*)
    complex(wp), intent(inout) :: y(n)
    integer :: i, k, k4
    complex(wp) :: c1, c2, c3, c4

    k4 = ncols - mod(ncols, 4)
    do k = 1, k4, 4
      c1 = coef(k)
      c2 = coef(k + 1)
      c3 = coef(k + 2)
      c4 = coef(k + 3)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, k) + c2 * a(i, k + 1) + c3 * a(i, k + 2) + &
               c4 * a(i, k + 3)
      end do
    end do
    do k = k4 + 1, ncols
      c1 = coef(k)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, k)
      end do
    end do
  end subroutine complex_add_leading_columns

  pure subroutine complex_add_columns(n, a, lda, ncols, cols, coef, y)
  !Subroutine complex_add_columns adds a linear combination of selected columns
  !of a to y,
  !
  !        y(1:n) = y(1:n) + sum_k coef(k)*a(1:n,cols(k)) ,
  !
  !four columns at a time.
  !
  !  Input parameters:
  !    n     - Number of rows.
  !    a     - Matrix with leading dimension lda.
  !    lda   - Leading dimension of a.
  !    ncols - Number of selected columns.
  !    cols  - Column indices, length ncols.
  !    coef  - Coefficients, length ncols.
  !
  !  Input/output parameter:
  !    y     - Accumulated vector of length n.
    integer, intent(in) :: n, lda, ncols, cols(*)
    complex(wp), intent(in) :: a(lda, *), coef(*)
    complex(wp), intent(inout) :: y(n)
    integer :: i, j1, j2, j3, j4, k, k4
    complex(wp) :: c1, c2, c3, c4

    k4 = ncols - mod(ncols, 4)
    do k = 1, k4, 4
      j1 = cols(k)
      j2 = cols(k + 1)
      j3 = cols(k + 2)
      j4 = cols(k + 3)
      c1 = coef(k)
      c2 = coef(k + 1)
      c3 = coef(k + 2)
      c4 = coef(k + 3)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, j1) + c2 * a(i, j2) + c3 * a(i, j3) + &
               c4 * a(i, j4)
      end do
    end do
    do k = k4 + 1, ncols
      j1 = cols(k)
      c1 = coef(k)
      do i = 1, n
        y(i) = y(i) + c1 * a(i, j1)
      end do
    end do
  end subroutine complex_add_columns

  pure subroutine complex_lu_factor(n, a, lda, piv, relative_tolerance, mag, &
                                    singular)
  !Subroutine complex_lu_factor computes the LU factorization with partial
  !pivoting of a small dense complex matrix in place, P*A = L*U, with unit
  !lower triangular L stored below the diagonal and the row interchanges
  !applied to whole rows (LAPACK convention).
  !
  !A pivot is rejected, and the matrix reported singular, when its magnitude
  !is zero, subnormal, or not finite. With a positive relative_tolerance a
  !pivot is also rejected when it has been lost to cancellation: the routine
  !accumulates for every element the magnitude sum
  !
  !        mag(i,j) = |a(i,j)| + sum_k |l(i,k)|*|u(k,j)|
  !
  !of the terms that were combined to produce it, and rejects a pivot whose
  !magnitude does not exceed relative_tolerance*mag. Such a pivot retains
  !fewer than about log10(1/relative_tolerance) of its digits relative to the
  !working precision. A pivot that is small without cancellation, as occurs
  !for badly scaled but well determined matrices, is accepted. Passing zero as
  !relative_tolerance disables the test and mag is not referenced.
  !
  !  Input parameters:
  !    n   - Order of the matrix; zero is allowed.
  !    lda - Leading dimension of a and mag.
  !    relative_tolerance - Relative pivot threshold (non-negative).
  !
  !  Input/output parameter:
  !    a   - On entry the matrix, on exit its factors. Undefined content is
  !          left behind when singular is true.
  !
  !  Output parameters:
  !    piv - Pivot row of every elimination step, length n.
  !    mag - Real work array (lda by n) of the magnitude sums.
  !    singular - True when a pivot was rejected.
    integer, intent(in) :: n, lda
    complex(wp), intent(inout) :: a(lda, *)
    integer, intent(out) :: piv(*)
    real(wp), intent(in) :: relative_tolerance
    real(wp), intent(inout) :: mag(lda, *)
    logical, intent(out) :: singular
    integer :: i, ip, j, k
    complex(wp) :: factor, swap
    real(wp) :: factor_magnitude, mag_swap, pivot_magnitude
    logical :: track

    singular = .false.
    track = (relative_tolerance > 0.0_wp)
    if (track) then
      do j = 1, n
        do i = 1, n
          mag(i, j) = abs(a(i, j))
        end do
      end do
    end if
    do k = 1, n
      ip = k
      pivot_magnitude = abs(a(k, k))
      do i = k + 1, n
        if (abs(a(i, k)) > pivot_magnitude) then
          ip = i
          pivot_magnitude = abs(a(i, k))
        end if
      end do
      piv(k) = ip
      !The negated form also rejects NaN.
      if (.not. (pivot_magnitude > tiny(1.0_wp) .and. &
                 pivot_magnitude <= huge(1.0_wp))) then
        singular = .true.
        return
      end if
      if (track) then
        if (.not. (pivot_magnitude > relative_tolerance * mag(ip, k))) then
          singular = .true.
          return
        end if
      end if
      if (ip /= k) then
        do j = 1, n
          swap = a(k, j)
          a(k, j) = a(ip, j)
          a(ip, j) = swap
        end do
        if (track) then
          do j = 1, n
            mag_swap = mag(k, j)
            mag(k, j) = mag(ip, j)
            mag(ip, j) = mag_swap
          end do
        end if
      end if
      do i = k + 1, n
        a(i, k) = a(i, k) / a(k, k)
      end do
      !The moduli of the multipliers are kept in the part of column k of mag
      !below the diagonal, which is not referenced again; the number of
      !complex moduli evaluated is then of second order in n.
      if (track) then
        do i = k + 1, n
          mag(i, k) = abs(a(i, k))
        end do
      end if
      do j = k + 1, n
        factor = a(k, j)
        do i = k + 1, n
          a(i, j) = a(i, j) - a(i, k) * factor
        end do
        if (track) then
          factor_magnitude = abs(factor)
          do i = k + 1, n
            mag(i, j) = mag(i, j) + mag(i, k) * factor_magnitude
          end do
        end if
      end do
    end do
  end subroutine complex_lu_factor

  pure subroutine complex_lu_solve(n, a, lda, piv, b)
  !Subroutine complex_lu_solve solves A*x = b in place with the factors
  !computed by complex_lu_factor.
  !
  !  Input parameters:
  !    n   - Order of the matrix; zero is allowed.
  !    a   - LU factors with leading dimension lda.
  !    lda - Leading dimension of a.
  !    piv - Pivot rows.
  !
  !  Input/output parameter:
  !    b   - Right-hand side on entry, solution on exit (length n).
    integer, intent(in) :: n, lda
    complex(wp), intent(in) :: a(lda, *)
    integer, intent(in) :: piv(*)
    complex(wp), intent(inout) :: b(*)
    integer :: i, k
    complex(wp) :: swap, value

    do k = 1, n
      if (piv(k) /= k) then
        swap = b(k)
        b(k) = b(piv(k))
        b(piv(k)) = swap
      end if
    end do
    do k = 1, n
      value = b(k)
      do i = k + 1, n
        b(i) = b(i) - a(i, k) * value
      end do
    end do
    do k = n, 1, -1
      b(k) = b(k) / a(k, k)
      value = b(k)
      do i = 1, k - 1
        b(i) = b(i) - a(i, k) * value
      end do
    end do
  end subroutine complex_lu_solve

  pure function complex_norm_squared(n, x) result(norm_squared)
  !Function complex_norm_squared computes the Hermitian Euclidean inner product
  !
  !              x^H*x = sum(real(x(i))^2+imag(x(i))^2)
  !
  !over the first n elements of a complex vector. The mathematically real
  !quantity is accumulated directly in real(wp). The result is the squared
  !norm; the square root is not taken.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Complex vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - The real value x^H*x for x(1:n).
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(*)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + real(x(i), wp)**2 + aimag(x(i))**2
    end do
  end function complex_norm_squared

  pure function complex_direction_difference(n, x, alpha, y, x_norm_squared) &
      result(relative_difference)
  !Function complex_direction_difference computes the scale-independent change
  !of direction between two complex vectors,
  !
  !                 ||x-alpha*y||_2 / ||x||_2 .
  !
  !In inverse iteration alpha=(y^H*x)/(y^H*y), so alpha*y is the component of
  !x parallel to the previous iterate y. The result measures only the remaining
  !directional change and is unaffected by a complex rescaling (in particular
  !a change of phase) of an eigenvector. Each squared magnitude is accumulated
  !in real(wp).
  !
  !  Input parameters:
  !    n              - Number of elements included in the calculation.
  !    x              - New complex iterate, containing at least n elements.
  !    alpha          - Complex projection coefficient multiplying y.
  !    y              - Previous complex iterate, containing at least n
  !                     elements.
  !    x_norm_squared - Precomputed positive real value x^H*x.
  !
  !  Result:
  !    relative_difference - Relative Euclidean norm shown above.
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(*), alpha, y(*)
    real(wp), intent(in) :: x_norm_squared
    real(wp) :: relative_difference
    complex(wp) :: difference
    real(wp) :: difference_norm_squared
    integer :: i

    difference_norm_squared = 0.0_wp
    do i = 1, n
      difference = x(i) - alpha * y(i)
      difference_norm_squared = difference_norm_squared + &
                                real(difference, wp)**2 + &
                                aimag(difference)**2
    end do
    relative_difference = sqrt(difference_norm_squared / x_norm_squared)
  end function complex_direction_difference

  pure function complex_diagonal_is_real(n, column, idx) result(is_real)
  !Function complex_diagonal_is_real tests whether the diagonal element of a
  !supplied conceptual column of a Hermitian matrix is real to within
  !roundoff, that is, whether
  !
  !        abs(imag(column(idx))) <= QR_DIAGONAL_TOLERANCE*scale,
  !        scale = max(1, max_i(abs(real(column(i)))+abs(imag(column(i))))) .
  !
  !The scale uses the sum of the absolute real and imaginary parts, which
  !bounds the modulus from above without a square root.
  !
  !  Input parameters:
  !    n      - Number of elements of the column.
  !    column - Complex column of length n.
  !    idx    - Position of the diagonal element, 1 <= idx <= n.
  !
  !  Result:
  !    is_real - True when the test above holds; false also when the imaginary
  !              part is not a number.
    integer, intent(in) :: n, idx
    complex(wp), intent(in) :: column(n)
    logical :: is_real
    integer :: i
    real(wp) :: scale

    scale = 1.0_wp
    do i = 1, n
      scale = max(scale, abs(real(column(i), wp)) + abs(aimag(column(i))))
    end do
    is_real = (abs(aimag(column(idx))) <= QR_DIAGONAL_TOLERANCE * scale)
  end function complex_diagonal_is_real

  pure function complex_reciprocal(z) result(w)
  !Function complex_reciprocal returns 1/z for a nonzero complex z by Smith's
  !formula, in which the real and imaginary parts are never squared. The
  !result is therefore accurate also when abs(z)**2 would overflow or
  !underflow, which the intrinsic complex division does not guarantee with
  !every compiler. It is used by the generation of Householder reflectors.
  !
  !  Input parameter:
  !    z - Nonzero complex number.
  !
  !  Result:
  !    w - 1/z.
    complex(wp), intent(in) :: z
    complex(wp) :: w
    real(wp) :: denominator, ratio, zi, zr

    zr = real(z, wp)
    zi = aimag(z)
    if (abs(zr) >= abs(zi)) then
      ratio = zi / zr
      denominator = zr + zi * ratio
      w = cmplx(1.0_wp / denominator, -ratio / denominator, kind=wp)
    else
      ratio = zr / zi
      denominator = zi + zr * ratio
      w = cmplx(ratio / denominator, -1.0_wp / denominator, kind=wp)
    end if
  end function complex_reciprocal

  pure function complex_max_abs_real_or_imag(n, x) result(max_component)
  !Function complex_max_abs_real_or_imag returns
  !
  !       max_i( max(abs(real(x(i))),abs(imag(x(i)))) )
  !
  !for the first n elements of x. This component norm is used to scale complex
  !inverse iterates. It differs from maxval(abs(x)), which uses the Euclidean
  !modulus of each complex element.
  !
  !  Input parameters:
  !    n - Number of vector elements to inspect.
  !    x - Complex vector containing at least n elements.
  !
  !  Result:
  !    max_component - Largest magnitude of any real or imaginary component.
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(*)
    real(wp) :: max_component
    integer :: i

    max_component = 0.0_wp
    do i = 1, n
      max_component = max(max_component, abs(real(x(i), wp)), &
                          abs(aimag(x(i))))
    end do
  end function complex_max_abs_real_or_imag

end module qrlinalg
