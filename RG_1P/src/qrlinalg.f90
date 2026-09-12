!Module qrlinalg provides serial QR factorizations and inverse iteration for
!the generalized symmetric and Hermitian eigenvalue problems
!
!                 H*x = lambda*S*x .
!
!A real state represents a symmetric problem and a complex state represents a
!Hermitian problem. For a real shift sigma, both states store the factorization
!
!                 H - sigma*S = Q*R ,
!
!where Q is explicit and orthogonal or unitary and R is upper triangular. The
!explicit representation is required by the structural QR-update routines and
!also permits inverse iteration without retaining a separate shifted matrix.
!
!The caller owns H, S, and all vectors passed through the public interface. A
!state does not retain pointers to caller arrays. It owns Q, R, the Householder
!coefficients used while constructing Q, and all factorization, update, and
!solve workspace. After initialization, factorization and inverse iteration do
!not allocate memory.
!
!The working kind wp is selected when the library is compiled. The bundled
!BLAS, LAPACK, qrupdate sources, interfaces in this module, and public arrays
!all use that same kind. No MPI object or process-dependent state is stored.
module qrlinalg
  use iso_fortran_env, only: int64
  use wp_def, only: wp
  use qrupdate, only: qr1up, qrinc, qrdec, qrinr, qrder
  implicit none
  private

  integer, parameter, public :: QR_SUCCESS = 0
    !The requested operation completed successfully.
  integer, parameter, public :: QR_ERR_INVALID_ARGUMENT = 1
    !A scalar control, index, or mathematical input property is invalid.
  integer, parameter, public :: QR_ERR_ALLOCATION = 2
    !Initialization could not allocate all required state storage.
  integer, parameter, public :: QR_ERR_NOT_IMPLEMENTED = 3
    !The requested numerical operation is not implemented in this version.
  integer, parameter, public :: QR_ERR_FACTORIZATION = 4
    !A bundled LAPACK factorization or workspace query reported an error.
  integer, parameter, public :: QR_ERR_SINGULAR = 5
    !The stored shifted factorization is singular or numerically unusable.
  integer, parameter, public :: QR_ERR_NO_CONVERGENCE = 6
    !Inverse iteration reached max_iter before satisfying the tolerance.
  integer, parameter, public :: QR_ERR_INVALID_STATE = 7
    !The state is uninitialized, invalid, or missing required owned storage.
  integer, parameter, public :: QR_ERR_DIMENSION_MISMATCH = 8
    !A caller-owned array has an incompatible or empty extent.
  integer, parameter, public :: QR_ERR_CAPACITY_EXCEEDED = 9
    !The requested active order is larger than the initialized capacity.
  integer, parameter, public :: QR_ERR_ZERO_INITIAL_VECTOR = 10
    !The supplied inverse-iteration starting vector is numerically zero.
  integer, parameter, public :: QR_ERR_NONPOSITIVE_OVERLAP = 11
    !The final vector has a non-positive or numerically zero overlap norm.
  integer, parameter :: QR_STATUS_MESSAGE_LENGTH = 96
    !Fixed result length used by qr_status_message; returned text is padded.

  public :: wp
  public :: qr_real_state
  public :: qr_complex_state
  public :: qr_status_message

  !Interfaces to the external BLAS and LAPACK routines supplied in
  !src/qrupdate. Explicit interfaces allow the compiler to verify argument
  !types and working kinds at every call site. The routine names retain their
  !traditional D/Z prefixes, but their scalar kind is the compile-time wp;
  !there is no assumption that a D routine always uses eight-byte reals.
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

    ! Multiply a real symmetric matrix by a vector using one stored triangle.
    subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda, incx, incy
      real(wp), intent(in) :: alpha, beta
      real(wp), intent(in) :: a(lda, *), x(*)
      real(wp), intent(inout) :: y(*)
    end subroutine dsymv

    ! Multiply a complex Hermitian matrix by a vector using one triangle.
    subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda, incx, incy
      complex(wp), intent(in) :: alpha, beta
      complex(wp), intent(in) :: a(lda, *), x(*)
      complex(wp), intent(inout) :: y(*)
    end subroutine zhemv

    ! General real matrix-vector product, used for Q and R applications.
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      real(wp), intent(in) :: alpha, beta
      real(wp), intent(in) :: a(lda, *), x(*)
      real(wp), intent(inout) :: y(*)
    end subroutine dgemv

    ! General complex matrix-vector product, including conjugate transpose.
    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      complex(wp), intent(in) :: alpha, beta
      complex(wp), intent(in) :: a(lda, *), x(*)
      complex(wp), intent(inout) :: y(*)
    end subroutine zgemv

    ! Solve an upper-triangular real system in place.
    subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      import wp
      character(len=1), intent(in) :: uplo, trans, diag
      integer, intent(in) :: n, lda, incx
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
    end subroutine dtrsv

    ! Solve an upper-triangular complex system in place.
    subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
      import wp
      character(len=1), intent(in) :: uplo, trans, diag
      integer, intent(in) :: n, lda, incx
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
    end subroutine ztrsv
  end interface

  type, public :: qr_real_state
#ifndef QRLINALG_TESTING
    private
#endif
    !The components are private in a normal library build. QRLINALG_TESTING is
    !defined only for the white-box test build, where direct access is needed
    !to verify factor and workspace invariants without adding copying accessors
    !to the public numerical interface.
    ! Current active order of the valid factors. Zero means no active factors.
    integer :: n = 0
    ! Largest order that fits in the allocated matrices and work arrays.
    integer :: capacity = 0
    ! Shift represented by Q*R = H - shift*S when valid is true.
    real(wp) :: shift = 0.0_wp
    ! True only when Q, R, n, and shift describe a successful factorization.
    logical :: valid = .false.
    ! Explicit orthogonal factor, stored with leading dimension capacity.
    real(wp), allocatable :: q(:,:)
    ! Explicit upper-triangular factor, also full leading-dimension storage.
    real(wp), allocatable :: r(:,:)
    ! Householder scalar factors produced temporarily by DGEQRF.
    real(wp), allocatable :: tau(:)
    ! Shared optimally sized workspace for DGEQRF and DORGQR.
    real(wp), allocatable :: factor_work(:)
    ! Scratch for destructive qrupdate arguments and residual-action vectors.
    real(wp), allocatable :: update_work(:)
    ! Two capacity-length vectors used by inverse iteration and quotients.
    real(wp), allocatable :: solve_work(:)
    ! Count of successful append, replacement, and deletion operations.
    integer(int64) :: structural_updates = 0_int64
    ! Successful structural updates since the last fresh factorization.
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => real_initialize
    procedure :: clear => clear_real_state
    procedure :: factorize_fresh => real_factorize_fresh
    procedure :: replace_symmetric => real_replace_symmetric
    procedure :: append_symmetric => real_append_symmetric
    procedure :: delete_symmetric => real_delete_symmetric
    procedure :: solve => real_solve
    procedure :: factorization_residual => real_factorization_residual
    procedure :: is_valid => real_is_valid
    procedure :: order => real_order
    procedure :: get_capacity => real_get_capacity
    procedure :: get_shift => real_get_shift
    procedure :: get_update_count => real_get_update_count
    procedure :: get_updates_since_fresh => real_get_updates_since_fresh
  end type qr_real_state

  type, public :: qr_complex_state
#ifndef QRLINALG_TESTING
    private
#endif
    !The component-access policy and metadata definitions are identical to
    !qr_real_state. Matrix and vector storage is complex; shift and the
    !rotation workspace used by complex qrupdate kernels remain real.
    integer :: n = 0
    integer :: capacity = 0
    real(wp) :: shift = 0.0_wp
    logical :: valid = .false.
    ! Explicit unitary Q and upper-triangular R.
    complex(wp), allocatable :: q(:,:)
    complex(wp), allocatable :: r(:,:)
    ! Householder scalars and optimally sized complex LAPACK workspace.
    complex(wp), allocatable :: tau(:)
    complex(wp), allocatable :: factor_work(:)
    ! Destructive qrupdate copies and residual-action vector scratch.
    complex(wp), allocatable :: update_work(:)
    ! Two-vector inverse-iteration and shifted-quotient workspace.
    complex(wp), allocatable :: solve_work(:)
    ! Real cosines and other real workspace required by complex qrupdate.
    real(wp), allocatable :: real_work(:)
    integer(int64) :: structural_updates = 0_int64
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => complex_initialize
    procedure :: clear => clear_complex_state
    procedure :: factorize_fresh => complex_factorize_fresh
    procedure :: replace_symmetric => complex_replace_symmetric
    procedure :: append_symmetric => complex_append_symmetric
    procedure :: delete_symmetric => complex_delete_symmetric
    procedure :: solve => complex_solve
    procedure :: factorization_residual => complex_factorization_residual
    procedure :: is_valid => complex_is_valid
    procedure :: order => complex_order
    procedure :: get_capacity => complex_get_capacity
    procedure :: get_shift => complex_get_shift
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
    case default
      message = 'unrecognized qrlinalg status code'
    end select
  end function qr_status_message

  subroutine clear_real_state(self)
  !Subroutine clear_real_state releases every allocation owned by a real QR
  !state and restores the state to its default, uninitialized condition. It is
  !the implementation of the public clear() binding and is also used before
  !reinitialization and after a partial allocation failure. Calling clear() for
  !an already empty state is valid and has no effect beyond restoring default
  !metadata values.
  !
  !  Input/output parameter:
  !    self - The real QR state. On exit all allocatable components are
  !           unallocated; n and capacity are zero; valid is false; shift and
  !           both structural-update counters are zero.
  !
  !The routine allocates no storage and has no failure status. After clear(),
  !initialize must be called before another factorization can be constructed.
    class(qr_real_state), intent(inout) :: self

    if (allocated(self%q)) deallocate(self%q)
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%tau)) deallocate(self%tau)
    if (allocated(self%factor_work)) deallocate(self%factor_work)
    if (allocated(self%update_work)) deallocate(self%update_work)
    if (allocated(self%solve_work)) deallocate(self%solve_work)
    self%n = 0
    self%capacity = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_real_state

  subroutine clear_complex_state(self)
  !Subroutine clear_complex_state releases every allocation owned by a complex
  !QR state and restores the state to its default, uninitialized condition. In
  !addition to the complex arrays it releases the real workspace required by
  !complex plane rotations. It implements the public clear() binding and is
  !also used by initialization failure handling. Calling clear() for an empty
  !state is valid.
  !
  !  Input/output parameter:
  !    self - The complex QR state. On exit all allocatable components are
  !           unallocated and all metadata and counters have their default
  !           values.
  !
  !The routine allocates no storage and has no failure status. After clear(),
  !initialize must be called before another factorization can be constructed.
    class(qr_complex_state), intent(inout) :: self

    if (allocated(self%q)) deallocate(self%q)
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%tau)) deallocate(self%tau)
    if (allocated(self%factor_work)) deallocate(self%factor_work)
    if (allocated(self%update_work)) deallocate(self%update_work)
    if (allocated(self%solve_work)) deallocate(self%solve_work)
    if (allocated(self%real_work)) deallocate(self%real_work)
    self%n = 0
    self%capacity = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_complex_state

  subroutine real_initialize(self, capacity, info)
  !Subroutine real_initialize prepares a real QR state for symmetric matrices
  !of order not greater than capacity. Any factorization and counters previously
  !held by the state are discarded.
  !
  !The following storage is allocated:
  !  - full capacity by capacity arrays for the explicit factors Q and R;
  !  - capacity Householder coefficients in tau;
  !  - 4*capacity elements for structural QR updates and residual actions;
  !  - 2*capacity elements for inverse iteration and Rayleigh quotients;
  !  - one factorization workspace large enough for both DGEQRF and DORGQR.
  !
  !The size of factor_work is obtained by querying DGEQRF and DORGQR with
  !LWORK=-1 for a square matrix of order capacity. The larger recommended size is
  !allocated once and reused. Initialization does not form a QR factorization;
  !therefore n is zero and valid is false on successful exit.
  !
  !  Input parameter:
  !    capacity - Maximum matrix order supported by the state. It must be
  !               positive. Initialization reserves storage but does not set
  !               the active matrix order.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Existing allocations and factorization
  !            metadata are destroyed before capacity is validated.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when all storage is ready;
  !            QR_ERR_INVALID_ARGUMENT when capacity is not positive;
  !            QR_ERR_ALLOCATION when an allocation fails;
  !            QR_ERR_FACTORIZATION when a LAPACK workspace query fails.
  !
  !If initialization fails, self is returned in the empty state described for
  !clear_real_state; no partial allocation remains owned by the object.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: capacity
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    real(wp) :: work_query(1)

    !Discard all previous storage before constructing the new state.
    call clear_real_state(self)
    if (capacity <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    !Allocate every array whose extent follows directly from capacity.
    !factor_work is allocated after the two workspace queries below.
    allocate(self%q(capacity, capacity), self%r(capacity, capacity), &
             self%tau(capacity), self%update_work(4 * capacity), &
             self%solve_work(2 * capacity), &
             stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = 0.0_wp
    self%r = 0.0_wp
    self%tau = 0.0_wp
    self%update_work = 0.0_wp
    self%solve_work = 0.0_wp

    !Query the workspace recommended for the compact Householder
    !factorization. LWORK=-1 performs no factorization and returns the
    !recommendation in WORK(1).
    call dgeqrf(capacity, capacity, self%q, capacity, self%tau, work_query, -1, &
                lapack_info)
    if (lapack_info /= 0) then
      call clear_real_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    factor_lwork = max(capacity, ceiling(work_query(1)))

    call dorgqr(capacity, capacity, capacity, self%q, capacity, self%tau, &
                work_query, -1, lapack_info)
    if (lapack_info /= 0) then
      call clear_real_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    generate_q_lwork = max(capacity, ceiling(work_query(1)))

    optimal_lwork = max(factor_lwork, generate_q_lwork)
    allocate(self%factor_work(optimal_lwork), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if
    self%factor_work = 0.0_wp

    !Only capacity is committed here. The active order, represented shift, and
    !validity flag are committed by a successful factorize_fresh operation.
    self%capacity = capacity
    info = QR_SUCCESS
  end subroutine real_initialize

  subroutine complex_initialize(self, capacity, info)
  !Subroutine complex_initialize prepares a complex QR state for Hermitian
  !matrices of order not greater than capacity. It follows the allocation and
  !failure semantics of real_initialize, with the following differences:
  !  - Q, R, tau, factor_work, update_work, and solve_work are complex(wp);
  !  - update_work also supplies three vector slices to the factorization
  !    residual operation;
  !  - real_work contains capacity real(wp) elements required by the complex
  !    qrupdate rotation kernels;
  !  - ZGEQRF and ZUNGQR supply the two factorization-work recommendations.
  !
  !  Input parameter:
  !    capacity - Maximum matrix order supported by the state; must be
  !               positive. It is a storage reservation, not the active order.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Its old factors, work arrays, metadata,
  !            and counters are discarded.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, QR_ERR_ALLOCATION, or
  !            QR_ERR_FACTORIZATION, with the meanings documented for
  !            real_initialize.
  !
  !On failure self is empty. On success self%capacity=capacity, n=0, and
  !valid=false.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: capacity
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    complex(wp) :: work_query(1)

    !Discard all previous complex and real state storage.
    call clear_complex_state(self)
    if (capacity <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    allocate(self%q(capacity, capacity), self%r(capacity, capacity), &
             self%tau(capacity), self%update_work(4 * capacity), &
             self%solve_work(2 * capacity), &
             self%real_work(capacity), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%r = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%tau = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%update_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%solve_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%real_work = 0.0_wp

    !In a complex LAPACK workspace query the recommended integer workspace
    !length is returned in the real part of WORK(1). Query both stages and use
    !the larger recommendation.
    call zgeqrf(capacity, capacity, self%q, capacity, self%tau, work_query, -1, &
                lapack_info)
    if (lapack_info /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    factor_lwork = max(capacity, ceiling(real(work_query(1), wp)))

    call zungqr(capacity, capacity, capacity, self%q, capacity, self%tau, &
                work_query, -1, lapack_info)
    if (lapack_info /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    generate_q_lwork = max(capacity, ceiling(real(work_query(1), wp)))

    optimal_lwork = max(factor_lwork, generate_q_lwork)
    allocate(self%factor_work(optimal_lwork), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if
    self%factor_work = cmplx(0.0_wp, 0.0_wp, kind=wp)

    self%capacity = capacity
    info = QR_SUCCESS
  end subroutine complex_initialize

  subroutine real_factorize_fresh(self, h, s, shift, info, active_order)
  !Subroutine real_factorize_fresh constructs a complete QR factorization of
  !the shifted real symmetric matrix
  !
  !                     M = H - shift*S = Q*R .
  !
  !The matrix M is formed directly in the state-owned Q array. DGEQRF replaces
  !M by the upper-triangular factor R and a compact set of Householder vectors.
  !The upper triangle is copied to the state-owned R array before DORGQR
  !replaces the compact representation by the explicit orthogonal matrix Q.
  !The strict lower triangle of R is explicitly set to zero.
  !
  !No temporary matrix of order n is created, and neither H nor S is modified
  !or retained. The state must have been initialized with capacity at least n.
  !When active_order is absent, H and S must be square and have equal order,
  !which preserves the original interface contract. When active_order is
  !present, H and S may have larger physical storage and need only have at
  !least n rows and n columns. In both forms, only the lower triangles of the
  !leading n by n principal blocks are referenced. The opposite triangle of M
  !is generated by symmetry, so upper-triangular and inactive input entries
  !need not be initialized.
  !
  !  Input parameters:
  !    h     - Real array whose leading n by n lower triangle contains the
  !            symmetric Hamiltonian matrix H. The upper triangle and entries
  !            outside the active principal block are ignored.
  !    s     - Real array whose leading n by n lower triangle contains the
  !            symmetric overlap matrix S. The upper triangle and entries
  !            outside the active principal block are ignored.
  !    shift - The real shift represented by the factorization.
  !    active_order - Optional active order n. It must be positive, must not
  !                   exceed state capacity, and must fit within both extents
  !                   of H and S. If absent, n is the common square order of
  !                   H and S.
  !
  !  Input/output parameter:
  !    self  - An initialized real QR state. On successful exit:
  !              self%n = n,
  !              self%shift = shift,
  !              self%valid = true,
  !              self%Q*self%R = H-shift*S to working precision,
  !              self%updates_since_fresh = 0.
  !            structural_updates is a lifetime counter and is not reset.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when both LAPACK stages succeed;
  !            QR_ERR_INVALID_STATE when self has not been initialized or its
  !            required state-owned factorization storage is unavailable;
  !            QR_ERR_INVALID_ARGUMENT when active_order is present but is
  !            not positive;
  !            QR_ERR_DIMENSION_MISMATCH when the original no-active-order
  !            matrix rules are violated, or when an explicit active block
  !            does not fit within H or S;
  !            QR_ERR_CAPACITY_EXCEEDED when n is greater than the initialized
  !            capacity;
  !            QR_ERR_FACTORIZATION when DGEQRF or DORGQR reports an error.
  !
  !All argument checks are completed before existing factors are overwritten.
  !An argument error therefore preserves the previous factorization. A LAPACK
  !error occurs after the state buffers have been modified; in that case n is
  !returned as zero and valid is false so that partial factors cannot be used.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer, intent(in), optional :: active_order
    integer :: i, j, lapack_info, matrix_n

    !Check initialization, matrix shapes, active order, and capacity before
    !overwriting any component of a previously valid factorization.
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
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

    !The Q buffer is about to be overwritten. Mark the factors invalid until
    !both the factorization and explicit-Q generation have completed.
    self%valid = .false.
    self%n = 0

    !Read each stored lower-triangular column contiguously and reflect every
    !off-diagonal value into the opposite triangle. LAPACK requires the full
    !shifted matrix, but callers need define only the canonical lower triangle.
    do j = 1, matrix_n
      do i = j, matrix_n
        self%q(i, j) = h(i, j) - shift * s(i, j)
        if (i /= j) self%q(j, i) = self%q(i, j)
      end do
    end do

    !Compute the compact Householder representation. On exit from DGEQRF the
    !upper triangle contains R; the strict lower triangle and tau describe Q.
    call dgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Copy R before DORGQR destroys the compact reflector representation. Clear
    !the strict lower triangle so all entries in the active R block have a
    !defined triangular meaning, including after refactorization at smaller n.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = 0.0_wp
      end do
    end do

    !Expand the Householder vectors and tau into the explicit orthogonal matrix
    !Q required by inverse iteration and the qrupdate kernels.
    call dorgqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Commit the active factorization only after every numerical stage succeeds.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine real_factorize_fresh

  subroutine complex_factorize_fresh(self, h, s, shift, info, active_order)
  !Subroutine complex_factorize_fresh constructs a complete QR factorization
  !of the shifted complex Hermitian matrix
  !
  !                     M = H - shift*S = Q*R ,
  !
  !where Q is unitary and R is upper triangular. ZGEQRF produces a compact
  !Householder representation, R is copied from its upper triangle, and ZUNGQR
  !generates explicit Q. The routine allocates no memory and does not modify or
  !retain H and S.
  !
  !Only the lower triangles of the active leading blocks of H and S, including
  !their diagonals, are read. active_order has the same optional leading-block
  !meaning and backward-compatible validation rules as in
  !real_factorize_fresh. The upper triangle of the shifted matrix is generated
  !by conjugating its lower triangle. Each active physical diagonal must be
  !real within
  !
  !       100*epsilon(1.0_wp)*max(1,maxval(abs(lower triangle))) .
  !
  !Roundoff-sized imaginary diagonal parts are accepted and discarded before
  !M is formed. H and S are checked separately so cancellation cannot make two
  !non-Hermitian inputs acceptable.
  !
  !  Input parameters:
  !    h     - Complex array whose leading n by n lower triangle contains the
  !            Hermitian Hamiltonian matrix H. All other entries are ignored.
  !    s     - Complex array whose leading n by n lower triangle contains the
  !            Hermitian overlap matrix S. All other entries are ignored.
  !    shift - The real shift represented by the factorization.
  !    active_order - Optional positive active order n. If absent, H and S
  !                   must retain the original equal-square shape. If present,
  !                   both arrays need only contain the leading n by n block.
  !
  !  Input/output parameter:
  !    self  - An initialized complex QR state with capacity at least n. The
  !            successful state transitions and counter rules are identical to
  !            real_factorize_fresh, with unitary Q in place of orthogonal Q.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS, QR_ERR_INVALID_STATE,
  !            QR_ERR_DIMENSION_MISMATCH, QR_ERR_CAPACITY_EXCEEDED,
  !            QR_ERR_INVALID_ARGUMENT, or QR_ERR_FACTORIZATION under the
  !            conditions documented for real_factorize_fresh.
  !
  !A diagonal that is not real within the tolerance above returns
  !QR_ERR_INVALID_ARGUMENT and preserves existing factors. A failure after
  !ZGEQRF begins leaves n=0 and valid=false.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer, intent(in), optional :: active_order
    integer :: i, j, lapack_info, matrix_n
    complex(wp) :: shifted_value
    real(wp) :: h_diagonal_tolerance, h_scale
    real(wp) :: s_diagonal_tolerance, s_scale

    !Validate all state and dimension requirements before modifying Q or R.
    info = QR_ERR_INVALID_STATE
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
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

    !Determine diagonal tolerances from the same lower triangles that define
    !the matrices. Reading the unused upper triangles here would defeat the
    !single-triangle ownership contract.
    h_scale = 1.0_wp
    s_scale = 1.0_wp
    do j = 1, matrix_n
      do i = j, matrix_n
        h_scale = max(h_scale, abs(h(i,j)))
        s_scale = max(s_scale, abs(s(i,j)))
      end do
    end do
    h_diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * h_scale
    s_diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * s_scale
    info = QR_ERR_INVALID_ARGUMENT
    do i = 1, matrix_n
      if (abs(aimag(h(i,i))) > h_diagonal_tolerance) return
      if (abs(aimag(s(i,i))) > s_diagonal_tolerance) return
    end do

    self%valid = .false.
    self%n = 0

    !Construct an exactly Hermitian LAPACK input from the canonical lower
    !triangle. Diagonal imaginary roundoff is removed explicitly; each strict
    !upper entry is the conjugate of its lower counterpart.
    do j = 1, matrix_n
      self%q(j,j) = cmplx(real(h(j,j), wp) - shift * real(s(j,j), wp), &
                          0.0_wp, kind=wp)
      do i = j + 1, matrix_n
        shifted_value = h(i,j) - &
                        cmplx(shift, 0.0_wp, kind=wp) * s(i,j)
        self%q(i,j) = shifted_value
        self%q(j,i) = conjg(shifted_value)
      end do
    end do

    !Compute R and the compact unitary Householder representation of Q.
    call zgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Copy the active upper triangle to R and define its lower triangle as zero
    !before the reflector storage is replaced by explicit Q.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      end do
    end do

    !Generate the explicit unitary matrix Q from the reflectors and tau.
    call zungqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Publish the active order, shift, and validity only after ZUNGQR succeeds.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine complex_factorize_fresh

  subroutine real_replace_symmetric(self, idx, delta_h, delta_s, info)
  !Subroutine real_replace_symmetric replaces one row and the corresponding
  !column of a real symmetric problem by updating the stored QR factors. If
  !delta_h and delta_s denote the changes in the physical H and S columns, the
  !change represented by the QR state is
  !
  !                 d = delta_h - self%shift*delta_s .
  !
  !A symmetric row-and-column replacement can be expressed as
  !
  !                 d*e_idx^T + e_idx*(d-d(idx)*e_idx)^T .
  !
  !The subtraction of d(idx)*e_idx prevents the diagonal change from being
  !applied twice. The operation applies these two rank-one terms sequentially
  !with qr1up. qr1up overwrites its u and v arguments, so both calls use
  !state-owned copies and leave delta_h and delta_s unchanged.
  !
  !  Input parameters:
  !    idx     - One-based index of the replaced row and column.
  !    delta_h - Change in H(:,idx), including the diagonal element.
  !    delta_s - Change in S(:,idx), including the diagonal element.
  !
  !  Input/output parameter:
  !    self    - A valid QR state of active order n. On success Q and R
  !              represent the updated shifted matrix, valid, n, capacity, and
  !              shift are preserved, and both structural-update counters are
  !              incremented by one.
  !
  !  Output parameter:
  !    info    - QR_SUCCESS when both rank-one updates are applied;
  !              QR_ERR_INVALID_STATE when self has no valid factorization or
  !              required state-owned factor or workspace storage is absent;
  !              QR_ERR_INVALID_ARGUMENT when idx is outside 1:n;
  !              QR_ERR_DIMENSION_MISMATCH when either change vector has
  !              length other than n.
  !
  !All validation precedes modification of Q or R, so every error status
  !preserves the complete state. The validated qr1up calls have no numerical
  !failure return. No allocation is performed.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info
    integer :: i, matrix_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    matrix_n = self%n
    if (size(self%update_work) < 4 * matrix_n) return
    if (idx < 1 .or. idx > matrix_n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(delta_h) /= matrix_n .or. size(delta_s) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if

    !First apply d*e_idx^T. The update vectors occupy the first two workspace
    !blocks and qr1up uses the remaining two blocks as rotation workspace.
    do i = 1, matrix_n
      self%update_work(i) = delta_h(i) - self%shift * delta_s(i)
      self%update_work(matrix_n + i) = 0.0_wp
    end do
    self%update_work(matrix_n + idx) = 1.0_wp
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:4 * matrix_n))

    !Then apply e_idx*(d-d(idx)*e_idx)^T. Reconstruct d because qr1up is
    !permitted to overwrite both vectors supplied to the first call.
    do i = 1, matrix_n
      self%update_work(i) = 0.0_wp
      self%update_work(matrix_n + i) = delta_h(i) - &
                                           self%shift * delta_s(i)
    end do
    self%update_work(idx) = 1.0_wp
    self%update_work(matrix_n + idx) = 0.0_wp
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:4 * matrix_n))

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_replace_symmetric

  subroutine complex_replace_symmetric(self, idx, delta_h, delta_s, info)
  !Subroutine complex_replace_symmetric is the Hermitian counterpart of
  !real_replace_symmetric. For d=delta_h-self%shift*delta_s, the represented
  !Hermitian change is
  !
  !                 d*e_idx^H + e_idx*(d-d(idx)*e_idx)^H .
  !
  !The diagonal of a Hermitian matrix is real. The represented diagonal change
  !is therefore rejected when its imaginary part exceeds
  !
  !       100*epsilon(1.0_wp)*max(1,maxval(abs(d))).
  !
  !An accepted roundoff-sized imaginary part is discarded before either
  !rank-one update. qr1up interprets its complex update as u*v^H, so the two
  !terms reconstruct the conjugate row without modifying caller arrays.
  !
  !  Input parameters:
  !    idx     - One-based index of the replaced row and column.
  !    delta_h - Change in H(:,idx), including its nominally real diagonal.
  !    delta_s - Change in S(:,idx), including its nominally real diagonal.
  !
  !  Input/output parameter:
  !    self    - A valid complex QR state. Successful state and counter changes
  !              are the same as for real_replace_symmetric.
  !
  !  Output parameter:
  !    info    - QR_SUCCESS when the Hermitian replacement is complete;
  !              QR_ERR_INVALID_STATE, QR_ERR_DIMENSION_MISMATCH, or
  !              QR_ERR_INVALID_ARGUMENT under the real-routine validation
  !              conditions. QR_ERR_INVALID_ARGUMENT also reports a
  !              represented diagonal change that is not real within the
  !              tolerance above.
  !
  !Every rejection occurs before Q or R is modified. The routine allocates no
  !memory and increments each counter once, rather than once per rank-one term.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info
    complex(wp) :: change
    real(wp) :: diagonal_tolerance, update_scale
    integer :: i, matrix_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    if (.not. allocated(self%real_work)) return
    matrix_n = self%n
    if (size(self%update_work) < 3 * matrix_n) return
    if (size(self%real_work) < matrix_n) return
    if (idx < 1 .or. idx > matrix_n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(delta_h) /= matrix_n .or. size(delta_s) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if

    !Determine the scale and validate the diagonal before either qr1up call.
    update_scale = 1.0_wp
    do i = 1, matrix_n
      change = delta_h(i) - &
               cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
      update_scale = max(update_scale, abs(change))
    end do
    change = delta_h(idx) - &
             cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(idx)
    diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * update_scale
    if (abs(aimag(change)) > diagonal_tolerance) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    !Apply d*e_idx^H. The first three complex workspace blocks contain u, v,
    !and the qr1up work vector; real_work stores the rotation cosines.
    do i = 1, matrix_n
      self%update_work(i) = delta_h(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
      self%update_work(matrix_n + i) = &
        cmplx(0.0_wp, 0.0_wp, kind=wp)
    end do
    self%update_work(idx) = &
      cmplx(real(self%update_work(idx), wp), 0.0_wp, kind=wp)
    self%update_work(matrix_n + idx) = &
      cmplx(1.0_wp, 0.0_wp, kind=wp)
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), &
               self%real_work(1:matrix_n))

    !Apply e_idx*(d-d(idx)*e_idx)^H. Re-form d after the destructive first
    !call and set its diagonal component to exact complex zero.
    do i = 1, matrix_n
      self%update_work(i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      self%update_work(matrix_n + i) = delta_h(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
    end do
    self%update_work(idx) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    self%update_work(matrix_n + idx) = &
      cmplx(0.0_wp, 0.0_wp, kind=wp)
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), &
               self%real_work(1:matrix_n))

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_replace_symmetric

  subroutine real_append_symmetric(self, h_column, s_column, info)
  !Subroutine real_append_symmetric increases the active real symmetric
  !problem from order n to n+1 without recomputing a fresh factorization.
  !h_column and s_column contain the complete new physical columns, including
  !their diagonal elements. The column of the represented shifted matrix is
  !
  !             shifted_column = h_column-self%shift*s_column.
  !
  !qrinc first inserts shifted_column(1:n) as column n+1 of the existing
  !n-row factorization. This produces factors for an n by n+1 matrix. qrinr
  !then inserts the symmetric row
  !
  !             [ shifted_column(1:n), shifted_column(n+1) ]
  !
  !at row n+1 and restores a square factorization of order n+1. qrinr modifies
  !its row-vector argument, so the complete row is formed in state-owned
  !workspace and both caller arrays remain unchanged.
  !
  !  Input parameters:
  !    h_column - New H column of length self%n+1. Elements 1:self%n
  !               determine the symmetric off-diagonal row by transposition;
  !               the final element is the new diagonal.
  !    s_column - Corresponding new S column with the same extent and storage
  !               convention.
  !
  !  Input/output parameter:
  !    self     - A valid QR state. On success n increases by one, Q and R
  !               represent the expanded shifted matrix, valid, capacity, and
  !               shift are preserved, and both structural-update counters
  !               increase by one.
  !
  !  Output parameter:
  !    info     - QR_SUCCESS when both insertions complete;
  !               QR_ERR_INVALID_STATE when self has no valid factors or
  !               required state-owned factor/work storage is absent or too
  !               small;
  !               QR_ERR_CAPACITY_EXCEEDED when active order has reached the
  !               initialized capacity;
  !               QR_ERR_DIMENSION_MISMATCH when either input length is not
  !               self%n+1.
  !
  !All recoverable failures are detected before qrinc modifies the factors and
  !therefore preserve the complete state. The validated qrinc and qrinr calls
  !have no numerical failure result. The operation allocates no memory.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info
    integer :: i, new_n, old_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    old_n = self%n
    if (old_n >= self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    new_n = old_n + 1
    if (size(h_column) /= new_n .or. size(s_column) /= new_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    info = QR_ERR_INVALID_STATE
    if (size(self%q, 1) < self%capacity .or. &
        size(self%q, 2) < self%capacity) return
    if (size(self%r, 1) < self%capacity .or. &
        size(self%r, 2) < self%capacity) return
    if (size(self%update_work) < 2 * new_n) return

    !Insert only the off-diagonal part of the new column. The diagonal belongs
    !to the bottom row, which does not exist until qrinr is called.
    do i = 1, old_n
      self%update_work(i) = h_column(i) - self%shift * s_column(i)
    end do
    call qrinc(old_n, old_n, old_n, self%q, self%capacity, self%r, &
               self%capacity, new_n, self%update_work(1:old_n), &
               self%update_work(new_n:2 * old_n))

    !Form the complete symmetric bottom row in workspace disjoint from the
    !rotation cosines written by qrinr.
    do i = 1, new_n
      self%update_work(i) = h_column(i) - self%shift * s_column(i)
    end do
    call qrinr(old_n, new_n, self%q, self%capacity, self%r, &
               self%capacity, new_n, self%update_work(1:new_n), &
               self%update_work(new_n + 1:2 * new_n))

    self%n = new_n
    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_append_symmetric

  subroutine complex_append_symmetric(self, h_column, s_column, info)
  !Subroutine complex_append_symmetric is the Hermitian counterpart of
  !real_append_symmetric. The supplied columns determine the new shifted
  !column and, by conjugation, the new bottom row. Because H and S are each
  !Hermitian, both supplied diagonal elements must be real. They are accepted
  !only when
  !
  !  abs(aimag(column(n+1))) <=
  !      100*epsilon(1.0_wp)*max(1,maxval(abs(column))).
  !
  !Accepted roundoff-sized imaginary diagonal parts are discarded separately
  !before H-shift*S is formed. qrinc inserts the off-diagonal shifted column;
  !qrinr inserts its conjugate-transposed row and the real shifted diagonal.
  !
  !  Input parameters:
  !    h_column - New complex H column of length self%n+1, including the real
  !               diagonal element.
  !    s_column - New complex S column of length self%n+1, including the real
  !               diagonal element.
  !
  !  Input/output parameter:
  !    self     - The complex QR state whose order is to be increased. State
  !               changes on success are the same as for the real routine.
  !
  !  Output parameter:
  !    info     - QR_SUCCESS when the Hermitian append completes;
  !               QR_ERR_INVALID_STATE, QR_ERR_CAPACITY_EXCEEDED, or
  !               QR_ERR_DIMENSION_MISMATCH under the corresponding validation
  !               conditions of real_append_symmetric;
  !               QR_ERR_INVALID_ARGUMENT when either physical diagonal has an
  !               imaginary part larger than its tolerance above.
  !
  !All rejection paths precede factor modification and preserve the complete
  !state. Caller arrays are not modified, and no allocation is performed.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info
    complex(wp) :: shifted_diagonal
    real(wp) :: h_diagonal_tolerance, h_scale
    real(wp) :: s_diagonal_tolerance, s_scale
    integer :: i, new_n, old_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    if (.not. allocated(self%real_work)) return
    old_n = self%n
    if (old_n >= self%capacity) then
      info = QR_ERR_CAPACITY_EXCEEDED
      return
    end if
    new_n = old_n + 1
    if (size(h_column) /= new_n .or. size(s_column) /= new_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    info = QR_ERR_INVALID_STATE
    if (size(self%q, 1) < self%capacity .or. &
        size(self%q, 2) < self%capacity) return
    if (size(self%r, 1) < self%capacity .or. &
        size(self%r, 2) < self%capacity) return
    if (size(self%update_work) < new_n) return
    if (size(self%real_work) < old_n) return

    !Hermitian diagonals are properties of H and S separately; cancellation
    !between their imaginary parts must not make two invalid inputs acceptable.
    h_scale = max(1.0_wp, maxval(abs(h_column)))
    s_scale = max(1.0_wp, maxval(abs(s_column)))
    h_diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * h_scale
    s_diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * s_scale
    if (abs(aimag(h_column(new_n))) > h_diagonal_tolerance .or. &
        abs(aimag(s_column(new_n))) > s_diagonal_tolerance) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    do i = 1, old_n
      self%update_work(i) = h_column(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * s_column(i)
    end do
    call qrinc(old_n, old_n, old_n, self%q, self%capacity, self%r, &
               self%capacity, new_n, self%update_work(1:old_n), &
               self%real_work(1:old_n))

    !The first n row elements are conjugates of the supplied shifted column.
    !Construct the diagonal from separately real H and S values so the expanded
    !matrix is exactly Hermitian rather than Hermitian only within tolerance.
    do i = 1, old_n
      self%update_work(i) = conjg(h_column(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * s_column(i))
    end do
    shifted_diagonal = cmplx(real(h_column(new_n), wp) - &
      self%shift * real(s_column(new_n), wp), 0.0_wp, kind=wp)
    self%update_work(new_n) = shifted_diagonal
    call qrinr(old_n, new_n, self%q, self%capacity, self%r, &
               self%capacity, new_n, self%update_work(1:new_n), &
               self%real_work(1:old_n))

    self%n = new_n
    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_append_symmetric

  subroutine real_delete_symmetric(self, idx, info)
  !Subroutine real_delete_symmetric removes row idx and column idx from an
  !active real symmetric factorization. If P deletes component idx, the new
  !represented matrix is the principal submatrix
  !
  !                       M_new = P^T*M*P .
  !
  !qrdec first removes column idx from the n by n factorization, leaving an
  !n by n-1 factorization. qrder then removes row idx and restores square
  !factors of order n-1. Deleting the only active row and column is rejected:
  !the library has no valid order-zero factorization, and append_symmetric
  !requires an existing valid factorization.
  !
  !  Input parameter:
  !    idx  - One-based row and column index in the active range 1:self%n.
  !
  !  Input/output parameter:
  !    self - A valid QR state of order at least two. On success n decreases
  !           by one, valid, capacity, and shift are preserved, Q and R
  !           represent the requested principal submatrix, and both update
  !           counters increase by one.
  !
  !  Output parameter:
  !    info - QR_SUCCESS when both deletions complete;
  !           QR_ERR_INVALID_STATE when self has no valid factorization or
  !           required state-owned factor/work storage is absent or too small;
  !           QR_ERR_INVALID_ARGUMENT when self has order one or idx is outside
  !           1:self%n.
  !
  !All recoverable failures precede qrdec and preserve the complete state.
  !The validated qrdec and qrder kernels have no numerical failure result.
  !After success the now-inactive trailing row and column of Q and R are
  !cleared, preventing stale factor data from being exposed by later capacity
  !growth or white-box inspection. No allocation is performed.
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info
    integer :: i, new_n, old_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    old_n = self%n
    if (old_n < 2 .or. idx < 1 .or. idx > old_n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(self%q, 1) < self%capacity .or. &
        size(self%q, 2) < self%capacity) return
    if (size(self%r, 1) < self%capacity .or. &
        size(self%r, 2) < self%capacity) return
    if (size(self%update_work) < 2 * old_n) return
    new_n = old_n - 1

    !Remove the selected column while retaining all n rows, then remove the
    !row with the same physical index from the rectangular factorization.
    call qrdec(old_n, old_n, old_n, self%q, self%capacity, self%r, &
               self%capacity, idx, self%update_work(1:old_n))
    call qrder(old_n, new_n, self%q, self%capacity, self%r, &
               self%capacity, idx, self%update_work(1:2 * old_n))

    do i = 1, self%capacity
      self%q(old_n,i) = 0.0_wp
      self%q(i,old_n) = 0.0_wp
      self%r(old_n,i) = 0.0_wp
      self%r(i,old_n) = 0.0_wp
    end do
    self%n = new_n
    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_delete_symmetric

  subroutine complex_delete_symmetric(self, idx, info)
  !Subroutine complex_delete_symmetric is the Hermitian counterpart of
  !real_delete_symmetric. Deletion introduces no numerical values and therefore
  !requires no diagonal-reality check or caller vector. qrdec and qrder preserve
  !the complex arithmetic and unitary Q of the retained principal submatrix.
  !
  !  Input parameter:
  !    idx  - One-based row and column index in the active range 1:self%n.
  !
  !  Input/output parameter:
  !    self - Complex QR state to update. Successful metadata, counter, and
  !           inactive-storage changes are the same as for the real routine.
  !
  !  Output parameter:
  !    info - QR_SUCCESS, QR_ERR_INVALID_STATE, or QR_ERR_INVALID_ARGUMENT
  !           under the validation conditions documented for
  !           real_delete_symmetric.
  !
  !All validation is completed before factor storage is modified. The routine
  !allocates no memory and rejects deletion from an order-one state.
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info
    integer :: i, new_n, old_n

    info = QR_ERR_INVALID_STATE
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    if (.not. allocated(self%real_work)) return
    old_n = self%n
    if (old_n < 2 .or. idx < 1 .or. idx > old_n) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (size(self%q, 1) < self%capacity .or. &
        size(self%q, 2) < self%capacity) return
    if (size(self%r, 1) < self%capacity .or. &
        size(self%r, 2) < self%capacity) return
    if (size(self%update_work) < old_n) return
    if (size(self%real_work) < old_n) return
    new_n = old_n - 1

    call qrdec(old_n, old_n, old_n, self%q, self%capacity, self%r, &
               self%capacity, idx, self%real_work(1:old_n))
    call qrder(old_n, new_n, self%q, self%capacity, self%r, &
               self%capacity, idx, self%update_work(1:old_n), &
               self%real_work(1:old_n))

    do i = 1, self%capacity
      self%q(old_n,i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      self%q(i,old_n) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      self%r(old_n,i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      self%r(i,old_n) = cmplx(0.0_wp, 0.0_wp, kind=wp)
    end do
    self%n = new_n
    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_delete_symmetric

  subroutine real_solve(self, s, v_initial, x, lambda, tol, max_iter, &
                        norm_mode, rel_acc, num_iter, info)
  !Subroutine real_solve finds one eigenvalue and its eigenvector for the real
  !generalized symmetric eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration. On entry self must contain valid factors
  !
  !                        H - shift*S = Q*R .
  !
  !The desired eigenvalue should be closer to shift than every other
  !eigenvalue. The convergence rate is governed principally by the ratio of
  !the distance from shift to the desired eigenvalue and the distance from
  !shift to the next closest eigenvalue. shift must not coincide with an
  !eigenvalue, because H-shift*S is then singular.
  !
  !A nondegenerate eigenvalue is assumed. If several eigenvalues are separated
  !only at the scale of working-precision roundoff, inverse iteration may
  !return a vector in their joint invariant subspace and may reach max_iter
  !without satisfying the requested directional accuracy. In that case the
  !best vector and Rayleigh quotient obtained are still returned with
  !QR_ERR_NO_CONVERGENCE.
  !
  !For a current vector v, one iteration solves
  !
  !             (H-shift*S)*x = S*v .
  !
  !Using the stored factors, the operation is performed as
  !
  !             w = S*v,
  !             y = Q^T*w,
  !             R*x = y.
  !
  !The new vector is divided by max(abs(x)). Its change of direction relative
  !to v is measured by
  !
  !             alpha   = (x^T*v)/(v^T*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0, iteration stops when rel_acc<=tol. For tol<=0, iteration
  !continues until rel_acc begins to increase and the current value is not
  !larger than abs(tol). The latter rule requests the most accurate attainable
  !direction subject to the floor abs(tol), and necessarily performs at least
  !one additional comparison iteration.
  !
  !After iteration, the eigenvalue is obtained from the shifted Rayleigh
  !quotient
  !
  ! lambda = shift + (x^T*(H-shift*S)*x)/(x^T*S*x).
  !
  !The state does not store H-shift*S separately. Its action on x is evaluated
  !as Q*(R*x). S may have physical extents larger than n; DSYMV uses its actual
  !first extent as the leading dimension and references only the lower triangle
  !of the leading n by n block. The upper triangle and inactive storage may
  !contain unrelated values.
  !
  !  Input parameters:
  !    s         - Real array with at least n rows and n columns. The lower
  !                triangle of its leading n by n block contains the symmetric
  !                overlap matrix and must be defined. S is expected to be
  !                positive definite.
  !    v_initial - A nonzero starting approximation of length n. A vector with
  !                a substantial component in the desired eigendirection
  !                generally converges faster. The array is not modified.
  !    tol       - Directional convergence tolerance. Its sign selects the
  !                stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^T*S*x = 1;
  !                  1: x^T*x = 1;
  !                other: max(abs(x)) = 1.
  !
  !  Input/output parameter:
  !    self      - A valid real QR state of active order n. Q, R, n, shift,
  !                validity, and structural-update counters are unchanged.
  !                solve_work is overwritten and remains internal scratch.
  !
  !  Output parameters:
  !    x         - The final eigenvector approximation in the normalization
  !                selected by norm_mode. It has length n.
  !    lambda    - The Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - The directional difference from the final iteration. This
  !                is a convergence estimate, not a rigorously bounded error.
  !    num_iter  - Number of inverse iterations performed.
  !    info      - QR_SUCCESS when the stopping criterion is satisfied;
  !                QR_ERR_INVALID_STATE when self has no valid factorization
  !                or required state-owned factor/work storage is absent;
  !                QR_ERR_DIMENSION_MISMATCH when S, v_initial, or x has an
  !                extent inconsistent with the active order;
  !                QR_ERR_INVALID_ARGUMENT when max_iter is not positive;
  !                QR_ERR_ZERO_INITIAL_VECTOR when v_initial is numerically
  !                zero in the working precision;
  !                QR_ERR_NONPOSITIVE_OVERLAP when the final x^T*S*x is not
  !                greater than the precision-scaled zero threshold;
  !                QR_ERR_SINGULAR when R cannot be used safely or an
  !                iteration produces a numerically zero vector;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. In this
  !                case x, lambda, rel_acc, and num_iter describe the best
  !                approximation reached.
  !
  !No allocation is performed. On an error detected before the first
  !iteration, x and lambda are zero, num_iter is zero, and rel_acc is huge.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: s(:,:)
    real(wp), intent(in) :: v_initial(:)
    real(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter, info
    integer :: matrix_n
    real(wp) :: coefficient, current_norm_squared, eigenvector_norm_squared
    real(wp) :: max_component, norm_of_diff, norm_of_diff_previous
    real(wp) :: overlap_norm_squared, shifted_numerator
    logical :: not_converged

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    info = QR_ERR_INVALID_STATE

    !Verify the factorization state, array dimensions, iteration limit, and
    !starting-vector norm before entering a BLAS routine.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v_initial) /= matrix_n .or. size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (max_iter <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (real_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    !DGEQRF can complete for a rank-deficient matrix. Test the diagonal of R
    !explicitly before DTRSV performs divisions during inverse iteration.
    if (real_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !Workspace layout during iteration:
    !  solve_work(1:n)       contains the current vector v;
    !  solve_work(n+1:2*n)   contains Q^T*S*v and triangular-solve scratch.
    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform inverse iterations until the direction criterion is satisfied or
    !the caller-supplied iteration limit is exhausted.
    do while (not_converged .and. num_iter < max_iter)
      !Form the right-hand side w=S*v from the lower triangle of S.
      call dsymv('L', matrix_n, 1.0_wp, s, size(s, 1), &
                 self%solve_work(1:matrix_n), 1, 0.0_wp, x, 1)
      !Transform w by Q^T. The result is the right-hand side of R*x=y.
      call dgemv('T', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
                 x, 1, 0.0_wp, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      !Solve the upper-triangular system in place to obtain the new iterate.
      call dtrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      !Scale the solution so its largest absolute component is one. Scaling
      !does not change the eigendirection and limits growth in repeated solves.
      max_component = maxval(abs(x))
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / max_component

      !Compute the least-direction-change coefficient and the relative norm of
      !the component of x not parallel to the previous iterate v.
      current_norm_squared = real_norm_squared( &
                               matrix_n, self%solve_work(1:matrix_n))
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = dot_product(x, self%solve_work(1:matrix_n)) / &
                    current_norm_squared
      eigenvector_norm_squared = real_norm_squared(matrix_n, x)
      norm_of_diff = real_direction_difference( &
                       matrix_n, x, coefficient, &
                       self%solve_work(1:matrix_n), &
                       eigenvector_norm_squared)

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        !For a non-positive tolerance, accept only after the direction error
        !has passed its minimum and remains within abs(tol).
        if (norm_of_diff > norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        self%solve_work(1:matrix_n) = x
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

    !Compute x^T*S*x for the Rayleigh quotient and, for norm_mode=0, final
    !normalization. A non-positive result violates the positive-definite S
    !precondition or indicates unusable numerical data.
    call dsymv('L', matrix_n, 1.0_wp, s, size(s, 1), x, 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = dot_product(x, &
                                       self%solve_work(1:matrix_n))
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if

    !Evaluate M*x as Q*(R*x), form x^T*M*x, and add the stored shift to the
    !quotient. This requires two matrix-vector products but no stored M matrix.
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%r, self%capacity, &
               x, 1, 0.0_wp, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1)
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    shifted_numerator = dot_product(x, self%solve_work(1:matrix_n))
    lambda = self%shift + shifted_numerator / overlap_norm_squared

    !Apply the normalization requested by the caller. No scaling is required
    !for other norm_mode values because each inverse iterate already has unit
    !largest-component magnitude.
    select case (norm_mode)
    case (0) !Normalize so that x^T*S*x=1.
      x = x / sqrt(overlap_norm_squared)
    case (1) !Normalize so that x^T*x=1.
      eigenvector_norm_squared = real_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / sqrt(eigenvector_norm_squared)
    end select
  end subroutine real_solve

  subroutine complex_solve(self, s, v_initial, x, lambda, tol, max_iter, &
                           norm_mode, rel_acc, num_iter, info)
  !Subroutine complex_solve finds one eigenvalue and its eigenvector for the
  !complex generalized Hermitian eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration. On entry self must contain valid factors
  !
  !                        H - shift*S = Q*R ,
  !
  !where Q is unitary and R is upper triangular. The desired eigenvalue must be
  !closer to shift than the remaining eigenvalues for ordinary nondegenerate
  !inverse iteration to select it. shift must not equal an eigenvalue.
  !Pathologically close or degenerate eigenvalues may yield a vector in their
  !joint invariant subspace and QR_ERR_NO_CONVERGENCE rather than a uniquely
  !determined eigenvector.
  !
  !For a current vector v, the iteration solves
  !
  !             (H-shift*S)*x = S*v
  !
  !through the three operations
  !
  !             w = S*v,
  !             y = Q^H*w,
  !             R*x = y.
  !
  !The vector is divided by the largest magnitude among all of its real and
  !imaginary components. This is not the same as division by max(abs(x)); it
  !preserves the component-scaling convention of the Hermitian inverse
  !iteration interface. Directional convergence is estimated by
  !
  !             alpha   = (v^H*x)/(v^H*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0 the first rel_acc<=tol terminates iteration. For tol<=0, the
  !routine waits until rel_acc begins to increase and is no greater than
  !abs(tol), thereby seeking the smallest attainable direction change subject
  !to the requested floor.
  !
  !The final real eigenvalue approximation is the shifted Hermitian Rayleigh
  !quotient
  !
  ! lambda = shift + real(x^H*(H-shift*S)*x)/(x^H*S*x).
  !
  !The product (H-shift*S)*x is evaluated as Q*(R*x). S may have physical
  !extents larger than n; ZHEMV uses size(S,1) as its leading dimension and
  !references only the lower triangle of the leading n by n block. Diagonal
  !elements of that block are assumed real and S is expected to be positive
  !definite. Upper-triangular and inactive entries are ignored.
  !
  !  Input parameters:
  !    s         - Complex array with at least n rows and n columns. The lower
  !                triangle and real diagonal of its leading n by n Hermitian
  !                block must be defined.
  !    v_initial - A nonzero complex starting approximation of length n. It is
  !                not modified.
  !    tol       - Real directional convergence tolerance. Its sign selects
  !                the stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^H*S*x = 1;
  !                  1: x^H*x = 1;
  !                other: the largest magnitude among every real and imaginary
  !                       component of x is one.
  !
  !  Input/output parameter:
  !    self      - A valid complex QR state of active order n. Numerical
  !                factors and public metadata are unchanged; solve_work is
  !                overwritten as private scratch storage.
  !
  !  Output parameters:
  !    x         - Final complex eigenvector approximation in the requested
  !                normalization.
  !    lambda    - Real Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - Real direction-change estimate from the final iteration.
  !    num_iter  - Number of inverse iterations performed.
  !    info      - QR_SUCCESS when convergence is detected;
  !                QR_ERR_INVALID_STATE when self has no valid factorization
  !                or required state-owned factor/work storage is absent;
  !                QR_ERR_DIMENSION_MISMATCH when S, v_initial, or x has an
  !                extent inconsistent with the active order;
  !                QR_ERR_INVALID_ARGUMENT when max_iter is not positive;
  !                QR_ERR_ZERO_INITIAL_VECTOR when v_initial is numerically
  !                zero in the working precision;
  !                QR_ERR_NONPOSITIVE_OVERLAP when the final x^H*S*x is not
  !                greater than the precision-scaled zero threshold;
  !                QR_ERR_SINGULAR for an unusable R or zero iterate;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. The latter
  !                status still returns the final eigenpair approximation.
  !
  !No allocation is performed. Before-iteration errors return zero x and
  !lambda, zero num_iter, and huge rel_acc.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: s(:,:)
    complex(wp), intent(in) :: v_initial(:)
    complex(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter, info
    integer :: matrix_n
    real(wp) :: current_norm_squared, eigenvector_norm_squared
    real(wp) :: max_component, norm_of_diff, norm_of_diff_previous
    real(wp) :: overlap_norm_squared, shifted_numerator
    complex(wp) :: coefficient
    complex(wp), parameter :: complex_one = (1.0_wp, 0.0_wp)
    complex(wp), parameter :: complex_zero = (0.0_wp, 0.0_wp)
    logical :: not_converged

    x = cmplx(0.0_wp, 0.0_wp, kind=wp)
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    info = QR_ERR_INVALID_STATE

    !Verify all state, dimension, iteration-limit, and starting-vector
    !requirements before using the stored factors or calling BLAS.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v_initial) /= matrix_n .or. size(x) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (max_iter <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
    if (complex_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    !ZGEQRF does not use a positive INFO value to report rank deficiency. The
    !diagonal of R is therefore tested explicitly before ZTRSV is called.
    if (complex_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !Workspace layout during iteration:
    !  solve_work(1:n)       contains the current vector v;
    !  solve_work(n+1:2*n)   contains Q^H*S*v and solve scratch.
    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform Hermitian inverse iterations until convergence or max_iter.
    do while (not_converged .and. num_iter < max_iter)
      !Form w=S*v from the stored lower triangle of the Hermitian matrix S.
      call zhemv('L', matrix_n, complex_one, s, size(s, 1), &
                 self%solve_work(1:matrix_n), 1, complex_zero, x, 1)
      !Apply Q^H to w to obtain the right-hand side of R*x=y.
      call zgemv('C', matrix_n, matrix_n, complex_one, self%q, &
                 self%capacity, x, 1, complex_zero, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      !Solve the complex upper-triangular system for the new iterate.
      call ztrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      !Scale with the largest real or imaginary component. This bounds both
      !parts of every element without changing the complex eigendirection.
      max_component = complex_max_abs_real_or_imag(matrix_n, x)
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(max_component, 0.0_wp, kind=wp)

      !Compute alpha=(v^H*x)/(v^H*v) and the prescribed relative difference
      !between x and alpha*v. Conjugating the previous iterate, rather than
      !the new iterate, makes this comparison invariant under the arbitrary
      !complex phase of an eigenvector.
      current_norm_squared = complex_norm_squared( &
                               matrix_n, self%solve_work(1:matrix_n))
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = dot_product(self%solve_work(1:matrix_n), x) / &
                    cmplx(current_norm_squared, 0.0_wp, kind=wp)
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      norm_of_diff = complex_direction_difference( &
                       matrix_n, x, coefficient, &
                       self%solve_work(1:matrix_n), &
                       eigenvector_norm_squared)

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        !A non-positive tolerance accepts a result only after the direction
        !error turns upward while remaining within abs(tol).
        if (norm_of_diff > norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        self%solve_work(1:matrix_n) = x
      end if
    end do

    rel_acc = norm_of_diff
    !A nonconverged final iterate remains a valid approximation and is carried
    !through the Rayleigh-quotient and normalization calculations below.
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    !Compute the Hermitian quadratic form x^H*S*x. Its real value is used both
    !in the Rayleigh quotient and, for norm_mode=0, in final normalization.
    call zhemv('L', matrix_n, complex_one, s, size(s, 1), x, 1, &
               complex_zero, self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = real( &
                               dot_product(x, &
                                 self%solve_work(1:matrix_n)), wp)
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_NONPOSITIVE_OVERLAP
      return
    end if

    !Evaluate M*x=Q*(R*x), form real(x^H*M*x), and add the stored shift. The
    !imaginary roundoff part of the Hermitian quadratic form is discarded.
    call zgemv('N', matrix_n, matrix_n, complex_one, self%r, &
               self%capacity, x, 1, complex_zero, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1)
    call zgemv('N', matrix_n, matrix_n, complex_one, self%q, &
               self%capacity, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1, &
               complex_zero, self%solve_work(1:matrix_n), 1)
    shifted_numerator = real( &
                            dot_product(x, &
                              self%solve_work(1:matrix_n)), wp)
    lambda = self%shift + shifted_numerator / overlap_norm_squared

    !Apply S or Euclidean normalization when requested. For every other mode,
    !retain the real/imaginary component scaling established in the iteration.
    select case (norm_mode)
    case (0) !Normalize so that x^H*S*x=1.
      x = x / cmplx(sqrt(overlap_norm_squared), 0.0_wp, kind=wp)
    case (1) !Normalize so that x^H*x=1.
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(sqrt(eigenvector_norm_squared), 0.0_wp, kind=wp)
    end select
  end subroutine complex_solve

  subroutine real_factorization_residual(self, h, s, v, &
                                         absolute_residual, &
                                         relative_residual, info)
  !Subroutine real_factorization_residual measures the action error between a
  !caller-owned real symmetric shifted matrix and the QR factors stored by the
  !state. For the active order n, represented shift sigma, and probe vector v,
  !it forms
  !
  !       physical_action = (H-sigma*S)*v,
  !       factor_action   = Q*(R*v),
  !       difference      = physical_action-factor_action.
  !
  !The returned residuals are
  !
  !       absolute_residual = ||difference||_2,
  !
  !       relative_residual = ||difference||_2 /
  !         (||physical_action||_2+||factor_action||_2+tiny(1.0_wp)).
  !
  !The positive safe-zero term makes an exactly zero comparison return a zero
  !relative residual rather than an undefined 0/0. The value is an action
  !residual for the caller-selected v, not a complete matrix norm or a bound on
  !factorization error for every direction. The library does not interpret the
  !value or trigger fresh factorization; refresh policy belongs to the caller.
  !
  !  Input parameters:
  !    h - Real array with at least n rows and n columns. The lower triangle of
  !        its leading n by n block contains H. Its physical first extent is
  !        used as the BLAS leading dimension. The upper triangle and inactive
  !        rows and columns are ignored.
  !    s - Real array with at least n rows and n columns, with the corresponding
  !        lower triangle of S. It may have a leading dimension different from
  !        H and is not required to have the same inactive extents.
  !    v - Nonzero real probe vector with exactly n elements. It is not
  !        modified.
  !
  !  Input/output parameter:
  !    self - Valid real QR state. Three n-element slices of update_work are
  !           overwritten as scratch. Q, R, order, capacity, shift, validity,
  !           and both structural-update counters are unchanged.
  !
  !  Output parameters:
  !    absolute_residual - Euclidean norm of difference. It is initialized to
  !                        zero before validation and remains zero on error.
  !    relative_residual - Symmetrically scaled action error defined above. It
  !                        is initialized to zero and remains zero on error.
  !    info              - QR_SUCCESS when both actions are evaluated;
  !                        QR_ERR_INVALID_STATE when no valid factors or
  !                        sufficient owned workspace are available;
  !                        QR_ERR_DIMENSION_MISMATCH when H or S cannot contain
  !                        the active block or V does not have length n;
  !                        QR_ERR_ZERO_INITIAL_VECTOR when V is numerically
  !                        zero according to the library norm threshold.
  !
  !The routine performs no allocation and no input/output. It reads only the
  !active lower triangles of H and S and does not modify caller arrays.
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in), contiguous :: h(:,:), s(:,:)
    real(wp), intent(in) :: v(:)
    real(wp), intent(out) :: absolute_residual, relative_residual
    integer, intent(out) :: info
    integer :: i, matrix_n
    real(wp) :: factor_norm, physical_norm

    absolute_residual = 0.0_wp
    relative_residual = 0.0_wp
    info = QR_ERR_INVALID_STATE

    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    matrix_n = self%n
    if (size(self%update_work) < 3 * matrix_n) return
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (real_norm_squared(matrix_n, v) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    !Form the physical shifted action from the two canonical lower triangles.
    !Using the descriptor extents as LDA preserves padded leading dimensions.
    call dsymv('L', matrix_n, 1.0_wp, h, size(h, 1), v, 1, 0.0_wp, &
               self%update_work(1:matrix_n), 1)
    call dsymv('L', matrix_n, -self%shift, s, size(s, 1), v, 1, 1.0_wp, &
               self%update_work(1:matrix_n), 1)

    !The second slice holds R*v and the third holds Q*(R*v), leaving the
    !physical action intact until its norm has been evaluated.
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%r, self%capacity, &
               v, 1, 0.0_wp, &
               self%update_work(matrix_n + 1:2 * matrix_n), 1)
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
               self%update_work(matrix_n + 1:2 * matrix_n), 1, 0.0_wp, &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), 1)

    physical_norm = sqrt(real_norm_squared( &
      matrix_n, self%update_work(1:matrix_n)))
    factor_norm = sqrt(real_norm_squared( &
      matrix_n, self%update_work(2 * matrix_n + 1:3 * matrix_n)))
    do i = 1, matrix_n
      self%update_work(i) = self%update_work(i) - &
                            self%update_work(2 * matrix_n + i)
    end do
    absolute_residual = sqrt(real_norm_squared( &
      matrix_n, self%update_work(1:matrix_n)))
    relative_residual = absolute_residual / &
      (physical_norm + factor_norm + tiny(1.0_wp))
    info = QR_SUCCESS
  end subroutine real_factorization_residual

  subroutine complex_factorization_residual(self, h, s, v, &
                                            absolute_residual, &
                                            relative_residual, info)
  !Subroutine complex_factorization_residual measures the action error between
  !a caller-owned Hermitian shifted matrix and the unitary QR factors stored by
  !a complex state. It evaluates the same physical_action, factor_action,
  !absolute-residual, and symmetrically scaled relative-residual equations as
  !real_factorization_residual, with complex matrix-vector products and
  !Euclidean norms induced by conjugate inner products.
  !
  !The lower triangles of the leading active blocks define H and S. ZHEMV uses
  !their separate physical first extents as leading dimensions, reconstructs
  !the conjugate upper action implicitly, and treats Hermitian diagonal values
  !as real. Upper-triangular entries and storage outside the active principal
  !blocks are ignored. The routine is a measurement only and never performs or
  !requests automatic refactorization.
  !
  !  Input parameters:
  !    h - Complex array with at least n rows and n columns whose leading lower
  !        triangle contains the Hermitian H matrix.
  !    s - Complex array with at least n rows and n columns whose leading lower
  !        triangle contains the Hermitian S matrix. H and S may have different
  !        physical leading dimensions and inactive extents.
  !    v - Nonzero complex probe vector with exactly n elements. It is not
  !        modified.
  !
  !  Input/output parameter:
  !    self - Valid complex QR state. Three n-element slices of update_work are
  !           used as scratch. Q, R, all public metadata, and update counters
  !           remain unchanged.
  !
  !  Output parameters:
  !    absolute_residual - Euclidean norm of the complex action difference;
  !                        initialized to zero on every entry path.
  !    relative_residual - Absolute residual divided by the sum of both action
  !                        norms and tiny(1.0_wp); initialized to zero.
  !    info              - QR_SUCCESS, QR_ERR_INVALID_STATE,
  !                        QR_ERR_DIMENSION_MISMATCH, or
  !                        QR_ERR_ZERO_INITIAL_VECTOR under the corresponding
  !                        conditions documented for the real routine.
  !
  !No allocation or input/output occurs. H, S, v, Q, R, order, shift,
  !validity, and both counters are unchanged on successful and failed calls.
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)
    complex(wp), intent(in) :: v(:)
    real(wp), intent(out) :: absolute_residual, relative_residual
    integer, intent(out) :: info
    integer :: i, matrix_n
    real(wp) :: factor_norm, physical_norm
    complex(wp), parameter :: complex_one = (1.0_wp, 0.0_wp)
    complex(wp), parameter :: complex_zero = (0.0_wp, 0.0_wp)

    absolute_residual = 0.0_wp
    relative_residual = 0.0_wp
    info = QR_ERR_INVALID_STATE

    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    matrix_n = self%n
    if (size(self%update_work) < 3 * matrix_n) return
    if (size(h, 1) < matrix_n .or. size(h, 2) < matrix_n .or. &
        size(s, 1) < matrix_n .or. size(s, 2) < matrix_n .or. &
        size(v) /= matrix_n) then
      info = QR_ERR_DIMENSION_MISMATCH
      return
    end if
    if (complex_norm_squared(matrix_n, v) <= tiny(1.0_wp)) then
      info = QR_ERR_ZERO_INITIAL_VECTOR
      return
    end if

    call zhemv('L', matrix_n, complex_one, h, size(h, 1), v, 1, &
               complex_zero, self%update_work(1:matrix_n), 1)
    call zhemv('L', matrix_n, &
               cmplx(-self%shift, 0.0_wp, kind=wp), &
               s, size(s, 1), v, 1, complex_one, &
               self%update_work(1:matrix_n), 1)

    call zgemv('N', matrix_n, matrix_n, complex_one, self%r, &
               self%capacity, v, 1, complex_zero, &
               self%update_work(matrix_n + 1:2 * matrix_n), 1)
    call zgemv('N', matrix_n, matrix_n, complex_one, self%q, &
               self%capacity, &
               self%update_work(matrix_n + 1:2 * matrix_n), 1, &
               complex_zero, &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), 1)

    physical_norm = sqrt(complex_norm_squared( &
      matrix_n, self%update_work(1:matrix_n)))
    factor_norm = sqrt(complex_norm_squared( &
      matrix_n, self%update_work(2 * matrix_n + 1:3 * matrix_n)))
    do i = 1, matrix_n
      self%update_work(i) = self%update_work(i) - &
                            self%update_work(2 * matrix_n + i)
    end do
    absolute_residual = sqrt(complex_norm_squared( &
      matrix_n, self%update_work(1:matrix_n)))
    relative_residual = absolute_residual / &
      (physical_norm + factor_norm + tiny(1.0_wp))
    info = QR_SUCCESS
  end subroutine complex_factorization_residual

  pure function real_is_valid(self) result(factors_are_valid)
  !Function real_is_valid reports whether a real state contains a complete QR
  !factorization that may be used by solve and the structural-update methods.
  !Initialization alone reserves storage but does not make the state valid.
  !The result becomes true only after factorize_fresh succeeds and remains true
  !across successful or rejected structural updates. A factorization failure
  !after numerical storage has been overwritten makes the result false.
  !
  !  Input parameter:
  !    self - Real QR state whose factorization validity is queried. The state
  !           and its internal workspace are not modified.
  !
  !  Result:
  !    factors_are_valid - True when the active Q and R factors are complete
  !                        and usable as a represented state; false otherwise.
    class(qr_real_state), intent(in) :: self
    logical :: factors_are_valid

    factors_are_valid = self%valid
  end function real_is_valid

  pure function complex_is_valid(self) result(factors_are_valid)
  !Function complex_is_valid reports whether a complex state contains a
  !complete QR factorization that may be used by solve and the Hermitian
  !structural-update methods. Its lifecycle semantics are identical to
  !real_is_valid. The query performs no allocation and does not modify self.
  !
  !  Input parameter:
  !    self - Complex QR state whose factorization validity is queried.
  !
  !  Result:
  !    factors_are_valid - True only when the active unitary Q and triangular R
  !                        constitute a valid represented factorization.
    class(qr_complex_state), intent(in) :: self
    logical :: factors_are_valid

    factors_are_valid = self%valid
  end function complex_is_valid

  pure function real_order(self) result(active_order)
  !Function real_order returns the active matrix order represented by a real
  !QR state. The result is zero before the first successful factorization,
  !after reinitialization, and after a numerical failure invalidates factors.
  !Successful append and deletion operations respectively increase and
  !decrease the result; replacement leaves it unchanged.
  !
  !  Input parameter:
  !    self - Real QR state whose active order is queried. It is not modified.
  !
  !  Result:
  !    active_order - Order of the represented factors, or zero when no valid
  !                   active factorization exists.
    class(qr_real_state), intent(in) :: self
    integer :: active_order

    active_order = self%n
  end function real_order

  pure function complex_order(self) result(active_order)
  !Function complex_order returns the active matrix order represented by a
  !complex QR state. The result follows the same initialization,
  !factorization, append, replacement, deletion, and invalidation rules as
  !real_order. The query performs no allocation and does not modify self.
  !
  !  Input parameter:
  !    self - Complex QR state whose active order is queried.
  !
  !  Result:
  !    active_order - Order of the represented factors, or zero when no valid
  !                   active factorization exists.
    class(qr_complex_state), intent(in) :: self
    integer :: active_order

    active_order = self%n
  end function complex_order

  pure function real_get_capacity(self) result(allocated_capacity)
  !Function real_get_capacity returns the largest matrix order supported by
  !the storage currently owned by a real QR state. Capacity is independent of
  !the active order and remains available after factorization and structural
  !updates. The result is zero for a default state and after failed or invalid
  !initialization has returned the state to its empty condition.
  !
  !  Input parameter:
  !    self - Real QR state whose allocated capacity is queried.
  !
  !  Result:
  !    allocated_capacity - Maximum order accepted without reinitialization,
  !                         or zero when the state owns no operating storage.
    class(qr_real_state), intent(in) :: self
    integer :: allocated_capacity

    allocated_capacity = self%capacity
  end function real_get_capacity

  pure function complex_get_capacity(self) result(allocated_capacity)
  !Function complex_get_capacity returns the largest matrix order supported by
  !the storage currently owned by a complex QR state. Its meaning and lifecycle
  !are identical to real_get_capacity. The query performs no allocation and
  !does not modify any complex or real workspace owned by self.
  !
  !  Input parameter:
  !    self - Complex QR state whose allocated capacity is queried.
  !
  !  Result:
  !    allocated_capacity - Maximum order accepted without reinitialization,
  !                         or zero when the state owns no operating storage.
    class(qr_complex_state), intent(in) :: self
    integer :: allocated_capacity

    allocated_capacity = self%capacity
  end function complex_get_capacity

  pure function real_get_shift(self) result(represented_shift)
  !Function real_get_shift returns the shift stored by a real QR state. When
  !is_valid() is true, the returned value is the sigma in
  !
  !                         H - sigma*S = Q*R .
  !
  !Before a successful factorization the stored value is zero. If is_valid()
  !is false after a numerical factorization failure, no shift is represented
  !and the returned stored value must not be used as factorization metadata.
  !A caller requiring an authoritative shift must first test is_valid().
  !
  !  Input parameter:
  !    self - Real QR state whose stored shift is queried. It is not modified.
  !
  !  Result:
  !    represented_shift - Stored real shift in the compile-time working kind.
    class(qr_real_state), intent(in) :: self
    real(wp) :: represented_shift

    represented_shift = self%shift
  end function real_get_shift

  pure function complex_get_shift(self) result(represented_shift)
  !Function complex_get_shift returns the real shift stored by a complex QR
  !state. When is_valid() is true, Q and R represent H-sigma*S at this shift,
  !where the matrices and factors are complex but the generalized Hermitian
  !eigenvalue and shift are real. As for real_get_shift, the result is not
  !factorization metadata while is_valid() is false.
  !
  !  Input parameter:
  !    self - Complex QR state whose stored real shift is queried.
  !
  !  Result:
  !    represented_shift - Stored real shift in the compile-time working kind.
    class(qr_complex_state), intent(in) :: self
    real(wp) :: represented_shift

    represented_shift = self%shift
  end function complex_get_shift

  pure function real_get_update_count(self) result(update_count)
  !Function real_get_update_count returns the lifetime number of successful
  !replace_symmetric, append_symmetric, and delete_symmetric calls recorded by
  !a real state. Each complete public structural operation contributes one,
  !irrespective of the number of internal QR rank-one operations. Rejected
  !operations and fresh factorizations do not increase the result.
  !Reinitialization starts a new state lifetime and resets the count to zero.
  !
  !  Input parameter:
  !    self - Real QR state whose lifetime update count is queried.
  !
  !  Result:
  !    update_count - Successful structural operations in the current state
  !                   lifetime, returned as a 64-bit integer.
    class(qr_real_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%structural_updates
  end function real_get_update_count

  pure function complex_get_update_count(self) result(update_count)
  !Function complex_get_update_count returns the lifetime number of successful
  !Hermitian replacement, append, and deletion operations recorded by a
  !complex state. Counting and reset semantics are identical to
  !real_get_update_count. The query does not modify factors or workspace.
  !
  !  Input parameter:
  !    self - Complex QR state whose lifetime update count is queried.
  !
  !  Result:
  !    update_count - Successful structural operations in the current state
  !                   lifetime, returned as a 64-bit integer.
    class(qr_complex_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%structural_updates
  end function complex_get_update_count

  pure function real_get_updates_since_fresh(self) result(update_count)
  !Function real_get_updates_since_fresh returns the number of successful
  !structural operations applied after the most recent successful fresh QR
  !factorization. A successful factorize_fresh resets this count to zero while
  !preserving the lifetime count returned by get_update_count(). Applications
  !may use this value as one input to an external refactorization policy; the
  !library does not impose a numerical-drift threshold itself.
  !
  !  Input parameter:
  !    self - Real QR state whose refresh counter is queried.
  !
  !  Result:
  !    update_count - Successful structural updates since fresh factorization,
  !                   returned as a 64-bit integer.
    class(qr_real_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%updates_since_fresh
  end function real_get_updates_since_fresh

  pure function complex_get_updates_since_fresh(self) result(update_count)
  !Function complex_get_updates_since_fresh returns the number of successful
  !Hermitian structural operations applied after the most recent successful
  !fresh QR factorization. The reset, rejection, and application-policy
  !semantics are identical to real_get_updates_since_fresh. The query performs
  !no allocation and does not modify self.
  !
  !  Input parameter:
  !    self - Complex QR state whose refresh counter is queried.
  !
  !  Result:
  !    update_count - Successful structural updates since fresh factorization,
  !                   returned as a 64-bit integer.
    class(qr_complex_state), intent(in) :: self
    integer(int64) :: update_count

    update_count = self%updates_since_fresh
  end function complex_get_updates_since_fresh

  function real_norm_squared(n, x) result(norm_squared)
  !Function real_norm_squared computes the real Euclidean inner product
  !
  !                         x^T*x = sum(x(i)^2)
  !
  !over the first n elements of x. The result is the squared norm; the square
  !root is not taken. The reduction is performed in wp and does not allocate an
  !array temporary.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Real vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - x^T*x for x(1:n).
    integer, intent(in) :: n
    real(wp), intent(in) :: x(:)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + x(i) * x(i)
    end do
  end function real_norm_squared

  function complex_norm_squared(n, x) result(norm_squared)
  !Function complex_norm_squared computes the Hermitian Euclidean inner
  !product
  !
  !              x^H*x = sum(real(x(i))^2+imag(x(i))^2)
  !
  !over the first n elements of a complex vector. The mathematically real
  !quantity is accumulated directly in real(wp), avoiding a complex reduction
  !and discarding no computed imaginary part.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Complex vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - The real value x^H*x for x(1:n).
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(:)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + real(x(i), wp)**2 + aimag(x(i))**2
    end do
  end function complex_norm_squared

  function real_direction_difference(n, x, alpha, y, x_norm_squared) &
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
    real(wp), intent(in) :: x(:), alpha, y(:), x_norm_squared
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

  function complex_direction_difference(n, x, alpha, y, x_norm_squared) &
      result(relative_difference)
  !Function complex_direction_difference computes the convergence measure
  !used to compare two successive complex inverse iterates,
  !
  !                 ||x-alpha*y||_2 / ||x||_2 .
  !
  !For Hermitian inverse iteration alpha=(y^H*x)/(y^H*y). Each squared
  !magnitude abs(x(i)-alpha*y(i))^2 is accumulated in real(wp).
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
    complex(wp), intent(in) :: x(:), alpha, y(:)
    real(wp), intent(in) :: x_norm_squared
    real(wp) :: relative_difference, difference_norm_squared
    complex(wp) :: difference
    integer :: i

    difference_norm_squared = 0.0_wp
    do i = 1, n
      difference = x(i) - alpha * y(i)
      difference_norm_squared = difference_norm_squared + abs(difference)**2
    end do
    relative_difference = sqrt(difference_norm_squared / x_norm_squared)
  end function complex_direction_difference

  function complex_max_abs_real_or_imag(n, x) result(max_component)
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
    complex(wp), intent(in) :: x(:)
    real(wp) :: max_component
    integer :: i

    max_component = 0.0_wp
    do i = 1, n
      max_component = max(max_component, abs(real(x(i), wp)), &
                          abs(aimag(x(i))))
    end do
  end function complex_max_abs_real_or_imag

  function real_upper_factor_is_singular(r, n) result(is_singular)
  !Function real_upper_factor_is_singular tests whether a real
  !upper-triangular factor can be used safely by an unguarded triangular solve.
  !The active factor scale is
  !
  !                     scale = max(abs(R(1:n,1:n)))
  !
  !and a diagonal element is considered unusable when
  !
  !        abs(R(i,i)) <= max(tiny(1.0_wp),epsilon(1.0_wp)*scale).
  !
  !This test detects exact rank deficiency as well as a diagonal that is lost
  !at the relative resolution of the stored factor. xGEQRF does not report
  !rank deficiency through INFO, so the test is required before DTRSV.
  !
  !  Input parameters:
  !    r - Real array containing the active upper-triangular factor.
  !    n - Active order of the factor.
  !
  !  Result:
  !    is_singular - True when at least one active diagonal element satisfies
  !                  the threshold above; false otherwise.
    real(wp), intent(in) :: r(:,:)
    integer, intent(in) :: n
    logical :: is_singular
    real(wp) :: factor_scale, threshold
    integer :: i

    factor_scale = maxval(abs(r(1:n,1:n)))
    threshold = max(tiny(1.0_wp), epsilon(1.0_wp) * factor_scale)
    is_singular = .false.
    do i = 1, n
      if (abs(r(i,i)) <= threshold) then
        is_singular = .true.
        return
      end if
    end do
  end function real_upper_factor_is_singular

  function complex_upper_factor_is_singular(r, n) result(is_singular)
  !Function complex_upper_factor_is_singular tests whether a complex
  !upper-triangular factor can be used safely by ZTRSV. The factor scale is the
  !largest complex modulus in R(1:n,1:n), and the diagonal threshold is
  !
  !        max(tiny(1.0_wp),epsilon(1.0_wp)*factor_scale).
  !
  !  Input parameters:
  !    r - Complex array containing the active upper-triangular factor.
  !    n - Active order of the factor.
  !
  !  Result:
  !    is_singular - True when a diagonal modulus is not greater than the
  !                  threshold; false otherwise.
    complex(wp), intent(in) :: r(:,:)
    integer, intent(in) :: n
    logical :: is_singular
    real(wp) :: factor_scale, threshold
    integer :: i

    factor_scale = maxval(abs(r(1:n,1:n)))
    threshold = max(tiny(1.0_wp), epsilon(1.0_wp) * factor_scale)
    is_singular = .false.
    do i = 1, n
      if (abs(r(i,i)) <= threshold) then
        is_singular = .true.
        return
      end if
    end do
  end function complex_upper_factor_is_singular

end module qrlinalg
