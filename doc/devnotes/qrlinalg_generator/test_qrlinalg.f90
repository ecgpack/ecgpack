!Program test_qrlinalg is the standalone MPI unit test and timing driver of the
!redesigned module qrlinalg (distributed base factorization with batched
!updates).
!
!  mpirun -np P ./test_qrlinalg [unit] [trace=FILE] [ref=FILE] [quick]
!                                [real|complex|both]
!                                [factor=auto|serial|distributed]
!  mpirun -np P ./test_qrlinalg timing [n=2850] [pending=32] [trials=10]
!                                [complex] [factor=...]
!  mpirun -np P ./test_qrlinalg factortime [n=2850] [reps=3] [complex]
!                                [skipserial]
!
!factor= selects the algorithm of the fresh factorization for the whole run
!(qr_set_factorization); the default is the automatic choice of the module.
!The environment variable ECG_QR_FACTOR overrides it inside the module. The
!direct factorization tests of the unit mode always exercise both algorithms.
!The factortime mode times factorize_fresh with both algorithms and reports
!what the automatic choice selects.
!
!The unit mode runs the suite of the real symmetric state (qr_real_state) and
!then the equivalent suite of the complex Hermitian state (qr_complex_state);
!real or complex selects one of them. The timing mode times the real state
!unless complex is given.
!
!Unit mode builds replicated physical matrices H and S (lower triangles
!authoritative, upper triangles filled with garbage), drives the module
!through random sequences of begin_active / replace_symmetric /
!append_symmetric / truncate / flush operations, and after every step compares
!the bordered solve with a dense reference (Gaussian elimination with partial
!pivoting on the full current matrix) and, every few steps, the eigenpair
!returned by solve with a dense reference inverse iteration. All ranks run the
!same deterministic sequence; rank zero prints. With trace=FILE rank zero
!writes the checked values, and with ref=FILE it compares them with the trace
!of another process count.
!
!The exit status is nonzero (MPI_Abort) when a check fails.
program test_qrlinalg
  use iso_fortran_env, only: int64
  use wp_def
  use qrlinalg
  implicit none

  integer, parameter :: MODEL_RITZ = 1, MODEL_GAUSS = 2
  real(wp), parameter :: EPS = epsilon(1.0_wp)
  !Kind of the reference evaluation of the shifted Rayleigh quotient: a real
  !kind with at least 30 decimal digits and more digits than wp, or wp itself
  !where no such kind exists (quadruple precision).
  integer, parameter :: HP_WIDER = selected_real_kind(max(precision(1.0_wp) + 1, 30))
  integer, parameter :: hp = merge(HP_WIDER, wp, HP_WIDER > 0)

  !MPI
  integer :: nprocs, myrank, ierr
  !Bookkeeping of the checks
  integer :: nfail = 0, nchecks = 0
  character(len=256) :: trace_file = ' ', ref_file = ' '
  logical :: quick = .false.
  !Trace of checked values (value, absolute tolerance)
  integer :: ntrace = 0
  real(wp), allocatable :: trace_val(:), trace_tol(:)
  !Random generator state (Park-Miller minimal standard)
  integer(int64) :: rng_state = 1_int64

  !Current problem
  integer :: model = MODEL_RITZ
  integer :: cap = 0, dmodes = 0, n = 0
  real(wp) :: shift = 0.0_wp
  real(wp), allocatable :: uvec(:,:), dval(:), expo(:), udelta(:), uener(:)
  !Parameters of the test problems: decades of the random norm factor of a
  !function, relative size of its unique part (small values make the basis
  !nearly linearly dependent), and smallest ratio of two Gaussian exponents.
  real(wp) :: scale_decades = 0.5_wp, dependence = 0.5_wp, min_ratio = 1.3_wp
  real(wp), allocatable :: hmat(:,:), smat(:,:)
  real(wp), allocatable :: hcol(:), scol(:)
  !Dense reference
  real(wp), allocatable :: aref(:,:), alu(:,:)
  integer, allocatable :: apiv(:)
  real(wp) :: anorm = 0.0_wp, ainv_norm = 0.0_wp, ainv2_norm = 0.0_wp
  real(wp) :: ref_res = 0.0_wp, ref_floor = 0.0_wp, ref_absform = 0.0_wp
  real(wp) :: ref_sabs = 0.0_wp
  real(wp), allocatable :: bvec(:), xvec(:), xref(:), wrk1(:), wrk2(:), wrk3(:)
  real(wp), allocatable :: wrk4(:), wrk5(:)
  real(wp), allocatable :: bcast_buf(:)
  real(wp) :: lam_ref = 0.0_wp
  !Module state
  type(qr_real_state) :: st
  !Complex Hermitian suite: problem data, dense reference, work arrays and
  !module state (see the header of the complex suite below). zdval, zudelta,
  !zuener, zrw4 and zrw5 are real.
  complex(wp), allocatable :: zuvec(:,:), zexpo(:)
  real(wp), allocatable :: zdval(:), zudelta(:), zuener(:)
  complex(wp), allocatable :: zhmat(:,:), zsmat(:,:), zhcol(:), zscol(:)
  complex(wp), allocatable :: zaref(:,:), zalu(:,:)
  integer, allocatable :: zapiv(:)
  complex(wp), allocatable :: zbvec(:), zxvec(:), zxref(:), zwrk1(:), zwrk2(:), zwrk3(:)
  complex(wp), allocatable :: zwrk4(:), zwrk5(:), zbcast_buf(:)
  real(wp), allocatable :: zrw4(:), zrw5(:)
  type(qr_complex_state) :: zst
  !Selected suite: 'both', 'real', or 'complex'.
  character(len=16) :: suite = 'both'
  !Statistics of one sequence
  real(wp) :: worst_res, worst_fwd, worst_lam, worst_vec, worst_eigres
  real(wp) :: worst_facres, worst_cond, worst_cond2, worst_norm
  !shifted_quotient: worst |lambda-reference eigenvalue|/(eps*absform) and
  !worst ratio of the difference from the quotient evaluated by the driver
  !in the kind hp to the rounding error bound of the evaluation in wp.
  real(wp) :: worst_shift_lam = 0.0_wp, worst_shift_ref = 0.0_wp
  real(wp) :: plain_lam, plain_vec, plain_eigres, plain_eigres_abs, worst_eigres_abs
  integer :: stat_polish_max, stat_polish_total
  integer :: stat_singular_flush, stat_resid_flush
  !allow_flush: the checks may apply the flush protocol of the application
  !(false in the stress cases, where every solve must succeed as it is).
  !flushed_in_check: a fresh factorization happened; the active set is gone.
  logical :: allow_flush = .true., flushed_in_check = .false.
  !ref_converged: the dense reference inverse iteration reached its stopping
  !rule. current_shift_rel: relative shift distance of the running sequence.
  logical :: ref_converged = .true.
  real(wp) :: current_shift_rel = 1.0e-3_wp
  integer :: stat_retarget
  real(wp) :: worst_facres_fresh
  integer :: stat_flush_requested, stat_flush_total, stat_max_pending
  integer :: stat_max_excluded, stat_below_base, stat_batch_full, stat_iter_max
  integer :: stat_solves, stat_eigs
  !Tolerance factors (multiples of the working precision)
  real(wp) :: tol_res_factor, tol_fwd_factor, tol_lam_factor, tol_vec_factor
  real(wp) :: tol_eigres_factor, tol_plain_factor

  character(len=256) :: arg
  character(len=16) :: mode
  integer :: iarg, timing_n, timing_pending, timing_trials, timing_reps
  !Requested factorization algorithm (factor=), and the algorithm forced by
  !the environment variable ECG_QR_FACTOR (-1 when it does not name one).
  integer :: factor_choice = QR_FACTOR_AUTO, env_factor = -1
  character(len=32) :: env_text
  logical :: skip_serial = .false.
  integer :: env_length, env_status
  !Direct factorization tests: bounds, as multiples of n*eps (factors) and of
  !n*eps*cond (solves), and the worst observed ratios (backward error,
  !columnwise backward error, orthogonality, solve; per algorithm).
  real(wp), parameter :: FACTOR_TOL = 4.0_wp, FACTOR_SOLVE_TOL = 1.0_wp
  real(wp) :: fstat(4, 2) = 0.0_wp, fstat_diff = 0.0_wp

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

  mode = 'unit'
  timing_n = 2850
  timing_pending = 32
  timing_trials = 10
  timing_reps = 3
  do iarg = 1, command_argument_count()
    call get_command_argument(iarg, arg)
    if (trim(arg) == 'unit' .or. trim(arg) == 'timing' .or. trim(arg) == 'factortime') then
      mode = trim(arg)
    else if (trim(arg) == 'factor=auto') then
      factor_choice = QR_FACTOR_AUTO
    else if (trim(arg) == 'factor=serial') then
      factor_choice = QR_FACTOR_SERIAL
    else if (trim(arg) == 'factor=distributed') then
      factor_choice = QR_FACTOR_DISTRIBUTED
    else if (trim(arg) == 'skipserial') then
      skip_serial = .true.
    else if (arg(1:5) == 'reps=') then
      read(arg(6:), *) timing_reps
    else if (trim(arg) == 'quick') then
      quick = .true.
    else if (trim(arg) == 'real' .or. trim(arg) == 'complex' .or. trim(arg) == 'both') then
      suite = trim(arg)
    else if (arg(1:6) == 'trace=') then
      trace_file = arg(7:)
    else if (arg(1:4) == 'ref=') then
      ref_file = arg(5:)
    else if (arg(1:2) == 'n=') then
      read(arg(3:), *) timing_n
    else if (arg(1:8) == 'pending=') then
      read(arg(9:), *) timing_pending
    else if (arg(1:7) == 'trials=') then
      read(arg(8:), *) timing_trials
    else
      if (myrank == 0) write(*, '(2a)') 'unknown argument: ', trim(arg)
      call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    end if
  end do

  !The module reads ECG_QR_FACTOR on rank zero only; the test does the same,
  !to know which algorithm a forced request really runs.
  if (myrank == 0) then
    call get_environment_variable('ECG_QR_FACTOR', env_text, env_length, env_status)
    if (env_status == 0 .and. env_length > 0) then
      do iarg = 1, len(env_text)
        if (env_text(iarg:iarg) >= 'A' .and. env_text(iarg:iarg) <= 'Z') &
          env_text(iarg:iarg) = achar(iachar(env_text(iarg:iarg)) + 32)
      end do
      if (trim(env_text) == 'serial') env_factor = QR_FACTOR_SERIAL
      if (trim(env_text) == 'distributed') env_factor = QR_FACTOR_DISTRIBUTED
      if (trim(env_text) == 'auto') env_factor = QR_FACTOR_AUTO
    end if
  end if
  if (nprocs > 1) call MPI_Bcast(env_factor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call qr_set_factorization(factor_choice)

  if (myrank == 0) then
    write(*, '(a,i0,a,i0,a,i0)') 'test_qrlinalg: wp=', wp, ' digits=', &
      precision(1.0_wp), ' processes=', nprocs
    write(*, '(a,i0,a,i0,a)') 'factorization request: ', factor_choice, &
      ' ECG_QR_FACTOR code: ', env_factor, ' (0 auto, 1 serial, 2 distributed, -1 not set)'
  end if

  if (trim(mode) == 'factortime') then
    if (trim(suite) == 'complex') then
      call zrun_factor_timing(timing_n, timing_reps)
    else
      call run_factor_timing(timing_n, timing_reps)
    end if
  else if (trim(mode) == 'timing') then
    if (trim(suite) == 'complex') then
      call zrun_timing(timing_n, timing_pending, timing_trials)
    else
      call run_timing(timing_n, timing_pending, timing_trials)
    end if
  else
    allocate(trace_val(200000), trace_tol(200000))
    if (myrank == 0) write(*, '(2a)') 'suite: ', trim(suite)
    if (trim(suite) /= 'complex') call run_unit_tests()
    if (trim(suite) /= 'real') call zrun_unit_tests()
    call finish_trace()
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, nfail, 1, MPI_INTEGER, MPI_MAX, &
                         MPI_COMM_WORLD, ierr)
    end if
    if (myrank == 0) then
      write(*, '(a,i0,a,i0)') 'TOTAL checks=', nchecks, ' failures=', nfail
      if (nfail == 0) then
        write(*, '(a)') 'RESULT: ALL TESTS PASSED'
      else
        write(*, '(a)') 'RESULT: FAILED'
      end if
    end if
    if (nfail /= 0) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call MPI_Finalize(ierr)

contains

  !=========================================================================
  ! Small utilities
  !=========================================================================

  function rnd() result(r)
    !Uniform pseudo-random number in (0,1); identical on all ranks.
    real(wp) :: r
    rng_state = mod(rng_state * 16807_int64, 2147483647_int64)
    r = real(rng_state, wp) / 2147483647.0_wp
  end function rnd

  function irnd(lo, hi) result(i)
    !Uniform pseudo-random integer in lo..hi.
    integer, intent(in) :: lo, hi
    integer :: i
    i = lo + int(rnd() * real(hi - lo + 1, wp))
    if (i > hi) i = hi
  end function irnd

  subroutine fail(message)
    character(len=*), intent(in) :: message
    nfail = nfail + 1
    if (myrank == 0 .and. nfail <= 60) write(*, '(2a)') '  FAIL: ', trim(message)
  end subroutine fail

  subroutine expect(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    nchecks = nchecks + 1
    if (.not. condition) call fail(message)
  end subroutine expect

  subroutine expect_info(info, expected, message)
    !Checks a status value and that it is identical on all ranks.
    integer, intent(in) :: info, expected
    character(len=*), intent(in) :: message
    character(len=200) :: text
    integer :: lohi(2)
    nchecks = nchecks + 1
    if (info /= expected) then
      write(text, '(a,a,i0,a,i0,3a)') trim(message), ': info=', info, &
        ' expected=', expected, ' (', trim(qr_status_message(info)), ')'
      call fail(text)
    end if
    if (nprocs > 1) then
      lohi(1) = info
      lohi(2) = -info
      call MPI_Allreduce(MPI_IN_PLACE, lohi, 2, MPI_INTEGER, MPI_MAX, &
                         MPI_COMM_WORLD, ierr)
      if (lohi(1) /= -lohi(2)) call fail(trim(message) // ': info differs between ranks')
    end if
  end subroutine expect_info

  subroutine expect_replicated(k, x, message)
    !Checks that a real vector is bitwise identical on all ranks.
    integer, intent(in) :: k
    real(wp), intent(in) :: x(k)
    character(len=*), intent(in) :: message
    integer :: i, nbad
    if (nprocs == 1) return
    nchecks = nchecks + 1
    bcast_buf(1:k) = x(1:k)
    call MPI_Bcast(bcast_buf, k, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    nbad = 0
    do i = 1, k
      if (.not. (bcast_buf(i) == x(i))) nbad = nbad + 1
    end do
    call MPI_Allreduce(MPI_IN_PLACE, nbad, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (nbad /= 0) call fail(trim(message) // ': result differs between ranks')
  end subroutine expect_replicated

  subroutine add_trace(value, tolerance)
    real(wp), intent(in) :: value, tolerance
    if (.not. allocated(trace_val)) return
    if (ntrace >= size(trace_val)) return
    ntrace = ntrace + 1
    trace_val(ntrace) = value
    trace_tol(ntrace) = tolerance
  end subroutine add_trace

  subroutine finish_trace()
    !Writes the trace and compares it with a reference trace (rank zero).
    integer :: i, k, nbad, ios
    real(wp) :: v, t, worst
    character(len=200) :: text
    if (myrank /= 0) return
    if (len_trim(trace_file) > 0) then
      open(21, file=trim(trace_file), status='replace', action='write')
      write(21, '(i0)') ntrace
      do i = 1, ntrace
        write(21, '(2es48.38e4)') trace_val(i), trace_tol(i)
      end do
      close(21)
    end if
    if (len_trim(ref_file) > 0) then
      open(22, file=trim(ref_file), status='old', action='read', iostat=ios)
      nchecks = nchecks + 1
      if (ios /= 0) then
        call fail('cannot open reference trace ' // trim(ref_file))
        return
      end if
      read(22, *) k
      if (k /= ntrace) then
        write(text, '(a,i0,a,i0)') 'trace length ', ntrace, ' differs from reference ', k
        call fail(text)
        close(22)
        return
      end if
      nbad = 0
      worst = 0.0_wp
      do i = 1, ntrace
        read(22, *) v, t
        t = 2.0_wp * max(t, trace_tol(i))
        if (.not. (abs(v - trace_val(i)) <= t)) nbad = nbad + 1
        if (t > 0.0_wp) worst = max(worst, abs(v - trace_val(i)) / t)
      end do
      close(22)
      write(*, '(a,i0,a,es9.2,a,i0)') 'cross-process trace comparison: values=', &
        ntrace, ' worst |diff|/tol=', worst, ' out of tolerance=', nbad
      if (nbad /= 0) call fail('trace differs from the reference process count')
    end if
  end subroutine finish_trace

  !=========================================================================
  ! Test problems
  !=========================================================================

  subroutine problem_allocate(model_in, capacity, modes)
    integer, intent(in) :: model_in, capacity, modes
    integer :: i, j
    call problem_free()
    model = model_in
    cap = capacity
    dmodes = modes
    allocate(hmat(cap, cap), smat(cap, cap), hcol(cap), scol(cap))
    allocate(aref(cap, cap), alu(cap, cap), apiv(cap))
    allocate(bvec(cap), xvec(cap), xref(cap), wrk1(cap), wrk2(cap), wrk3(cap))
    allocate(wrk4(cap), wrk5(cap))
    allocate(bcast_buf(cap))
    allocate(expo(cap), udelta(cap), uener(cap))
    allocate(uvec(max(1, dmodes), cap), dval(max(1, dmodes)))
    udelta = 0.0_wp
    uener = 0.0_wp
    !Garbage everywhere; the lower triangles are overwritten by the functions.
    do j = 1, cap
      do i = 1, cap
        hmat(i, j) = 7.7e33_wp * real(1 - 2 * mod(i + j, 2), wp)
        smat(i, j) = -3.3e31_wp * real(1 - 2 * mod(i, 2), wp)
      end do
    end do
    expo = 0.0_wp
    uvec = 0.0_wp
    dval = 0.0_wp
  end subroutine problem_allocate

  subroutine problem_free()
    if (allocated(hmat)) deallocate(hmat, smat, hcol, scol, aref, alu, apiv, &
                                    bvec, xvec, xref, wrk1, wrk2, wrk3, &
                                    wrk4, wrk5, bcast_buf, expo, uvec, dval, &
                                    udelta, uener)
  end subroutine problem_free

  subroutine new_function(k, norder)
    !Draws new parameters for function k; norder is the number of functions
    !considered to exist (used by the Gaussian model to keep exponents apart).
    !
    !Ritz model: function k is a vector with a "common" part u_k in a space of
    !dmodes modes with energies dval (one bound mode at -1, the others
    !positive and log-spaced up to dmax) and a "unique" part of amplitude
    !udelta(k) on a mode of its own with energy uener(k):
    !    S_ij = u_i.u_j + delta_ij*udelta(i)**2,
    !    H_ij = u_i.(dval*u_j) + delta_ij*udelta(i)**2*uener(i).
    !With more functions than common modes the lowest generalized eigenvalue
    !is close to -1 and well separated, and the size of the unique parts
    !(parameter dependence) controls how nearly dependent the basis is. A
    !random overall factor spreads the norms over 2*scale_decades decades.
    integer, intent(in) :: k, norder
    integer :: i, j, attempt
    real(wp) :: scale
    logical :: too_close
    if (model == MODEL_GAUSS) then
      !Exponent distributed log-uniformly over eight decades, not closer to
      !an existing exponent than the ratio min_ratio.
      do attempt = 1, 200
        expo(k) = 10.0_wp**(-2.0_wp + 8.0_wp * rnd())
        too_close = .false.
        do j = 1, norder
          if (j == k) cycle
          if (abs(log(expo(k) / expo(j))) < log(min_ratio)) too_close = .true.
        end do
        if (.not. too_close) exit
      end do
    else
      scale = 10.0_wp**(scale_decades * (2.0_wp * rnd() - 1.0_wp))
      uvec(1, k) = scale * (0.3_wp + 0.7_wp * rnd())
      if (rnd() < 0.5_wp) uvec(1, k) = -uvec(1, k)
      do i = 2, dmodes
        uvec(i, k) = scale * (2.0_wp * rnd() - 1.0_wp) / sqrt(dval(i))
      end do
      udelta(k) = scale * dependence * (0.5_wp + rnd())
      uener(k) = 2.0_wp + 10.0_wp * rnd()
    end if
  end subroutine new_function

  subroutine perturb_function(k, amount)
    !Small relative change of the parameters of function k (optimizer step).
    integer, intent(in) :: k
    real(wp), intent(in) :: amount
    integer :: i
    if (model == MODEL_GAUSS) then
      expo(k) = expo(k) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
    else
      do i = 1, dmodes
        uvec(i, k) = uvec(i, k) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
      end do
      udelta(k) = udelta(k) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
    end if
  end subroutine perturb_function

  subroutine function_column(k, norder)
    !Computes the full columns hcol(1:norder), scol(1:norder) of function k.
    integer, intent(in) :: k, norder
    integer :: i, j
    real(wp) :: ai, aj, hs, ss, sij, pi
    if (model == MODEL_GAUSS) then
      !Normalized s-type Gaussians exp(-a*r**2) for the hydrogen atom:
      !S = (2*sqrt(ai*aj)/(ai+aj))**1.5, T = 3*ai*aj/(ai+aj)*S,
      !V = -2*pi/(ai+aj)*(4*ai*aj/pi**2)**0.75.
      pi = 4.0_wp * atan(1.0_wp)
      ai = expo(k)
      do j = 1, norder
        aj = expo(j)
        sij = (2.0_wp * sqrt(ai * aj) / (ai + aj))**1.5_wp
        if (j == k) sij = 1.0_wp
        scol(j) = sij
        hcol(j) = 3.0_wp * ai * aj / (ai + aj) * sij - &
                  2.0_wp * pi / (ai + aj) * (4.0_wp * ai * aj / pi**2)**0.75_wp
      end do
    else
      do j = 1, norder
        hs = 0.0_wp
        ss = 0.0_wp
        do i = 1, dmodes
          ss = ss + uvec(i, k) * uvec(i, j)
          hs = hs + uvec(i, k) * dval(i) * uvec(i, j)
        end do
        if (j == k) then
          ss = ss + udelta(k)**2
          hs = hs + udelta(k)**2 * uener(k)
        end if
        hcol(j) = hs
        scol(j) = ss
      end do
    end if
  end subroutine function_column

  subroutine store_column(k, norder)
    !Writes hcol, scol into row and column k of the physical lower triangles.
    integer, intent(in) :: k, norder
    integer :: j
    do j = 1, norder
      hmat(max(j, k), min(j, k)) = hcol(j)
      smat(max(j, k), min(j, k)) = scol(j)
    end do
  end subroutine store_column

  subroutine problem_build(norder, dmax)
    !Creates norder functions and the physical matrices.
    integer, intent(in) :: norder
    real(wp), intent(in) :: dmax
    integer :: k
    if (model == MODEL_RITZ) then
      dval(1) = -1.0_wp
      do k = 2, dmodes
        dval(k) = 0.5_wp * (dmax / 0.5_wp)**(real(k - 2, wp) / real(max(1, dmodes - 2), wp))
      end do
    end if
    do k = 1, norder
      call new_function(k, k - 1)
    end do
    do k = 1, norder
      call function_column(k, k)
      call store_column(k, k)
    end do
    n = norder
  end subroutine problem_build

  !=========================================================================
  ! Dense reference
  !=========================================================================

  subroutine ref_symv(norder, a, x, y)
    !y = A*x using the lower triangle of a only.
    integer, intent(in) :: norder
    real(wp), intent(in) :: a(:,:), x(:)
    real(wp), intent(out) :: y(:)
    integer :: i, j
    y(1:norder) = 0.0_wp
    do j = 1, norder
      y(j) = y(j) + a(j, j) * x(j)
      do i = j + 1, norder
        y(i) = y(i) + a(i, j) * x(j)
        y(j) = y(j) + a(i, j) * x(i)
      end do
    end do
  end subroutine ref_symv

  subroutine ref_build()
    !Full current shifted matrix, its LU factorization with partial pivoting,
    !its Frobenius norm, and an estimate of the 2-norm of its inverse.
    integer :: i, j, k, ip, it
    real(wp) :: pmax, t, znorm, m11, m12, m22
    anorm = 0.0_wp
    do j = 1, n
      do i = j, n
        aref(i, j) = hmat(i, j) - shift * smat(i, j)
        aref(j, i) = aref(i, j)
      end do
    end do
    do j = 1, n
      do i = 1, n
        alu(i, j) = aref(i, j)
        anorm = anorm + aref(i, j)**2
      end do
    end do
    anorm = sqrt(anorm)
    do k = 1, n
      ip = k
      pmax = abs(alu(k, k))
      do i = k + 1, n
        if (abs(alu(i, k)) > pmax) then
          ip = i
          pmax = abs(alu(i, k))
        end if
      end do
      apiv(k) = ip
      if (ip /= k) then
        do j = 1, n
          t = alu(k, j)
          alu(k, j) = alu(ip, j)
          alu(ip, j) = t
        end do
      end if
      if (pmax > 0.0_wp) then
        do i = k + 1, n
          alu(i, k) = alu(i, k) / alu(k, k)
        end do
        do j = k + 1, n
          t = alu(k, j)
          do i = k + 1, n
            alu(i, j) = alu(i, j) - alu(i, k) * t
          end do
        end do
      end if
    end do
    !Two-vector subspace iteration on the inverse (the matrix is symmetric):
    !estimates of the largest and of the second largest singular value of
    !the inverse. The first gives the condition number of the shifted matrix;
    !the second governs the sensitivity of everything except the component
    !along the near-null vector (in particular of the eigenvector).
    do i = 1, n
      wrk1(i) = 1.0_wp + 0.37_wp * real(mod(7 * i, 11), wp)
      wrk2(i) = (1.0_wp + 0.21_wp * real(mod(5 * i, 13), wp)) * real(1 - 2 * mod(i, 2), wp)
    end do
    ainv_norm = 0.0_wp
    ainv2_norm = 0.0_wp
    do it = 1, 6
      znorm = sqrt(sum(wrk1(1:n)**2))
      wrk1(1:n) = wrk1(1:n) / znorm
      if (n > 1) then
        t = sum(wrk1(1:n) * wrk2(1:n))
        wrk2(1:n) = wrk2(1:n) - t * wrk1(1:n)
        znorm = sqrt(sum(wrk2(1:n)**2))
        wrk2(1:n) = wrk2(1:n) / znorm
      else
        wrk2(1) = 0.0_wp
      end if
      wrk4(1:n) = wrk1(1:n)
      wrk5(1:n) = wrk2(1:n)
      call ref_solve(wrk1)
      call ref_solve(wrk2)
      m11 = sum(wrk4(1:n) * wrk1(1:n))
      m12 = sum(wrk4(1:n) * wrk2(1:n))
      m22 = sum(wrk5(1:n) * wrk2(1:n))
      t = sqrt((0.5_wp * (m11 - m22))**2 + m12**2)
      ainv_norm = max(abs(0.5_wp * (m11 + m22) + t), abs(0.5_wp * (m11 + m22) - t))
      ainv2_norm = min(abs(0.5_wp * (m11 + m22) + t), abs(0.5_wp * (m11 + m22) - t))
    end do
    if (n == 1) ainv2_norm = ainv_norm
  end subroutine ref_build

  subroutine ref_solve(b)
    !Solves the current shifted system in place with the reference LU.
    real(wp), intent(inout) :: b(:)
    integer :: i, k
    real(wp) :: t
    do k = 1, n
      if (apiv(k) /= k) then
        t = b(k)
        b(k) = b(apiv(k))
        b(apiv(k)) = t
      end if
    end do
    do k = 1, n
      t = b(k)
      do i = k + 1, n
        b(i) = b(i) - alu(i, k) * t
      end do
    end do
    do k = n, 1, -1
      b(k) = b(k) / alu(k, k)
      t = b(k)
      do i = 1, k - 1
        b(i) = b(i) - alu(i, k) * t
      end do
    end do
  end subroutine ref_solve

  subroutine ref_abs_symv(norder, a, x, y)
    !y = |A|*|x| using the lower triangle of a only.
    integer, intent(in) :: norder
    real(wp), intent(in) :: a(:,:), x(:)
    real(wp), intent(out) :: y(:)
    integer :: i, j
    y(1:norder) = 0.0_wp
    do j = 1, norder
      y(j) = y(j) + abs(a(j, j) * x(j))
      do i = j + 1, norder
        y(i) = y(i) + abs(a(i, j) * x(j))
        y(j) = y(j) + abs(a(i, j) * x(i))
      end do
    end do
  end subroutine ref_abs_symv

  subroutine ref_eigen(x, lambda)
    !Dense reference inverse iteration with the current shift; S-normalized
    !eigenvector with positive largest component, physical Rayleigh quotient.
    !Also sets, for the reference eigenpair,
    !  ref_res     - relative eigenpair residual (definition of the module);
    !  ref_floor   - roundoff floor of that residual, eps*(|| |H||x| || +
    !                |lambda|*|| |S||x| ||)/(||H x|| + |lambda|*||S x||);
    !  ref_absform - |x|^T|H||x| + |lambda|*|x|^T|S||x|, the scale of the
    !                roundoff of the Rayleigh quotient of the S-normalized x;
    !  ref_sabs    - |x|^T|S||x|, the scale of the roundoff of x^T S x.
    real(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    integer :: i, it
    real(wp) :: change, change_prev, alpha, xs, vnorm2, hn, sn
    do i = 1, n
      wrk1(i) = 1.0_wp
    end do
    change_prev = huge(1.0_wp)
    ref_converged = .false.
    do it = 1, 120
      call ref_symv(n, smat, wrk1, x)
      call ref_solve(x)
      x(1:n) = x(1:n) / maxval(abs(x(1:n)))
      vnorm2 = sum(wrk1(1:n)**2)
      alpha = sum(x(1:n) * wrk1(1:n)) / vnorm2
      change = sqrt(sum((x(1:n) - alpha * wrk1(1:n))**2) / sum(x(1:n)**2))
      wrk1(1:n) = x(1:n)
      if ((change >= change_prev .and. change < 1.0e-3_wp) .or. change == 0.0_wp) then
        ref_converged = .true.
        exit
      end if
      change_prev = change
    end do
    call ref_symv(n, smat, x, wrk3)
    xs = sum(x(1:n) * wrk3(1:n))
    x(1:n) = x(1:n) / sqrt(xs)
    call fix_sign(x)
    call ref_symv(n, hmat, x, wrk2)
    call ref_symv(n, smat, x, wrk3)
    lambda = sum(x(1:n) * wrk2(1:n)) / sum(x(1:n) * wrk3(1:n))
    hn = sqrt(sum(wrk2(1:n)**2))
    sn = sqrt(sum(wrk3(1:n)**2))
    ref_res = sqrt(sum((wrk2(1:n) - lambda * wrk3(1:n))**2)) / (hn + abs(lambda) * sn)
    call ref_abs_symv(n, hmat, x, wrk4)
    call ref_abs_symv(n, smat, x, wrk5)
    ref_floor = EPS * (sqrt(sum(wrk4(1:n)**2)) + abs(lambda) * sqrt(sum(wrk5(1:n)**2))) / &
                (hn + abs(lambda) * sn)
    ref_sabs = sum(abs(x(1:n)) * wrk5(1:n))
    ref_absform = sum(abs(x(1:n)) * wrk4(1:n)) + abs(lambda) * ref_sabs
  end subroutine ref_eigen

  subroutine fix_sign(x)
    !Makes the component of largest magnitude positive.
    real(wp), intent(inout) :: x(:)
    integer :: i, imax
    imax = 1
    do i = 2, n
      if (abs(x(i)) > abs(x(imax))) imax = i
    end do
    if (x(imax) < 0.0_wp) x(1:n) = -x(1:n)
  end subroutine fix_sign

  !=========================================================================
  ! Checks against the reference
  !=========================================================================

  subroutine reset_stats()
    worst_res = 0.0_wp
    worst_fwd = 0.0_wp
    worst_lam = 0.0_wp
    worst_vec = 0.0_wp
    worst_eigres = 0.0_wp
    worst_eigres_abs = 0.0_wp
    worst_facres = 0.0_wp
    worst_cond = 0.0_wp
    worst_cond2 = 0.0_wp
    worst_norm = 0.0_wp
    worst_shift_lam = 0.0_wp
    worst_shift_ref = 0.0_wp
    plain_lam = 0.0_wp
    plain_vec = 0.0_wp
    plain_eigres = 0.0_wp
    plain_eigres_abs = 0.0_wp
    stat_flush_requested = 0
    stat_flush_total = 0
    stat_max_pending = 0
    stat_max_excluded = 0
    stat_below_base = 0
    stat_batch_full = 0
    stat_iter_max = 0
    stat_solves = 0
    stat_eigs = 0
    stat_polish_max = 0
    stat_polish_total = 0
    stat_singular_flush = 0
    stat_resid_flush = 0
    stat_retarget = 0
    worst_facres_fresh = 0.0_wp
  end subroutine reset_stats

  subroutine check_linear(label)
    !Bordered solve with a random right-hand side against the dense
    !reference solve of the current matrix.
    character(len=*), intent(in) :: label
    integer :: i, info
    real(wp) :: res, fwd, xnorm, bnorm, cond, tol_fwd
    character(len=200) :: text

    stat_max_pending = max(stat_max_pending, st%pending_count())
    stat_max_excluded = max(stat_max_excluded, st%excluded_count())
    if (st%order() < st%base_order()) stat_below_base = stat_below_base + 1
    call expect(st%order() == n, trim(label) // ': order() differs from the test order')

    do i = 1, n
      bvec(i) = 2.0_wp * rnd() - 1.0_wp
    end do
    call st%solve_shifted(bvec(1:n), xvec(1:n), info)
    if (info == QR_ERR_SINGULAR .and. allow_flush .and. st%pending_count() > 0) then
      !The committed block failed the pivot test: flush, as the application.
      call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': solve_shifted')
      stat_singular_flush = stat_singular_flush + 1
      call do_flush(0.0_wp, .false.)
      call st%solve_shifted(bvec(1:n), xvec(1:n), info)
    end if
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve_shifted')
    if (info /= QR_SUCCESS) return
    call expect_replicated(n, xvec, trim(label) // ': solve_shifted')
    xref(1:n) = bvec(1:n)
    call ref_solve(xref)

    xnorm = sqrt(sum(xvec(1:n)**2))
    bnorm = sqrt(sum(bvec(1:n)**2))
    wrk1(1:n) = matmul(aref(1:n, 1:n), xvec(1:n)) - bvec(1:n)
    res = sqrt(sum(wrk1(1:n)**2)) / (anorm * xnorm + bnorm)
    fwd = sqrt(sum((xvec(1:n) - xref(1:n))**2)) / sqrt(sum(xref(1:n)**2))
    cond = anorm * ainv_norm
    tol_fwd = tol_fwd_factor * EPS * cond
    worst_res = max(worst_res, res / EPS)
    worst_fwd = max(worst_fwd, fwd / (EPS * cond))
    worst_cond = max(worst_cond, cond)
    stat_solves = stat_solves + 1
    nchecks = nchecks + 2
    if (.not. (res <= tol_res_factor * EPS)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': solve residual/eps=', &
        res / EPS, ' cond=', cond
      call fail(text)
    end if
    if (.not. (fwd <= tol_fwd)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': solve forward error/(eps*cond)=', &
        fwd / (EPS * cond), ' cond=', cond
      call fail(text)
    end if
    !Trace: a norm-like functional of the solution. Its sensitivity to
    !roundoff is bounded by the forward error bound.
    call add_trace(xnorm, tol_fwd * xnorm)
  end subroutine check_linear

  subroutine check_eigen(label)
    !Eigenpair of solve against the dense reference inverse iteration. The
    !plain solve (no residual correction) is checked against loose bounds and
    !its errors are reported; the solve with residual correction is checked
    !against bounds derived from the conditioning of the problem:
    !  eigenvalue   |lambda-ref| <= tol_lam_factor*eps*ref_absform
    !               (roundoff of the Rayleigh quotient of an S-normalized x),
    !  eigenvector  ||x-ref||/||ref|| <= tol_vec_factor*eps*cond2, where cond2
    !               is the condition number of the shifted matrix on the
    !               complement of its near-null vector,
    !  residual     eig_residual <= tol_eigres_factor*max(ref_floor,ref_res),
    !  norm         |x^T S x - 1| <= tol_lam_factor*eps*ref_sabs.
    character(len=*), intent(in) :: label
    integer :: i, info, num_iter, npol
    real(wp) :: lambda, rel_acc, eig_res, dl, dv, xsx, cond2, res_scale
    character(len=200) :: text

    call ref_eigen(xref, lam_ref)
    if (.not. ref_converged .and. allow_flush) then
      !The shift of the last scheduled flush has become stale: the wanted
      !eigenvalue has drifted so far that even the dense reference iteration
      !converges too slowly to define a reference eigenvector. Retarget the
      !shift with a fresh factorization. The decision depends only on the
      !dense reference and is therefore the same for every process count.
      stat_retarget = stat_retarget + 1
      call do_flush(current_shift_rel, .true.)
      call ref_eigen(xref, lam_ref)
    end if
    nchecks = nchecks + 1
    if (.not. ref_converged) then
      call fail(trim(label) // ': the dense reference inverse iteration did not converge')
      return
    end if
    cond2 = anorm * ainv2_norm
    res_scale = max(ref_floor, ref_res)
    worst_cond2 = max(worst_cond2, cond2)
    stat_eigs = stat_eigs + 1

    !Plain solve, exactly the iteration of the design specification.
    do i = 1, n
      wrk1(i) = 1.0_wp
    end do
    call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                  rel_acc, num_iter, eig_res, info)
    call expect_replicated(n, xvec, trim(label) // ': solve')
    if (allow_flush .and. st%pending_count() > 0 .and. &
        (info == QR_ERR_NO_CONVERGENCE .or. info == QR_ERR_SINGULAR .or. &
         (info == QR_SUCCESS .and. .not. (eig_res <= tol_eigres_factor * res_scale)))) then
      !Retry policy of the application (DESIGN.md section 5): the solve did
      !not converge or the physical residual is too large while functions
      !are pending. The failure was reported by the module; a fresh
      !factorization with the same shift must give the right eigenpair.
      nchecks = nchecks + 1
      stat_resid_flush = stat_resid_flush + 1
      call do_flush(0.0_wp, .false.)
      do i = 1, n
        wrk1(i) = 1.0_wp
      end do
      call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                    rel_acc, num_iter, eig_res, info)
      call expect_replicated(n, xvec, trim(label) // ': solve after flush')
    end if
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve')
    if (info /= QR_SUCCESS .and. info /= QR_ERR_NO_CONVERGENCE) return
    call fix_sign(xvec)
    dl = abs(lambda - lam_ref)
    dv = sqrt(sum((xvec(1:n) - xref(1:n))**2)) / sqrt(sum(xref(1:n)**2))
    plain_lam = max(plain_lam, dl / (EPS * ref_absform))
    plain_vec = max(plain_vec, dv / (EPS * cond2))
    plain_eigres = max(plain_eigres, eig_res / res_scale)
    plain_eigres_abs = max(plain_eigres_abs, eig_res)
    stat_iter_max = max(stat_iter_max, num_iter)
    nchecks = nchecks + 3
    if (.not. (dl <= tol_plain_factor * tol_lam_factor * EPS * ref_absform)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve |lambda-ref|/(eps*absform)=', &
        dl / (EPS * ref_absform)
      call fail(text)
    end if
    if (.not. (dv <= tol_plain_factor * tol_vec_factor * EPS * cond2)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve eigenvector difference/(eps*cond2)=', &
        dv / (EPS * cond2)
      call fail(text)
    end if
    if (.not. (eig_res <= tol_plain_factor * tol_eigres_factor * res_scale)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve eigenpair residual/floor=', &
        eig_res / res_scale
      call fail(text)
    end if

    !Solve with residual correction.
    do i = 1, n
      wrk1(i) = 1.0_wp
    end do
    call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                  rel_acc, num_iter, eig_res, info, 0.0_wp, npol)
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve with correction')
    if (info /= QR_SUCCESS .and. info /= QR_ERR_NO_CONVERGENCE) return
    call expect_replicated(n, xvec, trim(label) // ': solve with correction')
    call ref_symv(n, smat, xvec, wrk2)
    xsx = sum(xvec(1:n) * wrk2(1:n))
    call fix_sign(xvec)
    dl = abs(lambda - lam_ref)
    dv = sqrt(sum((xvec(1:n) - xref(1:n))**2)) / sqrt(sum(xref(1:n)**2))
    worst_lam = max(worst_lam, dl / (EPS * ref_absform))
    worst_vec = max(worst_vec, dv / (EPS * cond2))
    worst_eigres = max(worst_eigres, eig_res / res_scale)
    worst_eigres_abs = max(worst_eigres_abs, eig_res)
    worst_norm = max(worst_norm, abs(xsx - 1.0_wp) / (EPS * ref_sabs))
    stat_polish_max = max(stat_polish_max, npol)
    stat_polish_total = stat_polish_total + npol
    nchecks = nchecks + 4
    if (.not. (dl <= tol_lam_factor * EPS * ref_absform)) then
      write(text, '(2a,es10.3,a,es24.16)') trim(label), ': |lambda-ref|/(eps*absform)=', &
        dl / (EPS * ref_absform), ' ref=', lam_ref
      call fail(text)
    end if
    if (.not. (dv <= tol_vec_factor * EPS * cond2)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': eigenvector difference/(eps*cond2)=', &
        dv / (EPS * cond2), ' cond2=', cond2
      call fail(text)
    end if
    if (.not. (eig_res <= tol_eigres_factor * res_scale)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': eigenpair residual/floor=', &
        eig_res / res_scale, ' floor=', res_scale
      call fail(text)
    end if
    if (.not. (abs(xsx - 1.0_wp) <= tol_lam_factor * EPS * ref_sabs)) then
      write(text, '(2a,es10.3)') trim(label), ': (x^T S x - 1)/(eps*sabs)=', &
        (xsx - 1.0_wp) / (EPS * ref_sabs)
      call fail(text)
    end if
    call add_trace(lambda, 2.0_wp * tol_lam_factor * EPS * ref_absform)
    call check_shifted_quotient(label)
  end subroutine check_eigen

  subroutine ref_shifted_quotient(x, sigma, lambda, bound)
    !Reference value of the shifted Rayleigh quotient
    !
    !        lambda = sigma + F/G,  F = x^T*(H-sigma*S)*x,  G = x^T*S*x,
    !
    !evaluated in the kind hp by an independent formula: the full symmetric
    !matrix is rebuilt from the lower triangle and summed row by row, where
    !the module sums the columns of the lower triangle and doubles the part
    !below the diagonal. bound is the rounding error bound of an evaluation in
    !the kind wp by any order of summation,
    !
    !        bound = (n+8)*eps*(Fabs + |F/G|*Gabs)/G + eps*|lambda| ,
    !        Fabs = |x|^T*(|H|+|sigma|*|S|)*|x| ,  Gabs = |x|^T*|S|*|x| :
    !
    !every term of F carries three rounded operations and at most n+P+2
    !additions follow, every term of G one operation less, the unit roundoff
    !is eps/2, and the final sum sigma + F/G is rounded once more.
    real(wp), intent(in) :: x(:), sigma
    real(wp), intent(out) :: lambda, bound
    integer :: i, j
    real(hp) :: f, g, fa, ga, row_f, row_g, row_fa, row_ga, hij, sij, sg, xj

    sg = real(sigma, hp)
    f = 0.0_hp
    g = 0.0_hp
    fa = 0.0_hp
    ga = 0.0_hp
    do i = 1, n
      row_f = 0.0_hp
      row_g = 0.0_hp
      row_fa = 0.0_hp
      row_ga = 0.0_hp
      do j = 1, n
        if (j <= i) then
          hij = real(hmat(i, j), hp)
          sij = real(smat(i, j), hp)
        else
          hij = real(hmat(j, i), hp)
          sij = real(smat(j, i), hp)
        end if
        xj = real(x(j), hp)
        row_f = row_f + (hij - sg * sij) * xj
        row_g = row_g + sij * xj
        row_fa = row_fa + (abs(hij) + abs(sg) * abs(sij)) * abs(xj)
        row_ga = row_ga + abs(sij) * abs(xj)
      end do
      f = f + real(x(i), hp) * row_f
      g = g + real(x(i), hp) * row_g
      fa = fa + abs(real(x(i), hp)) * row_fa
      ga = ga + abs(real(x(i), hp)) * row_ga
    end do
    lambda = real(sg + f / g, wp)
    bound = real(real(n + 8, hp) * real(EPS, hp) * (fa + abs(f / g) * ga) / g + &
                 real(EPS, hp) * abs(sg + f / g), wp)
    !Without a wider kind the reference carries the same rounding errors.
    if (hp == wp) bound = 2.0_wp * bound
  end subroutine ref_shifted_quotient

  subroutine expect_same_scalar(value, message)
    !Checks that a real scalar is bitwise identical on all ranks.
    real(wp), intent(in) :: value
    character(len=*), intent(in) :: message
    real(wp) :: buf(1)
    integer :: nbad
    if (nprocs == 1) return
    nchecks = nchecks + 1
    buf(1) = value
    call MPI_Bcast(buf, 1, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    nbad = 0
    if (.not. (buf(1) == value)) nbad = 1
    call MPI_Allreduce(MPI_IN_PLACE, nbad, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (nbad /= 0) call fail(trim(message) // ': result differs between ranks')
  end subroutine expect_same_scalar

  subroutine check_shifted_quotient(label)
    !shifted_quotient for the eigenvector left in xvec by check_eigen and for
    !a fixed vector that is far from every eigenvector (F is then not small,
    !so that a wrong weight of a matrix element cannot hide in the bound).
    !Each value must be identical on all ranks and agree with the quotient
    !of ref_shifted_quotient to within its bound, which does not depend on
    !the number of processes; the value of the eigenvector must also agree
    !with the reference eigenvalue as the eigenvalue of solve does. Both
    !values enter the trace that is compared between process counts.
    character(len=*), intent(in) :: label
    integer :: i

    call shifted_quotient_case(label, xvec, 1)
    do i = 1, n
      wrk1(i) = 1.0_wp + 0.5_wp * real(1 - 2 * mod(i, 2), wp) * &
                real(i, wp) / real(n, wp)
    end do
    call shifted_quotient_case(label, wrk1, 2)
  end subroutine check_shifted_quotient

  subroutine shifted_quotient_case(label, x, pass)
    !One vector of check_shifted_quotient; pass 1 is the eigenvector.
    character(len=*), intent(in) :: label
    real(wp), intent(in) :: x(:)
    integer, intent(in) :: pass
    integer :: info
    real(wp) :: lam_shift, lam_hp, bound, dl, sigma
    character(len=200) :: text

    sigma = st%get_shift()
    call st%shifted_quotient(hmat, smat, x(1:n), lam_shift, info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': shifted_quotient')
    if (info /= QR_SUCCESS) return
    call expect_same_scalar(lam_shift, trim(label) // ': shifted_quotient')
    call ref_shifted_quotient(x, sigma, lam_hp, bound)
    nchecks = nchecks + 1
    dl = abs(lam_shift - lam_hp)
    worst_shift_ref = max(worst_shift_ref, dl / bound)
    if (.not. (dl <= bound)) then
      write(text, '(2a,i0,a,es10.3,a,es24.16)') trim(label), ': shifted quotient, vector ', pass, &
        ': |lambda-driver|/bound=', dl / bound, ' driver=', lam_hp
      call fail(text)
    end if
    if (pass == 1) then
      nchecks = nchecks + 1
      dl = abs(lam_shift - lam_ref)
      worst_shift_lam = max(worst_shift_lam, dl / (EPS * ref_absform))
      if (.not. (dl <= tol_lam_factor * EPS * ref_absform)) then
        write(text, '(2a,es10.3)') trim(label), ': shifted quotient |lambda-ref|/(eps*absform)=', &
          dl / (EPS * ref_absform)
        call fail(text)
      end if
      call add_trace(lam_shift, 2.0_wp * tol_lam_factor * EPS * ref_absform)
    else
      call add_trace(lam_shift, 2.0_wp * bound)
    end if
  end subroutine shifted_quotient_case

  subroutine check_factor_residual(label)
    character(len=*), intent(in) :: label
    integer :: i, info
    real(wp) :: abs_res, rel_res, fwd, bound
    character(len=200) :: text
    do i = 1, n
      wrk1(i) = 1.0_wp + real(mod(17 * i, 23), wp) / 23.0_wp
      if (mod(i, 2) == 0) wrk1(i) = -wrk1(i)
    end do
    call st%factorization_residual(hmat, smat, wrk1(1:n), abs_res, rel_res, info, fwd)
    call expect_info(info, QR_SUCCESS, trim(label) // ': factorization_residual')
    !With an empty border the probe measures the fresh factorization and
    !must be at the level of the working precision (the bound of the
    !application is 1000*eps*n). With pending functions the probe solution
    !is of moderate size while the cancelling terms of block elimination are
    !not, so the probe shows the full growth factor ||A_base^{-1}*dA||; it is
    !bounded here by 1000 times the bound of the random-right-hand-side
    !solve and reported.
    if (info /= QR_SUCCESS) return
    nchecks = nchecks + 1
    if (st%pending_count() == 0) then
      worst_facres_fresh = max(worst_facres_fresh, rel_res / EPS)
      bound = 1000.0_wp * real(n, wp)
    else
      worst_facres = max(worst_facres, rel_res / EPS)
      bound = 1000.0_wp * tol_res_factor
    end if
    if (.not. (rel_res <= bound * EPS)) then
      write(text, '(2a,es10.3,a,es10.3,a,i0)') trim(label), ': factorization_residual/eps=', &
        rel_res / EPS, ' forward=', fwd, ' pending=', st%pending_count()
      call fail(text)
    end if
  end subroutine check_factor_residual

  !=========================================================================
  ! Driving the module
  !=========================================================================

  subroutine do_flush(shift_rel, retarget)
    !Fresh factorization of the current physical matrices. With retarget the
    !shift is placed shift_rel (relative) below the current reference
    !eigenvalue (the dense reference must be current); otherwise the shift is
    !kept. Only flushes that are scheduled by the random sequence itself
    !retarget the shift: flushes requested by the bookkeeping or by a
    !numerical safeguard may happen at different steps for different process
    !counts, and must not change the represented matrix.
    real(wp), intent(in) :: shift_rel
    logical, intent(in) :: retarget
    integer :: info
    if (retarget) then
      call ref_eigen(xref, lam_ref)
      shift = lam_ref - shift_rel * max(abs(lam_ref), 1.0e-3_wp)
    end if
    call st%factorize_fresh(hmat, smat, shift, info, n)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(st%pending_count() == 0 .and. st%excluded_count() == 0 .and. &
                st%base_order() == n .and. st%order() == n .and. st%is_valid(), &
                'state after factorize_fresh')
    stat_flush_total = stat_flush_total + 1
    flushed_in_check = .true.
    if (retarget) call ref_build()
  end subroutine do_flush

  subroutine announce(set, shift_rel, label)
    !begin_active with the flush protocol of the application: a fresh
    !factorization when the border is full (need_flush) or when the committed
    !block fails the pivot test (QR_ERR_SINGULAR).
    integer, intent(in) :: set(:)
    real(wp), intent(in) :: shift_rel
    character(len=*), intent(in) :: label
    logical :: need_flush
    integer :: info
    call st%begin_active(set, need_flush, info)
    if (info == QR_ERR_SINGULAR .and. allow_flush) then
      call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': begin_active')
      stat_singular_flush = stat_singular_flush + 1
      need_flush = .true.
    else
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active')
      if (need_flush) stat_flush_requested = stat_flush_requested + 1
    end if
    if (need_flush) then
      call do_flush(shift_rel, .false.)
      call st%begin_active(set, need_flush, info)
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active after flush')
      call expect(.not. need_flush, trim(label) // ': need_flush after a flush')
    end if
  end subroutine announce

  subroutine replace_function(k, label)
    !New column of function k in the physical matrices and in the module.
    integer, intent(in) :: k
    character(len=*), intent(in) :: label
    integer :: info
    call function_column(k, n)
    call store_column(k, n)
    call st%replace_symmetric(k, hcol(1:n), scol(1:n), info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': replace_symmetric')
  end subroutine replace_function

  subroutine append_function(label)
    !Appends function n+1 (parameters already drawn).
    character(len=*), intent(in) :: label
    integer :: info
    call function_column(n + 1, n + 1)
    call store_column(n + 1, n + 1)
    call st%append_symmetric(hcol(1:n + 1), scol(1:n + 1), info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': append_symmetric')
    n = n + 1
  end subroutine append_function

  subroutine run_sequence(name, model_in, n_init, capacity, max_pending, nsteps, &
                          shift_rel, dmax, seed)
    !Random sequence of operations with a check after every step.
    character(len=*), intent(in) :: name
    integer, intent(in) :: model_in, n_init, capacity, max_pending, nsteps, seed
    real(wp), intent(in) :: shift_rel, dmax
    integer :: step, info, nact, act(8), set(8), nset, i, j, k, op, t, base, nmin
    integer :: nfail_before
    logical :: need_flush, found
    real(wp) :: r
    character(len=80) :: label

    nfail_before = nfail
    allow_flush = .true.
    current_shift_rel = shift_rel
    rng_state = int(seed, int64)
    call reset_stats()
    call problem_allocate(model_in, capacity, max(2, n_init / 3))
    call problem_build(n_init, dmax)
    nmin = max(3, n_init / 2)
    call st%initialize(capacity, max_pending, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    !Initial shift well below the spectrum, then a shift close to the lowest
    !eigenvalue found with it.
    shift = -2.0_wp * maxval(abs(dval)) - 10.0_wp
    if (model == MODEL_GAUSS) shift = -5.0_wp
    call ref_build()
    call do_flush(shift_rel, .true.)
    call check_linear('fresh')
    call check_eigen('fresh')
    call check_factor_residual('fresh')
    nact = 0

    do step = 1, nsteps
      write(label, '(a,a,i0)') trim(name), ' step ', step
      r = rnd()
      op = 1
      if (r > 0.22_wp) op = 2
      if (r > 0.55_wp) op = 3
      if (r > 0.70_wp) op = 4
      if (r > 0.80_wp) op = 5
      if (r > 0.86_wp) op = 6
      if (r > 0.95_wp) op = 7
      if (op == 2 .and. nact == 0) op = 1

      select case (op)
      case (1)
        !New active set of one or three existing functions: ascending,
        !descending, or scattered, optionally sharing members with the
        !previous active set.
        nset = 1
        if (n >= 6 .and. rnd() < 0.6_wp) nset = 3
        t = irnd(1, 3)
        if (nset == 3 .and. t == 1) then
          k = irnd(1, n - 2)
          set(1:3) = [k, k + 1, k + 2]
        else if (nset == 3 .and. t == 2) then
          k = irnd(1, n - 2)
          set(1:3) = [k + 2, k + 1, k]
        else
          i = 0
          do while (i < nset)
            k = irnd(1, n)
            found = .false.
            do j = 1, i
              if (set(j) == k) found = .true.
            end do
            if (.not. found) then
              i = i + 1
              set(i) = k
            end if
          end do
        end if
        !Shared intersection with the previous set.
        if (nact > 0 .and. rnd() < 0.5_wp) then
          k = act(irnd(1, nact))
          found = .false.
          do j = 1, nset
            if (set(j) == k) found = .true.
          end do
          if (.not. found .and. k <= n) set(irnd(1, nset)) = k
        end if
        call announce(set(1:nset), shift_rel, label)
        !The matrix is unchanged by the announcement.
        flushed_in_check = .false.
        call check_linear(trim(label) // ' (announced)')
        if (flushed_in_check) call announce(set(1:nset), shift_rel, label)
        do i = 1, nset
          call new_function(set(i), n)
          call replace_function(set(i), label)
        end do
        nact = nset
        act(1:nset) = set(1:nset)
      case (2)
        !Trial: new columns for one or all members of the active set.
        if (rnd() < 0.5_wp) then
          k = act(irnd(1, nact))
          call perturb_function(k, 0.05_wp)
          call replace_function(k, label)
        else
          do i = 1, nact
            call new_function(act(i), n)
          end do
          do i = 1, nact
            call replace_function(act(i), label)
          end do
        end if
      case (3)
        !Append one or two functions, sometimes together with an existing
        !active function.
        t = irnd(1, 2)
        if (n + t > capacity) t = capacity - n
        if (t <= 0) then
          call st%truncate(n - 1, info)
          if (info == QR_ERR_BATCH_FULL) then
            stat_batch_full = stat_batch_full + 1
            n = n - 1
            call do_flush(shift_rel, .false.)
          else
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate at capacity')
            n = n - 1
          end if
          nact = 0
          call st%begin_active(set(1:0), need_flush, info)
          call expect_info(info, QR_SUCCESS, trim(label) // ': commit')
        else
          nset = 0
          if (rnd() < 0.3_wp) then
            nset = 1
            set(1) = irnd(1, n)
          end if
          do i = 1, t
            set(nset + i) = n + t + 1 - i
          end do
          call announce(set(1:nset + t), shift_rel, label)
          if (nset == 1) then
            call new_function(set(1), n)
            call replace_function(set(1), label)
          end if
          nact = nset
          act(1:nset) = set(1:nset)
          do i = 1, t
            call new_function(n + 1, n)
            call append_function(label)
            nact = nact + 1
            act(nact) = n
          end do
        end if
      case (4)
        !Truncate one to three trailing functions, possibly below the base.
        t = irnd(1, 3)
        if (n - t < nmin) t = n - nmin
        if (t > 0) then
          call st%truncate(n - t, info)
          if (info == QR_ERR_BATCH_FULL) then
            stat_batch_full = stat_batch_full + 1
            call expect(st%order() == n, trim(label) // ': order changed by a refused truncate')
            n = n - t
            call do_flush(shift_rel, .false.)
            nact = 0
          else
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate')
            n = n - t
          end if
          !Forget active indices that no longer exist.
          j = 0
          do i = 1, nact
            if (act(i) <= n) then
              j = j + 1
              act(j) = act(i)
            end if
          end do
          nact = j
        end if
      case (5)
        !Basis-enlargement pattern: several candidates for the same slot.
        if (n < capacity) then
          base = n
          do i = 1, 4
            call st%truncate(base, info)
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate to base')
            n = base
            set(1) = base + 1
            call announce(set(1:1), shift_rel, label)
            call new_function(n + 1, n)
            call append_function(label)
            if (i < 4) then
              call ref_build()
              call check_linear(trim(label) // ' (candidate)')
            end if
          end do
          nact = 1
          act(1) = n
        end if
      case (6)
        !Flush at random.
        call do_flush(shift_rel, .true.)
        call check_factor_residual(trim(label) // ' (fresh)')
        nact = 0
      case (7)
        !Commit the tail by announcing the empty set.
        call st%begin_active(set(1:0), need_flush, info)
        call expect_info(info, QR_SUCCESS, trim(label) // ': commit')
        call expect(.not. need_flush .and. st%active_count() == 0, trim(label) // ': commit state')
        nact = 0
      end select

      call ref_build()
      flushed_in_check = .false.
      call check_linear(label)
      if (mod(step, 4) == 0) call check_eigen(label)
      if (mod(step, 25) == 0) call check_factor_residual(label)
      !A flush inside the checks emptied the active set.
      if (flushed_in_check) nact = 0
    end do
    !Numerically triggered flushes must stay exceptional.
    call expect(10 * (stat_resid_flush + stat_singular_flush) <= stat_eigs, &
                trim(name) // ': too many numerically triggered flushes')

    call report_sequence(name, nfail - nfail_before)
    call st%clear()
  end subroutine run_sequence

  subroutine report_sequence(name, failures)
    character(len=*), intent(in) :: name
    integer, intent(in) :: failures
    if (myrank /= 0) return
    write(*, '(a)') 'sequence ' // trim(name)
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '   solves=', stat_solves, &
      ' eigenpairs=', stat_eigs, ' flushes=', stat_flush_total, &
      ' (requested by begin_active: ', stat_flush_requested, &
      ', truncate batch-full: ', stat_batch_full, ') max pending=', &
      stat_max_pending, ' max excluded=', stat_max_excluded
    write(*, '(a,i0,a,i0,a,es9.2,a,es9.2)') '   checks with order below base=', stat_below_base, &
      ' max iterations=', stat_iter_max, ' max cond(A)=', worst_cond, ' max cond2(A)=', worst_cond2
    write(*, '(a,es9.2,a,es9.2)') '   linear solve: residual/eps=', worst_res, &
      ' forward/(eps*cond)=', worst_fwd
    write(*, '(a,es9.2,a,es9.2)') '   factorization_residual/eps: empty border=', worst_facres_fresh, &
      ' with pending functions=', worst_facres
    write(*, '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)') '   plain solve:     lambda/(eps*absform)=', plain_lam, &
      ' vector/(eps*cond2)=', plain_vec, ' residual/floor=', plain_eigres, ' residual=', plain_eigres_abs
    write(*, '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)') '   corrected solve: lambda/(eps*absform)=', worst_lam, &
      ' vector/(eps*cond2)=', worst_vec, ' residual/floor=', worst_eigres, ' residual=', worst_eigres_abs
    write(*, '(a,es9.2,a,es9.2)') '   shifted quotient: lambda/(eps*absform)=', worst_shift_lam, &
      ' difference from the driver/bound=', worst_shift_ref
    write(*, '(a,es9.2,a,i0,a,i0)') '   corrected solve: norm error/(eps*sabs)=', worst_norm, &
      ' max correction steps=', stat_polish_max, ' total correction steps=', stat_polish_total
    write(*, '(a,i0,a,i0,a,i0)') '   numerically triggered flushes: pivot test=', stat_singular_flush, &
      ' residual/no convergence=', stat_resid_flush, ' ; shift retargeted (stale shift)=', stat_retarget
    write(*, '(a,i0)') '   failures=', failures
  end subroutine report_sequence

  subroutine run_stress(name, model_in, n_init, max_pending, npending, shift_rel, &
                        dmax, seed)
    !Near-singular stress: the shift lies shift_rel (relative) below the
    !lowest eigenvalue and npending functions are replaced by slightly
    !changed ones in active sets of three without any flush, so that the
    !shifted matrix stays nearly singular throughout.
    character(len=*), intent(in) :: name
    integer, intent(in) :: model_in, n_init, max_pending, npending, seed
    real(wp), intent(in) :: shift_rel, dmax
    integer :: info, i, k, set(3), done, nfail_before, stride
    logical :: need_flush_dummy
    character(len=80) :: label

    nfail_before = nfail
    allow_flush = .false.
    current_shift_rel = shift_rel
    rng_state = int(seed, int64)
    call reset_stats()
    call problem_allocate(model_in, n_init + 4, max(2, n_init / 3))
    call problem_build(n_init, dmax)
    call st%initialize(n_init + 4, max_pending, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    shift = -2.0_wp * maxval(abs(dval)) - 10.0_wp
    if (model == MODEL_GAUSS) shift = -5.0_wp
    call ref_build()
    !Two passes: the first shift is only a rough guide to the eigenvalue.
    call do_flush(1.0e-2_wp, .true.)
    call do_flush(shift_rel, .true.)
    call check_linear('fresh')
    call check_eigen('fresh')
    done = 0
    stride = max(1, (n - 3) / max(1, npending / 3 + 1))
    k = 1
    do while (done < npending)
      write(label, '(a,a,i0)') trim(name), ' pending ', done + 3
      set = [k + 2, k, k + 1]
      k = k + max(3, stride)
      if (k + 2 > n) k = 1
      call st%begin_active(set, need_flush_dummy, info)
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active')
      call expect(.not. need_flush_dummy, trim(label) // ': unexpected need_flush')
      do i = 1, 3
        call perturb_function(set(i), 2.0e-3_wp)
        call replace_function(set(i), label)
      end do
      done = done + 3
      call ref_build()
      call check_linear(label)
      call check_eigen(label)
    end do
    call expect(st%pending_count() >= npending, trim(name) // ': pending count')
    call check_factor_residual(name)
    !Trials on the last active set with everything else pending.
    do i = 1, 6
      write(label, '(a,a,i0)') trim(name), ' trial ', i
      call perturb_function(set(1 + mod(i, 3)), 1.0e-3_wp)
      call replace_function(set(1 + mod(i, 3)), label)
      call ref_build()
      call check_linear(label)
      call check_eigen(label)
    end do
    if (myrank == 0) then
      write(*, '(a,es9.2,a,es9.2)') 'stress ' // trim(name) // ': relative distance shift-eigenvalue=', &
        abs(lam_ref - shift) / abs(lam_ref), ' lambda=', lam_ref
    end if
    call report_sequence(name, nfail - nfail_before)
    call st%clear()
  end subroutine run_stress

  subroutine run_error_paths()
    !Error paths and state preservation.
    integer :: info, num_iter, nfail_before, i, j
    logical :: need_flush
    real(wp) :: lambda, rel_acc, eig_res, abs_res, rel_res
    real(wp) :: h4(4, 4), s4(4, 4), h2(2, 2), s2(2, 2), v4(4), x4(4)
    type(qr_real_state) :: e

    nfail_before = nfail
    rng_state = 4711_int64
    call problem_allocate(MODEL_RITZ, 10, 3)
    call problem_build(8, 100.0_wp)

    !Uninitialized state.
    call e%factorize_fresh(hmat, smat, -50.0_wp, info, 8)
    call expect_info(info, QR_ERR_INVALID_STATE, 'factorize_fresh on uninitialized state')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'begin_active on uninitialized state')
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'solve on uninitialized state')
    call e%shifted_quotient(hmat, smat, wrk1(1:8), lambda, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'shifted_quotient on uninitialized state')
    call expect(.not. e%is_valid() .and. e%order() == 0, 'queries on uninitialized state')

    !Initialization arguments.
    call e%initialize(0, 4, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'initialize capacity=0')
    call e%initialize(10, 0, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'initialize max_pending=0')
    call e%initialize(10, 99, info)
    call expect_info(info, QR_SUCCESS, 'initialize with max_pending above capacity')
    call expect(e%max_pending() == 10 .and. e%get_capacity() == 10, 'max_pending clamp')
    call e%initialize(10, 3, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    call expect(e%max_pending() == 3 .and. .not. e%is_valid(), 'state after initialize')
    call e%replace_symmetric(1, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'replace before factorization')

    !Factorization arguments.
    call e%factorize_fresh(hmat, smat, -50.0_wp, info, 11)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'factorize_fresh order above capacity')
    call e%factorize_fresh(hmat, smat, -50.0_wp, info, 0)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'factorize_fresh order zero')
    call e%factorize_fresh(hmat(1:8, 1:7), smat(1:8, 1:7), -50.0_wp, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorize_fresh non-square')
    call e%factorize_fresh(hmat(1:7, 1:7), smat, -50.0_wp, info, 8)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorize_fresh too small H')
    shift = -50.0_wp
    call e%factorize_fresh(hmat, smat, shift, info, 8)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(e%is_valid() .and. e%order() == 8 .and. e%base_order() == 8 .and. &
                e%get_shift() == shift, 'queries after factorize_fresh')

    !begin_active arguments.
    call e%begin_active([0], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active index 0')
    call e%begin_active([3, 3], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active repeated index')
    call e%begin_active([10], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active gap in appended range')
    call e%begin_active([9, 10, 11], need_flush, info)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'begin_active append beyond capacity')
    call e%begin_active([1, 2, 3, 4], need_flush, info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'begin_active set larger than max_pending')
    call expect(e%pending_count() == 0 .and. e%active_count() == 0, 'state after refused begin_active')

    !replace_symmetric arguments.
    call function_column(2, 8)
    call e%replace_symmetric(2, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace of an unannounced function')
    call e%begin_active([3, 1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [3,1,2]')
    call expect(.not. need_flush .and. e%pending_count() == 3 .and. e%active_count() == 3 &
                .and. e%excluded_count() == 3, 'state after begin_active')
    call e%replace_symmetric(0, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace index 0')
    call e%replace_symmetric(9, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace index above order')
    call e%replace_symmetric(2, hcol(1:7), scol(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'replace with short column')
    call e%replace_symmetric(5, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace of a clean function')
    call e%replace_symmetric(2, hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_SUCCESS, 'replace with an unchanged column')
    call expect(e%get_update_count() == 1_int64 .and. e%get_updates_since_fresh() == 1_int64, &
                'update counters')

    !Batch full: need_flush and unchanged state.
    call e%begin_active([4], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active on a full border')
    call expect(need_flush, 'need_flush on a full border')
    call expect(e%pending_count() == 3 .and. e%active_count() == 3, 'state after need_flush')
    call e%begin_active([2, 3, 1], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active same set')
    call expect(.not. need_flush .and. e%active_count() == 3, 'same set is a no-op')
    call e%begin_active([9], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active append on a full border')
    call expect(need_flush, 'need_flush for an append on a full border')
    call function_column(9, 9)
    call new_function(9, 8)
    call function_column(9, 9)
    call e%append_symmetric(hcol(1:9), scol(1:9), info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'unannounced append on a full border')
    call expect(e%order() == 8, 'order after refused append')
    call e%append_symmetric(hcol(1:8), scol(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'append with short column')

    !truncate / delete arguments.
    call e%truncate(0, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'truncate to zero')
    call e%truncate(9, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'truncate above order')
    call e%delete_symmetric(3, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'delete interior function')
    call e%truncate(7, info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'truncate with full D0')
    call expect(e%order() == 8, 'order after refused truncate')
    call e%truncate(8, info)
    call expect_info(info, QR_SUCCESS, 'truncate no-op')

    !solve arguments.
    wrk1(1:8) = 0.0_wp
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_ZERO_INITIAL_VECTOR, 'solve with zero start vector')
    call expect(lambda == 0.0_wp .and. num_iter == 0 .and. all(xvec(1:8) == 0.0_wp), &
                'outputs after zero start vector')
    wrk1(1:8) = 1.0_wp
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e-8_wp, 0, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'solve with max_iter=0')
    call e%solve(hmat, smat, wrk1(1:7), xvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve with short start vector')
    call e%solve(hmat(1:7, 1:7), smat, wrk1(1:8), xvec(1:8), lambda, 1.0e-8_wp, 10, 0, &
                 rel_acc, num_iter, eig_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve with small H')
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e-30_wp, 1, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_NO_CONVERGENCE, 'solve with max_iter=1')
    call expect(num_iter == 1 .and. maxval(abs(xvec(1:8))) > 0.0_wp .and. &
                eig_res < 1.0_wp, 'approximation after no convergence')
    !Normalization modes with the far shift (slow but sure convergence).
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e3_wp * EPS, 5000, 1, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_SUCCESS, 'solve norm_mode=1')
    call expect(abs(sum(xvec(1:8)**2) - 1.0_wp) < 100.0_wp * EPS, 'norm_mode=1 normalization')
    call e%solve(hmat, smat, wrk1(1:8), xvec(1:8), lambda, 1.0e3_wp * EPS, 5000, 2, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_SUCCESS, 'solve norm_mode=2')
    call expect(abs(maxval(abs(xvec(1:8))) - 1.0_wp) < 100.0_wp * EPS, 'norm_mode=2 normalization')
    call e%factorization_residual(hmat, smat, wrk1(1:7), abs_res, rel_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorization_residual short probe')
    wrk2(1:8) = 0.0_wp
    call e%factorization_residual(hmat, smat, wrk2(1:8), abs_res, rel_res, info)
    call expect_info(info, QR_ERR_ZERO_INITIAL_VECTOR, 'factorization_residual zero probe')
    call e%solve_shifted(wrk1(1:7), xvec(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve_shifted short right-hand side')
    lambda = 1.0_wp
    call e%shifted_quotient(hmat, smat, wrk1(1:7), lambda, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'shifted_quotient short vector')
    call expect(lambda == 0.0_wp, 'shifted_quotient returns zero on a dimension error')
    call e%shifted_quotient(hmat(1:7, 1:8), smat, wrk1(1:8), lambda, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'shifted_quotient short H')
    wrk2(1:8) = 0.0_wp
    lambda = 1.0_wp
    call e%shifted_quotient(hmat, smat, wrk2(1:8), lambda, info)
    call expect_info(info, QR_ERR_NONPOSITIVE_OVERLAP, 'shifted_quotient zero vector')
    call expect(lambda == 0.0_wp, 'shifted_quotient returns zero for a zero vector')

    !Append at capacity.
    call e%factorize_fresh(hmat, smat, shift, info, 8)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(e%get_updates_since_fresh() == 0_int64 .and. e%get_update_count() == 1_int64, &
                'counters after factorize_fresh')
    n = 8
    do i = 9, 10
      call new_function(i, i - 1)
      call function_column(i, i)
      call store_column(i, i)
      call e%append_symmetric(hcol(1:i), scol(1:i), info)
      call expect_info(info, QR_SUCCESS, 'unannounced append')
    end do
    call expect(e%order() == 10 .and. e%pending_count() == 2 .and. e%active_count() == 2, &
                'state after unannounced appends')
    call e%append_symmetric(hcol(1:10), scol(1:10), info)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'append at capacity')
    call e%delete_symmetric(10, info)
    call expect_info(info, QR_SUCCESS, 'delete last function')
    call expect(e%order() == 9 .and. e%pending_count() == 1, 'state after delete')

    !Exactly singular base matrix: H = shift*S.
    do j = 1, 4
      do i = 1, 4
        s4(i, j) = 0.0_wp
        h4(i, j) = 0.0_wp
      end do
      s4(j, j) = 1.0_wp
      h4(j, j) = 2.0_wp
    end do
    v4 = 1.0_wp
    call e%initialize(4, 2, info)
    call e%factorize_fresh(h4, s4, 2.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh of a singular matrix')
    call e%solve(h4, s4, v4, x4, lambda, 1.0e-8_wp, 10, 0, rel_acc, num_iter, eig_res, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve with singular R')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular R')

    !Singular small matrix ZDD: A is nonsingular but A without function 1 is
    !singular, A = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1].
    h4 = 0.0_wp
    h4(2, 1) = 1.0_wp
    h4(3, 3) = 1.0_wp
    h4(4, 4) = 1.0_wp
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh (ZDD case)')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted before the singular border')
    call expect(maxval(abs(x4 - 1.0_wp)) < 100.0_wp * EPS, 'solution of the permutation system')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular committed block')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular committed block')
    !Columns can still be stored; only the solve needs the singular block.
    call e%replace_symmetric(1, h4(:, 1), s4(:, 1), info)
    call expect_info(info, QR_SUCCESS, 'replace with singular committed block')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted after the replace')
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'recovery by factorize_fresh')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after recovery')

    !Singular committed Schur complement Scc: commit a replacement that makes
    !the matrix without the next active function singular. A = diag(1,1,1,1);
    !replace column 1 by e2 (A stays nonsingular with column 2 = e1 later).
    do j = 1, 4
      do i = 1, 4
        h4(i, j) = 0.0_wp
      end do
      h4(j, j) = 1.0_wp
    end do
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call e%begin_active([1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [1,2]')
    !New matrix: rows/columns 1 and 2 form [0 1; 1 0].
    h4(1, 1) = 0.0_wp
    h4(2, 1) = 1.0_wp
    h4(2, 2) = 0.0_wp
    v4 = [0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    x4 = 0.0_wp
    call e%replace_symmetric(1, v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'replace 1')
    v4 = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    call e%replace_symmetric(2, v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'replace 2')
    v4 = 1.0_wp
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted with the 2x2 swap')
    call expect(maxval(abs(x4 - 1.0_wp)) < 100.0_wp * EPS, 'solution with the 2x2 swap')
    !Committing function 1 alone leaves Scc = (0): singular.
    call e%begin_active([2], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular Scc')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular Scc')
    !A different partition is usable again.
    call e%begin_active([1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active back to [1,2]')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after repartition')

    !Singular active Schur complement Sigma: the represented matrix itself
    !becomes exactly singular, A = [1 1; 1 1].
    !The base is the identity, so that every stored quantity is exact.
    h2 = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
    s2 = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
    call e%initialize(2, 2, info)
    call e%factorize_fresh(h2, s2, 0.0_wp, info)
    call e%begin_active([2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [2] (Sigma case)')
    call e%replace_symmetric(2, [1.0_wp, 1.0_wp], [0.0_wp, 0.0_wp], info)
    call expect_info(info, QR_SUCCESS, 'replace giving a singular matrix')
    call e%solve_shifted([1.0_wp, 2.0_wp], x4(1:2), info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular Sigma')
    call e%replace_symmetric(2, [1.0_wp, 3.0_wp], [0.0_wp, 0.0_wp], info)
    call expect_info(info, QR_SUCCESS, 'replace giving a regular matrix')
    call e%solve_shifted([2.0_wp, 4.0_wp], x4(1:2), info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after the regular replace')
    call expect(maxval(abs(x4(1:2) - 1.0_wp)) < 100.0_wp * EPS, 'solution after the regular replace')

    call e%clear()
    call expect(.not. e%is_valid() .and. e%get_capacity() == 0 .and. e%max_pending() == 0, &
                'state after clear')
    call expect(trim(qr_status_message(QR_ERR_BATCH_FULL)) /= trim(qr_status_message(-1)), &
                'status message for QR_ERR_BATCH_FULL')
    if (myrank == 0) write(*, '(a,i0)') 'error paths: failures=', nfail - nfail_before
  end subroutine run_error_paths

  subroutine run_unit_tests()
    integer :: steps
    logical :: high_precision
    !Tolerances, as multiples of the working precision eps.
    !
    !Eigenpairs (the quantity the application uses) are held to bounds that
    !follow from the conditioning of each individual problem, see
    !check_eigen: 100*eps*(roundoff scale of the Rayleigh quotient) for the
    !eigenvalue, 100*eps*cond2 for the eigenvector, and 100 times the larger
    !of the roundoff floor of the residual and the residual of the dense
    !reference eigenvector for the eigenpair residual.
    !
    !The raw bordered linear solve is not backward stable in the usual sense:
    !block elimination on the nearly singular base matrix has a growth factor
    !(the size of the cancelling terms Zb*yE and Gs*y relative to the
    !solution) that is proportional to the norm of the inverse of the base
    !matrix and to the accumulated change of the pending rows applied to the
    !solution. Its residual ||A*x-b||/(||A||*||x||+||b||) is bounded by
    !eps*growth and its forward error by eps*growth*cond. Three classes are
    !distinguished:
    !  - random sequences of the Ritz model, where pending functions are
    !    replaced by unrelated ones: growth allowed up to 1e6;
    !  - the Gaussian sequence D, where functions with matrix elements up to
    !    1e6 are replaced by unrelated ones at a shift 1e-5 from the
    !    eigenvalue: growth allowed up to 1e9;
    !  - the near-singular stress cases, where 30 or more functions carry
    !    small accumulated changes (the regime of an optimizer): growth
    !    allowed up to 1e6.
    !The observed values are printed. Inverse iteration is insensitive to
    !this growth because all iterations after the first are done in
    !residual-correction form.
    tol_res_factor = 1.0e6_wp
    tol_fwd_factor = 1.0e6_wp
    tol_lam_factor = 1.0e2_wp
    tol_vec_factor = 1.0e2_wp
    tol_eigres_factor = 100.0_wp
    tol_plain_factor = 1.0_wp
    high_precision = (precision(1.0_wp) > 15)

    call run_error_paths()
#ifdef QRLINALG_TESTING
    call run_factor_tests()
#endif

    steps = 320
    if (precision(1.0_wp) > 20) steps = 100
    if (quick) steps = 40
    !Well-conditioned random problems: the algebra and the bookkeeping.
    scale_decades = 0.5_wp
    dependence = 0.5_wp
    call run_sequence('A ritz n0=120 cap=160 mmax=16', MODEL_RITZ, 120, 160, 16, &
                      steps, 1.0e-3_wp, 1.0e3_wp, 20241)
    !Small border: frequent flushes requested by begin_active.
    call run_sequence('B ritz n0=61 cap=75 mmax=6', MODEL_RITZ, 61, 75, 6, &
                      steps, 1.0e-4_wp, 1.0e5_wp, 20242)
    !Fewer functions than processes when run with 14 ranks.
    call run_sequence('C ritz n0=9 cap=13 mmax=4', MODEL_RITZ, 9, 13, 4, &
                      steps, 1.0e-4_wp, 1.0e2_wp, 20243)
    !Hydrogen atom in s-type Gaussians: norms up to 1e6, nearly dependent
    !basis, shift 1e-5 below the ground state. The smallest ratio of two
    !exponents is chosen such that eps*cond(A) stays below about 1e-2; beyond
    !that even a fresh Householder QR factorization of this strongly
    !row-graded matrix is too inaccurate for inverse iteration to reach the
    !roundoff floor of the residual (observed with ratio 1.15 in extended
    !precision: cond(A) = 2e17, residual 100 times the floor on fresh
    !factors, with and without a border).
    min_ratio = 1.3_wp
    if (high_precision) min_ratio = 1.2_wp
    tol_res_factor = 1.0e9_wp
    tol_fwd_factor = 1.0e9_wp
    call run_sequence('D gauss n0=24 cap=32 mmax=8', MODEL_GAUSS, 24, 32, 8, &
                      steps, 1.0e-5_wp, 0.0_wp, 20244)
    tol_res_factor = 1.0e6_wp
    tol_fwd_factor = 1.0e6_wp
    !Larger problem with badly scaled, nearly dependent functions.
    scale_decades = 1.0_wp
    dependence = 0.1_wp
    if (precision(1.0_wp) <= 20 .and. .not. quick) then
      call run_sequence('E ritz n0=257 cap=300 mmax=40', MODEL_RITZ, 257, 300, 40, &
                        160, 1.0e-5_wp, 1.0e4_wp, 20245)
    else
      call run_sequence('E ritz n0=101 cap=120 mmax=24', MODEL_RITZ, 101, 120, 24, &
                        40, 1.0e-5_wp, 1.0e4_wp, 20245)
    end if

    !Near-singular stress cases: shift 1e-6 (relative) below the lowest
    !eigenvalue, norms of the shifted matrix about 1e7, nearly dependent
    !basis, 30 or more pending functions without any flush.
    scale_decades = 1.5_wp
    dependence = 0.05_wp
    if (high_precision) dependence = 0.01_wp
    if (precision(1.0_wp) <= 20 .and. .not. quick) then
      call run_stress('S1 ritz n=300 shift 1e-6 pending 36', MODEL_RITZ, 300, 48, 36, &
                      1.0e-6_wp, 1.0e4_wp, 20246)
    else
      call run_stress('S1 ritz n=120 shift 1e-6 pending 33', MODEL_RITZ, 120, 40, 33, &
                      1.0e-6_wp, 1.0e4_wp, 20246)
    end if
    min_ratio = 1.3_wp
    if (high_precision) min_ratio = 1.15_wp
    call run_stress('S2 gauss n=40 shift 1e-6 pending 30', MODEL_GAUSS, 40, 36, 30, &
                    1.0e-6_wp, 0.0_wp, 20247)
    call problem_free()
  end subroutine run_unit_tests

  !=========================================================================
  ! Timing mode
  !=========================================================================

  subroutine timing_matrices(order)
    !Synthetic dense matrices for the timing mode: S = I + small random
    !symmetric part (positive definite); H = diagonal of positive energies
    !(a few of them large) + small random symmetric part - 8*q*q^T with a
    !dense vector q, which creates one well separated bound state whose
    !eigenvector involves all functions.
    integer, intent(in) :: order
    integer :: i, j
    real(wp) :: rs
    rs = 1.0_wp / sqrt(real(order, wp))
    do i = 1, order
      wrk1(i) = rs * (0.5_wp + rnd())
      if (rnd() < 0.5_wp) wrk1(i) = -wrk1(i)
      wrk2(i) = 0.5_wp + 4.5_wp * rnd()
      if (mod(i, 50) == 0) wrk2(i) = 10.0_wp**(1.0_wp + 5.0_wp * rnd())
    end do
    do j = 1, order
      do i = j, order
        smat(i, j) = 0.2_wp * rs * (2.0_wp * rnd() - 1.0_wp)
        hmat(i, j) = 0.3_wp * rs * (2.0_wp * rnd() - 1.0_wp) - 8.0_wp * wrk1(i) * wrk1(j)
      end do
      smat(j, j) = smat(j, j) + 1.0_wp
      hmat(j, j) = hmat(j, j) + wrk2(j)
    end do
    n = order
  end subroutine timing_matrices

  subroutine timing_perturb(k, amount)
    !New columns hcol, scol of function k: the current ones with a relative
    !random change; the physical lower triangles are updated.
    integer, intent(in) :: k
    real(wp), intent(in) :: amount
    integer :: j
    do j = 1, n
      hcol(j) = hmat(max(j, k), min(j, k)) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
      scol(j) = smat(max(j, k), min(j, k)) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
    end do
    call store_column(k, n)
  end subroutine timing_perturb

  function wall() result(t)
    !Wall-clock time after synchronizing the ranks.
    real(kind(1.0d0)) :: t
    if (nprocs > 1) call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t = MPI_Wtime()
  end function wall

  subroutine time_trials(label, k, ntrials, amount, tol)
    !Trials of the active function k: one replace_symmetric followed by solve
    !started from the previous eigenvector.
    character(len=*), intent(in) :: label
    integer, intent(in) :: k, ntrials
    real(wp), intent(in) :: amount, tol
    integer :: it, info, num_iter, iters, i
    real(wp) :: lambda, rel_acc, eig_res, worst
    real(kind(1.0d0)) :: t0, t1, t2, t_replace, t_solve, t_min
    t_replace = 0.0d0
    t_solve = 0.0d0
    t_min = huge(1.0d0)
    iters = 0
    worst = 0.0_wp
    lambda = 0.0_wp
    do it = 1, ntrials
      call timing_perturb(k, amount)
      t0 = wall()
      call st%replace_symmetric(k, hcol(1:n), scol(1:n), info)
      t1 = wall()
      if (info /= QR_SUCCESS) call fail(label // ': replace_symmetric')
      do i = 1, n
        wrk1(i) = xvec(i)
      end do
      call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, tol, 50, 0, rel_acc, &
                    num_iter, eig_res, info)
      t2 = wall()
      if (info /= QR_SUCCESS) call fail(label // ': solve')
      t_replace = t_replace + (t1 - t0)
      t_solve = t_solve + (t2 - t1)
      t_min = min(t_min, t2 - t0)
      iters = iters + num_iter
      worst = max(worst, eig_res)
    end do
    lam_ref = lambda
    if (myrank == 0) then
      write(*, '(a,a,i0,a,i0)') trim(label), ': pending=', st%pending_count(), &
        ' excluded=', st%excluded_count()
      write(*, '(a,f9.3,a,f9.3,a,f9.3,a,f9.3,a)') '   trial mean ', &
        1.0d3 * (t_replace + t_solve) / ntrials, ' ms (replace ', &
        1.0d3 * t_replace / ntrials, ' ms + solve ', 1.0d3 * t_solve / ntrials, &
        ' ms), fastest trial ', 1.0d3 * t_min, ' ms'
      write(*, '(a,f6.2,a,es9.2,a,es24.16)') '   mean iterations ', &
        real(iters) / real(ntrials), ' worst eigen-residual ', worst, ' lambda ', lambda
    end if
  end subroutine time_trials

  subroutine run_timing(order, npending, ntrials)
    !Wall times of the building blocks and of full trials for a synthetic
    !problem of the given order.
    integer, intent(in) :: order, npending, ntrials
    integer :: info, num_iter, i, j, k, rep, set(1), stride
    logical :: need_flush
    real(wp) :: lambda, rel_acc, eig_res, tol, lam_border, abs_res, rel_res, fwd
    real(kind(1.0d0)) :: t0, t1, t_sum, t_begin
    integer, parameter :: NREP = 10

    tol = 1.0e-12_wp
    rng_state = 987654321_int64
    model = MODEL_RITZ
    cap = order
    allocate(hmat(order, order), smat(order, order), hcol(order), scol(order))
    allocate(bvec(order), xvec(order), xref(order), wrk1(order), wrk2(order))
    do j = 1, order
      do i = 1, order
        hmat(i, j) = 7.7e33_wp
        smat(i, j) = -3.3e31_wp
      end do
    end do
    call timing_matrices(order)
    call st%initialize(order, max(npending + 8, 40), info)
    if (info /= QR_SUCCESS) call fail('timing: initialize')
    if (myrank == 0) write(*, '(a,i0,a,i0,a,i0)') 'timing: n=', order, ' max_pending=', &
      st%max_pending(), ' trials=', ntrials

    !Rough shift below the spectrum, then a shift 1e-6 (relative) below the
    !lowest eigenvalue, as the application chooses it.
    t0 = wall()
    call st%factorize_fresh(hmat, smat, -20.0_wp, info)
    t1 = wall()
    if (info /= QR_SUCCESS) call fail('timing: factorize_fresh')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'factorize_fresh (far shift)      ', t1 - t0, ' s'
    wrk1(1:n) = 1.0_wp
    call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, 1.0e-10_wp, 500, 0, rel_acc, &
                  num_iter, eig_res, info)
    if (info /= QR_SUCCESS) call fail('timing: first solve')
    if (myrank == 0) write(*, '(a,es24.16,a,i0,a,es9.2)') '   lowest eigenvalue ', lambda, &
      ' iterations ', num_iter, ' residual ', eig_res
    shift = lambda - 1.0e-6_wp * abs(lambda)
    t0 = wall()
    call st%factorize_fresh(hmat, smat, shift, info)
    t1 = wall()
    if (info /= QR_SUCCESS) call fail('timing: factorize_fresh')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'factorize_fresh (near shift)     ', t1 - t0, ' s'

    !Base solve and symmetric product.
    do i = 1, n
      bvec(i) = 2.0_wp * rnd() - 1.0_wp
    end do
    t_sum = 0.0d0
    do rep = 1, NREP
      t0 = wall()
      call st%solve_shifted(bvec(1:n), xref(1:n), info)
      t1 = wall()
      t_sum = t_sum + (t1 - t0)
    end do
    if (myrank == 0) write(*, '(a,f10.3,a)') 'base solve (no border)           ', &
      1.0d3 * t_sum / NREP, ' ms'
    t_sum = 0.0d0
    do rep = 1, NREP
      t0 = wall()
      call st%symmetric_product(smat, bvec(1:n), xref(1:n), info)
      t1 = wall()
      t_sum = t_sum + (t1 - t0)
    end do
    if (myrank == 0) write(*, '(a,f10.3,a)') 'symmetric product                ', &
      1.0d3 * t_sum / NREP, ' ms'
    call st%factorization_residual(hmat, smat, bvec(1:n), abs_res, rel_res, info, fwd)
    if (myrank == 0) write(*, '(a,es9.2,a,es9.2)') '   factorization_residual ', rel_res, &
      ' forward error ', fwd

    !Eigenvector with the near shift as the start of the trials.
    wrk1(1:n) = xvec(1:n)
    t0 = wall()
    call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, tol, 50, 0, rel_acc, &
                  num_iter, eig_res, info)
    t1 = wall()
    if (myrank == 0) write(*, '(a,f10.3,a,i0,a,es9.2)') 'solve after fresh factorization  ', &
      1.0d3 * (t1 - t0), ' ms, iterations ', num_iter, ' residual ', eig_res

    !Trials with one pending function.
    k = n / 2
    set(1) = k
    t0 = wall()
    call st%begin_active(set, need_flush, info)
    t1 = wall()
    if (info /= QR_SUCCESS .or. need_flush) call fail('timing: begin_active')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'begin_active (new base function) ', &
      1.0d3 * (t1 - t0), ' ms'
    call time_trials('trials, 1 pending function', k, ntrials, 1.0e-3_wp, tol)

    !Make npending-1 further functions pending, one optimized function at a
    !time, each with a few accepted changes.
    stride = max(1, n / (npending + 1))
    t_begin = 0.0d0
    do j = 1, npending - 1
      k = 1 + mod(j * stride, n)
      if (k == n / 2) k = k + 1
      set(1) = k
      t0 = wall()
      call st%begin_active(set, need_flush, info)
      t1 = wall()
      t_begin = t_begin + (t1 - t0)
      if (info /= QR_SUCCESS .or. need_flush) call fail('timing: begin_active (pending)')
      call timing_perturb(k, 1.0e-3_wp)
      call st%replace_symmetric(k, hcol(1:n), scol(1:n), info)
      if (info /= QR_SUCCESS) call fail('timing: replace (pending)')
    end do
    if (npending > 1) then
      if (myrank == 0) write(*, '(a,f10.3,a)') 'begin_active (mean while filling)', &
        1.0d3 * t_begin / (npending - 1), ' ms'
      !Bring the eigenvector up to date, then time trials of the last one.
      wrk1(1:n) = xvec(1:n)
      call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, tol, 100, 0, rel_acc, &
                    num_iter, eig_res, info)
      if (myrank == 0) write(*, '(a,i0,a,es9.2)') '   solve after filling: iterations ', &
        num_iter, ' residual ', eig_res
      call time_trials('trials, full border', k, ntrials, 1.0e-3_wp, tol)
      t_sum = 0.0d0
      do rep = 1, NREP
        t0 = wall()
        call st%solve_shifted(bvec(1:n), xref(1:n), info)
        t1 = wall()
        t_sum = t_sum + (t1 - t0)
      end do
      if (myrank == 0) write(*, '(a,f10.3,a)') 'bordered solve (full border)     ', &
        1.0d3 * t_sum / NREP, ' ms'
    end if

    !Agreement of the bordered representation with a fresh factorization.
    lam_border = lam_ref
    call st%factorize_fresh(hmat, smat, shift, info)
    wrk1(1:n) = xvec(1:n)
    call st%solve(hmat, smat, wrk1(1:n), xref(1:n), lambda, tol, 100, 0, rel_acc, &
                  num_iter, eig_res, info)
    if (myrank == 0) then
      write(*, '(a,es24.16,a,es9.2)') 'fresh factorization of the final matrices: lambda ', &
        lambda, ' residual ', eig_res
      write(*, '(a,es9.2,a,es9.2)') '   |lambda(border)-lambda(fresh)|/|lambda| ', &
        abs(lam_border - lambda) / abs(lambda), '  max|x(border)-x(fresh)| ', &
        maxval(abs(xvec(1:n) - xref(1:n)))
      if (nfail == 0) then
        write(*, '(a)') 'TIMING RUN COMPLETED'
      else
        write(*, '(a,i0)') 'TIMING RUN FAILED: ', nfail
      end if
    end if
    call st%clear()
  end subroutine run_timing

  !=========================================================================
  !=========================================================================
  ! Complex Hermitian suite (type qr_complex_state)
  !
  ! The procedures below mirror the real ones above; their names and the
  ! names of their data carry the prefix z. The physical matrices follow the
  ! convention of CG_0S: the stored element (i,j), i>=j, of the lower triangle
  ! is <phi_i|O|phi_j> (the row function is the bra), the conceptual element
  ! (i,j) with i<j is the complex conjugate of the stored element (j,i), and
  ! the columns zhcol, zscol passed to the module are conceptual columns:
  ! element j of the column of function k is <phi_j|O|phi_k>. The upper
  ! triangles hold garbage. The dense reference builds the full Hermitian
  ! matrix from the lower triangles independently of the module.
  !=========================================================================
  !=========================================================================

  function zstart(i) result(z)
    !Element i of the fixed complex starting vector of the eigen-solves.
    integer, intent(in) :: i
    complex(wp) :: z
    z = cmplx(1.0_wp, 0.25_wp * real(1 - 2 * mod(i, 2), wp), kind=wp)
  end function zstart

  function zrnd() result(z)
    !Complex pseudo-random number with real and imaginary parts in (-1,1).
    complex(wp) :: z
    real(wp) :: re, im
    re = 2.0_wp * rnd() - 1.0_wp
    im = 2.0_wp * rnd() - 1.0_wp
    z = cmplx(re, im, kind=wp)
  end function zrnd

  subroutine zexpect_replicated(k, x, message)
    !Checks that a complex vector is bitwise identical on all ranks.
    integer, intent(in) :: k
    complex(wp), intent(in) :: x(k)
    character(len=*), intent(in) :: message
    integer :: i, nbad
    if (nprocs == 1) return
    nchecks = nchecks + 1
    zbcast_buf(1:k) = x(1:k)
    call MPI_Bcast(zbcast_buf, 2 * k, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    nbad = 0
    do i = 1, k
      if (.not. (zbcast_buf(i) == x(i))) nbad = nbad + 1
    end do
    call MPI_Allreduce(MPI_IN_PLACE, nbad, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (nbad /= 0) call fail(trim(message) // ': result differs between ranks')
  end subroutine zexpect_replicated

  !=========================================================================
  ! Complex test problems
  !=========================================================================

  subroutine zproblem_allocate(model_in, capacity, modes)
    integer, intent(in) :: model_in, capacity, modes
    integer :: i, j
    call zproblem_free()
    model = model_in
    cap = capacity
    dmodes = modes
    allocate(zhmat(cap, cap), zsmat(cap, cap), zhcol(cap), zscol(cap))
    allocate(zaref(cap, cap), zalu(cap, cap), zapiv(cap))
    allocate(zbvec(cap), zxvec(cap), zxref(cap), zwrk1(cap), zwrk2(cap), zwrk3(cap))
    allocate(zwrk4(cap), zwrk5(cap), zrw4(cap), zrw5(cap))
    allocate(zbcast_buf(cap))
    allocate(zexpo(cap), zudelta(cap), zuener(cap))
    allocate(zuvec(max(1, dmodes), cap), zdval(max(1, dmodes)))
    zudelta = 0.0_wp
    zuener = 0.0_wp
    !Complex garbage everywhere; the lower triangles are overwritten by the
    !functions.
    do j = 1, cap
      do i = 1, cap
        zhmat(i, j) = cmplx(7.7e33_wp * real(1 - 2 * mod(i + j, 2), wp), &
                            -5.5e32_wp * real(1 - 2 * mod(j, 2), wp), kind=wp)
        zsmat(i, j) = cmplx(-3.3e31_wp * real(1 - 2 * mod(i, 2), wp), 9.9e30_wp, kind=wp)
      end do
    end do
    zexpo = 0.0_wp
    zuvec = 0.0_wp
    zdval = 0.0_wp
  end subroutine zproblem_allocate

  subroutine zproblem_free()
    if (allocated(zhmat)) deallocate(zhmat, zsmat, zhcol, zscol, zaref, zalu, zapiv, &
                                     zbvec, zxvec, zxref, zwrk1, zwrk2, zwrk3, &
                                     zwrk4, zwrk5, zrw4, zrw5, zbcast_buf, zexpo, &
                                     zuvec, zdval, zudelta, zuener)
  end subroutine zproblem_free

  subroutine znew_function(k, norder)
    !Draws new parameters for function k; norder is the number of functions
    !considered to exist.
    !
    !Ritz model: as in the real suite, with complex vectors u_k,
    !    S_ij = u_i^H.u_j + delta_ij*udelta(i)**2,
    !    H_ij = u_i^H.(dval*u_j) + delta_ij*udelta(i)**2*uener(i),
    !where the component on the bound mode has a random phase and the other
    !components have independent real and imaginary parts. H and S are
    !Hermitian with genuinely complex off-diagonal elements, and S is positive
    !definite.
    !
    !Gaussian model: hydrogen atom in complex s-type Gaussians exp(-a*r**2)
    !with a = alpha*(1+i*t): alpha is distributed log-uniformly over eight
    !decades, not closer to the real part of an existing exponent than the
    !ratio min_ratio, and t is uniform in (-0.6,0.6).
    integer, intent(in) :: k, norder
    integer :: i, j, attempt
    real(wp) :: scale, alpha, theta, re, im
    logical :: too_close
    if (model == MODEL_GAUSS) then
      alpha = 1.0_wp
      do attempt = 1, 200
        alpha = 10.0_wp**(-2.0_wp + 8.0_wp * rnd())
        too_close = .false.
        do j = 1, norder
          if (j == k) cycle
          if (abs(log(alpha / real(zexpo(j), wp))) < log(min_ratio)) too_close = .true.
        end do
        if (.not. too_close) exit
      end do
      zexpo(k) = cmplx(alpha, alpha * 0.6_wp * (2.0_wp * rnd() - 1.0_wp), kind=wp)
    else
      scale = 10.0_wp**(scale_decades * (2.0_wp * rnd() - 1.0_wp))
      theta = 8.0_wp * atan(1.0_wp) * rnd()
      zuvec(1, k) = scale * (0.3_wp + 0.7_wp * rnd()) * cmplx(cos(theta), sin(theta), kind=wp)
      do i = 2, dmodes
        re = 2.0_wp * rnd() - 1.0_wp
        im = 2.0_wp * rnd() - 1.0_wp
        zuvec(i, k) = scale * cmplx(re, im, kind=wp) / sqrt(2.0_wp * zdval(i))
      end do
      zudelta(k) = scale * dependence * (0.5_wp + rnd())
      zuener(k) = 2.0_wp + 10.0_wp * rnd()
    end if
  end subroutine znew_function

  subroutine zperturb_function(k, amount)
    !Small relative change of the parameters of function k (optimizer step).
    integer, intent(in) :: k
    real(wp), intent(in) :: amount
    integer :: i
    real(wp) :: re, im
    if (model == MODEL_GAUSS) then
      re = real(zexpo(k), wp) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
      im = aimag(zexpo(k)) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
      zexpo(k) = cmplx(re, im, kind=wp)
    else
      do i = 1, dmodes
        re = 1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp)
        im = amount * (2.0_wp * rnd() - 1.0_wp)
        zuvec(i, k) = zuvec(i, k) * cmplx(re, im, kind=wp)
      end do
      zudelta(k) = zudelta(k) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
    end if
  end subroutine zperturb_function

  subroutine zfunction_column(k, norder)
    !Computes the full conceptual columns zhcol(1:norder), zscol(1:norder) of
    !function k: element j is <phi_j|O|phi_k>, the function j being the bra.
    integer, intent(in) :: k, norder
    integer :: i, j
    real(wp) :: pi, nn
    complex(wp) :: ak, aj, p, hs, ss, sjk
    if (model == MODEL_GAUSS) then
      !Normalized complex s-type Gaussians exp(-a*r**2) for the hydrogen atom.
      !With p = conjg(aj)+ak and N = (4*Re(aj)*Re(ak)/pi**2)**0.75:
      !S = N*(pi/p)**1.5, T = 3*conjg(aj)*ak/p*S, V = -2*pi/p*N.
      pi = 4.0_wp * atan(1.0_wp)
      ak = zexpo(k)
      do j = 1, norder
        aj = zexpo(j)
        p = conjg(aj) + ak
        nn = (4.0_wp * real(aj, wp) * real(ak, wp) / pi**2)**0.75_wp
        sjk = nn * (pi / p)**1.5_wp
        if (j == k) sjk = 1.0_wp
        zscol(j) = sjk
        zhcol(j) = 3.0_wp * conjg(aj) * ak / p * sjk - 2.0_wp * pi / p * nn
      end do
    else
      do j = 1, norder
        hs = 0.0_wp
        ss = 0.0_wp
        do i = 1, dmodes
          ss = ss + conjg(zuvec(i, j)) * zuvec(i, k)
          hs = hs + conjg(zuvec(i, j)) * zdval(i) * zuvec(i, k)
        end do
        if (j == k) then
          ss = ss + zudelta(k)**2
          hs = hs + zudelta(k)**2 * zuener(k)
        end if
        zhcol(j) = hs
        zscol(j) = ss
      end do
    end if
    !The diagonal elements are real.
    zhcol(k) = cmplx(real(zhcol(k), wp), 0.0_wp, kind=wp)
    zscol(k) = cmplx(real(zscol(k), wp), 0.0_wp, kind=wp)
  end subroutine zfunction_column

  subroutine zstore_column(k, norder)
    !Writes the conceptual columns zhcol, zscol of function k into row and
    !column k of the physical lower triangles: element j goes to the stored
    !element (j,k) when j >= k, and its complex conjugate to the stored
    !element (k,j) when j < k.
    integer, intent(in) :: k, norder
    integer :: j
    do j = 1, norder
      if (j >= k) then
        zhmat(j, k) = zhcol(j)
        zsmat(j, k) = zscol(j)
      else
        zhmat(k, j) = conjg(zhcol(j))
        zsmat(k, j) = conjg(zscol(j))
      end if
    end do
  end subroutine zstore_column

  subroutine zproblem_build(norder, dmax)
    !Creates norder functions and the physical matrices.
    integer, intent(in) :: norder
    real(wp), intent(in) :: dmax
    integer :: k
    if (model == MODEL_RITZ) then
      zdval(1) = -1.0_wp
      do k = 2, dmodes
        zdval(k) = 0.5_wp * (dmax / 0.5_wp)**(real(k - 2, wp) / real(max(1, dmodes - 2), wp))
      end do
    end if
    do k = 1, norder
      call znew_function(k, k - 1)
    end do
    do k = 1, norder
      call zfunction_column(k, k)
      call zstore_column(k, k)
    end do
    n = norder
  end subroutine zproblem_build

  !=========================================================================
  ! Dense complex reference
  !=========================================================================

  subroutine zref_symv(norder, a, x, y)
    !y = A*x for the Hermitian matrix whose lower triangle is stored in a.
    integer, intent(in) :: norder
    complex(wp), intent(in) :: a(:,:), x(:)
    complex(wp), intent(out) :: y(:)
    integer :: i, j
    y(1:norder) = 0.0_wp
    do j = 1, norder
      y(j) = y(j) + real(a(j, j), wp) * x(j)
      do i = j + 1, norder
        y(i) = y(i) + a(i, j) * x(j)
        y(j) = y(j) + conjg(a(i, j)) * x(i)
      end do
    end do
  end subroutine zref_symv

  subroutine zref_build()
    !Full current shifted Hermitian matrix, its LU factorization with partial
    !pivoting, its Frobenius norm, and estimates of the two largest singular
    !values of its inverse (see ref_build).
    integer :: i, j, k, ip, it
    real(wp) :: pmax, znorm, m11, m22, tr
    complex(wp) :: t, m12
    anorm = 0.0_wp
    do j = 1, n
      zaref(j, j) = cmplx(real(zhmat(j, j), wp) - shift * real(zsmat(j, j), wp), 0.0_wp, kind=wp)
      do i = j + 1, n
        zaref(i, j) = zhmat(i, j) - shift * zsmat(i, j)
        zaref(j, i) = conjg(zaref(i, j))
      end do
    end do
    do j = 1, n
      do i = 1, n
        zalu(i, j) = zaref(i, j)
        anorm = anorm + real(zaref(i, j), wp)**2 + aimag(zaref(i, j))**2
      end do
    end do
    anorm = sqrt(anorm)
    do k = 1, n
      ip = k
      pmax = abs(zalu(k, k))
      do i = k + 1, n
        if (abs(zalu(i, k)) > pmax) then
          ip = i
          pmax = abs(zalu(i, k))
        end if
      end do
      zapiv(k) = ip
      if (ip /= k) then
        do j = 1, n
          t = zalu(k, j)
          zalu(k, j) = zalu(ip, j)
          zalu(ip, j) = t
        end do
      end if
      if (pmax > 0.0_wp) then
        do i = k + 1, n
          zalu(i, k) = zalu(i, k) / zalu(k, k)
        end do
        do j = k + 1, n
          t = zalu(k, j)
          do i = k + 1, n
            zalu(i, j) = zalu(i, j) - zalu(i, k) * t
          end do
        end do
      end if
    end do
    !Two-vector subspace iteration on the inverse (the matrix is Hermitian).
    do i = 1, n
      zwrk1(i) = cmplx(1.0_wp + 0.37_wp * real(mod(7 * i, 11), wp), &
                       0.3_wp * real(mod(3 * i, 7), wp) - 0.8_wp, kind=wp)
      zwrk2(i) = cmplx(1.0_wp + 0.21_wp * real(mod(5 * i, 13), wp), 0.4_wp, kind=wp) * &
                 real(1 - 2 * mod(i, 2), wp)
    end do
    ainv_norm = 0.0_wp
    ainv2_norm = 0.0_wp
    do it = 1, 6
      znorm = sqrt(sum(abs(zwrk1(1:n))**2))
      zwrk1(1:n) = zwrk1(1:n) / znorm
      if (n > 1) then
        t = sum(conjg(zwrk1(1:n)) * zwrk2(1:n))
        zwrk2(1:n) = zwrk2(1:n) - t * zwrk1(1:n)
        znorm = sqrt(sum(abs(zwrk2(1:n))**2))
        zwrk2(1:n) = zwrk2(1:n) / znorm
      else
        zwrk2(1) = 0.0_wp
      end if
      zwrk4(1:n) = zwrk1(1:n)
      zwrk5(1:n) = zwrk2(1:n)
      call zref_solve(zwrk1)
      call zref_solve(zwrk2)
      m11 = real(sum(conjg(zwrk4(1:n)) * zwrk1(1:n)), wp)
      m12 = sum(conjg(zwrk4(1:n)) * zwrk2(1:n))
      m22 = real(sum(conjg(zwrk5(1:n)) * zwrk2(1:n)), wp)
      tr = sqrt((0.5_wp * (m11 - m22))**2 + abs(m12)**2)
      ainv_norm = max(abs(0.5_wp * (m11 + m22) + tr), abs(0.5_wp * (m11 + m22) - tr))
      ainv2_norm = min(abs(0.5_wp * (m11 + m22) + tr), abs(0.5_wp * (m11 + m22) - tr))
    end do
    if (n == 1) ainv2_norm = ainv_norm
  end subroutine zref_build

  subroutine zref_solve(b)
    !Solves the current shifted system in place with the reference LU.
    complex(wp), intent(inout) :: b(:)
    integer :: i, k
    complex(wp) :: t
    do k = 1, n
      if (zapiv(k) /= k) then
        t = b(k)
        b(k) = b(zapiv(k))
        b(zapiv(k)) = t
      end if
    end do
    do k = 1, n
      t = b(k)
      do i = k + 1, n
        b(i) = b(i) - zalu(i, k) * t
      end do
    end do
    do k = n, 1, -1
      b(k) = b(k) / zalu(k, k)
      t = b(k)
      do i = 1, k - 1
        b(i) = b(i) - zalu(i, k) * t
      end do
    end do
  end subroutine zref_solve

  subroutine zref_abs_symv(norder, a, x, y)
    !y = |A|*|x| (moduli) using the lower triangle of a only.
    integer, intent(in) :: norder
    complex(wp), intent(in) :: a(:,:), x(:)
    real(wp), intent(out) :: y(:)
    integer :: i, j
    y(1:norder) = 0.0_wp
    do j = 1, norder
      y(j) = y(j) + abs(real(a(j, j), wp)) * abs(x(j))
      do i = j + 1, norder
        y(i) = y(i) + abs(a(i, j)) * abs(x(j))
        y(j) = y(j) + abs(a(i, j)) * abs(x(i))
      end do
    end do
  end subroutine zref_abs_symv

  subroutine zref_eigen(x, lambda)
    !Dense reference inverse iteration with the current shift; S-normalized
    !eigenvector (x^H S x = 1) whose component of largest modulus is real and
    !positive, physical Rayleigh quotient. Sets ref_res, ref_floor,
    !ref_absform and ref_sabs as ref_eigen does, with moduli in place of
    !absolute values.
    complex(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    integer :: i, it
    real(wp) :: change, change_prev, xs, vnorm2, hn, sn
    complex(wp) :: alpha
    do i = 1, n
      zwrk1(i) = zstart(i)
    end do
    change_prev = huge(1.0_wp)
    ref_converged = .false.
    do it = 1, 120
      call zref_symv(n, zsmat, zwrk1, x)
      call zref_solve(x)
      x(1:n) = x(1:n) / maxval(abs(x(1:n)))
      vnorm2 = sum(abs(zwrk1(1:n))**2)
      alpha = sum(conjg(zwrk1(1:n)) * x(1:n)) / vnorm2
      change = sqrt(sum(abs(x(1:n) - alpha * zwrk1(1:n))**2) / sum(abs(x(1:n))**2))
      zwrk1(1:n) = x(1:n)
      if ((change >= change_prev .and. change < 1.0e-3_wp) .or. change == 0.0_wp) then
        ref_converged = .true.
        exit
      end if
      change_prev = change
    end do
    call zref_symv(n, zsmat, x, zwrk3)
    xs = real(sum(conjg(x(1:n)) * zwrk3(1:n)), wp)
    x(1:n) = x(1:n) / sqrt(xs)
    call zfix_phase(x)
    call zref_symv(n, zhmat, x, zwrk2)
    call zref_symv(n, zsmat, x, zwrk3)
    lambda = real(sum(conjg(x(1:n)) * zwrk2(1:n)), wp) / &
             real(sum(conjg(x(1:n)) * zwrk3(1:n)), wp)
    hn = sqrt(sum(abs(zwrk2(1:n))**2))
    sn = sqrt(sum(abs(zwrk3(1:n))**2))
    ref_res = sqrt(sum(abs(zwrk2(1:n) - lambda * zwrk3(1:n))**2)) / (hn + abs(lambda) * sn)
    call zref_abs_symv(n, zhmat, x, zrw4)
    call zref_abs_symv(n, zsmat, x, zrw5)
    ref_floor = EPS * (sqrt(sum(zrw4(1:n)**2)) + abs(lambda) * sqrt(sum(zrw5(1:n)**2))) / &
                (hn + abs(lambda) * sn)
    ref_sabs = sum(abs(x(1:n)) * zrw5(1:n))
    ref_absform = sum(abs(x(1:n)) * zrw4(1:n)) + abs(lambda) * ref_sabs
  end subroutine zref_eigen

  subroutine zfix_phase(x)
    !Makes the component of largest modulus real and positive.
    complex(wp), intent(inout) :: x(:)
    integer :: i, imax
    imax = 1
    do i = 2, n
      if (abs(x(i)) > abs(x(imax))) imax = i
    end do
    if (abs(x(imax)) > 0.0_wp) x(1:n) = x(1:n) * (conjg(x(imax)) / abs(x(imax)))
  end subroutine zfix_phase

  subroutine zalign_phase(x, ref)
    !Removes the arbitrary complex phase of an eigenvector: x is multiplied by
    !the unit factor that makes ref^H x real and positive.
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(in) :: ref(:)
    complex(wp) :: c
    c = sum(conjg(ref(1:n)) * x(1:n))
    if (abs(c) > 0.0_wp) x(1:n) = x(1:n) * (conjg(c) / abs(c))
  end subroutine zalign_phase

  !=========================================================================
  ! Complex checks against the reference
  !=========================================================================

  subroutine zcheck_linear(label)
    !Bordered solve with a random complex right-hand side against the dense
    !reference solve of the current matrix.
    character(len=*), intent(in) :: label
    integer :: i, info
    real(wp) :: res, fwd, xnorm, bnorm, cond, tol_fwd
    character(len=200) :: text

    stat_max_pending = max(stat_max_pending, zst%pending_count())
    stat_max_excluded = max(stat_max_excluded, zst%excluded_count())
    if (zst%order() < zst%base_order()) stat_below_base = stat_below_base + 1
    call expect(zst%order() == n, trim(label) // ': order() differs from the test order')

    do i = 1, n
      zbvec(i) = zrnd()
    end do
    call zst%solve_shifted(zbvec(1:n), zxvec(1:n), info)
    if (info == QR_ERR_SINGULAR .and. allow_flush .and. zst%pending_count() > 0) then
      !The committed block failed the pivot test: flush, as the application.
      call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': solve_shifted')
      stat_singular_flush = stat_singular_flush + 1
      call zdo_flush(0.0_wp, .false.)
      call zst%solve_shifted(zbvec(1:n), zxvec(1:n), info)
    end if
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve_shifted')
    if (info /= QR_SUCCESS) return
    call zexpect_replicated(n, zxvec, trim(label) // ': solve_shifted')
    zxref(1:n) = zbvec(1:n)
    call zref_solve(zxref)

    xnorm = sqrt(sum(abs(zxvec(1:n))**2))
    bnorm = sqrt(sum(abs(zbvec(1:n))**2))
    zwrk1(1:n) = matmul(zaref(1:n, 1:n), zxvec(1:n)) - zbvec(1:n)
    res = sqrt(sum(abs(zwrk1(1:n))**2)) / (anorm * xnorm + bnorm)
    fwd = sqrt(sum(abs(zxvec(1:n) - zxref(1:n))**2)) / sqrt(sum(abs(zxref(1:n))**2))
    cond = anorm * ainv_norm
    tol_fwd = tol_fwd_factor * EPS * cond
    worst_res = max(worst_res, res / EPS)
    worst_fwd = max(worst_fwd, fwd / (EPS * cond))
    worst_cond = max(worst_cond, cond)
    stat_solves = stat_solves + 1
    nchecks = nchecks + 2
    if (.not. (res <= tol_res_factor * EPS)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': solve residual/eps=', &
        res / EPS, ' cond=', cond
      call fail(text)
    end if
    if (.not. (fwd <= tol_fwd)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': solve forward error/(eps*cond)=', &
        fwd / (EPS * cond), ' cond=', cond
      call fail(text)
    end if
    !Trace: a norm-like functional of the solution. Its sensitivity to
    !roundoff is bounded by the forward error bound.
    call add_trace(xnorm, tol_fwd * xnorm)
  end subroutine zcheck_linear

  subroutine zcheck_eigen(label)
    !Eigenpair of solve against the dense reference inverse iteration (the
    !arbitrary complex phase of the eigenvector is removed by aligning it with
    !the reference vector). The
    !plain solve (no residual correction) is checked against loose bounds and
    !its errors are reported; the solve with residual correction is checked
    !against bounds derived from the conditioning of the problem:
    !  eigenvalue   |lambda-ref| <= tol_lam_factor*eps*ref_absform
    !               (roundoff of the Rayleigh quotient of an S-normalized x),
    !  eigenvector  ||x-ref||/||ref|| <= tol_vec_factor*eps*cond2, where cond2
    !               is the condition number of the shifted matrix on the
    !               complement of its near-null vector,
    !  residual     eig_residual <= tol_eigres_factor*max(ref_floor,ref_res),
    !  norm         |x^H S x - 1| <= tol_lam_factor*eps*ref_sabs.
    character(len=*), intent(in) :: label
    integer :: i, info, num_iter, npol
    real(wp) :: lambda, rel_acc, eig_res, dl, dv, xsx, cond2, res_scale
    character(len=200) :: text

    call zref_eigen(zxref, lam_ref)
    if (.not. ref_converged .and. allow_flush) then
      !The shift of the last scheduled flush has become stale: the wanted
      !eigenvalue has drifted so far that even the dense reference iteration
      !converges too slowly to define a reference eigenvector. Retarget the
      !shift with a fresh factorization. The decision depends only on the
      !dense reference and is therefore the same for every process count.
      stat_retarget = stat_retarget + 1
      call zdo_flush(current_shift_rel, .true.)
      call zref_eigen(zxref, lam_ref)
    end if
    nchecks = nchecks + 1
    if (.not. ref_converged) then
      call fail(trim(label) // ': the dense reference inverse iteration did not converge')
      return
    end if
    !The product of the Frobenius norm of a matrix and the second largest
    !singular value of its inverse is at least one. The six-step estimate of
    !zref_build can fall below that bound when the shifted matrix is
    !indefinite (stale shift) and the order is small; it is then replaced by
    !the bound.
    cond2 = max(1.0_wp, anorm * ainv2_norm)
    res_scale = max(ref_floor, ref_res)
    worst_cond2 = max(worst_cond2, cond2)
    stat_eigs = stat_eigs + 1

    !Plain solve, exactly the iteration of the design specification.
    do i = 1, n
      zwrk1(i) = zstart(i)
    end do
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                  rel_acc, num_iter, eig_res, info)
    call zexpect_replicated(n, zxvec, trim(label) // ': solve')
    if (allow_flush .and. zst%pending_count() > 0 .and. &
        (info == QR_ERR_NO_CONVERGENCE .or. info == QR_ERR_SINGULAR .or. &
         (info == QR_SUCCESS .and. .not. (eig_res <= tol_eigres_factor * res_scale)))) then
      !Retry policy of the application (DESIGN.md section 5): the solve did
      !not converge or the physical residual is too large while functions
      !are pending. The failure was reported by the module; a fresh
      !factorization with the same shift must give the right eigenpair.
      nchecks = nchecks + 1
      stat_resid_flush = stat_resid_flush + 1
      call zdo_flush(0.0_wp, .false.)
      do i = 1, n
        zwrk1(i) = zstart(i)
      end do
      call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                    rel_acc, num_iter, eig_res, info)
      call zexpect_replicated(n, zxvec, trim(label) // ': solve after flush')
    end if
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve')
    if (info /= QR_SUCCESS .and. info /= QR_ERR_NO_CONVERGENCE) return
    call zalign_phase(zxvec, zxref)
    dl = abs(lambda - lam_ref)
    dv = sqrt(sum(abs(zxvec(1:n) - zxref(1:n))**2)) / sqrt(sum(abs(zxref(1:n))**2))
    plain_lam = max(plain_lam, dl / (EPS * ref_absform))
    plain_vec = max(plain_vec, dv / (EPS * cond2))
    plain_eigres = max(plain_eigres, eig_res / res_scale)
    plain_eigres_abs = max(plain_eigres_abs, eig_res)
    stat_iter_max = max(stat_iter_max, num_iter)
    nchecks = nchecks + 3
    if (.not. (dl <= tol_plain_factor * tol_lam_factor * EPS * ref_absform)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve |lambda-ref|/(eps*absform)=', &
        dl / (EPS * ref_absform)
      call fail(text)
    end if
    if (.not. (dv <= tol_plain_factor * tol_vec_factor * EPS * cond2)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve eigenvector difference/(eps*cond2)=', &
        dv / (EPS * cond2)
      call fail(text)
    end if
    if (.not. (eig_res <= tol_plain_factor * tol_eigres_factor * res_scale)) then
      write(text, '(2a,es10.3)') trim(label), ': plain solve eigenpair residual/floor=', &
        eig_res / res_scale
      call fail(text)
    end if

    !Solve with residual correction.
    do i = 1, n
      zwrk1(i) = zstart(i)
    end do
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, -sqrt(EPS), 300, 0, &
                  rel_acc, num_iter, eig_res, info, 0.0_wp, npol)
    call expect_info(info, QR_SUCCESS, trim(label) // ': solve with correction')
    if (info /= QR_SUCCESS .and. info /= QR_ERR_NO_CONVERGENCE) return
    call zexpect_replicated(n, zxvec, trim(label) // ': solve with correction')
    call zref_symv(n, zsmat, zxvec, zwrk2)
    xsx = real(sum(conjg(zxvec(1:n)) * zwrk2(1:n)), wp)
    call zalign_phase(zxvec, zxref)
    dl = abs(lambda - lam_ref)
    dv = sqrt(sum(abs(zxvec(1:n) - zxref(1:n))**2)) / sqrt(sum(abs(zxref(1:n))**2))
    worst_lam = max(worst_lam, dl / (EPS * ref_absform))
    worst_vec = max(worst_vec, dv / (EPS * cond2))
    worst_eigres = max(worst_eigres, eig_res / res_scale)
    worst_eigres_abs = max(worst_eigres_abs, eig_res)
    worst_norm = max(worst_norm, abs(xsx - 1.0_wp) / (EPS * ref_sabs))
    stat_polish_max = max(stat_polish_max, npol)
    stat_polish_total = stat_polish_total + npol
    nchecks = nchecks + 4
    if (.not. (dl <= tol_lam_factor * EPS * ref_absform)) then
      write(text, '(2a,es10.3,a,es24.16)') trim(label), ': |lambda-ref|/(eps*absform)=', &
        dl / (EPS * ref_absform), ' ref=', lam_ref
      call fail(text)
    end if
    if (.not. (dv <= tol_vec_factor * EPS * cond2)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': eigenvector difference/(eps*cond2)=', &
        dv / (EPS * cond2), ' cond2=', cond2
      call fail(text)
    end if
    if (.not. (eig_res <= tol_eigres_factor * res_scale)) then
      write(text, '(2a,es10.3,a,es10.3)') trim(label), ': eigenpair residual/floor=', &
        eig_res / res_scale, ' floor=', res_scale
      call fail(text)
    end if
    if (.not. (abs(xsx - 1.0_wp) <= tol_lam_factor * EPS * ref_sabs)) then
      write(text, '(2a,es10.3)') trim(label), ': (x^H S x - 1)/(eps*sabs)=', &
        (xsx - 1.0_wp) / (EPS * ref_sabs)
      call fail(text)
    end if
    call add_trace(lambda, 2.0_wp * tol_lam_factor * EPS * ref_absform)
    call zcheck_shifted_quotient(label)
  end subroutine zcheck_eigen

  subroutine zref_shifted_quotient(x, sigma, lambda, bound)
    !Complex counterpart of ref_shifted_quotient: F = x^H*(H-sigma*S)*x and
    !G = x^H*S*x for the full Hermitian matrices rebuilt from the stored lower
    !triangles (element (i,j), i<j, is the conjugate of the stored (j,i); the
    !imaginary parts of the diagonals are ignored), summed row by row in the
    !kind hp with separate real and imaginary parts. The bound uses moduli
    !and is doubled: a complex product carries up to twice the relative error
    !of a real one.
    complex(wp), intent(in) :: x(:)
    real(wp), intent(in) :: sigma
    real(wp), intent(out) :: lambda, bound
    integer :: i, j
    real(hp) :: f, g, fa, ga, row_fr, row_fi, row_gr, row_gi, row_fa, row_ga
    real(hp) :: hr, hi, sr, si, ar, ai, sg, xr, xi, xm

    sg = real(sigma, hp)
    f = 0.0_hp
    g = 0.0_hp
    fa = 0.0_hp
    ga = 0.0_hp
    do i = 1, n
      row_fr = 0.0_hp
      row_fi = 0.0_hp
      row_gr = 0.0_hp
      row_gi = 0.0_hp
      row_fa = 0.0_hp
      row_ga = 0.0_hp
      do j = 1, n
        if (j < i) then
          hr = real(real(zhmat(i, j), wp), hp)
          hi = real(aimag(zhmat(i, j)), hp)
          sr = real(real(zsmat(i, j), wp), hp)
          si = real(aimag(zsmat(i, j)), hp)
        else if (j == i) then
          hr = real(real(zhmat(i, i), wp), hp)
          hi = 0.0_hp
          sr = real(real(zsmat(i, i), wp), hp)
          si = 0.0_hp
        else
          hr = real(real(zhmat(j, i), wp), hp)
          hi = -real(aimag(zhmat(j, i)), hp)
          sr = real(real(zsmat(j, i), wp), hp)
          si = -real(aimag(zsmat(j, i)), hp)
        end if
        xr = real(real(x(j), wp), hp)
        xi = real(aimag(x(j)), hp)
        xm = sqrt(xr * xr + xi * xi)
        ar = hr - sg * sr
        ai = hi - sg * si
        row_fr = row_fr + ar * xr - ai * xi
        row_fi = row_fi + ar * xi + ai * xr
        row_gr = row_gr + sr * xr - si * xi
        row_gi = row_gi + sr * xi + si * xr
        row_fa = row_fa + (sqrt(hr * hr + hi * hi) + abs(sg) * sqrt(sr * sr + si * si)) * xm
        row_ga = row_ga + sqrt(sr * sr + si * si) * xm
      end do
      !Real part of conjg(x(i)) times the row sum.
      xr = real(real(x(i), wp), hp)
      xi = real(aimag(x(i)), hp)
      xm = sqrt(xr * xr + xi * xi)
      f = f + xr * row_fr + xi * row_fi
      g = g + xr * row_gr + xi * row_gi
      fa = fa + xm * row_fa
      ga = ga + xm * row_ga
    end do
    lambda = real(sg + f / g, wp)
    bound = real(2.0_hp * real(n + 8, hp) * real(EPS, hp) * (fa + abs(f / g) * ga) / g + &
                 real(EPS, hp) * abs(sg + f / g), wp)
    !Without a wider kind the reference carries the same rounding errors.
    if (hp == wp) bound = 2.0_wp * bound
  end subroutine zref_shifted_quotient

  subroutine zcheck_shifted_quotient(label)
    !Complex counterpart of check_shifted_quotient: the eigenvector left in
    !zxvec by zcheck_eigen and a fixed complex vector with varying phases that
    !is far from every eigenvector.
    character(len=*), intent(in) :: label
    integer :: i

    call zshifted_quotient_case(label, zxvec, 1)
    do i = 1, n
      zwrk1(i) = cmplx(1.0_wp + 0.5_wp * real(1 - 2 * mod(i, 2), wp) * &
                       real(i, wp) / real(n, wp), &
                       0.75_wp * real(mod(i, 3) - 1, wp) + 0.1_wp, kind=wp)
    end do
    call zshifted_quotient_case(label, zwrk1, 2)
  end subroutine zcheck_shifted_quotient

  subroutine zshifted_quotient_case(label, x, pass)
    !One vector of zcheck_shifted_quotient; pass 1 is the eigenvector.
    character(len=*), intent(in) :: label
    complex(wp), intent(in) :: x(:)
    integer, intent(in) :: pass
    integer :: info
    real(wp) :: lam_shift, lam_hp, bound, dl, sigma
    character(len=200) :: text

    sigma = zst%get_shift()
    call zst%shifted_quotient(zhmat, zsmat, x(1:n), lam_shift, info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': shifted_quotient')
    if (info /= QR_SUCCESS) return
    call expect_same_scalar(lam_shift, trim(label) // ': shifted_quotient')
    call zref_shifted_quotient(x, sigma, lam_hp, bound)
    nchecks = nchecks + 1
    dl = abs(lam_shift - lam_hp)
    worst_shift_ref = max(worst_shift_ref, dl / bound)
    if (.not. (dl <= bound)) then
      write(text, '(2a,i0,a,es10.3,a,es24.16)') trim(label), ': shifted quotient, vector ', pass, &
        ': |lambda-driver|/bound=', dl / bound, ' driver=', lam_hp
      call fail(text)
    end if
    if (pass == 1) then
      nchecks = nchecks + 1
      dl = abs(lam_shift - lam_ref)
      worst_shift_lam = max(worst_shift_lam, dl / (EPS * ref_absform))
      if (.not. (dl <= tol_lam_factor * EPS * ref_absform)) then
        write(text, '(2a,es10.3)') trim(label), ': shifted quotient |lambda-ref|/(eps*absform)=', &
          dl / (EPS * ref_absform)
        call fail(text)
      end if
      call add_trace(lam_shift, 2.0_wp * tol_lam_factor * EPS * ref_absform)
    else
      call add_trace(lam_shift, 2.0_wp * bound)
    end if
  end subroutine zshifted_quotient_case

  subroutine zcheck_factor_residual(label)
    character(len=*), intent(in) :: label
    integer :: i, info
    real(wp) :: abs_res, rel_res, fwd, bound
    character(len=200) :: text
    do i = 1, n
      zwrk1(i) = 1.0_wp + real(mod(17 * i, 23), wp) / 23.0_wp
      if (mod(i, 2) == 0) zwrk1(i) = -zwrk1(i)
      zwrk1(i) = zwrk1(i) * cmplx(1.0_wp, real(mod(5 * i, 7), wp) / 7.0_wp - 0.4_wp, kind=wp)
    end do
    call zst%factorization_residual(zhmat, zsmat, zwrk1(1:n), abs_res, rel_res, info, fwd)
    call expect_info(info, QR_SUCCESS, trim(label) // ': factorization_residual')
    !With an empty border the probe measures the fresh factorization and
    !must be at the level of the working precision (the bound of the
    !application is 1000*eps*n). With pending functions the probe solution
    !is of moderate size while the cancelling terms of block elimination are
    !not, so the probe shows the full growth factor ||A_base^{-1}*dA||; it is
    !bounded here by 1000 times the bound of the random-right-hand-side
    !solve and reported.
    if (info /= QR_SUCCESS) return
    nchecks = nchecks + 1
    if (zst%pending_count() == 0) then
      worst_facres_fresh = max(worst_facres_fresh, rel_res / EPS)
      bound = 1000.0_wp * real(n, wp)
    else
      worst_facres = max(worst_facres, rel_res / EPS)
      bound = 1000.0_wp * tol_res_factor
    end if
    if (.not. (rel_res <= bound * EPS)) then
      write(text, '(2a,es10.3,a,es10.3,a,i0)') trim(label), ': factorization_residual/eps=', &
        rel_res / EPS, ' forward=', fwd, ' pending=', zst%pending_count()
      call fail(text)
    end if
  end subroutine zcheck_factor_residual

  !=========================================================================
  ! Driving the complex state
  !=========================================================================

  subroutine zdo_flush(shift_rel, retarget)
    !Fresh factorization of the current physical matrices. With retarget the
    !shift is placed shift_rel (relative) below the current reference
    !eigenvalue (the dense reference must be current); otherwise the shift is
    !kept. Only flushes that are scheduled by the random sequence itself
    !retarget the shift: flushes requested by the bookkeeping or by a
    !numerical safeguard may happen at different steps for different process
    !counts, and must not change the represented matrix.
    real(wp), intent(in) :: shift_rel
    logical, intent(in) :: retarget
    integer :: info
    if (retarget) then
      call zref_eigen(zxref, lam_ref)
      shift = lam_ref - shift_rel * max(abs(lam_ref), 1.0e-3_wp)
    end if
    call zst%factorize_fresh(zhmat, zsmat, shift, info, n)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(zst%pending_count() == 0 .and. zst%excluded_count() == 0 .and. &
                zst%base_order() == n .and. zst%order() == n .and. zst%is_valid(), &
                'state after factorize_fresh')
    stat_flush_total = stat_flush_total + 1
    flushed_in_check = .true.
    if (retarget) call zref_build()
  end subroutine zdo_flush

  subroutine zannounce(set, shift_rel, label)
    !begin_active with the flush protocol of the application: a fresh
    !factorization when the border is full (need_flush) or when the committed
    !block fails the pivot test (QR_ERR_SINGULAR).
    integer, intent(in) :: set(:)
    real(wp), intent(in) :: shift_rel
    character(len=*), intent(in) :: label
    logical :: need_flush
    integer :: info
    call zst%begin_active(set, need_flush, info)
    if (info == QR_ERR_SINGULAR .and. allow_flush) then
      call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': begin_active')
      stat_singular_flush = stat_singular_flush + 1
      need_flush = .true.
    else
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active')
      if (need_flush) stat_flush_requested = stat_flush_requested + 1
    end if
    if (need_flush) then
      call zdo_flush(shift_rel, .false.)
      call zst%begin_active(set, need_flush, info)
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active after flush')
      call expect(.not. need_flush, trim(label) // ': need_flush after a flush')
    end if
  end subroutine zannounce

  subroutine zreplace_function(k, label)
    !New column of function k in the physical matrices and in the module.
    integer, intent(in) :: k
    character(len=*), intent(in) :: label
    integer :: info
    call zfunction_column(k, n)
    call zstore_column(k, n)
    call zst%replace_symmetric(k, zhcol(1:n), zscol(1:n), info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': replace_symmetric')
  end subroutine zreplace_function

  subroutine zappend_function(label)
    !Appends function n+1 (parameters already drawn).
    character(len=*), intent(in) :: label
    integer :: info
    call zfunction_column(n + 1, n + 1)
    call zstore_column(n + 1, n + 1)
    call zst%append_symmetric(zhcol(1:n + 1), zscol(1:n + 1), info)
    call expect_info(info, QR_SUCCESS, trim(label) // ': append_symmetric')
    n = n + 1
  end subroutine zappend_function

  subroutine zrun_sequence(name, model_in, n_init, capacity, max_pending, nsteps, &
                          shift_rel, dmax, seed)
    !Random sequence of operations with a check after every step.
    character(len=*), intent(in) :: name
    integer, intent(in) :: model_in, n_init, capacity, max_pending, nsteps, seed
    real(wp), intent(in) :: shift_rel, dmax
    integer :: step, info, nact, act(8), set(8), nset, i, j, k, op, t, base, nmin
    integer :: nfail_before
    logical :: need_flush, found
    real(wp) :: r
    character(len=80) :: label

    nfail_before = nfail
    allow_flush = .true.
    current_shift_rel = shift_rel
    rng_state = int(seed, int64)
    call reset_stats()
    call zproblem_allocate(model_in, capacity, max(2, n_init / 3))
    call zproblem_build(n_init, dmax)
    nmin = max(3, n_init / 2)
    call zst%initialize(capacity, max_pending, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    !Initial shift well below the spectrum, then a shift close to the lowest
    !eigenvalue found with it.
    shift = -2.0_wp * maxval(abs(zdval)) - 10.0_wp
    if (model == MODEL_GAUSS) shift = -5.0_wp
    call zref_build()
    call zdo_flush(shift_rel, .true.)
    call zcheck_linear('fresh')
    call zcheck_eigen('fresh')
    call zcheck_factor_residual('fresh')
    nact = 0

    do step = 1, nsteps
      write(label, '(a,a,i0)') trim(name), ' step ', step
      r = rnd()
      op = 1
      if (r > 0.22_wp) op = 2
      if (r > 0.55_wp) op = 3
      if (r > 0.70_wp) op = 4
      if (r > 0.80_wp) op = 5
      if (r > 0.86_wp) op = 6
      if (r > 0.95_wp) op = 7
      if (op == 2 .and. nact == 0) op = 1

      select case (op)
      case (1)
        !New active set of one or three existing functions: ascending,
        !descending, or scattered, optionally sharing members with the
        !previous active set.
        nset = 1
        if (n >= 6 .and. rnd() < 0.6_wp) nset = 3
        t = irnd(1, 3)
        if (nset == 3 .and. t == 1) then
          k = irnd(1, n - 2)
          set(1:3) = [k, k + 1, k + 2]
        else if (nset == 3 .and. t == 2) then
          k = irnd(1, n - 2)
          set(1:3) = [k + 2, k + 1, k]
        else
          i = 0
          do while (i < nset)
            k = irnd(1, n)
            found = .false.
            do j = 1, i
              if (set(j) == k) found = .true.
            end do
            if (.not. found) then
              i = i + 1
              set(i) = k
            end if
          end do
        end if
        !Shared intersection with the previous set.
        if (nact > 0 .and. rnd() < 0.5_wp) then
          k = act(irnd(1, nact))
          found = .false.
          do j = 1, nset
            if (set(j) == k) found = .true.
          end do
          if (.not. found .and. k <= n) set(irnd(1, nset)) = k
        end if
        call zannounce(set(1:nset), shift_rel, label)
        !The matrix is unchanged by the announcement.
        flushed_in_check = .false.
        call zcheck_linear(trim(label) // ' (announced)')
        if (flushed_in_check) call zannounce(set(1:nset), shift_rel, label)
        do i = 1, nset
          call znew_function(set(i), n)
          call zreplace_function(set(i), label)
        end do
        nact = nset
        act(1:nset) = set(1:nset)
      case (2)
        !Trial: new columns for one or all members of the active set.
        if (rnd() < 0.5_wp) then
          k = act(irnd(1, nact))
          call zperturb_function(k, 0.05_wp)
          call zreplace_function(k, label)
        else
          do i = 1, nact
            call znew_function(act(i), n)
          end do
          do i = 1, nact
            call zreplace_function(act(i), label)
          end do
        end if
      case (3)
        !Append one or two functions, sometimes together with an existing
        !active function.
        t = irnd(1, 2)
        if (n + t > capacity) t = capacity - n
        if (t <= 0) then
          call zst%truncate(n - 1, info)
          if (info == QR_ERR_BATCH_FULL) then
            stat_batch_full = stat_batch_full + 1
            n = n - 1
            call zdo_flush(shift_rel, .false.)
          else
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate at capacity')
            n = n - 1
          end if
          nact = 0
          call zst%begin_active(set(1:0), need_flush, info)
          call expect_info(info, QR_SUCCESS, trim(label) // ': commit')
        else
          nset = 0
          if (rnd() < 0.3_wp) then
            nset = 1
            set(1) = irnd(1, n)
          end if
          do i = 1, t
            set(nset + i) = n + t + 1 - i
          end do
          call zannounce(set(1:nset + t), shift_rel, label)
          if (nset == 1) then
            call znew_function(set(1), n)
            call zreplace_function(set(1), label)
          end if
          nact = nset
          act(1:nset) = set(1:nset)
          do i = 1, t
            call znew_function(n + 1, n)
            call zappend_function(label)
            nact = nact + 1
            act(nact) = n
          end do
        end if
      case (4)
        !Truncate one to three trailing functions, possibly below the base.
        t = irnd(1, 3)
        if (n - t < nmin) t = n - nmin
        if (t > 0) then
          call zst%truncate(n - t, info)
          if (info == QR_ERR_BATCH_FULL) then
            stat_batch_full = stat_batch_full + 1
            call expect(zst%order() == n, trim(label) // ': order changed by a refused truncate')
            n = n - t
            call zdo_flush(shift_rel, .false.)
            nact = 0
          else
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate')
            n = n - t
          end if
          !Forget active indices that no longer exist.
          j = 0
          do i = 1, nact
            if (act(i) <= n) then
              j = j + 1
              act(j) = act(i)
            end if
          end do
          nact = j
        end if
      case (5)
        !Basis-enlargement pattern: several candidates for the same slot.
        if (n < capacity) then
          base = n
          do i = 1, 4
            call zst%truncate(base, info)
            call expect_info(info, QR_SUCCESS, trim(label) // ': truncate to base')
            n = base
            set(1) = base + 1
            call zannounce(set(1:1), shift_rel, label)
            call znew_function(n + 1, n)
            call zappend_function(label)
            if (i < 4) then
              call zref_build()
              call zcheck_linear(trim(label) // ' (candidate)')
            end if
          end do
          nact = 1
          act(1) = n
        end if
      case (6)
        !Flush at random.
        call zdo_flush(shift_rel, .true.)
        call zcheck_factor_residual(trim(label) // ' (fresh)')
        nact = 0
      case (7)
        !Commit the tail by announcing the empty set.
        call zst%begin_active(set(1:0), need_flush, info)
        call expect_info(info, QR_SUCCESS, trim(label) // ': commit')
        call expect(.not. need_flush .and. zst%active_count() == 0, trim(label) // ': commit state')
        nact = 0
      end select

      call zref_build()
      flushed_in_check = .false.
      call zcheck_linear(label)
      if (mod(step, 4) == 0) call zcheck_eigen(label)
      if (mod(step, 25) == 0) call zcheck_factor_residual(label)
      !A flush inside the checks emptied the active set.
      if (flushed_in_check) nact = 0
    end do
    !Numerically triggered flushes must stay exceptional.
    call expect(10 * (stat_resid_flush + stat_singular_flush) <= stat_eigs, &
                trim(name) // ': too many numerically triggered flushes')

    call report_sequence(name, nfail - nfail_before)
    call zst%clear()
  end subroutine zrun_sequence

  subroutine zrun_stress(name, model_in, n_init, max_pending, npending, shift_rel, &
                        dmax, seed)
    !Near-singular stress: the shift lies shift_rel (relative) below the
    !lowest eigenvalue and npending functions are replaced by slightly
    !changed ones in active sets of three without any flush, so that the
    !shifted matrix stays nearly singular throughout.
    character(len=*), intent(in) :: name
    integer, intent(in) :: model_in, n_init, max_pending, npending, seed
    real(wp), intent(in) :: shift_rel, dmax
    integer :: info, i, k, set(3), done, nfail_before, stride
    logical :: need_flush_dummy
    character(len=80) :: label

    nfail_before = nfail
    allow_flush = .false.
    current_shift_rel = shift_rel
    rng_state = int(seed, int64)
    call reset_stats()
    call zproblem_allocate(model_in, n_init + 4, max(2, n_init / 3))
    call zproblem_build(n_init, dmax)
    call zst%initialize(n_init + 4, max_pending, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    shift = -2.0_wp * maxval(abs(zdval)) - 10.0_wp
    if (model == MODEL_GAUSS) shift = -5.0_wp
    call zref_build()
    !Two passes: the first shift is only a rough guide to the eigenvalue.
    call zdo_flush(1.0e-2_wp, .true.)
    call zdo_flush(shift_rel, .true.)
    call zcheck_linear('fresh')
    call zcheck_eigen('fresh')
    done = 0
    stride = max(1, (n - 3) / max(1, npending / 3 + 1))
    k = 1
    do while (done < npending)
      write(label, '(a,a,i0)') trim(name), ' pending ', done + 3
      set = [k + 2, k, k + 1]
      k = k + max(3, stride)
      if (k + 2 > n) k = 1
      call zst%begin_active(set, need_flush_dummy, info)
      call expect_info(info, QR_SUCCESS, trim(label) // ': begin_active')
      call expect(.not. need_flush_dummy, trim(label) // ': unexpected need_flush')
      do i = 1, 3
        call zperturb_function(set(i), 2.0e-3_wp)
        call zreplace_function(set(i), label)
      end do
      done = done + 3
      call zref_build()
      call zcheck_linear(label)
      call zcheck_eigen(label)
    end do
    call expect(zst%pending_count() >= npending, trim(name) // ': pending count')
    call zcheck_factor_residual(name)
    !Trials on the last active set with everything else pending.
    do i = 1, 6
      write(label, '(a,a,i0)') trim(name), ' trial ', i
      call zperturb_function(set(1 + mod(i, 3)), 1.0e-3_wp)
      call zreplace_function(set(1 + mod(i, 3)), label)
      call zref_build()
      call zcheck_linear(label)
      call zcheck_eigen(label)
    end do
    if (myrank == 0) then
      write(*, '(a,es9.2,a,es9.2)') 'stress ' // trim(name) // ': relative distance shift-eigenvalue=', &
        abs(lam_ref - shift) / abs(lam_ref), ' lambda=', lam_ref
    end if
    call report_sequence(name, nfail - nfail_before)
    call zst%clear()
  end subroutine zrun_stress

  subroutine zrun_error_paths()
    !Error paths and state preservation.
    integer :: info, num_iter, nfail_before, i, j
    logical :: need_flush
    real(wp) :: lambda, rel_acc, eig_res, abs_res, rel_res
    complex(wp) :: h4(4, 4), s4(4, 4), h2(2, 2), s2(2, 2), v4(4), x4(4)
    complex(wp) :: c2(2), z2(2), ci
    type(qr_complex_state) :: e

    nfail_before = nfail
    ci = cmplx(0.0_wp, 1.0_wp, kind=wp)
    rng_state = 4711_int64
    call zproblem_allocate(MODEL_RITZ, 10, 3)
    call zproblem_build(8, 100.0_wp)

    !Uninitialized state.
    call e%factorize_fresh(zhmat, zsmat, -50.0_wp, info, 8)
    call expect_info(info, QR_ERR_INVALID_STATE, 'factorize_fresh on uninitialized state')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'begin_active on uninitialized state')
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'solve on uninitialized state')
    call e%shifted_quotient(zhmat, zsmat, zwrk1(1:8), lambda, info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'shifted_quotient on uninitialized state')
    call expect(.not. e%is_valid() .and. e%order() == 0, 'queries on uninitialized state')

    !Initialization arguments.
    call e%initialize(0, 4, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'initialize capacity=0')
    call e%initialize(10, 0, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'initialize max_pending=0')
    call e%initialize(10, 99, info)
    call expect_info(info, QR_SUCCESS, 'initialize with max_pending above capacity')
    call expect(e%max_pending() == 10 .and. e%get_capacity() == 10, 'max_pending clamp')
    call e%initialize(10, 3, info)
    call expect_info(info, QR_SUCCESS, 'initialize')
    call expect(e%max_pending() == 3 .and. .not. e%is_valid(), 'state after initialize')
    call e%replace_symmetric(1, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_STATE, 'replace before factorization')

    !Factorization arguments.
    call e%factorize_fresh(zhmat, zsmat, -50.0_wp, info, 11)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'factorize_fresh order above capacity')
    call e%factorize_fresh(zhmat, zsmat, -50.0_wp, info, 0)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'factorize_fresh order zero')
    call e%factorize_fresh(zhmat(1:8, 1:7), zsmat(1:8, 1:7), -50.0_wp, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorize_fresh non-square')
    call e%factorize_fresh(zhmat(1:7, 1:7), zsmat, -50.0_wp, info, 8)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorize_fresh too small H')
    shift = -50.0_wp
    call e%factorize_fresh(zhmat, zsmat, shift, info, 8)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(e%is_valid() .and. e%order() == 8 .and. e%base_order() == 8 .and. &
                e%get_shift() == shift, 'queries after factorize_fresh')

    !begin_active arguments.
    call e%begin_active([0], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active index 0')
    call e%begin_active([3, 3], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active repeated index')
    call e%begin_active([10], need_flush, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'begin_active gap in appended range')
    call e%begin_active([9, 10, 11], need_flush, info)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'begin_active append beyond capacity')
    call e%begin_active([1, 2, 3, 4], need_flush, info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'begin_active set larger than max_pending')
    call expect(e%pending_count() == 0 .and. e%active_count() == 0, 'state after refused begin_active')

    !replace_symmetric arguments.
    call zfunction_column(2, 8)
    call e%replace_symmetric(2, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace of an unannounced function')
    call e%begin_active([3, 1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [3,1,2]')
    call expect(.not. need_flush .and. e%pending_count() == 3 .and. e%active_count() == 3 &
                .and. e%excluded_count() == 3, 'state after begin_active')
    call e%replace_symmetric(0, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace index 0')
    call e%replace_symmetric(9, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace index above order')
    call e%replace_symmetric(2, zhcol(1:7), zscol(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'replace with short column')
    call e%replace_symmetric(5, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'replace of a clean function')
    call e%replace_symmetric(2, zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_SUCCESS, 'replace with an unchanged column')
    call expect(e%get_update_count() == 1_int64 .and. e%get_updates_since_fresh() == 1_int64, &
                'update counters')

    !Batch full: need_flush and unchanged state.
    call e%begin_active([4], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active on a full border')
    call expect(need_flush, 'need_flush on a full border')
    call expect(e%pending_count() == 3 .and. e%active_count() == 3, 'state after need_flush')
    call e%begin_active([2, 3, 1], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active same set')
    call expect(.not. need_flush .and. e%active_count() == 3, 'same set is a no-op')
    call e%begin_active([9], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active append on a full border')
    call expect(need_flush, 'need_flush for an append on a full border')
    call zfunction_column(9, 9)
    call znew_function(9, 8)
    call zfunction_column(9, 9)
    call e%append_symmetric(zhcol(1:9), zscol(1:9), info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'unannounced append on a full border')
    call expect(e%order() == 8, 'order after refused append')
    call e%append_symmetric(zhcol(1:8), zscol(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'append with short column')

    !truncate / delete arguments.
    call e%truncate(0, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'truncate to zero')
    call e%truncate(9, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'truncate above order')
    call e%delete_symmetric(3, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'delete interior function')
    call e%truncate(7, info)
    call expect_info(info, QR_ERR_BATCH_FULL, 'truncate with full D0')
    call expect(e%order() == 8, 'order after refused truncate')
    call e%truncate(8, info)
    call expect_info(info, QR_SUCCESS, 'truncate no-op')

    !solve arguments.
    zwrk1(1:8) = 0.0_wp
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_ZERO_INITIAL_VECTOR, 'solve with zero start vector')
    call expect(lambda == 0.0_wp .and. num_iter == 0 .and. all(zxvec(1:8) == 0.0_wp), &
                'outputs after zero start vector')
    zwrk1(1:8) = 1.0_wp
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e-8_wp, 0, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'solve with max_iter=0')
    call e%solve(zhmat, zsmat, zwrk1(1:7), zxvec(1:8), lambda, 1.0e-8_wp, 10, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve with short start vector')
    call e%solve(zhmat(1:7, 1:7), zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e-8_wp, 10, 0, &
                 rel_acc, num_iter, eig_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve with small H')
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e-30_wp, 1, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_ERR_NO_CONVERGENCE, 'solve with max_iter=1')
    call expect(num_iter == 1 .and. maxval(abs(zxvec(1:8))) > 0.0_wp .and. &
                eig_res < 1.0_wp, 'approximation after no convergence')
    !Normalization modes with the far shift (slow but sure convergence).
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e3_wp * EPS, 5000, 1, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_SUCCESS, 'solve norm_mode=1')
    call expect(abs(sum(abs(zxvec(1:8))**2) - 1.0_wp) < 100.0_wp * EPS, 'norm_mode=1 normalization')
    call e%solve(zhmat, zsmat, zwrk1(1:8), zxvec(1:8), lambda, 1.0e3_wp * EPS, 5000, 2, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_SUCCESS, 'solve norm_mode=2')
    call expect(abs(maxval(max(abs(real(zxvec(1:8), wp)), abs(aimag(zxvec(1:8))))) - 1.0_wp) < &
                100.0_wp * EPS, 'norm_mode=2 normalization')
    call e%factorization_residual(zhmat, zsmat, zwrk1(1:7), abs_res, rel_res, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'factorization_residual short probe')
    zwrk2(1:8) = 0.0_wp
    call e%factorization_residual(zhmat, zsmat, zwrk2(1:8), abs_res, rel_res, info)
    call expect_info(info, QR_ERR_ZERO_INITIAL_VECTOR, 'factorization_residual zero probe')
    call e%solve_shifted(zwrk1(1:7), zxvec(1:8), info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'solve_shifted short right-hand side')
    lambda = 1.0_wp
    call e%shifted_quotient(zhmat, zsmat, zwrk1(1:7), lambda, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'shifted_quotient short vector')
    call expect(lambda == 0.0_wp, 'shifted_quotient returns zero on a dimension error')
    call e%shifted_quotient(zhmat(1:7, 1:8), zsmat, zwrk1(1:8), lambda, info)
    call expect_info(info, QR_ERR_DIMENSION_MISMATCH, 'shifted_quotient short H')
    zwrk2(1:8) = 0.0_wp
    lambda = 1.0_wp
    call e%shifted_quotient(zhmat, zsmat, zwrk2(1:8), lambda, info)
    call expect_info(info, QR_ERR_NONPOSITIVE_OVERLAP, 'shifted_quotient zero vector')
    call expect(lambda == 0.0_wp, 'shifted_quotient returns zero for a zero vector')

    !Append at capacity.
    call e%factorize_fresh(zhmat, zsmat, shift, info, 8)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh')
    call expect(e%get_updates_since_fresh() == 0_int64 .and. e%get_update_count() == 1_int64, &
                'counters after factorize_fresh')
    n = 8
    do i = 9, 10
      call znew_function(i, i - 1)
      call zfunction_column(i, i)
      call zstore_column(i, i)
      call e%append_symmetric(zhcol(1:i), zscol(1:i), info)
      call expect_info(info, QR_SUCCESS, 'unannounced append')
    end do
    call expect(e%order() == 10 .and. e%pending_count() == 2 .and. e%active_count() == 2, &
                'state after unannounced appends')
    call e%append_symmetric(zhcol(1:10), zscol(1:10), info)
    call expect_info(info, QR_ERR_CAPACITY_EXCEEDED, 'append at capacity')
    call e%delete_symmetric(10, info)
    call expect_info(info, QR_SUCCESS, 'delete last function')
    call expect(e%order() == 9 .and. e%pending_count() == 1, 'state after delete')

    !Exactly singular base matrix: H = shift*S.
    do j = 1, 4
      do i = 1, 4
        s4(i, j) = 0.0_wp
        h4(i, j) = 0.0_wp
      end do
      s4(j, j) = 1.0_wp
      h4(j, j) = 2.0_wp
    end do
    v4 = 1.0_wp
    call e%initialize(4, 2, info)
    call e%factorize_fresh(h4, s4, 2.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh of a singular matrix')
    call e%solve(h4, s4, v4, x4, lambda, 1.0e-8_wp, 10, 0, rel_acc, num_iter, eig_res, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve with singular R')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular R')

    !Singular small matrix ZDD: A is nonsingular but A without function 1 is
    !singular, A = [0 -i 0 0; i 0 0 0; 0 0 1 0; 0 0 0 1]. The right-hand side
    !is A times the vector of ones.
    h4 = 0.0_wp
    h4(2, 1) = ci
    h4(3, 3) = 1.0_wp
    h4(4, 4) = 1.0_wp
    v4 = [-ci, ci, cmplx(1.0_wp, 0.0_wp, kind=wp), cmplx(1.0_wp, 0.0_wp, kind=wp)]
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'factorize_fresh (ZDD case)')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted before the singular border')
    call expect(maxval(abs(x4 - 1.0_wp)) < 100.0_wp * EPS, 'solution of the permutation system')
    call e%begin_active([1], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular committed block')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular committed block')
    !Columns can still be stored; only the solve needs the singular block.
    call e%replace_symmetric(1, h4(:, 1), s4(:, 1), info)
    call expect_info(info, QR_SUCCESS, 'replace with singular committed block')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted after the replace')
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call expect_info(info, QR_SUCCESS, 'recovery by factorize_fresh')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after recovery')

    !Singular committed Schur complement Scc: commit a replacement that makes
    !the matrix without the next active function singular. A = diag(1,1,1,1);
    !replace column 1 by e2 (A stays nonsingular with column 2 = e1 later).
    do j = 1, 4
      do i = 1, 4
        h4(i, j) = 0.0_wp
      end do
      h4(j, j) = 1.0_wp
    end do
    call e%factorize_fresh(h4, s4, 0.0_wp, info)
    call e%begin_active([1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [1,2]')
    !New matrix: rows/columns 1 and 2 form [0 -i; i 0]. The conceptual column
    !of function 1 is (0,i,0,0) and that of function 2 is (-i,0,0,0).
    h4(1, 1) = 0.0_wp
    h4(2, 1) = ci
    h4(2, 2) = 0.0_wp
    v4 = 0.0_wp
    v4(2) = ci
    x4 = 0.0_wp
    call e%replace_symmetric(1, v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'replace 1')
    v4 = 0.0_wp
    v4(1) = -ci
    call e%replace_symmetric(2, v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'replace 2')
    v4 = [-ci, ci, cmplx(1.0_wp, 0.0_wp, kind=wp), cmplx(1.0_wp, 0.0_wp, kind=wp)]
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted with the 2x2 swap')
    call expect(maxval(abs(x4 - 1.0_wp)) < 100.0_wp * EPS, 'solution with the 2x2 swap')
    !Committing function 1 alone leaves Scc = (0): singular.
    call e%begin_active([2], need_flush, info)
    call expect_info(info, QR_ERR_SINGULAR, 'begin_active with singular Scc')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular Scc')
    !A different partition is usable again.
    call e%begin_active([1, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active back to [1,2]')
    call e%solve_shifted(v4, x4, info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after repartition')

    !Singular active Schur complement Sigma: the represented matrix itself
    !becomes exactly singular, A = [1 i; -i 1] (conceptual column 2 = (i,1)).
    !The base is the identity, so that every stored quantity is exact.
    h2 = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
    s2 = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2, 2])
    call e%initialize(2, 2, info)
    call e%factorize_fresh(h2, s2, 0.0_wp, info)
    call e%begin_active([2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'begin_active [2] (Sigma case)')
    z2 = 0.0_wp
    c2 = [ci, cmplx(1.0_wp, 0.0_wp, kind=wp)]
    call e%replace_symmetric(2, c2, z2, info)
    call expect_info(info, QR_SUCCESS, 'replace giving a singular matrix')
    c2 = [cmplx(1.0_wp, 0.0_wp, kind=wp), cmplx(2.0_wp, 0.0_wp, kind=wp)]
    call e%solve_shifted(c2, x4(1:2), info)
    call expect_info(info, QR_ERR_SINGULAR, 'solve_shifted with singular Sigma')
    !A = [1 i; -i 3]; the right-hand side is A times the vector of ones.
    c2 = [ci, cmplx(3.0_wp, 0.0_wp, kind=wp)]
    call e%replace_symmetric(2, c2, z2, info)
    call expect_info(info, QR_SUCCESS, 'replace giving a regular matrix')
    c2 = [1.0_wp + ci, 3.0_wp - ci]
    call e%solve_shifted(c2, x4(1:2), info)
    call expect_info(info, QR_SUCCESS, 'solve_shifted after the regular replace')
    call expect(maxval(abs(x4(1:2) - 1.0_wp)) < 100.0_wp * EPS, 'solution after the regular replace')

    call e%clear()
    call expect(.not. e%is_valid() .and. e%get_capacity() == 0 .and. e%max_pending() == 0, &
                'state after clear')
    call expect(trim(qr_status_message(QR_ERR_BATCH_FULL)) /= trim(qr_status_message(-1)), &
                'status message for QR_ERR_BATCH_FULL')
    if (myrank == 0) write(*, '(a,i0)') 'complex error paths: failures=', nfail - nfail_before
  end subroutine zrun_error_paths

  !=========================================================================
  ! Explicit tests of the Hermitian conventions
  !=========================================================================

  subroutine zconv_full(m, nord, sigma, h, s, hf, sf, af)
    !Full Hermitian H, S and shifted matrix of order nord from the lower
    !triangles of h and s (helper of zrun_conventions).
    integer, intent(in) :: m, nord
    real(wp), intent(in) :: sigma
    complex(wp), intent(in) :: h(m, m), s(m, m)
    complex(wp), intent(inout) :: hf(m, m), sf(m, m), af(m, m)
    integer :: ii, jj
    do jj = 1, nord
      do ii = jj, nord
        hf(ii, jj) = h(ii, jj)
        sf(ii, jj) = s(ii, jj)
        hf(jj, ii) = conjg(h(ii, jj))
        sf(jj, ii) = conjg(s(ii, jj))
      end do
    end do
    af(1:nord, 1:nord) = hf(1:nord, 1:nord) - sigma * sf(1:nord, 1:nord)
  end subroutine zconv_full

  subroutine zconv_store(m, k, norder, hc, sc, h, s)
    !Stores the conceptual columns hc, sc of function k in the lower
    !triangles of h and s (helper of zrun_conventions).
    integer, intent(in) :: m, k, norder
    complex(wp), intent(in) :: hc(m), sc(m)
    complex(wp), intent(inout) :: h(m, m), s(m, m)
    integer :: jj
    do jj = 1, norder
      if (jj >= k) then
        h(jj, k) = hc(jj)
        s(jj, k) = sc(jj)
      else
        h(k, jj) = conjg(hc(jj))
        s(k, jj) = conjg(sc(jj))
      end if
    end do
  end subroutine zconv_store

  subroutine zconv_solve(e, m, nord, af, xt, tol, text)
    !solve_shifted with the right-hand side af*xt must return xt (helper of
    !zrun_conventions).
    type(qr_complex_state), intent(inout) :: e
    integer, intent(in) :: m, nord
    complex(wp), intent(in) :: af(m, m), xt(m)
    real(wp), intent(in) :: tol
    character(len=*), intent(in) :: text
    integer :: status
    complex(wp) :: b(m), x(m)
    b(1:nord) = matmul(af(1:nord, 1:nord), xt(1:nord))
    call e%solve_shifted(b(1:nord), x(1:nord), status)
    call expect_info(status, QR_SUCCESS, text // ': solve_shifted')
    call expect(maxval(abs(x(1:nord) - xt(1:nord))) <= tol, text)
  end subroutine zconv_solve

  subroutine zrun_conventions()
    !Small hand-made problem that pins down the conventions of the complex
    !state independently of the random suites: which triangle is read, what a
    !supplied column means, how the diagonal is treated, and what
    !symmetric_product computes. af is the full shifted Hermitian matrix kept
    !by the test; every solve is compared with a right-hand side af*xt.
    integer, parameter :: m = 5
    integer :: info, i, j, nfail_before, nord
    logical :: need_flush
    real(wp) :: sigma, tol, lambda, rel_acc, eig_res
    integer :: num_iter
    complex(wp) :: h(m, m), s(m, m), hbad(m, m), af(m, m), hf(m, m), sf(m, m)
    complex(wp) :: xt(m), b(m), x(m), y(m), hc(m), sc(m), ci, v(m)
    type(qr_complex_state) :: e

    nfail_before = nfail
    ci = cmplx(0.0_wp, 1.0_wp, kind=wp)
    sigma = 0.5_wp
    tol = 1.0e4_wp * EPS
    !Garbage everywhere, then the lower triangles of the leading 4 by 4
    !blocks: stored (i,j), i>=j, is the element (i,j).
    do j = 1, m
      do i = 1, m
        h(i, j) = cmplx(7.0e20_wp, -3.0e20_wp, kind=wp)
        s(i, j) = cmplx(-5.0e20_wp, 2.0e20_wp, kind=wp)
      end do
    end do
    h(1, 1) = 2.0_wp
    h(2, 1) = cmplx(0.5_wp, 0.25_wp, kind=wp)
    h(3, 1) = cmplx(0.1_wp, -0.3_wp, kind=wp)
    h(4, 1) = cmplx(0.0_wp, 0.3_wp, kind=wp)
    h(2, 2) = 3.0_wp
    h(3, 2) = cmplx(-0.2_wp, 0.4_wp, kind=wp)
    h(4, 2) = 0.6_wp
    h(3, 3) = 4.0_wp
    h(4, 3) = cmplx(-0.1_wp, -0.2_wp, kind=wp)
    h(4, 4) = 5.0_wp
    s(1, 1) = 1.0_wp
    s(2, 1) = cmplx(0.1_wp, 0.05_wp, kind=wp)
    s(3, 1) = cmplx(0.0_wp, 0.02_wp, kind=wp)
    s(4, 1) = 0.04_wp
    s(2, 2) = 1.0_wp
    s(3, 2) = cmplx(-0.03_wp, 0.01_wp, kind=wp)
    s(4, 2) = cmplx(0.01_wp, -0.02_wp, kind=wp)
    s(3, 3) = 1.0_wp
    s(4, 3) = cmplx(0.0_wp, 0.03_wp, kind=wp)
    s(4, 4) = 1.0_wp
    nord = 4
    call zconv_full(m, nord, sigma, h, s, hf, sf, af)
    xt = [cmplx(1.0_wp, 1.0_wp, kind=wp), cmplx(2.0_wp, -1.0_wp, kind=wp), &
          cmplx(-1.0_wp, 0.5_wp, kind=wp), cmplx(0.3_wp, 0.0_wp, kind=wp), &
          cmplx(-0.7_wp, 0.2_wp, kind=wp)]

    call e%initialize(m, 3, info)
    call expect_info(info, QR_SUCCESS, 'conventions: initialize')
    call e%factorize_fresh(h, s, sigma, info, nord)
    call expect_info(info, QR_SUCCESS, 'conventions: factorize_fresh')
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: lower triangle is <phi_i|O|phi_j>, upper is its conjugate')

    !symmetric_product is the Hermitian product; the imaginary part of the
    !diagonal and the upper triangle are ignored.
    hbad = h
    do i = 1, nord
      hbad(i, i) = h(i, i) + 0.7_wp * ci
    end do
    call e%symmetric_product(hbad(1:nord, 1:nord), xt(1:nord), y(1:nord), info)
    call expect_info(info, QR_SUCCESS, 'conventions: symmetric_product')
    b(1:nord) = matmul(hf(1:nord, 1:nord), xt(1:nord))
    call expect(maxval(abs(y(1:nord) - b(1:nord))) <= tol, &
                'conventions: symmetric_product is the Hermitian product of the lower triangle')

    !shifted_quotient reads the same Hermitian matrices: garbage in the upper
    !triangles and imaginary parts of the diagonal of H do not enter, and the
    !result is sigma + x^H*(H-sigma*S)*x/(x^H*S*x) for the full matrices.
    call e%shifted_quotient(hbad, s, xt(1:nord), lambda, info)
    call expect_info(info, QR_SUCCESS, 'conventions: shifted_quotient')
    b(1:nord) = matmul(hf(1:nord, 1:nord), xt(1:nord))
    y(1:nord) = matmul(sf(1:nord, 1:nord), xt(1:nord))
    call expect(abs(lambda - real(sum(conjg(xt(1:nord)) * b(1:nord)), wp) / &
                    real(sum(conjg(xt(1:nord)) * y(1:nord)), wp)) <= tol, &
                'conventions: shifted_quotient is the Hermitian quotient of the lower triangles')

    !A diagonal that is not real is refused and the state is preserved.
    call e%factorize_fresh(hbad, s, sigma, info, nord)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'conventions: factorize_fresh with complex diagonal of H')
    hbad = s
    hbad(3, 3) = s(3, 3) - 0.01_wp * ci
    call e%factorize_fresh(h, hbad, sigma, info, nord)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'conventions: factorize_fresh with complex diagonal of S')
    call expect(e%is_valid() .and. e%order() == nord, 'conventions: state after refused factorize_fresh')
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: factors after refused factorize_fresh')
    !A roundoff-sized imaginary part of the diagonal is discarded.
    hbad = h
    do i = 1, nord
      hbad(i, i) = h(i, i) + 20.0_wp * EPS * ci
    end do
    call e%factorize_fresh(hbad, s, sigma, info, nord)
    call expect_info(info, QR_SUCCESS, 'conventions: factorize_fresh with roundoff on the diagonal')
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: roundoff on the diagonal is discarded')

    !replace_symmetric: element i of the column of function k is
    !<phi_i|O|phi_k>, for i below and above k alike.
    call e%begin_active([2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'conventions: begin_active [2]')
    hc(1:4) = [cmplx(-0.4_wp, 0.7_wp, kind=wp), cmplx(2.5_wp, 0.0_wp, kind=wp), &
               cmplx(0.3_wp, 0.2_wp, kind=wp), cmplx(-0.6_wp, -0.1_wp, kind=wp)]
    sc(1:4) = [cmplx(0.05_wp, -0.08_wp, kind=wp), cmplx(1.1_wp, 0.0_wp, kind=wp), &
               cmplx(0.02_wp, 0.03_wp, kind=wp), cmplx(-0.01_wp, 0.04_wp, kind=wp)]
    v(1:4) = hc(1:4)
    v(2) = hc(2) + 0.3_wp * ci
    call e%replace_symmetric(2, v(1:4), sc(1:4), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'conventions: replace with complex diagonal of H')
    v(1:4) = sc(1:4)
    v(2) = sc(2) - 0.3_wp * ci
    call e%replace_symmetric(2, hc(1:4), v(1:4), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'conventions: replace with complex diagonal of S')
    call expect(e%get_updates_since_fresh() == 0_int64, 'conventions: counters after refused replace')
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: matrix after refused replace')
    v(1:4) = hc(1:4)
    v(2) = hc(2) + 20.0_wp * EPS * ci
    call e%replace_symmetric(2, v(1:4), sc(1:4), info)
    call expect_info(info, QR_SUCCESS, 'conventions: replace_symmetric')
    call zconv_store(m, 2, 4, hc, sc, h, s)
    call zconv_full(m, nord, sigma, h, s, hf, sf, af)
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: replaced column is the conceptual column')

    !append_symmetric with the same convention.
    call e%begin_active([5, 2], need_flush, info)
    call expect_info(info, QR_SUCCESS, 'conventions: begin_active [5,2]')
    hc(1:5) = [cmplx(0.2_wp, -0.5_wp, kind=wp), cmplx(-0.3_wp, 0.1_wp, kind=wp), &
               cmplx(0.4_wp, 0.4_wp, kind=wp), cmplx(0.0_wp, -0.2_wp, kind=wp), &
               cmplx(6.0_wp, 0.0_wp, kind=wp)]
    sc(1:5) = [cmplx(0.03_wp, 0.02_wp, kind=wp), cmplx(0.0_wp, -0.05_wp, kind=wp), &
               cmplx(-0.02_wp, 0.0_wp, kind=wp), cmplx(0.01_wp, 0.01_wp, kind=wp), &
               cmplx(1.0_wp, 0.0_wp, kind=wp)]
    v(1:5) = hc(1:5)
    v(5) = hc(5) - 0.2_wp * ci
    call e%append_symmetric(v(1:5), sc(1:5), info)
    call expect_info(info, QR_ERR_INVALID_ARGUMENT, 'conventions: append with complex diagonal')
    call expect(e%order() == 4, 'conventions: order after refused append')
    call e%append_symmetric(hc(1:5), sc(1:5), info)
    call expect_info(info, QR_SUCCESS, 'conventions: append_symmetric')
    nord = 5
    call zconv_store(m, 5, 5, hc, sc, h, s)
    call zconv_full(m, nord, sigma, h, s, hf, sf, af)
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: appended column is the conceptual column')
    !Second trial of function 2 with the appended function present.
    hc(1:5) = [cmplx(0.1_wp, 0.9_wp, kind=wp), cmplx(2.0_wp, 0.0_wp, kind=wp), &
               cmplx(-0.3_wp, 0.5_wp, kind=wp), cmplx(0.2_wp, 0.2_wp, kind=wp), &
               cmplx(0.6_wp, -0.3_wp, kind=wp)]
    sc(1:5) = [cmplx(-0.02_wp, 0.06_wp, kind=wp), cmplx(0.9_wp, 0.0_wp, kind=wp), &
               cmplx(0.04_wp, -0.01_wp, kind=wp), cmplx(0.0_wp, 0.02_wp, kind=wp), &
               cmplx(0.03_wp, 0.05_wp, kind=wp)]
    call e%replace_symmetric(2, hc(1:5), sc(1:5), info)
    call expect_info(info, QR_SUCCESS, 'conventions: replace with an appended function')
    call zconv_store(m, 2, 5, hc, sc, h, s)
    call zconv_full(m, nord, sigma, h, s, hf, sf, af)
    call zconv_solve(e, m, nord, af, xt, tol, 'conventions: column shared by a base and an appended function')

    !Eigenpair: real eigenvalue, x^H S x = 1, Hermitian residual.
    do i = 1, nord
      v(i) = zstart(i)
    end do
    call e%factorize_fresh(h, s, -3.0_wp, info, nord)
    call expect_info(info, QR_SUCCESS, 'conventions: factorize_fresh (eigenpair)')
    call e%solve(h, s, v(1:nord), x(1:nord), lambda, 1.0e3_wp * EPS, 5000, 0, rel_acc, &
                 num_iter, eig_res, info)
    call expect_info(info, QR_SUCCESS, 'conventions: solve')
    b(1:nord) = matmul(sf(1:nord, 1:nord), x(1:nord))
    y(1:nord) = matmul(hf(1:nord, 1:nord), x(1:nord))
    call expect(abs(real(sum(conjg(x(1:nord)) * b(1:nord)), wp) - 1.0_wp) <= tol, &
                'conventions: x^H S x = 1')
    call expect(maxval(abs(y(1:nord) - lambda * b(1:nord))) <= 1.0e4_wp * tol, &
                'conventions: H x = lambda S x')
    call expect(abs(lambda - real(sum(conjg(x(1:nord)) * y(1:nord)), wp)) <= tol * 10.0_wp, &
                'conventions: lambda = x^H H x')
    call e%clear()
    if (myrank == 0) write(*, '(a,i0)') 'complex conventions: failures=', nfail - nfail_before

  end subroutine zrun_conventions

  subroutine zrun_unit_tests()
    !Complex Hermitian counterpart of run_unit_tests: the same kinds of
    !sequences, sizes and tolerances, with other seeds. The names of the
    !sequences carry the prefix Z.
    integer :: steps
    logical :: high_precision
    tol_res_factor = 1.0e6_wp
    tol_fwd_factor = 1.0e6_wp
    tol_lam_factor = 1.0e2_wp
    tol_vec_factor = 1.0e2_wp
    tol_eigres_factor = 100.0_wp
    tol_plain_factor = 1.0_wp
    high_precision = (precision(1.0_wp) > 15)

    call zrun_conventions()
    call zrun_error_paths()
#ifdef QRLINALG_TESTING
    call zrun_factor_tests()
#endif

    steps = 320
    if (precision(1.0_wp) > 20) steps = 100
    if (quick) steps = 40
    !Well-conditioned random problems: the algebra and the bookkeeping.
    scale_decades = 0.5_wp
    dependence = 0.5_wp
    call zrun_sequence('ZA ritz n0=120 cap=160 mmax=16', MODEL_RITZ, 120, 160, 16, &
                       steps, 1.0e-3_wp, 1.0e3_wp, 30241)
    !Small border: frequent flushes requested by begin_active.
    call zrun_sequence('ZB ritz n0=61 cap=75 mmax=6', MODEL_RITZ, 61, 75, 6, &
                       steps, 1.0e-4_wp, 1.0e5_wp, 30242)
    !Fewer functions than processes when run with 14 ranks.
    call zrun_sequence('ZC ritz n0=9 cap=13 mmax=4', MODEL_RITZ, 9, 13, 4, &
                       steps, 1.0e-4_wp, 1.0e2_wp, 30243)
    !Hydrogen atom in complex s-type Gaussians (see run_unit_tests).
    min_ratio = 1.3_wp
    if (high_precision) min_ratio = 1.2_wp
    tol_res_factor = 1.0e9_wp
    tol_fwd_factor = 1.0e9_wp
    call zrun_sequence('ZD gauss n0=24 cap=32 mmax=8', MODEL_GAUSS, 24, 32, 8, &
                       steps, 1.0e-5_wp, 0.0_wp, 30244)
    tol_res_factor = 1.0e6_wp
    tol_fwd_factor = 1.0e6_wp
    !Larger problem with badly scaled, nearly dependent functions.
    scale_decades = 1.0_wp
    dependence = 0.1_wp
    if (precision(1.0_wp) <= 20 .and. .not. quick) then
      call zrun_sequence('ZE ritz n0=257 cap=300 mmax=40', MODEL_RITZ, 257, 300, 40, &
                         160, 1.0e-5_wp, 1.0e4_wp, 30245)
    else
      call zrun_sequence('ZE ritz n0=101 cap=120 mmax=24', MODEL_RITZ, 101, 120, 24, &
                         40, 1.0e-5_wp, 1.0e4_wp, 30245)
    end if

    !Near-singular stress cases (see run_unit_tests).
    scale_decades = 1.5_wp
    dependence = 0.05_wp
    if (high_precision) dependence = 0.01_wp
    if (precision(1.0_wp) <= 20 .and. .not. quick) then
      call zrun_stress('ZS1 ritz n=300 shift 1e-6 pending 36', MODEL_RITZ, 300, 48, 36, &
                       1.0e-6_wp, 1.0e4_wp, 30246)
    else
      call zrun_stress('ZS1 ritz n=120 shift 1e-6 pending 33', MODEL_RITZ, 120, 40, 33, &
                       1.0e-6_wp, 1.0e4_wp, 30246)
    end if
    min_ratio = 1.3_wp
    if (high_precision) min_ratio = 1.15_wp
    call zrun_stress('ZS2 gauss n=40 shift 1e-6 pending 30', MODEL_GAUSS, 40, 36, 30, &
                     1.0e-6_wp, 0.0_wp, 30247)
    call zproblem_free()
  end subroutine zrun_unit_tests

  !=========================================================================
  ! Complex timing mode
  !=========================================================================

  subroutine ztiming_matrices(order)
    !Synthetic dense Hermitian matrices for the timing mode: S = I + small
    !random Hermitian part (positive definite); H = diagonal of positive
    !energies (a few of them large) + small random Hermitian part - 8*q*q^H
    !with a dense complex vector q, which creates one well separated bound
    !state whose eigenvector involves all functions.
    integer, intent(in) :: order
    integer :: i, j
    real(wp) :: rs, theta
    rs = 1.0_wp / sqrt(real(order, wp))
    do i = 1, order
      theta = 8.0_wp * atan(1.0_wp) * rnd()
      zwrk1(i) = rs * (0.5_wp + rnd()) * cmplx(cos(theta), sin(theta), kind=wp)
      zrw4(i) = 0.5_wp + 4.5_wp * rnd()
      if (mod(i, 50) == 0) zrw4(i) = 10.0_wp**(1.0_wp + 5.0_wp * rnd())
    end do
    do j = 1, order
      do i = j, order
        zsmat(i, j) = 0.2_wp * rs * zrnd()
        zhmat(i, j) = 0.3_wp * rs * zrnd() - 8.0_wp * zwrk1(i) * conjg(zwrk1(j))
      end do
      zsmat(j, j) = cmplx(real(zsmat(j, j), wp) + 1.0_wp, 0.0_wp, kind=wp)
      zhmat(j, j) = cmplx(real(zhmat(j, j), wp) + zrw4(j), 0.0_wp, kind=wp)
    end do
    n = order
  end subroutine ztiming_matrices

  subroutine ztiming_perturb(k, amount)
    !New conceptual columns zhcol, zscol of function k: the current ones with
    !a relative random change (real factors, so that the diagonal stays
    !real); the physical lower triangles are updated.
    integer, intent(in) :: k
    real(wp), intent(in) :: amount
    integer :: j
    do j = 1, n
      if (j >= k) then
        zhcol(j) = zhmat(j, k)
        zscol(j) = zsmat(j, k)
      else
        zhcol(j) = conjg(zhmat(k, j))
        zscol(j) = conjg(zsmat(k, j))
      end if
      zhcol(j) = zhcol(j) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
      zscol(j) = zscol(j) * (1.0_wp + amount * (2.0_wp * rnd() - 1.0_wp))
    end do
    call zstore_column(k, n)
  end subroutine ztiming_perturb

  subroutine ztime_trials(label, k, ntrials, amount, tol)
    !Trials of the active function k: one replace_symmetric followed by solve
    !started from the previous eigenvector.
    character(len=*), intent(in) :: label
    integer, intent(in) :: k, ntrials
    real(wp), intent(in) :: amount, tol
    integer :: it, info, num_iter, iters, i
    real(wp) :: lambda, rel_acc, eig_res, worst
    real(kind(1.0d0)) :: t0, t1, t2, t_replace, t_solve, t_min
    t_replace = 0.0d0
    t_solve = 0.0d0
    t_min = huge(1.0d0)
    iters = 0
    worst = 0.0_wp
    lambda = 0.0_wp
    do it = 1, ntrials
      call ztiming_perturb(k, amount)
      t0 = wall()
      call zst%replace_symmetric(k, zhcol(1:n), zscol(1:n), info)
      t1 = wall()
      if (info /= QR_SUCCESS) call fail(label // ': replace_symmetric')
      do i = 1, n
        zwrk1(i) = zxvec(i)
      end do
      call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, tol, 50, 0, rel_acc, &
                    num_iter, eig_res, info)
      t2 = wall()
      if (info /= QR_SUCCESS) call fail(label // ': solve')
      t_replace = t_replace + (t1 - t0)
      t_solve = t_solve + (t2 - t1)
      t_min = min(t_min, t2 - t0)
      iters = iters + num_iter
      worst = max(worst, eig_res)
    end do
    lam_ref = lambda
    if (myrank == 0) then
      write(*, '(a,a,i0,a,i0)') trim(label), ': pending=', zst%pending_count(), &
        ' excluded=', zst%excluded_count()
      write(*, '(a,f9.3,a,f9.3,a,f9.3,a,f9.3,a)') '   trial mean ', &
        1.0d3 * (t_replace + t_solve) / ntrials, ' ms (replace ', &
        1.0d3 * t_replace / ntrials, ' ms + solve ', 1.0d3 * t_solve / ntrials, &
        ' ms), fastest trial ', 1.0d3 * t_min, ' ms'
      write(*, '(a,f6.2,a,es9.2,a,es24.16)') '   mean iterations ', &
        real(iters) / real(ntrials), ' worst eigen-residual ', worst, ' lambda ', lambda
    end if
  end subroutine ztime_trials

  subroutine zrun_timing(order, npending, ntrials)
    !Wall times of the building blocks and of full trials for a synthetic
    !complex Hermitian problem of the given order.
    integer, intent(in) :: order, npending, ntrials
    integer :: info, num_iter, i, j, k, rep, set(1), stride
    logical :: need_flush
    real(wp) :: lambda, rel_acc, eig_res, tol, lam_border, abs_res, rel_res, fwd
    real(kind(1.0d0)) :: t0, t1, t_sum, t_begin
    integer, parameter :: NREP = 10

    tol = 1.0e-12_wp
    rng_state = 987654321_int64
    model = MODEL_RITZ
    cap = order
    allocate(zhmat(order, order), zsmat(order, order), zhcol(order), zscol(order))
    allocate(zbvec(order), zxvec(order), zxref(order), zwrk1(order), zwrk2(order))
    allocate(zrw4(order))
    do j = 1, order
      do i = 1, order
        zhmat(i, j) = 7.7e33_wp
        zsmat(i, j) = -3.3e31_wp
      end do
    end do
    call ztiming_matrices(order)
    call zst%initialize(order, max(npending + 8, 40), info)
    if (info /= QR_SUCCESS) call fail('timing: initialize')
    if (myrank == 0) write(*, '(a,i0,a,i0,a,i0)') 'timing (complex Hermitian): n=', order, ' max_pending=', &
      zst%max_pending(), ' trials=', ntrials

    !Rough shift below the spectrum, then a shift 1e-6 (relative) below the
    !lowest eigenvalue, as the application chooses it.
    t0 = wall()
    call zst%factorize_fresh(zhmat, zsmat, -20.0_wp, info)
    t1 = wall()
    if (info /= QR_SUCCESS) call fail('timing: factorize_fresh')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'factorize_fresh (far shift)      ', t1 - t0, ' s'
    zwrk1(1:n) = 1.0_wp
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, 1.0e-10_wp, 500, 0, rel_acc, &
                  num_iter, eig_res, info)
    if (info /= QR_SUCCESS) call fail('timing: first solve')
    if (myrank == 0) write(*, '(a,es24.16,a,i0,a,es9.2)') '   lowest eigenvalue ', lambda, &
      ' iterations ', num_iter, ' residual ', eig_res
    shift = lambda - 1.0e-6_wp * abs(lambda)
    t0 = wall()
    call zst%factorize_fresh(zhmat, zsmat, shift, info)
    t1 = wall()
    if (info /= QR_SUCCESS) call fail('timing: factorize_fresh')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'factorize_fresh (near shift)     ', t1 - t0, ' s'

    !Base solve and symmetric product.
    do i = 1, n
      zbvec(i) = zrnd()
    end do
    t_sum = 0.0d0
    do rep = 1, NREP
      t0 = wall()
      call zst%solve_shifted(zbvec(1:n), zxref(1:n), info)
      t1 = wall()
      t_sum = t_sum + (t1 - t0)
    end do
    if (myrank == 0) write(*, '(a,f10.3,a)') 'base solve (no border)           ', &
      1.0d3 * t_sum / NREP, ' ms'
    t_sum = 0.0d0
    do rep = 1, NREP
      t0 = wall()
      call zst%symmetric_product(zsmat, zbvec(1:n), zxref(1:n), info)
      t1 = wall()
      t_sum = t_sum + (t1 - t0)
    end do
    if (myrank == 0) write(*, '(a,f10.3,a)') 'Hermitian product                ', &
      1.0d3 * t_sum / NREP, ' ms'
    call zst%factorization_residual(zhmat, zsmat, zbvec(1:n), abs_res, rel_res, info, fwd)
    if (myrank == 0) write(*, '(a,es9.2,a,es9.2)') '   factorization_residual ', rel_res, &
      ' forward error ', fwd

    !Eigenvector with the near shift as the start of the trials.
    zwrk1(1:n) = zxvec(1:n)
    t0 = wall()
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, tol, 50, 0, rel_acc, &
                  num_iter, eig_res, info)
    t1 = wall()
    if (myrank == 0) write(*, '(a,f10.3,a,i0,a,es9.2)') 'solve after fresh factorization  ', &
      1.0d3 * (t1 - t0), ' ms, iterations ', num_iter, ' residual ', eig_res

    !Trials with one pending function.
    k = n / 2
    set(1) = k
    t0 = wall()
    call zst%begin_active(set, need_flush, info)
    t1 = wall()
    if (info /= QR_SUCCESS .or. need_flush) call fail('timing: begin_active')
    if (myrank == 0) write(*, '(a,f10.3,a)') 'begin_active (new base function) ', &
      1.0d3 * (t1 - t0), ' ms'
    call ztime_trials('trials, 1 pending function', k, ntrials, 1.0e-3_wp, tol)

    !Make npending-1 further functions pending, one optimized function at a
    !time, each with a few accepted changes.
    stride = max(1, n / (npending + 1))
    t_begin = 0.0d0
    do j = 1, npending - 1
      k = 1 + mod(j * stride, n)
      if (k == n / 2) k = k + 1
      set(1) = k
      t0 = wall()
      call zst%begin_active(set, need_flush, info)
      t1 = wall()
      t_begin = t_begin + (t1 - t0)
      if (info /= QR_SUCCESS .or. need_flush) call fail('timing: begin_active (pending)')
      call ztiming_perturb(k, 1.0e-3_wp)
      call zst%replace_symmetric(k, zhcol(1:n), zscol(1:n), info)
      if (info /= QR_SUCCESS) call fail('timing: replace (pending)')
    end do
    if (npending > 1) then
      if (myrank == 0) write(*, '(a,f10.3,a)') 'begin_active (mean while filling)', &
        1.0d3 * t_begin / (npending - 1), ' ms'
      !Bring the eigenvector up to date, then time trials of the last one.
      zwrk1(1:n) = zxvec(1:n)
      call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, tol, 100, 0, rel_acc, &
                    num_iter, eig_res, info)
      if (myrank == 0) write(*, '(a,i0,a,es9.2)') '   solve after filling: iterations ', &
        num_iter, ' residual ', eig_res
      call ztime_trials('trials, full border', k, ntrials, 1.0e-3_wp, tol)
      t_sum = 0.0d0
      do rep = 1, NREP
        t0 = wall()
        call zst%solve_shifted(zbvec(1:n), zxref(1:n), info)
        t1 = wall()
        t_sum = t_sum + (t1 - t0)
      end do
      if (myrank == 0) write(*, '(a,f10.3,a)') 'bordered solve (full border)     ', &
        1.0d3 * t_sum / NREP, ' ms'
    end if

    !Agreement of the bordered representation with a fresh factorization.
    lam_border = lam_ref
    call zst%factorize_fresh(zhmat, zsmat, shift, info)
    zwrk1(1:n) = zxvec(1:n)
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxref(1:n), lambda, tol, 100, 0, rel_acc, &
                  num_iter, eig_res, info)
    !The complex phase of an eigenvector is arbitrary (and its sign changes
    !when the eigenvalue has moved across the shift): compare aligned vectors.
    call zalign_phase(zxref, zxvec)
    if (myrank == 0) then
      write(*, '(a,es24.16,a,es9.2)') 'fresh factorization of the final matrices: lambda ', &
        lambda, ' residual ', eig_res
      write(*, '(a,es9.2,a,es9.2)') '   |lambda(border)-lambda(fresh)|/|lambda| ', &
        abs(lam_border - lambda) / abs(lambda), '  max|x(border)-x(fresh)| ', &
        maxval(abs(zxvec(1:n) - zxref(1:n)))
      if (nfail == 0) then
        write(*, '(a)') 'TIMING RUN COMPLETED'
      else
        write(*, '(a,i0)') 'TIMING RUN FAILED: ', nfail
      end if
    end if
    call zst%clear()
  end subroutine zrun_timing

#ifdef QRLINALG_TESTING
  !=========================================================================
  !=========================================================================
  ! Direct tests of the fresh factorization (white-box: the components of the
  ! state types are public in a QRLINALG_TESTING build)
  !
  ! Every case factorizes the same matrix with the serial algorithm (LAPACK on
  ! rank zero) and with the distributed blocked Householder algorithm,
  ! gathers the distributed factors, and checks for each algorithm
  !   - the backward error ||A - Q*R||_F/||A||_F and the largest columnwise
  !     backward error max_c ||(A - Q*R)(:,c)||/||A(:,c)||,
  !   - the orthogonality ||Q^H*Q - I||_F,
  !     all three against FACTOR_TOL*n*eps,
  !   - that R is exactly upper triangular (and has a real diagonal in the
  !     complex case), that the replicated diagonal blocks rd equal the
  !     diagonal blocks of R, and the values of nloc, nblk and r_singular,
  !   - base solves against the dense reference, forward error bounded by
  !     FACTOR_SOLVE_TOL*n*eps*cond(A),
  ! and that the base solves of the two algorithms agree within the same
  ! conditioning-aware bound.
  !
  ! The test matrix is D*A0 with a diagonal matrix D of powers of two (the
  ! scaling is exact), where A0 = H0 - shift*S0 has elements of order one and
  ! D is either a multiple of the identity or constant on the two diagonal
  ! blocks of a block-diagonal A0. The dense reference factorizes A0. Variants:
  !   1  random symmetric (Hermitian) matrix;
  !   2  prescribed singular values from 1 to 10**(-0.8*precision), that is,
  !      condition number 1e12 in double precision;
  !   3  variant 1 times 2**(0.8*maxexponent): squares of elements overflow;
  !   4  variant 1 times 2**(0.8*minexponent): squares of elements underflow;
  !   5  second diagonal block times 2**(minexponent/2-20): the squared norms
  !      of its columns underflow, the columns themselves are far above the
  !      safe minimum;
  !   6  second diagonal block times 2**(minexponent+digits/2): its columns
  !      are below the safe minimum of the reflector generation;
  !   7  variant 1 with one zero row and column (exactly singular R);
  !   8  the zero matrix.
  !=========================================================================
  !=========================================================================

  subroutine factor_matrix(variant, norder, dscale)
    !Fills the lower triangles of hmat and smat with D*H0 and D*S0, aref with
    !the full matrix D*A0, and dscale with the diagonal of D. For the
    !nonsingular variants the reference LU factorization of A0, anorm and
    !ainv_norm are set by ref_build.
    integer, intent(in) :: variant, norder
    real(wp), intent(out) :: dscale(:)
    integer :: emax, emin, i, j, k, n1, ndig, rep
    real(wp) :: kappa, sigma, t, unorm

    n = norder
    shift = 0.37_wp
    n1 = norder
    if (variant == 5 .or. variant == 6) n1 = max(1, norder / 2 - 2)
    do j = 1, n
      do i = j, n
        smat(i, j) = 0.2_wp * (2.0_wp * rnd() - 1.0_wp)
        hmat(i, j) = 2.0_wp * rnd() - 1.0_wp
        if (i > n1 .and. j <= n1) then
          smat(i, j) = 0.0_wp
          hmat(i, j) = 0.0_wp
        end if
      end do
      smat(j, j) = smat(j, j) + 1.0_wp
    end do
    if (variant == 2) then
      !A0 = U*diag(sigma)*U^T with U a product of three random reflectors.
      kappa = 10.0_wp**(0.8_wp * real(precision(1.0_wp), wp))
      aref(1:n, 1:n) = 0.0_wp
      do k = 1, n
        sigma = kappa**(-real(k - 1, wp) / real(max(1, n - 1), wp))
        if (mod(k, 3) == 0) sigma = -sigma
        aref(k, k) = sigma
      end do
      do rep = 1, 3
        do i = 1, n
          wrk1(i) = 2.0_wp * rnd() - 1.0_wp
        end do
        unorm = sqrt(sum(wrk1(1:n)**2))
        wrk1(1:n) = wrk1(1:n) / unorm
        wrk2(1:n) = matmul(aref(1:n, 1:n), wrk1(1:n))
        t = sum(wrk1(1:n) * wrk2(1:n))
        do j = 1, n
          do i = 1, n
            aref(i, j) = aref(i, j) - 2.0_wp * wrk1(i) * wrk2(j) - &
                         2.0_wp * wrk2(i) * wrk1(j) + 4.0_wp * t * wrk1(i) * wrk1(j)
          end do
        end do
      end do
      do j = 1, n
        do i = j, n
          hmat(i, j) = 0.5_wp * (aref(i, j) + aref(j, i)) + shift * smat(i, j)
        end do
      end do
    else if (variant == 7) then
      k = (2 * n) / 5 + 1
      do j = 1, n
        do i = j, n
          if (i == k .or. j == k) then
            hmat(i, j) = 0.0_wp
            smat(i, j) = 0.0_wp
          end if
        end do
      end do
    else if (variant == 8) then
      do j = 1, n
        do i = j, n
          hmat(i, j) = 0.0_wp
          smat(i, j) = 0.0_wp
        end do
      end do
    end if

    if (variant <= 6) then
      call ref_build()
    else
      do j = 1, n
        do i = j, n
          aref(i, j) = hmat(i, j) - shift * smat(i, j)
          aref(j, i) = aref(i, j)
        end do
      end do
      anorm = 0.0_wp
      ainv_norm = 0.0_wp
    end if

    dscale(1:n) = 1.0_wp
    emax = maxexponent(1.0_wp)
    emin = minexponent(1.0_wp)
    ndig = digits(1.0_wp)
    select case (variant)
    case (3)
      dscale(1:n) = scale(1.0_wp, (4 * emax) / 5)
    case (4)
      dscale(1:n) = scale(1.0_wp, (4 * emin) / 5)
    case (5)
      dscale(n1 + 1:n) = scale(1.0_wp, emin / 2 - 20)
    case (6)
      dscale(n1 + 1:n) = scale(1.0_wp, emin + ndig / 2)
    end select
    do j = 1, n
      do i = j, n
        hmat(i, j) = dscale(i) * hmat(i, j)
        smat(i, j) = dscale(i) * smat(i, j)
      end do
      do i = 1, n
        aref(i, j) = dscale(i) * aref(i, j)
      end do
    end do
  end subroutine factor_matrix

  subroutine factor_check(label, alg, expect_singular, ratios)
    !Gathers the distributed factors of st and checks them against aref.
    !ratios: backward error, columnwise backward error and orthogonality,
    !each divided by n*eps.
    character(len=*), intent(in) :: label
    integer, intent(in) :: alg
    logical, intent(in) :: expect_singular
    real(wp), intent(out) :: ratios(3)
    real(wp), allocatable :: qf(:,:), rf(:,:)
    integer :: a, b, blk, c, i, ii, jj, k, l, nbad, nexp, nn
    real(wp) :: amax, colmax, colnorm2, d, err2, orth2, rescol, sa, sg, t
    real(wp) :: tot_a2, tot_e2
    character(len=200) :: text

    nn = st%n0
    call expect(nn == n, trim(label) // ': base order')
    allocate(qf(nn, nn), rf(nn, nn))
    qf = 0.0_wp
    rf = 0.0_wp
    do l = 1, st%nloc
      i = st%rank + 1 + (l - 1) * st%nprocs
      qf(i, 1:nn) = st%qt(1:nn, l)
      rf(i, 1:nn) = st%rt(1:nn, l)
    end do
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, qf, nn * nn, MPI_WP, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, rf, nn * nn, MPI_WP, MPI_SUM, MPI_COMM_WORLD, ierr)
    end if
    nexp = 0
    do i = 1, nn
      if (mod(i - 1, nprocs) == myrank) nexp = nexp + 1
    end do
    call expect(st%nloc == nexp, trim(label) // ': nloc')
    call expect(st%nblk == (nn + nprocs - 1) / nprocs, trim(label) // ': nblk')
    call expect(st%r_singular .eqv. expect_singular, trim(label) // ': r_singular')

    !R is exactly upper triangular; rd holds its diagonal blocks.
    nbad = 0
    do c = 1, nn
      do i = c + 1, nn
        if (.not. (rf(i, c) == 0.0_wp)) nbad = nbad + 1
      end do
    end do
    call expect(nbad == 0, trim(label) // ': R is not upper triangular')
    nbad = 0
    do blk = 1, st%nblk
      do jj = 1, min(nprocs, nn - (blk - 1) * nprocs)
        do ii = 1, jj
          if (.not. (st%rd(ii, jj, blk) == rf((blk - 1) * nprocs + ii, (blk - 1) * nprocs + jj))) &
            nbad = nbad + 1
        end do
      end do
    end do
    call expect(nbad == 0, trim(label) // ': rd differs from the diagonal blocks of R')

    !Backward errors. Differences are multiplied by powers of two before they
    !are squared, so that no square overflows or underflows.
    amax = maxval(abs(aref(1:nn, 1:nn)))
    sg = 1.0_wp
    if (amax > 0.0_wp) sg = scale(1.0_wp, -exponent(amax))
    tot_a2 = 0.0_wp
    tot_e2 = 0.0_wp
    rescol = 0.0_wp
    do c = 1, nn
      colmax = maxval(abs(aref(1:nn, c)))
      sa = 1.0_wp
      if (colmax > 0.0_wp) sa = scale(1.0_wp, -exponent(colmax))
      err2 = 0.0_wp
      colnorm2 = 0.0_wp
      do i = 1, nn
        t = 0.0_wp
        do k = 1, c
          t = t + qf(i, k) * rf(k, c)
        end do
        d = aref(i, c) - t
        err2 = err2 + (d * sa)**2
        colnorm2 = colnorm2 + (aref(i, c) * sa)**2
        tot_e2 = tot_e2 + (d * sg)**2
        tot_a2 = tot_a2 + (aref(i, c) * sg)**2
      end do
      if (colnorm2 > 0.0_wp) then
        rescol = max(rescol, sqrt(err2 / colnorm2))
      else if (err2 > 0.0_wp) then
        rescol = huge(1.0_wp)
      end if
    end do
    ratios(1) = 0.0_wp
    if (tot_a2 > 0.0_wp) then
      ratios(1) = sqrt(tot_e2 / tot_a2) / (real(nn, wp) * EPS)
    else if (tot_e2 > 0.0_wp) then
      ratios(1) = huge(1.0_wp)
    end if
    ratios(2) = rescol / (real(nn, wp) * EPS)
    orth2 = 0.0_wp
    do b = 1, nn
      do a = 1, nn
        t = 0.0_wp
        do i = 1, nn
          t = t + qf(i, a) * qf(i, b)
        end do
        if (a == b) t = t - 1.0_wp
        orth2 = orth2 + t * t
      end do
    end do
    ratios(3) = sqrt(orth2) / (real(nn, wp) * EPS)
    do k = 1, 3
      fstat(k, alg) = max(fstat(k, alg), ratios(k))
      nchecks = nchecks + 1
      if (.not. (ratios(k) <= FACTOR_TOL)) then
        write(text, '(2a,i0,a,es10.3)') trim(label), ': factorization error ', k, &
          ' /(n*eps)=', ratios(k)
        call fail(text)
      end if
    end do
    deallocate(qf, rf)
  end subroutine factor_check

  subroutine factor_case(name, variant, norder, seed)
    !One matrix, both algorithms; see the header of this section.
    character(len=*), intent(in) :: name
    integer, intent(in) :: variant, norder, seed
    integer, parameter :: NRHS = 2
    real(wp), allocatable :: dscale(:), b0(:,:), xs(:,:,:), xr(:,:)
    integer :: alg, algs(2), i, info, irhs
    real(wp) :: cond, diff, fwd(2), ratios(3, 2), xnorm
    character(len=11) :: algname(2)
    character(len=200) :: label, text

    algs(1) = QR_FACTOR_SERIAL
    algs(2) = QR_FACTOR_DISTRIBUTED
    algname(1) = 'serial'
    algname(2) = 'distributed'
    call problem_allocate(MODEL_RITZ, norder, 1)
    rng_state = int(seed, int64)
    allocate(dscale(norder), b0(norder, NRHS), xs(norder, NRHS, 2), xr(norder, NRHS))
    call factor_matrix(variant, norder, dscale)
    cond = anorm * ainv_norm
    do irhs = 1, NRHS
      do i = 1, n
        b0(i, irhs) = 2.0_wp * rnd() - 1.0_wp
      end do
      if (variant <= 6) then
        xr(1:n, irhs) = b0(1:n, irhs)
        call ref_solve(xr(:, irhs))
      end if
    end do
    call st%initialize(norder, 2, info)
    call expect_info(info, QR_SUCCESS, trim(name) // ': initialize')
    fwd = 0.0_wp
    diff = 0.0_wp
    ratios = 0.0_wp
    do alg = 1, 2
      write(label, '(4a,i0)') trim(name), ' ', trim(algname(alg)), ' n=', norder
      call qr_set_factorization(algs(alg))
      call st%factorize_fresh(hmat, smat, shift, info, n)
      call expect_info(info, QR_SUCCESS, trim(label) // ': factorize_fresh')
      if (info /= QR_SUCCESS) cycle
      if (env_factor < 0) then
        call expect(st%factorization_kind() == algs(alg), trim(label) // ': factorization_kind')
      else if (env_factor /= QR_FACTOR_AUTO) then
        !The environment variable overrides the request.
        call expect(st%factorization_kind() == env_factor, trim(label) // &
                    ': factorization_kind with ECG_QR_FACTOR')
      end if
      call factor_check(label, alg, variant >= 7, ratios(:, alg))
      do irhs = 1, NRHS
        bvec(1:n) = dscale(1:n) * b0(1:n, irhs)
        call st%solve_shifted(bvec(1:n), xvec(1:n), info)
        if (variant >= 7) then
          call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': solve_shifted')
          cycle
        end if
        call expect_info(info, QR_SUCCESS, trim(label) // ': solve_shifted')
        call expect_replicated(n, xvec, trim(label) // ': solve_shifted')
        xs(1:n, irhs, alg) = xvec(1:n)
        xnorm = sqrt(sum(xr(1:n, irhs)**2))
        fwd(alg) = max(fwd(alg), sqrt(sum((xvec(1:n) - xr(1:n, irhs))**2)) / &
                                 (xnorm * EPS * cond * real(n, wp)))
        if (alg == 2) then
          diff = max(diff, sqrt(sum((xs(1:n, irhs, 1) - xs(1:n, irhs, 2))**2)) / &
                           (xnorm * EPS * cond * real(n, wp)))
          call add_trace(sqrt(sum(xvec(1:n)**2)), &
                         FACTOR_SOLVE_TOL * real(n, wp) * EPS * cond * xnorm)
        end if
      end do
      if (variant <= 6) then
        fstat(4, alg) = max(fstat(4, alg), fwd(alg))
        nchecks = nchecks + 1
        if (.not. (fwd(alg) <= FACTOR_SOLVE_TOL)) then
          write(text, '(2a,es10.3)') trim(label), ': solve forward error/(n*eps*cond)=', fwd(alg)
          call fail(text)
        end if
      end if
    end do
    if (variant <= 6) then
      fstat_diff = max(fstat_diff, diff)
      nchecks = nchecks + 1
      if (.not. (diff <= FACTOR_SOLVE_TOL)) then
        write(text, '(2a,i0,a,es10.3)') trim(name), ' n=', norder, &
          ': solves of the two algorithms differ, /(n*eps*cond)=', diff
        call fail(text)
      end if
    end if
    if (myrank == 0) then
      write(*, '(2a,i0,a,es9.2)') trim(name), ' n=', norder, ' cond=', cond
      do alg = 1, 2
        write(*, '(3x,a11,4(a,es9.2))') algname(alg), ': res/(n*eps)=', ratios(1, alg), &
          ' col/(n*eps)=', ratios(2, alg), ' orth/(n*eps)=', ratios(3, alg), &
          ' solve/(n*eps*cond)=', fwd(alg)
      end do
      if (variant <= 6) write(*, '(3x,a,es9.2)') 'serial against distributed solve: diff/(n*eps*cond)=', diff
    end if
    call qr_set_factorization(factor_choice)
    call st%clear()
    deallocate(dscale, b0, xs, xr)
  end subroutine factor_case

  subroutine run_factor_tests()
    !Direct factorization tests of the real state and a report of what the
    !automatic choice selects.
    integer, parameter :: NSIZE = 18
    integer :: sizes(NSIZE), i, info, k, kind_used, nrun
    real(wp), allocatable :: dscale(:)
    character(len=120) :: line

    sizes = (/ 1, 2, 3, 5, 9, 13, 14, 15, 31, 32, 33, 47, 64, 65, 97, 130, 200, 301 /)
    fstat = 0.0_wp
    fstat_diff = 0.0_wp
    if (myrank == 0) write(*, '(a)') 'factorization tests (real symmetric)'
    do i = 1, NSIZE
      if (quick .and. mod(i, 3) /= 0) cycle
      if (precision(1.0_wp) > 20 .and. sizes(i) > 200) cycle
      call factor_case('random', 1, sizes(i), 41000 + i)
    end do
    call factor_case('ill-conditioned', 2, 60, 42001)
    if (.not. quick) call factor_case('ill-conditioned', 2, 130, 42002)
    call factor_case('scaled up (squares overflow)', 3, 45, 43001)
    call factor_case('scaled down (squares underflow)', 4, 45, 43002)
    call factor_case('graded blocks (small column norms)', 5, 50, 43003)
    call factor_case('graded blocks (below safe minimum)', 6, 50, 43004)
    call factor_case('zero row and column', 7, 40, 43005)
    call factor_case('zero matrix', 8, 6, 43006)
    if (myrank == 0) then
      write(*, '(a)') 'factorization tests (real symmetric): worst ratios'
      write(*, '(a,4es10.2)') '   serial:      res col orth solve ', fstat(1:4, 1)
      write(*, '(a,4es10.2)') '   distributed: res col orth solve ', fstat(1:4, 2)
      write(*, '(a,es10.2)') '   serial against distributed solve ', fstat_diff
    end if

    !What the automatic choice selects on this machine with this number of
    !processes (information only).
    call qr_set_factorization(QR_FACTOR_AUTO)
    line = ' '
    nrun = 0
    do i = 1, 8
      k = 8 * 2**(i - 1)
      if (k > 512 .or. (quick .and. k > 128)) exit
      call problem_allocate(MODEL_RITZ, k, 1)
      allocate(dscale(k))
      call factor_matrix(1, k, dscale)
      call st%initialize(k, 2, info)
      call st%factorize_fresh(hmat, smat, shift, info)
      call expect_info(info, QR_SUCCESS, 'automatic choice: factorize_fresh')
      kind_used = st%factorization_kind()
      nrun = nrun + 1
      if (kind_used == QR_FACTOR_DISTRIBUTED) then
        write(line(len_trim(line) + 1:), '(1x,i0,a)') k, ':distributed'
      else
        write(line(len_trim(line) + 1:), '(1x,i0,a)') k, ':serial'
      end if
      call st%clear()
      deallocate(dscale)
    end do
    if (myrank == 0) then
      write(*, '(2a)') 'automatic choice (real) by order:', trim(line)
      write(*, '(a,4es11.3)') '   timing model (serial n^3, n^2; distributed n^3, n): ', real_factor_model
    end if
    call qr_set_factorization(factor_choice)
    call problem_free()
  end subroutine run_factor_tests
#endif

#ifdef QRLINALG_TESTING
  !=========================================================================
  ! Direct tests of the fresh factorization, complex Hermitian state (see the
  ! header of the real tests; the matrices are complex Hermitian, and the
  ! diagonal of R must be real)
  !=========================================================================

  subroutine zfactor_matrix(variant, norder, dscale)
    !Fills the lower triangles of zhmat and zsmat with D*H0 and D*S0, zaref with
    !the full matrix D*A0, and dscale with the diagonal of D. For the
    !nonsingular variants the reference LU factorization of A0, anorm and
    !ainv_norm are set by zref_build.
    integer, intent(in) :: variant, norder
    real(wp), intent(out) :: dscale(:)
    integer :: emax, emin, i, j, k, n1, ndig, rep
    real(wp) :: kappa, sigma, theta, unorm
    complex(wp) :: t

    n = norder
    shift = 0.37_wp
    n1 = norder
    if (variant == 5 .or. variant == 6) n1 = max(1, norder / 2 - 2)
    do j = 1, n
      do i = j, n
        zsmat(i, j) = 0.2_wp * zrnd()
        zhmat(i, j) = zrnd()
        if (i > n1 .and. j <= n1) then
          zsmat(i, j) = 0.0_wp
          zhmat(i, j) = 0.0_wp
        end if
      end do
      zsmat(j, j) = cmplx(real(zsmat(j, j), wp) + 1.0_wp, 0.0_wp, kind=wp)
      zhmat(j, j) = cmplx(real(zhmat(j, j), wp), 0.0_wp, kind=wp)
    end do
    if (variant == 2) then
      !A0 = P*U*diag(sigma)*U^T*P^H with U a product of three random real
      !reflectors and P a diagonal matrix of random phases.
      kappa = 10.0_wp**(0.8_wp * real(precision(1.0_wp), wp))
      zaref(1:n, 1:n) = 0.0_wp
      do k = 1, n
        sigma = kappa**(-real(k - 1, wp) / real(max(1, n - 1), wp))
        if (mod(k, 3) == 0) sigma = -sigma
        zaref(k, k) = sigma
      end do
      do rep = 1, 3
        do i = 1, n
          zwrk1(i) = cmplx(2.0_wp * rnd() - 1.0_wp, 0.0_wp, kind=wp)
        end do
        unorm = sqrt(sum(abs(zwrk1(1:n))**2))
        zwrk1(1:n) = zwrk1(1:n) / unorm
        zwrk2(1:n) = matmul(zaref(1:n, 1:n), zwrk1(1:n))
        t = sum(zwrk1(1:n) * zwrk2(1:n))
        do j = 1, n
          do i = 1, n
            zaref(i, j) = zaref(i, j) - 2.0_wp * zwrk1(i) * zwrk2(j) - &
                         2.0_wp * zwrk2(i) * zwrk1(j) + 4.0_wp * t * zwrk1(i) * zwrk1(j)
          end do
        end do
      end do
      do i = 1, n
        theta = 8.0_wp * atan(1.0_wp) * rnd()
        zwrk1(i) = cmplx(cos(theta), sin(theta), kind=wp)
      end do
      do j = 1, n
        do i = j, n
          zhmat(i, j) = 0.5_wp * (zaref(i, j) + zaref(j, i)) * zwrk1(i) * conjg(zwrk1(j)) + &
                        shift * zsmat(i, j)
        end do
        zhmat(j, j) = cmplx(real(zhmat(j, j), wp), 0.0_wp, kind=wp)
      end do
    else if (variant == 7) then
      k = (2 * n) / 5 + 1
      do j = 1, n
        do i = j, n
          if (i == k .or. j == k) then
            zhmat(i, j) = 0.0_wp
            zsmat(i, j) = 0.0_wp
          end if
        end do
      end do
    else if (variant == 8) then
      do j = 1, n
        do i = j, n
          zhmat(i, j) = 0.0_wp
          zsmat(i, j) = 0.0_wp
        end do
      end do
    end if

    if (variant <= 6) then
      call zref_build()
    else
      do j = 1, n
        zaref(j, j) = cmplx(real(zhmat(j, j), wp) - shift * real(zsmat(j, j), wp), 0.0_wp, kind=wp)
        do i = j + 1, n
          zaref(i, j) = zhmat(i, j) - shift * zsmat(i, j)
          zaref(j, i) = conjg(zaref(i, j))
        end do
      end do
      anorm = 0.0_wp
      ainv_norm = 0.0_wp
    end if

    dscale(1:n) = 1.0_wp
    emax = maxexponent(1.0_wp)
    emin = minexponent(1.0_wp)
    ndig = digits(1.0_wp)
    select case (variant)
    case (3)
      dscale(1:n) = scale(1.0_wp, (4 * emax) / 5)
    case (4)
      dscale(1:n) = scale(1.0_wp, (4 * emin) / 5)
    case (5)
      dscale(n1 + 1:n) = scale(1.0_wp, emin / 2 - 20)
    case (6)
      dscale(n1 + 1:n) = scale(1.0_wp, emin + ndig / 2)
    end select
    do j = 1, n
      do i = j, n
        zhmat(i, j) = dscale(i) * zhmat(i, j)
        zsmat(i, j) = dscale(i) * zsmat(i, j)
      end do
      do i = 1, n
        zaref(i, j) = dscale(i) * zaref(i, j)
      end do
    end do
  end subroutine zfactor_matrix

  subroutine zfactor_check(label, alg, expect_singular, ratios)
    !Gathers the distributed factors of st and checks them against zaref.
    !ratios: backward error, columnwise backward error and orthogonality,
    !each divided by n*eps.
    character(len=*), intent(in) :: label
    integer, intent(in) :: alg
    logical, intent(in) :: expect_singular
    real(wp), intent(out) :: ratios(3)
    complex(wp), allocatable :: qf(:,:), rf(:,:)
    complex(wp) :: d, t
    integer :: a, b, blk, c, i, ii, jj, k, l, nbad, nexp, nn
    real(wp) :: amax, colmax, colnorm2, err2, orth2, rescol, sa, sg
    real(wp) :: tot_a2, tot_e2
    character(len=200) :: text

    nn = zst%n0
    call expect(nn == n, trim(label) // ': base order')
    allocate(qf(nn, nn), rf(nn, nn))
    qf = 0.0_wp
    rf = 0.0_wp
    do l = 1, zst%nloc
      i = zst%rank + 1 + (l - 1) * zst%nprocs
      !The state stores the complex conjugates of the rows of Q.
      qf(i, 1:nn) = conjg(zst%qt(1:nn, l))
      rf(i, 1:nn) = zst%rt(1:nn, l)
    end do
    if (nprocs > 1) then
      call MPI_Allreduce(MPI_IN_PLACE, qf, 2 * nn * nn, MPI_WP, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, rf, 2 * nn * nn, MPI_WP, MPI_SUM, MPI_COMM_WORLD, ierr)
    end if
    nexp = 0
    do i = 1, nn
      if (mod(i - 1, nprocs) == myrank) nexp = nexp + 1
    end do
    call expect(zst%nloc == nexp, trim(label) // ': nloc')
    call expect(zst%nblk == (nn + nprocs - 1) / nprocs, trim(label) // ': nblk')
    call expect(zst%r_singular .eqv. expect_singular, trim(label) // ': r_singular')

    !R is exactly upper triangular; rd holds its diagonal blocks.
    nbad = 0
    do c = 1, nn
      do i = c + 1, nn
        if (.not. (rf(i, c) == 0.0_wp)) nbad = nbad + 1
      end do
    end do
    call expect(nbad == 0, trim(label) // ': R is not upper triangular')
    nbad = 0
    do c = 1, nn
      if (.not. (aimag(rf(c, c)) == 0.0_wp)) nbad = nbad + 1
    end do
    call expect(nbad == 0, trim(label) // ': the diagonal of R is not real')
    nbad = 0
    do blk = 1, zst%nblk
      do jj = 1, min(nprocs, nn - (blk - 1) * nprocs)
        do ii = 1, jj
          if (.not. (zst%rd(ii, jj, blk) == rf((blk - 1) * nprocs + ii, (blk - 1) * nprocs + jj))) &
            nbad = nbad + 1
        end do
      end do
    end do
    call expect(nbad == 0, trim(label) // ': rd differs from the diagonal blocks of R')

    !Backward errors. Differences are multiplied by powers of two before they
    !are squared, so that no square overflows or underflows.
    amax = maxval(abs(zaref(1:nn, 1:nn)))
    sg = 1.0_wp
    if (amax > 0.0_wp) sg = scale(1.0_wp, -exponent(amax))
    tot_a2 = 0.0_wp
    tot_e2 = 0.0_wp
    rescol = 0.0_wp
    do c = 1, nn
      colmax = maxval(abs(zaref(1:nn, c)))
      sa = 1.0_wp
      if (colmax > 0.0_wp) sa = scale(1.0_wp, -exponent(colmax))
      err2 = 0.0_wp
      colnorm2 = 0.0_wp
      do i = 1, nn
        t = 0.0_wp
        do k = 1, c
          t = t + qf(i, k) * rf(k, c)
        end do
        d = zaref(i, c) - t
        err2 = err2 + real(d * sa, wp)**2 + aimag(d * sa)**2
        colnorm2 = colnorm2 + real(zaref(i, c) * sa, wp)**2 + aimag(zaref(i, c) * sa)**2
        tot_e2 = tot_e2 + real(d * sg, wp)**2 + aimag(d * sg)**2
        tot_a2 = tot_a2 + real(zaref(i, c) * sg, wp)**2 + aimag(zaref(i, c) * sg)**2
      end do
      if (colnorm2 > 0.0_wp) then
        rescol = max(rescol, sqrt(err2 / colnorm2))
      else if (err2 > 0.0_wp) then
        rescol = huge(1.0_wp)
      end if
    end do
    ratios(1) = 0.0_wp
    if (tot_a2 > 0.0_wp) then
      ratios(1) = sqrt(tot_e2 / tot_a2) / (real(nn, wp) * EPS)
    else if (tot_e2 > 0.0_wp) then
      ratios(1) = huge(1.0_wp)
    end if
    ratios(2) = rescol / (real(nn, wp) * EPS)
    orth2 = 0.0_wp
    do b = 1, nn
      do a = 1, nn
        t = 0.0_wp
        do i = 1, nn
          t = t + conjg(qf(i, a)) * qf(i, b)
        end do
        if (a == b) t = t - 1.0_wp
        orth2 = orth2 + real(t, wp)**2 + aimag(t)**2
      end do
    end do
    ratios(3) = sqrt(orth2) / (real(nn, wp) * EPS)
    do k = 1, 3
      fstat(k, alg) = max(fstat(k, alg), ratios(k))
      nchecks = nchecks + 1
      if (.not. (ratios(k) <= FACTOR_TOL)) then
        write(text, '(2a,i0,a,es10.3)') trim(label), ': factorization error ', k, &
          ' /(n*eps)=', ratios(k)
        call fail(text)
      end if
    end do
    deallocate(qf, rf)
  end subroutine zfactor_check

  subroutine zfactor_case(name, variant, norder, seed)
    !One matrix, both algorithms; see the header of this section.
    character(len=*), intent(in) :: name
    integer, intent(in) :: variant, norder, seed
    integer, parameter :: NRHS = 2
    real(wp), allocatable :: dscale(:)
    complex(wp), allocatable :: b0(:,:), xs(:,:,:), xr(:,:)
    integer :: alg, algs(2), i, info, irhs
    real(wp) :: cond, diff, fwd(2), ratios(3, 2), xnorm
    character(len=11) :: algname(2)
    character(len=200) :: label, text

    algs(1) = QR_FACTOR_SERIAL
    algs(2) = QR_FACTOR_DISTRIBUTED
    algname(1) = 'serial'
    algname(2) = 'distributed'
    call zproblem_allocate(MODEL_RITZ, norder, 1)
    rng_state = int(seed, int64)
    allocate(dscale(norder), b0(norder, NRHS), xs(norder, NRHS, 2), xr(norder, NRHS))
    call zfactor_matrix(variant, norder, dscale)
    cond = anorm * ainv_norm
    do irhs = 1, NRHS
      do i = 1, n
        b0(i, irhs) = zrnd()
      end do
      if (variant <= 6) then
        xr(1:n, irhs) = b0(1:n, irhs)
        call zref_solve(xr(:, irhs))
      end if
    end do
    call zst%initialize(norder, 2, info)
    call expect_info(info, QR_SUCCESS, trim(name) // ': initialize')
    fwd = 0.0_wp
    diff = 0.0_wp
    ratios = 0.0_wp
    do alg = 1, 2
      write(label, '(4a,i0)') trim(name), ' ', trim(algname(alg)), ' n=', norder
      call qr_set_factorization(algs(alg))
      call zst%factorize_fresh(zhmat, zsmat, shift, info, n)
      call expect_info(info, QR_SUCCESS, trim(label) // ': factorize_fresh')
      if (info /= QR_SUCCESS) cycle
      if (env_factor < 0) then
        call expect(zst%factorization_kind() == algs(alg), trim(label) // ': factorization_kind')
      else if (env_factor /= QR_FACTOR_AUTO) then
        !The environment variable overrides the request.
        call expect(zst%factorization_kind() == env_factor, trim(label) // &
                    ': factorization_kind with ECG_QR_FACTOR')
      end if
      call zfactor_check(label, alg, variant >= 7, ratios(:, alg))
      do irhs = 1, NRHS
        zbvec(1:n) = dscale(1:n) * b0(1:n, irhs)
        call zst%solve_shifted(zbvec(1:n), zxvec(1:n), info)
        if (variant >= 7) then
          call expect_info(info, QR_ERR_SINGULAR, trim(label) // ': solve_shifted')
          cycle
        end if
        call expect_info(info, QR_SUCCESS, trim(label) // ': solve_shifted')
        call zexpect_replicated(n, zxvec, trim(label) // ': solve_shifted')
        xs(1:n, irhs, alg) = zxvec(1:n)
        xnorm = sqrt(sum(abs(xr(1:n, irhs))**2))
        fwd(alg) = max(fwd(alg), sqrt(sum(abs(zxvec(1:n) - xr(1:n, irhs))**2)) / &
                                 (xnorm * EPS * cond * real(n, wp)))
        if (alg == 2) then
          diff = max(diff, sqrt(sum(abs(xs(1:n, irhs, 1) - xs(1:n, irhs, 2))**2)) / &
                           (xnorm * EPS * cond * real(n, wp)))
          call add_trace(sqrt(sum(abs(zxvec(1:n))**2)), &
                         FACTOR_SOLVE_TOL * real(n, wp) * EPS * cond * xnorm)
        end if
      end do
      if (variant <= 6) then
        fstat(4, alg) = max(fstat(4, alg), fwd(alg))
        nchecks = nchecks + 1
        if (.not. (fwd(alg) <= FACTOR_SOLVE_TOL)) then
          write(text, '(2a,es10.3)') trim(label), ': solve forward error/(n*eps*cond)=', fwd(alg)
          call fail(text)
        end if
      end if
    end do
    if (variant <= 6) then
      fstat_diff = max(fstat_diff, diff)
      nchecks = nchecks + 1
      if (.not. (diff <= FACTOR_SOLVE_TOL)) then
        write(text, '(2a,i0,a,es10.3)') trim(name), ' n=', norder, &
          ': solves of the two algorithms differ, /(n*eps*cond)=', diff
        call fail(text)
      end if
    end if
    if (myrank == 0) then
      write(*, '(2a,i0,a,es9.2)') trim(name), ' n=', norder, ' cond=', cond
      do alg = 1, 2
        write(*, '(3x,a11,4(a,es9.2))') algname(alg), ': res/(n*eps)=', ratios(1, alg), &
          ' col/(n*eps)=', ratios(2, alg), ' orth/(n*eps)=', ratios(3, alg), &
          ' solve/(n*eps*cond)=', fwd(alg)
      end do
      if (variant <= 6) write(*, '(3x,a,es9.2)') 'serial against distributed solve: diff/(n*eps*cond)=', diff
    end if
    call qr_set_factorization(factor_choice)
    call zst%clear()
    deallocate(dscale, b0, xs, xr)
  end subroutine zfactor_case

  subroutine zrun_factor_tests()
    !Direct factorization tests of the complex state and a report of what the
    !automatic choice selects.
    integer, parameter :: NSIZE = 18
    integer :: sizes(NSIZE), i, info, k, kind_used, nrun
    real(wp), allocatable :: dscale(:)
    character(len=120) :: line

    sizes = (/ 1, 2, 3, 5, 9, 13, 14, 15, 31, 32, 33, 47, 64, 65, 97, 130, 200, 301 /)
    fstat = 0.0_wp
    fstat_diff = 0.0_wp
    if (myrank == 0) write(*, '(a)') 'factorization tests (complex Hermitian)'
    do i = 1, NSIZE
      if (quick .and. mod(i, 3) /= 0) cycle
      if (precision(1.0_wp) > 20 .and. sizes(i) > 200) cycle
      call zfactor_case('random', 1, sizes(i), 51000 + i)
    end do
    call zfactor_case('ill-conditioned', 2, 60, 52001)
    if (.not. quick) call zfactor_case('ill-conditioned', 2, 130, 52002)
    call zfactor_case('scaled up (squares overflow)', 3, 45, 53001)
    call zfactor_case('scaled down (squares underflow)', 4, 45, 53002)
    call zfactor_case('graded blocks (small column norms)', 5, 50, 53003)
    call zfactor_case('graded blocks (below safe minimum)', 6, 50, 53004)
    call zfactor_case('zero row and column', 7, 40, 53005)
    call zfactor_case('zero matrix', 8, 6, 53006)
    if (myrank == 0) then
      write(*, '(a)') 'factorization tests (complex Hermitian): worst ratios'
      write(*, '(a,4es10.2)') '   serial:      res col orth solve ', fstat(1:4, 1)
      write(*, '(a,4es10.2)') '   distributed: res col orth solve ', fstat(1:4, 2)
      write(*, '(a,es10.2)') '   serial against distributed solve ', fstat_diff
    end if

    !What the automatic choice selects on this machine with this number of
    !processes (information only).
    call qr_set_factorization(QR_FACTOR_AUTO)
    line = ' '
    nrun = 0
    do i = 1, 8
      k = 8 * 2**(i - 1)
      if (k > 512 .or. (quick .and. k > 128)) exit
      call zproblem_allocate(MODEL_RITZ, k, 1)
      allocate(dscale(k))
      call zfactor_matrix(1, k, dscale)
      call zst%initialize(k, 2, info)
      call zst%factorize_fresh(zhmat, zsmat, shift, info)
      call expect_info(info, QR_SUCCESS, 'automatic choice: factorize_fresh')
      kind_used = zst%factorization_kind()
      nrun = nrun + 1
      if (kind_used == QR_FACTOR_DISTRIBUTED) then
        write(line(len_trim(line) + 1:), '(1x,i0,a)') k, ':distributed'
      else
        write(line(len_trim(line) + 1:), '(1x,i0,a)') k, ':serial'
      end if
      call zst%clear()
      deallocate(dscale)
    end do
    if (myrank == 0) then
      write(*, '(2a)') 'automatic choice (complex) by order:', trim(line)
      write(*, '(a,4es11.3)') '   timing model (serial n^3, n^2; distributed n^3, n): ', complex_factor_model
    end if
    call qr_set_factorization(factor_choice)
    call zproblem_free()
  end subroutine zrun_factor_tests
#endif

  !=========================================================================
  ! Timing of the fresh factorization (mode factortime)
  !=========================================================================

  subroutine run_factor_timing(order, reps)
    !Wall times of factorize_fresh with the serial and with the distributed
    !algorithm for the synthetic problem of the timing mode, the agreement of
    !the two base solves, and what the automatic choice selects.
    integer, intent(in) :: order, reps
    integer :: alg, algs(2), i, info, j, nrep, num_iter, rep
    real(wp) :: abs_res, eig_res, fwd, lambda, rel_acc, rel_res(2), xdiff
    real(kind(1.0d0)) :: t0, t1, tbest(2)
    character(len=11) :: algname(2)

    algs(1) = QR_FACTOR_SERIAL
    algs(2) = QR_FACTOR_DISTRIBUTED
    algname(1) = 'serial'
    algname(2) = 'distributed'
    rng_state = 987654321_int64
    model = MODEL_RITZ
    cap = order
    allocate(hmat(order, order), smat(order, order), hcol(order), scol(order))
    allocate(bvec(order), xvec(order), xref(order), wrk1(order), wrk2(order))
    do j = 1, order
      do i = 1, order
        hmat(i, j) = 7.7e33_wp
        smat(i, j) = -3.3e31_wp
      end do
    end do
    call timing_matrices(order)
    call st%initialize(order, 8, info)
    if (info /= QR_SUCCESS) call fail('factortime: initialize')
    if (myrank == 0) write(*, '(a,i0,a,i0)') 'factortime (real symmetric): n=', order, &
      ' processes=', nprocs

    !Shift 1e-6 (relative) below the lowest eigenvalue, as the application
    !chooses it; the preliminary factorization uses the distributed algorithm.
    call qr_set_factorization(QR_FACTOR_DISTRIBUTED)
    call st%factorize_fresh(hmat, smat, -20.0_wp, info)
    if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh')
    wrk1(1:n) = 1.0_wp
    call st%solve(hmat, smat, wrk1(1:n), xvec(1:n), lambda, 1.0e-10_wp, 500, 0, rel_acc, &
                  num_iter, eig_res, info)
    if (info /= QR_SUCCESS) call fail('factortime: first solve')
    shift = lambda - 1.0e-6_wp * abs(lambda)
    do i = 1, n
      bvec(i) = 2.0_wp * rnd() - 1.0_wp
    end do

    do alg = 1, 2
      if (alg == 1 .and. skip_serial) cycle
      call qr_set_factorization(algs(alg))
      nrep = reps
      if (alg == 1) nrep = 1
      tbest(alg) = huge(1.0d0)
      do rep = 1, nrep
        t0 = wall()
        call st%factorize_fresh(hmat, smat, shift, info)
        t1 = wall()
        if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh')
        tbest(alg) = min(tbest(alg), t1 - t0)
        if (myrank == 0) write(*, '(3a,i0,a,f10.4,a)') 'factorize_fresh ', algname(alg), &
          ' run ', rep, ': ', t1 - t0, ' s'
      end do
      if (alg == 1) then
        call st%solve_shifted(bvec(1:n), xref(1:n), info)
      else
        call st%solve_shifted(bvec(1:n), xvec(1:n), info)
      end if
      if (info /= QR_SUCCESS) call fail('factortime: solve_shifted')
      call st%factorization_residual(hmat, smat, bvec(1:n), abs_res, rel_res(alg), info, fwd)
    end do
    if (skip_serial) then
      tbest(1) = 0.0d0
      rel_res(1) = 0.0_wp
      xref(1:n) = xvec(1:n)
    end if
    xdiff = maxval(abs(xvec(1:n) - xref(1:n))) / maxval(abs(xref(1:n)))
    if (myrank == 0) then
      write(*, '(a,i0,a,i0,a,f10.4,a,f10.4,a,f8.2)') 'FACTOR n=', order, ' P=', nprocs, &
        ' serial ', tbest(1), ' s distributed ', tbest(2), ' s speedup ', tbest(1) / tbest(2)
      write(*, '(a,2es10.2,a,es10.2)') '   factorization_residual (serial, distributed) ', &
        rel_res, '  max|x_serial - x_distributed|/max|x| ', xdiff
    end if

    !Automatic choice: the first call includes the calibration.
    call qr_set_factorization(QR_FACTOR_AUTO)
    do rep = 1, 2
      t0 = wall()
      call st%factorize_fresh(hmat, smat, shift, info)
      t1 = wall()
      if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh (automatic)')
      if (myrank == 0) then
        if (st%factorization_kind() == QR_FACTOR_DISTRIBUTED) then
          write(*, '(a,i0,a,f10.4,a)') 'automatic choice, call ', rep, ': distributed, ', t1 - t0, ' s'
        else
          write(*, '(a,i0,a,f10.4,a)') 'automatic choice, call ', rep, ': serial, ', t1 - t0, ' s'
        end if
      end if
    end do
#ifdef QRLINALG_TESTING
    if (myrank == 0 .and. real_factor_model(1) >= 0.0d0) then
      write(*, '(a,4es11.3)') '   timing model (serial n^3, n^2; distributed n^3, n): ', real_factor_model
      write(*, '(a,f10.4,a,f10.4,a)') '   model prediction for this order: serial ', &
        real_factor_model(1) * dble(order)**3 + real_factor_model(2) * dble(order)**2, &
        ' s, distributed ', &
        real_factor_model(3) * dble(order)**3 + real_factor_model(4) * dble(order), ' s'
    end if
#endif
    if (myrank == 0) then
      if (nfail == 0) then
        write(*, '(a)') 'TIMING RUN COMPLETED'
      else
        write(*, '(a,i0)') 'TIMING RUN FAILED: ', nfail
      end if
    end if
    call st%clear()
    deallocate(hmat, smat, hcol, scol, bvec, xvec, xref, wrk1, wrk2)
  end subroutine run_factor_timing

  subroutine zrun_factor_timing(order, reps)
    !Wall times of factorize_fresh with the serial and with the distributed
    !algorithm for the synthetic complex problem of the timing mode, the agreement of
    !the two base solves, and what the automatic choice selects.
    integer, intent(in) :: order, reps
    integer :: alg, algs(2), i, info, j, nrep, num_iter, rep
    real(wp) :: abs_res, eig_res, fwd, lambda, rel_acc, rel_res(2), xdiff
    real(kind(1.0d0)) :: t0, t1, tbest(2)
    character(len=11) :: algname(2)

    algs(1) = QR_FACTOR_SERIAL
    algs(2) = QR_FACTOR_DISTRIBUTED
    algname(1) = 'serial'
    algname(2) = 'distributed'
    rng_state = 987654321_int64
    model = MODEL_RITZ
    cap = order
    allocate(zhmat(order, order), zsmat(order, order), zhcol(order), zscol(order))
    allocate(zbvec(order), zxvec(order), zxref(order), zwrk1(order), zwrk2(order))
    allocate(zrw4(order))
    do j = 1, order
      do i = 1, order
        zhmat(i, j) = 7.7e33_wp
        zsmat(i, j) = -3.3e31_wp
      end do
    end do
    call ztiming_matrices(order)
    call zst%initialize(order, 8, info)
    if (info /= QR_SUCCESS) call fail('factortime: initialize')
    if (myrank == 0) write(*, '(a,i0,a,i0)') 'factortime (complex Hermitian): n=', order, &
      ' processes=', nprocs

    !Shift 1e-6 (relative) below the lowest eigenvalue, as the application
    !chooses it; the preliminary factorization uses the distributed algorithm.
    call qr_set_factorization(QR_FACTOR_DISTRIBUTED)
    call zst%factorize_fresh(zhmat, zsmat, -20.0_wp, info)
    if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh')
    zwrk1(1:n) = 1.0_wp
    call zst%solve(zhmat, zsmat, zwrk1(1:n), zxvec(1:n), lambda, 1.0e-10_wp, 500, 0, rel_acc, &
                  num_iter, eig_res, info)
    if (info /= QR_SUCCESS) call fail('factortime: first solve')
    shift = lambda - 1.0e-6_wp * abs(lambda)
    do i = 1, n
      zbvec(i) = zrnd()
    end do

    do alg = 1, 2
      if (alg == 1 .and. skip_serial) cycle
      call qr_set_factorization(algs(alg))
      nrep = reps
      if (alg == 1) nrep = 1
      tbest(alg) = huge(1.0d0)
      do rep = 1, nrep
        t0 = wall()
        call zst%factorize_fresh(zhmat, zsmat, shift, info)
        t1 = wall()
        if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh')
        tbest(alg) = min(tbest(alg), t1 - t0)
        if (myrank == 0) write(*, '(3a,i0,a,f10.4,a)') 'factorize_fresh ', algname(alg), &
          ' run ', rep, ': ', t1 - t0, ' s'
      end do
      if (alg == 1) then
        call zst%solve_shifted(zbvec(1:n), zxref(1:n), info)
      else
        call zst%solve_shifted(zbvec(1:n), zxvec(1:n), info)
      end if
      if (info /= QR_SUCCESS) call fail('factortime: solve_shifted')
      call zst%factorization_residual(zhmat, zsmat, zbvec(1:n), abs_res, rel_res(alg), info, fwd)
    end do
    if (skip_serial) then
      tbest(1) = 0.0d0
      rel_res(1) = 0.0_wp
      zxref(1:n) = zxvec(1:n)
    end if
    xdiff = maxval(abs(zxvec(1:n) - zxref(1:n))) / maxval(abs(zxref(1:n)))
    if (myrank == 0) then
      write(*, '(a,i0,a,i0,a,f10.4,a,f10.4,a,f8.2)') 'FACTOR n=', order, ' P=', nprocs, &
        ' serial ', tbest(1), ' s distributed ', tbest(2), ' s speedup ', tbest(1) / tbest(2)
      write(*, '(a,2es10.2,a,es10.2)') '   factorization_residual (serial, distributed) ', &
        rel_res, '  max|x_serial - x_distributed|/max|x| ', xdiff
    end if

    !Automatic choice: the first call includes the calibration.
    call qr_set_factorization(QR_FACTOR_AUTO)
    do rep = 1, 2
      t0 = wall()
      call zst%factorize_fresh(zhmat, zsmat, shift, info)
      t1 = wall()
      if (info /= QR_SUCCESS) call fail('factortime: factorize_fresh (automatic)')
      if (myrank == 0) then
        if (zst%factorization_kind() == QR_FACTOR_DISTRIBUTED) then
          write(*, '(a,i0,a,f10.4,a)') 'automatic choice, call ', rep, ': distributed, ', t1 - t0, ' s'
        else
          write(*, '(a,i0,a,f10.4,a)') 'automatic choice, call ', rep, ': serial, ', t1 - t0, ' s'
        end if
      end if
    end do
#ifdef QRLINALG_TESTING
    if (myrank == 0 .and. complex_factor_model(1) >= 0.0d0) then
      write(*, '(a,4es11.3)') '   timing model (serial n^3, n^2; distributed n^3, n): ', complex_factor_model
      write(*, '(a,f10.4,a,f10.4,a)') '   model prediction for this order: serial ', &
        complex_factor_model(1) * dble(order)**3 + complex_factor_model(2) * dble(order)**2, &
        ' s, distributed ', &
        complex_factor_model(3) * dble(order)**3 + complex_factor_model(4) * dble(order), ' s'
    end if
#endif
    if (myrank == 0) then
      if (nfail == 0) then
        write(*, '(a)') 'TIMING RUN COMPLETED'
      else
        write(*, '(a,i0)') 'TIMING RUN FAILED: ', nfail
      end if
    end if
    call zst%clear()
    deallocate(zhmat, zsmat, zhcol, zscol, zbvec, zxvec, zxref, zwrk1, zwrk2, zrw4)
  end subroutine zrun_factor_timing

end program test_qrlinalg
