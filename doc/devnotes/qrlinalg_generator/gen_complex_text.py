# Text blocks inserted by gen_complex.py into the declaration part of module
# qrlinalg. This file is executed by gen_complex.py; it is not a program.

HEADER_TEXT = '''\
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
'''

CONSTANT_TEXT = '''\
  real(wp), parameter :: QR_DIAGONAL_TOLERANCE = 100.0_wp * epsilon(1.0_wp)
    !Complex state only: largest accepted imaginary part of a diagonal element
    !of H or S, relative to the scale of the supplied matrix or column.
'''

MODEL_TEXT = '''\
  real(kind(1.0d0)), save :: complex_factor_model(4) = -1.0d0
    !The same timing models for the complex state (complex_factor_select).
'''

INTERFACE_TEXT = '''\
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
'''

TAIL_BANNER = '''\
  !===========================================================================
  !===========================================================================
  ! Complex Hermitian twin: type qr_complex_state
  !
  ! Everything below mirrors the real procedures above one to one (same
  ! order, same structure, names complex_* instead of real_*). The conventions
  ! that differ from the real case are listed in the module header.
  !===========================================================================
  !===========================================================================
'''

_COMPLEX_COMPONENTS = ('qt|rt|rd|bt|gs|zb|ab|cbb|bg|bz|sc|scc_lu|sig_lu|xca|'
                       'wn|wsm|wloc|wblk|wred')

TYPE_EDITS = [
    (re.compile(r'^    real\(wp\), allocatable :: (%s)\(' % _COMPLEX_COMPONENTS,
                re.M), r'    complex(wp), allocatable :: \1(', 19),
    ('    ! qt(1:n0,l) is row i of Q and rt(1:n0,l) is row i of R, where\n'
     '    ! i = rank+1+(l-1)*P. Leading dimension capacity.\n',
     '    ! qt(1:n0,l) is the complex conjugate of row i of Q (column i of Q^H)\n'
     '    ! and rt(1:n0,l) is row i of R, where i = rank+1+(l-1)*P. Leading\n'
     '    ! dimension capacity.\n', 1),
    ('    ! cbb: current border-border values (slot by slot); bg = Bt^T*Gs (slot by\n'
     '    ! slot); bz = Bt^T*Zb (slot by D0 position).\n',
     '    ! cbb: current border-border values (slot by slot), a Hermitian matrix\n'
     '    ! whose element (s1,s2) is the matrix element (bidx of s1, bidx of s2);\n'
     '    ! bg = Bt^T*Gs (slot by slot); bz = Bt^T*Zb (slot by D0 position).\n', 1),
    ('    ! (tail order); scc_mag is the magnitude bookkeeping of the pivot test;\n',
     '    ! (tail order); scc_mag is the real magnitude bookkeeping of the pivot test;\n', 1),
]
