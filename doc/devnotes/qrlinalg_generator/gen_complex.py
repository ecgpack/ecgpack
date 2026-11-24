#!/usr/bin/env python3
"""Generate the complex Hermitian twin of module qrlinalg from its real part.

    python3 gen_complex.py --check MODULE
    python3 gen_complex.py --regenerate MODULE OUTPUT
    python3 gen_complex.py --extract MODULE REAL_OUTPUT
    python3 gen_complex.py REAL_SOURCE OUTPUT

MODULE is an installed src/qrlinalg.f90, which holds the real state type and
the generated complex one. REAL_SOURCE is a qrlinalg.f90 that contains only
the real state type. A generated module is its real source with pure
insertions: no line of the real source is changed or removed, so

    diff REAL_SOURCE OUTPUT | grep -c '^<'      prints 0.

Inserted are a header paragraph, two named entities of the declaration part,
two public statements, the interfaces of ZGEQRF and ZUNGQR, type
qr_complex_state, and, after the last real procedure, a banner and the
complex_* procedures. Every complex procedure is produced from its real
counterpart by

  1. procedure-specific substitutions (types of the declarations, conjugate
     transposes, Hermitian fills, MPI counts); each one must match an exact
     number of times, otherwise the script stops. A change of the real
     source that touches such a place is therefore noticed, not silently
     mistranslated;
  2. generic renaming (real_* -> complex_*, DGEQRF -> ZGEQRF, ...);
  3. generic rewording of comments (symmetric -> Hermitian, ^T -> ^H, ...).

Because the insertions are found again by their position and, for the
declaration part, by their exact text (gen_complex_text.py), the real source
can be recovered from a module, whether or not its complex part is up to
date:

  --extract     writes the real source contained in MODULE to REAL_OUTPUT.
  --regenerate  extracts the real source of MODULE and writes the module
                generated from it to OUTPUT (which may be MODULE itself).
                This is the normal way to work: edit the real part of
                src/qrlinalg.f90, regenerate, review the diff of the complex
                part.
  --check       extracts, generates in memory, and compares with MODULE.
                Exit status 0 when MODULE is byte for byte what the generator
                produces from its own real part, 1 otherwise.

The last form generates from a separate real source, as the module was
developed. See README.md in this directory.
"""
import re
import sys

USAGE = ('usage: gen_complex.py --check MODULE\n'
         '       gen_complex.py --regenerate MODULE OUTPUT\n'
         '       gen_complex.py --extract MODULE REAL_OUTPUT\n'
         '       gen_complex.py REAL_SOURCE OUTPUT')

ARGS = sys.argv[1:]
MODE = 'generate'
if ARGS and ARGS[0] in ('--check', '--regenerate', '--extract'):
    MODE = ARGS.pop(0)[2:]
if len(ARGS) != (1 if MODE == 'check' else 2) or ARGS[0].startswith('-'):
    raise SystemExit(USAGE)

# --------------------------------------------------------------------------
# Text inserted into the declaration part
# --------------------------------------------------------------------------

exec(open(__file__.replace('gen_complex.py', 'gen_complex_text.py')).read())


def text_lines(text):
    return text.rstrip('\n').split('\n')


def count(lines):
    """Number of lines of a file given as the list text.split('\\n')."""
    return len(lines) - (1 if lines and lines[-1] == '' else 0)


# --------------------------------------------------------------------------
# Recovery of the real source from a generated module
# --------------------------------------------------------------------------

def extract_real(lines):
    """Removes from the lines of a module everything that the assembly at the
    end of this script inserts and returns the lines of the real source. The
    blocks of the declaration part must be present with the exact text of
    gen_complex_text.py; the complex type and the complex procedures are
    delimited by their first and last lines, so their content may be out of
    date."""
    lines = list(lines)

    def locate(block, what):
        hits = [i for i in range(len(lines) - len(block) + 1)
                if lines[i:i + len(block)] == block]
        if len(hits) != 1:
            raise SystemExit('cannot extract the real source: %s found %d times '
                             '(expected once). If gen_complex_text.py was '
                             'changed, extract with the unchanged script first.'
                             % (what, len(hits)))
        return hits[0]

    def remove(block, what, before=None, after=None):
        i = locate(block, what)
        if before is not None and lines[i + len(block)] != before:
            raise SystemExit('cannot extract the real source: %s is not followed '
                             'by %r' % (what, before))
        if after is not None and lines[i - 1] != after:
            raise SystemExit('cannot extract the real source: %s is not preceded '
                             'by %r' % (what, after))
        del lines[i:i + len(block)]

    remove(text_lines(HEADER_TEXT), 'the header paragraph', before='module qrlinalg')
    remove(text_lines(CONSTANT_TEXT), 'the constant of the complex state')
    remove(text_lines(MODEL_TEXT), 'the timing model of the complex state')
    remove(['  public :: qr_complex_state'], 'public :: qr_complex_state',
           after='  public :: qr_real_state')
    remove(['  public :: complex_factor_model'], 'public :: complex_factor_model',
           after='  public :: real_factor_model')
    remove([''] + text_lines(INTERFACE_TEXT), 'the ZGEQRF/ZUNGQR interfaces',
           after='  end interface')
    i0 = locate(['  end type qr_real_state', '', '  type :: qr_complex_state'],
                'the start of type qr_complex_state') + 1
    i1 = locate(['  end type qr_complex_state'], 'the end of type qr_complex_state')
    del lines[i0:i1 + 1]
    i0 = locate(text_lines(TAIL_BANNER), 'the banner of the complex procedures')
    i1 = locate(['end module qrlinalg'], 'end module qrlinalg')
    del lines[i0:i1]
    return lines


if MODE == 'generate':
    REAL, OUT = ARGS
    src = open(REAL).read().split('\n')
else:
    MODULE = ARGS[0]
    OUT = ARGS[1] if MODE != 'check' else None
    module_lines = open(MODULE).read().split('\n')
    src = extract_real(module_lines)
    if MODE == 'extract':
        open(OUT, 'w', newline='\n').write('\n'.join(src))
        print('wrote %s: %d real source lines (of %d module lines)'
              % (OUT, count(src), count(module_lines)))
        # The generation below still runs, without writing anything: it
        # reports whether the complex part of MODULE is up to date.


def find(line_text, start=0):
    for i in range(start, len(src)):
        if src[i] == line_text:
            return i
    raise SystemExit('marker line not found: %r' % line_text)


# --------------------------------------------------------------------------
# Generic transformations
# --------------------------------------------------------------------------

def is_comment(line):
    return line.lstrip().startswith('!')


RENAMES = [
    (re.compile(r'\bclear_real_state\b'), 'clear_complex_state'),
    (re.compile(r'\bqr_real_state\b'), 'qr_complex_state'),
    (re.compile(r'\breal_dot\b'), 'complex_dotc'),
    (re.compile(r'\breal_(?=[a-z])'), 'complex_'),
    (re.compile(r'\bdgeqrf\b'), 'zgeqrf'),
    (re.compile(r'\bdorgqr\b'), 'zungqr'),
    (re.compile(r'\bDGEQRF\b'), 'ZGEQRF'),
    (re.compile(r'\bDORGQR\b'), 'ZUNGQR'),
    (re.compile(r'\bDLARFG\b'), 'ZLARFG'),
    (re.compile(r'\bDLARFT\b'), 'ZLARFT'),
]

COMMENT_WORDS = [
    (re.compile(r'\^T\b'), '^H'),
    (re.compile(r'\breal QR\b'), 'complex QR'),
    (re.compile(r'\bmirrored row contribution'), 'conjugated row contribution'),
    (re.compile(r'\breal generalized symmetric\b'), 'complex generalized Hermitian'),
    (re.compile(r'\breal symmetric\b'), 'complex Hermitian'),
    (re.compile(r'(?<![\w])symmetric(?![\w])'), 'Hermitian'),
    (re.compile(r'\borthogonal\b'), 'unitary'),
    (re.compile(r'\bReal array\b'), 'Complex array'),
    (re.compile(r'\bReal vector\b'), 'Complex vector'),
    (re.compile(r'\bcompact real QR\b'), 'compact complex QR'),
]


def rename_line(line):
    """Returns the renamed line and the list of (column, growth) of the
    substitutions, for the realignment of continuation lines."""
    shifts = []
    for pattern, new in RENAMES:
        pos = 0
        while True:
            m = pattern.search(line, pos)
            if not m:
                break
            repl = m.expand(new)
            line = line[:m.start()] + repl + line[m.end():]
            if len(repl) != len(m.group(0)):
                shifts.append((m.start(), len(repl) - len(m.group(0))))
            pos = m.start() + len(repl)
    return line, shifts


def generic(lines):
    """Renames identifiers, rewords comments, and keeps continuation lines
    aligned with the (longer) first line of their statement."""
    out = []
    pending = []          # shifts of the current continued statement
    for line in lines:
        if is_comment(line):
            new, _ = rename_line(line)
            for pattern, word in COMMENT_WORDS:
                new = pattern.sub(word, new)
            #'^TT' protects a plain transpose from the rule above.
            out.append(new.replace('^TT', '^T'))
            continue
        indent = len(line) - len(line.lstrip())
        grow = sum(d for col, d in pending if indent > col)
        new, shifts = rename_line(line)
        if grow > 0 and new.strip():
            new = ' ' * grow + new
        shifts = [(col + grow, d) for col, d in shifts]
        if line.rstrip().endswith('&'):
            pending = pending + shifts
        else:
            pending = []
        out.append(new)
    return wrap_calls(wrap_comments(out))


def wrap_comments(lines, width=79):
    """Moves the words that exceed the width (after the renaming made a
    comment line longer) to the following line of the same paragraph."""
    out = list(lines)
    k = 0
    while k < len(out):
        line = out[k]
        if not is_comment(line) or len(line) <= width:
            k += 1
            continue
        m = re.match(r'^(\s*!\s*)(.*)$', line)
        prefix, words = m.group(1), m.group(2).split(' ')
        keep = []
        while words and len(prefix + ' '.join(keep + [words[0]])) <= width:
            keep.append(words.pop(0))
        if not keep or not words:
            raise SystemExit('cannot wrap comment line: %r' % line)
        out[k] = prefix + ' '.join(keep)
        rest = ' '.join(words)
        nxt = out[k + 1] if k + 1 < len(out) else ''
        m2 = re.match(r'^(\s*!\s*)(\S.*)$', nxt) if is_comment(nxt) else None
        item = re.match(r'^(\s*!\s+)(\S+\s+- )', line)
        if m2 and len(m2.group(1)) >= len(prefix) and not (
                item and len(m2.group(1)) == len(prefix)):
            out[k + 1] = m2.group(1) + rest + ' ' + m2.group(2)
        else:
            hang = prefix
            if item:
                hang = item.group(1) + ' ' * len(item.group(2))
            out.insert(k + 1, hang + rest)
        k += 1
    return out


def split_arguments(text):
    """Splits an argument list at its top-level commas."""
    args, depth, cur = [], 0, ''
    for ch in text:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            args.append(cur.strip())
            cur = ''
        else:
            cur += ch
    args.append(cur.strip())
    return args


def wrap_calls(lines, width=80):
    """Refills the argument list of a call (or of a subroutine statement)
    that the renaming pushed beyond the width."""
    out = []
    k = 0
    head = re.compile(r'^(\s*(?:call |pure subroutine |subroutine )\w+\()(.*)$')
    while k < len(lines):
        line = lines[k]
        j = k
        while lines[j].rstrip().endswith('&') and not is_comment(lines[j]):
            j += 1
        stmt = lines[k:j + 1]
        m = head.match(line)
        if is_comment(line) or not m or max(len(s) for s in stmt) <= width:
            out += stmt
            k = j + 1
            continue
        text = m.group(2).rstrip().rstrip('&')
        for s in stmt[1:]:
            text += ' ' + s.strip().rstrip('&')
        text = text.rstrip()
        if not text.endswith(')'):
            raise SystemExit('cannot wrap statement: %r' % line)
        args = split_arguments(text[:-1])
        first, pad = m.group(1), ' ' * len(m.group(1))
        cur = first
        for n_arg, arg in enumerate(args):
            piece = arg + (', ' if n_arg < len(args) - 1 else ')')
            tail_len = 1 if n_arg < len(args) - 1 else 0
            if len(cur) + len(piece.rstrip()) + 1 + tail_len > width and cur.strip() != first.strip():
                out.append(cur.rstrip() + ' &')
                cur = pad
            cur += piece
        out.append(cur.rstrip())
        k = j + 1
    return out


def apply_edits(name, body, edits):
    """edits: list of (old, new, count); count None means at least one."""
    for old, new, count in edits:
        if isinstance(old, str):
            found = body.count(old)
            if found == 0 or (count is not None and found != count):
                raise SystemExit('%s: expected %s occurrence(s) of %r, found %d'
                                 % (name, count, old, found))
            body = body.replace(old, new)
        else:
            body, found = old.subn(new, body)
            if found == 0 or (count is not None and found != count):
                raise SystemExit('%s: expected %s match(es) of %r, found %d'
                                 % (name, count, old.pattern, found))
    return body


# --------------------------------------------------------------------------
# Procedure-specific substitutions (old text is the REAL source text)
# --------------------------------------------------------------------------

DIAG_CHECK_FRESH = '''\
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

'''

DIAG_CHECK_COLUMN = '''\
    !The diagonal elements of H and of S must be real to within roundoff.
    if (.not. (complex_diagonal_is_real(%(n)s, h_col, %(d)s) .and. &
               complex_diagonal_is_real(%(n)s, s_col, %(d)s))) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if
'''

DIAG_REAL = '''\
    !The diagonal element is real; a roundoff-sized imaginary part of the
    !supplied values is discarded.
    self%%wn(%(d)s, 3) = cmplx(real(h_col(%(d)s), wp) - &
                           self%%shift * real(s_col(%(d)s), wp), 0.0_wp, kind=wp)
'''

SQUARED = ('difference * difference',
           'real(difference, wp)**2 + aimag(difference)**2')

EDITS = {}

# Rewording of comments that the generic rules cannot derive. Applied after
# the procedure-specific code edits and before the generic pass; every entry
# must match exactly once in the named procedure.
COMMENT_EDITS = {}

EDITS['real_factorize_fresh'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n', 1),
    ('    integer :: blk, ii, matrix_n, np\n',
     '    integer :: blk, i, ii, j, matrix_n, np\n'
     '    real(wp) :: h_scale, s_scale\n', 1),
    ('    if (matrix_n > self%capacity) then\n'
     '      info = QR_ERR_CAPACITY_EXCEEDED\n'
     '      return\n'
     '    end if\n'
     '\n'
     '    np = self%nprocs\n',
     '    if (matrix_n > self%capacity) then\n'
     '      info = QR_ERR_CAPACITY_EXCEEDED\n'
     '      return\n'
     '    end if\n'
     '\n' + DIAG_CHECK_FRESH +
     '    np = self%nprocs\n', 1),
]

EDITS['real_factor_serial'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n', 1),
    ('    real(wp), allocatable :: a(:,:), tau(:), work(:), row(:)\n'
     '    real(wp) :: work_query(1)\n',
     '    complex(wp), allocatable :: a(:,:), tau(:), work(:), row(:)\n'
     '    complex(wp) :: shifted_value, work_query(1)\n', 1),
    ('          lwork = max(matrix_n, ceiling(work_query(1)))\n',
     '          lwork = max(matrix_n, ceiling(real(work_query(1), wp)))\n', 1),
    ('          lwork = max(lwork, ceiling(work_query(1)))\n',
     '          lwork = max(lwork, ceiling(real(work_query(1), wp)))\n', 1),
    ('    !Read each stored lower-triangular column contiguously and reflect every\n'
     '    !off-diagonal value into the opposite triangle. LAPACK requires the full\n'
     '    !shifted matrix, but callers need define only the canonical lower triangle.\n',
     '    !Read each stored lower-triangular column contiguously and reflect the\n'
     '    !complex conjugate of every off-diagonal value into the opposite triangle.\n'
     '    !LAPACK requires the full shifted matrix, but callers need define only the\n'
     '    !canonical lower triangle. The diagonal is formed from the real parts, so\n'
     '    !that the factorized matrix is exactly Hermitian.\n', 1),
    ('      do j = 1, matrix_n\n'
     '        do i = j, matrix_n\n'
     '          a(i, j) = h(i, j) - shift * s(i, j)\n'
     '          if (i /= j) a(j, i) = a(i, j)\n'
     '        end do\n'
     '      end do\n',
     '      do j = 1, matrix_n\n'
     '        a(j, j) = cmplx(real(h(j, j), wp) - shift * real(s(j, j), wp), &\n'
     '                        0.0_wp, kind=wp)\n'
     '        do i = j + 1, matrix_n\n'
     '          shifted_value = h(i, j) - shift * s(i, j)\n'
     '          a(i, j) = shifted_value\n'
     '          a(j, i) = conjg(shifted_value)\n'
     '        end do\n'
     '      end do\n', 1),
    ('          call MPI_Send(row, matrix_n, MPI_WP, dest, QR_TAG_R, self%comm, &\n',
     '          call MPI_Send(row, 2 * matrix_n, MPI_WP, dest, QR_TAG_R, &\n'
     '                        self%comm, &\n', 1),
    ('        call MPI_Recv(self%rt(1, l), matrix_n, MPI_WP, 0, QR_TAG_R, &\n',
     '        call MPI_Recv(self%rt(1, l), 2 * matrix_n, MPI_WP, 0, QR_TAG_R, &\n', 1),
    ('      call MPI_Bcast(self%rd, np * np * self%nblk, MPI_WP, 0, self%comm, &\n'
     '                     mpi_err)\n',
     '      call MPI_Bcast(self%rd, 2 * np * np * self%nblk, MPI_WP, 0, &\n'
     '                     self%comm, mpi_err)\n', 1),
    ('    !Expand the Householder vectors and tau into the explicit orthogonal\n'
     '    !matrix Q and distribute its rows.\n',
     '    !Expand the Householder vectors and tau into the explicit unitary matrix\n'
     '    !Q and distribute its rows. A rank stores the complex conjugate of each of\n'
     '    !its rows, that is, a column of Q^H.\n', 1),
    ('            self%qt(j, blk) = a(i, j)\n',
     '            self%qt(j, blk) = conjg(a(i, j))\n', 1),
    ('          do j = 1, matrix_n\n'
     '            row(j) = a(i, j)\n'
     '          end do\n'
     '          call MPI_Send(row, matrix_n, MPI_WP, dest, QR_TAG_Q, self%comm, &\n',
     '          do j = 1, matrix_n\n'
     '            row(j) = conjg(a(i, j))\n'
     '          end do\n'
     '          call MPI_Send(row, 2 * matrix_n, MPI_WP, dest, QR_TAG_Q, &\n'
     '                        self%comm, &\n', 1),
    ('        call MPI_Recv(self%qt(1, l), matrix_n, MPI_WP, 0, QR_TAG_Q, &\n',
     '        call MPI_Recv(self%qt(1, l), 2 * matrix_n, MPI_WP, 0, QR_TAG_Q, &\n', 1),
]

# ---- Distributed (stage B) factorization ----------------------------------

EDITS['real_factor_calibrate'] = [
    ('    real(wp), allocatable :: a(:,:)\n',
     '    complex(wp), allocatable :: a(:,:)\n', 1),
    ('          a(i, j) = 1.0_wp / real(1 + abs(i - j), wp)\n',
     '          a(i, j) = cmplx(1.0_wp / real(1 + abs(i - j), wp), &\n'
     '                          0.25_wp * real(i - j, wp) / &\n'
     '                          real((1 + abs(i - j))**2, wp), kind=wp)\n', 1),
]

COMMENT_EDITS['real_factor_calibrate'] = [
    ('      !Synthetic diagonally dominant matrix; every leading block is well\n'
     '      !conditioned. It serves as H, and the shift is zero.\n',
     '      !Synthetic Hermitian matrix with a dominant diagonal; every leading\n'
     '      !block is well conditioned. It serves as H, and the shift is zero.\n', 1),
]

EDITS['real_factor_distributed'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n', 1),
    ('    real(wp), allocatable :: pan(:,:), red(:), tau(:), tsave(:,:), wbuf(:)\n'
     '    real(wp), allocatable :: work(:), ymat(:,:)\n',
     '    complex(wp), allocatable :: pan(:,:), red(:), tau(:), tsave(:,:), wbuf(:)\n'
     '    complex(wp), allocatable :: work(:), ymat(:,:)\n', 1),
    ('      self%rt(i, l) = h(i, i) - shift * s(i, i)\n',
     '      self%rt(i, l) = cmplx(real(h(i, i), wp) - shift * real(s(i, i), wp), &\n'
     '                            0.0_wp, kind=wp)\n', 1),
    ('        self%rt(j, l) = h(j, i) - shift * s(j, i)\n',
     '        self%rt(j, l) = conjg(h(j, i) - shift * s(j, i))\n', 1),
    ('        call MPI_Allreduce(MPI_IN_PLACE, wbuf, (m + pb) * pb, MPI_WP, &\n'
     '                           MPI_SUM, self%comm, mpi_err)\n',
     '        call MPI_Allreduce(MPI_IN_PLACE, wbuf, 2 * (m + pb) * pb, MPI_WP, &\n'
     '                           MPI_SUM, self%comm, mpi_err)\n', 1),
    ('        call MPI_Allreduce(MPI_IN_PLACE, wbuf, m * pb, MPI_WP, MPI_SUM, &\n'
     '                           self%comm, mpi_err)\n',
     '        call MPI_Allreduce(MPI_IN_PLACE, wbuf, 2 * m * pb, MPI_WP, MPI_SUM, &\n'
     '                           self%comm, mpi_err)\n', 1),
    ('      call MPI_Allreduce(MPI_IN_PLACE, self%rd, np * np * nblk, MPI_WP, &\n'
     '                         MPI_SUM, self%comm, mpi_err)\n',
     '      call MPI_Allreduce(MPI_IN_PLACE, self%rd, 2 * np * np * nblk, MPI_WP, &\n'
     '                         MPI_SUM, self%comm, mpi_err)\n', 1),
    ('    !Clear the reflector vectors stored below the diagonal of R, so that\n',
     '    !The reflectors were applied to Q itself. A rank stores the complex\n'
     '    !conjugate of each of its rows, that is, a column of Q^H.\n'
     '    do l = 1, nloc\n'
     '      do j = 1, n\n'
     '        self%qt(j, l) = conjg(self%qt(j, l))\n'
     '      end do\n'
     '    end do\n'
     '\n'
     '    !Clear the reflector vectors stored below the diagonal of R, so that\n', 1),
]

COMMENT_EDITS['real_factor_distributed'] = [
    ('  !lower triangles of H and S; there is no n by n temporary on any rank.\n',
     '  !lower triangles of H and S (the part of a row to the right of the diagonal\n'
     '  !is the complex conjugate of a stored column, and the diagonal is formed\n'
     '  !from the real parts); there is no n by n temporary on any rank.\n', 1),
    ('  !with the same two kernels and one reduction per panel. Finally the part of\n'
     '  !rt below the diagonal is cleared and the diagonal blocks of R are\n'
     '  !collected in rd by one reduction.\n',
     '  !with the same two kernels and one reduction per panel. Finally qt is\n'
     '  !replaced by its complex conjugate (a rank stores columns of Q^H), the part\n'
     '  !of rt below the diagonal is cleared, and the diagonal blocks of R are\n'
     '  !collected in rd by one reduction. The diagonal of R is real, as in ZGEQRF.\n', 1),
    ('    !right is column i of the lower triangles.\n',
     '    !right is the complex conjugate of column i of the lower triangles.\n', 1),
]

EDITS['real_panel_factor'] = [
    ('    real(wp), intent(inout) :: pan(ldp, *)\n'
     '    real(wp), intent(out) :: tau(pb)\n'
     '    real(wp), intent(inout) :: red(2 * pb)\n'
     '    real(wp) :: alpha, factor2, vl\n',
     '    complex(wp), intent(inout) :: pan(ldp, *)\n'
     '    complex(wp), intent(out) :: tau(pb)\n'
     '    complex(wp), intent(inout) :: red(2 * pb)\n'
     '    complex(wp) :: alpha, factor2, vl\n', 1),
    ('xnorm2 = red(1)\n', 'xnorm2 = real(red(1), wp)\n', 2),
    ('        red(c) = tau(jj) * (red(nc + c) + (factor2 * red(c)) * factor3)\n',
     '        red(c) = conjg(tau(jj)) * (red(nc + c) + &\n'
     '                                   (conjg(factor2) * red(c)) * factor3)\n', 1),
    ('        call MPI_Allreduce(MPI_IN_PLACE, red, 2 * nc, MPI_WP, MPI_SUM, comm, &\n'
     '                           mpi_err)\n',
     '        call MPI_Allreduce(MPI_IN_PLACE, red, 4 * nc, MPI_WP, MPI_SUM, comm, &\n'
     '                           mpi_err)\n', 1),
]

COMMENT_EDITS['real_panel_factor'] = [
    ('  !     z(c) = sum_i a(i,j)*a(i,c) for the panel columns c >= j (z(j) is the\n',
     '  !     z(c) = sum_i conjg(a(i,j))*a(i,c) for the panel columns c >= j (z(j) is the\n', 1),
    ('  !     norm (real_reflector): the new diagonal element beta, tau, and the\n',
     '  !     norm (real_reflector): the new real diagonal element beta, tau, and the\n', 1),
    ('      !red(c) becomes tau*(v^T*a(:,c)), the row that the reflector subtracts\n',
     '      !red(c) becomes conjg(tau)*(v^T*a(:,c)), the row that the adjoint of the\n'
     '      !reflector subtracts\n', 1),
]

EDITS['real_panel_column'] = [
    ('    real(wp), intent(in) :: x(ldx, *)\n', '    complex(wp), intent(in) :: x(ldx, *)\n', 1),
    ('    real(wp), intent(out) :: z(nc)\n', '    complex(wp), intent(out) :: z(nc)\n', 1),
    ('    real(wp) :: xs\n', '    complex(wp) :: xs\n', 1),
    ('      xs = x(1, i) * sc\n', '      xs = conjg(x(1, i)) * sc\n', 1),
    ('      z(1) = z(1) + xs * xs\n',
     '      z(1) = z(1) + (real(xs, wp)**2 + aimag(xs)**2)\n', 1),
]

COMMENT_EDITS['real_panel_column'] = [
    ('  !        z(1) = sum_i (sc*x(1,i))**2 ,\n'
     '  !        z(c) = sum_i sc*x(1,i)*x(c,i) ,   c = 2..nc ,\n',
     '  !        z(1) = sum_i abs(sc*x(1,i))**2 ,\n'
     '  !        z(c) = sum_i sc*conjg(x(1,i))*x(c,i) ,   c = 2..nc ,\n', 1),
]

EDITS['real_scale_code'] = [
    ('    real(wp), intent(in) :: x(*)\n', '    complex(wp), intent(in) :: x(*)\n', 1),
    ('      magnitude = abs(x(k))\n',
     '      magnitude = max(abs(real(x(k), wp)), abs(aimag(x(k))))\n', 1),
]

COMMENT_EDITS['real_scale_code'] = [
    ('  !        otherwise the binary exponent of the largest absolute value.\n',
     '  !        otherwise the binary exponent of the largest absolute value of a\n'
     '  !        real or imaginary part.\n', 1),
]

EDITS['real_reflector'] = [
    ('    real(wp), intent(in) :: alpha\n', '    complex(wp), intent(in) :: alpha\n', 1),
    ('    real(wp), intent(out) :: tau, factor2\n',
     '    complex(wp), intent(out) :: tau, factor2\n', 1),
    ('    real(wp) :: alpha_scaled\n', '    complex(wp) :: alpha_scaled\n', 1),
    ('    if (xnorm == 0.0_wp) then\n'
     '      beta = alpha\n',
     '    if (xnorm == 0.0_wp .and. aimag(alpha) == 0.0_wp) then\n'
     '      beta = real(alpha, wp)\n', 1),
    ('      beta = -sign(norm, alpha)\n',
     '      beta = -sign(norm, real(alpha, wp))\n', 1),
    ('      beta_scaled = -sign(norm, alpha_scaled)\n',
     '      beta_scaled = -sign(norm, real(alpha_scaled, wp))\n', 1),
    # No complex modulus and no division by a complex number is left to the
    # compiler: not every compiler evaluates them without squaring the parts
    # (which overflows or underflows for very large or very small arguments).
    ('    norm = hypot(abs(alpha), xnorm)\n',
     '    norm = hypot(hypot(real(alpha, wp), aimag(alpha)), xnorm)\n', 1),
    ('      norm = hypot(abs(alpha_scaled), xnorm_scaled)\n',
     '      norm = hypot(hypot(real(alpha_scaled, wp), aimag(alpha_scaled)), &\n'
     '                   xnorm_scaled)\n', 1),
    ('      tau = (beta - alpha) / beta\n',
     '      tau = cmplx((beta - real(alpha, wp)) / beta, -aimag(alpha) / beta, &\n'
     '                  kind=wp)\n', 1),
    ('      tau = (beta_scaled - alpha_scaled) / beta_scaled\n',
     '      tau = cmplx((beta_scaled - real(alpha_scaled, wp)) / beta_scaled, &\n'
     '                  -aimag(alpha_scaled) / beta_scaled, kind=wp)\n', 1),
    ('      factor2 = 1.0_wp / (alpha - beta)\n',
     '      factor2 = complex_reciprocal(alpha - beta)\n', 1),
    ('      factor2 = 1.0_wp / (alpha_scaled - beta_scaled)\n',
     '      factor2 = complex_reciprocal(alpha_scaled - beta_scaled)\n', 1),
]

COMMENT_EDITS['real_reflector'] = [
    ('  !with the conventions of the LAPACK routine DLARFG: beta has the sign\n'
     '  !opposite to alpha, tau = (beta-alpha)/beta, v2 = x/(alpha-beta), and\n'
     '  !H = I (tau = 0) when x is zero.\n',
     '  !with the conventions of the LAPACK routine DLARFG: beta is real and has\n'
     '  !the sign opposite to the real part of alpha, tau = (beta-alpha)/beta is\n'
     '  !complex, v2 = x/(alpha-beta), and H = I (tau = 0) when x is zero and alpha\n'
     '  !is real. The diagonal of R is therefore real.\n', 1),
    ('  !        v^T*a(:,c) = a(j,c) + (factor2*z(c))*factor3 ,\n'
     '  !\n'
     '  !where z(c) = sum_i sc*x(i)*a(i,c). When alpha and x are both so small\n',
     '  !        v^T*a(:,c) = a(j,c) + (conjg(factor2)*z(c))*factor3 ,\n'
     '  !\n'
     '  !where z(c) = sum_i sc*conjg(x(i))*a(i,c). When alpha and x are both so small\n', 1),
    ('  !    beta   - The new diagonal element.\n',
     '  !    beta   - The new diagonal element (real).\n', 1),
    ('  !that beta is below QR_SAFE_MINIMUM, the reflector is generated from the\n'
     '  !scaled quantities, as DLARFG does.\n',
     '  !that beta is below QR_SAFE_MINIMUM, the reflector is generated from the\n'
     '  !scaled quantities, as DLARFG does. Moduli are formed with HYPOT and the\n'
     '  !reciprocal of alpha-beta with real_reciprocal, so that no real or\n'
     '  !imaginary part is squared.\n', 1),
]

EDITS['real_block_product'] = [
    ('    real(wp), intent(in) :: b(ldb, *), v(ldv, *)\n'
     '    real(wp), intent(inout) :: w(ldw, *)\n'
     '    real(wp), intent(inout) :: work(*)\n',
     '    complex(wp), intent(in) :: b(ldb, *), v(ldv, *)\n'
     '    complex(wp), intent(inout) :: w(ldw, *)\n'
     '    complex(wp), intent(inout) :: work(*)\n', 1),
    ('    real(wp) :: c1, c2, c3, c4\n', '    complex(wp) :: c1, c2, c3, c4\n', 1),
    ('            work(l + (jj - 1) * nl) = v(jj, l0 + l - 1)\n',
     '            work(l + (jj - 1) * nl) = conjg(v(jj, l0 + l - 1))\n', 1),
    (re.compile(r'(c\d) = v\(jj, (l[^)]*)\)\n'), r'\1 = conjg(v(jj, \2))\n', 5),
]

COMMENT_EDITS['real_block_product'] = [
    ('  !        w(i,jj) = w(i,jj) + sum_l b(i,l)*v(jj,l) ,\n',
     '  !        w(i,jj) = w(i,jj) + sum_l b(i,l)*conjg(v(jj,l)) ,\n', 1),
]

EDITS['real_block_update'] = [
    ('    real(wp), intent(in) :: y(ldy, *), w(ldw, *)\n'
     '    real(wp), intent(inout) :: b(ldb, *)\n'
     '    real(wp), intent(inout) :: work(*)\n',
     '    complex(wp), intent(in) :: y(ldy, *), w(ldw, *)\n'
     '    complex(wp), intent(inout) :: b(ldb, *)\n'
     '    complex(wp), intent(inout) :: work(*)\n', 1),
    ('    real(wp) :: c1, c2, c3, c4, c5, c6, c7, c8\n',
     '    complex(wp) :: c1, c2, c3, c4, c5, c6, c7, c8\n', 1),
    ('    logical, parameter :: EIGHT_COLUMNS = .false.\n',
     '    logical, parameter :: EIGHT_COLUMNS = .true.\n', 1),
]

COMMENT_EDITS['real_block_update'] = [
    ('  !in real_block_product, and four columns of w at a time (eight when\n'
     '  !EIGHT_COLUMNS is set; four were measured to be slightly faster for real\n'
     '  !data in double precision).\n',
     '  !in real_block_product, and eight columns of w at a time (EIGHT_COLUMNS;\n'
     '  !for complex data in double precision eight were measured to be clearly\n'
     '  !faster than four).\n', 1),
]

EDITS['real_dot_multiply'] = [
    ('    real(wp), intent(in) :: a(lda, *), b(ldb, *)\n',
     '    complex(wp), intent(in) :: a(lda, *), b(ldb, *)\n', 1),
    ('    real(wp), intent(inout) :: c(ldc, *)\n',
     '    complex(wp), intent(inout) :: c(ldc, *)\n', 1),
    ('    real(wp) :: t1, t2, t3, t4, x\n', '    complex(wp) :: t1, t2, t3, t4, x\n', 1),
]

COMMENT_EDITS['real_dot_multiply'] = [
    ('  !Subroutine real_dot_multiply adds a multiple of the product of the\n'
     '  !transpose of one matrix with another matrix to a third,\n',
     '  !Subroutine real_dot_multiply adds a multiple of the product of the plain\n'
     '  !(unconjugated) transpose of one matrix with another matrix to a third,\n', 1),
]

EDITS['real_block_gram'] = [
    ('    real(wp), intent(in) :: v(ldv, *)\n'
     '    real(wp), intent(inout) :: g(ldg, *)\n',
     '    complex(wp), intent(in) :: v(ldv, *)\n'
     '    complex(wp), intent(inout) :: g(ldg, *)\n', 1),
    ('    real(wp) :: vb\n', '    complex(wp) :: vb\n', 1),
    ('          g(a, b) = g(a, b) + v(a, l) * vb\n',
     '          g(a, b) = g(a, b) + conjg(v(a, l)) * vb\n', 1),
]

COMMENT_EDITS['real_block_gram'] = [
    ('  !        g(a,b) = g(a,b) + sum_l v(a,l)*v(b,l) ,   a < b ,\n',
     '  !        g(a,b) = g(a,b) + sum_l conjg(v(a,l))*v(b,l) ,   a < b ,\n', 1),
]

EDITS['real_block_tfactor'] = [
    ('    real(wp), intent(in) :: tau(k), g(ldg, *)\n'
     '    real(wp), intent(out) :: t(ldt, *)\n',
     '    complex(wp), intent(in) :: tau(k), g(ldg, *)\n'
     '    complex(wp), intent(out) :: t(ldt, *)\n', 1),
    ('    real(wp) :: tsum\n', '    complex(wp) :: tsum\n', 1),
]

EDITS['real_block_vt'] = [
    ('    real(wp), intent(in) :: v(ldv, *), t(ldt, *)\n',
     '    complex(wp), intent(in) :: v(ldv, *), t(ldt, *)\n', 1),
    ('    real(wp), intent(out) :: y(ldy, *)\n', '    complex(wp), intent(out) :: y(ldy, *)\n', 1),
    ('    real(wp) :: tsum\n', '    complex(wp) :: tsum\n', 1),
    ('            tsum = tsum + t(jj, a) * v(a, l)\n',
     '            tsum = tsum + conjg(t(jj, a)) * v(a, l)\n', 1),
]

COMMENT_EDITS['real_block_vt'] = [
    ('  !        Y = V*T^T   (adjoint true:  y(jj,l) = sum_{a>=jj} t(jj,a)*v(a,l))\n',
     '  !        Y = V*T^T   (adjoint true:  y(jj,l) = sum_{a>=jj} conjg(t(jj,a))*v(a,l))\n', 1),
]

EDITS['real_replace_symmetric'] = [
    ('    real(wp), intent(in) :: h_col(:), s_col(:)\n',
     '    complex(wp), intent(in) :: h_col(:), s_col(:)\n', 1),
    ('    if (size(h_col) /= self%n .or. size(s_col) /= self%n) then\n'
     '      info = QR_ERR_DIMENSION_MISMATCH\n'
     '      return\n'
     '    end if\n',
     '    if (size(h_col) /= self%n .or. size(s_col) /= self%n) then\n'
     '      info = QR_ERR_DIMENSION_MISMATCH\n'
     '      return\n'
     '    end if\n' + DIAG_CHECK_COLUMN % {'n': 'self%n', 'd': 'idx'}, 1),
    ('      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)\n'
     '    end do\n',
     '      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)\n'
     '    end do\n' + DIAG_REAL % {'d': 'idx'}, 1),
]

EDITS['real_append_symmetric'] = [
    ('    real(wp), intent(in) :: h_col(:), s_col(:)\n',
     '    complex(wp), intent(in) :: h_col(:), s_col(:)\n', 1),
    ('    if (size(h_col) /= new_n .or. size(s_col) /= new_n) then\n'
     '      info = QR_ERR_DIMENSION_MISMATCH\n'
     '      return\n'
     '    end if\n',
     '    if (size(h_col) /= new_n .or. size(s_col) /= new_n) then\n'
     '      info = QR_ERR_DIMENSION_MISMATCH\n'
     '      return\n'
     '    end if\n' + DIAG_CHECK_COLUMN % {'n': 'new_n', 'd': 'new_n'}, 1),
    ('      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)\n'
     '    end do\n',
     '      self%wn(i, 3) = h_col(i) - self%shift * s_col(i)\n'
     '    end do\n' + DIAG_REAL % {'d': 'new_n'}, 1),
]

EDITS['real_solve'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    real(wp), intent(in) :: v_initial(:)\n'
     '    real(wp), intent(out) :: x(:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    complex(wp), intent(in) :: v_initial(:)\n'
     '    complex(wp), intent(out) :: x(:)\n', 1),
]

EDITS['real_inverse_iteration'] = [
    ('    real(wp), intent(in) :: h(ldh, *), s(lds, *)\n'
     '    real(wp), intent(inout) :: v(n), w(n), y(n)\n'
     '    real(wp), intent(out) :: x(n)\n',
     '    complex(wp), intent(in) :: h(ldh, *), s(lds, *)\n'
     '    complex(wp), intent(inout) :: v(n), w(n), y(n)\n'
     '    complex(wp), intent(out) :: x(n)\n', 1),
    ('    real(wp) :: coefficient, current_norm_squared, eigenvector_norm_squared\n',
     '    complex(wp) :: coefficient\n'
     '    real(wp) :: current_norm_squared, eigenvector_norm_squared, quotient\n', 1),
    ('        overlap_norm_squared = real_dot(n, v, w)\n',
     '        overlap_norm_squared = real(complex_dotc(n, v, w), wp)\n', 1),
    ('          coefficient = real_dot(n, v, y) / overlap_norm_squared\n'
     '          do i = 1, n\n'
     '            y(i) = y(i) - coefficient * w(i)\n',
     '          quotient = real(complex_dotc(n, v, y), wp) / overlap_norm_squared\n'
     '          do i = 1, n\n'
     '            y(i) = y(i) - quotient * w(i)\n', 1),
    ('      !Scale the solution so its largest absolute component is one. Scaling\n'
     '      !does not change the eigendirection and limits growth in repeated solves.\n'
     '      max_component = 0.0_wp\n'
     '      do i = 1, n\n'
     '        max_component = max(max_component, abs(x(i)))\n'
     '      end do\n',
     '      !Scale the solution with its largest real or imaginary component. This\n'
     '      !bounds both parts of every element without changing the complex\n'
     '      !eigendirection and limits growth in repeated solves.\n'
     '      max_component = complex_max_abs_real_or_imag(n, x)\n', 1),
    ('      !Compute the least-direction-change coefficient and the relative norm of\n'
     '      !the component of x not parallel to the previous iterate v.\n',
     '      !Compute the least-direction-change coefficient (v^H*x)/(v^H*v) and the\n'
     '      !relative norm of the component of x not parallel to the previous\n'
     '      !iterate v. Conjugating the previous iterate makes the comparison\n'
     '      !invariant under the arbitrary complex phase of an eigenvector.\n', 1),
    ('      coefficient = real_dot(n, x, v) / current_norm_squared\n',
     '      coefficient = complex_dotc(n, v, x) / current_norm_squared\n', 1),
    ('        max_component = 0.0_wp\n'
     '        do i = 1, n\n'
     '          y(i) = x(i) - y(i)\n'
     '          max_component = max(max_component, abs(y(i)))\n'
     '        end do\n',
     '        do i = 1, n\n'
     '          y(i) = x(i) - y(i)\n'
     '        end do\n'
     '        max_component = complex_max_abs_real_or_imag(n, y)\n', 1),
    ('      coefficient = sqrt(', '      quotient = sqrt(', 2),
    ('        x(i) = x(i) / coefficient\n', '        x(i) = x(i) / quotient\n', 2),
    ('    !for other norm_mode values because each iterate already has unit\n'
     '    !largest-component magnitude.\n',
     '    !for other norm_mode values because the largest real or imaginary\n'
     '    !component of each iterate already has unit magnitude.\n', 1),
]

EDITS['real_physical_quotient'] = [
    ('    real(wp), intent(in) :: h(ldh, *), s(lds, *), x(n)\n'
     '    real(wp), intent(out) :: hx(n), sx(n)\n',
     '    complex(wp), intent(in) :: h(ldh, *), s(lds, *), x(n)\n'
     '    complex(wp), intent(out) :: hx(n), sx(n)\n', 1),
    ('    real(wp) :: difference, residual_norm_squared\n',
     '    complex(wp) :: difference\n'
     '    real(wp) :: residual_norm_squared\n', 1),
    ('    overlap = real_dot(n, x, sx)\n',
     '    overlap = real(complex_dotc(n, x, sx), wp)\n', 1),
    ('    lambda = real_dot(n, x, hx) / overlap\n',
     '    lambda = real(complex_dotc(n, x, hx), wp) / overlap\n', 1),
    ('      residual_norm_squared = residual_norm_squared + difference * difference\n',
     '      residual_norm_squared = residual_norm_squared + &\n'
     '                              real(difference, wp)**2 + aimag(difference)**2\n', 1),
]

EDITS['real_solve_shifted'] = [
    ('    real(wp), intent(in) :: b(:)\n'
     '    real(wp), intent(out) :: x(:)\n',
     '    complex(wp), intent(in) :: b(:)\n'
     '    complex(wp), intent(out) :: x(:)\n', 1),
]

EDITS['real_symmetric_product'] = [
    ('    real(wp), intent(in), contiguous :: a(:,:)\n'
     '    real(wp), intent(in) :: x(:)\n'
     '    real(wp), intent(out) :: y(:)\n',
     '    complex(wp), intent(in), contiguous :: a(:,:)\n'
     '    complex(wp), intent(in) :: x(:)\n'
     '    complex(wp), intent(out) :: y(:)\n', 1),
]

EDITS['real_shifted_quotient'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    real(wp), intent(in) :: x(:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    complex(wp), intent(in) :: x(:)\n', 1),
]

# The quadratic forms are real: a column contributes the real diagonal term
# and twice the real part of x(j) times the sum over the stored elements below
# the diagonal, each multiplied by the conjugate of its row component.
EDITS['real_shifted_column'] = [
    ('    real(wp), intent(in) :: h_col(0:m), s_col(0:m), shift, x(0:m)\n',
     '    complex(wp), intent(in) :: h_col(0:m), s_col(0:m), x(0:m)\n'
     '    real(wp), intent(in) :: shift\n', 1),
    ('    real(wp) :: form_below, overlap_below\n',
     '    complex(wp) :: form_below, overlap_below\n'
     '    real(wp) :: x_squared\n', 1),
    ('      form_below = form_below + (h_col(i) - shift * s_col(i)) * x(i)\n'
     '      overlap_below = overlap_below + s_col(i) * x(i)\n',
     '      form_below = form_below + (h_col(i) - shift * s_col(i)) * conjg(x(i))\n'
     '      overlap_below = overlap_below + s_col(i) * conjg(x(i))\n', 1),
    ('    form = x(0) * ((h_col(0) - shift * s_col(0)) * x(0) + 2.0_wp * form_below)\n'
     '    overlap = x(0) * (s_col(0) * x(0) + 2.0_wp * overlap_below)\n',
     '    !The diagonal elements are real; the imaginary parts of the supplied\n'
     '    !values are ignored, as in the products with H and S.\n'
     '    x_squared = real(x(0), wp)**2 + aimag(x(0))**2\n'
     '    form = (real(h_col(0), wp) - shift * real(s_col(0), wp)) * x_squared + &\n'
     '           2.0_wp * real(x(0) * form_below, wp)\n'
     '    overlap = real(s_col(0), wp) * x_squared + &\n'
     '              2.0_wp * real(x(0) * overlap_below, wp)\n', 1),
]

EDITS['real_factorization_residual'] = [
    ('    real(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    real(wp), intent(in) :: v(:)\n',
     '    complex(wp), intent(in), contiguous :: h(:,:), s(:,:)\n'
     '    complex(wp), intent(in) :: v(:)\n', 1),
    ('    real(wp) :: action_norm, difference, difference_norm_squared\n',
     '    complex(wp) :: difference\n'
     '    real(wp) :: action_norm, difference_norm_squared\n', 1),
    ('                                  difference * difference\n',
     '                                  real(difference, wp)**2 + &\n'
     '                                  aimag(difference)**2\n', 1),
    ('                                difference * difference\n',
     '                                real(difference, wp)**2 + &\n'
     '                                aimag(difference)**2\n', 1),
]

EDITS['real_join_base_slot'] = [
    ('    real(wp) :: current, entry\n', '    complex(wp) :: current, entry\n', 1),
    ('      self%cbb(slot, slot2) = current\n'
     '      self%cbb(slot2, slot) = current\n',
     '      self%cbb(slot, slot2) = current\n'
     '      self%cbb(slot2, slot) = conjg(current)\n', 1),
    ('      !Current value of the matrix element (idx,idx2): the stored column of a\n'
     '      !member of D0 is a difference from its base column.\n',
     '      !Current value of the matrix element (idx,idx2): the stored column of a\n'
     '      !member of D0 is a difference from its base column. The element\n'
     '      !(idx2,idx) of the Hermitian border-border block is its conjugate.\n', 1),
]

EDITS['real_install_column'] = [
    ('    real(wp), intent(in) :: a(self%n)\n', '    complex(wp), intent(in) :: a(self%n)\n', 1),
    ('    real(wp) :: current, value\n', '    complex(wp) :: current, value\n', 1),
    ('          current = self%ab(idx, q)\n'
     '          if (current == 0.0_wp) current = 1.0_wp\n'
     '        else\n'
     '          current = a(idx)\n'
     '        end if\n',
     '          current = cmplx(real(self%ab(idx, q), wp), 0.0_wp, kind=wp)\n'
     '          if (current == 0.0_wp) current = 1.0_wp\n'
     '        else\n'
     '          current = cmplx(real(a(idx), wp), 0.0_wp, kind=wp)\n'
     '        end if\n', 1),
    ('      self%cbb(slot, slot2) = current\n'
     '      self%cbb(slot2, slot) = current\n',
     '      !a(idx2) is the element (idx2,idx) of the conceptual column idx; the\n'
     '      !element (idx,idx2) of the Hermitian matrix is its conjugate.\n'
     '      self%cbb(slot2, slot) = current\n'
     '      self%cbb(slot, slot2) = conjg(current)\n', 1),
]

EDITS['real_fill_schur'] = [
    ('    real(wp) :: base12, base21, own\n',
     '    complex(wp) :: base12, base21\n'
     '    real(wp) :: own\n', 1),
]

EDITS['real_rebuild_small'] = [
    ('        call MPI_Allreduce(MPI_IN_PLACE, self%wred, count, MPI_WP, MPI_SUM, &\n'
     '                           self%comm, mpi_err)\n',
     '        call MPI_Allreduce(MPI_IN_PLACE, self%wred, 2 * count, MPI_WP, &\n'
     '                           MPI_SUM, self%comm, mpi_err)\n', 1),
]

EDITS['real_refresh_small'] = [
    ('    real(wp) :: tsum\n', '    complex(wp) :: tsum\n', 1),
]

EDITS['real_bordered_solve'] = [
    ('    real(wp), intent(in) :: r(self%n)\n'
     '    real(wp), intent(out) :: x(self%n)\n',
     '    complex(wp), intent(in) :: r(self%n)\n'
     '    complex(wp), intent(out) :: x(self%n)\n', 1),
    ('    real(wp) :: tsum\n', '    complex(wp) :: tsum\n', 1),
]

EDITS['real_base_solve'] = [
    ('    real(wp), intent(in) :: b(self%n0)\n'
     '    real(wp), intent(out) :: x(self%n0)\n',
     '    complex(wp), intent(in) :: b(self%n0)\n'
     '    complex(wp), intent(out) :: x(self%n0)\n', 1),
    ('    real(wp) :: tsum\n', '    complex(wp) :: tsum\n', 1),
    ('    !y = Q^T*b, accumulated in x.\n',
     '    !y = Q^H*b, accumulated in x. The local columns of qt are columns of Q^H,\n'
     '    !so no conjugation is needed here.\n', 1),
    ('      call MPI_Allreduce(MPI_IN_PLACE, x, n0, MPI_WP, MPI_SUM, self%comm, &\n'
     '                         mpi_err)\n',
     '      call MPI_Allreduce(MPI_IN_PLACE, x, 2 * n0, MPI_WP, MPI_SUM, &\n'
     '                         self%comm, mpi_err)\n', 1),
    ('          tsum = tsum - real_dot(n0 - i, self%rt(i + 1:n0, i), x(i + 1:n0))\n',
     '          tsum = tsum - complex_dotu(n0 - i, self%rt(i + 1:n0, i), &\n'
     '                                     x(i + 1:n0))\n', 1),
    ('            self%wblk(self%rank + 1) = real_dot(n0 - ie, &\n'
     '                                                self%rt(ie + 1:n0, blk), &\n'
     '                                                x(ie + 1:n0))\n',
     '            self%wblk(self%rank + 1) = complex_dotu(n0 - ie, &\n'
     '                                                    self%rt(ie + 1:n0, blk), &\n'
     '                                                    x(ie + 1:n0))\n', 1),
    ('          call MPI_Allreduce(MPI_IN_PLACE, self%wblk, nb, MPI_WP, MPI_SUM, &\n'
     '                             self%comm, mpi_err)\n',
     '          call MPI_Allreduce(MPI_IN_PLACE, self%wblk, 2 * nb, MPI_WP, &\n'
     '                             MPI_SUM, self%comm, mpi_err)\n', 1),
]

EDITS['real_base_column'] = [
    ('    real(wp), intent(out) :: column(self%n0), work(self%n0)\n',
     '    complex(wp), intent(out) :: column(self%n0), work(self%n0)\n', 1),
    ('      call MPI_Allreduce(MPI_IN_PLACE, work, n0, MPI_WP, MPI_SUM, self%comm, &\n'
     '                         mpi_err)\n',
     '      call MPI_Allreduce(MPI_IN_PLACE, work, 2 * n0, MPI_WP, MPI_SUM, &\n'
     '                         self%comm, mpi_err)\n', 1),
    ('      call MPI_Allreduce(MPI_IN_PLACE, column, n0, MPI_WP, MPI_SUM, &\n',
     '      call MPI_Allreduce(MPI_IN_PLACE, column, 2 * n0, MPI_WP, MPI_SUM, &\n', 1),
    ('      column(i) = real_dot(idx, self%qt(1:idx, l), work(1:idx))\n',
     '      !Row i of Q is the conjugate of the stored column of Q^H.\n'
     '      column(i) = real_dot(idx, self%qt(1:idx, l), work(1:idx))\n', 1),
]

EDITS['real_sym_product'] = [
    ('    real(wp), intent(in) :: a(lda, *), x(n)\n'
     '    real(wp), intent(out) :: y(n)\n',
     '    complex(wp), intent(in) :: a(lda, *), x(n)\n'
     '    complex(wp), intent(out) :: y(n)\n', 1),
    ('    real(wp) :: t1, t2, t3, t4, u1, u2, u3\n',
     '    complex(wp) :: t1, t2, t3, t4, u1, u2, u3\n', 1),
    (re.compile(r'a\((j\d?), \1\) \* x\(\1\)'), r'real(a(\1, \1), wp) * x(\1)', 8),
    ('      call MPI_Allreduce(MPI_IN_PLACE, y, n, MPI_WP, MPI_SUM, comm, mpi_err)\n',
     '      call MPI_Allreduce(MPI_IN_PLACE, y, 2 * n, MPI_WP, MPI_SUM, comm, &\n'
     '                         mpi_err)\n', 1),
]

EDITS['real_symv_column'] = [
    ('    real(wp), intent(in) :: a(n), xj, x(n)\n'
     '    real(wp), intent(inout) :: y(n)\n'
     '    real(wp), intent(out) :: tsum\n',
     '    complex(wp), intent(in) :: a(n), xj, x(n)\n'
     '    complex(wp), intent(inout) :: y(n)\n'
     '    complex(wp), intent(out) :: tsum\n', 1),
    ('    real(wp) :: t1, t2, t3, t4\n', '    complex(wp) :: t1, t2, t3, t4\n', 1),
    (re.compile(r'(t\d) = \1 \+ a\(([^)]*)\) \* x\(([^)]*)\)'),
     r'\1 = \1 + conjg(a(\2)) * x(\3)', 5),
]

EDITS['real_symv_columns4'] = [
    ('    real(wp), intent(in) :: a1(n), a2(n), a3(n), a4(n), x1, x2, x3, x4, x(n)\n'
     '    real(wp), intent(inout) :: y(n)\n'
     '    real(wp), intent(out) :: t1, t2, t3, t4\n',
     '    complex(wp), intent(in) :: a1(n), a2(n), a3(n), a4(n), x1, x2, x3, x4\n'
     '    complex(wp), intent(in) :: x(n)\n'
     '    complex(wp), intent(inout) :: y(n)\n'
     '    complex(wp), intent(out) :: t1, t2, t3, t4\n', 1),
    ('    real(wp) :: xi\n', '    complex(wp) :: xi\n', 1),
    (re.compile(r'(t\d) = \1 \+ (a\d\(i\)) \* xi'), r'\1 = \1 + conjg(\2) * xi', 4),
]

EDITS['real_symv_columns2'] = [
    ('    real(wp), intent(in) :: a1(n), a2(n), x1, x2, x(n)\n'
     '    real(wp), intent(inout) :: y(n)\n'
     '    real(wp), intent(out) :: t1, t2\n',
     '    complex(wp), intent(in) :: a1(n), a2(n), x1, x2, x(n)\n'
     '    complex(wp), intent(inout) :: y(n)\n'
     '    complex(wp), intent(out) :: t1, t2\n', 1),
    ('    real(wp) :: s1, s2\n', '    complex(wp) :: s1, s2\n', 1),
    (re.compile(r'([ts]\d) = \1 \+ (a\d\([^)]*\)) \* (x\([^)]*\))'),
     r'\1 = \1 + conjg(\2) * \3', 6),
]

EDITS['real_dot'] = [
    ('    real(wp), intent(in) :: x(n), y(n)\n'
     '    real(wp) :: dot\n',
     '    complex(wp), intent(in) :: x(n), y(n)\n'
     '    complex(wp) :: dot\n', 1),
    ('    real(wp) :: t1, t2, t3, t4\n', '    complex(wp) :: t1, t2, t3, t4\n', 1),
    (re.compile(r'(t\d) = \1 \+ x\(([^)]*)\) \* y\('), r'\1 = \1 + conjg(x(\2)) * y(', 5),
    ('  !Function real_dot computes the real inner product x^T*y of the first n\n'
     '  !elements with four independent partial sums.\n',
     '  !Function real_dot computes the Hermitian inner product x^H*y (the first\n'
     '  !argument is conjugated) of the first n elements with four independent\n'
     '  !partial sums.\n', 1),
    ('  !    x - Real vector of length n.\n', '  !    x - Complex vector of length n; conjugated.\n', 1),
]

# complex_dotu is a second copy of real_dot without the conjugation.
EDITS_DOTU = [
    ('real_dot', 'complex_dotu', 3),
    ('    real(wp), intent(in) :: x(n), y(n)\n'
     '    real(wp) :: dot\n',
     '    complex(wp), intent(in) :: x(n), y(n)\n'
     '    complex(wp) :: dot\n', 1),
    ('    real(wp) :: t1, t2, t3, t4\n', '    complex(wp) :: t1, t2, t3, t4\n', 1),
    ('  !Function complex_dotu computes the real inner product x^T*y of the first n\n'
     '  !elements with four independent partial sums.\n',
     '  !Function complex_dotu computes the unconjugated product x^TT*y of the\n'
     '  !first n elements of two complex vectors with four independent partial\n'
     '  !sums. It is the kernel of the back substitution with the rows of R.\n', 1),
    ('  !    dot - x^T*y.\n', '  !    dot - x^TT*y (no conjugation).\n', 1),
]

EDITS['real_axpy'] = [
    ('    real(wp), intent(in) :: alpha, x(n)\n'
     '    real(wp), intent(inout) :: y(n)\n',
     '    complex(wp), intent(in) :: alpha, x(n)\n'
     '    complex(wp), intent(inout) :: y(n)\n', 1),
]

for _name in ('real_add_leading_columns', 'real_add_columns'):
    EDITS[_name] = [
        ('    real(wp), intent(in) :: a(lda, *), coef(*)\n'
         '    real(wp), intent(inout) :: y(n)\n',
         '    complex(wp), intent(in) :: a(lda, *), coef(*)\n'
         '    complex(wp), intent(inout) :: y(n)\n', 1),
        ('    real(wp) :: c1, c2, c3, c4\n', '    complex(wp) :: c1, c2, c3, c4\n', 1),
    ]

EDITS['real_lu_factor'] = [
    ('    real(wp), intent(inout) :: a(lda, *)\n',
     '    complex(wp), intent(inout) :: a(lda, *)\n', 1),
    ('    real(wp) :: factor, pivot_magnitude, swap\n',
     '    complex(wp) :: factor, swap\n'
     '    real(wp) :: factor_magnitude, mag_swap, pivot_magnitude\n', 1),
    ('            swap = mag(k, j)\n'
     '            mag(k, j) = mag(ip, j)\n'
     '            mag(ip, j) = swap\n',
     '            mag_swap = mag(k, j)\n'
     '            mag(k, j) = mag(ip, j)\n'
     '            mag(ip, j) = mag_swap\n', 1),
    ('      do i = k + 1, n\n'
     '        a(i, k) = a(i, k) / a(k, k)\n'
     '      end do\n',
     '      do i = k + 1, n\n'
     '        a(i, k) = a(i, k) / a(k, k)\n'
     '      end do\n'
     '      !The moduli of the multipliers are kept in the part of column k of mag\n'
     '      !below the diagonal, which is not referenced again; the number of\n'
     '      !complex moduli evaluated is then of second order in n.\n'
     '      if (track) then\n'
     '        do i = k + 1, n\n'
     '          mag(i, k) = abs(a(i, k))\n'
     '        end do\n'
     '      end if\n', 1),
    ('        if (track) then\n'
     '          do i = k + 1, n\n'
     '            mag(i, j) = mag(i, j) + abs(a(i, k) * factor)\n'
     '          end do\n'
     '        end if\n',
     '        if (track) then\n'
     '          factor_magnitude = abs(factor)\n'
     '          do i = k + 1, n\n'
     '            mag(i, j) = mag(i, j) + mag(i, k) * factor_magnitude\n'
     '          end do\n'
     '        end if\n', 1),
    ('  !pivoting of a small dense real matrix in place',
     '  !pivoting of a small dense complex matrix in place', 1),
    ('  !    mag - Work array (lda by n) of the magnitude sums.\n',
     '  !    mag - Real work array (lda by n) of the magnitude sums.\n', 1),
]

EDITS['real_lu_solve'] = [
    ('    real(wp), intent(in) :: a(lda, *)\n', '    complex(wp), intent(in) :: a(lda, *)\n', 1),
    ('    real(wp), intent(inout) :: b(*)\n', '    complex(wp), intent(inout) :: b(*)\n', 1),
    ('    real(wp) :: swap, value\n', '    complex(wp) :: swap, value\n', 1),
]

EDITS['real_norm_squared'] = [
    ('    real(wp), intent(in) :: x(*)\n', '    complex(wp), intent(in) :: x(*)\n', 1),
    ('      norm_squared = norm_squared + x(i) * x(i)\n',
     '      norm_squared = norm_squared + real(x(i), wp)**2 + aimag(x(i))**2\n', 1),
    ('  !Function real_norm_squared computes the real Euclidean inner product\n'
     '  !\n'
     '  !                         x^T*x = sum(x(i)^2)\n'
     '  !\n'
     '  !over the first n elements of x. The result is the squared norm; the square\n'
     '  !root is not taken.\n',
     '  !Function real_norm_squared computes the Hermitian Euclidean inner product\n'
     '  !\n'
     '  !              x^H*x = sum(real(x(i))^2+imag(x(i))^2)\n'
     '  !\n'
     '  !over the first n elements of a complex vector. The mathematically real\n'
     '  !quantity is accumulated directly in real(wp). The result is the squared\n'
     '  !norm; the square root is not taken.\n', 1),
    ('  !    norm_squared - x^T*x for x(1:n).\n',
     '  !    norm_squared - The real value x^H*x for x(1:n).\n', 1),
]

EDITS['real_direction_difference'] = [
    ('    real(wp), intent(in) :: x(*), alpha, y(*), x_norm_squared\n',
     '    complex(wp), intent(in) :: x(*), alpha, y(*)\n'
     '    real(wp), intent(in) :: x_norm_squared\n', 1),
    ('    real(wp) :: difference, difference_norm_squared\n',
     '    complex(wp) :: difference\n'
     '    real(wp) :: difference_norm_squared\n', 1),
    ('                                difference * difference\n',
     '                                real(difference, wp)**2 + &\n'
     '                                aimag(difference)**2\n', 1),
    ('  !direction between two real vectors,\n',
     '  !direction between two complex vectors,\n', 1),
    ('  !In inverse iteration alpha=(x^T*y)/(y^T*y), so alpha*y is the component of\n'
     '  !x parallel to the previous iterate y. The result measures only the remaining\n'
     '  !directional change and is unaffected by real rescaling of an eigenvector.\n',
     '  !In inverse iteration alpha=(y^H*x)/(y^H*y), so alpha*y is the component of\n'
     '  !x parallel to the previous iterate y. The result measures only the remaining\n'
     '  !directional change and is unaffected by a complex rescaling (in particular\n'
     '  !a change of phase) of an eigenvector. Each squared magnitude is accumulated\n'
     '  !in real(wp).\n', 1),
    ('  !    x              - New real iterate, containing at least n elements.\n'
     '  !    alpha          - Scalar projection coefficient multiplying y.\n'
     '  !    y              - Previous real iterate, containing at least n elements.\n'
     '  !    x_norm_squared - Precomputed positive value x^T*x.\n',
     '  !    x              - New complex iterate, containing at least n elements.\n'
     '  !    alpha          - Complex projection coefficient multiplying y.\n'
     '  !    y              - Previous complex iterate, containing at least n\n'
     '  !                     elements.\n'
     '  !    x_norm_squared - Precomputed positive real value x^H*x.\n', 1),
]


COMMENT_EDITS['real_factor_serial'] = [
    ('  !explicit Q, ships the rows of Q, and releases its temporaries (one n by n\n'
     '  !array, tau, one row buffer, and the LAPACK workspace). The other ranks only\n'
     '  !receive. With one process no message is sent.\n',
     '  !explicit Q, ships the complex conjugates of the rows of Q, and releases its\n'
     '  !temporaries (one n by n array, tau, one row buffer, and the LAPACK\n'
     '  !workspace). The other ranks only receive. With one process no message is\n'
     '  !sent. The diagonal of the shifted matrix is formed from the real parts of\n'
     '  !the diagonals of H and S.\n', 1),
]

COMMENT_EDITS['real_factorize_fresh'] = [
    ('  !columns. In both forms, only the lower triangles of the leading n by n\n'
     '  !principal blocks are referenced.\n',
     '  !columns. In both forms, only the lower triangles of the leading n by n\n'
     '  !principal blocks are referenced; the upper triangle of the shifted matrix\n'
     '  !is generated by conjugating its lower triangle. Each diagonal element of H\n'
     '  !and of S must be real to within\n'
     '  !\n'
     '  !   QR_DIAGONAL_TOLERANCE*max(1,max(abs(Re)+abs(Im)) over the triangle) .\n'
     '  !\n'
     '  !H and S are tested separately, so that cancellation cannot make two\n'
     '  !non-Hermitian inputs acceptable. Roundoff-sized imaginary parts of the\n'
     '  !diagonals are accepted and discarded.\n', 1),
    ('  !            QR_ERR_INVALID_ARGUMENT when active_order is present but is\n'
     '  !            not positive;\n',
     '  !            QR_ERR_INVALID_ARGUMENT when active_order is present but is\n'
     '  !            not positive, or when a diagonal element of H or of S is\n'
     '  !            not real within the tolerance above;\n', 1),
]

COMMENT_EDITS['real_replace_symmetric'] = [
    ('  !matrix becomes the symmetric matrix whose row and column idx equal\n'
     '  !\n'
     '  !                 a = h_col - shift*s_col ,\n'
     '  !\n'
     '  !all other entries being unchanged. Nothing in Q or R changes. The routine\n'
     '  !is collective.\n',
     '  !matrix becomes the Hermitian matrix whose column idx equals\n'
     '  !\n'
     '  !                 a = h_col - shift*s_col ,\n'
     '  !\n'
     '  !and whose row idx is the complex conjugate of a, all other entries being\n'
     '  !unchanged. The arguments are conceptual columns: h_col(i) is the matrix\n'
     '  !element <phi_i|H|phi_idx> with the row function as the bra, whether i is\n'
     '  !smaller or larger than idx. The diagonal elements h_col(idx) and\n'
     '  !s_col(idx) must be real to within QR_DIAGONAL_TOLERANCE times\n'
     '  !max(1,max(abs(Re)+abs(Im))) of their column; their imaginary parts are\n'
     '  !discarded. Nothing in Q or R changes. The routine is collective.\n', 1),
    ('  !    h_col - New column idx of H, length order().\n'
     '  !    s_col - New column idx of S, length order().\n',
     '  !    h_col - New conceptual column idx of H, length order().\n'
     '  !    s_col - New conceptual column idx of S, length order().\n', 1),
    ('  !            QR_ERR_INVALID_ARGUMENT when idx is out of range or is not an\n'
     '  !            active border entry;\n',
     '  !            QR_ERR_INVALID_ARGUMENT when idx is out of range, or when the\n'
     '  !            diagonal element of a column is not real within the\n'
     '  !            tolerance above, or when idx is not an active border entry;\n', 1),
]

COMMENT_EDITS['real_append_symmetric'] = [
    ('  !    h_col - Column order()+1 of H, length order()+1 (the last element is\n'
     '  !            the new diagonal element).\n'
     '  !    s_col - Column order()+1 of S, length order()+1.\n',
     '  !    h_col - Conceptual column order()+1 of H, length order()+1: element i\n'
     '  !            is <phi_i|H|phi_new> with the row function as the bra. The\n'
     '  !            last element is the new diagonal element; it must be real\n'
     '  !            to within QR_DIAGONAL_TOLERANCE times\n'
     '  !            max(1,max(abs(Re)+abs(Im))) of the column, and its imaginary\n'
     '  !            part is discarded.\n'
     '  !    s_col - Conceptual column order()+1 of S, length order()+1, with the\n'
     '  !            same convention.\n', 1),
    ('  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from\n'
     '  !            order()+1;\n',
     '  !            QR_ERR_DIMENSION_MISMATCH when a column length differs from\n'
     '  !            order()+1;\n'
     '  !            QR_ERR_INVALID_ARGUMENT when the new diagonal element of H or\n'
     '  !            of S is not real within the tolerance above;\n', 1),
]

COMMENT_EDITS['real_solve'] = [
    ('  !Subroutine real_solve finds one eigenvalue and its eigenvector for the real\n'
     '  !generalized symmetric eigenvalue problem\n',
     '  !Subroutine real_solve finds one eigenvalue and its eigenvector for the\n'
     '  !complex generalized Hermitian eigenvalue problem\n', 1),
    ('  !at the cost of one additional symmetric product. The reason is numerical.\n',
     '  !at the cost of one additional Hermitian product; mu is real (the real\n'
     '  !parts of the two computed Hermitian forms are used). The reason is numerical.\n', 1),
    ('  !physical pair (H,S). The new vector is divided by max(abs(x)). Its change\n'
     '  !of direction relative to v is measured by\n'
     '  !\n'
     '  !             alpha   = (x^T*v)/(v^T*v),\n',
     '  !physical pair (H,S). The new vector is divided by its largest absolute\n'
     '  !real or imaginary component. Its change of direction relative to v is\n'
     '  !measured, independently of the complex phase of the iterates, by\n'
     '  !\n'
     '  !             alpha   = (v^T*x)/(v^T*v),\n', 1),
    ('  !             lambda       = (x^T*hs)/(x^T*ss),\n',
     '  !             lambda       = Re(x^T*hs)/Re(x^T*ss),\n', 1),
    ('  !                other: max(abs(x)) = 1.\n',
     '  !                other: largest absolute real or imaginary component of\n'
     '  !                       x equal to one.\n'
     '  !                The complex phase of x is not normalized.\n', 1),
    ('  !    lambda    - The Rayleigh-quotient eigenvalue approximation.\n',
     '  !    lambda    - The real Rayleigh-quotient eigenvalue approximation.\n', 1),
]

COMMENT_EDITS['real_physical_quotient'] = [
    ('  !        hx = H*x,  sx = S*x,  overlap = x^T*sx,  lambda = (x^T*hx)/overlap,\n',
     '  !        hx = H*x,  sx = S*x,  overlap = Re(x^T*sx),\n'
     '  !        lambda = Re(x^T*hx)/overlap,\n', 1),
    ('  !    overlap  - x^T*S*x.\n', '  !    overlap  - Real part of x^T*S*x.\n', 1),
]

COMMENT_EDITS['real_symmetric_product'] = [
    ('  !Subroutine real_symmetric_product computes y = A*x for a replicated real\n'
     '  !symmetric matrix of which only the lower triangle is referenced, with the\n',
     '  !Subroutine real_symmetric_product computes y = A*x for a replicated complex\n'
     '  !Hermitian matrix of which only the lower triangle is referenced (the\n'
     '  !strictly upper triangle is taken to be the conjugate transpose of the\n'
     '  !strictly lower one, and the imaginary parts of the diagonal elements are\n'
     '  !ignored), with the\n', 1),
    ('  !residual checks, and to the timing driver. The routine is collective.\n',
     '  !residual checks, and to the timing driver. The routine is collective. The\n'
     '  !binding keeps the name symmetric_product of the real state.\n', 1),
]

COMMENT_EDITS['real_shifted_quotient'] = [
    ('  !this way from the stored shifted matrix. For a 2850-function basis in\n'
     '  !extended precision the plain quotient was measured to be wrong by up to 15\n'
     '  !units of the last place and the shifted quotient by at most one half,\n'
     '  !which matters to an optimizer that resolves energy differences of that\n'
     '  !size.\n',
     '  !this way from the stored shifted matrix. For a 2850-function real basis in\n'
     '  !extended precision the plain quotient of the real state was measured to be\n'
     '  !wrong by up to 15 units of the last place and the shifted quotient by at\n'
     '  !most one half, which matters to an optimizer that resolves energy\n'
     '  !differences of that size. Both quadratic forms are real: the strictly\n'
     '  !upper triangles are taken to be the conjugate transposes of the strictly\n'
     '  !lower ones and the imaginary parts of the diagonal elements are ignored,\n'
     '  !as in symmetric_product.\n', 1),
    ('    !columns owned by this rank. A column contributes its diagonal term and\n'
     '    !twice its part below the diagonal.\n',
     '    !columns owned by this rank. A column contributes its real diagonal term\n'
     '    !and twice the real part of its part below the diagonal.\n', 1),
]

COMMENT_EDITS['real_shifted_column'] = [
    ('  !the lower triangles to the quadratic forms x^T*(H-shift*S)*x and x^T*S*x:\n'
     '  !the diagonal term plus twice the part below the diagonal.\n',
     '  !the lower triangles to the real quadratic forms x^T*(H-shift*S)*x and\n'
     '  !x^T*S*x. With the stored elements a(i,j), i>=j, of a Hermitian matrix,\n'
     '  !column j contributes\n'
     '  !\n'
     '  !        Re(a(j,j))*abs(x(j))**2 + 2*Re(x(j)*sum_{i>j} a(i,j)*conj(x(i))) ,\n'
     '  !\n'
     '  !the second term being the element (i,j) and its conjugate (j,i) taken\n'
     '  !together.\n', 1),
    ('  !    shift - The shift.\n', '  !    shift - The real shift.\n', 1),
]

COMMENT_EDITS['real_factorization_residual'] = [
    ('  !    v - Nonzero real probe vector with exactly n elements. A fixed dense\n',
     '  !    v - Nonzero complex probe vector with exactly n elements. A fixed dense\n', 1),
]

COMMENT_EDITS['real_install_column'] = [
    ('  !values on the border indices go to row and column of the border-border\n'
     '  !block, Gs is obtained by one distributed base solve, and the row and column\n',
     '  !values on the border indices go to the column, and their complex\n'
     '  !conjugates to the row, of the Hermitian border-border block, Gs is\n'
     '  !obtained by one distributed base solve, and the row and column\n', 1),
    ('  !decoupling column: zero off the diagonal and the base diagonal element on\n'
     '  !the diagonal (one if that element is zero).\n',
     '  !decoupling column: zero off the diagonal and the real part of the base\n'
     '  !diagonal element on the diagonal (one if that value is zero).\n', 1),
    ('  !    a     - Shifted column of length n (current order); not referenced\n'
     '  !            when dummy is true.\n',
     '  !    a     - Shifted conceptual column of length n (current order): a(i)\n'
     '  !            is the element (i,idx). Not referenced when dummy is true.\n', 1),
]

COMMENT_EDITS['real_fill_schur'] = [
    ('  !  Input parameter:\n'
     '  !    p    - Border position of the entry.\n',
     '  !BG and BZ hold inner products that conjugate the border columns Bt, and C\n'
     '  !is Hermitian. The Schur complement is treated as a general complex matrix:\n'
     '  !the base elements recovered from the factors are Hermitian only to\n'
     '  !roundoff.\n'
     '  !\n'
     '  !  Input parameter:\n'
     '  !    p    - Border position of the entry.\n', 1),
]

COMMENT_EDITS['real_base_solve'] = [
    ('  !  1. y = Q^T*b: every rank accumulates b(i) times its rows of Q (contiguous\n'
     '  !     columns of qt), followed by one reduction of n0 values.\n',
     '  !  1. y = Q^T*b: every rank accumulates b(i) times the complex conjugates\n'
     '  !     of its rows of Q (the contiguous columns of qt, which stores columns\n'
     '  !     of Q^T), followed by one reduction of n0 complex values.\n', 1),
]

COMMENT_EDITS['real_base_column'] = [
    ('  !entries of column idx of R (one reduction) and then the products of its\n'
     '  !rows of Q with that column (one reduction).\n',
     '  !entries of column idx of R (one reduction) and then the products of its\n'
     '  !rows of Q, the conjugates of the stored columns of qt, with that column\n'
     '  !(one reduction).\n', 1),
]

COMMENT_EDITS['real_sym_product'] = [
    ('  !for 80-bit extended precision (QR_PRODUCT_PAIRS), where pairs keep all\n'
     '  !running values in the eight registers of the x87 unit and were measured\n'
     '  !to be faster.\n',
     '  !for 80-bit extended precision (QR_PRODUCT_PAIRS), where the pairs of the\n'
     '  !real kernel are kept: for complex data no significant difference between\n'
     '  !pairs and groups of four was measured, and single columns were slower.\n', 1),
    ('  !          triangle holds the matrix. The upper triangle is not referenced.\n',
     '  !          triangle holds the matrix. The upper triangle is not referenced,\n'
     '  !          and the imaginary parts of the diagonal elements are ignored.\n', 1),
]

EXTRA_PROCEDURES = '''\
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
'''

# --------------------------------------------------------------------------
# Assemble
# --------------------------------------------------------------------------

i_module = find('module qrlinalg')
i_pivot = find('    !Relative pivot threshold of the LU factorization of Scc.')
i_public = find('  public :: qr_real_state')
i_model = find('    !see real_factor_select; a negative first element means not calibrated.')
i_model_public = find('  public :: real_factor_model')
i_iface = find('  end interface')
i_type0 = find('  type :: qr_real_state')
i_type1 = find('  end type qr_real_state')
i_tail0 = find('  subroutine clear_real_state(self)')
i_end = find('end module qrlinalg')

# Type definition.
type_text = apply_edits('type', '\n'.join(src[i_type0:i_type1 + 1]) + '\n', TYPE_EDITS)
type_lines = generic(type_text.rstrip('\n').split('\n'))

# Procedures.
start_re = re.compile(r'^  (?:pure )?(?:subroutine|function) (\w+)')
units = []
i = i_tail0
gap = []
while i < i_end:
    m = start_re.match(src[i])
    if m:
        name = m.group(1)
        j = i
        while not re.match(r'^  end (?:subroutine|function) %s$' % name, src[j]):
            j += 1
        units.append(('gap', None, gap))
        units.append(('proc', name, src[i:j + 1]))
        gap = []
        i = j + 1
    else:
        gap.append(src[i])
        i += 1
trailing_gap = gap

used = set()
tail = []
skip_gap = False
for kind, name, lines in units:
    if kind == 'gap':
        if not skip_gap:
            tail += generic(lines)
        skip_gap = False
        continue
    if name.startswith('qr_'):
        # Type independent; shared with the real part.
        skip_gap = True
        continue
    body = '\n'.join(lines) + '\n'
    if name in EDITS:
        body = apply_edits(name, body, EDITS[name])
        used.add(name)
    if name in COMMENT_EDITS:
        body = apply_edits(name + ' (comments)', body, COMMENT_EDITS[name])
        used.add(name + '#c')
    tail += generic(body.rstrip('\n').split('\n'))
    if name == 'real_dot':
        body = apply_edits('complex_dotu', '\n'.join(lines) + '\n', EDITS_DOTU)
        tail += [''] + generic(body.rstrip('\n').split('\n'))
missing = set(EDITS) - used
if missing:
    raise SystemExit('edits for unknown procedures: %s' % sorted(missing))
missing = set(k + '#c' for k in COMMENT_EDITS) - used
if missing:
    raise SystemExit('comment edits for unknown procedures: %s' % sorted(missing))
while tail and tail[0] == '':
    tail.pop(0)
tail += [''] + EXTRA_PROCEDURES.rstrip('\n').split('\n') + ['']

out = []
out += src[:i_module] + text_lines(HEADER_TEXT)
out += src[i_module:i_pivot + 1] + text_lines(CONSTANT_TEXT)
out += src[i_pivot + 1:i_model + 1] + text_lines(MODEL_TEXT)
out += src[i_model + 1:i_public + 1] + ['  public :: qr_complex_state']
out += src[i_public + 1:i_model_public + 1] + ['  public :: complex_factor_model']
out += src[i_model_public + 1:i_iface + 1] + [''] + text_lines(INTERFACE_TEXT)
out += src[i_iface + 1:i_type1 + 1] + [''] + type_lines
out += src[i_type1 + 1:i_end] + text_lines(TAIL_BANNER) + [''] + tail
out += src[i_end:]

long_lines = [k + 1 for k, line in enumerate(out) if len(line) > 96]
if long_lines:
    print('lines longer than 96 characters:', long_lines)

if MODE in ('generate', 'regenerate'):
    open(OUT, 'w', newline='\n').write('\n'.join(out))
    print('wrote %s: %d lines (%d real source lines)' % (OUT, count(out), count(src)))
    if MODE == 'regenerate':
        if out == module_lines:
            print('%s was up to date: the output is identical to it' % MODULE)
        else:
            print('the complex part of %s was not up to date: the output differs '
                  'from it' % MODULE)
elif out == module_lines:
    print('%s: the complex part is up to date (%d lines, %d real source lines)'
          % (MODULE, count(out), count(src)))
else:
    print('%s: the complex part is NOT what the generator produces from the real '
          'part; run --regenerate and review the diff' % MODULE)
    if MODE == 'check':
        raise SystemExit(1)
