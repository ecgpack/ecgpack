# Q-method design and maintenance guide

## Purpose and status

This document is the design reference for generalized-eigenvalue method `Q` in
ECGPACK. It describes the numerical contract, the canonical matrix layout, the
QR state (a distributed base factorization plus a batch of deferred updates),
the MPI rules, the supported BBOP commands, the measured cost, and the
procedure for porting the method to another ECGPACK variant.

All seven energy codes (`RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `PG_0S`, `PG_1P`,
and `CG_0S`) implement the design described here, with the same
`src/qrlinalg.f90`. `RG_0S` is the reference implementation: the Q helper
routines of `RG_1P`, `RG_2D`, and `RG_2P` are its text line for line, those
of `PG_0S` and `PG_1P` are its text in the source style of those files, and
`CG_0S` applies the Hermitian rules of section 11 to them and uses the
complex state type. Statements about `RG_0S` routines therefore hold for the
other codes unless section 11 says otherwise. Section 16 is the procedure by
which the design was carried from `RG_0S` to the other six codes; it remains
the procedure for a later change of the helper layer.

The previous implementation ("previous" throughout this document) held
explicit Q and R factors on MPI rank zero and updated them by Givens
rotations for every trial point, using the bundled `src/qrupdate/` library.
No energy code calls that library any more (section 13).

Why the redesign: on a 2850-function, four-particle basis in extended
precision with 14 MPI ranks the previous implementation needed 527 ms per
energy evaluation (method I: 20 ms), because every trial applied a full
symmetric rank-two update to dense factors and all QR work ran on one rank.
The present design needs 25 ms on the same problem (22 to 29 ms, depending
on how many functions are pending), and its fresh factorization takes 2.4 s
on the 14 ranks instead of 40 s on one (section 12).

The Q algorithms are not implemented in the off-diagonal matrix-element codes.

The following commands accept `Q` in the energy codes:

| Command | Q behavior |
| --- | --- |
| `BASIS_ENL` | Appends candidate functions in canonical order as deferred columns. |
| `OPT_CYCLE` | Optimizes arbitrary canonical functions through deferred column replacements. |
| `FULL_OPT1` | Optimizes an arbitrary contiguous function range without permuting the physical basis. |
| `ELIM_LCFN` | Deletes the coefficient-selected canonical functions and factorizes the survivors. |
| `ELIM_LND1` | Deletes the overlap-selected canonical functions and factorizes the survivors. |
| `SEPR_LND1` | Replaces perturbed canonical columns selected by the overlap rule. |
| `SEPR_FLCF` | Replaces perturbed canonical columns selected by the coefficient rule. |
| `EXPC_VALS` | Constructs and solves the Q eigenproblem before evaluating operators. |
| `DENSITIES` | Uses the Q eigenpair before computing position-space distributions. |
| `MOMT_DENS` | Uses the Q eigenpair in the momentum-density path. |
| `SAVE_HSWF` | Saves physical unshifted H/S, the eigenvector, and/or the wave function. |

`SAVE_FILE` has no solver argument. It requires no QR operation because Q
keeps basis arrays in their user-visible canonical order.

## 1. Numerical semantics

Method Q solves the same physical generalized symmetric eigenproblem as the
existing methods:

```text
H c = E S c .
```

Its QR state represents the shifted matrix

```text
A = H - sigma*S ,
```

where `sigma` is the current inverse-iteration shift. The shift remains fixed
between fresh factorizations. Any shift change requires a fresh factorization.

Q is a shift-targeted inverse-iteration method. It follows method I's
root-near-shift semantics:

- `CURRENT_ENERGY`, `INVITPARAMETER`, and `EIGVAL_TOLERANCE` control the
  target and convergence; `EIGVAL_TOLERANCE` bounds the relative change of the
  eigenvector in one direct inverse iteration, exactly as in method I
  (section 9);
- `WHICH_EIGENVALUE` is not a root-selection guarantee for Q; and
- an inappropriate shift can converge to a different nearby root.

Method G remains the safe reference when a calculation must select a root by
ordered eigenvalue index. Comparisons between G and Q are meaningful only
when both methods target the same physical root.

The returned coefficient vector is S-normalized:

```text
c^T S c = 1.
```

This normalization is required by the existing ECGPACK energy-gradient and
observable formulas.

The energy is the shifted Rayleigh quotient of the returned vector,

```text
E = sigma + c^T (H - sigma*S) c / (c^T S c) ,
```

with every element of `H - sigma*S` formed before it is multiplied by the
coefficients, which is how method I obtains its energy from its stored
shifted matrix. This is not a detail: section 9 explains why the plain
quotient is not accurate enough for the optimizer.

### Order-one special case

For one basis function, `SolveQ` uses the exact result

```text
E = H(1,1)/S(1,1),
c(1) = 1/sqrt(S(1,1)).
```

This avoids unnecessary inverse iteration and prevents catastrophic
cancellation when a new calculation still carries the conventional enormous
`CURRENT_ENERGY` sentinel. A one-function problem never uses its factors:
`ApplyQTrial` stores the column and factorizes fresh, and `FactorizeQFresh`
does not probe it.

### Shift retargeting in BASIS_ENL

The shift of a command is `CURRENT_ENERGY*INVITPARAMETER`, set by `main`
before every BBOP step. `BASIS_ENL Q` alone may move it during the command,
at a function boundary: after a block has been accepted and optimized,
`BasisEnlQ` calls `RetargetQStaleShift` with the energy of the represented
problem. The routine sets the shift to that energy times `INVITPARAMETER`
and factorizes fresh when

- the shift is **provisional** (the command started from an empty basis, see
  below), or
- the shift is **stale**: the solves of the block just finished (candidates,
  the discrete index or premultiplier scan of the codes that have one, and
  optimizer evaluations, counted by `Glob_InvItTempCounter1` and
  `Glob_InvItTempCounter2`, which `BasisEnlQ` resets for every block) needed
  more than `Q_STALE_SHIFT_PASSES` = 4 passes of `SolveQ` on average.

Inverse iteration converges per pass by the ratio of the distances from the
shift to the wanted and to the nearest other eigenvalue. Near a well placed
shift a solve needs 1.1 to 1.9 passes; many passes mean that the energy has
moved far from the shift, which happens while a basis grows from a few
functions and the energy drops by much more than the distance
`|E|*|INVITPARAMETER-1|` that the shift keeps from it. Section 12 gives the
measurements.

The decision uses replicated integer counters only, so every rank takes the
same branch; the call is collective. A shift that already has the target
value is left alone. Nothing else needs to know: `Glob_ApproxEnergy` is read
by the Q helpers (which compare it with the shift of the QR state) and by the
method-I routines, `main` derives it anew before the next BBOP step, the
input/output file holds `CURRENT_ENERGY` and `INVITPARAMETER` but not the
shift, and the swap file holds unshifted matrices. `OPT_CYCLE Q` and
`FULL_OPT1 Q` keep their shift for the whole command.

### Starting from an empty basis

An empty input basis carries a placeholder `CURRENT_ENERGY` (conventionally
1E30), so `main` cannot derive a meaningful shift. `BASIS_ENL Q` handles this
as follows.

1. While the first block is being selected and optimized, the shift is
   *provisional* (`Q_Workspace%ShiftIsProvisional`). `AppendQTrial` detects a
   remote shift (magnitude above 1E6 times the largest diagonal energy) and
   sets the shift of every first-block candidate from that candidate's own
   lowest normalized diagonal energy times `INVITPARAMETER`. Each such
   candidate is factorized fresh anyway. For the usual `Kstep=1` start the
   one-function problem is solved in closed form and the shift is not used.
2. After the first block has been accepted and optimized,
   `RetargetQStaleShift` ends the provisional phase: the shift becomes the
   accepted energy times `INVITPARAMETER`.
3. The following blocks lower the energy by large steps, so the stale-shift
   rule retargets after most of the first functions, and after about every
   second one while the basis has a few tens of functions (section 12). It
   stops once a function no longer moves the energy by a sizable fraction of
   the distance to the next eigenvalue.

Method I has no counterpart: it keeps the shift derived from the placeholder
and cannot start from an empty basis (it returns energies such as 0.0 for the
first functions). For a nonempty basis the RG and CG codes keep the shift of
method I for the whole command; `PG_0S` and `PG_1P` re-anchor it without a
test before every block but the first of `BASIS_ENL` and before every window
of `OPT_CYCLE` (`RefreshINVITShift` in their `matform.f90`), which is why
the comment of `RetargetQStaleShift` on method I differs in those two codes.

Limits that remain: the shift is retargeted between functions, not within
one. The candidates of the second function are still solved with the shift
of the one-function energy (6 passes per solve for He); with `Kstep>1` the
first block targets the root nearest its lowest diagonal energy, which is
the ground state only if the block is well conditioned; and if an excited
state comes closer to the shift than the ground state within one function,
the solve can follow the wrong root. That case is reported as `EC0128` (no
candidate improves the energy), not as a wrong energy, because every other
root lies above the accepted energy.

## 2. Canonical basis and matrix layout

Q never permutes the physical basis merely to make active functions trailing
rows or columns. Canonical basis index `i` always identifies the same entries
in every basis-indexed object:

```text
Glob_NonlinParam(:,i)   nonlinear parameters
Glob_FuncNum(i)         function identity and history order
Glob_c(i)               normalized linear coefficient
Glob_diagS(i)           raw diagonal overlap <phi_i|phi_i>
basis metadata(i)       angular prefactor index or indices, when present
Glob_H and Glob_S       physical normalized matrices
qrlinalg state          the identical canonical row/column order
```

This stable-order rule is the central Q invariant. The G implementation may
move active functions into a trailing block for Cholesky-based work; that
permutation must not be copied into Q. If a user-visible operation genuinely
reorders a Q basis, all basis-indexed arrays must be permuted consistently and
the QR state must be rebuilt.

### Authoritative matrix triangle

The lower triangles of `Glob_H` and `Glob_S`, including their diagonals, are
authoritative. They contain normalized, unshifted physical matrix elements:

```text
Glob_H(i,i) = <phi_i|H|phi_i> / Glob_diagS(i)
Glob_S(i,i) = 1

Glob_H(i,j) = <phi_i|H|phi_j> /
              sqrt(Glob_diagS(i)*Glob_diagS(j)),  i > j

Glob_S(i,j) = <phi_i|phi_j> /
              sqrt(Glob_diagS(i)*Glob_diagS(j)),  i > j
```

The upper triangles are scratch or unspecified storage and must not be read as
physical Q matrices. Access a conceptual symmetric element as

```text
A(max(i,j),min(i,j)).
```

`QCanonicalMatrixElement`, `GatherQCanonicalColumn`, and
`StoreQCanonicalColumn` centralize this rule. qrlinalg reads only the lower
triangles, on every rank.

### Why H and S remain unshifted

`Glob_H` must never be overwritten with `H-sigma*S` in a Q calculation.
Keeping physical H and S separate is required for:

- changing the inverse-iteration shift;
- evaluating physical eigenpair residuals;
- saving H and S through `SAVE_HSWF`;
- expectation values and overlap penalties;
- reconstructing fresh factors after a flush or a failure; and
- comparing Q against G or I without undoing an in-place transformation.

Only the private qrlinalg state represents the shifted matrix.

## 3. The QR state: base factorization plus batch

`src/qrlinalg.f90` represents `A = H - sigma*S` of the current order `n` by

1. a **base factorization** `A_base = Q*R` of order `n0`, computed at a
   snapshot and never modified afterwards; and
2. a **batch** (in the module: border) of at most `max_pending` entries, one
   for every *pending* function, that is, a function whose current row and
   column differ from the snapshot or that was appended after it, and one for
   every trailing base slot deleted by `truncate`.

Replacing or appending a function costs one solve with the base factors and
O(batch*n) local work. Nothing in Q or R changes. The only O(n^3) operation
is a fresh factorization ("flush").

### 3.1 Bordered system

Let D0 be the set of base slots (indices not above `n0`) that own a batch
entry, and I the clean set of the other base indices. On the base indices the
current matrix is a symmetric low-rank modification of `A_base`, and the
appended functions border it. Both are written as one augmented system on the
full, unmodified base matrix:

```text
[ A_base  E    D    B  ] [x0]   [r0]     E   unit vectors e_k, k in D0
[ E^T     0   -I    0  ] [yE] = [0 ]     D   (current - base) column of each D0
[ D^T    -I    DD   K  ] [yD]   [0 ]         member on the rows of I, zero on D0
[ B^T     0    K^T  CA ] [xA]   [rA]     B   current column of each appended
                                             function on the rows of I
DD  (current - base) values of the D0 x D0 block
K   current values of the D0 x appended block
CA  current values of the appended x appended block
```

Eliminating `yE` and `yD` shows that `(x0, xA)` solves the current system. A
deleted trailing slot is a member of D0 whose current row and column are zero
except for the diagonal, which decouples it. Stored, replicated on all ranks:

```text
Zb(:,q)  = A_base^{-1} e_{D0(q)}        per member of D0
Ab(:,q)  = A_base e_{D0(q)}             base column, recovered as Q*(R*e)
Bt(:,s)  column of D or of B            per batch entry
Gs(:,s)  = A_base^{-1} Bt(:,s)          per batch entry
C, BG = Bt^T Gs, BZ = Bt^T Zb           small matrices
Sc = Caug - [E D B]^T A_base^{-1} [E D B]   Schur complement
```

A solve `A x = r` is one base solve `g = A_base^{-1} r0`, a small system with
`Sc`, and `x0 = g - Zb*yE - Gs*[yD; xA]`.

Why the base matrix is modified and not restricted: excluding the rows and
columns of D0 from the base and bordering the clean block was implemented
first and is algebraically equivalent, but the inverse of the clean block is
then a difference of terms that exceed the result by the ratio of the
near-singularities of the two matrices, whatever the size of the pending
changes (bordered-solve backward errors of 1E5 to 1E7 epsilon in the module's
unit tests). In the scheme above the cancelling terms are proportional to the
accumulated *change* of the pending rows, which is small for optimizer steps.

### 3.2 Active tail

The last entries of the batch are the **active tail**: the functions being
changed now. All other entries and all unit-vector unknowns are committed.
`LU(Scc)` of the committed block is computed when the committed set changes;
the tail quantities are recomputed lazily, at O(batch^2) cost, at the first
solve after a tail column changed.

The module never guesses the tail. workproc announces the active set with
`begin_active(indices)`. That call commits the previous tail if the set
differs, moves the announced functions to the end of the batch (a clean base
function joins D0 with a zero modification: one base solve, one column
recovery, one distributed rebuild of the small matrices), and is a no-op when
the set is unchanged. When the batch cannot take the set it returns
`need_flush` and changes nothing.

`BASIS_ENL` evaluates hundreds of candidates for the same appended slot:
`truncate` drops the previous candidate without touching the committed set,
`begin_active` announces the slot again, `append_symmetric` supplies the
candidate. No O(batch^3) or O(batch^2*n) work happens per candidate. This
also holds when the slot is a base slot because a flush included the
candidate: the deleted slot keeps its tail entry.

### 3.3 Inverse iteration in the module

The first iteration of a `solve` call is a direct inverse iteration,
`x = A^{-1} S v`. Later iterations of the same call, and the optional steps
enabled by `residual_goal`, use the residual-correction form

```text
mu = (v^T H v)/(v^T S v),   x = v - A^{-1} (H v - mu*S v) ,
```

which gives the same direction in exact arithmetic. The correction form makes
the attainable residual independent of the error growth of block elimination
over the batch. Its weakness, found on the 2850-function problem, is the
subject of section 9: workproc drives the iteration itself, one `solve` call
per iteration.

After the iteration `solve` forms `H x` and `S x` once and returns the plain
Rayleigh quotient and the relative physical residual.

### 3.4 Distribution

With P ranks of `MPI_COMM_WORLD`:

- Rows of Q and rows of R are distributed cyclically: row `i` lives on rank
  `mod(i-1,P)`, stored transposed so that a row is contiguous.
- The P by P diagonal blocks of R are replicated.
- A base solve `x = R^{-1} Q^T b` is a local product, one reduction of `n0`
  values, and a blocked back substitution with one reduction of at most P
  values per block.
- Symmetric products with the replicated lower triangles assign the columns
  cyclically; one reduction of `n` values.
- All O(batch*n) and smaller work is executed identically on every rank on
  replicated data.
- With one process no MPI routine is called with `wp` data, so the serial
  path works in quadruple precision with gfortran.

**Every procedure of the state type is collective.** All ranks call it with
identical replicated arguments and obtain identical results and status
values. Rank-local failures (allocation) are reduced before any rank branches.

### 3.5 Fresh factorization

Two algorithms compute the base factors and fill the same distributed layout.

- **Serial** (`real_factor_serial`): rank zero forms the shifted matrix in
  one temporary n by n array, calls the bundled `DGEQRF`, ships the rows of
  R, calls `DORGQR`, and ships the rows of Q. Its cost does not decrease with
  the number of ranks.
- **Distributed** (`real_factor_distributed`): a blocked Householder QR
  factorization on the distributed rows. Every rank forms its own rows of
  `H - sigma*S` from the replicated lower triangles, so no rank holds an n by
  n temporary (the work arrays take about `2*(1+1/P)*32*n` elements per
  rank). The columns are processed in panels of `QR_FACTOR_BLOCK` = 32: a
  panel is factorized with one small reduction per column, the reflectors of
  the panel are applied to the remaining columns in compact WY form
  (`A := A - (V*T^T)*(V^T*A)`, one reduction of `V^T*A` per panel), and Q is
  accumulated afterwards by applying the block reflectors in reverse order to
  the identity in the same layout. The reflectors follow the conventions and
  safeguards of `DLARFG`. Everything a rank branches on comes from a
  reduction or from replicated data. With one process no MPI routine is
  called. In 80-bit extended precision the two block kernels accumulate inner
  products in registers (`QR_FACTOR_DOT_KERNELS`); this shortened the
  factorization of order 2850 from 3.06 s to 2.44 s. The other precisions
  update whole columns, which the compiler vectorizes.

**Selection.** `factorize_fresh` asks `real_factor_select`:

1. the environment variable `ECG_QR_FACTOR` (`serial`, `distributed`, or
   `auto`; case-insensitive; read once per run on rank zero and broadcast),
   if it names one of the three;
2. otherwise the value set with `qr_set_factorization` (workproc does not
   call it, so this is `QR_FACTOR_AUTO`);
3. automatic choice: with one process the serial algorithm; with several,
   the faster one according to two timing models, `t_serial = a*n^3 + b*n^2`
   and `t_distributed = c*n^3 + d*n`, whose coefficients are measured once
   per run at the first fresh factorization (`real_factor_calibrate`: both
   algorithms on synthetic matrices of orders 128 and 256, or 64 and 128 in
   quadruple precision; rank zero fits and broadcasts the coefficients, so
   all ranks decide identically).

The calibration costs a few factorizations of order 256 (a fraction of a
second) per run. Because it is a timing, the order at which the automatic
choice changes from the serial to the distributed algorithm can differ from
run to run, and the two algorithms differ by roundoff; a calculation that
must be reproducible to the last bit across runs should set `ECG_QR_FACTOR`.
`factorization_kind()` reports which algorithm produced the current factors.

Measured (muon, 14 ranks, one node; module timing driver, synthetic
matrices): order 2850 in extended precision 40.5 s serial and 2.44 s
distributed; order 2850 in double precision 9.1 s and 0.50 s; order 2000 in
double precision with 2 and 4 ranks 2.9 s and 3.1 s serial against 0.86 s and
0.44 s distributed. The automatic choice was the distributed algorithm for
every order from 1000 upwards in these runs; in the unit tests (2 to 14
ranks) it changed from the serial to the distributed algorithm between
orders 32 and 128, not always at the same order.

## 4. qrlinalg contract

`src/qrlinalg.f90` is the API authority; read the header of the procedure you
call. A `qr_real_state` owns the distributed factors, the batch, and its
workspace. It does not own H, S, the coefficient vector, or ECG parameters,
and it depends only on `wp_def`.

The same file holds the complex Hermitian twin, type `qr_complex_state`, with
the same procedures (names `complex_*`), for `CG_0S`; the module header lists
its conventions (section 11 summarizes them). The file is byte-identical in
all seven energy codes: a real code simply does not use the complex type. The
complex part is generated from the real part, see "Maintaining the module" in
section 13.

| Operation | Contract |
| --- | --- |
| `initialize(capacity,max_pending,info)` | Allocates the state on every rank; no factorization. `max_pending` above the capacity is reduced to it. |
| `factorize_fresh(H,S,shift,info,active_order=N)` | Fresh base factorization of the leading block; empties the batch. Physical leading dimensions may exceed N. The algorithm is selected as described in section 3.5. |
| `qr_set_factorization(choice[,info])` (module procedure) | Requests `QR_FACTOR_AUTO` (default), `QR_FACTOR_SERIAL`, or `QR_FACTOR_DISTRIBUTED` for all states; overridden by `ECG_QR_FACTOR`. Must be called with the same value on all ranks. Not used by workproc. |
| `begin_active(indices,need_flush,info)` | Announces the active set (existing functions and/or the range `order()+1..`). `need_flush`: the batch cannot take it, state unchanged. `QR_ERR_SINGULAR`: the committed block cannot be factorized for this partition. `QR_ERR_BATCH_FULL`: the set is larger than the whole batch. |
| `replace_symmetric(index,h_column,s_column,info)` | Supplies the NEW complete conceptual column of an announced function (not a difference). |
| `append_symmetric(h_column,s_column,info)` | Adds function `order()+1` with its complete column. The order grows only here. |
| `truncate(target_order,info)` | Removes trailing functions. `QR_ERR_BATCH_FULL` (state unchanged) when the batch cannot take the removed base slots. |
| `delete_symmetric(index,info)` | Last index only; interior deletion is a flush. Not used by workproc. |
| `solve(H,S,v,x,lambda,tol,max_iter,norm_mode,rel_acc,num_iter,eig_residual,info[,residual_goal,num_polish])` | Inverse iteration (section 3.3). Reads the physical lower triangles, which must already hold every trial column. Success means the direction test was met; `QR_ERR_NO_CONVERGENCE` still returns the best vector. |
| `shifted_quotient(H,S,x,lambda,info)` | Shifted Rayleigh quotient with the elementwise shifted matrix (section 9). |
| `factorization_residual(H,S,v,abs,rel,info[,forward_error])` | Solves `A x = A v` with the state and returns the relative residual. Meaningful as a check of the factors only while nothing is pending. |
| `solve_shifted`, `symmetric_product` | One bordered solve; one distributed product. Tests and diagnostics. |
| queries | `is_valid`, `order`, `base_order`, `get_capacity`, `get_shift`, `pending_count`, `active_count`, `excluded_count`, `max_pending`, `get_update_count`, `get_updates_since_fresh`, `factorization_kind`. |
| `clear()` | Releases the state. No communication. |

`pending_count()` counts batch entries, including deleted trailing base slots.

The state components are private in production. Do not define
`QRLINALG_TESTING` in an ECGPACK build.

### Error boundary

The procedures validate indices, dimensions, and state before changing
anything; an argument error leaves the state unchanged. The shifted matrix is
deliberately nearly singular, so a small diagonal element of R or a small
pivot of the active Schur complement is accepted; only zero, subnormal, and
non-finite values are rejected. The committed block is tested with a relative
pivot threshold that detects pivots lost to cancellation. Whether an
eigenpair is usable is decided by the physical residual.

## 5. Q workspace and ownership

`QMethodWorkspace` is defined in each energy code's `workproc.f90`. It is
kept there because optimizer scheduling, trial acceptance, physical matrices,
and QR recovery form one lifecycle.

### MPI rules

- `Glob_H`, `Glob_S`, `Glob_diagS`, coefficients, nonlinear parameters, basis
  metadata, relationship flags, and the Q statistics are replicated.
- Matrix-element work remains distributed with the established MPI routing.
- The qrlinalg state exists on every rank. **No Q helper guards a qrlinalg
  call by the rank number, and none broadcasts its results.** A helper that
  did would deadlock or desynchronize the ranks.
- No rank may return early from a collective Q path on the basis of a local
  condition. All conditions that the helpers branch on are replicated.

### Relationship flags

The workspace tracks relationships that qrlinalg cannot infer:

```text
MatricesAreCanonical
FactorsMatchMatrices
MatrixParametersAreStored
TrialIsReady
TrialHasDerivatives
AcceptedPointIsStored
ShiftIsProvisional
```

`FactorsMatchMatrices` means that the QR state (base plus batch) represents
the leading `MatrixOrder` block of the physical matrices.

`MatrixParam` records the active nonlinear parameters represented by H/S.
`PreviousParam` belongs to the transaction rollback snapshot.
`AcceptedParam` and `AcceptedEnergy` belong to the optimizer's best accepted
point. These generations must not be conflated.

### Active-map layout

Optimizer variable block `a` is mapped explicitly:

```text
ActiveFunction(a) = canonical basis index
ActivePosition(i) = a, or zero when canonical function i is inactive
```

The active list must contain unique in-range indices. It may be ascending,
descending, contiguous, or noncontiguous. `SetQActiveFunctions` sets only
these maps; it may be called before the first factorization. The set is
announced to qrlinalg by `AnnounceQActiveFunctions` at the start of every
`ApplyQTrial` and `AppendQTrial`.

### Batch size

`QBatchSize(Capacity,MaxActive)` returns the `max_pending` that workproc
requests (P is the number of MPI ranks):

```text
default = min(max(32, Capacity/64), 256)
limit   = max(default, Capacity/(4*sqrt(P)))
batch   = max(default, MaxActive+8)   if MaxActive+8 <= limit
        = default                     otherwise
batch   = min(batch, Capacity)
```

**Default.** Every pending function adds O(N) work, replicated on all ranks,
to each energy evaluation, and a full batch costs one fresh factorization, so
the best batch balances the two and grows in proportion to the order. With a
distributed factorization a flush is cheap and the best batch is small: for
the 2850-function problem of section 12 the cost per evaluation is flat for
batches from 16 to 64 and 20% higher at 178, the default of the time when
the factorization was serial (`Capacity/16`). With one rank the factorization
is the serial one and a larger batch would be slightly better, but a batch of
44 then costs only about 5% (section 12), so one formula serves both.

**Active sets and direct trials.** All functions that are active at the same
time must fit into the batch together with some room for the functions
optimized before them, hence `MaxActive+8`. A block above the limit (a large
`FULL_OPT1` range) does not get a batch of its own size. Such a trial is
applied *directly*: the columns are stored and the matrices are factorized
fresh. A trial that changes M functions through the batch costs M base
solves, which become cheaper with every rank, plus work on the batch that
grows as M^2*N and is replicated on every rank; a direct trial costs one
factorization, which the ranks share. The limit therefore falls as
1/sqrt(P): 712 functions for the 2850-function problem with one rank and 190
with 14 ranks, where the two paths were measured to cost the same near 210
(section 12). The limit also keeps the batch storage (four Capacity-by-batch
arrays on every rank) below that of the physical matrices. The same direct
path serves a separation command whose number of selected functions exceeds
the batch chosen for one active function (`EnsureQActiveCapacity` does not
change the batch).

## 6. Trial assembly and derivatives

`AssembleQTrial` stages complete normalized target columns before it modifies
either the physical matrices or the QR state. For each unordered matrix
pair touching at least one active function, it:

1. assigns the pair through the existing MPI work distribution;
2. calls the local basis-qualified H/S matrix-element routine;
3. computes all active raw diagonal overlaps needed for normalization;
4. places the normalized physical element into every affected staged column;
5. places derivatives in the optimizer block associated with the active
   endpoint; and
6. combines distributed results before a transaction begins.

When both endpoints are active, the shared matrix element is evaluated once
and copied into both staged columns, so the two columns agree at their
intersection. Derivatives with respect to both endpoint parameter blocks may
be needed: the matrix-element call must request `DActive` for the first
active endpoint and `DOther` only when the other endpoint is also active.
Never assume a trailing active block when indexing `Glob_D`.

Overlap penalties use the same active map. Each unordered pair touching an
active function is counted exactly once.

The matrix-element routine names and argument order are basis-specific. See
section 11 before copying this code to another variant.

## 7. Replacement transaction

`ApplyQTrial` changes an existing fixed-order basis atomically at the ECGPACK
level:

1. Validate method, order, active indices, staged columns, and allocations.
2. Verify that qrlinalg is valid and has the same order, capacity, and shift
   as the physical workspace.
3. Snapshot all original active columns, raw norms, and represented nonlinear
   parameters.
4. Announce the active set (`AnnounceQActiveFunctions`). On `need_flush` or
   `QR_ERR_SINGULAR` the current physical matrices are factorized fresh and
   the set is announced again.
5. For each active function call `replace_symmetric` with its complete staged
   column.
6. Only after every column has been accepted, store the staged columns and
   diagonal norms in `Glob_H`, `Glob_S`, and `Glob_diagS` on every rank. The
   solve that follows reads these matrices.
7. Mark the factor and physical generations as matched and update
   `MatrixParam`.

A one-function problem and an active set larger than the batch replace steps
4 to 6 by "store the columns, factorize fresh" (a direct trial).

If a qrlinalg call fails unexpectedly, restore all physical snapshots and
represented parameters, invalidate the staged trial, and factorize fresh.
Never continue from a state whose factor/matrix relationship cannot be proven.

DRMNG evaluates many points that are not accepted. The physical matrices and
the QR state follow the most recently evaluated point; with the batch this
costs nothing extra, because a new trial overwrites the tail entry of the
previous one. The optimizer's best point remains separate. When accepting
`x_best` or restoring `x_init`, reconstruct that complete point through the
same transaction; restoring parameters alone is invalid.

## 8. Growth and deletion transactions

### Basis growth

`BASIS_ENL Q` keeps the accepted prefix canonical and represents a candidate
block as a consecutive tentative suffix.

- From order zero, stage the complete first block, store it physically, and
  factorize fresh because qrlinalg has no valid order-zero state. The shift
  rules of section 1 apply.
- With a nonempty accepted prefix, `AppendQTrial` announces the suffix and
  appends the staged functions one at a time with `append_symmetric`; each
  physical column is stored after its append succeeded.
- Competing candidates at the same tentative order replace the previous
  suffix: `EvaluateQAppendedTrial` first calls `TrimQFactors`, which calls
  `truncate`.
- `Glob_CurrBasisSize` remains the accepted size while
  `Q_Workspace%MatrixOrder` may include the tentative suffix.
- When the batch is full, the accepted prefix (not the candidate) is
  factorized fresh. `truncate` reporting `QR_ERR_BATCH_FULL` is the same
  ordinary flush at the target order.
- `TrimQFactors` accepts a QR state whose order is anywhere from the target
  order to the physical order: a candidate that could not be assembled leaves
  the state at the accepted order.
- On an append failure the accepted leading block is exposed again and
  factorized fresh.

Discrete angular metadata is not part of `MatrixParam`. Every metadata-only
candidate change must therefore rebuild the candidate suffix with
`EvaluateQAppendedTrial`. After scanning choices, explicitly reconstruct the
winning metadata before continuous optimization.

### Basis deletion

The batched state has no interior deletion. `DeleteQMaskedFunctions` checks
the QR state against the workspace, compacts the lower-triangular physical
matrices and every basis-indexed array in survivor order, and factorizes the
survivors fresh. No matrix element is recalculated. A cleanup command deletes
once, so one O(n^3) factorization is the appropriate cost. Code that can
delete every function must handle the order-zero case outside qrlinalg.

### Separation

The separation commands perturb existing nonlinear parameters without
changing canonical indices. Their Q path stages and replaces the selected
columns through `ApplyQTrial`.

## 9. Solve, residuals, and flush policy

### What was measured

Three properties of a large, nearly linearly dependent basis shaped
`SolveQ`. They were measured on the 2850-function problem (extended
precision, epsilon = 1.1E-19, `EIGVAL_TOLERANCE` 1E-12, shift 4E-6 below the
eigenvalue) with instrumented scratch builds.

1. **The direction test fails in the residual-correction form.** The
   correction `A^{-1}(H v - mu S v)` carries roundoff noise along the
   near-null space of S. The physical residual reaches 2E-19 after three
   iterations, but the Euclidean change of direction between iterates stays
   near 5E-11 and never meets the tolerance; a call with the normal iteration
   limit then always runs to the limit (60 iterations and 450 ms per
   evaluation in the first integration). A direct iteration does not have
   this defect: `S v` annihilates such components of `v`, and the errors of
   the solve are small beside its solution, which is dominated by the wanted
   eigenvector. Method I and the previous Q implementation iterate directly
   and meet the test after 1.7 to 1.8 iterations.
2. **The plain Rayleigh quotient is too noisy for DRMNG.** Compared with a
   quadruple-precision quotient of the same vector, `(c^T H c)/(c^T S c)`
   was wrong by up to 15 units of the last place (typically 3 to 10), the
   shifted quotient with separately formed `H c` and `S c` by up to 10
   (typically 2 to 4), and the shifted quotient with elementwise
   `H - sigma*S` by at most one half. With the plain quotient DRMNG
   stopped early (16 105 evaluations for 151 functions where method I made
   19 925) and the energy fell 9E-14 behind. With the elementwise form the
   evaluation counts and energies follow method I (section 12).
3. **One direct iteration is usually enough for the residual, not for the
   rule of method I.** After one direct iteration from the previous
   eigenvector the residual was below 1E-12 in 99% of the solves and below
   1E-15 in 76%. A residual goal of 1E-12 would need 1.01 solves per
   evaluation. The direction test of method I needs 1.7.

### SolveQ

`SolveQ` calls `solve` with an iteration limit of one (a *pass*: one direct
inverse iteration plus the physical quotient and residual) and iterates
itself.

```text
tolerance = EIGVAL_TOLERANCE          direction test of one direct iteration;
                                      solve returns success exactly when met
goal      = max(|EIGVAL_TOLERANCE|, 100*epsilon), not above the threshold
threshold = max(1000*epsilon*n, 100*|EIGVAL_TOLERANCE|)
```

1. **Direct passes**, up to the normal iteration limit. They end when the
   direction test is met (the stopping rule of method I), when the residual
   is at roundoff level (10*epsilon), or when a pass no longer lowers the
   residual. The eigenpair is accepted if the direction test was met and the
   residual does not exceed the threshold, or if the residual meets the goal.
2. Otherwise **passes with `residual_goal`** follow, each one direct
   iteration plus residual-correction steps of qrlinalg, until the goal is
   met or the residual no longer decreases. The corrections recover the
   accuracy that block elimination over a large batch can cost a direct
   iteration. The eigenpair reached is accepted when its residual does not
   exceed the threshold. These solves are counted as "polished".
3. If the eigenpair is still not acceptable, or qrlinalg reported an error,
   the matrices are **factorized fresh** when anything is pending (a
   residual-triggered flush), and passes as in 2 are repeated with the larger
   limit `max(120, normal limit, 4*n)`. The result is final. The active set
   is not announced again here; these passes use the fresh factors alone and
   the next trial announces it.

A final eigenpair that is not acceptable is reported as
`QR_ERR_NO_CONVERGENCE`, or as the qrlinalg error when no vector was
returned. `Glob_c` keeps the last returned vector. An accepted solve ends
with `shifted_quotient`, whose result is the energy.

Cost: a pass is one solve with the QR state and three symmetric products.
The alternative of item 3 above (accept on a residual goal alone) would save
about 40% of the passes. It changes the meaning of `EIGVAL_TOLERANCE`, and
its effect on the optimizer has not been validated with the shifted quotient,
so it was not adopted.

### Physical generalized-eigenpair residual

```text
relative = ||H*c - E*S*c|| / (||H*c|| + |E|*||S*c|| + tiny)
```

`solve` returns it from the two products it forms for the quotient. It is
computed from H and S themselves and is independent of the QR state, so it
also detects a batch that has become too inaccurate. It is the only monitor
of the batch.

### Factorization probe

`ProbeQFactorization` maps a fixed dense vector `v` to `b = A v` with the
physical matrices, solves `A x = b` with the QR state, and stores the relative
residual of `x` in `Q_Workspace%LastFactorResidual`. `FactorizeQFresh` runs it
once on every fresh factorization of order above one, the only point at which
it measures the factors: with pending functions the same quantity measures
the error growth of block elimination, which is not a defect of the
representation. A result above `max(1000*epsilon*n, 10*EIGVAL_TOLERANCE)` is
reported as an invalid state; the factors then still represent the matrices
as well as a fresh factorization can, and the caller decides what to do with
the failed evaluation. The probe costs one solve and four products.

### When a fresh factorization happens

- at the start of every Q command, and for every direct trial;
- **batch full**: `begin_active` returns `need_flush`, or `truncate` returns
  `QR_ERR_BATCH_FULL`;
- **singular batch**: `begin_active` returns `QR_ERR_SINGULAR`;
- **residual**: stage 3 of `SolveQ`;
- **shift retarget**: `RetargetQStaleShift` in `BASIS_ENL` (section 1);
- deletion, and every recovery from a failed transaction.

There is no update-count limit and no periodic probe any more: the base
factors do not degrade, and the batch is bounded by its size.

The residual-triggered flush remains the last resort of `SolveQ` although a
flush has become cheap: on the 2850-function problem a pass costs about 9 ms
(one base solve and three products) and a distributed factorization 2.4 s,
so a few polishing passes are still far cheaper, and no run of section 12
needed a residual-triggered flush (none in about 50 000 evaluations, and at
most 3 polished solves per hundred functions).

`PrintQStatistics` prints one line at the end of `BASIS_ENL`, `OPT_CYCLE`,
and `FULL_OPT1` with Q:

```text
QR factorizations: 4 (batch full 3, singular batch 0, residual 0, direct trials 0, shift retargets 0); polished solves: 0
```

Frequent residual-triggered factorizations mean that the batch loses accuracy
for this basis and shift; many polished solves mean that direct iterations
alone do not reach the goal (slow convergence, or an `EIGVAL_TOLERANCE` below
what roundoff permits). The same happens when the energy has moved far from
the shift within one `OPT_CYCLE` or `FULL_OPT1` command: direct passes then
stagnate at a residual of about epsilon*|E-shift|/|E0-shift|, where E0 is the
energy the shift was derived from (1E-12 for a displacement of 0.05 hartree
at `INVITPARAMETER` 1.000001 in double precision), and almost every solve is
polished. The energies stay correct.

## 10. BBOP driver behavior

The drivers did not change with the redesign beyond the statistics line and,
in `BasisEnlQ`, the call of `RetargetQStaleShift` and the handling of a block
that every attempt rejected (below).

### OPT_CYCLE Q

The driver preserves G's validation, scheduling, DRMNG settings, counters,
acceptance rules, history, and save frequency, but it never physically moves
the optimized functions. The active map reproduces the optimizer-variable
order expected by the legacy G schedule. Every optimized function stays in
the batch after its optimization, so the batch fills at one entry per
function and a fresh factorization follows every `max_pending` functions.

### FULL_OPT1 Q

The active map is the requested `InitFunc:FinalFunc` range, and Hessian I/O
uses that same variable order. A trial changing M functions costs M base
solves plus O(M^2*n) work on the batch, or one fresh factorization when M
exceeds the batch (section 5).

### BASIS_ENL Q

Candidates are appended as tail entries (section 3.2); the accepted function
stays in the batch. Shift retargeting and the start from an empty basis are
described in section 1.

A block rejected by the overlap or by the linear-coefficient threshold is
generated again, up to `Glob_BadOverlapOrLinCoeffLim` = 10 attempts. When all
attempts are rejected, `BasisEnlQ` keeps the last block and continues, with
the message `No acceptable candidate block after 10 attempts; the last block
is kept`, exactly as `BasisEnlG` and `BasisEnlI` do; like them it records the
energy before the block as the current energy. (The previous `BasisEnlQ`
stopped with `EC0130` here. There was no Q-specific reason for it: the QR
state, the physical matrices, and the eigenvector describe the kept block
consistently.) `EC0130` remains for the one case that is Q-specific: the
final energy evaluation of the last attempt failed. After a failed Q
transaction the physical matrices may represent other nonlinear parameters
than `Glob_NonlinParam` holds, and `Glob_c` is not an eigenvector, so the
block cannot be recorded. `RetargetQStaleShift` is called with the energy of
the represented problem (`Evalue`), not with the recorded one.

The codes with discrete metadata scan it before the continuous optimization
(section 11), and the PG codes add a shape test to the acceptance tests: a
block in which the Young operator nearly annihilates a new function is drawn
again without limit and outside the attempt budget, as in their `BasisEnlG`
and `BasisEnlI`. A block that is kept after ten rejected attempts has
therefore always passed the shape test.

### Output and observable commands

`SAVE_HSWF Q` exports canonical, normalized, unshifted physical H and S. It
does not export `H-sigma*S` or private QR storage. Expectation-value and
density commands build canonical Q state, factorize, solve once, and then use
the existing operator paths with the S-normalized coefficient vector.

### Swap state

Swap files persist physical matrix information and enough metadata to restore
canonical H/S. The QR state is never treated as portable persistent state; it
is reconstructed after a restore. Any swap format change must be tested for G
and I as well as Q because the routines are shared.

## 11. Basis-specific adaptations

The canonical layout and the helper layer are common, but matrix-element
interfaces and discrete metadata are not:

| Code | Metadata for function `i` | H/S matrix-element routine |
| --- | --- | --- |
| `RG_0S` | none | `MatrixElementsHS_RG_0S` |
| `RG_1P` | `Glob_ZIndex(i)` | `MatrixElementsHS_RG_1P` |
| `RG_2D` | `Glob_Index(i,1:2)` | `MatrixElementsHS_RG_2D` |
| `RG_2P` | `Glob_Index(i,1:2)` | `MatrixElementsHS_RG_2P` |
| `PG_0S` | `Glob_PWR(i)` | `MatrixElementsOpt` |
| `PG_1P` | `Glob_PWR(i)`, `Glob_X_index(i)` | `MatrixElements` |
| `CG_0S` | none | `MatrixElementsHS_CG_0S` |

Consult the actual routine declaration before adapting a call. In particular,
RG_2D and RG_2P use different argument ordering, and RG_2D returns additional
`Tkl` and `Vkl` values. The PG routines take the premultiplier power (and, in
`PG_1P`, the index of the angular prefactor) of both functions, and `PG_0S`
passes precomputed matrices of the active function. The active function must
be passed as the derivative endpoint associated with `DActive`.

Discrete metadata must be generated, broadcast, copied during acceptance,
saved, restored, printed, optimized when supported, and compacted during
deletion:

- RG_1P carries one Z index;
- RG_2D carries two indices and may restrict their relationship with
  `VECTOR_COUPLING_SCHEME`;
- RG_2P carries two distinct indices;
- PG_0S carries the even power of the premultiplier `r_1^m`; and
- PG_1P carries that power and the index of the angular prefactor `z_i`
  (fixed by `FIXED_INDEX` when that keyword is given).

### Discrete scans in BasisEnlQ

With `OptimizationType` 1 the five codes with metadata scan it for the
selected block before DRMNG starts: `RG_1P` tries every other Z index,
`RG_2D` and `RG_2P` every allowed value of each index, `PG_0S` every even
power up to `Glob_MaxPowerAllowed` = 250 (125 alternatives per function), and
`PG_1P` every combination of prefactor index and power. Each alternative is
one `EvaluateQAppendedTrial` (section 8), and the winner is reconstructed
before the continuous optimization. The scans interact with the rest of
`BasisEnlQ` as follows.

- They use `Evalue` as a scratch variable. The `EnergyQAM` call after the
  optimization sets it again, and nothing after that call changes it, so the
  energy passed to `RetargetQStaleShift` is that of the represented problem
  in every code.
- Their solves are counted by `Glob_InvItTempCounter1` and
  `Glob_InvItTempCounter2` together with the candidate and optimizer solves
  of the block. They are solves of the same kind with the same shift, so they
  neither dilute nor inflate the average on which the stale-shift rule
  decides. A solve that fails counts once in the first counter and with all
  its passes (up to the limits of the three stages of `SolveQ`) in the second,
  for a scan alternative exactly as for a random candidate; a block with
  several failed solves can therefore trigger a retarget. That costs one
  fresh factorization and moves the shift to where `main` would put it for
  the next command.
- A block kept after ten rejected attempts keeps the metadata that its last
  scan selected; the physical matrices and the QR state were built with it.

Measured with the sample inputs (double precision, 1 and 4 ranks; averages
per function as printed by `BasisEnlQ`): in bases of more than 30 functions
the average was 2.1 to 2.2 passes per solve in `RG_1P` and `RG_2D` (9000
functions each, none above 4), against 2.0 in `RG_0S`; in the 10- to
20-function PG cases 1.5 to 1.7. Values above 4 occurred only while a basis
grew from a few functions, which is the case the rule is for.

### RG_2D vector-coupling schemes

```text
scheme 0: no relationship restriction
scheme 1: the two indices must differ
scheme 2: the two indices must be equal
```

`GenerateTrialParam` establishes the invariant for random candidates.
`BasisEnlQ` must preserve it during discrete optimization: scheme 1 excludes
the other index from each one-index scan; scheme 2 varies both indices
together as one discrete choice; and if either scheme-2 index is fixed, both
are effectively fixed.

### CG_0S: Hermitian matrices

In `CG_0S` the basis functions are complex, so H and S are Hermitian and the
linear coefficients are complex. The workproc conventions are:

- The authoritative lower triangle stores `<phi_i|O|phi_j>` for `i>=j`, with
  function `i` (the row index) being the bra. The conceptual element `(i,j)`
  with `i<j` is the complex conjugate of the stored element `(j,i)`. Every
  staged column is a conceptual column: element `i` of the column of function
  `k` is `<phi_i|O|phi_k>`.
- `MatrixElementsHS_CG_0S` is called with the active function as the bra, so
  `AssembleQTrial` stores the complex conjugate of the computed element in the
  column of the active function. Diagonal elements are real.
- In `Glob_D` the differentiated function is always the bra. The diagonal
  derivative is twice the real part of the bra derivative.
- The energy gradient is `2*Re[conjg(c_k)*W] - |c_k|^2*(...)`, the overlap
  penalty depends on `|S_ij|^2`, and the S-normalization is `c^H S c = 1`.

`qrlinalg.f90` contains the complex twin of the distributed, batched state
(`qr_complex_state`, same procedures and the conventions above: lower
triangle with the row function as the bra, real diagonals, conceptual
columns). `shifted_quotient` of the complex state returns the real quotient
`sigma + c^H (H - sigma*S) c / (c^H S c)` from the stored lower triangles.

A conjugation error in the `CG_0S` helper layer does not show in the
energies: `SolveQ` detects the damaged representation through the physical
residual and repairs it with fresh factorizations or polishing passes. It
shows in the statistics line (residual-triggered factorizations, polished
solves) and in a direct probe of the QR state while columns are pending
(`factorization_residual` with a complex probe vector: 1E-12 or below when
correct, of order 0.1 to 1 with a wrong conjugation). Validate `CG_0S` on
bases with genuinely complex parameters; the He sample basis is real to
1E-7.

## 12. Memory and performance model

For matrix order N, capacity C, batch size B, active block M, and P ranks:

- physical H and S: 2*C^2 values on every rank, as in the other methods;
- base factors: 2*C^2/P values per rank, plus P*C for the diagonal blocks;
- batch: four C-by-B arrays and matrices of order 2*B on every rank;
- staged current/previous active columns: O(M*C);
- during a fresh factorization: about 2*(1+1/P)*32*N values per rank with
  the distributed algorithm; N^2 values on rank zero with the serial one.

Operation counts:

- fresh factorization: O(N^3), shared by the ranks (distributed algorithm)
  or on rank zero (serial algorithm);
- base solve: O(N^2/P) per rank plus N/P small reductions;
- symmetric product: O(N^2/P) per rank plus one reduction;
- replacing or appending one function: one base solve and O(B*N) work that is
  replicated on every rank;
- one pass of `SolveQ`: one base solve, O(B*N) replicated work, three
  products;
- `begin_active` for a new base function: one base solve, one column
  recovery, O(B^2*N/P) per rank.

### Measured: 2850 functions, extended precision, 14 ranks

Machine: one exclusive node of muon (2 x 7 cores of Xeon E5-2690 v4),
gfortran 14.3, Open MPI 5.0.8, `-O3 -march=native`, bundled BLAS/LAPACK.
Input: four particles, `OPT_CYCLE` over a 2850-function basis, one function
at a time, at most 100 energy evaluations per function, `EIGVAL_TOLERANCE`
1E-12, default batch (44). An "evaluation" is one energy or
energy-plus-gradient evaluation including its matrix elements. The I and Q
columns were measured one after the other on the same node and cover the
same 108 functions; the previous implementation was measured earlier on the
same input.

| Quantity | Method I | Q, this design | Q, previous implementation |
| --- | ---: | ---: | ---: |
| Start-up (matrix elements, factorization, first solve) | 3.1 s | 3.6 s | not recorded |
| Time per evaluation, first 20 functions | 19.2 ms | 22.8 ms | 527 ms |
| Time per evaluation, average over the 108 functions, flushes included | 19.4 ms | 25.0 ms | 527 ms |
| Direct inverse iterations (passes) per solve | 1.78 | 1.56 | 1.7 |
| Evaluations made for the 108 functions | 13 937 | 13 995 | |
| Seconds per optimized function, flushes included | 2.50 | 3.24 | 75 |
| Fresh factorizations after the first | | 2 (batch full) | |
| Cost of a fresh factorization (distributed) | | 2.45 s | |
| Residual-triggered factorizations, polished solves | | 0, 3 | |

The energies after the same functions agree to 4.1E-14 hartree (3E-15
relative) or better over the 108 functions; two runs of method I from
different builds differ by as much. A fit of the Q column gives 22.0 ms plus
0.15 ms per pending function (0.06 to 0.13 ms in three other runs of 75
functions).

**Batch size.** Runs of about 75 functions with other batch sizes, flushes
included: 24.97 ms per evaluation with a batch of 16 (4 batch-full flushes),
24.58 ms with 32 (2 flushes), 25.63 ms with 64 (1 flush); with the former
default of 178 and the serial factorization, 29.5 ms (no flush within 98
functions, and 41 s to start). A function during which a flush happens takes
2.4 to 2.7 s longer. No residual-triggered flush occurred in any run.

**One rank.** The same input with one rank (on a node shared with other
jobs): 49 s to start (matrix elements, and 41.2 s for the serial
factorization), 133 ms per evaluation and 17.5 s per function over 55
functions including one batch-full flush of 41.2 s (122 ms per evaluation
over the first 20 functions), 1.55 passes per solve, energies within 4.9E-14
hartree of the 14-rank method-I run. With one process
the automatic choice is always the serial algorithm; the distributed one
was measured to be faster there as well (order 2850 on the 2.4 GHz login
node: 38.0 s against 47.0 s in extended precision, 6.2 s against 10.6 s in
double precision) and can be selected with `ECG_QR_FACTOR=distributed`.

**Double precision.** The same input (a basis optimized in extended
precision) run with the precision-8 binary and 14 ranks: 12.7 ms per
evaluation with Q over 101 functions (11.8 ms over the 68 functions that the
I run also covered) against 11.4 ms with I, start-up 1.7 s against 2.4 s,
energies within 5.6E-13 hartree.

**Large active blocks** (`FULL_OPT1` over the first M functions of the same
basis; seconds per trial for the transaction and the solve, without the
matrix elements, which take 0.15 to 0.55 s for an energy and 0.6 to 2.2 s
for a gradient with 14 ranks):

| M | 14 ranks, batch | 14 ranks, direct | 1 rank, batch | 1 rank, direct |
| ---: | ---: | ---: | ---: | ---: |
| 178 | 1.95 | 2.50 | | |
| 356 | 6.2 | 2.52 | 10.1 | 48.0 |
| 700 | 24.1 | 2.57 | 25.0 | 46.5 |

The limit of section 5 sends M above 182 to the direct path with 14 ranks
and M above 704 with one rank. The one-rank columns were measured on the
login node (Xeon E5-2640 v4) with two such runs at a time.

**`BASIS_ENL`** on the same basis (functions 2851 to 2854, 500 candidates
per function, fixed seed): 39 s per function with Q (37 s with I, measured
before the distributed factorization), 1.90 to 1.92 passes per solve for
every function, no shift retarget, no flush.

Components (same machine, synthetic matrices of the same order, from the
module's timing driver): fresh factorization 2.44 s distributed and 40.5 s
serial; base solve 4.8 ms; symmetric product 1.3 ms; `begin_active` for a
new base function 6 to 7 ms.

Observations that matter for tuning:

- The time per evaluation grows by roughly 0.1 ms per pending function. The
  O(B*N) work is replicated on all ranks and the node is bound by memory
  bandwidth. Distributing that work is the most promising optimization of
  the module, and it would also move the limit for direct trials upwards.
- The base solve is dominated by its N/P small reductions, not by arithmetic.
- A pass costs one solve and three products. Accepting a solve as soon as
  its residual meets a goal of 1E-12 would need about one pass per solve
  (section 9).
- Method I remains 1.3 times faster per evaluation on this input (it was 1.5
  times faster before the factorization was distributed and the batch
  reduced).

### Measured: start from an empty basis

Double precision, He (three particles) and Li (four particles) sample
headers, `RANDOM_SEED` fixed, 500 candidates and up to 200 evaluations per
function, thresholds 0.95 and 3.0.

Passes per solve by function, He sample command list with every solver
letter switched to Q (one rank):

| Functions | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| shift fixed after the first function (before) | 6.3 | 10.6 | 14.7 | 15.5 | 3.0 | 4.6 | 4.6 | 4.9 | 4.9 |
| retargeted when stale | 6.2 | 4.7 | 4.3 | 3.9 | 3.0 | 4.6 | 2.9 | 4.1 | 2.7 |

(Functions 6 to 10 belong to a second `BASIS_ENL` command, whose shift
`main` derives from the energy of five functions.) The all-Q run reproduces
the all-G run with the same seed to the last digit with one rank
(-2.900872393343730 against -2.900872393343731); with four ranks the two
differ by 1.6E-8, as two G runs with one and four ranks do (3E-9). Six
shift retargets, and 26 polished solves in the first command where the fixed
shift needed 1202.

One command from the empty basis to 40 functions (`BASIS_ENL Q 1 40`): 14 to
16 retargets, 2.1 to 6.5 passes per solve. He: G and Q end at
-2.9023309706 and -2.9023309702 (one rank); 21 of the 40 blocks are kept
after ten rejected attempts, in G as in Q. Li: -7.4735622454 (G) and
-7.4735622859 (Q) with one rank, -7.4735623058 and -7.4735622783 with four.
In small bases a function lowers the energy by a sizable fraction of the
distance to the next eigenvalue, so 3 to 5 passes per solve are normal there
and the rule retargets after about every second function, at the cost of a
tiny factorization.

`FULL_OPT1` becomes cubic per trial when M exceeds the limit for direct
trials. Matrix element generation, MPI communication, and optimizer behavior
can dominate the observed wall time of small problems.

OpenMP is process-local. Every MPI rank may create `OMP_NUM_THREADS` threads.
Only the serial fresh factorization on rank zero (through the `DGEMM` calls
of `DGEQRF` and `DORGQR` in the bundled BLAS) uses them; the distributed
factorization, solves, products, and batch updates are parallelized over MPI
ranks only. The calibration of the automatic choice times both algorithms as
they run, threads included.

## 13. Build and third-party integration

`qrlinalg.f90` depends only on `wp_def`, on MPI, and on four
LAPACK routines of the bundled (or external) provider: `DGEQRF` and `DORGQR`
for the serial factorization of the real state, `ZGEQRF` and `ZUNGQR` for
that of the complex state. All four are in `src/LAPACK.f` of every energy
code, and the distributed factorization calls no library routine.

```text
wp_def -> qrlinalg -> workproc
```

It requires preprocessing: `-cpp` for gfortran, `-fpp` for ifort/ifx, or
`-Mpreprocess` for nvfortran (the only preprocessor symbol is
`QRLINALG_TESTING`, which must stay undefined). The fixed-form bundled
BLAS/LAPACK sources also require the compiler-specific unlimited-line-length
option. The module uses the Fortran 2003 intrinsic
`get_environment_variable` (for `ECG_QR_FACTOR`) and `MPI_Wtime` (for the
calibration of the automatic choice); it contains no OpenMP directive.

### Maintaining the module

The real part of `qrlinalg.f90` is the source of truth. The complex twin is
generated from it by `gen_complex.py` (with `gen_complex_text.py`), which
applies procedure-specific substitutions that must each match an exact
number of times, then generic renaming, and writes the full module as the
real text plus pure insertions. A change of a real procedure therefore either
translates mechanically or stops the generator. The generator, the standalone
MPI unit test (`test_qrlinalg.f90`), and a script that builds and runs it are
in `doc/devnotes/qrlinalg_generator/`; its `README.md` is the procedure. The
real source is not kept as a separate file: the generator extracts it from
the installed module (`--regenerate`, `--extract`) and can verify that a
module is what it produces from its own real part (`--check`).

The rule: edit the real part of `RG_0S/src/qrlinalg.f90`, regenerate, run the
unit test, and copy the file to the other six codes. Do not edit the complex
part by hand, and keep the file byte-identical in all seven energy codes.

No energy code calls `src/qrupdate/` any more. In all seven codes its sources
and Makefile rules are still present, `qrlinalg.o` no longer depends on them,
and its objects are still compiled and linked:

```text
wp_def -> qrupdate_linalg -> qrupdate_error -> qrupdate_real, qrupdate_complex
       -> qrupdate                              (compiled and linked, unused)
```

Whether the bundled library is removed from the codes (sources, Makefile
rules, the `qr1up_unit` addition, and its section of
`THIRD-PARTY-NOTICES.md`) is an open decision for the maintainer.

Only one BLAS/LAPACK provider may be linked. With `LINALG=netlib`, the build
uses `src/BLAS.f` and `src/LAPACK.f`. Optimized precision-8 providers resolve
the calls externally. Precisions 10 and 16 require the bundled provider.

The bundled QR-update kernels come from Martin Köhler's `qrupdate-ng`:

<https://gitlab.mpi-magdeburg.mpg.de/koehlerm/qrupdate-ng>

The imported revision is `3fa4f77d6259c00cf0c36632287063f1277a3fe0`, with one
ECGPACK addition (`qr1up_unit`) that is marked in the sources, in
`src/qrupdate/README.md`, and in `THIRD-PARTY-NOTICES.md`. Preserve the source
headers, per-directory `COPYING` files, and notices for as long as the code is
bundled.

## 14. Source map

| File | Routines or responsibility |
| --- | --- |
| `src/matform.f90` | `StoreHS` selects canonical normalized Q storage. |
| `src/workproc.f90` | `QBatchSize`, `PrepareQWorkspace`, `EnsureQActiveCapacity`, `ClearQWorkspace`. |
| `src/workproc.f90` | `SetQActiveFunctions`, `CaptureQMatrixParameters`, `AnnounceQActiveFunctions`. |
| `src/workproc.f90` | `QCanonicalMatrixElement`, `GatherQCanonicalColumn`, `StoreQCanonicalColumn`. |
| `src/workproc.f90` | `AssembleQTrial`, `ApplyQTrial`, `AppendQTrial`, `EvaluateQAppendedTrial`, `TrimQFactors`. |
| `src/workproc.f90` | `FactorizeQFresh`, `ProbeQFactorization`, `SolveQ`, `RetargetQStaleShift`, `PrintQStatistics`. |
| `src/workproc.f90` | `ComputeQOverlapPenalty`, `ComputeQOverlapPenaltyAndAddGradient`, `GetQOverlapStatistics`. |
| `src/workproc.f90` | `EnergyQA`, `EnergyQAM`, `EnergyQB`. |
| `src/workproc.f90` | `BasisEnlQ`, `OptCycleQ`, `FullOpt1Q`. |
| `src/workproc.f90` | `DeleteQMaskedFunctions`, `SolveEliminationGSEP`, and Q branches in cleanup/output routines. |
| `src/main.f90` | BBOP dispatch for method Q. |
| `src/qrlinalg.f90` | Real and complex QR states: serial and distributed base factorization, batch, bordered solves, inverse iteration, shifted quotient. |
| `src/BLAS.f`, `src/LAPACK.f` | Bundled Netlib provider (`DGEQRF`, `DORGQR`, `ZGEQRF`, `ZUNGQR` for the serial fresh factorization). |
| `Makefile` | Dependency order, preprocessing, provider selection, OpenMP build isolation. |
| `doc/devnotes/qrlinalg_generator/` | Generator of the complex half of `qrlinalg.f90`, unit test of the module, and their instructions. |

Line numbers are intentionally omitted. Search by routine name because the
duplicated workproc files evolve independently. The PG codes write the same
names in their own source style (upper-case keywords).

## 15. Validation requirements

Numerical changes to Q are incomplete until they pass tests that exercise both
the physical eigenproblem and the batch lifecycle. Use disposable copies of
sample inputs; energy-code runs rewrite `inout.txt`.

### Module

The standalone MPI unit test of `qrlinalg.f90` covers the real and the
complex state: random symmetric and Hermitian problems, sequences of
replace/append/truncate/commit with interleaved solves, every solve compared
with a dense reference; stress cases with a full batch; direct tests of both
factorization algorithms (backward error, orthogonality, solves); the
shifted quotient against a quotient evaluated by the driver in a wider real
kind, within the rounding-error bound of the working precision; error paths
and the conventions of the complex state. It must pass at 1, 2, 3, 5, and 14
processes in precisions 8 and 10, at one process in precision 16, and in
strict debug builds, each with the serial, the distributed, and the
automatic factorization, with identical results on all ranks and agreement
between process counts within the stated tolerances. The driver and the
script that builds and runs it are in `doc/devnotes/qrlinalg_generator/`.

### Required build coverage

- strict debug build at precision 8;
- release build at precisions 8, 10, and 16;
- bundled `LINALG=netlib` path; and
- at least one optimized precision-8 BLAS/LAPACK provider when available.

### Required runtime coverage

- order one, order two, active order smaller than capacity, and growth to
  capacity, including a start from an empty basis;
- active first, middle, last, contiguous, and noncontiguous functions;
- two active functions sharing a matrix intersection;
- an active block larger than the batch (direct trials);
- enough optimized or appended functions to fill the batch (flush);
- repeated rejected trials followed by reconstruction of `x_best` and
  restoration of `x_init`;
- G/Q and I/Q energy and expectation-value agreement on the same fixed basis;
- the self-check: the energy that a Q optimization reports must be what
  `EXPC_VALS G` computes for the basis it wrote;
- all four cleanup commands, `DENSITIES`, `MOMT_DENS`, `SAVE_HSWF`;
- one and several MPI ranks, with equal results to roundoff;
- precision 16 with one rank; and
- a representative large basis for timing, compared with method I on the
  same input: time per evaluation, iterations per solve, evaluations per
  function, energies after the same functions, and the statistics line.

For indexed bases, also verify metadata after generation, optimization,
deletion, separation, save, and restart. RG_2D tests must cover coupling
schemes 1 and 2; RG_2P tests must confirm distinct indices.

Stochastic energies are not reference values. Compare methods on the exact
same saved physical basis and report tolerances.

## 16. Porting procedure

These are the edits by which the design was carried from `RG_0S` to the
other six energy codes. They are kept for two purposes: as the record of what
a ported code must contain, and as the procedure for a later change. A change
of the Q layer is made in `RG_0S` first and then applied to the siblings unit
by unit, as far as the steps below name the unit. "Verbatim" means the
`RG_0S` text of the named unit, comments included.

Source style. `RG_1P`, `RG_2D`, and `RG_2P` are written like `RG_0S`, and
their verbatim units are identical to it line for line, so a change can be
applied to them as a patch
(`git diff <old> <new> -- RG_0S/src/workproc.f90 | sed 's#RG_0S/#RG_1P/#g' | patch -p1`).
`PG_0S` and `PG_1P` are written in another style (upper-case keywords, `::`
in every declaration, blanks around operators and after commas, comments
indented with the code, `!====` banners). There every unit taken from
`RG_0S` is rewritten in the style of the file, and nothing else in it
changes: after undoing the style, the routine must be identical to the
`RG_0S` routine. `CG_0S` is written like `RG_0S` and differs by the Hermitian
rules of section 11.

1. **Module.** Copy `RG_0S/src/qrlinalg.f90` unchanged. The file is
   byte-identical in every ported code: the real codes use `qr_real_state`,
   `CG_0S` uses `qr_complex_state`. Nothing else is needed for it: the four
   LAPACK routines it references (`DGEQRF`, `DORGQR`, `ZGEQRF`, `ZUNGQR`) are
   in the bundled `src/LAPACK.f` of all seven codes, and every Makefile
   already compiles `qrlinalg.f90` with `$(QR_PREPROCESS_FLAGS)`.
2. **Makefile.** Make `qrlinalg.o` depend on `wp_def_$(PREC).o` instead of
   `qrupdate.o`. Leave the `qrupdate` objects and link line alone (section
   13). The PG Makefiles describe the dependencies in a comment headed
   "Module compile order"; change its `qrupdate -> qrlinalg` chain to
   `wp_def -> qrlinalg -> workproc` there as well.
3. **Module head of `workproc.f90`.** Replace by the `RG_0S` text:
   - the `use qrlinalg` list (it gains `QR_ERR_SINGULAR` and
     `QR_ERR_BATCH_FULL`; `CG_0S` keeps `qr_complex_state` in place of
     `qr_real_state`);
   - the comment block on the Q method that follows `implicit none` and
     precedes the `Q_METHOD_*` constants (the comment above the `use`
     statements is not affected; in the PG codes this block is the banner
     "State of the QR solution method ('Q')", which follows the module
     variables of method I);
   - the `Q_METHOD_*` constants, followed by the new constant
     `Q_STALE_SHIFT_PASSES` with its comment; and
   - the `QMethodWorkspace` type.

   Remove the `Q_FACTOR_PROBE_INTERVAL` parameter. Remove the
   `iso_fortran_env` import only if nothing else in the file uses it:
   `PG_1P` keeps `use iso_fortran_env, only: int64`, which another routine
   needs.
4. **Helper routines, verbatim:** `QBatchSize` (new), `PrepareQWorkspace`,
   `ClearQWorkspace`, `TrimQFactors`, `FactorizeQFresh`,
   `AnnounceQActiveFunctions` (new), `AppendQTrial`, `ApplyQTrial`,
   `ProbeQFactorization`, `SolveQ`, `RetargetQStaleShift` (new),
   `PrintQStatistics` (new), and `DeleteQMaskedFunctions`. They contain
   nothing basis-specific. Place the new routines as in `RG_0S`: `QBatchSize`
   directly before `PrepareQWorkspace`, `AnnounceQActiveFunctions` between
   `FactorizeQFresh` and `AppendQTrial`, and `RetargetQStaleShift` followed
   by `PrintQStatistics` directly after `SolveQ`. In `CG_0S` apply the
   Hermitian rules of section 11 to these routines (complex columns, vectors,
   and coefficient array; real norms, parameters, energies, and residuals;
   the one-function closed form uses the real parts of the diagonals).
5. **Comment-only changes:** the header comments of `EnsureQActiveCapacity`,
   `SetQActiveFunctions`, `GatherQCanonicalColumn`, `StoreQCanonicalColumn`,
   `EvaluateQAppendedTrial`, and `SolveEliminationGSEP`, and the comment in
   the Q branch of `EliminateLittleContribFunc` above the call of
   `DeleteQMaskedFunctions` ("factorize the compacted matrices", no longer
   "delete the matching factor rows and columns"). The code of these
   routines is unchanged.
6. **`OptCycleQ` and `FullOpt1Q`:** add `call PrintQStatistics()` immediately
   before `call ClearQWorkspace()`. No other change.
7. **`BasisEnlQ`:**
   1. add `call PrintQStatistics()` immediately before
      `call ClearQWorkspace()`;
   2. replace the block after the acceptance loop that stops with `EC0130`
      whenever `IsOverlapBad.or.IsAnyLinCoeffBad.or.IsEnergyBad` by the
      `RG_0S` text: stop with `EC0130` only `if (IsEnergyBad)`, otherwise
      print the "No acceptable candidate block after ... attempts; the last
      block is kept" line when a threshold flag is still set, and continue.
      The other `EC0130` messages of `RG_2D` and `RG_2P` ("optimized indices
      cannot be reconstructed") and the additional tests of the PG codes
      (shape test, which redraws without limit) are not part of this edit;
   3. directly after that block, before the block is reported and recorded,
      add

      ```fortran
      call RetargetQStaleShift(Evalue,ErrCode)
      ```

      with an abort on failure: `EC0131` in the RG and CG codes, `EC0132` in
      the PG codes, whose `BasisEnlQ` already uses `EC0131` for an
      unsupported `OptimizationType`. `Evalue` must be the value returned by
      the last successful `EnergyQAM` of the acceptance loop, that is, the
      energy of the represented problem; check in each code that nothing
      between that call and this point overwrites it (the discrete scans use
      `Evalue` as scratch, but they run before that call; the shape test of
      the PG codes does not touch it). `Glob_InvItTempCounter1` and
      `Glob_InvItTempCounter2` must be reset once per block before the
      candidate loop (they already are, for the printed average) and must
      not be reset between the solves of the block and this call. The solves
      of the discrete scans stay in the counters (section 11);
   4. take the header comment of the `RG_0S` routine (attempt limit, shift
      retargeting). In the PG codes merge it into the existing banner
      comment of the routine, and correct the comment of the acceptance loop
      there, which described the stop with `EC0130`.
8. **Elimination commands:** update the two "Deleting selected rows and
   columns" messages (`EliminateLittleContribFunc`, `EliminateLinDepFunc`)
   as in `RG_0S`.
9. Leave untouched everything basis-specific: `AssembleQTrial` and the
   overlap helpers (matrix-element interface, metadata), the discrete index
   scans of `BasisEnlQ`, metadata compaction in the cleanup commands, the
   PG-only commands, and `EnergyQA/QAM/QB`. `DeleteQMaskedFunctions`
   compacts `Glob_H`, `Glob_S`, `Glob_diagS`, and `Glob_c` only; metadata
   compaction stays in the callers, as before.
10. **Check the result** with `grep`: no remaining reference to
    `RetargetQProvisionalShift`, `Q_FACTOR_PROBE_INTERVAL`, rank guards
    around `Q_Workspace%Factors%` calls, or broadcasts of their results;
    every `Q_Workspace%Factors%` call inside a helper routine of step 4 or 5;
    `diff` of every verbatim routine and of the module head against `RG_0S`
    empty for the RG codes, empty after the style is undone for the PG codes,
    and reduced to the Hermitian changes for `CG_0S`. The error codes of a
    code are not unique across its routines (`EC0125` to `EC0130` occur in
    more than one), so check a new code within its routine.
11. **Build and validate:** `make -j` debug and release, precisions 8, 10,
    and 16 (and `OPENMP=yes` where supported); then the fixed-basis G/I/Q
    comparison, the optimizer self-checks, the cleanup commands, and a
    `BASIS_ENL Q` start from an empty basis against `G` with the same
    `RANDOM_SEED`, at 1 and 4 ranks.
12. **User documentation:** `doc/input_file_format.md` and
    `doc/compilation_and_execution.md` describe `ECG_QR_FACTOR`, the batch,
    the shift retargeting, and the attempt limit for all seven energy codes
    alike. A change of behavior goes there for all of them at once, once all
    of them carry it.

Do not copy an entire `workproc.f90` over a sibling. Use targeted
replacements so basis-specific physics, file formats, and newer local
semantics remain intact.

What a porter must know about behavior:

- The automatic choice of the factorization algorithm is a timing made once
  per run (section 3.5). Results obtained with the two algorithms differ by
  roundoff, so two runs of the same input can differ in the last digits
  unless `ECG_QR_FACTOR` is set.
- The batch is now small (section 5), so batch-full flushes are frequent and
  cheap. A code whose trial evaluation is cheap relative to a flush (few
  particles, large basis, one rank) is the unfavorable case; measure it.
- `RetargetQStaleShift` changes `Glob_ApproxEnergy` inside `BASIS_ENL Q`. In
  every code confirm that, in Q mode, nothing reads `Glob_ApproxEnergy`
  except the Q helpers and `main` (which resets it before each BBOP step).

## 17. Production invariant checklist

Before approving a Q change, confirm all of the following:

- physical basis order is canonical and stable;
- lower H/S triangles are normalized and unshifted;
- upper H/S triangles are never treated as authoritative;
- every qrlinalg call is made on all ranks with replicated arguments, without
  rank guards and without broadcasts of its results;
- qrlinalg order, capacity, and shift match the physical workspace before a
  transaction;
- the active set is announced before `replace_symmetric` or
  `append_symmetric`, and `need_flush`, `QR_ERR_SINGULAR`, and
  `QR_ERR_BATCH_FULL` are handled by a fresh factorization;
- staged trials are complete and mutually consistent before the first
  qrlinalg call, and the physical columns are stored before the solve;
- failed or uncertain transactions recover through proven physical state and
  a fresh factorization;
- convergence is decided by direct passes and the physical residual, and the
  energy is the shifted quotient;
- the shift changes only together with a fresh factorization;
- saved H/S and wave functions remain in canonical user-visible order;
- basis-specific metadata follows every accepted, deleted, or saved function;
  and
- performance claims state the batch size, the rank count, and the number of
  pending functions.

If any relationship cannot be demonstrated from current state and flags,
factorize fresh. A fresh O(N^3) factorization is preferable to reusing an
unproven state.
