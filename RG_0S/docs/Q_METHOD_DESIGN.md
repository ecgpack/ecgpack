# Q-method implementation design and porting log

## Status and scope

This document is the implementation record for adding the `Q` generalized
eigenvalue method to `RG_0S`. It is also intended to be a recipe for porting
the same method to the other ECGPACK variants. Update it after every reviewable
chunk with the exact files changed, validation performed, open questions, and
any decision that a later port must reproduce.

The first three chunks define code ownership and canonical layout, integrate
qrlinalg and its build dependencies, and implement non-derivative canonical Q
storage plus arbitrary-index staged column assembly. Q remains absent from
`main.f90` dispatch, so no Q BBOP is usable yet and existing G/I input behavior
is unchanged.

The G implementations are the behavioral references:

- `BasisEnlG` for random candidates, optimization, acceptance, history, and
  growth policy;
- `OptCycleG` for cyclic scheduling, restart history, failure policy, and
  saving;
- `FullOpt1G` for simultaneous optimization, overlap penalties, Hessian I/O,
  and periodic saves;
- `EnergyGA`, `EnergyGAM`, and `EnergyGB` for normalized matrix elements,
  energy/gradient conventions, and MPI behavior.

Use I only to compare shifted inverse-iteration semantics. Do not copy I's
destructive storage of `H-shift*S`, its full-S layout, or its factorization
ordering into Q.

## 1. qrlinalg API findings

`src/qrlinalg.f90` is the API authority. Its real state represents

```text
H - sigma*S = Q*R
```

for one fixed shift `sigma`. It owns explicit dense Q and R factors plus all
factorization, update, and inverse-iteration workspace. It does not own or
retain H, S, or caller vectors, and it contains no MPI state.

The relevant public operations are:

| Operation | Required caller data | Effect and important restriction |
| --- | --- | --- |
| `initialize(capacity,info)` | Positive maximum order | Allocates two `capacity*capacity` factor arrays and workspace. It creates no valid factorization. Reinitializing destroys the previous state. |
| `factorize_fresh(H,S,shift,info,active_order=n)` | Capacity-sized H/S; lower triangle of the leading active block | Forms a fresh QR factorization. The optional active order avoids packed sections when capacity exceeds n. |
| `replace_symmetric(idx,delta_h,delta_s,info)` | Full conceptual length-n column changes | Replaces an arbitrary one-based row/column through two rank-one updates. It does not require the changed index to be last. |
| `append_symmetric(h_column,s_column,info)` | Full length-(n+1) physical columns including diagonals | Appends only at the end and requires an already valid nonempty state with spare capacity. The first basis function must use fresh factorization. |
| `delete_symmetric(idx,info)` | One valid principal index | Deletes an arbitrary row/column and preserves survivor order. It refuses to delete the only remaining function. |
| `solve(S,v_initial,x,lambda,tol,max_iter,norm_mode,rel_acc,num_iter,info)` | S with extents at least n and two distinct length-n vectors | Finds the root closest to the stored shift by inverse iteration. `norm_mode=0` returns `x^T*S*x=1`, which is required by the existing gradient formula. |
| `factorization_residual(H,S,v,absolute,relative,info)` | Capacity-sized physical H/S and a nonzero length-n probe | Compares `(H-shift*S)*v` with `Q*(R*v)` without allocation; it detects factor/physical-matrix drift but is not the eigenpair residual. |
| metadata and `clear()` | None | Query validity/order/capacity/shift/update counts or release the state. State components are private in production. |

Structural update routines validate ordinary argument errors before changing
the factors. The underlying qrupdate kernels do not return a numerical status.
A structurally successful update can still produce a singular shifted problem;
`solve` detects that later. `QR_ERR_NO_CONVERGENCE` returns an approximate
eigenpair, but Q optimizer code must still treat it as a failed evaluation
unless a separate acceptance policy is explicitly designed.

The shift must remain fixed while updates are reused. Changing
`Glob_ApproxEnergy` requires `factorize_fresh`; never silently retarget existing
factors. Initially Q should follow I's root-near-shift semantics. It cannot
promise G's `WHICH_EIGENVALUE` index semantics, especially for excited states.

### Completed qrlinalg API refinement

`BasisEnlQ` allocates H/S to `Kstop` while the active order grows from a smaller
value. Passing `Glob_S(1:n,1:n)` to the current assumed-shape `solve` interface
can create a packed temporary because the section is not contiguous when
`n<Kstop`. That copy can occur on every inverse-iteration solve. Fresh
factorization has the same shape issue, although it runs less frequently.

The imported library now provides the preferred solution: `factorize_fresh`
has optional `active_order`, while factorization, solve, and residual checks
use physical leading dimensions and ignore inactive storage. Exact-size calls
remain source-compatible. ECGPACK therefore does not need a second packed S
matrix or a quadratic packing copy on every solve.

## 2. Canonical basis and matrix layout

QR itself imposes no preferred basis-function permutation. It only requires
that factors, physical matrices, coefficients, parameters, and derivatives
refer to the same order. The “functions being optimized are the last rows”
rule comes from current G assembly and gradient code, not from qrlinalg.

Q therefore uses a stable physical basis order:

```text
canonical index i
    Glob_NonlinParam(:,i)   parameters of basis function i
    Glob_FuncNum(i)         existing function identity/history value
    Glob_c(i)               coefficient of normalized function i
    Glob_diagS(i)           raw diagonal overlap <phi_i|phi_i>
    Glob_H(i,j), i>=j       normalized, unshifted physical H
    Glob_S(i,j), i>=j       normalized physical S
    QR factors              same order, representing H-sigma*S
```

The lower triangles, including diagonals, are authoritative:

```text
Glob_H(i,i) = <phi_i|H|phi_i> / Glob_diagS(i)
Glob_S(i,i) = 1
Glob_H(i,j) = <phi_i|H|phi_j> /
              sqrt(Glob_diagS(i)*Glob_diagS(j)), i>j
Glob_S(i,j) = <phi_i|phi_j> /
              sqrt(Glob_diagS(i)*Glob_diagS(j)), i>j
```

The upper triangles are unspecified and must not be read by Q. A conceptual
element is always `A(max(i,j),min(i,j))`. qrlinalg replacement needs a full
column, so `GatherQCanonicalColumn` performs an O(n) symmetric gather.
`StoreQCanonicalColumn` writes both halves of a conceptual column into their
proper locations in the lower triangle while leaving the upper triangle
untouched.

This layout differs deliberately from both existing methods:

- G keeps the normalized H diagonal in `Glob_diagH`, leaves the S diagonal in
  `Glob_diagS`, and reconstructs matrix diagonals immediately before DSYGVX.
  DSYGVX may overwrite its chosen upper triangles.
- I stores the already shifted H matrix and mirrors S into both triangles.

For Q, H/S must always remain physical and unshifted. `Glob_diagH` can be kept
as a synchronized compatibility cache while shared swap/output routines still
need it, but Q must regard `Glob_H(i,i)` as authoritative. `Glob_diagS` retains
its existing raw-norm meaning; it is not a copy of `Glob_S(i,i)`.

### Why the G permutations must not be copied

`OptCycleG` currently:

1. reverses `FuncBegin:FuncEnd` around lines 5002-5005;
2. moves that range to the end around lines 5007-5011;
3. repeatedly rotates trailing blocks around lines 5079-5100;
4. sorts and reverses everything back around lines 5326-5336.

`FullOpt1G` similarly moves an interior range to the end around lines
6121-6137 and moves it back near the end of the routine. These operations let
the existing APIs express active functions as the suffix
`Glob_nfru+1:Glob_nfa`.

Calling those permutations in Q would require applying the identical symmetric
permutation to QR factors. qrlinalg has no public permutation operation, so a
fresh O(n^3) factorization would be required after every scheduling movement.
That would defeat updates. Q must instead keep the basis fixed and schedule
through explicit active indices.

Do not refactor G merely to make Q possible. Implement and validate the indexed
Q path first. A later shared refactor may let G use it after numerical and
performance equivalence have been demonstrated.

## 3. Active-index and derivative layout

`Q_Workspace%ActiveFunction(a)` maps optimizer block `a` to canonical function
index `i`. `ActivePosition(i)` is the inverse, or zero for inactive functions.
The order of `ActiveFunction` defines all of the following:

```text
x((a-1)*Glob_npt+1:a*Glob_npt)
gradient((a-1)*Glob_npt+1:a*Glob_npt)
DRMNG scaling and Hessian block a
Glob_D(:,a,j)
```

`Glob_D(:,a,j)` keeps the existing G scaled-derivative convention, but its
second index means active position rather than `i-Glob_nfru`. Its third index
is the canonical partner function `j`. The energy-gradient contraction must
therefore use:

```text
i = ActiveFunction(a)
coefficient of differentiated function = Glob_c(i)
diagonal derivative                    = Glob_D(:,a,i)
partner derivative                     = Glob_D(:,a,j)
```

Do not reinterpret `Glob_D` as a derivative of the already normalized H/S
matrices. `StoreHSD` scales raw derivatives and the current gradient applies
the diagonal normalization correction separately. Preserve that algebra until
it is rederived and finite-difference tested.

### Indexed assembly algorithm

The current `ComputeMatElem(Nmin,Nmax)` and
`ComputeMatElemAndDeriv(Nmin,Nmax)` are suffix-oriented. They only visit pairs
whose higher index lies in `Nmin:Nmax`, which misses `(inactive later,
active earlier)` pairs for an arbitrary active list.

The Q assembly routine must:

1. validate the active map and copy all trial parameters into their canonical
   `Glob_NonlinParam` columns;
2. calculate every active diagonal first and stage its raw overlap norm;
3. enumerate each unordered off-diagonal pair `(i,j)`, `i>j`, exactly once
   when either endpoint is active;
4. request the derivative for endpoint i only if i is active and for endpoint
   j only if j is active;
5. normalize with the staged raw norm for active endpoints and existing
   `Glob_diagS` for inactive endpoints;
6. store results in full staged columns without modifying canonical H/S;
7. preserve the existing symmetry-term distribution and MPI reductions;
8. preserve the no-gradient fast path in `MatrixElementsHS_RG_0S`.

The essential correctness test is one active function in the middle of the
basis. Its staged column must contain interactions with both earlier and later
functions. A second test must use two noncontiguous active functions and check
their shared intersection.

The indexed assembly is kept in `workproc.f90`, as requested for this port,
because its active map and transaction buffers are workproc-owned. Only the
small Q diagonal/off-diagonal storage branch belongs in `matform.f90`'s
`StoreHS`. This leaves both suffix-oriented G/I assembly loops unchanged and
does not create a separate qmethods module.

## 4. Memory ownership

Physical matrices and basis data remain replicated on all MPI ranks, matching
the existing assembly and distributed gradient calculations. qrlinalg is
serial, so its state and solve vectors should exist only on rank zero.

`QMethodWorkspace` in `workproc.f90` owns:

| Data | Shape and ranks | Purpose |
| --- | --- | --- |
| `ActiveFunction`, `ActivePosition` | max-active and capacity; all ranks | Stable optimizer/canonical mappings. |
| `PreviousH/S`, `TrialH/S` | `(capacity,max_active)`; all ranks | Complete conceptual columns for transaction recovery and final staged columns. |
| `PreviousDiagS`, `TrialDiagS` | max-active; all ranks | Raw norms corresponding to old and trial columns. |
| `MatrixParam`, `PreviousParam`, `AcceptedParam` | `(Glob_npt,max_active)`; all ranks | Parameters represented by current matrices, pre-transaction recovery parameters, and the explicitly accepted optimizer point. They must remain distinct from the next trial in `Glob_NonlinParam`. |
| `qr_real_state` | two capacity-square factors; rank zero only after lifecycle implementation | Private factorization and qrlinalg workspace. The type is already owned by `QMethodWorkspace`; non-root instances remain empty. |
| `InitialVector`, `SolvedVector`, `DeltaH/S` | capacity; rank zero only | Nonaliasing inverse-iteration vectors and one sequential replacement delta. |

The workspace does not duplicate physical H/S. It also does not stage a second
derivative tensor: `Glob_D(2*Glob_npt,max_active,matrix_order)` describes the
currently assembled trial and can be overwritten on the next trial. If failure
recovery later requires retaining derivatives, prefer recomputation over an
automatic second tensor because derivative storage may dominate memory.

Approximate aggregate matrix storage for scalar byte size `b`, capacity `C`,
active order `n`, active count `m`, and process count `P` is:

```text
replicated physical H/S       2*b*C*C*P
root qrlinalg Q/R             2*b*C*C
replicated transaction cols   4*b*C*m*P
replicated derivatives        2*b*Glob_npt*m*n*P
```

This estimate excludes small vectors, optimizer Hessians, matrix-element
buffers, and existing global data.

### Architecture gate result

The post-import audit found no architectural blocker for `OptCycleQ`.

- One replacement is two rank-one QR updates and costs O(n^2). A fixed-size
  cyclic block therefore keeps one energy evaluation O(n^2) in linear algebra.
- Indexed matrix-element assembly visits `m*n-m*(m-1)/2` physical pairs and
  performs two bulk reductions of the staged `(capacity,max_active)` arrays.
  It does not introduce a collective per matrix element.
- Stable canonical order removes all cyclic scheduling permutations and the
  O(n^3) refactorizations they would force.
- Physical H/S remain replicated exactly as in the existing program. The two
  extra dense Q/R arrays live only on rank zero. For `FULL_OPT1`, transaction
  and derivative storage is O(n^2), but the existing derivative tensor and
  DRMNG Hessian are already the dominant memory consumers.
- A dense full-basis change still costs O(n^3), whether expressed as n column
  replacements or one fresh factorization. This is a FULL_OPT1 policy choice,
  not an OPT_CYCLE bottleneck.

The audit did identify and correct one state-model omission: current matrix
parameters cannot live only in `Glob_NonlinParam`, because DRMNG overwrites
that array with the next trial before H/S changes. `MatrixParam` and
`MatrixParametersAreStored` now preserve the represented generation.

The refresh mechanism is also sufficient, but its final constants belong to
the factor-lifecycle benchmark chunk. The initial policy to test is:

1. factorize fresh after initialization, restart, any shift change, an
   uncertain rollback, or an order/metadata mismatch;
2. compute the physical eigenpair residual after every solve;
3. periodically use a dense deterministic probe for
   `factorization_residual`, not only the eigenvector, which could miss drift
   in another direction;
4. use both a precision-scaled residual threshold and an O(n) update-count
   guard so the amortized refresh cost remains O(n^2) per replacement.

The qrlinalg order-eight stress test through 10,000 reversible cycles showed
smooth, approximately accumulated-roundoff growth rather than an abrupt
stability boundary. For real wp=8 moderate replacements, the reported factor
residual was about `7e-13` after 20,000 public replacements. This supports a
measured caller policy but does not justify using 20,000 as a universal limit;
ECG matrices have different orders and conditioning.

## 5. Trial transaction and accepted-point rules

The parameters, physical matrices, coefficients, and QR factors must always be
treated as related generations. qrlinalg's `is_valid()` cannot detect edits to
caller-owned H/S. The architectural audit therefore requires these independent
workproc states:

```text
MatricesAreCanonical       lower H/S and raw diagS obey Q storage
MatrixParametersAreStored  MatrixParam identifies the active parameters in H/S
TrialIsReady               all final normalized trial columns are staged
TrialHasDerivatives        Glob_D corresponds to that same staged trial
FactorsMatchMatrices       root Q/R represents the current physical H/S and shift
```

No routine may infer one flag from another. In particular, a successful matrix
assembly does not imply factor validity, and changing `Glob_NonlinParam` to a
DRMNG trial does not invalidate `MatrixParam`, which records the physical
matrices' preceding generation.

For a trial changing `m` functions:

1. Require a captured `MatrixParam`, then save current conceptual columns, raw
   norms, and represented parameters in `Previous*`.
2. Assemble every final target column using all new parameters and all new raw
   norms. Do not build one target against partially old parameters.
3. Preflight every index, extent, factor order, capacity, shift, and MPI-side
   invariant before changing the first factor.
4. For active column a, gather the progressively current physical column and
   compute `delta=target-current`.
5. On rank zero call `replace_symmetric`, broadcast status, then commit that
   column to physical H/S on every rank.
6. At an intersection with a column committed earlier, the later delta must be
   zero. Computing all deltas from the original matrix double-counts shared
   intersections.
7. Solve only after all m columns and factors represent the same final trial.

DRMNG may request many rejected points. It is valid for factors and physical
matrices to follow the latest evaluated point as long as they remain matched;
the next trial can update directly from it. Separately retain the best or
accepted parameter vector. When the driver decides to accept `x_best`, or to
restore `x_init` after overlap/coefficient rejection, reassemble that point and
make matrices and factors follow it as one transaction. Restoring only
nonlinear parameters is incorrect.

All ordinary qrlinalg replacement errors are caught before its first rank-one
update, so a complete preflight makes a mid-column recoverable failure
impossible. If a later multi-column failure or solve failure requires rollback,
restore physical columns and parameters and apply inverse column changes in
reverse order. Immediately check factor drift; fall back to fresh factorization
when the check fails. A partially modified state whose relationship cannot be
proved must never be reused.

## 6. Mapping each BBOP from G to Q

### OPT_CYCLE

Start here after low-level assembly/update/solve validation because it changes
a fixed-size basis and exercises repeated replacement naturally.

Keep G's input checks, DRMNG settings, counters, acceptance rules, saving
frequency, and `Glob_History` meaning. Replace all physical permutations with:

```text
CurrFunc = next original/canonical function from history
NumActive = min(FuncEnd-CurrFunc+1,NumOfFuncToOpt)
ActiveFunction(a) = CurrFunc+NumActive-a, a=1:NumActive
```

The descending optimizer map reproduces the block order created when G first
reverses `FuncBegin:FuncEnd` and then moves it to the trailing matrix block.
It is independent of the ascending physical matrix order. Verify the map
against printed `Glob_FuncNum` values from G for partial ranges, overlapping
blocks, final short blocks, `NumOfFuncToShift>1`, and restart in the middle of
a cycle. Q output no longer needs `SaveResults(Sort='yes')` merely to undo
internal optimization order, because canonical order never changed.

### BASIS_ENL

Keep G's generation and acceptance logic. The accepted prefix stays canonical;
new functions occupy a tentative suffix.

- For an empty basis, assemble the first nonempty block and factorize fresh;
  qrlinalg cannot append to order zero.
- When a block grows an existing accepted basis, append each new function in
  order after all target raw norms and columns are available.
- Later random candidates at the same tentative order are replacements, not
  appends.
- `Glob_CurrBasisSize` remains the accepted size. Track the represented
  tentative order separately in `Q_Workspace%MatrixOrder`.
- If a candidate block is rejected, restore the accepted prefix/order and
  factors. Deleting tentative suffix functions in descending order is valid;
  fresh factorization is the initial recovery fallback.
- qrlinalg capacity is `Kstop`, matching G's H/S allocation.

### FULL_OPT1

Use `ActiveFunction=InitFunc:FinalFunc` and do not move the range to the end.
The x/gradient/Hessian order remains the user's ascending range, which must be
recorded consistently in Hessian I/O and restarts.

Replacing m columns costs O(m*n^2). When m is proportional to n, that is cubic
and may be slower or less stable than fresh factorization. The first correct
implementation may explicitly factorize fresh for dense trials. Benchmark
fresh versus sequential updates at representative n and m before choosing a
crossover. Document the selected rule and do not claim that a full-basis
change becomes quadratic through row/column updates.

### Other BBOP actions

These routines must be audited after the three optimizer drivers:

| BBOP | Existing implementation and Q work |
| --- | --- |
| `ELIM_LCFN` | `EliminateLittleContribFunc`: retain its coefficient-based selection policy; apply principal deletions in descending index order and compact every basis-indexed array consistently. |
| `ELIM_LND1` | `EliminateLinDepFunc`: retain overlap-mask semantics; use principal deletion rather than rebuilding unaffected elements. Define behavior if every function would be deleted because qrlinalg has no order-zero state. |
| `SEPR_LND1` | `SeparateLinDepFunc`: active indices are the selected mask; assemble replacements and solve in canonical order. |
| `SEPR_FLCF` | `SeparateFuncLargeCoeff`: same replacement path for coefficient-selected indices. |
| `EXPC_VALS`, `DENSITIES`, `MOMT_DENS` | `ExpectationValues`: build/solve Q once, retain normalized coefficients and existing operator calculations. Q returns one shift-selected root, unlike G's index selection. |
| `SAVE_HSWF` | `SaveHSWF`: solve through Q but export physical unshifted H/S and coefficients in canonical order. |
| `SAVE_FILE` | No eigenproblem. Confirm canonical parameter/history order and do not create Q factor state. |

Elimination and separation currently dispatch only G and reject I. Adding Q is
new behavior, not an alias to an existing iterative branch.

## 7. Required routine changes

The following is the working implementation map. Line numbers are intentionally
omitted because this file will grow; search by routine name.

| File/routines | Change |
| --- | --- |
| `workproc.f90`: Q workspace/helpers | Add root factor lifecycle after qrlinalg build integration; keep active mapping, transaction, solve, acceptance, and BBOP orchestration here. |
| `matform.f90`: `StoreHS` | Q storage is implemented: physical normalized H/S, including diagonals, go directly into the lower triangles and raw S norm remains in `Glob_diagS`. |
| `workproc.f90`: `AssembleQTrial` | Non-derivative arbitrary-index staging is implemented. Add indexed derivative placement later without changing the suffix-oriented G/I loops. |
| `workproc.f90`: `EnergyGA/GAM/GB` | Use as behavioral references. Fill `EnergyQA/QAM/QB`; do not change G until Q equivalence is demonstrated. |
| `workproc.f90`: overlap penalty routines | Add active-index versions or generalize carefully. Current loops and gradient offsets assume a suffix. Count each unordered pair once. |
| `workproc.f90`: swap routines | Treat Q like canonical unshifted G data but with physical diagonals already in H/S. Serialize through upper-triangle packing without changing the authoritative lower triangle. Rebuild factors after restart; never serialize private Q/R. |
| `workproc.f90`: permutation/sort helpers | Q optimization scheduling must not call them. If a user-visible operation truly reorders the basis, permute all basis-indexed arrays then factorize fresh. |
| `workproc.f90`: `GetOverlapStatistics` | Add an active-index form for noncontiguous sets; the current `Nmin:Nmax` form remains suitable for contiguous FULL_OPT1. |
| `workproc.f90`: `BasisEnlQ`, `OptCycleQ`, `FullOpt1Q` | Fill one driver per review chunk from its G counterpart. |
| `main.f90` | Add Q cases only after the corresponding driver is operational. Unsupported BBOP/method combinations must produce an explicit error. |
| `globvars.f90` | Update comments for Q physical matrix invariants if Q storage is selected. Avoid moving workproc-only transaction state into globals. |
| `Makefile` | Compile qrupdate modules and qrlinalg in dependency order with preprocessing; select one BLAS/LAPACK provider; make workproc depend on qrlinalg. |
| shared documentation | Document Q as shift-targeted and list exactly which BBOP actions are enabled. |

Input parsing and `SaveResults` already preserve the one-character method token;
they need no format change for Q. Validation/dispatch does need explicit Q
coverage when routines become usable.

## 8. Build, provider, and qrlinalg API integration

The second implementation chunk imported `src/qrlinalg.f90` from qrlinalg
commit `18c587d` (`Support active matrix storage and QR drift checks`). The
vendored file and the upstream file had SHA-256
`c8f38cafd6f9c8a2ebeef586bd30186a0de2eceb72d9bad9f19a1fe7a1d50e0c` at
import time. A recursive comparison showed that the already vendored
`src/qrupdate/` directory was unchanged, so it was not recopied.

The Makefile now compiles the dependencies in this order:

```text
wp_def
  -> qrupdate_linalg
  -> qrupdate_error
  -> qrupdate_real, qrupdate_complex
  -> qrupdate
  -> qrlinalg
  -> workproc
```

`qrlinalg.f90` contains preprocessing around private testing access. The
Makefile therefore selects `-cpp`, `-fpp`, or `-Mpreprocess` for gfortran,
ifort/ifx, or nvfortran respectively. Production does not define
`QRLINALG_TESTING`. The fixed-form qrupdate aggregates also use the matching
compiler option that removes the 72-column source limit.

Only one BLAS/LAPACK provider is linked. A definition-name audit found 36 BLAS
and 81 LAPACK routines in the old ECGPACK aggregates. The qrupdate aggregates
contain all of them, plus `DROT` in BLAS and `DGEQR2`, `DGEQRF`, `ZGEQR2`,
`ZGEQRF`, `ZLARTG`, and `ZROT` in LAPACK. Both sets use `wp_def` and retain the
ECGPACK selectable-kind modifications. Consequently `LINALG=netlib` now
compiles `src/qrupdate/BLAS.f` and `src/qrupdate/LAPACK.f` as the existing
`BLAS.o` and `LAPACK.o` targets. The old aggregates remain in the source tree
but are not linked. Optimized precision-8 builds continue to exclude bundled
BLAS/LAPACK and resolve both ECGPACK and qrlinalg calls from the selected
external provider.

The imported API removes the two storage blockers identified in the first
chunk:

- `factorize_fresh(H,S,shift,info,active_order=n)` accepts matrices whose
  physical leading dimensions are at least `n`, and reads only the lower
  triangle of the active leading principal block;
- `solve(S,...)` accepts the same capacity-sized overlap storage and passes
  its physical leading dimension to BLAS;
- `factorization_residual(H,S,v,absolute,relative,info)` checks
  `(H-shift*S)*v-Q*(R*v)` without allocating and without reading inactive or
  upper-triangle storage;
- `get_updates_since_fresh()` exposes the structural-update count needed by
  the caller's refresh policy.

This distinction must remain explicit: `factorization_residual` monitors QR
factor drift, while the physical eigenpair residual
`||H*c-E*S*c||` monitors the inverse-iteration result. `workproc` must evaluate
both at the checkpoints defined in the lifecycle chunk.

The vendored qrupdate code states GPL-3.0-or-later licensing. Preserve its
notices and confirm project-wide distribution compatibility before release;
update `THIRD-PARTY-NOTICES.md` when integration becomes distributable.

## 9. Reviewable implementation sequence

Each chunk stops for review and is committed only on explicit request.

1. **Completed in this working tree:** API study, canonical stable-order
   decision, memory/code scaffold in workproc, routine map, and porting record.
2. **Completed in this working tree:** qrlinalg shape/residual API import,
   build/provider integration, all-precision compilation, analytical residual
   test, and an existing I-method sample check.
3. **Completed in this working tree:** Q canonical `StoreHS`, explicit matrix-
   parameter generation state, and staged arbitrary-index non-derivative H/S
   columns compared with full G/Q recomputation under MPI.
4. **Root factor lifecycle and solve:** fresh factorization, status broadcasts,
   S-normalized coefficients, residual checks, and single-function behavior.
5. **Sequential replacement transactions:** repeated trials, intersecting
   active columns, restoration after failed solves, and refresh diagnostics.
6. **Indexed gradients and penalties:** compare with G suffix gradients and
   finite differences for interior/noncontiguous active sets under MPI.
7. **OPT_CYCLE Q:** preserve G scheduling/history semantics without physical
   permutations; test restart and final short blocks.
8. **BASIS_ENL Q:** append lifecycle, repeated candidate replacement, accepted
   order versus tentative order, and rejection recovery.
9. **FULL_OPT1 Q:** dense-change policy, Hessian ordering/restart, penalties,
   and periodic saves.
10. **Remaining BBOP actions and shared docs:** enable one operation only after
    its explicit Q path is validated.
11. **Port to sibling real codes:** repeat the checklist below for each basis
    variant; treat CG_0S separately because Hermitian conjugation and complex
    gradients require their own derivation.

## 10. Validation checklist for every implementation and port

At minimum, record exact commands and results for:

- syntax/build at precision 8 and, with the bundled generic provider,
  precisions 10 and 16;
- one and multiple MPI processes;
- order 1, order 2, `n<capacity`, and growth to capacity;
- active first, middle, last, contiguous, and noncontiguous functions;
- replacement of two active functions sharing one matrix intersection;
- append, arbitrary deletion, and stable survivor order;
- repeated trial points followed by restoration of `x_best` and `x_init`;
- singular or nearly singular shifted systems and non-convergence status;
- physical residual `||H*c-E*S*c||` and `c^T*S*c`;
- fresh versus updated factor results within expected roundoff;
- G/Q energy and analytic-gradient agreement where both target the same root;
- finite-difference gradients and overlap penalties;
- swap-file restart, optimizer-history restart, saved basis order, and Hessian
  variable order;
- representative performance and memory at realistic basis sizes.

For RG_1P, RG_2D, and RG_2P, identify the local basis-qualified matrix-element
routine and verify its diagonal/derivative conventions before copying indexed
assembly. Do not assume `Glob_npt`, extra basis metadata, or normalization
details are identical. Preserve each code's no-gradient fast path and MPI
symmetry routing.

## 11. First-chunk change log

Changed:

- `src/workproc.f90`: added `QMethodWorkspace`, allocation/cleanup, validated
  active-index maps, canonical symmetric column gather/store helpers, and
  fail-closed stubs for trial assembly, QR lifecycle, Q energy wrappers, and
  the three requested Q BBOP drivers.
- `docs/Q_METHOD_DESIGN.md`: replaced the disposable prior design with this
  canonical layout, routine map, implementation sequence, and porting log.
- removed disposable `src/qcontext.f90`, `src/qmethods.f90`,
  `docs/Q_CODE_STRUCTURE.md`, and `docs/ECGPACK_INTEGRATION.md`; their useful
  information was consolidated here.

Not changed:

- no G/I numerical code or permutations;
- no matform storage/assembly behavior;
- no qrlinalg/qrupdate implementation;
- no Makefile, main dispatch, or input behavior;
- no generated build files and no commit.

Validation performed:

```bash
make debug COMPILER=gfortran MACHINE=linux-generic PREC=8 \
  LINALG=netlib EXEFILE=ecg
```

The existing executable compiled and linked successfully. The compiler emitted
only the pre-existing legacy-format warning in `workproc.f90`.

A temporary debug caller under `/tmp` linked against the existing objects and
successfully checked:

- allocation for `MatrixOrder=3`, `Capacity=5`, and `MaxActive=2`;
- active list `(1,3)` and inverse map `(1,0,2)`;
- rejection of duplicate active indices without changing the previous map;
- gathering conceptual column 2 from a lower-only 3-by-3 matrix;
- storing column 2 only in canonical lower-triangle locations;
- deterministic workspace cleanup.

The generated debug objects, modules, and executable were removed with
`make cleaner`. No numerical QR test was run because qrlinalg is intentionally
not connected in this chunk.

## 12. Second-chunk change log

Changed:

- imported the updated `src/qrlinalg.f90` from upstream commit `18c587d`;
- verified that the already present `src/qrupdate/` tree was identical to the
  upstream tree and left it unchanged;
- updated `Makefile` with compiler-specific preprocessing/fixed-form flags,
  qrupdate/qrlinalg dependency order, and one-provider Netlib selection;
- added `qr_real_state` ownership to `QMethodWorkspace` and mapped temporary Q
  status constants to qrlinalg's public status values.

Provider audit:

- the qrupdate BLAS routine set is a strict superset of the original ECGPACK
  set: 37 versus 36 definitions, with only `DROT` added;
- the qrupdate LAPACK set is a strict superset: 87 versus 81 definitions, with
  only the real/complex QR factorization and rotation prerequisites added;
- an explicit link excluding bundled objects and using `-llapack -lblas`
  succeeded, confirming the optimized-provider path resolves the new symbols.

Validation performed:

```bash
make debug COMPILER=gfortran MACHINE=linux-generic PREC=8 \
  LINALG=netlib EXEFILE=ecg
make release COMPILER=gfortran MACHINE=linux-generic PREC=10 \
  LINALG=netlib EXEFILE=ecg_qr_wp10
make release COMPILER=gfortran MACHINE=linux-generic PREC=16 \
  LINALG=netlib EXEFILE=ecg_qr_wp16
```

All three builds compiled and linked. The upstream
`test_factorization_residual` suite was compiled against this project's wp=8
objects and passed, including capacity-sized matrices, inactive NaN sentinels,
replacement/append/delete lifecycle, and real/complex residual checks.

The existing 100-function Li expectation-value sample reached and solved its
I-method generalized eigenproblem with energy
`-7.4773765161376833`. It later stopped at the pre-existing malformed format
`'(a,i1,i1.1x)'` in `ExpectationValues`; that unrelated debug-runtime failure
was not changed in this Q chunk.

Not changed:

- no Q matrix assembly, factor initialization, solve, update, or BBOP dispatch;
- no qrlinalg or qrupdate algorithm was edited after import;
- no commit was created.

## 13. Third-chunk change log and architecture gate

Changed:

- `matform.f90`: `StoreHS` now has a Q branch that writes normalized,
  unshifted physical H/S into the canonical lower triangles, including
  `H(i,i)` and exact `S(i,i)=1`, while preserving raw `Glob_diagS`;
- `workproc.f90`: added `MatrixParam` and explicit parameter/trial generation
  flags, plus `CaptureQMatrixParameters`;
- `workproc.f90`: implemented the non-derivative `AssembleQTrial` path for any
  distinct active indices. It enumerates every touched unordered pair once,
  keeps the matrix-element kernel's no-gradient fast path, performs two bulk
  MPI reductions, normalizes only after all active raw norms are available,
  and leaves physical H/S unchanged.

Architecture conclusion:

- the stable-order design supports O(n^2) QR linear algebra per OPT_CYCLE
  energy evaluation when the active block size is fixed;
- no required qrlinalg operation assumes an active function is last;
- capacity-sized matrices no longer cause packed solve/factorization copies;
- transaction state, MPI ownership, swap representation, order-zero growth,
  root-selection semantics, and dense FULL_OPT1 cost all have explicit paths;
- no unresolved architectural bottleneck remains before factor lifecycle and
  replacement implementation. Refresh thresholds still require
  representative ECG benchmarks, as planned.

Validation used the 100-function Li basis and active canonical functions
`(1,50,100)`. Each active parameter block was changed by a different small
factor. A temporary caller compared staged columns with a complete Q matrix
recomputation and also compared an unmodified full Q matrix with G's physical
matrix convention.

```text
1 MPI rank: G/Q relative difference       0.0000e+00
            staged/full relative difference 2.3581e-14
            physical staging mutation       0.0000e+00
2 MPI ranks: G/Q relative difference      0.0000e+00
             staged/full relative difference 3.1443e-14
             physical staging mutation      0.0000e+00
```

The one-rank check was repeated with workspace capacity `n+7`; it produced the
same result, exercising nonpacked transaction storage. A two-rank check with
the descending active map `(100,50,1)`, which matches OPT_CYCLE scheduling,
also passed with staged/full relative difference `1.9654e-15` and zero physical
matrix mutation. The normal wp=8 debug build passed with only the pre-existing
format warning.

Not changed:

- indexed derivatives and overlap-penalty gradients remain fail-closed;
- `ApplyQTrial`, `FactorizeQFresh`, `SolveQ`, and all Q energy/BBOP drivers
  remain stubs;
- Q is not accepted by `main.f90`, so user-visible behavior is unchanged;
- no commit was created.
