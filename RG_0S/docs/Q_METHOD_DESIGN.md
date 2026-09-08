# Q-method implementation design and porting log

## Status and scope

This document is the implementation record for adding the `Q` generalized
eigenvalue method to `RG_0S`. It is also intended to be a recipe for porting
the same method to the other ECGPACK variants. Update it after every reviewable
chunk with the exact files changed, validation performed, open questions, and
any decision that a later port must reproduce.

The first committed foundation (`8561b61`) defines code ownership and canonical
layout, integrates qrlinalg and its build dependencies, and implements
non-derivative canonical Q storage plus arbitrary-index staged column assembly.
The second committed chunk (`4b07fa4`) completes `OPT_CYCLE Q`: factor
lifecycle, replacement transactions, solve/residual monitoring, indexed
derivatives, energy/gradient wrappers, swap handling, driver, and `main.f90`
dispatch are implemented and validated.

The current review chunk implements the remaining Q-facing work in `RG_0S`:
`BASIS_ENL`, `FULL_OPT1`, indexed overlap penalties, `EXPC_VALS`, `DENSITIES`,
`MOMT_DENS`, `SAVE_HSWF`, and all four elimination/separation commands. The
optimizer drivers exercise append and replacement updates. The one-shot
elimination/separation routines use qrlinalg but deliberately rebuild and
factorize the resulting matrix, preserving their established full-recompute
behavior. Sections 15 and 16 record the exact implementation and validation.

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
eigenpair. `SolveQ` accepts that approximation only if the independent
physical generalized-eigenpair residual satisfies the requested tolerance;
otherwise it refreshes the factors and retries once.

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

These routines were audited after the three optimizer drivers. The resulting
implementation choices are:

| BBOP | Existing implementation and Q work |
| --- | --- |
| `ELIM_LCFN` | `EliminateLittleContribFunc`: retain its coefficient-based selection policy, compact the selected nonlinear parameters, rebuild canonical H/S for the compacted basis, and factorize fresh. A future performance-only refinement may pair `delete_symmetric` with in-place matrix compaction. |
| `ELIM_LND1` | `EliminateLinDepFunc`: retain overlap-mask semantics, rebuild canonical H/S for the survivors, and factorize fresh. Define behavior if every function would be deleted before attempting qrlinalg order zero. |
| `SEPR_LND1` | `SeparateLinDepFunc`: retain the selected mask and random perturbation policy, then rebuild canonical H/S and factorize fresh. A replacement-only path is a possible later optimization. |
| `SEPR_FLCF` | `SeparateFuncLargeCoeff`: use the same full-rebuild policy for coefficient-selected functions. |
| `EXPC_VALS`, `DENSITIES`, `MOMT_DENS` | `ExpectationValues`: build/solve Q once, retain normalized coefficients and existing operator calculations. Q returns one shift-selected root, unlike G's index selection. |
| `SAVE_HSWF` | `SaveHSWF`: solve through Q but export physical unshifted H/S and coefficients in canonical order. |
| `SAVE_FILE` | No eigenproblem. Confirm canonical parameter/history order and do not create Q factor state. |

Elimination and separation now dispatch G and Q and continue to reject I. Q is
new behavior, not an alias to the existing iterative branch.

## 7. Required routine changes

The following is the working implementation map. Line numbers are intentionally
omitted because this file will grow; search by routine name.

| File/routines | Change |
| --- | --- |
| `workproc.f90`: Q workspace/helpers | Add root factor lifecycle after qrlinalg build integration; keep active mapping, transaction, solve, acceptance, and BBOP orchestration here. |
| `matform.f90`: `StoreHS` | Q storage is implemented: physical normalized H/S, including diagonals, go directly into the lower triangles and raw S norm remains in `Glob_diagS`. |
| `workproc.f90`: `AssembleQTrial` | Arbitrary-index H/S staging and derivative placement are implemented without changing the suffix-oriented G/I loops. |
| `workproc.f90`: `EnergyGA/GAM/GB` | Behavioral references for implemented `EnergyQA/QAM/QB`; G is unchanged apart from an independent uninitialized-best-point correction in `BasisEnlG`. |
| `workproc.f90`: overlap penalty routines | Indexed Q energy, gradient, and statistics variants enumerate each unordered pair touching the active map exactly once. |
| `workproc.f90`: swap routines | Q uses canonical unshifted data with physical diagonals already in H/S. Serialization borrows the unused upper triangle without changing authoritative lower storage; private Q/R is always rebuilt. |
| `workproc.f90`: permutation/sort helpers | Q optimization scheduling must not call them. If a user-visible operation truly reorders the basis, permute all basis-indexed arrays then factorize fresh. |
| `workproc.f90`: `GetOverlapStatistics` | `GetQOverlapStatistics` covers arbitrary active maps; the old contiguous routine remains unchanged for G/I. |
| `workproc.f90`: `BasisEnlQ`, `OptCycleQ`, `FullOpt1Q` | All three drivers are implemented with stable canonical basis order. |
| `main.f90` | Q dispatch is enabled for every method-taking BBOP after its implementation; I remains explicitly unsupported by the four legacy elimination/separation commands. |
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

1. **Completed in foundation commit `8561b61`:** API study, canonical stable-order
   decision, memory/code scaffold in workproc, routine map, and porting record.
2. **Completed in foundation commit `8561b61`:** qrlinalg shape/residual API import,
   build/provider integration, all-precision compilation, analytical residual
   test, and an existing I-method sample check.
3. **Completed in foundation commit `8561b61`:** Q canonical `StoreHS`, explicit matrix-
   parameter generation state, and staged arbitrary-index non-derivative H/S
   columns compared with full G/Q recomputation under MPI.
4. **Completed for OPT_CYCLE Q:** root factor lifecycle and solve: fresh factorization, status broadcasts,
   S-normalized coefficients, residual checks, and single-function behavior.
5. **Completed for OPT_CYCLE Q:** sequential replacement transactions: repeated trials, intersecting
   active columns, restoration after failed solves, and refresh diagnostics.
6. **Indexed gradients completed for OPT_CYCLE Q:** compare with finite
   differences for interior and multi-function active sets under MPI. Indexed
   overlap-penalty gradients remain part of `FULL_OPT1 Q`.
7. **Completed:** `OPT_CYCLE Q` preserves G scheduling/history semantics
   without physical permutations and supports final short blocks and swap
   restart.
8. **Completed:** `BASIS_ENL Q` append lifecycle, repeated candidate
   replacement, accepted order versus tentative order, and rejection recovery.
9. **Completed:** `FULL_OPT1 Q`, canonical Hessian ordering/restart, indexed
   penalties, and periodic saves.
10. **Completed in code:** remaining method-taking BBOP actions use Q. Runtime
    validation is complete for `SAVE_HSWF`, `EXPC_VALS`, `ELIM_LCFN`, and
    `ELIM_LND1`; the two separation commands await replay because the local
    MPI execution approval service became unavailable after compilation.
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

## 14. OPT_CYCLE Q implementation log

This section records the complete implementation after foundation commit
`8561b61`. The statements in sections 11--13 are historical chunk boundaries;
their “not changed” lists describe those earlier reviews, not current code.

### Step 1: factor ownership and fresh construction

`PrepareQWorkspace` now allocates the four capacity-length vectors used for
inverse iteration and column changes. The dense `qr_real_state` itself remains
root-owned. `FactorizeQFresh` checks the canonical matrix relationship and
physical extents on every rank, initializes root storage only when capacity
changes, calls

```fortran
call Factors%factorize_fresh(Glob_H,Glob_S,Glob_ApproxEnergy,info, &
                             active_order=MatrixOrder)
```

on rank zero, broadcasts the status, and changes `FactorsMatchMatrices` only
after collective success. This exact call form is important for future
`BASIS_ENL Q`: passing capacity-sized matrices with `active_order` avoids a
packed temporary when active order is smaller than capacity.

### Step 2: transactional symmetric replacements

`ApplyQTrial` preflights factor validity, order, capacity, and shift before the
first update. It then snapshots all original active physical columns, raw
overlap diagonals, and represented nonlinear parameters. Columns are applied
in active-map order:

1. gather the current conceptual physical column from the lower triangle;
2. subtract it from the staged final column;
3. call `replace_symmetric` on rank zero;
4. broadcast its status;
5. commit that final column to canonical H/S on every rank;
6. update the corresponding raw `Glob_diagS` and parameter generation.

The current column is gathered again at each step. Consequently, an
active--active intersection installed by an earlier replacement has zero delta
when the later endpoint is processed and is not applied twice. Although all
ordinary qrlinalg failures are eliminated by preflight, the defensive failure
path restores every physical snapshot and parameter block and constructs fresh
factors. It never reuses a factor state whose relationship is uncertain.

### Step 3: solve and two distinct residual monitors

`SolveQ` uses the preceding `Glob_c` as the inverse-iteration initial vector
and falls back to an all-ones vector only when it is numerically zero. It calls
qrlinalg with `norm_mode=0`, so the returned coefficient vector obeys
`c^T*S*c=1`, and broadcasts status, energy, coefficients, direction accuracy,
and iteration count.

Every successful solve evaluates two different diagnostics:

- `ComputeQEigenpairResidual` computes
  `||H*c-E*S*c||/(||H*c||+|E|*||S*c||+tiny)` from canonical H/S;
- `factorization_residual` compares the physical shifted action with `Q*(R*v)`
  using a fixed, alternating dense probe vector.

Do not use the eigenvector as the factor probe. Because the shift is close to
the desired eigenvalue, `(H-shift*S)*c` is nearly zero and makes the relative
factor-action ratio artificially sensitive. On the 100-function Li test, an
eigenvector probe reported about `1.6e-10` immediately after a fresh
factorization, while the deterministic probe reported about `3.9e-16` and the
physical eigenpair residual was about `2.6e-16`.

A factor or physical residual beyond its precision-scaled tolerance causes one
fresh-factorization and solve retry. qrlinalg's internal
`QR_ERR_NO_CONVERGENCE` approximation is also checked against the physical
residual and is accepted if that stronger condition passes. On a genuine
physical-residual failure, retry iterations are raised to
`max(120,Glob_MaxIterForGSEPIIS,4*n)` because an append can rotate the desired
eigenvector substantially. Continued physical-residual failure is reported as
`QR_ERR_NO_CONVERGENCE`; continued factor mismatch is reported as invalid
state. Independently, a fresh factorization is forced after
`max(64,8*n)` replacements. Because this guard is O(n), its O(n^3) refresh cost
amortizes to O(n^2) per replacement. `LastEigenpairResidual`,
`LastFactorResidual`, and `FreshFactorizations` preserve the latest diagnostics
in `Q_Workspace` for tests and future reporting.

### Step 4: arbitrary-index derivatives

`AssembleQTrial(.true.)` uses the same unordered-pair enumeration as the
non-derivative path. For active position `a` and canonical partner `j`, it
stores the G-compatible scaled raw derivatives in `Glob_D(:,a,j)`. When both
endpoints are active, the matrix-element kernel returns both endpoint
derivatives and the second one is placed through `ActivePosition`; no second
matrix-element evaluation is made.

Local derivative contributions are combined with one in-place MPI all-reduce
over `Glob_D`. This avoids another `2*npt*max_active*capacity` receive tensor.
Active raw diagonal norms are reduced before normalization. Off-diagonal
derivatives are multiplied by `1/sqrt(Sii*Sjj)`; diagonal derivatives are
multiplied by `2/Sii`, matching `StoreHSD` and accounting for both identical
endpoints.

`EnergyQB` contracts this tensor with canonical coefficients. Its diagonal
correction uses `Glob_D(:,a,ActiveFunction(a))`, so no suffix assumption
remains. `OPT_CYCLE` disables the smooth overlap penalty exactly as G does;
the Q energy wrappers deliberately return `QR_ERR_NOT_IMPLEMENTED` if a future
caller enables that penalty before the indexed `FULL_OPT1 Q` implementation.

### Step 5: cyclic driver without matrix permutations

`OptCycleQ` retains the input validation, DRMNG configuration, best-point
selection, evaluation limits, overlap/coefficient rejection, history fields,
saving frequency, and warning policy of `OptCycleG`. It removes all calls to
`ReverseFuncOrder`, `PermuteFunctions*`, `PermuteMatrixElements*`, and
`SortBasisFuncAndMatElem`.

For a step beginning at canonical `CurrFunc`, it constructs

```text
ActiveFunction(a) = CurrFunc+nfo-a
```

This descending order reproduces the optimizer block order created by G's
initial reversal while leaving every physical index fixed. Before DRMNG begins,
the driver captures the matrix parameter generation. Every energy or gradient
request copies the requested blocks into their mapped canonical functions and
executes a complete assemble/apply/solve transaction. At termination it
re-evaluates `x_best`; on failed solve, excessive overlap, or excessive linear
coefficient, it transactionally re-evaluates `x_init` rather than restoring
parameters alone.

Overlap rejection enumerates every unordered pair touching an active function.
This explicitly includes inactive canonical indices greater than the active
index, a case that G handles implicitly by moving active functions to the end.
Results are saved with `Sort='no'` because canonical order is already the
user-visible order.

### Step 6: swap format and dispatch

The on-disk swap representation is unchanged: H occupies the lower triangle
and diagonal, normalized S occupies the upper triangle, and raw `diagS` follows
the packed matrix. For Q, restore copies packed S into only its canonical lower
triangle and sets `Sii=1`; it does not shift H or populate an authoritative
upper triangle. Store borrows only the unused upper H triangle, so live lower
H/S and QR factors remain consistent.

`main.f90` dispatches `OPT_CYCLE Q` to `OptCycleQ` and reports a returned Q
status. At the historical boundary recorded by this section, `BASIS_ENL Q`
and `FULL_OPT1 Q` were still unsupported. They are implemented in the next
review chunk documented in sections 15 and 16.

The accepted input line has the same fields as G and I, with method `Q`:

```text
OPT_CYCLE Q K FuncBegin FuncEnd NumToOpt NumToShift NumCycles MaxEnergyEval \
              OverlapThreshold LinCoeffThreshold SavingFreq
```

### Validation performed

Strict and optimized wp=8 builds completed with only the pre-existing legacy
format warning in `ExpectationValues`. Optimized wp=10 and wp=16 builds also
completed successfully:

```bash
make debug   COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib EXEFILE=ecg
make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib EXEFILE=ecg
make release COMPILER=gfortran MACHINE=linux-generic PREC=10 LINALG=netlib EXEFILE=ecg_q_wp10
make release COMPILER=gfortran MACHINE=linux-generic PREC=16 LINALG=netlib EXEFILE=ecg_q_wp16
```

A disposable harness used the 100-function Li basis and canonical function 50.
All six `EnergyQB` components were compared with central differences through
the full staged-replacement and QR-solve path. The maximum scaled discrepancy
was `7.39e-5` at a finite-difference step of `1e-5*max(1,abs(parameter))`.
The same harness with descending, noncontiguous active functions `(100,50)`
validated all twelve gradient entries, including the shared active--active
matrix element, to the same scaled tolerance on both one and two MPI ranks.
The two-rank maximum scaled discrepancy was `2.05e-5`. The tiny function-100
gradient components were at the energy finite-difference resolution; their
absolute discrepancies were at most about `3e-11`.

End-to-end release runs on the same basis completed for:

- one MPI rank, functions 99--100 in separate one-function steps, with two
  energy and one gradient evaluation per step;
- two MPI ranks with the same schedule;
- two MPI ranks with one descending two-function block `(100,99)`, exercising
  both endpoint derivatives and sequential intersection handling;
- two consecutive Q BBOP steps, where the first wrote the swap file and the
  second restored it and reproduced the same initial energy.

No sample input was modified. All calculation inputs and gradient harnesses
used for these checks were disposable files under `/tmp`.

## 15. BASIS_ENL Q and FULL_OPT1 Q implementation log

This is the first uncommitted review chunk after `4b07fa4`. It completes the
two optimizer drivers that were intentionally left fail-closed in the
OPT_CYCLE chunk.

### Step 1: suffix order transactions for basis growth

`TrimQFactors(TargetOrder,...)` and `AppendQTrial(...)` separate accepted basis
order from tentative trial order.

- Rejected suffixes are removed from the factor state in descending order by
  repeated `delete_symmetric(last)`. Their obsolete physical slots remain
  outside `MatrixOrder` and are overwritten by the next candidate.
- An accepted nonempty prefix is grown one column at a time with
  `append_symmetric`. Each canonical physical column is committed only after
  the corresponding factor update succeeds.
- qrlinalg has no valid order-zero factorization. The first candidate block is
  therefore staged completely, committed to physical H/S, and factorized
  fresh. Later candidates at the same order use replacements.
- A failed structural update restores the accepted physical generation and
  constructs fresh factors rather than reusing an uncertain partial state.
- When a calculation starts from the conventional enormous energy sentinel,
  the first exact order-one energy retargets the fixed QR shift before order
  two is attempted. This avoids loss of the physical energy in
  `shift+(lambda-shift)`.

`EvaluateQAppendedTrial(BaseOrder,TargetOrder,...)` composes trim, suffix map,
staging, append, and solve. Keeping this lifecycle below the BBOP driver makes
the accepted/tentative transition explicit and reusable in sibling ports.

### Step 2: BasisEnlQ driver

`BasisEnlQ` follows `BasisEnlG` for random generation, best-candidate choice,
optional DRMNG refinement, overlap/coefficient acceptance, history, saving,
and retry limits. The differences are deliberate:

- H/S capacity is `Kstop`, while `Q_Workspace%MatrixOrder` is only the current
  accepted or tentative order.
- New functions always occupy a canonical suffix; existing functions are
  never permuted.
- Random candidates use delete/append order transactions. Nonlinear optimizer
  trials use ordinary indexed replacement transactions.
- If the configured number of candidate attempts is exhausted, Q returns a
  failure instead of accepting the final rejected candidate. This prevents a
  matrix/factor generation rejected by overlap, coefficients, or solver
  status from leaking into the saved basis.
- Q allocates neither `Glob_diagH` nor DSYGVX/GSEPIIS factor work arrays.

While auditing the reference driver, `BasisEnlG` was given an initialization
for `x_best` in `OptimizationType=0`; without it, that path could copy an
undefined candidate vector. No other G numerical path was changed.

### Step 3: indexed overlap penalty

`ComputeQOverlapPenalty`, `ComputeQOverlapPenaltyAndAddGradient`, and
`GetQOverlapStatistics` operate on the arbitrary canonical active map. A pair
is included exactly once when either endpoint is active. For normalized
overlap `S(i,j)`, the derivative contribution for an active endpoint uses the
existing scaled derivative convention:

```text
dS_normalized = D_cross - 0.5*S(i,j)*D_diagonal
```

An active--active pair adds the appropriate derivative to both optimizer
blocks but adds only one energy penalty. The gradient is accumulated before
the final MPI all-reduce in `EnergyQB`, so it follows the existing distributed
gradient ownership model.

### Step 4: FullOpt1Q driver

`FullOpt1Q` uses the ascending canonical map `InitFunc:FinalFunc`. This order is
also the x vector, gradient, DRMNG Hessian, and Hessian-file order; there is no
temporary suffix permutation to reverse during saving or restart.

The initial physical matrix is restored or assembled once and factorized
fresh. Each DRMNG trial is then staged and committed through the replacement
transaction. A full-basis trial still has cubic asymptotic work because it
contains O(n) symmetric column updates of O(n^2) each. This implementation
prioritizes one correct transaction path for arbitrary subsets; a future
benchmark may choose fresh factorization above a measured active-density
crossover without changing matrix layout or driver semantics.

Periodic basis saves, Hessian saves/restarts, overlap penalties, physical
energy reporting, history updates, and final accepted-point restoration match
the G driver. Output uses `Sort='no'` because Q never changes canonical order.

### Validation performed

The strict wp=8 debug build passed after each driver was added. Disposable
end-to-end tests included:

- basis growth from order zero to one, exercising exact order-one solution and
  shift retargeting;
- restart growth from one to three on one and two MPI ranks, exercising suffix
  append, rejection trimming, and active order smaller than capacity;
- a three-function `FULL_OPT1 Q` run optimizing interior canonical function 2
  on one and two ranks, with two energy and one gradient evaluation;
- indexed overlap penalty and gradient for noncontiguous active functions
  `(100,50)` in the 100-function Li basis. All twelve analytic entries matched
  central finite differences; the maximum relative discrepancy was
  `6.7315e-10`, and factor/eigenpair residuals were of order `1e-16`.

The append tests motivated the residual retry iteration floor described in
section 14. A qrlinalg direction-change timeout was not automatically treated
as success: the independent physical residual remained the acceptance test.

## 16. Remaining Q BBOP and output-method implementation log

### Step 1: SAVE_HSWF

`SaveHSWF` accepts `Q`. Matrix assembly leaves canonical lower H/S untouched.
When matrix text is requested, `QCanonicalMatrixElement` materializes each
conceptual symmetric value only at the output boundary, so no full upper
triangle is created in memory. When an eigenvector or wave function is
requested, the routine constructs fresh factors, calls `SolveQ`, and exports
the S-normalized canonical coefficient vector.

### Step 2: EXPC_VALS, DENSITIES, and MOMT_DENS

All three commands enter `ExpectationValues`; its Q branch performs one fresh
canonical factorization and solve before the existing operator kernels run.
The density commands therefore need no separate eigensolver implementation.
As everywhere else, Q returns the root closest to the shift, not an eigenvalue
selected by ordinal index.

The first complete debug execution exposed the historical malformed format
`'(a,i1,i1.1x)'` while writing `expvals.txt`. It was corrected to
`'(a,i1,i1,1x)'`. This is independent of Q but was necessary for every debug
expectation-value calculation to finish.

### Step 3: elimination and separation

`SolveEliminationGSEP` is the common local boundary for the four legacy
save-and-stop routines. G retains its DSYGVX upper-triangle materialization.
Q consumes the canonical lower triangles, constructs fresh qrlinalg factors,
and returns an S-normalized coefficient vector. The routines take an optional
method argument so all existing three-argument/four-argument G call sites keep
their source behavior; `main.f90` passes `Q` explicitly on the new paths.

After elimination or random separation, these routines deliberately retain
their historical full `ComputeMatElem(1,cbs)` reconstruction and then call the
common solver. This means their one-time linear-algebra step remains O(n^3),
like DSYGVX. It is not an iterative optimizer bottleneck. Principal deletion
for elimination and indexed replacement for separation are valid future
performance refinements, but require careful in-place compaction and separate
tests; they are not prerequisites for correct Q semantics.

For Q overlap reporting, the lower element `Glob_S(j,i)`, `j>i`, is used
directly. The former print statement read `Glob_S(i,j)` only because G had
materialized the upper triangle; leaving that read in Q would observe
unspecified storage.

The same audit corrected two structural edge cases shared with G. `ELIM_LCFN`
now stops without saving if its threshold would remove the entire basis,
because neither a generalized eigenproblem nor a qrlinalg order-zero state
exists. `ELIM_LND1` now sets the new order from the number of uniquely masked
functions, not the number of offending pairs; one function may participate in
several printed pairs but must be removed only once.

### Validation performed and replay still required

The wp=8 strict debug build compiles and links all new dispatch paths.
Disposable one-rank executions produced:

- `SAVE_HSWF Q 3 ...`: energy `-2.2487037512529300`, full H and S text,
  eigenvector, and wave-function files;
- `EXPC_VALS Q 3`: the same energy, `S=1.0000000000000002`, all expectation
  values, symmetrized values, and output file completed;
- `ELIM_LCFN Q 3 0.1 ...`: function 3 with coefficient
  `0.089552928023` was removed; the rebuilt order-two energy was
  `-2.2109294026962369`, and the reduced basis was saved;
- `ELIM_LND1 Q 3 0.2 ...`: pair `(2,3)` with overlap
  `0.32073776366828` was detected; the rebuilt order-two result had the same
  energy and was saved.

The elimination routines intentionally call `MPI_Abort` after saving, exactly
as their G implementations do, so exit status 1 is expected after the success
messages and output file.

The environment approval service became unavailable before the two final MPI
runs could start. After review, replay these disposable tests (the final exit
status is again expected to be 1):

```text
SEPR_LND1 Q 3 0.2 0.1 separated-lnd-q.txt
SEPR_FLCF Q 3 0.7 0.1 separated-flcf-q.txt
```

Optimized Netlib builds compiled and linked at wp=8, wp=10, and wp=16. The
build directory was cleaned between precisions so module/object reuse could
not mix working kinds. After the final elimination edge-case correction, the
strict wp=8 debug build was repeated successfully. No sample input or tracked
generated output was modified.
