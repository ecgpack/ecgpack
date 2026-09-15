# Q-method design and maintenance guide

## Purpose and status

This document is the production design reference for generalized-eigenvalue
method `Q` in ECGPACK. It describes the numerical contract, canonical matrix
layout, QR-update lifecycle, MPI ownership, supported BBOP commands,
basis-specific adaptations, performance characteristics, and the procedure
for porting the method to another ECGPACK variant.

Method `Q` is implemented in the four real-ECG energy codes:

- `RG_0S`;
- `RG_1P`;
- `RG_2D`; and
- `RG_2P`.

The implementations are intentionally local to each code because ECGPACK's
energy variants have separate source trees and basis-specific matrix-element
interfaces. `RG_0S` remains the reference implementation. The Q algorithms
are not currently implemented in `CG_0S` or the off-diagonal matrix-element
codes.

The following commands accept `Q` in the real-ECG energy codes:

| Command | Q behavior |
| --- | --- |
| `BASIS_ENL` | Appends candidate functions in canonical order and updates or rebuilds the QR factors. |
| `OPT_CYCLE` | Optimizes arbitrary canonical functions through symmetric column replacements. |
| `FULL_OPT1` | Optimizes an arbitrary contiguous function range without permuting the physical basis. |
| `ELIM_LCFN` | Deletes the coefficient-selected canonical functions and updates the QR state. |
| `ELIM_LND1` | Deletes the overlap-selected canonical functions and updates the QR state. |
| `SEPR_LND1` | Replaces perturbed canonical columns selected by the overlap rule. |
| `SEPR_FLCF` | Replaces perturbed canonical columns selected by the coefficient rule. |
| `EXPC_VALS` | Constructs and solves the Q eigenproblem before evaluating operators. |
| `DENSITIES` | Uses the Q eigenpair before computing position-space distributions. |
| `MOMT_DENS` | Uses the Q eigenpair in the RG_0S momentum-density path. |
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
A = H - sigma*S = Q*R,
```

where `sigma` is the current inverse-iteration shift. The shift remains fixed
while structural QR updates are reused. Any shift change requires a fresh
factorization.

Q is a shift-targeted inverse-iteration method. It follows method I's
root-near-shift semantics:

- `CURRENT_ENERGY`, `INVITPARAMETER`, and `EIGVAL_TOLERANCE` control the
  target and convergence;
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

### Order-one special case

For one basis function, `SolveQ` uses the exact result

```text
E = H(1,1)/S(1,1),
c(1) = 1/sqrt(S(1,1)).
```

This avoids unnecessary inverse iteration and prevents catastrophic
cancellation when a new calculation still carries the conventional enormous
`CURRENT_ENERGY` sentinel. Before growing beyond order one, basis enlargement
retargets the shift to the first accepted physical energy.

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
qrlinalg factors        the identical canonical row/column order
```

This stable-order rule is the central Q invariant. The G implementation may
move active functions into a trailing block for Cholesky-based work; that
permutation must not be copied into Q. If a user-visible operation genuinely
reorders a Q basis, all basis-indexed arrays must be permuted consistently and
the QR factors must be rebuilt.

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
`StoreQCanonicalColumn` centralize this rule. QR replacement operations need
a full conceptual column, so gathering or storing one symmetric column costs
O(N) while retaining one authoritative physical copy.

### Why H and S remain unshifted

`Glob_H` must never be overwritten with `H-sigma*S` in a Q calculation.
Keeping physical H and S separate is required for:

- changing the inverse-iteration shift;
- evaluating physical eigenpair residuals;
- saving H and S through `SAVE_HSWF`;
- expectation values and overlap penalties;
- reconstructing fresh factors after drift or failure; and
- comparing Q against G or I without undoing an in-place transformation.

Only the private qrlinalg factors represent the shifted matrix.

## 3. qrlinalg contract

`src/qrlinalg.f90` is the API authority. A `qr_real_state` owns dense Q and R
factors and its private factorization/iteration workspace. It does not own H,
S, the coefficient vector, ECG nonlinear parameters, or MPI state.

The production code uses these operations:

| Operation | Contract |
| --- | --- |
| `initialize(capacity,info)` | Allocates capacity-sized state but creates no valid factorization. Reinitialization discards previous factors. |
| `factorize_fresh(H,S,shift,info,active_order=N)` | Reads the lower triangle of the active leading block and constructs factors for `H-shift*S`. Physical leading dimensions may exceed N. |
| `replace_symmetric(index,delta_h,delta_s,info)` | Applies a symmetric row/column change at any one-based canonical index. Both delta vectors have conceptual length N. |
| `append_symmetric(h_column,s_column,info)` | Appends the last row/column to an existing nonempty state. The first function must be factorized fresh. |
| `delete_symmetric(index,info)` | Deletes an arbitrary canonical row/column and preserves survivor order. It cannot represent order zero. |
| `solve(S,v_initial,x,lambda,tol,max_iter,norm_mode,rel_acc,num_iter,info)` | Performs inverse iteration using distinct input/output vectors. ECGPACK passes `norm_mode=0` for S normalization. |
| `factorization_residual(H,S,v,absolute,relative,info)` | Checks the factor action against `(H-shift*S)*v` without changing state. This is not the physical eigenpair residual. |
| metadata queries | Report validity, order, capacity, shift, and structural updates since the last fresh factorization. |
| `clear()` | Releases the private state. |

The active-order and physical-leading-dimension API is performance-critical.
It allows `BASIS_ENL Q` to allocate matrices to final capacity while solving a
smaller active leading block without creating a packed O(N^2) temporary on
every solve.

The qrlinalg state is private in production. Do not define
`QRLINALG_TESTING` in an ECGPACK build.

### Update error boundary

The qrlinalg wrappers validate indices, dimensions, state, and capacity before
calling the underlying QR-update kernels. The kernels themselves do not
return a numerical status. A structurally successful update can still leave a
singular or ill-conditioned shifted problem; `solve` and the independent
residual checks detect that condition later.

`QR_ERR_NO_CONVERGENCE` may carry an approximate eigenpair. ECGPACK accepts it
only when the independently computed physical residual meets the Q tolerance.

## 4. Q workspace and ownership

`QMethodWorkspace` is defined in each real energy code's `workproc.f90`. It is
kept there because optimizer scheduling, trial acceptance, physical matrices,
and QR recovery form one lifecycle. Moving only some of this state into a
separate context module would obscure ownership without reducing coupling.

### MPI ownership

- `Glob_H`, `Glob_S`, `Glob_diagS`, coefficients, nonlinear parameters, basis
  metadata, and relationship flags are replicated as in existing ECGPACK
  calculations.
- Matrix-element work remains distributed with the established MPI routing.
- The serial `qr_real_state` and its dense Q/R arrays are meaningful only on
  rank zero.
- Every qrlinalg status and every returned eigenpair is broadcast before
  collective code continues.
- No rank may return early from a collective Q path on the basis of a local
  error that other ranks have not received.

### Relationship flags

The workspace tracks relationships that qrlinalg cannot infer from its private
validity flag:

```text
MatricesAreCanonical
FactorsMatchMatrices
MatrixParametersAreStored
TrialIsReady
TrialHasDerivatives
AcceptedPointIsStored
```

These flags are independent. For example, valid factors do not prove that a
caller has not edited `Glob_H`; and writing a new DRMNG point into
`Glob_NonlinParam` does not change which parameter generation the current
physical matrices represent.

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
descending, contiguous, or noncontiguous. The map defines the order of
optimizer variables and derivatives independently of physical matrix order.

## 5. Trial assembly and derivatives

`AssembleQTrial` stages complete normalized target columns before it modifies
either the physical matrices or the QR factors. For each unordered matrix
pair touching at least one active function, it:

1. assigns the pair through the existing MPI work distribution;
2. calls the local basis-qualified H/S matrix-element routine;
3. computes all active raw diagonal overlaps needed for normalization;
4. places the normalized physical element into every affected staged column;
5. places derivatives in the optimizer block associated with the active
   endpoint; and
6. combines distributed results before a transaction begins.

When both endpoints are active, the shared matrix element is evaluated once,
but derivatives with respect to both endpoint parameter blocks may be needed.
The matrix-element call must request `DActive` for the first active endpoint
and `DOther` only when the other endpoint is also active. Never assume a
trailing active block when indexing `Glob_D`.

Overlap penalties use the same active map. Each unordered pair touching an
active function is counted exactly once; its energy and gradient contribution
is assigned to the correct active block.

The matrix-element routine names and argument order are basis-specific. See
section 10 before copying this code to another variant.

## 6. Replacement transaction

`ApplyQTrial` changes an existing fixed-order basis atomically at the ECGPACK
level:

1. Validate method, order, capacity, active indices, staged columns, and every
   required allocation.
2. On rank zero verify that qrlinalg is valid and has the same order, capacity,
   and shift as the physical workspace.
3. Snapshot all original active columns, raw norms, and represented nonlinear
   parameters before the first update.
4. For each active function in active-map order, gather the progressively
   current physical column and compute `delta=target-current`.
5. Apply `replace_symmetric` on rank zero and broadcast its status.
6. Only after success, store the target column and diagonal norm on every
   rank.
7. After all columns succeed, mark the factor and physical generations as
   matched and update `MatrixParam`.

The progressively current gather in step 4 is essential. Once the first of
two active columns installs their shared intersection, the later column must
see that new value and apply zero delta at the intersection. Computing every
delta from the original matrix would update the same symmetric element twice.

If an unexpected update failure occurs, restore all physical snapshots and
represented parameters, invalidate the staged trial, and factorize fresh.
Never continue from a state whose factor/matrix relationship cannot be proven.

DRMNG evaluates many points that are not accepted. It is valid for physical
matrices and factors to follow the most recently evaluated point, because the
next trial can update from it. The optimizer's best point remains separate.
When accepting `x_best` or restoring `x_init`, reconstruct that complete point
through the same transaction; restoring parameters alone is invalid.

## 7. Growth and deletion transactions

### Basis growth

`BASIS_ENL Q` keeps the accepted prefix canonical and represents a candidate
block as a consecutive tentative suffix.

- From order zero, stage the complete first block, store it physically, and
  factorize fresh because qrlinalg has no valid order-zero factorization.
- With a nonempty accepted prefix, append staged functions one at a time with
  `append_symmetric`.
- Competing candidates at the same tentative order replace the previous
  suffix; they are not appended repeatedly.
- `Glob_CurrBasisSize` remains the accepted size while
  `Q_Workspace%MatrixOrder` may include the tentative suffix.
- Rejection trims suffix factors and exposes the untouched accepted leading
  block.

Discrete angular metadata is not part of `MatrixParam`. Every metadata-only
candidate change must therefore rebuild the candidate suffix with
`EvaluateQAppendedTrial`. After scanning choices, explicitly reconstruct the
winning metadata before continuous optimization. Parameter equality alone
cannot prove that H/S belongs to the winning discrete indices.

### Basis deletion

`DeleteQMaskedFunctions` deletes selected indices from the QR state in
descending canonical order, compacts the lower-triangular physical matrices,
and compacts every basis-indexed array in the same survivor order. Descending
factor deletion prevents earlier canonical indices from changing before they
are processed.

If an optimized delete path fails or its state relationship becomes
uncertain, rebuild the survivor matrices and factorize fresh. Code that can
delete every function must handle the order-zero case outside qrlinalg.

### Separation

The separation commands perturb existing nonlinear parameters without
changing canonical indices. Their Q path stages and replaces the selected
columns. Basis-specific discrete metadata remains attached to the same
functions.

## 8. Solve, residuals, and refresh policy

Two residuals monitor different numerical relationships and neither may be
removed:

### Physical generalized-eigenpair residual

```text
r_eig = H*c - E*S*c

relative = ||r_eig|| /
           (||H*c|| + |E|*||S*c|| + tiny)
```

`ComputeQEigenpairResidual` evaluates the two symmetric products with
`DSYMV('L',...)`, reading only the authoritative lower triangles. It is called
on rank zero and must use BLAS directly rather than an MPI-collective wrapper.

### Factorization residual

```text
r_fac = (H-sigma*S)*v - Q*(R*v)
```

qrlinalg evaluates this action residual with a deterministic dense probe. The
eigenvector is intentionally not used as the probe because it is nearly a
null vector of the shifted matrix and can make harmless factorization error
look disproportionately large.

### Current refresh thresholds

Before solving, factors are rebuilt when structural updates since the last
fresh factorization reach

```text
max(64,8*N).
```

After a solve, a fresh factorization and one retry are requested when either

```text
factor residual > max(1000*epsilon*N, 10*EIGVAL_TOLERANCE)
```

or

```text
eigenpair residual > max(1000*epsilon*N, 100*EIGVAL_TOLERANCE).
```

The retry iteration limit is

```text
max(120,default_limit,4*N).
```

The O(N) update-count guard keeps an O(N^3) refresh amortized to O(N^2) per
replacement. Residual monitoring can trigger an earlier refresh when matrix
conditioning or accumulated rotations require it.

A fresh factorization is also mandatory after initialization, swap restore,
shift change, factor/matrix metadata mismatch, or uncertain rollback.

## 9. BBOP driver behavior

### OPT_CYCLE Q

The driver preserves G's validation, scheduling, DRMNG settings, counters,
acceptance rules, history, and save frequency, but it never physically moves
the optimized functions.

The active map reproduces the optimizer-variable order expected by the legacy
G schedule. It must be verified for partial ranges, final short blocks,
overlapping blocks, shifts larger than one, and restart in the middle of a
cycle. Q does not require `SaveResults(Sort='yes')` to undo an internal
permutation because no such permutation occurs.

### FULL_OPT1 Q

The active map is the requested `InitFunc:FinalFunc` range, and Hessian I/O
uses that same variable order. A trial changing M functions requires O(M*N^2)
sequential symmetric replacements. When M is proportional to N, this remains
an O(N^3) operation; QR updates do not make full-basis changes quadratic.

### Output and observable commands

`SAVE_HSWF Q` exports canonical, normalized, unshifted physical H and S. It
does not export `H-sigma*S` or private Q/R storage. Expectation-value and
density commands build or restore canonical Q state, solve once, and then use
the existing operator paths with the S-normalized coefficient vector.

### Swap state

Swap files persist physical matrix information and enough metadata to restore
canonical H/S. Private Q/R factors are never treated as portable persistent
state; they are reconstructed after a restore. Any swap format change must be
tested for G and I as well as Q because the routines are shared.

## 10. Basis-specific adaptations

The canonical layout and qrlinalg lifecycle are common, but matrix-element
interfaces and discrete metadata are not:

| Code | Metadata for function `i` | H/S matrix-element routine |
| --- | --- | --- |
| `RG_0S` | none | `MatrixElementsHS_RG_0S` |
| `RG_1P` | `Glob_ZIndex(i)` | `MatrixElementsHS_RG_1P` |
| `RG_2D` | `Glob_Index(i,1:2)` | `MatrixElementsHS_RG_2D` |
| `RG_2P` | `Glob_Index(i,1:2)` | `MatrixElementsHS_RG_2P` |

Consult the actual routine declaration before adapting a call. In particular,
RG_2D and RG_2P use different argument ordering, and RG_2D returns additional
`Tkl` and `Vkl` values. The active function must be passed as the derivative
endpoint associated with `DActive`.

Discrete metadata must be generated, broadcast, copied during acceptance,
saved, restored, printed, optimized when supported, and compacted during
deletion:

- RG_1P carries one Z index;
- RG_2D carries two indices and may restrict their relationship with
  `VECTOR_COUPLING_SCHEME`; and
- RG_2P carries two distinct indices.

### RG_2D vector-coupling schemes

The canonical ecgpack RG_2D code supports:

```text
scheme 0: no relationship restriction
scheme 1: the two indices must differ
scheme 2: the two indices must be equal
```

`GenerateTrialParam` establishes the invariant for random candidates.
`BasisEnlQ` must preserve it during discrete optimization:

- scheme 1 excludes the other index from each one-index scan;
- scheme 2 varies both indices together as one discrete choice; and
- if either scheme-2 index is fixed, both are effectively fixed.

Do not replace RG_2D's newer `globvars.f90`, input/output support, candidate
generation, or G/I index scans with files from an older branch while porting Q.

## 11. Memory and performance model

For matrix order N and active block size M:

- physical H and S remain O(N^2) replicated storage, as in existing methods;
- rank zero additionally owns dense Q and R, approximately 2*N^2 values;
- staged current/previous active columns require O(M*N) storage;
- one symmetric replacement costs O(N^2);
- one inverse-iteration step costs O(N^2);
- both residual checks cost O(N^2);
- an M-function trial costs O(M*N^2), excluding matrix-element evaluation;
- fresh factorization costs O(N^3); and
- append/delete operations use QR structural updates, with fresh rebuild as
  the recovery path.

For `OPT_CYCLE` with fixed small M, the energy-solve/update portion is
quadratic per trial. `FULL_OPT1` becomes cubic when M scales with N. Matrix
element generation, MPI communication, optimizer behavior, and the selected
BLAS implementation can dominate the observed wall time.

OpenMP is process-local. Every MPI rank may create `OMP_NUM_THREADS` threads,
so approximate CPU demand is

```text
MPI ranks * OMP_NUM_THREADS.
```

Use explicit rank placement and thread affinity on production machines and
avoid unintentional oversubscription.

## 12. Build and third-party integration

The real-code Makefiles compile dependencies in this order:

```text
wp_def
  -> qrupdate_linalg
  -> qrupdate_error
  -> qrupdate_real, qrupdate_complex
  -> qrupdate
  -> qrlinalg
  -> workproc
```

`qrlinalg.f90` requires preprocessing: `-cpp` for gfortran, `-fpp` for
ifort/ifx, or `-Mpreprocess` for nvfortran. Fixed-form qrupdate BLAS/LAPACK
sources also require the compiler-specific unlimited-line-length option.

Only one BLAS/LAPACK provider may be linked. With `LINALG=netlib`, the build
uses the qrupdate directory's aggregate `BLAS.f` and `LAPACK.f`, which contain
the routines required by both ECGPACK and qrlinalg. Optimized precision-8
providers resolve both call sets externally. Precisions 10 and 16 require the
bundled generic provider.

Serial objects use `debug/` or `release/`. `OPENMP=yes` uses isolated
`debug-omp/` or `release-omp/` trees so incompatible objects cannot be reused.

Example direct build:

```bash
make -C RG_0S release COMPILER=gfortran MACHINE=linux-generic \
  PREC=8 LINALG=netlib OPENMP=yes EXEFILE=ecg
```

Example batch build:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release \
  code=RG_0S nparticles=4 precision=8 linalg=netlib openmp=yes
```

The bundled QR-update kernels come from Martin Köhler's `qrupdate-ng`:

<https://gitlab.mpi-magdeburg.mpg.de/koehlerm/qrupdate-ng>

The imported revision is
`3fa4f77d6259c00cf0c36632287063f1277a3fe0`. Preserve the source headers,
per-directory `COPYING` files, README provenance, and
`THIRD-PARTY-NOTICES.md` entry when updating or copying the implementation.

## 13. Source map

The primary production entry points are:

| File | Routines or responsibility |
| --- | --- |
| `src/matform.f90` | `StoreHS` selects canonical normalized Q storage. |
| `src/workproc.f90` | `PrepareQWorkspace`, `EnsureQActiveCapacity`, `ClearQWorkspace`. |
| `src/workproc.f90` | `SetQActiveFunctions`, `CaptureQMatrixParameters`. |
| `src/workproc.f90` | `QCanonicalMatrixElement`, `GatherQCanonicalColumn`, `StoreQCanonicalColumn`. |
| `src/workproc.f90` | `AssembleQTrial`, `ApplyQTrial`, `AppendQTrial`, `EvaluateQAppendedTrial`, `TrimQFactors`. |
| `src/workproc.f90` | `FactorizeQFresh`, `SolveQ`, `ComputeQEigenpairResidual`. |
| `src/workproc.f90` | `ComputeQOverlapPenalty`, `ComputeQOverlapPenaltyAndAddGradient`, `GetQOverlapStatistics`. |
| `src/workproc.f90` | `EnergyQA`, `EnergyQAM`, `EnergyQB`. |
| `src/workproc.f90` | `BasisEnlQ`, `OptCycleQ`, `FullOpt1Q`. |
| `src/workproc.f90` | `DeleteQMaskedFunctions` and Q branches in cleanup/output routines. |
| `src/main.f90` | BBOP dispatch for method Q. |
| `src/qrlinalg.f90` | QR state, structural updates, inverse iteration, and factor residual. |
| `src/qrupdate/` | Generic real/complex QR-update kernels and required Netlib aggregates. |
| `Makefile` | Dependency order, preprocessing, provider selection, and OpenMP build isolation. |

Line numbers are intentionally omitted. Search by routine name because the
four duplicated workproc files evolve independently.

## 14. Validation requirements

Numerical changes to Q are incomplete until they pass tests that exercise both
the physical eigenproblem and the update lifecycle. Use disposable copies of
sample inputs; energy-code runs rewrite `inout.txt`.

### Required build coverage

- strict debug build at precision 8;
- release build at precision 8;
- bundled `LINALG=netlib` path;
- at least one optimized precision-8 BLAS/LAPACK provider when available;
- OpenMP build for changes touching compiler flags or qrupdate; and
- precisions 10 and 16 when changing generic numerical interfaces.

Fortran modules in these legacy Makefiles are safest to build serially unless
all module dependencies have been verified. Do not assume `make -j` is valid.

### Required runtime coverage

- order one, order two, active order smaller than capacity, and growth to
  capacity;
- active first, middle, last, contiguous, and noncontiguous functions;
- two active functions sharing a matrix intersection;
- repeated rejected trials followed by reconstruction of `x_best` and
  restoration of `x_init`;
- append, arbitrary deletion, and stable survivor order;
- singular or nearly singular shifted systems and nonconvergence status;
- physical residual and S normalization;
- fresh-factor versus updated-factor agreement;
- G/Q or I/Q energy and gradient agreement when targeting the same root;
- finite-difference gradient and overlap-penalty checks;
- swap restart and chained Q BBOP steps;
- all four cleanup commands;
- one and multiple MPI ranks; and
- representative basis sizes for memory and timing.

For indexed bases, also verify metadata after generation, optimization,
deletion, separation, save, and restart. RG_2D tests must cover coupling
schemes 1 and 2; RG_2P tests must confirm distinct indices.

### Canonical-port baseline

The canonical ecgpack port was validated with strict-debug and release/OpenMP
builds of all four real variants. Disposable three-function bases constructed
from each variant's own sample header gave fresh I/Q energy agreement within
approximately 1E-15. Separately assembled overlap matrices were byte-identical;
Hamiltonian entries agreed within approximately 1.1E-14. RG_2D scheme-1
functions retained unequal indices and scheme-2 functions retained equal
indices through optimization type 1.

These stochastic energies are not reference values. Future validation should
compare methods on the exact same saved physical basis and report tolerances,
not expect a separately generated random basis to reproduce an old energy.

## 15. Porting procedure

Use the following sequence when adding Q to another ECGPACK variant or a
newer divergent branch:

1. Identify a common pre-Q baseline and compare the target against it.
2. Inventory target-only changes in `main.f90`, `matform.f90`, `workproc.f90`,
   `globvars.f90`, the Makefile, and input/output documentation.
3. Study the target's basis-qualified matrix-element interface, normalization,
   derivative endpoint order, and discrete metadata.
4. Transfer `qrlinalg.f90` and `qrupdate/` unchanged unless the numerical
   library itself is deliberately being upgraded.
5. Add Makefile dependency order, preprocessing, one-provider linking, and
   isolated OpenMP build directories.
6. Add Q storage to `StoreHS` without changing G/I storage.
7. Port the workspace, canonical column helpers, arbitrary-index assembly,
   factor lifecycle, solve, residuals, and overlap helpers.
8. Port drivers one at a time: `OPT_CYCLE`, `BASIS_ENL`, `FULL_OPT1`, then
   output/observable and cleanup commands.
9. Adapt generation, acceptance, printing, save/restore, and compaction for
   every basis-specific metadata array.
10. Merge target-only semantics explicitly. A conflict-free patch does not
    prove that a newly added Q path enforces them.
11. Compile after each reviewable driver and run a fixed-basis comparison
    before stochastic optimization tests.
12. Update this guide, public input/build manuals, agent guidance, and
    third-party provenance when interfaces or behavior change.

Do not copy an entire final `workproc.f90` over a sibling. Use a three-way
merge or targeted patches so basis-specific physics, file formats, and newer
local semantics remain intact. Treat `CG_0S` as a separate numerical design:
Hermitian storage, conjugation, complex QR updates, and complex gradients need
their own derivation and validation.

## 16. Production invariant checklist

Before approving a Q change, confirm all of the following:

- physical basis order is canonical and stable;
- lower H/S triangles are normalized and unshifted;
- upper H/S triangles are never treated as authoritative;
- qrlinalg order, capacity, shift, and basis order match physical matrices;
- matrix parameters describe the generation currently stored in H/S;
- staged trials are complete before the first QR update;
- active-active intersections are changed exactly once;
- every factor status and eigenpair is synchronized across MPI ranks;
- failed or uncertain transactions recover through proven physical state and
  fresh factorization;
- physical and factorization residuals are both monitored;
- saved H/S and wave functions remain in canonical user-visible order;
- basis-specific metadata follows every accepted, deleted, or saved function;
  and
- performance claims distinguish fixed-small-block updates from full-basis
  optimization.

If any relationship cannot be demonstrated from current state and flags,
invalidate the factors and rebuild. A fresh O(N^3) factorization is preferable
to reusing an unproven O(N^2) update state.
