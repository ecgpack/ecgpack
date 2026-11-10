# Quadruple-precision modified Bessel functions

Made with the assistance of Claude Code.

## Purpose and scope

This note documents the working-kind adaptation of the SLATEC routines for
computing the modified Bessel function of the second kind,

\[
K_\nu(x), \qquad x>0,\quad \nu\geq 0,
\]

and its exponentially scaled form

\[
\widehat K_\nu(x)=e^xK_\nu(x).
\]

The implementation is the module `bessel` in
[`RG_0S/src/bessel.f90`](../../../RG_0S/src/bessel.f90), with its coefficient
tables in the module `bessel_coeff` in
[`RG_0S/src/bessel_coeff.f90`](../../../RG_0S/src/bessel_coeff.f90). Identical
copies of both files are compiled into `RG_1P`, `RG_2D`, and `RG_2P`. The
public entry point is `DBESK(X,FNU,KODE,N,Y,NZ[,F1MinusOne])`. `KODE=1` returns the unscaled function and
`KODE=2` returns the scaled function. The output sequence is

\[
Y_i=K_{\mathrm{FNU}+i-1}(X)
\]

or its scaled counterpart. Scaling is essential at large `X`, where
\(K_\nu(X)\) can underflow even though \(e^XK_\nu(X)\) remains representable.

The implementation is based on the SLATEC 4.1 routines `DBESK`, `DBSKNU`, and
`DASYIK` (author D. E. Amos, SNLA; public domain), available from
<https://netlib.org/slatec/>; provenance and the list of modifications are in
the repository's `THIRD-PARTY-NOTICES.md`. The original fixed-form sources are
not bundled with ECGPACK. A direct working-kind port of them retains the
approximately 14-decimal-digit accuracy limit stated in the original `DBSKNU`
prologue. The ECGPACK version extends the coefficient data and revises the
stopping and region-selection criteria so that extended and quadruple precision
can be used meaningfully.

The implementation is intended for the ECGPACK precision choices:

| `wp` | Typical representation with gfortran | Decimal precision |
|---:|---|---:|
| 8 | IEEE binary64 | 15 digits |
| 10 | x87 extended precision | 18 digits |
| 16 | IEEE binary128 | 33 digits |

Compiler support for these kind numbers is platform-dependent. In particular,
`wp=10` is primarily a GNU/x86 choice. Inside ECGPACK the kind is the `wp`
of the build (`PREC`), taken from `src/wp_def_<PREC>.f90` through `globvars`.

## Files and call structure

The relevant files are:

- [`RG_0S/src/bessel_coeff.f90`](../../../RG_0S/src/bessel_coeff.f90)
  (module `bessel_coeff`): regenerated Debye and reciprocal-gamma coefficient
  tables;
- [`RG_0S/src/bessel.f90`](../../../RG_0S/src/bessel.f90) (module `bessel`):
  - `DBESK`: public sequence driver, range checks, large-order routing, and
    forward order recurrence;
  - `DBSKNU`: small-order starting values from a power series, Miller
    algorithm, or large-argument expansion;
  - `DASYIK`: uniform large-order asymptotic expansion;
- [`test_bessel_k01_wp.f90`](test_bessel_k01_wp.f90): comparison program for
  orders zero and one, built against the module `bessel`;
- [`generate_bessel_k01_mpmath_reference.py`](generate_bessel_k01_mpmath_reference.py):
  generator for independent mpmath reference values; and
- [`bessel_k01_mpmath_reference.dat`](bessel_k01_mpmath_reference.dat):
  committed reference table.

The numerical routing is approximately:

```text
DBESK
  ├─ sufficiently large order ───────────────> DASYIK
  └─ reduced order -1/2 <= dnu < 1/2 ───────> DBSKNU
       ├─ x <= 2 ────────────────────────────> convergent power series
       ├─ 2 < x <= X2(wp) ──────────────────> Miller backward algorithm
       └─ x > X2(wp) ───────────────────────> large-x asymptotic series

Higher orders are formed by forward recurrence.
```

For \(K_\nu\), forward recurrence is numerically appropriate because the
desired solution is the dominant solution in the increasing-order direction.
In the notation used here,

\[
K_{\nu+1}(x)=K_{\nu-1}(x)+\frac{2\nu}{x}K_\nu(x).
\]

This is the recurrence implemented by the updates in `DBESK` and `DBSKNU`.

## Why changing the real kind was insufficient

The original double-precision routine contains several accuracy limits that do
not disappear when every variable is changed to `real(wp)`:

1. the eight-entry reciprocal-gamma table contains roughly 14 significant
   digits;
2. the 65-entry Debye table contains only \(U_1\) through \(U_{10}\), again
   rounded to roughly 14 digits;
3. the stopping tolerance is bounded below by \(10^{-15}\);
4. the large-argument expansion is limited to 30 terms and begins at the
   double-precision boundary `X2=17`;
5. the Miller workspace contains only 160 entries; and
6. the large-order routing thresholds 35 and 70 were selected for the
   original double-precision error target.

Promoting those values to `real(wp)` preserves already-rounded numbers; it
does not reconstruct the discarded digits. All of these restrictions therefore
had to be addressed separately.

## Uniform large-order expansion and the Debye table

### Published expansion

For positive real \(z\), define

\[
p=(1+z^2)^{-1/2}
\]

and

\[
\eta=\sqrt{1+z^2}+ln\!\left(\frac{z}{1+\sqrt{1+z^2}}\right).
\]

The large-order uniform expansion used by `DASYIK` is

\[
K_\nu(\nu z)
\sim
\sqrt{\frac{\pi}{2\nu}}
\frac{e^{-\nu\eta}}{(1+z^2)^{1/4}}
\sum_{k=0}^{\infty}(-1)^k\frac{U_k(p)}{\nu^k}.
\]

These are equations 10.41.4, 10.41.7, and 10.41.8 of the
[NIST Digital Library of Mathematical Functions, §10.41](https://dlmf.nist.gov/10.41).

The Debye polynomials start with \(U_0(p)=1\) and obey

\[
U_{k+1}(p)=
\frac{1}{2}p^2(1-p^2)U'_k(p)
+\frac{1}{8}\int_0^p(1-5t^2)U_k(t)\,dt.
\]

This is [DLMF equation 10.41.9](https://dlmf.nist.gov/10.41.E9).
The first polynomials provide useful checks on a coefficient generator:

\[
U_1(p)=\frac{3p-5p^3}{24},
\]

\[
U_2(p)=\frac{81p^2-462p^4+385p^6}{1152},
\]

\[
U_3(p)=
\frac{30375p^3-369603p^5+765765p^7-425425p^9}{414720}.
\]

They agree with [DLMF equation 10.41.10](https://dlmf.nist.gov/10.41.E10)
and with the first entries of the original SLATEC table.

### How the extended table was calculated

The recurrence above contains only differentiation, integration, addition,
multiplication, and rational constants. Starting from \(U_0=1\), every
coefficient is therefore rational. The extended table was generated as
follows:

1. represent each polynomial with exact rational coefficients;
2. apply the DLMF recurrence symbolically;
3. expand and collect powers of \(p\);
4. repeat through \(U_{24}\);
5. convert the exact rational values to decimal using 120-digit arithmetic;
6. write 50 significant digits for each Fortran literal and append `_wp`.

In pseudocode, the generation step is

```python
U = [1]
for k in range(24):
    next_U = (
        Rational(1, 2) * p**2 * (1 - p**2) * diff(U[k], p)
        + Rational(1, 8)
        * integrate((1 - 5*t**2) * U[k].subs(p, t), (t, 0, p))
    )
    U.append(expand(next_U))
```

Exact rational arithmetic is important here. Computing the recurrence in
binary64 and merely printing additional decimal places would preserve
binary64 rounding error rather than produce quadruple-precision coefficients.

Each \(U_k\) contains \(k+1\) powers,

\[
U_k(p)=c_{k,0}p^{3k}+c_{k,1}p^{3k-2}+\cdots+c_{k,k}p^k.
\]

The coefficients are stored from the highest power to the lowest power. The
flattened table therefore contains

\[
\sum_{k=1}^{24}(k+1)=324
\]

entries. The first entry for order \(k\) is at the one-based position

\[
L_k=1+\frac{(k-1)(k+2)}{2},
\]

which is the expression used in `DASYIK`.

### How `DASYIK` evaluates the table

With `RA=sqrt(1+z*z)`, the code has \(p=1/\mathrm{RA}\). It evaluates the
polynomial in \(p^2\) by Horner's rule. For the \(K\) function,

```fortran
T = -p/FN
AP = AP*T
AK = AP*S1
```

so `AP` supplies \((-p/\nu)^k\), while `S1` supplies the remaining polynomial
in \(p^2\). Their product is exactly

\[
(-1)^k\frac{U_k(p)}{\nu^k}.
\]

The leading factor in the code is also the DLMF leading factor because

\[
\sqrt{|T|}\sqrt{\frac{\pi}{2}}
=
\sqrt{\frac{\pi p}{2\nu}}
=
\sqrt{\frac{\pi}{2\nu}}(1+z^2)^{-1/4}.
\]

The summation now stops relative to `epsilon(1.0_wp)` rather than a fixed
\(10^{-15}\) floor.

## Reciprocal-gamma table for the small-order power series

### Why this table is needed

`DBSKNU` reduces the initial order to

\[
-\frac12\leq\mu<\frac12.
\]

For small \(x\), its power series needs

\[
G_1(\mu)=
\frac{\Gamma(1-\mu)^{-1}-\Gamma(1+\mu)^{-1}}{2\mu}.
\]

Direct evaluation loses digits near \(\mu=0\), because its numerator is the
difference of two nearly equal values. The code therefore uses a Taylor series
when \(|\mu|\leq0.1\).

Define

\[
F(z)=\frac{1}{\Gamma(1+z)}=\sum_{n=0}^{\infty}a_nz^n.
\]

Then

\[
\frac{F(-\mu)-F(\mu)}{2\mu}
=-(a_1+a_3\mu^2+a_5\mu^4+\cdots),
\]

and hence

\[
G_1(\mu)=-(a_1+a_3\mu^2+a_5\mu^4+\cdots).
\]

This explains the code sequence

```fortran
S = Bessel_RecipGammaC(1)
...
S = S+Bessel_RecipGammaC(K)*DNU2**(K-1)
G1 = -S
```

and the table mapping

\[
\texttt{Bessel\_RecipGammaC}(j)=a_{2j-1}.
\]

For example, its first entries are

\[
a_1=\gamma=0.5772156649015328606\ldots,
\]

\[
a_3=-0.0420026350340952355\ldots,
\qquad
a_5=-0.0421977345555443367\ldots.
\]

### How the reciprocal-gamma coefficients were calculated

The starting point is the reciprocal-gamma Maclaurin series in
[DLMF equation 5.7.1](https://dlmf.nist.gov/5.7.E1) and its coefficient
recurrence in [DLMF equation 5.7.2](https://dlmf.nist.gov/5.7.E2). An
equivalent and convenient generation formula follows from the logarithmic
gamma series in [DLMF equation 5.7.3](https://dlmf.nist.gov/5.7.E3):

\[
\ln F(z)
=
\gamma z-
\sum_{n=2}^{\infty}\frac{(-1)^n\zeta(n)}{n}z^n.
\]

Write

\[
\ln F(z)=\sum_{n=1}^{\infty}b_nz^n,
\]

where

\[
b_1=\gamma,
\qquad
b_n=-\frac{(-1)^n\zeta(n)}{n},\quad n\geq2.
\]

Since \(F'(z)=F(z)(\ln F(z))'\), equating coefficients gives

\[
a_0=1,
\qquad
a_n=\frac1n\sum_{k=1}^{n}k\,b_k\,a_{n-k}.
\]

This recurrence was evaluated with mpmath at 120 decimal digits. The 20 odd
coefficients

\[
a_1,a_3,\ldots,a_{39}
\]

were emitted with 50 significant digits. At the branch endpoint
\(|\mu|=0.1\), successive terms receive factors \(10^{-2j}\), so the extended
table provides ample margin beyond the approximately 33 decimal digits of
binary128.

## Large-argument expansion

For large positive \(x\), `DBSKNU` uses

\[
K_\nu(x)
\sim
\sqrt{\frac{\pi}{2x}}e^{-x}
\sum_{k=0}^{\infty}\frac{a_k(\nu)}{x^k},
\]

from [DLMF equation 10.40.2](https://dlmf.nist.gov/10.40.E2), with

\[
a_k(\nu)=
\frac{\displaystyle\prod_{m=1}^{k}
\left(4\nu^2-(2m-1)^2\right)}{k!\,8^k}.
\]

Consequently the terms obey

\[
t_k=t_{k-1}
\frac{4\nu^2-(2k-1)^2}{8kx}.
\]

This is the recurrence implemented by `CK`, `FMU`, `SQK`, and `DK` in
`DBSKNU`. No stored coefficient table is required for this expansion.

The series is asymptotic rather than convergent. The revised loop therefore:

- permits up to 200 terms instead of 30;
- stops when the next term is small relative to the accumulated sum; and
- stops before adding a term whose magnitude has begun to increase.

The original fixed boundary `X2=17` was suitable for a roughly
double-precision target. The revised boundary is

\[
X_2=\max\!\left(17,-\frac12\ln\epsilon_{\mathrm{wp}}+2\right).
\]

The motivation is that the best attainable remainder near optimal truncation
has an exponential scale comparable to \(e^{-2x}\). The additional `+2` is a
guard margin. This boundary is an implementation choice validated numerically,
not a formula quoted directly from DLMF. It gives approximately:

| Kind | `X2` |
|---|---:|
| binary64 | 20.0 |
| x87 extended | 23.8 |
| binary128 | 40.8 |

Moving `X2` outward means the intermediate Miller algorithm must cover a
larger region and meet a smaller tolerance. Its work arrays were therefore
increased from 160 to 512 elements.

## Large-order routing thresholds

The original `DASYIK` documentation specifies order 35 or greater, and the
SLATEC driver uses thresholds 35 and 70 depending on whether one or two
starting values are required. Those values target double precision with the
original ten Debye polynomials.

The revised driver retains `[35,70]` when the working kind has at most 15
decimal digits. For wider kinds it uses

\[
[N_1,N_2]=[3P,6P],
\]

where \(P=\texttt{precision(1.0\_wp)}\). This gives:

| `wp` | Decimal precision `P` | `NULIM` |
|---:|---:|---:|
| 8 | 15 | `[35,70]` |
| 10 | 18 | `[54,108]` |
| 16 | 33 | `[99,198]` |

These conservative thresholds keep the uniform asymptotic expansion in a
region where 24 Debye terms are sufficient for the requested tolerance. Like
the revised `X2`, they are precision-routing choices confirmed by numerical
comparison rather than universal mathematical constants.

## Other working-kind changes

The following changes are also required for genuine wider-precision
calculation:

- every floating-point variable and literal uses `real(wp)` and `_wp`;
- `epsilon(1.0_wp)`, `tiny(1.0_wp)`, and `huge(1.0_wp)` determine limits from
  the selected kind;
- the generic `gamma` intrinsic replaces the double-specific SLATEC
  dependency `DGAMMA`;
- `Glob_Pi` carries enough digits for binary128;
- the hard lower bound `max(D1MACH(3),1.0e-15)` is removed;
- SLATEC `XERMSG` calls are replaced by descriptive `error stop` statements;
  and
- integer-order starting values use `DBSKNU`, avoiding the separate
  double-precision `DBESK0`/`DBESK1` Chebyshev dependency chain.

The machine-constant routines in `RG_0S/src/X1MACH.f90` use Fortran inquiry
intrinsics, so their values follow the active `wp` rather than retaining
hard-coded binary64 constants.

## Verification against mpmath

The Python generator evaluates

```python
k0 = mp.besselk(0, x)
k1 = mp.besselk(1, x)
scaled_k0 = mp.exp(x) * k0
scaled_k1 = mp.exp(x) * k1
```

with a default `mp.dps` of 120 and writes 90 significant digits. Decimal input
strings are converted directly with `mp.mpf`, avoiding an intermediate
binary64 conversion.

The Fortran test reads the reference table and tests both `KODE=1` and
`KODE=2`. The present table spans \(10^{-40}\leq x\leq10^{32}\). When an
unscaled reference value is below the exponent range of the selected kind, a
zero or subnormal Fortran result is accepted; the scaled result is still
compared by relative error. The acceptance tolerance is

\[
1000\,\epsilon_{\mathrm{wp}}.
\]

The clean validation run produced:

| Kind | Comparisons | Maximum measured relative error |
|---|---:|---:|
| `wp=8` | 80 | \(5.67\times10^{-16}\) |
| `wp=10` | 80 | \(2.10\times10^{-19}\) |
| `wp=16` | 80 | \(1.12\times10^{-33}\) |

These tests cover orders zero and one across the power-series, Miller, and
large-argument regions. Additional noninteger-order and large-order checks were
also performed during development, including the order recurrence. The tests
provide strong regression evidence, but they are not a formal error proof over
the entire allowed \((x,\nu)\) domain.

## Building and running the standalone test

The test program is linked against the same sources that ECGPACK compiles.
`globvars` needs the C preprocessor and MPI modules, so use the MPI compiler
wrapper with `-cpp`. Choose the precision with `PREC` (`8`, `10`, or `16`).
From the repository root:

```bash
PREC=16
build_dir=/tmp/ecgpack_bessel_p$PREC
mkdir -p "$build_dir"

for f in RG_0S/src/wp_def_$PREC.f90 RG_0S/src/X1MACH.f90 \
         RG_0S/src/globvars.f90 RG_0S/src/bessel_coeff.f90 \
         RG_0S/src/bessel.f90 \
         doc/devnotes/mod_Bessel_high_precision/test_bessel_k01_wp.f90; do
  mpif90 -cpp -O2 -J"$build_dir" -I"$build_dir" \
    -c "$f" -o "$build_dir/$(basename "$f" .f90).o"
done

mpif90 "$build_dir"/{wp_def_$PREC,X1MACH,globvars,bessel_coeff,bessel,test_bessel_k01_wp}.o \
  -o "$build_dir/test_bessel_k01_wp"
```

Run the test with an explicit reference path:

```bash
"$build_dir/test_bessel_k01_wp" \
  doc/devnotes/mod_Bessel_high_precision/bessel_k01_mpmath_reference.dat
```

With gfortran this reports 80 passed, 0 failed for `PREC=8`, `10`, and `16`.
Generated objects, modules, and executables should remain outside tracked
source directories.

## Sources

1. NIST Digital Library of Mathematical Functions,
   [§10.41, Asymptotic Expansions for Large Order](https://dlmf.nist.gov/10.41),
   especially equations 10.41.4 and 10.41.7–10.41.10.
2. NIST Digital Library of Mathematical Functions,
   [§10.40, Asymptotic Expansions for Large Argument](https://dlmf.nist.gov/10.40),
   especially equation 10.40.2.
3. NIST Digital Library of Mathematical Functions,
   [§5.7, Gamma Function Series Expansions](https://dlmf.nist.gov/5.7),
   especially equations 5.7.1–5.7.3.
4. NIST Digital Library of Mathematical Functions,
   [§10.29, Modified-Bessel Recurrence Relations](https://dlmf.nist.gov/10.29),
   especially equation 10.29.1.
5. N. M. Temme, “On the numerical evaluation of the modified Bessel function
   of the third kind,” *Journal of Computational Physics* **19** (3), 324–337
   (1975), [doi:10.1016/0021-9991(75)90082-0](https://doi.org/10.1016/0021-9991(75)90082-0).
6. [SLATEC Common Mathematical Library, version 4.1](https://netlib.org/slatec/)
   and its [individual source directory](https://netlib.org/slatec/src/).

The formulas above come from DLMF and the algorithmic structure comes from
SLATEC and its cited Temme paper. The higher-order table sizes, guard digits,
precision-dependent routing thresholds, enlarged workspace, and test
tolerances are ECGPACK adaptation choices documented here so they can be
reproduced and revised independently of the historical binary64 constants.
