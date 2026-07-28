# Third-party notices

ECGPACK bundles some amount of source code written by others. Those
files remain under the terms of their original authors, which are
reproduced below. The terms in this file govern the listed files only;
all other files in this repository are covered by the ECGPACK license.

---

## Reference BLAS and LAPACK (University of Tennessee et al.)

**Files:** `BLAS.f` and `LAPACK.f` in each of `CG_0S/src/`, `RG_0S/src/`,
`RG_1P/src/`, `RG_2D/src/` and `RG_2P/src/` (five identical copies of
each file).

**Origin:** the netlib reference implementations of BLAS and LAPACK,
<https://www.netlib.org/lapack/>. The routine headers retained in the
sources identify them as version 3.4.0, with one Level 1 BLAS routine
carrying a version 3.1 header.

**Scope:** these are subsets, not the complete libraries. `BLAS.f`
contains 42 routines — the real double precision (`D`) and complex
(`Z`) Level 1, 2 and 3 routines used by the codes, plus `LSAME`,
`XERBLA`, `IDAMAX` and `DCABS1`. `LAPACK.f` contains 87 routines
supporting the symmetric and Hermitian generalized eigenvalue problem:
the `DSYEVX`/`DSYGVX`/`ZHEEVX`/`ZHEGVX` drivers, Cholesky factorization
(`DPOTRF`/`ZPOTRF`), reduction to tridiagonal form (`DSYTRD`/`ZHETRD`),
bisection and inverse iteration (`DSTEBZ`, `DSTEIN`), the QL/QR
iterations (`DSTEQR`, `DSTERF`) and their supporting auxiliaries.

These two files are compiled only when the build selects
`LINALG=netlib`. When an optimized library is chosen instead (`mkl`,
`openblas`, `lblas`, `aocl`) they are excluded from the build and the
external library is linked in their place.

**Modifications made for ECGPACK:** as with the other bundled sources,
the real kind has been made selectable. `DOUBLE PRECISION` declarations
were replaced by `REAL(wp)` (and the corresponding complex declarations
by `COMPLEX(wp)`), and the routines obtain the kind parameter `wp` from
module `wp_def`, so that the codes can be built in double (fp64),
extended (fp80) or quadruple (fp128) precision. The individual netlib
source files were also concatenated into the two aggregate files
`BLAS.f` and `LAPACK.f`, with separator comment lines inserted between
routines. Problems observed in these copies should be reported to the
ECGPACK authors rather than to the LAPACK developers.

### Copyright and license (BSD 3-Clause)

Reproduced verbatim from <https://www.netlib.org/lapack/LICENSE.txt>.
The `$COPYRIGHT$` and `$HEADER$` markers are present in the upstream
file and are retained here so that the text is unaltered.

```text
Copyright (c) 1992-2013 The University of Tennessee and The University
                        of Tennessee Research Foundation.  All rights
                        reserved.
Copyright (c) 2000-2013 The University of California Berkeley. All
                        rights reserved.
Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
                        reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

The copyright holders provide no reassurances that the source code
provided does not infringe any patent, copyright, or any other
intellectual property rights of third parties.  The copyright holders
disclaim any liability to any recipient for claims brought against
recipient by any third party for infringement of that parties
intellectual property rights.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

## PORT optimization routines (David M. Gay, Bell Laboratories)

**Files:** `CG_0S/src/dmng.f`, `RG_0S/src/dmng.f`, `RG_1P/src/dmng.f`,
`RG_2D/src/dmng.f`, `RG_2P/src/dmng.f` (five identical copies)

**Origin:** the unconstrained minimization routines of the PORT
Mathematical Subroutine Library, written by David M. Gay at Bell
Laboratories and obtained from the publicly available subset
distributed by Netlib, <https://www.netlib.org/port/>.

The Netlib PORT `readme` distinguishes a publicly available subset —
the PORT utility routines together with the NL2SOL variants and the
`mnf`/`mng`/`mnh` minimizers and their bounded variants — from the
complete PORT library, which is the thing covered by the separate
non-exclusive limited-use agreement described there. The routines
bundled here (`DMNG`/`DRMNG` and their supporting routines, the double
precision `mng`) belong to the publicly available subset.

**Algorithms documented in:**

1. Dennis, J. E., Jr., Gay, D. M., and Welsch, R. E. (1981),
   "An adaptive nonlinear least-squares algorithm",
   *ACM Trans. Math. Software* **7**(3), 348–368.
2. Gay, D. M. (1983), "Algorithm 611: Subroutines for unconstrained
   minimization using a model/trust-region approach",
   *ACM Trans. Math. Software* **9**(4), 503–524.

**Modifications made for ECGPACK:** the real kind has been made
selectable — `DOUBLE PRECISION` declarations were replaced by
`REAL(wp)` and each routine obtains the kind parameter `wp` from module
`wp_def`, so that the minimizer can be built in double (fp64), extended
(fp80), or quadruple (fp128) precision. The numerical algorithms are
otherwise unchanged. Problems observed in this copy should be reported
to the ECGPACK authors rather than to Netlib or the original author.

### Copyright and permission notice

```text
All software in the directory is Copyright Lucent Technologies.
Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this
entire notice is included in all copies of any software which is
or includes a copy or modification of this software and in all
copies of the supporting documentation for such software.

THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR
IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR LUCENT
TECHNOLOGIES MAKE ANY REPRESENTATION OR WARRANTY OF ANY KIND
CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS
FOR ANY PARTICULAR PURPOSE.
```

---

## SLATEC machine constants (`I1MACH` and `D1MACH`)

**Files:** `X1MACH.f90` in each of `CG_0S/src/`, `RG_0S/src/`,
`RG_1P/src/`, `RG_2D/src/` and `RG_2P/src/` (five identical copies).

**Origin:** obtained from Netlib, <https://www.netlib.org/slatec/>. The
file provides the machine constant enquiry functions `I1MACH` and
`D1MACH` of the SLATEC Common Mathematical Library, which are used by
the bundled PORT minimizer. They originate in the PORT framework of
Fox, Hall and Schryer at Bell Laboratories.

The copy bundled here is a Fortran 90 adaptation rather than the
original SLATEC source: the function bodies determine the machine
constants from Fortran 90 intrinsic enquiry functions (`RADIX`,
`DIGITS`, `MINEXPONENT`, `MAXEXPONENT`, `HUGE`, `TINY`, `EPSILON`,
`BIT_SIZE`) instead of the per-architecture `DATA` tables used by the
original fixed-form routines. The SLATEC prologue documentation has
been retained, and records that adaptation:

```text
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)
```

**Reference:** P. A. Fox, A. D. Hall and N. L. Schryer, "Framework for
a portable library", *ACM Trans. Math. Software* **4**(2) (June 1978),
177–188.

**Modifications made for ECGPACK:** as with the other bundled sources,
the real kind has been made selectable. `DOUBLE PRECISION` declarations
were replaced by `REAL(wp)` and the functions obtain the kind parameter
`wp` from module `wp_def`, so that the constants returned describe the
working precision actually in use, whether double (fp64), extended
(fp80) or quadruple (fp128).

### Copyright and license (public domain)

The SLATEC Common Mathematical Library is in the public domain, so no
license conditions attach to this file. From the SLATEC guide published
at <https://www.netlib.org/slatec/guide>:

```text
The Library is in the public domain and distributed by the Energy
Science and Technology Software Center.
```

and, on the requirements for routines accepted into the library:

```text
The SLATEC Library is intended to have no restriction on its
distribution; therefore, new routines must be in the public domain.
```
