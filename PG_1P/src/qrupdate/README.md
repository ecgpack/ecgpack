# Generic-precision QR updates

This directory vendors the generic QR-update implementation from Martin
Köhler's `qrupdate-ng` repository:

<https://gitlab.mpi-magdeburg.mpg.de/koehlerm/qrupdate-ng>

The imported revision is commit
`3fa4f77d6259c00cf0c36632287063f1277a3fe0`.

The public `qrupdate` module provides generic real/complex interfaces for
`qr1up`, `qrinc`, `qrdec`, `qrinr`, `qrder`, `qrshc`, and `gqvec`. Its
`real(wp)` and `complex(wp)` procedures use the same build-selected `wp_def`
module as the enclosing ECGPACK harness.

## Local modification

This copy differs from the imported revision by one addition made for ECGPACK
in November 2026: the generic interface `qr1up_unit` in `qrupdate.f90` and its
specific procedures `real_qr1up_unit` (`qrupdate_real.f90`) and
`complex_qr1up_unit` (`qrupdate_complex.f90`). They apply the rank-one update
`A + e_j*v^H` to a full QR factorization, where `e_j` is a unit vector. This
is `qr1up` with `u = e_j`, except that `Q^H*u` is read from row `j` of Q
instead of being accumulated by dot products. The upstream routines are
unchanged, and each affected file states the modification in its header.

The precision-independent bundled Netlib aggregates live one directory above
as `../BLAS.f` and `../LAPACK.f`. They provide the routines required by both
the QR updates and ECGPACK's existing linear-algebra paths. Keeping the single
provider in the conventional `src/` location avoids duplicate BLAS/LAPACK
sources. The aggregates intentionally use D/Z symbol names for build-selected
`real(wp)` and `complex(wp)` routines. Do not link a conventional fixed-double
BLAS/LAPACK into `wp=10` or `wp=16` builds.

The QR-update implementation is GPL-3.0-or-later; see `COPYING`. The imported
Netlib routines retain their upstream notices and licensing documentation in
the aggregate source files.
