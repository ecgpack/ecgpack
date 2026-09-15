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
