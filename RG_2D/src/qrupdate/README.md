# Generic-precision QR updates

This directory vendors the generic QR-update implementation from
`do1exist7/qrupdate-ng` at commit
`3fa4f77d6259c00cf0c36632287063f1277a3fe0`.

The public `qrupdate` module provides generic real/complex interfaces for
`qr1up`, `qrinc`, `qrdec`, `qrinr`, `qrder`, `qrshc`, and `gqvec`. Its
`real(wp)` and `complex(wp)` procedures use the same build-selected `wp_def`
module as the enclosing ECGPACK harness.

`BLAS.f` and `LAPACK.f` are the precision-independent bundled Netlib
aggregates required by both QR updates and the existing optional BLAS paths.
They intentionally use D/Z symbol names for build-selected `real(wp)` and
`complex(wp)` routines. Do not link a conventional fixed-double BLAS/LAPACK
into `wp=10` or `wp=16` builds.

The QR-update implementation is GPL-3.0-or-later; see `COPYING`. The imported
Netlib routines retain their upstream notices and licensing documentation in
the aggregate source files.
