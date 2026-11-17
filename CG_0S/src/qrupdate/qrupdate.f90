! Copyright (C) 2026 Martin Köhler <koehlerm(AT)mpi-magdeburg.mpg.de>
!
! This file is part of qrupdate-ng.
!
! qrupdate-ng is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version. See COPYING for the full license text.

module qrupdate
    use wp_def, only: wp
    use qrupdate_error, only: error_handler_if, qrupdate_set_error, &
        qrupdate_set_error_data
    use qrupdate_real, only: real_qr1up, real_qrinc, real_qrdec, real_qrinr, &
        real_qrder, real_qrshc, real_gqvec
    use qrupdate_complex, only: complex_qr1up, complex_qrinc, complex_qrdec, &
        complex_qrinr, complex_qrder, complex_qrshc, complex_gqvec
    implicit none
    private

    public :: wp
    public :: qr1up, qrinc, qrdec, qrinr, qrder, qrshc, gqvec
    public :: error_handler_if, qrupdate_set_error, qrupdate_set_error_data

    interface qr1up
        procedure :: real_qr1up
        procedure :: complex_qr1up
    end interface qr1up

    interface qrinc
        procedure :: real_qrinc
        procedure :: complex_qrinc
    end interface qrinc

    interface qrdec
        procedure :: real_qrdec
        procedure :: complex_qrdec
    end interface qrdec

    interface qrinr
        procedure :: real_qrinr
        procedure :: complex_qrinr
    end interface qrinr

    interface qrder
        procedure :: real_qrder
        procedure :: complex_qrder
    end interface qrder

    interface qrshc
        procedure :: real_qrshc
        procedure :: complex_qrshc
    end interface qrshc

    interface gqvec
        procedure :: real_gqvec
        procedure :: complex_gqvec
    end interface gqvec

end module qrupdate
