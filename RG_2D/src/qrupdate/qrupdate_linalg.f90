! Copyright (C) 2026 Martin Köhler <koehlerm(AT)mpi-magdeburg.mpg.de>
!
! This file is part of qrupdate-ng.
!
! qrupdate-ng is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your
! option) any later version.
!
! Explicit interfaces to the precision-independent D/Z routines bundled in
! linalg/.  The routine names follow reference BLAS/LAPACK convention while
! their argument kinds are selected by wp_def for the whole build.
module qrupdate_linalg
    use wp_def, only: wp
    implicit none
    private

    public :: daxpy, dcopy, ddot, dnrm2, dscal, drot
    public :: zaxpy, zcopy, zdotc, dznrm2, zdscal, zrot
    public :: dlamch, dlartg, zlartg
    public :: lsame, xerbla

    interface
        subroutine daxpy(n, da, dx, incx, dy, incy)
            import wp
            integer, intent(in) :: n, incx, incy
            real(wp), intent(in) :: da, dx(*)
            real(wp), intent(inout) :: dy(*)
        end subroutine daxpy

        subroutine dcopy(n, dx, incx, dy, incy)
            import wp
            integer, intent(in) :: n, incx, incy
            real(wp), intent(in) :: dx(*)
            real(wp), intent(inout) :: dy(*)
        end subroutine dcopy

        function ddot(n, dx, incx, dy, incy) result(res)
            import wp
            integer, intent(in) :: n, incx, incy
            real(wp), intent(in) :: dx(*), dy(*)
            real(wp) :: res
        end function ddot

        function dnrm2(n, x, incx) result(res)
            import wp
            integer, intent(in) :: n, incx
            real(wp), intent(in) :: x(*)
            real(wp) :: res
        end function dnrm2

        subroutine dscal(n, da, dx, incx)
            import wp
            integer, intent(in) :: n, incx
            real(wp), intent(in) :: da
            real(wp), intent(inout) :: dx(*)
        end subroutine dscal

        subroutine drot(n, dx, incx, dy, incy, c, s)
            import wp
            integer, intent(in) :: n, incx, incy
            real(wp), intent(inout) :: dx(*), dy(*)
            real(wp), intent(in) :: c, s
        end subroutine drot

        subroutine zaxpy(n, za, zx, incx, zy, incy)
            import wp
            integer, intent(in) :: n, incx, incy
            complex(wp), intent(in) :: za, zx(*)
            complex(wp), intent(inout) :: zy(*)
        end subroutine zaxpy

        subroutine zcopy(n, zx, incx, zy, incy)
            import wp
            integer, intent(in) :: n, incx, incy
            complex(wp), intent(in) :: zx(*)
            complex(wp), intent(inout) :: zy(*)
        end subroutine zcopy

        function zdotc(n, zx, incx, zy, incy) result(res)
            import wp
            integer, intent(in) :: n, incx, incy
            complex(wp), intent(in) :: zx(*), zy(*)
            complex(wp) :: res
        end function zdotc

        function dznrm2(n, x, incx) result(res)
            import wp
            integer, intent(in) :: n, incx
            complex(wp), intent(in) :: x(*)
            real(wp) :: res
        end function dznrm2

        subroutine zdscal(n, da, zx, incx)
            import wp
            integer, intent(in) :: n, incx
            real(wp), intent(in) :: da
            complex(wp), intent(inout) :: zx(*)
        end subroutine zdscal

        subroutine zrot(n, cx, incx, cy, incy, c, s)
            import wp
            integer, intent(in) :: n, incx, incy
            complex(wp), intent(inout) :: cx(*), cy(*)
            real(wp), intent(in) :: c
            complex(wp), intent(in) :: s
        end subroutine zrot

        function dlamch(cmach) result(res)
            import wp
            character, intent(in) :: cmach
            real(wp) :: res
        end function dlamch

        subroutine dlartg(f, g, c, s, r)
            import wp
            real(wp), intent(in) :: f, g
            real(wp), intent(out) :: c, s, r
        end subroutine dlartg

        subroutine zlartg(f, g, c, s, r)
            import wp
            complex(wp), intent(in) :: f, g
            real(wp), intent(out) :: c
            complex(wp), intent(out) :: s, r
        end subroutine zlartg

        subroutine xerbla(srname, info)
            character(len=*), intent(in) :: srname
            integer, intent(in) :: info
        end subroutine xerbla
    end interface

contains

    pure logical function lsame(ca, cb)
        character, intent(in) :: ca, cb
        integer :: inta, intb, zcode

        lsame = ca == cb
        if (lsame) return

        zcode = ichar('Z')
        inta = ichar(ca)
        intb = ichar(cb)

        if (zcode == 90 .or. zcode == 122) then
            if (inta >= 97 .and. inta <= 122) inta = inta - 32
            if (intb >= 97 .and. intb <= 122) intb = intb - 32
        else if (zcode == 233 .or. zcode == 169) then
            if ((inta >= 129 .and. inta <= 137) .or. &
                (inta >= 145 .and. inta <= 153) .or. &
                (inta >= 162 .and. inta <= 169)) inta = inta + 64
            if ((intb >= 129 .and. intb <= 137) .or. &
                (intb >= 145 .and. intb <= 153) .or. &
                (intb >= 162 .and. intb <= 169)) intb = intb + 64
        else if (zcode == 218 .or. zcode == 250) then
            if (inta >= 225 .and. inta <= 250) inta = inta - 32
            if (intb >= 225 .and. intb <= 250) intb = intb - 32
        end if

        lsame = inta == intb
    end function lsame

end module qrupdate_linalg
