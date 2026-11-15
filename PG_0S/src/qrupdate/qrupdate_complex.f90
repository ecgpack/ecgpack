! Copyright (C) 2008, 2009 VZLU Prague, a.s., Czech Republic
! Copyright (C) 2026 Martin Koehler <koehlerm(AT)mpi-magdeburg.mpg.de>
!
! This file is part of qrupdate-ng.
!
! qrupdate-ng is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your
! option) any later version.
!
! qrupdate-ng is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.

module qrupdate_complex
    use wp_def, only: wp
    use qrupdate_linalg, only: dlamch, dznrm2, lsame, zaxpy, zcopy, zdotc, &
                              zdscal, zlartg, zrot
    use qrupdate_error, only: qrupdate_xerror
    implicit none
    private

    public :: complex_qr1up
    public :: complex_qrinc
    public :: complex_qrdec
    public :: complex_qrinr
    public :: complex_qrder
    public :: complex_qrshc
    public :: complex_gqvec

contains

    subroutine complex_qrtv1(n, u, w)
        integer, intent(in) :: n
        complex(wp), intent(inout) :: u(*)
        real(wp), intent(out) :: w(*)
        complex(wp) :: rr, t
        integer :: i

        if (n <= 0) return

        rr = u(n)
        do i = n - 1, 1, -1
            call zlartg(u(i), rr, w(i), u(i + 1), t)
            rr = t
        end do
        u(1) = rr
    end subroutine complex_qrtv1

    subroutine complex_qrot(dir, m, n, q, ldq, c, s)
        character, intent(in) :: dir
        integer, intent(in) :: m, n, ldq
        complex(wp), intent(inout) :: q(ldq, *)
        real(wp), intent(in) :: c(*)
        complex(wp), intent(in) :: s(*)
        logical :: fwd
        integer :: info, i

        if (m == 0 .or. n == 0 .or. n == 1) return

        info = 0
        fwd = lsame(dir, 'F')
        if (.not. (fwd .or. lsame(dir, 'B'))) then
            info = 1
        else if (m < 0) then
            info = 2
        else if (n < 0) then
            info = 3
        else if (ldq < m) then
            info = 5
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQROT', info)
            return
        end if

        if (fwd) then
            do i = 1, n - 1
                call zrot(m, q(1, i), 1, q(1, i + 1), 1, c(i), conjg(s(i)))
            end do
        else
            do i = n - 1, 1, -1
                call zrot(m, q(1, i), 1, q(1, i + 1), 1, c(i), conjg(s(i)))
            end do
        end if
    end subroutine complex_qrot

    subroutine complex_qrqh(m, n, r, ldr, c, s)
        integer, intent(in) :: m, n, ldr
        complex(wp), intent(inout) :: r(ldr, *)
        real(wp), intent(in) :: c(*)
        complex(wp), intent(in) :: s(*)
        complex(wp) :: t
        integer :: info, i, ii, j

        if (m == 0 .or. m == 1 .or. n == 0) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (ldr < m) then
            info = 4
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRQH', info)
            return
        end if

        do i = 1, n
            ii = min(m - 1, i)
            t = r(ii + 1, i)
            do j = ii, 1, -1
                r(j + 1, i) = c(j)*t - conjg(s(j))*r(j, i)
                t = c(j)*r(j, i) + s(j)*t
            end do
            r(1, i) = t
        end do
    end subroutine complex_qrqh

    subroutine complex_qhqr(m, n, r, ldr, c, s)
        integer, intent(in) :: m, n, ldr
        complex(wp), intent(inout) :: r(ldr, *)
        real(wp), intent(out) :: c(*)
        complex(wp), intent(out) :: s(*)
        complex(wp) :: t
        integer :: info, i, ii, j

        if (m == 0 .or. m == 1 .or. n == 0) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (ldr < m) then
            info = 4
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQHQR', info)
            return
        end if

        do i = 1, n
            t = r(1, i)
            ii = min(m, i)
            do j = 1, ii - 1
                r(j, i) = c(j)*t + s(j)*r(j + 1, i)
                t = c(j)*r(j + 1, i) - conjg(s(j))*t
            end do
            if (ii < m) then
                call zlartg(t, r(ii + 1, i), c(i), s(i), r(ii, i))
                r(ii + 1, i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            else
                r(ii, i) = t
            end if
        end do
    end subroutine complex_qhqr

    subroutine complex_gqvec(m, n, q, ldq, u)
        integer, intent(in) :: m, n, ldq
        complex(wp), intent(in) :: q(ldq, *)
        complex(wp), intent(out) :: u(*)
        real(wp) :: rnorm
        complex(wp) :: rc
        integer :: info, i, j

        if (m == 0) return
        if (n == 0) then
            u(1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
            do i = 2, m
                u(i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            end do
            return
        end if

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (ldq < m) then
            info = 4
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZGQVEC', info)
            return
        end if

        j = 1
        rnorm = 0.0_wp
        do while (rnorm == 0.0_wp)
            do i = 1, m
                u(i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            end do
            u(j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
            do i = 1, n
                rc = zdotc(m, q(1, i), 1, u, 1)
                call zaxpy(m, -rc, q(1, i), 1, u, 1)
            end do
            rnorm = dznrm2(m, u, 1)
            if (rnorm == 0.0_wp) then
                j = j + 1
                if (j > m) then
                    stop 'fatal: impossible condition in ZGQVEC'
                end if
            end if
        end do
        call zdscal(m, 1.0_wp/rnorm, u, 1)
    end subroutine complex_gqvec

    subroutine complex_ch1up(n, r, ldr, u, w)
        integer, intent(in) :: n, ldr
        complex(wp), intent(inout) :: r(ldr, *), u(*)
        real(wp), intent(out) :: w(*)
        complex(wp) :: rr, ui, t
        integer :: i, j

        do i = 1, n
            ui = conjg(u(i))
            do j = 1, i - 1
                t = w(j)*r(j, i) + u(j)*ui
                ui = w(j)*ui - conjg(u(j))*r(j, i)
                r(j, i) = t
            end do
            call zlartg(r(i, i), ui, w(i), u(i), rr)
            r(i, i) = rr
        end do
    end subroutine complex_ch1up

    subroutine complex_axcpy(n, a, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(wp), intent(in) :: a
        complex(wp), intent(in) :: x(*)
        complex(wp), intent(inout) :: y(*)
        integer :: i, ix, iy

        if (n <= 0) return
        if (incx /= 1 .or. incy /= 1) then
            ix = 1
            iy = 1
            if (incx < 0) ix = (-n + 1)*incx + 1
            if (incy < 0) iy = (-n + 1)*incy + 1
            do i = 1, n
                y(iy) = y(iy) + a*conjg(x(ix))
                ix = ix + incx
                iy = iy + incy
            end do
        else
            do i = 1, n
                y(i) = y(i) + a*conjg(x(i))
            end do
        end if
    end subroutine complex_axcpy

    subroutine complex_qr1up(m, n, k, q, ldq, r, ldr, u, v, w, rw)
        integer, intent(in) :: m, n, k, ldq, ldr
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *), u(*), v(*)
        complex(wp), intent(out) :: w(*)
        real(wp), intent(out) :: rw(*)
        real(wp) :: ru, ruu
        integer :: info, i
        logical :: full

        if (k == 0 .or. n == 0) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (k /= m .and. (k /= n .or. n > m)) then
            info = 3
        else if (ldq < m) then
            info = 5
        else if (ldr < k) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQR1UP', info)
            return
        end if

        full = k == m
        ru = 1.0_wp
        if (.not. full) ru = dznrm2(m, u, 1)
        do i = 1, k
            w(i) = zdotc(m, q(1, i), 1, u, 1)
            if (.not. full) call zaxpy(m, -w(i), q(1, i), 1, u, 1)
        end do
        call complex_qrtv1(k, w, rw)
        call complex_qrqh(k, n, r, ldr, rw, w(2))
        call complex_qrot('B', m, k, q, ldq, rw, w(2))
        call complex_axcpy(n, w(1), v, 1, r(1, 1), ldr)
        call complex_qhqr(k, n, r, ldr, rw, w)
        call complex_qrot('F', m, min(k, n + 1), q, ldq, rw, w)
        if (full) return

        ruu = dznrm2(m, u, 1)
        ru = ru*dlamch('E')
        if (ruu <= ru) return
        call zdscal(n, ruu, v, 1)
        call zdscal(m, 1.0_wp/ruu, u, 1)
        call complex_ch1up(n, r, ldr, v, rw)
        do i = 1, n
            call zrot(m, q(1, i), 1, u, 1, rw(i), conjg(v(i)))
        end do
    end subroutine complex_qr1up

    subroutine complex_qrinc(m, n, k, q, ldq, r, ldr, j, x, rw)
        integer, intent(in) :: m, n, k, ldq, ldr, j
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *)
        complex(wp), intent(in) :: x(*)
        real(wp), intent(out) :: rw(*)
        real(wp) :: rx
        integer :: info, i, k1
        logical :: full

        if (m == 0) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (k /= m .and. (k /= n .or. n >= m)) then
            info = 3
        else if (ldq < m) then
            info = 5
        else if (ldr < min(m, k + 1)) then
            info = 7
        else if (j < 1 .or. j > n + 1) then
            info = 8
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRINC', info)
            return
        end if

        full = k == m
        do i = n, j, -1
            call zcopy(k, r(1, i), 1, r(1, i + 1), 1)
        end do
        if (full) then
            k1 = k
            do i = 1, k
                r(i, j) = zdotc(m, q(1, i), 1, x, 1)
            end do
        else
            k1 = k + 1
            do i = 1, n + 1
                r(k1, i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            end do
            call zcopy(m, x, 1, q(1, k1), 1)
            do i = 1, k
                r(i, j) = zdotc(m, q(1, i), 1, q(1, k1), 1)
                call zaxpy(m, -r(i, j), q(1, i), 1, q(1, k1), 1)
            end do
            rx = dznrm2(m, q(1, k1), 1)
            r(k1, j) = cmplx(rx, 0.0_wp, kind=wp)
            if (rx == 0.0_wp) then
                call complex_gqvec(m, k, q, ldq, q(1, k1))
            else
                call zdscal(m, 1.0_wp/rx, q(1, k1), 1)
            end if
        end if
        if (j > k) return

        call complex_qrtv1(k1 + 1 - j, r(j, j), rw)
        if (j <= n) then
            call complex_qrqh(k1 + 1 - j, n + 1 - j, r(j, j + 1), ldr, &
                              rw, r(j + 1, j))
        end if
        call complex_qrot('B', m, k1 + 1 - j, q(1, j), ldq, rw, r(j + 1, j))
        do i = j + 1, k1
            r(i, j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        end do
    end subroutine complex_qrinc

    subroutine complex_qrdec(m, n, k, q, ldq, r, ldr, j, rw)
        integer, intent(in) :: m, n, k, ldq, ldr, j
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *)
        real(wp), intent(out) :: rw(*)
        integer :: info, i

        if (m == 0 .or. n == 0 .or. j == n) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (k /= m .and. (k /= n .or. n >= m)) then
            info = 3
        else if (ldq < m) then
            info = 5
        else if (ldr < k) then
            info = 7
        else if (j < 1 .or. j > n) then
            info = 8
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRDEC', info)
            return
        end if

        do i = j, n - 1
            call zcopy(k, r(1, i + 1), 1, r(1, i), 1)
        end do
        if (j < k) then
            call complex_qhqr(k + 1 - j, n - j, r(j, j), ldr, rw, r(1, n))
            call complex_qrot('F', m, min(k, n) + 1 - j, q(1, j), ldq, &
                              rw, r(1, n))
        end if
    end subroutine complex_qrdec

    subroutine complex_qrinr(m, n, q, ldq, r, ldr, j, x, rw)
        integer, intent(in) :: m, n, j, ldq, ldr
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *), x(*)
        real(wp), intent(out) :: rw(*)
        integer :: info, i, k

        info = 0
        if (n < 0) then
            info = 2
        else if (j < 1 .or. j > m + 1) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRINR', info)
            return
        end if

        do i = m, 1, -1
            if (j > 1) call zcopy(j - 1, q(1, i), 1, q(1, i + 1), 1)
            q(j, i + 1) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            if (j <= m) then
                call zcopy(m + 1 - j, q(j, i), 1, q(j + 1, i + 1), 1)
            end if
        end do
        do i = 1, j - 1
            q(i, 1) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        end do
        q(j, 1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
        do i = j + 1, m + 1
            q(i, 1) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        end do
        do k = 1, n
            if (k < m) r(m + 1, k) = cmplx(0.0_wp, 0.0_wp, kind=wp)
            do i = min(m, k), 1, -1
                r(i + 1, k) = r(i, k)
            end do
            r(1, k) = x(k)
        end do
        call complex_qhqr(m + 1, n, r, ldr, rw, x)
        call complex_qrot('F', m + 1, min(m, n) + 1, q, ldq, rw, x)
    end subroutine complex_qrinr

    subroutine complex_qrder(m, n, q, ldq, r, ldr, j, w, rw)
        integer, intent(in) :: m, n, ldq, ldr, j
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *)
        complex(wp), intent(out) :: w(*)
        real(wp), intent(out) :: rw(*)
        integer :: info, i, k

        if (m == 1) return

        info = 0
        if (m < 1) then
            info = 1
        else if (j < 1 .or. j > m) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRDER', info)
            return
        end if

        do k = 1, m
            w(k) = conjg(q(j, k))
        end do
        call complex_qrtv1(m, w, rw)
        call complex_qrot('B', m, m, q, ldq, rw, w(2))
        do k = 1, m - 1
            if (j > 1) call zcopy(j - 1, q(1, k + 1), 1, q(1, k), 1)
            if (j < m) call zcopy(m - j, q(j + 1, k + 1), 1, q(j, k), 1)
        end do
        call complex_qrqh(m, n, r, ldr, rw, w(2))
        do k = 1, n
            do i = 1, m - 1
                r(i, k) = r(i + 1, k)
            end do
        end do
    end subroutine complex_qrder

    subroutine complex_qrshc(m, n, k, q, ldq, r, ldr, i, j, w, rw)
        integer, intent(in) :: m, n, k, ldq, ldr, i, j
        complex(wp), intent(inout) :: q(ldq, *), r(ldr, *)
        complex(wp), intent(out) :: w(*)
        real(wp), intent(out) :: rw(*)
        integer :: info, jj, kk, l

        if (m == 0 .or. n == 1) return

        info = 0
        if (m < 0) then
            info = 1
        else if (n < 0) then
            info = 2
        else if (k /= m .and. (k /= n .or. n > m)) then
            info = 3
        else if (i < 1 .or. i > n) then
            info = 6
        else if (j < 1 .or. j > n) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('ZQRSHC', info)
            return
        end if

        if (i < j) then
            call zcopy(k, r(1, i), 1, w, 1)
            do l = i, j - 1
                call zcopy(k, r(1, l + 1), 1, r(1, l), 1)
            end do
            call zcopy(k, w, 1, r(1, j), 1)
            if (i < k) then
                kk = min(k, j)
                call complex_qhqr(kk + 1 - i, n + 1 - i, r(i, i), ldr, rw, w)
                call complex_qrot('F', m, kk + 1 - i, q(1, i), ldq, rw, w)
            end if
        else if (j < i) then
            call zcopy(k, r(1, i), 1, w, 1)
            do l = i, j + 1, -1
                call zcopy(k, r(1, l - 1), 1, r(1, l), 1)
            end do
            call zcopy(k, w, 1, r(1, j), 1)
            if (j < k) then
                jj = min(j + 1, n)
                kk = min(k, i)
                call complex_qrtv1(kk + 1 - j, r(j, j), rw)
                call complex_qrqh(kk + 1 - j, n - j, r(j, jj), ldr, rw, &
                                  r(j + 1, j))
                call complex_qrot('B', m, kk + 1 - j, q(1, j), ldq, rw, &
                                  r(j + 1, j))
                do l = j + 1, kk
                    r(l, j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
                end do
            end if
        end if
    end subroutine complex_qrshc

end module qrupdate_complex
