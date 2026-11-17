! Copyright (C) 2026 Martin Köhler <koehlerm(AT)mpi-magdeburg.mpg.de>
!
! This file is part of qrupdate-ng.
!
! qrupdate-ng is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version. See COPYING for the full license text.

module qrupdate_real
    use wp_def, only: wp
    use qrupdate_linalg, only: daxpy, dcopy, ddot, dnrm2, dscal, drot, &
        dlartg, dlamch, lsame
    use qrupdate_error, only: qrupdate_xerror
    implicit none
    private

    public :: real_qr1up
    public :: real_qrinc
    public :: real_qrdec
    public :: real_qrinr
    public :: real_qrder
    public :: real_qrshc
    public :: real_gqvec

contains

    subroutine real_qrtv1(n, u, w)
        integer, intent(in) :: n
        real(wp), intent(inout) :: u(*)
        real(wp), intent(out) :: w(*)
        real(wp) :: rr, t
        integer :: i

        if (n <= 0) return

        rr = u(n)
        do i = n - 1, 1, -1
            call dlartg(u(i), rr, w(i), u(i + 1), t)
            rr = t
        end do
        u(1) = rr
    end subroutine real_qrtv1

    subroutine real_qrot(dir, m, n, Q, ldq, c, s)
        character, intent(in) :: dir
        integer, intent(in) :: m, n, ldq
        real(wp), intent(inout) :: Q(ldq, *)
        real(wp), intent(in) :: c(*), s(*)
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
            call qrupdate_xerror('DQROT', info)
            return
        end if

        if (fwd) then
            do i = 1, n - 1
                call drot(m, Q(1, i), 1, Q(1, i + 1), 1, c(i), s(i))
            end do
        else
            do i = n - 1, 1, -1
                call drot(m, Q(1, i), 1, Q(1, i + 1), 1, c(i), s(i))
            end do
        end if
    end subroutine real_qrot

    subroutine real_qrqh(m, n, R, ldr, c, s)
        integer, intent(in) :: m, n, ldr
        real(wp), intent(inout) :: R(ldr, *)
        real(wp), intent(in) :: c(*), s(*)
        real(wp) :: t
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
            call qrupdate_xerror('DQRQH', info)
            return
        end if

        do i = 1, n
            ii = min(m - 1, i)
            t = R(ii + 1, i)
            do j = ii, 1, -1
                R(j + 1, i) = c(j) * t - s(j) * R(j, i)
                t = c(j) * R(j, i) + s(j) * t
            end do
            R(1, i) = t
        end do
    end subroutine real_qrqh

    subroutine real_qhqr(m, n, R, ldr, c, s)
        integer, intent(in) :: m, n, ldr
        real(wp), intent(inout) :: R(ldr, *)
        real(wp), intent(out) :: c(*), s(*)
        real(wp) :: t
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
            call qrupdate_xerror('DQHQR', info)
            return
        end if

        do i = 1, n
            t = R(1, i)
            ii = min(m, i)
            do j = 1, ii - 1
                R(j, i) = c(j) * t + s(j) * R(j + 1, i)
                t = c(j) * R(j + 1, i) - s(j) * t
            end do
            if (ii < m) then
                call dlartg(t, R(ii + 1, i), c(i), s(i), R(ii, i))
                R(ii + 1, i) = 0.0_wp
            else
                R(ii, i) = t
            end if
        end do
    end subroutine real_qhqr

    subroutine real_gqvec(m, n, Q, ldq, u)
        integer, intent(in) :: m, n, ldq
        real(wp), intent(in) :: Q(ldq, *)
        real(wp), intent(out) :: u(*)
        real(wp) :: r
        integer :: info, i, j

        if (m == 0) return
        if (n == 0) then
            u(1) = 1.0_wp
            do i = 2, m
                u(i) = 0.0_wp
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
            call qrupdate_xerror('DGQVEC', info)
            return
        end if

        j = 1
        r = 0.0_wp
        do while (r == 0.0_wp)
            do i = 1, m
                u(i) = 0.0_wp
            end do
            u(j) = 1.0_wp
            do i = 1, n
                r = ddot(m, Q(1, i), 1, u, 1)
                call daxpy(m, -r, Q(1, i), 1, u, 1)
            end do
            r = dnrm2(m, u, 1)
            if (r == 0.0_wp) then
                j = j + 1
                if (j > m) then
                    stop 'fatal: impossible condition in DGQVEC'
                end if
            end if
        end do
        call dscal(m, 1.0_wp / r, u, 1)
    end subroutine real_gqvec

    subroutine real_ch1up(n, R, ldr, u, w)
        integer, intent(in) :: n, ldr
        real(wp), intent(inout) :: R(ldr, *), u(*)
        real(wp), intent(out) :: w(*)
        real(wp) :: rr, ui, t
        integer :: i, j

        do i = 1, n
            ui = u(i)
            do j = 1, i - 1
                t = w(j) * R(j, i) + u(j) * ui
                ui = w(j) * ui - u(j) * R(j, i)
                R(j, i) = t
            end do
            call dlartg(R(i, i), ui, w(i), u(i), rr)
            R(i, i) = rr
        end do
    end subroutine real_ch1up

    subroutine real_qr1up(m, n, k, Q, ldq, R, ldr, u, v, w)
        integer, intent(in) :: m, n, k, ldq, ldr
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *), u(*), v(*)
        real(wp), intent(out) :: w(*)
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
            call qrupdate_xerror('DQR1UP', info)
            return
        end if

        full = k == m
        ru = 1.0_wp
        if (.not. full) ru = dnrm2(m, u, 1)
        do i = 1, k
            w(i) = ddot(m, Q(1, i), 1, u, 1)
            if (.not. full) call daxpy(m, -w(i), Q(1, i), 1, u, 1)
        end do
        call real_qrtv1(k, w, w(k + 1))
        call real_qrqh(k, n, R, ldr, w(k + 1), w(2))
        call real_qrot('B', m, k, Q, ldq, w(k + 1), w(2))
        call daxpy(n, w(1), v, 1, R(1, 1), ldr)
        call real_qhqr(k, n, R, ldr, w(k + 1), w)
        call real_qrot('F', m, min(k, n + 1), Q, ldq, w(k + 1), w)
        if (full) return

        ruu = dnrm2(m, u, 1)
        ru = ru * dlamch('e')
        if (ruu <= ru) return
        call dscal(n, ruu, v, 1)
        call dscal(m, 1.0_wp / ruu, u, 1)
        call real_ch1up(n, R, ldr, v, w(k + 1))
        do i = 1, n
            call drot(m, Q(1, i), 1, u, 1, w(k + i), v(i))
        end do
    end subroutine real_qr1up

    subroutine real_qrinc(m, n, k, Q, ldq, R, ldr, j, x, w)
        integer, intent(in) :: m, n, k, ldq, ldr, j
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *)
        real(wp), intent(in) :: x(*)
        real(wp), intent(out) :: w(*)
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
            call qrupdate_xerror('DQRINC', info)
            return
        end if

        full = k == m
        do i = n, j, -1
            call dcopy(k, R(1, i), 1, R(1, i + 1), 1)
        end do
        if (full) then
            k1 = k
            do i = 1, k
                R(i, j) = ddot(m, Q(1, i), 1, x, 1)
            end do
        else
            k1 = k + 1
            do i = 1, n + 1
                R(k1, i) = 0.0_wp
            end do
            call dcopy(m, x, 1, Q(1, k1), 1)
            do i = 1, k
                R(i, j) = ddot(m, Q(1, i), 1, Q(1, k1), 1)
                call daxpy(m, -R(i, j), Q(1, i), 1, Q(1, k1), 1)
            end do
            rx = dnrm2(m, Q(1, k1), 1)
            R(k1, j) = rx
            if (rx == 0.0_wp) then
                call real_gqvec(m, k, Q, ldq, Q(1, k1))
            else
                call dscal(m, 1.0_wp / rx, Q(1, k1), 1)
            end if
        end if
        if (j > k) return

        call real_qrtv1(k1 + 1 - j, R(j, j), w)
        if (j <= n) then
            call real_qrqh(k1 + 1 - j, n + 1 - j, R(j, j + 1), ldr, &
                w, R(j + 1, j))
        end if
        call real_qrot('B', m, k1 + 1 - j, Q(1, j), ldq, w, R(j + 1, j))
        do i = j + 1, k1
            R(i, j) = 0.0_wp
        end do
    end subroutine real_qrinc

    subroutine real_qrdec(m, n, k, Q, ldq, R, ldr, j, w)
        integer, intent(in) :: m, n, k, ldq, ldr, j
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *)
        real(wp), intent(out) :: w(*)
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
            call qrupdate_xerror('DQRDEC', info)
            return
        end if

        do i = j, n - 1
            call dcopy(k, R(1, i + 1), 1, R(1, i), 1)
        end do
        if (j < k) then
            call real_qhqr(k + 1 - j, n - j, R(j, j), ldr, w, R(1, n))
            call real_qrot('F', m, min(k, n) + 1 - j, Q(1, j), ldq, &
                w, R(1, n))
        end if
    end subroutine real_qrdec

    subroutine real_qrinr(m, n, Q, ldq, R, ldr, j, x, w)
        integer, intent(in) :: m, n, ldq, ldr, j
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *), x(*)
        real(wp), intent(out) :: w(*)
        integer :: info, i, k

        info = 0
        if (n < 0) then
            info = 2
        else if (j < 1 .or. j > m + 1) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('DQRINR', info)
            return
        end if

        do i = m, 1, -1
            if (j > 1) call dcopy(j - 1, Q(1, i), 1, Q(1, i + 1), 1)
            Q(j, i + 1) = 0.0_wp
            if (j <= m) then
                call dcopy(m + 1 - j, Q(j, i), 1, Q(j + 1, i + 1), 1)
            end if
        end do
        do i = 1, j - 1
            Q(i, 1) = 0.0_wp
        end do
        Q(j, 1) = 1.0_wp
        do i = j + 1, m + 1
            Q(i, 1) = 0.0_wp
        end do
        do k = 1, n
            if (k < m) R(m + 1, k) = 0.0_wp
            do i = min(m, k), 1, -1
                R(i + 1, k) = R(i, k)
            end do
            R(1, k) = x(k)
        end do
        call real_qhqr(m + 1, n, R, ldr, w, x)
        call real_qrot('F', m + 1, min(m, n) + 1, Q, ldq, w, x)
    end subroutine real_qrinr

    subroutine real_qrder(m, n, Q, ldq, R, ldr, j, w)
        integer, intent(in) :: m, n, ldq, ldr, j
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *)
        real(wp), intent(out) :: w(*)
        integer :: info, i, k

        if (m == 1) return

        info = 0
        if (m < 1) then
            info = 1
        else if (j < 1 .or. j > m) then
            info = 7
        end if
        if (info /= 0) then
            call qrupdate_xerror('DQRDER', info)
            return
        end if

        call dcopy(m, Q(j, 1), ldq, w, 1)
        call real_qrtv1(m, w, w(m + 1))
        call real_qrot('B', m, m, Q, ldq, w(m + 1), w(2))
        do k = 1, m - 1
            if (j > 1) call dcopy(j - 1, Q(1, k + 1), 1, Q(1, k), 1)
            if (j < m) call dcopy(m - j, Q(j + 1, k + 1), 1, Q(j, k), 1)
        end do
        call real_qrqh(m, n, R, ldr, w(m + 1), w(2))
        do k = 1, n
            do i = 1, m - 1
                R(i, k) = R(i + 1, k)
            end do
        end do
    end subroutine real_qrder

    subroutine real_qrshc(m, n, k, Q, ldq, R, ldr, i, j, w)
        integer, intent(in) :: m, n, k, ldq, ldr, i, j
        real(wp), intent(inout) :: Q(ldq, *), R(ldr, *)
        real(wp), intent(out) :: w(*)
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
            call qrupdate_xerror('DQRSHC', info)
            return
        end if

        if (i < j) then
            call dcopy(k, R(1, i), 1, w, 1)
            do l = i, j - 1
                call dcopy(k, R(1, l + 1), 1, R(1, l), 1)
            end do
            call dcopy(k, w, 1, R(1, j), 1)
            if (i < k) then
                kk = min(k, j)
                call real_qhqr(kk + 1 - i, n + 1 - i, R(i, i), ldr, &
                    w(k + 1), w)
                call real_qrot('F', m, kk + 1 - i, Q(1, i), ldq, &
                    w(k + 1), w)
            end if
        else if (j < i) then
            call dcopy(k, R(1, i), 1, w, 1)
            do l = i, j + 1, -1
                call dcopy(k, R(1, l - 1), 1, R(1, l), 1)
            end do
            call dcopy(k, w, 1, R(1, j), 1)
            if (j < k) then
                jj = min(j + 1, n)
                kk = min(k, i)
                call real_qrtv1(kk + 1 - j, R(j, j), w(k + 1))
                call real_qrqh(kk + 1 - j, n - j, R(j, jj), ldr, &
                    w(k + 1), R(j + 1, j))
                call real_qrot('B', m, kk + 1 - j, Q(1, j), ldq, &
                    w(k + 1), R(j + 1, j))
                do l = j + 1, kk
                    R(l, j) = 0.0_wp
                end do
            end if
        end if
    end subroutine real_qrshc

end module qrupdate_real
