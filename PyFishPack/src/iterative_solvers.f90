module iterative_solvers

    use, intrinsic :: iso_c_binding, only: &
        c_double, &
        c_int, &
        c_ptr, &
        c_f_pointer

    implicit none

contains

    subroutine pyfp_sor_standard1d(n, s_ptr, a_ptr, b_ptr, f_ptr, delx, bcx, &
        optarg, undef, mxloop, tolerance, overflow, relerr, loops) &
        bind(C, name="pyfp_sor_standard1d")

        integer(c_int), value :: n, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, f_ptr
        real(c_double), value :: delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:), a(:), b(:), f(:)

        call c_f_pointer(s_ptr, s, [n])
        call c_f_pointer(a_ptr, a, [n])
        call c_f_pointer(b_ptr, b, [n])
        call c_f_pointer(f_ptr, f, [n])

        call sor_standard1d(n, s, a, b, f, delx, bcx, optarg, undef, &
            mxloop, tolerance, overflow, relerr, loops)

    end subroutine pyfp_sor_standard1d

    subroutine pyfp_sor_standard2d(ny, nx, s_ptr, a_ptr, b_ptr, ccoef_ptr, f_ptr, &
        dely, delx, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, &
        relerr, loops) bind(C, name="pyfp_sor_standard2d")

        integer(c_int), value :: ny, nx, bcy, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, ccoef_ptr, f_ptr
        real(c_double), value :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:, :), a(:, :), b(:, :), ccoef(:, :), f(:, :)

        call c_f_pointer(s_ptr, s, [ny, nx])
        call c_f_pointer(a_ptr, a, [ny, nx])
        call c_f_pointer(b_ptr, b, [ny, nx])
        call c_f_pointer(ccoef_ptr, ccoef, [ny, nx])
        call c_f_pointer(f_ptr, f, [ny, nx])

        call sor_standard2d(ny, nx, s, a, b, ccoef, f, dely, delx, bcy, bcx, &
            optarg, undef, mxloop, tolerance, overflow, relerr, loops)

    end subroutine pyfp_sor_standard2d

    subroutine pyfp_sor_general2d(ny, nx, s_ptr, a_ptr, b_ptr, ccoef_ptr, &
        d_ptr, e_ptr, fcoef_ptr, g_ptr, dely, delx, bcy, bcx, optarg, undef, &
        mxloop, tolerance, overflow, relerr, loops) bind(C, name="pyfp_sor_general2d")

        integer(c_int), value :: ny, nx, bcy, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, ccoef_ptr, d_ptr, e_ptr, fcoef_ptr, g_ptr
        real(c_double), value :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:, :), a(:, :), b(:, :), ccoef(:, :)
        real(c_double), pointer :: d(:, :), e(:, :), fcoef(:, :), g(:, :)

        call c_f_pointer(s_ptr, s, [ny, nx])
        call c_f_pointer(a_ptr, a, [ny, nx])
        call c_f_pointer(b_ptr, b, [ny, nx])
        call c_f_pointer(ccoef_ptr, ccoef, [ny, nx])
        call c_f_pointer(d_ptr, d, [ny, nx])
        call c_f_pointer(e_ptr, e, [ny, nx])
        call c_f_pointer(fcoef_ptr, fcoef, [ny, nx])
        call c_f_pointer(g_ptr, g, [ny, nx])

        call sor_general2d(ny, nx, s, a, b, ccoef, d, e, fcoef, g, dely, delx, &
            bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

    end subroutine pyfp_sor_general2d

    subroutine sor_standard1d(n, s, a, b, f, delx, bcx, optarg, undef, &
        mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: n, bcx, mxloop
        real(c_double), intent(inout) :: s(n)
        real(c_double), intent(in) :: a(n), b(n), f(n), delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, loop
        real(c_double) :: temp, norm, norm_prev, delx_sqr, denom

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        delx_sqr = delx * delx
        loop = 0

        do
            if (bcx == 1_c_int) then
                if (s(2) /= undef) s(1) = s(2)
                if (s(n - 1) /= undef) s(n) = s(n - 1)
            end if

            if (bcx == 2_c_int) then
                if (f(1) /= undef .and. a(1) /= undef .and. a(2) /= undef .and. b(1) /= undef) then
                    denom = (a(2) + a(1)) / delx_sqr - b(1)
                    temp = (a(2) * (s(2) - s(1)) - a(1) * (s(1) - s(n))) / delx_sqr &
                        + (b(1) * s(1) - f(1))
                    s(1) = s(1) + optarg * temp / denom
                end if
            end if

            do i = 2, n - 1
                if (f(i) /= undef .and. a(i) /= undef .and. a(i + 1) /= undef .and. b(i) /= undef) then
                    denom = (a(i + 1) + a(i)) / delx_sqr - b(i)
                    temp = (a(i + 1) * (s(i + 1) - s(i)) - a(i) * (s(i) - s(i - 1))) / delx_sqr &
                        + (b(i) * s(i) - f(i))
                    s(i) = s(i) + optarg * temp / denom
                end if
            end do

            if (bcx == 2_c_int) then
                if (f(n) /= undef .and. a(n) /= undef .and. a(1) /= undef .and. b(n) /= undef) then
                    denom = (a(1) + a(n)) / delx_sqr - b(n)
                    temp = (a(1) * (s(1) - s(n)) - a(n) * (s(n) - s(n - 1))) / delx_sqr &
                        + (b(n) * s(n) - f(n))
                    s(n) = s(n) + optarg * temp / denom
                end if
            end if

            norm = abs_norm1d(n, s, undef)
            if (norm /= norm .or. norm > 1.0e100_c_double) then
                overflow = 1_c_int
                exit
            end if

            relerr = abs(norm - norm_prev) / norm_prev
            loops = int(loop, c_int)
            if (relerr < tolerance .or. loop >= mxloop .or. norm == 0.0_c_double) exit

            norm_prev = norm
            loop = loop + 1
        end do

    end subroutine sor_standard1d

    subroutine sor_standard2d(ny, nx, s, a, b, ccoef, f, dely, delx, bcy, bcx, &
        optarg, undef, mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: ny, nx, bcy, bcx, mxloop
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx), f(ny, nx)
        real(c_double), intent(in) :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, j, loop
        real(c_double) :: norm, norm_prev, ratio_sqr, ratio_qtr, delx_sqr

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        delx_sqr = delx * delx
        ratio_sqr = (delx / dely) ** 2
        ratio_qtr = (delx / dely) / 4.0_c_double
        loop = 0

        do
            if (bcy == 1_c_int) call apply_extend_y(ny, nx, s, undef, bcx)
            if (bcx == 1_c_int) call apply_extend_x(ny, nx, s, undef)

            do j = 2, ny - 1
                if (bcx == 2_c_int) then
                    i = 1
                    call update_standard2d_point(ny, nx, s, a, b, ccoef, f, j, i, &
                        delx_sqr, ratio_sqr, ratio_qtr, optarg, undef)
                end if

                do i = 2, nx - 1
                    call update_standard2d_point(ny, nx, s, a, b, ccoef, f, j, i, &
                        delx_sqr, ratio_sqr, ratio_qtr, optarg, undef)
                end do

                if (bcx == 2_c_int) then
                    i = nx
                    call update_standard2d_point(ny, nx, s, a, b, ccoef, f, j, i, &
                        delx_sqr, ratio_sqr, ratio_qtr, optarg, undef)
                end if
            end do

            norm = abs_norm2d(ny, nx, s, undef)
            if (norm /= norm .or. norm > 1.0e100_c_double) then
                overflow = 1_c_int
                exit
            end if

            relerr = abs(norm - norm_prev) / norm_prev
            loops = int(loop, c_int)
            if (relerr < tolerance .or. loop >= mxloop) exit

            norm_prev = norm
            loop = loop + 1
        end do

    end subroutine sor_standard2d

    subroutine sor_general2d(ny, nx, s, a, b, ccoef, d, e, fcoef, g, dely, delx, &
        bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: ny, nx, bcy, bcx, mxloop
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx)
        real(c_double), intent(in) :: d(ny, nx), e(ny, nx), fcoef(ny, nx), g(ny, nx)
        real(c_double), intent(in) :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, j, loop
        real(c_double) :: norm, norm_prev, ratio, ratio_sqr, ratio_qtr, delx_sqr

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        delx_sqr = delx * delx
        ratio = delx / dely
        ratio_sqr = ratio * ratio
        ratio_qtr = ratio / 4.0_c_double
        loop = 0

        do
            if (bcy == 1_c_int) call apply_extend_y(ny, nx, s, undef, bcx)
            if (bcx == 1_c_int) call apply_extend_x(ny, nx, s, undef)

            do j = 2, ny - 1
                if (bcx == 2_c_int) then
                    i = 1
                    call update_general2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, g, j, i, &
                        delx, delx_sqr, ratio, ratio_sqr, ratio_qtr, optarg, undef)
                end if

                do i = 2, nx - 1
                    call update_general2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, g, j, i, &
                        delx, delx_sqr, ratio, ratio_sqr, ratio_qtr, optarg, undef)
                end do

                if (bcx == 2_c_int) then
                    i = nx
                    call update_general2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, g, j, i, &
                        delx, delx_sqr, ratio, ratio_sqr, ratio_qtr, optarg, undef)
                end if
            end do

            norm = abs_norm2d(ny, nx, s, undef)
            if (norm /= norm .or. norm > 1.0e100_c_double) then
                overflow = 1_c_int
                exit
            end if

            relerr = abs(norm - norm_prev) / norm_prev
            loops = int(loop, c_int)
            if (relerr < tolerance .or. loop >= mxloop) exit

            norm_prev = norm
            loop = loop + 1
        end do

    end subroutine sor_general2d

    subroutine update_standard2d_point(ny, nx, s, a, b, ccoef, f, j, i, &
        delx_sqr, ratio_sqr, ratio_qtr, optarg, undef)

        integer, intent(in) :: ny, nx, j, i
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx), f(ny, nx)
        real(c_double), intent(in) :: delx_sqr, ratio_sqr, ratio_qtr, optarg, undef

        integer :: im, ip
        real(c_double) :: temp, denom

        im = i - 1
        ip = i + 1
        if (i == 1) im = nx
        if (i == nx) ip = 1

        if (f(j, i) == undef) return
        if (a(j + 1, i) == undef .or. a(j, i) == undef) return
        if (b(j, ip) == undef .or. b(j, im) == undef .or. b(j + 1, i) == undef .or. b(j - 1, i) == undef) return
        if (ccoef(j, ip) == undef .or. ccoef(j, i) == undef) return

        temp = (a(j + 1, i) * (s(j + 1, i) - s(j, i)) - a(j, i) * (s(j, i) - s(j - 1, i))) &
            * ratio_sqr &
            + (b(j + 1, ip) * (s(j + 1, ip) - s(j + 1, im)) &
            - b(j - 1, i) * (s(j - 1, i) - s(j - 1, im))) * ratio_qtr &
            + (b(j, ip) * (s(j + 1, ip) - s(j - 1, ip)) &
            - b(j, im) * (s(j + 1, im) - s(j - 1, im))) * ratio_qtr &
            + ccoef(j, ip) * (s(j, ip) - s(j, i)) &
            - ccoef(j, i) * (s(j, i) - s(j, im)) &
            - f(j, i) * delx_sqr

        denom = (a(j + 1, i) + a(j, i)) * ratio_sqr + ccoef(j, ip) + ccoef(j, i)
        s(j, i) = s(j, i) + optarg * temp / denom

    end subroutine update_standard2d_point

    subroutine update_general2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, g, j, i, &
        delx, delx_sqr, ratio, ratio_sqr, ratio_qtr, optarg, undef)

        integer, intent(in) :: ny, nx, j, i
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx)
        real(c_double), intent(in) :: d(ny, nx), e(ny, nx), fcoef(ny, nx), g(ny, nx)
        real(c_double), intent(in) :: delx, delx_sqr, ratio, ratio_sqr, ratio_qtr, optarg, undef

        integer :: im, ip
        real(c_double) :: temp, denom

        im = i - 1
        ip = i + 1
        if (i == 1) im = nx
        if (i == nx) ip = 1

        if (g(j, i) == undef) return
        if (a(j, i) == undef .or. b(j, i) == undef .or. ccoef(j, i) == undef) return
        if (d(j, i) == undef .or. e(j, i) == undef .or. fcoef(j, i) == undef) return

        temp = a(j, i) * ((s(j + 1, i) - s(j, i)) - (s(j, i) - s(j - 1, i))) * ratio_sqr &
            + b(j, i) * ((s(j + 1, ip) - s(j - 1, ip)) - (s(j + 1, im) - s(j - 1, im))) * ratio_qtr &
            + ccoef(j, i) * ((s(j, ip) - s(j, i)) - (s(j, i) - s(j, im))) &
            + (d(j, i) * (s(j + 1, i) - s(j - 1, i)) * ratio &
            + e(j, i) * (s(j, ip) - s(j, im))) * delx / 2.0_c_double &
            + (fcoef(j, i) * s(j, i) - g(j, i)) * delx_sqr

        denom = (a(j, i) * ratio_sqr + ccoef(j, i)) * 2.0_c_double - fcoef(j, i) * delx_sqr
        s(j, i) = s(j, i) + optarg * temp / denom

    end subroutine update_general2d_point

    subroutine apply_extend_y(ny, nx, s, undef, bcx)
        integer, intent(in) :: ny, nx, bcx
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: undef
        integer :: i

        if (bcx == 2_c_int) then
            do i = 1, nx
                if (s(2, i) /= undef) s(1, i) = s(2, i)
                if (s(ny - 1, i) /= undef) s(ny, i) = s(ny - 1, i)
            end do
        else
            do i = 2, nx - 1
                if (s(2, i) /= undef) s(1, i) = s(2, i)
                if (s(ny - 1, i) /= undef) s(ny, i) = s(ny - 1, i)
            end do
            if (s(2, 2) /= undef) s(1, 1) = s(2, 2)
            if (s(2, nx - 1) /= undef) s(1, nx) = s(2, nx - 1)
            if (s(ny - 1, 2) /= undef) s(ny, 1) = s(ny - 1, 2)
            if (s(ny - 1, nx - 1) /= undef) s(ny, nx) = s(ny - 1, nx - 1)
        end if
    end subroutine apply_extend_y

    subroutine apply_extend_x(ny, nx, s, undef)
        integer, intent(in) :: ny, nx
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: undef
        integer :: j

        do j = 2, ny - 1
            if (s(j, 2) /= undef) s(j, 1) = s(j, 2)
            if (s(j, nx - 1) /= undef) s(j, nx) = s(j, nx - 1)
        end do
    end subroutine apply_extend_x

    real(c_double) function abs_norm1d(n, s, undef) result(norm)
        integer, intent(in) :: n
        real(c_double), intent(in) :: s(n), undef
        integer :: i, count

        norm = 0.0_c_double
        count = 0
        do i = 1, n
            if (s(i) /= undef) then
                norm = norm + abs(s(i))
                count = count + 1
            end if
        end do
        if (count > 0) then
            norm = norm / real(count, c_double)
        else
            norm = huge(1.0_c_double)
        end if
    end function abs_norm1d

    real(c_double) function abs_norm2d(ny, nx, s, undef) result(norm)
        integer, intent(in) :: ny, nx
        real(c_double), intent(in) :: s(ny, nx), undef
        integer :: i, j, count

        norm = 0.0_c_double
        count = 0
        do i = 1, nx
            do j = 1, ny
                if (s(j, i) /= undef) then
                    norm = norm + abs(s(j, i))
                    count = count + 1
                end if
            end do
        end do
        if (count > 0) then
            norm = norm / real(count, c_double)
        else
            norm = huge(1.0_c_double)
        end if
    end function abs_norm2d

end module iterative_solvers
