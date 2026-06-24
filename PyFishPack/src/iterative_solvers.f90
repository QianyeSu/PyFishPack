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

    subroutine pyfp_sor_standard3d(nz, ny, nx, s_ptr, a_ptr, b_ptr, ccoef_ptr, f_ptr, &
        delz, dely, delx, bcz, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, &
        relerr, loops) bind(C, name="pyfp_sor_standard3d")

        integer(c_int), value :: nz, ny, nx, bcz, bcy, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, ccoef_ptr, f_ptr
        real(c_double), value :: delz, dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:, :, :), a(:, :, :), b(:, :, :), ccoef(:, :, :), f(:, :, :)

        call c_f_pointer(s_ptr, s, [nz, ny, nx])
        call c_f_pointer(a_ptr, a, [nz, ny, nx])
        call c_f_pointer(b_ptr, b, [nz, ny, nx])
        call c_f_pointer(ccoef_ptr, ccoef, [nz, ny, nx])
        call c_f_pointer(f_ptr, f, [nz, ny, nx])

        call sor_standard3d(nz, ny, nx, s, a, b, ccoef, f, delz, dely, delx, &
            bcz, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

    end subroutine pyfp_sor_standard3d

    subroutine pyfp_sor_general3d(nz, ny, nx, s_ptr, a_ptr, b_ptr, ccoef_ptr, &
        d_ptr, e_ptr, fcoef_ptr, g_ptr, h_ptr, delz, dely, delx, bcz, bcy, bcx, &
        optarg, undef, mxloop, tolerance, overflow, relerr, loops) bind(C, name="pyfp_sor_general3d")

        integer(c_int), value :: nz, ny, nx, bcz, bcy, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, ccoef_ptr, d_ptr, e_ptr, fcoef_ptr, g_ptr, h_ptr
        real(c_double), value :: delz, dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:, :, :), a(:, :, :), b(:, :, :), ccoef(:, :, :)
        real(c_double), pointer :: d(:, :, :), e(:, :, :), fcoef(:, :, :), g(:, :, :), h(:, :, :)

        call c_f_pointer(s_ptr, s, [nz, ny, nx])
        call c_f_pointer(a_ptr, a, [nz, ny, nx])
        call c_f_pointer(b_ptr, b, [nz, ny, nx])
        call c_f_pointer(ccoef_ptr, ccoef, [nz, ny, nx])
        call c_f_pointer(d_ptr, d, [nz, ny, nx])
        call c_f_pointer(e_ptr, e, [nz, ny, nx])
        call c_f_pointer(fcoef_ptr, fcoef, [nz, ny, nx])
        call c_f_pointer(g_ptr, g, [nz, ny, nx])
        call c_f_pointer(h_ptr, h, [nz, ny, nx])

        call sor_general3d(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, &
            delz, dely, delx, bcz, bcy, bcx, optarg, undef, mxloop, tolerance, &
            overflow, relerr, loops)

    end subroutine pyfp_sor_general3d

    subroutine pyfp_sor_biharmonic2d(ny, nx, s_ptr, a_ptr, b_ptr, ccoef_ptr, &
        d_ptr, e_ptr, fcoef_ptr, g_ptr, h_ptr, icoef_ptr, jcoef_ptr, dely, delx, &
        bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops) &
        bind(C, name="pyfp_sor_biharmonic2d")

        integer(c_int), value :: ny, nx, bcy, bcx, mxloop
        type(c_ptr), value :: s_ptr, a_ptr, b_ptr, ccoef_ptr, d_ptr, e_ptr
        type(c_ptr), value :: fcoef_ptr, g_ptr, h_ptr, icoef_ptr, jcoef_ptr
        real(c_double), value :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        real(c_double), pointer :: s(:, :), a(:, :), b(:, :), ccoef(:, :)
        real(c_double), pointer :: d(:, :), e(:, :), fcoef(:, :), g(:, :)
        real(c_double), pointer :: h(:, :), icoef(:, :), jcoef(:, :)

        call c_f_pointer(s_ptr, s, [ny, nx])
        call c_f_pointer(a_ptr, a, [ny, nx])
        call c_f_pointer(b_ptr, b, [ny, nx])
        call c_f_pointer(ccoef_ptr, ccoef, [ny, nx])
        call c_f_pointer(d_ptr, d, [ny, nx])
        call c_f_pointer(e_ptr, e, [ny, nx])
        call c_f_pointer(fcoef_ptr, fcoef, [ny, nx])
        call c_f_pointer(g_ptr, g, [ny, nx])
        call c_f_pointer(h_ptr, h, [ny, nx])
        call c_f_pointer(icoef_ptr, icoef, [ny, nx])
        call c_f_pointer(jcoef_ptr, jcoef, [ny, nx])

        call sor_biharmonic2d(ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, icoef, jcoef, &
            dely, delx, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

    end subroutine pyfp_sor_biharmonic2d

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

    subroutine sor_standard3d(nz, ny, nx, s, a, b, ccoef, f, delz, dely, delx, &
        bcz, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: nz, ny, nx, bcz, bcy, bcx, mxloop
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: a(nz, ny, nx), b(nz, ny, nx), ccoef(nz, ny, nx), f(nz, ny, nx)
        real(c_double), intent(in) :: delz, dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, j, k, loop
        real(c_double) :: norm, norm_prev, ratio_z_sqr, ratio_y_sqr, delx_sqr

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        delx_sqr = delx * delx
        ratio_z_sqr = (delx / delz) ** 2
        ratio_y_sqr = (delx / dely) ** 2
        loop = 0

        do
            if (bcz == 1_c_int) call apply_extend_z3d(nz, ny, nx, s, undef)
            if (bcy == 1_c_int) call apply_extend_y3d(nz, ny, nx, s, undef, bcx)
            if (bcx == 1_c_int) call apply_extend_x3d(nz, ny, nx, s, undef)

            do k = 2, nz - 1
                do j = 2, ny - 1
                    if (bcx == 2_c_int) then
                        i = 1
                        call update_standard3d_point(nz, ny, nx, s, a, b, ccoef, f, k, j, i, &
                            delx_sqr, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end if

                    do i = 2, nx - 1
                        call update_standard3d_point(nz, ny, nx, s, a, b, ccoef, f, k, j, i, &
                            delx_sqr, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end do

                    if (bcx == 2_c_int) then
                        i = nx
                        call update_standard3d_point(nz, ny, nx, s, a, b, ccoef, f, k, j, i, &
                            delx_sqr, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end if
                end do
            end do

            norm = abs_norm3d(nz, ny, nx, s, undef)
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

    end subroutine sor_standard3d

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

    subroutine sor_general3d(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, &
        delz, dely, delx, bcz, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: nz, ny, nx, bcz, bcy, bcx, mxloop
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: a(nz, ny, nx), b(nz, ny, nx), ccoef(nz, ny, nx)
        real(c_double), intent(in) :: d(nz, ny, nx), e(nz, ny, nx), fcoef(nz, ny, nx)
        real(c_double), intent(in) :: g(nz, ny, nx), h(nz, ny, nx)
        real(c_double), intent(in) :: delz, dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, j, k, loop
        real(c_double) :: norm, norm_prev, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, delx_sqr

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        delx_sqr = delx * delx
        ratio_z = delx / delz
        ratio_y = delx / dely
        ratio_z_sqr = ratio_z * ratio_z
        ratio_y_sqr = ratio_y * ratio_y
        loop = 0

        do
            if (bcz == 1_c_int) call apply_extend_z3d(nz, ny, nx, s, undef)
            if (bcy == 1_c_int) call apply_extend_y3d(nz, ny, nx, s, undef, bcx)
            if (bcx == 1_c_int) call apply_extend_x3d(nz, ny, nx, s, undef)

            do k = 2, nz - 1
                do j = 2, ny - 1
                    if (bcx == 2_c_int) then
                        i = 1
                        call update_general3d_point(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, k, j, i, &
                            delx, delx_sqr, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end if

                    do i = 2, nx - 1
                        call update_general3d_point(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, k, j, i, &
                            delx, delx_sqr, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end do

                    if (bcx == 2_c_int) then
                        i = nx
                        call update_general3d_point(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, k, j, i, &
                            delx, delx_sqr, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, optarg, undef)
                    end if
                end do
            end do

            norm = abs_norm3d(nz, ny, nx, s, undef)
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

    end subroutine sor_general3d

    subroutine sor_biharmonic2d(ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, icoef, jcoef, &
        dely, delx, bcy, bcx, optarg, undef, mxloop, tolerance, overflow, relerr, loops)

        integer(c_int), intent(in) :: ny, nx, bcy, bcx, mxloop
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx)
        real(c_double), intent(in) :: d(ny, nx), e(ny, nx), fcoef(ny, nx)
        real(c_double), intent(in) :: g(ny, nx), h(ny, nx), icoef(ny, nx), jcoef(ny, nx)
        real(c_double), intent(in) :: dely, delx, optarg, undef, tolerance
        integer(c_int), intent(out) :: overflow, loops
        real(c_double), intent(out) :: relerr

        integer :: i, j, loop
        real(c_double) :: norm, norm_prev, ratio, ratio_sqr, ratio_qtr
        real(c_double) :: ratio_fourth, delx_sqr, delx_cubed, delx_fourth

        overflow = 0_c_int
        loops = 0_c_int
        relerr = 1.0_c_double
        norm_prev = huge(1.0_c_double)
        ratio = delx / dely
        ratio_sqr = ratio * ratio
        ratio_qtr = ratio / 4.0_c_double
        ratio_fourth = ratio_sqr * ratio_sqr
        delx_sqr = delx * delx
        delx_cubed = delx_sqr * delx
        delx_fourth = delx_sqr * delx_sqr
        loop = 0

        do
            if (bcy == 1_c_int) call apply_extend_biharmonic_y(ny, nx, s, undef, bcx)
            if (bcx == 1_c_int) call apply_extend_biharmonic_x(ny, nx, s, undef)

            do j = 3, ny - 2
                if (bcx == 2_c_int) then
                    i = 1
                    call update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
                        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, &
                        ratio, ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)
                    i = 2
                    call update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
                        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, &
                        ratio, ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)
                end if

                do i = 3, nx - 2
                    call update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
                        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, &
                        ratio, ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)
                end do

                if (bcx == 2_c_int) then
                    i = nx - 1
                    call update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
                        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, &
                        ratio, ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)
                    i = nx
                    call update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
                        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, &
                        ratio, ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)
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

    end subroutine sor_biharmonic2d

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

    subroutine update_standard3d_point(nz, ny, nx, s, a, b, ccoef, f, k, j, i, &
        delx_sqr, ratio_z_sqr, ratio_y_sqr, optarg, undef)

        integer, intent(in) :: nz, ny, nx, k, j, i
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: a(nz, ny, nx), b(nz, ny, nx), ccoef(nz, ny, nx), f(nz, ny, nx)
        real(c_double), intent(in) :: delx_sqr, ratio_z_sqr, ratio_y_sqr, optarg, undef

        integer :: im, ip
        real(c_double) :: temp, denom

        im = i - 1
        ip = i + 1
        if (i == 1) im = nx
        if (i == nx) ip = 1

        if (f(k, j, i) == undef) return
        if (a(k + 1, j, i) == undef .or. a(k, j, i) == undef) return
        if (b(k, j + 1, i) == undef .or. b(k, j, i) == undef) return
        if (ccoef(k, j, ip) == undef .or. ccoef(k, j, i) == undef) return

        temp = (a(k + 1, j, i) * (s(k + 1, j, i) - s(k, j, i)) &
            - a(k, j, i) * (s(k, j, i) - s(k - 1, j, i))) * ratio_z_sqr &
            + (b(k, j + 1, i) * (s(k, j + 1, i) - s(k, j, i)) &
            - b(k, j, i) * (s(k, j, i) - s(k, j - 1, i))) * ratio_y_sqr &
            + ccoef(k, j, ip) * (s(k, j, ip) - s(k, j, i)) &
            - ccoef(k, j, i) * (s(k, j, i) - s(k, j, im)) &
            - f(k, j, i) * delx_sqr

        denom = (a(k + 1, j, i) + a(k, j, i)) * ratio_z_sqr &
            + (b(k, j + 1, i) + b(k, j, i)) * ratio_y_sqr &
            + ccoef(k, j, ip) + ccoef(k, j, i)
        s(k, j, i) = s(k, j, i) + optarg * temp / denom

    end subroutine update_standard3d_point

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

    subroutine update_general3d_point(nz, ny, nx, s, a, b, ccoef, d, e, fcoef, g, h, k, j, i, &
        delx, delx_sqr, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, optarg, undef)

        integer, intent(in) :: nz, ny, nx, k, j, i
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: a(nz, ny, nx), b(nz, ny, nx), ccoef(nz, ny, nx)
        real(c_double), intent(in) :: d(nz, ny, nx), e(nz, ny, nx), fcoef(nz, ny, nx)
        real(c_double), intent(in) :: g(nz, ny, nx), h(nz, ny, nx)
        real(c_double), intent(in) :: delx, delx_sqr, ratio_z, ratio_y, ratio_z_sqr, ratio_y_sqr, optarg, undef

        integer :: im, ip
        real(c_double) :: temp, denom

        im = i - 1
        ip = i + 1
        if (i == 1) im = nx
        if (i == nx) ip = 1

        if (h(k, j, i) == undef .or. g(k, j, i) == undef) return
        if (a(k, j, i) == undef .or. b(k, j, i) == undef .or. ccoef(k, j, i) == undef) return
        if (d(k, j, i) == undef .or. e(k, j, i) == undef .or. fcoef(k, j, i) == undef) return

        temp = a(k, j, i) * ((s(k + 1, j, i) - s(k, j, i)) - (s(k, j, i) - s(k - 1, j, i))) &
            * ratio_z_sqr &
            + b(k, j, i) * ((s(k, j + 1, i) - s(k, j, i)) - (s(k, j, i) - s(k, j - 1, i))) &
            * ratio_y_sqr &
            + ccoef(k, j, i) * ((s(k, j, ip) - s(k, j, i)) - (s(k, j, i) - s(k, j, im))) &
            + (d(k, j, i) * (s(k + 1, j, i) - s(k - 1, j, i)) * ratio_z &
            + e(k, j, i) * (s(k, j + 1, i) - s(k, j - 1, i)) * ratio_y &
            + fcoef(k, j, i) * (s(k, j, ip) - s(k, j, im))) * delx / 2.0_c_double &
            + (g(k, j, i) * s(k, j, i) - h(k, j, i)) * delx_sqr

        denom = (a(k, j, i) * ratio_z_sqr + b(k, j, i) * ratio_y_sqr + ccoef(k, j, i)) &
            * 2.0_c_double - g(k, j, i) * delx_sqr
        s(k, j, i) = s(k, j, i) + optarg * temp / denom

    end subroutine update_general3d_point

    subroutine update_biharmonic2d_point(ny, nx, s, a, b, ccoef, d, e, fcoef, &
        g, h, icoef, jcoef, j, i, delx_sqr, delx_cubed, delx_fourth, ratio, &
        ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef)

        integer, intent(in) :: ny, nx, j, i
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: a(ny, nx), b(ny, nx), ccoef(ny, nx)
        real(c_double), intent(in) :: d(ny, nx), e(ny, nx), fcoef(ny, nx)
        real(c_double), intent(in) :: g(ny, nx), h(ny, nx), icoef(ny, nx), jcoef(ny, nx)
        real(c_double), intent(in) :: delx_sqr, delx_cubed, delx_fourth, ratio
        real(c_double), intent(in) :: ratio_sqr, ratio_qtr, ratio_fourth, optarg, undef

        integer :: im1, im2, ip1, ip2
        real(c_double) :: temp, denom

        im1 = i - 1
        im2 = i - 2
        ip1 = i + 1
        ip2 = i + 2
        if (i == 1) then
            im1 = nx
            im2 = nx - 1
        else if (i == 2) then
            im2 = nx
        end if
        if (i == nx) then
            ip1 = 1
            ip2 = 2
        else if (i == nx - 1) then
            ip2 = 1
        end if

        if (jcoef(j, i) == undef) return
        if (a(j, i) == undef .or. b(j, i) == undef .or. ccoef(j, i) == undef) return
        if (d(j, i) == undef .or. e(j, i) == undef .or. fcoef(j, i) == undef) return
        if (g(j, i) == undef .or. h(j, i) == undef .or. icoef(j, i) == undef) return

        temp = a(j, i) * (s(j + 2, i) - 4.0_c_double * s(j + 1, i) + 6.0_c_double * s(j, i) &
            - 4.0_c_double * s(j - 1, i) + s(j - 2, i)) * ratio_fourth &
            + b(j, i) * (s(j + 2, ip2) - 2.0_c_double * s(j + 2, i) + s(j + 2, im2) &
            - 2.0_c_double * s(j, ip2) + 4.0_c_double * s(j, i) - 2.0_c_double * s(j, im2) &
            + s(j - 2, ip2) - 2.0_c_double * s(j - 2, i) + s(j - 2, im2)) * ratio_sqr / 16.0_c_double &
            + ccoef(j, i) * (s(j, ip2) - 4.0_c_double * s(j, ip1) + 6.0_c_double * s(j, i) &
            - 4.0_c_double * s(j, im1) + s(j, im2)) &
            + d(j, i) * ((s(j + 1, i) - s(j, i)) - (s(j, i) - s(j - 1, i))) * ratio_sqr * delx_sqr &
            + e(j, i) * ((s(j + 1, ip1) - s(j - 1, ip1)) - (s(j + 1, im1) - s(j - 1, im1))) &
            * ratio_qtr * delx_sqr &
            + fcoef(j, i) * ((s(j, ip1) - s(j, i)) - (s(j, i) - s(j, im1))) * delx_sqr &
            + g(j, i) * (s(j + 1, i) - s(j - 1, i)) * delx_cubed * ratio / 2.0_c_double &
            + h(j, i) * (s(j, ip1) - s(j, im1)) * delx_cubed / 2.0_c_double &
            + (icoef(j, i) * s(j, i) - jcoef(j, i)) * delx_fourth

        denom = (a(j, i) * ratio_fourth + ccoef(j, i)) * 6.0_c_double &
            + b(j, i) * ratio_sqr / 4.0_c_double &
            - (d(j, i) * ratio_sqr + fcoef(j, i)) * 2.0_c_double * delx_sqr &
            + icoef(j, i) * delx_fourth
        s(j, i) = s(j, i) - optarg * temp / denom

    end subroutine update_biharmonic2d_point

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

    subroutine apply_extend_biharmonic_y(ny, nx, s, undef, bcx)
        integer, intent(in) :: ny, nx, bcx
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: undef
        integer :: i, istart, iend

        istart = 3
        iend = nx - 2
        if (bcx == 2_c_int) then
            istart = 1
            iend = nx
        end if

        do i = istart, iend
            if (s(3, i) /= undef) then
                s(1, i) = s(3, i)
                s(2, i) = s(3, i)
            end if
            if (s(ny - 2, i) /= undef) then
                s(ny, i) = s(ny - 2, i)
                s(ny - 1, i) = s(ny - 2, i)
            end if
        end do
    end subroutine apply_extend_biharmonic_y

    subroutine apply_extend_biharmonic_x(ny, nx, s, undef)
        integer, intent(in) :: ny, nx
        real(c_double), intent(inout) :: s(ny, nx)
        real(c_double), intent(in) :: undef
        integer :: j

        do j = 3, ny - 2
            if (s(j, 3) /= undef) then
                s(j, 1) = s(j, 3)
                s(j, 2) = s(j, 3)
            end if
            if (s(j, nx - 2) /= undef) then
                s(j, nx) = s(j, nx - 2)
                s(j, nx - 1) = s(j, nx - 2)
            end if
        end do
    end subroutine apply_extend_biharmonic_x

    subroutine apply_extend_z3d(nz, ny, nx, s, undef)
        integer, intent(in) :: nz, ny, nx
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: undef
        integer :: i, j

        do i = 1, nx
            do j = 1, ny
                if (s(2, j, i) /= undef) s(1, j, i) = s(2, j, i)
                if (s(nz - 1, j, i) /= undef) s(nz, j, i) = s(nz - 1, j, i)
            end do
        end do
    end subroutine apply_extend_z3d

    subroutine apply_extend_y3d(nz, ny, nx, s, undef, bcx)
        integer, intent(in) :: nz, ny, nx, bcx
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: undef
        integer :: i, k, istart, iend

        istart = 2
        iend = nx - 1
        if (bcx == 2_c_int) then
            istart = 1
            iend = nx
        end if

        do k = 2, nz - 1
            do i = istart, iend
                if (s(k, 2, i) /= undef) s(k, 1, i) = s(k, 2, i)
                if (s(k, ny - 1, i) /= undef) s(k, ny, i) = s(k, ny - 1, i)
            end do
        end do
    end subroutine apply_extend_y3d

    subroutine apply_extend_x3d(nz, ny, nx, s, undef)
        integer, intent(in) :: nz, ny, nx
        real(c_double), intent(inout) :: s(nz, ny, nx)
        real(c_double), intent(in) :: undef
        integer :: j, k

        do k = 2, nz - 1
            do j = 2, ny - 1
                if (s(k, j, 2) /= undef) s(k, j, 1) = s(k, j, 2)
                if (s(k, j, nx - 1) /= undef) s(k, j, nx) = s(k, j, nx - 1)
            end do
        end do
    end subroutine apply_extend_x3d

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

    real(c_double) function abs_norm3d(nz, ny, nx, s, undef) result(norm)
        integer, intent(in) :: nz, ny, nx
        real(c_double), intent(in) :: s(nz, ny, nx), undef
        integer :: i, j, k, count

        norm = 0.0_c_double
        count = 0
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    if (s(k, j, i) /= undef) then
                        norm = norm + abs(s(k, j, i))
                        count = count + 1
                    end if
                end do
            end do
        end do
        if (count > 0) then
            norm = norm / real(count, c_double)
        else
            norm = huge(1.0_c_double)
        end if
    end function abs_norm3d

end module iterative_solvers
