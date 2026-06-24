module fftpack_c_api

    use, intrinsic :: iso_c_binding, only: &
        c_double, &
        c_int, &
        c_ptr, &
        c_f_pointer
    use type_PeriodicFastFourierTransform, only: &
        rffti, rfftf, rfftb, &
        sinti, sint, costi, cost, &
        sinqi, sinqf, sinqb, cosqi, cosqf, cosqb, &
        cffti, cfftf, cfftb

    implicit none

contains

    subroutine pyfp_fftpack_real_transform(n, x_ptr, kind, direction) &
        bind(C, name="pyfp_fftpack_real_transform")

        integer(c_int), value :: n, kind, direction
        type(c_ptr), value :: x_ptr
        real(c_double), pointer :: x(:)
        real(c_double), allocatable :: wsave(:)

        call c_f_pointer(x_ptr, x, [n])

        select case (kind)
        case (0_c_int)
            allocate(wsave(2*n + 15))
            call rffti(n, wsave)
            if (direction == 0_c_int) then
                call rfftf(n, x, wsave)
            else
                call rfftb(n, x, wsave)
            end if
        case (1_c_int)
            allocate(wsave(3*n + 15))
            call sinti(n, wsave)
            call sint(n, x, wsave)
        case (2_c_int)
            allocate(wsave(3*n + 15))
            call costi(n, wsave)
            call cost(n, x, wsave)
        case (3_c_int)
            allocate(wsave(3*n + 15))
            call sinqi(n, wsave)
            if (direction == 0_c_int) then
                call sinqf(n, x, wsave)
            else
                call sinqb(n, x, wsave)
            end if
        case (4_c_int)
            allocate(wsave(3*n + 15))
            call cosqi(n, wsave)
            if (direction == 0_c_int) then
                call cosqf(n, x, wsave)
            else
                call cosqb(n, x, wsave)
            end if
        case default
            return
        end select

    end subroutine pyfp_fftpack_real_transform

    subroutine pyfp_fftpack_complex_transform(n, cdata_ptr, direction) &
        bind(C, name="pyfp_fftpack_complex_transform")

        integer(c_int), value :: n, direction
        type(c_ptr), value :: cdata_ptr
        real(c_double), pointer :: c(:)
        real(c_double), allocatable :: wsave(:)

        call c_f_pointer(cdata_ptr, c, [2*n])
        allocate(wsave(4*n + 15))
        call cffti(n, wsave)
        if (direction == 0_c_int) then
            call cfftf(n, c, wsave)
        else
            call cfftb(n, c, wsave)
        end if

    end subroutine pyfp_fftpack_complex_transform

end module fftpack_c_api
