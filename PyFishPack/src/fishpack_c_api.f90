module fishpack_c_api

    use, intrinsic :: iso_c_binding, only: &
        c_double, &
        c_int, &
        c_ptr, &
        c_f_pointer

    use centered_real_linear_systems_solver, only: &
        genbun

    use centered_helmholtz_solvers, only: &
        hwscrt, &
        hwsplr, &
        hwscyl, &
        hwsssp, &
        hwscsp

    use staggered_helmholtz_solvers, only: &
        hstcrt, &
        hstplr, &
        hstcyl, &
        hstssp, &
        hstcsp

    use staggered_real_linear_systems_solver, only: &
        poistg

    use three_dimensional_solvers, only: &
        pois3d, &
        hw3crt

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    implicit none

contains

    subroutine pyfp_genbun(nperod, n, mperod, m, idimy, a_ptr, b_ptr, ccoef_ptr, y_ptr, ierror) &
        bind(C, name="pyfp_genbun")

        integer(c_int), value      :: nperod, n, mperod, m, idimy
        type(c_ptr),    value      :: a_ptr, b_ptr, ccoef_ptr, y_ptr
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: a(:), b(:), c(:), y(:, :)

        call c_f_pointer(a_ptr, a, [m])
        call c_f_pointer(b_ptr, b, [m])
        call c_f_pointer(ccoef_ptr, c, [m])
        call c_f_pointer(y_ptr, y, [idimy, n])

        call genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)

    end subroutine pyfp_genbun

    subroutine pyfp_poistg(nperod, n, mperod, m, idimy, a_ptr, b_ptr, ccoef_ptr, y_ptr, ierror) &
        bind(C, name="pyfp_poistg")

        integer(c_int), value      :: nperod, n, mperod, m, idimy
        type(c_ptr),    value      :: a_ptr, b_ptr, ccoef_ptr, y_ptr
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: a(:), b(:), c(:), y(:, :)

        call c_f_pointer(a_ptr, a, [m])
        call c_f_pointer(b_ptr, b, [m])
        call c_f_pointer(ccoef_ptr, c, [m])
        call c_f_pointer(y_ptr, y, [idimy, n])

        call poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)

    end subroutine pyfp_poistg

    subroutine pyfp_pois3d(lperod, l, c1, mperod, m, c2, nperod, n, ldimf, mdimf, &
        a_ptr, b_ptr, ccoef_ptr, f_ptr, ierror) bind(C, name="pyfp_pois3d")

        integer(c_int), value      :: lperod, l, mperod, m, nperod, n, ldimf, mdimf
        real(c_double), value      :: c1, c2
        type(c_ptr),    value      :: a_ptr, b_ptr, ccoef_ptr, f_ptr
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: a(:), b(:), c(:), f(:, :, :)

        call c_f_pointer(a_ptr, a, [n])
        call c_f_pointer(b_ptr, b, [n])
        call c_f_pointer(ccoef_ptr, c, [n])
        call c_f_pointer(f_ptr, f, [ldimf, mdimf, n])

        call pois3d(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, ldimf, mdimf, f, ierror)

    end subroutine pyfp_pois3d

    subroutine pyfp_hwscrt(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hwscrt")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n + 1])

        call hwscrt(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hwscrt

    subroutine pyfp_hstcrt(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hstcrt")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n])

        call hstcrt(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hstcrt

    subroutine pyfp_hwsplr(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hwsplr")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n + 1])

        call hwsplr(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hwsplr

    subroutine pyfp_hwscyl(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hwscyl")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n + 1])

        call hwscyl(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hwscyl

    subroutine pyfp_hstplr(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hstplr")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n])

        call hstplr(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hstplr

    subroutine pyfp_hstcyl(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hstcyl")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n])

        call hstcyl(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hstcyl

    subroutine pyfp_hwsssp(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hwsssp")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n + 1])

        call hwsssp(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hwsssp

    subroutine pyfp_hstssp(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hstssp")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n])

        call hstssp(a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    end subroutine pyfp_hstssp

    subroutine pyfp_hwscsp(ts, tf, m, mbdcnd, bdts_ptr, len_bdts, bdtf_ptr, len_bdtf, &
        rs, rf, n, nbdcnd, bdrs_ptr, len_bdrs, bdrf_ptr, len_bdrf, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hwscsp")

        real(c_double), value      :: ts, tf, rs, rf, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bdts, len_bdtf, n, nbdcnd, len_bdrs, len_bdrf, idimf
        type(c_ptr),    value      :: bdts_ptr, bdtf_ptr, bdrs_ptr, bdrf_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        integer(c_int) :: intl
        real(c_double), pointer :: bdts(:), bdtf(:), bdrs(:), bdrf(:), f(:, :)
        type(FishpackWorkspace) :: workspace

        call c_f_pointer(bdts_ptr, bdts, [len_bdts])
        call c_f_pointer(bdtf_ptr, bdtf, [len_bdtf])
        call c_f_pointer(bdrs_ptr, bdrs, [len_bdrs])
        call c_f_pointer(bdrf_ptr, bdrf, [len_bdrf])
        call c_f_pointer(f_ptr, f, [idimf, n + 1])

        intl = 0
        call hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
            nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, workspace)
        call workspace%destroy()

    end subroutine pyfp_hwscsp

    subroutine pyfp_hstcsp(a_lower, b_upper, m, mbdcnd, bda_ptr, len_bda, bdb_ptr, len_bdb, &
        c_lower, d_upper, n, nbdcnd, bdc_ptr, len_bdc, bdd_ptr, len_bdd, elmbda, &
        idimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hstcsp")

        real(c_double), value      :: a_lower, b_upper, c_lower, d_upper, elmbda
        integer(c_int), value      :: m, mbdcnd, len_bda, len_bdb, n, nbdcnd, len_bdc, len_bdd, idimf
        type(c_ptr),    value      :: bda_ptr, bdb_ptr, bdc_ptr, bdd_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        integer(c_int) :: intl
        real(c_double), pointer :: bda(:), bdb(:), bdc(:), bdd(:), f(:, :)
        type(FishpackWorkspace) :: workspace

        call c_f_pointer(bda_ptr, bda, [len_bda])
        call c_f_pointer(bdb_ptr, bdb, [len_bdb])
        call c_f_pointer(bdc_ptr, bdc, [len_bdc])
        call c_f_pointer(bdd_ptr, bdd, [len_bdd])
        call c_f_pointer(f_ptr, f, [idimf, n])

        intl = 0
        call hstcsp(intl, a_lower, b_upper, m, mbdcnd, bda, bdb, c_lower, d_upper, n, &
            nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, workspace)
        call workspace%destroy()

    end subroutine pyfp_hstcsp

    subroutine pyfp_hw3crt(xs, xf, l, lbdcnd, bdxs_ptr, bdxs_dim1, bdxs_dim2, &
        bdxf_ptr, bdxf_dim1, bdxf_dim2, ys, yf, m, mbdcnd, bdys_ptr, bdys_dim1, &
        bdys_dim2, bdyf_ptr, bdyf_dim1, bdyf_dim2, zs, zf, n, nbdcnd, bdzs_ptr, &
        bdzs_dim1, bdzs_dim2, bdzf_ptr, bdzf_dim1, bdzf_dim2, elmbda, ldimf, &
        mdimf, f_ptr, pertrb, ierror) bind(C, name="pyfp_hw3crt")

        real(c_double), value      :: xs, xf, ys, yf, zs, zf, elmbda
        integer(c_int), value      :: l, lbdcnd, m, mbdcnd, n, nbdcnd, ldimf, mdimf
        integer(c_int), value      :: bdxs_dim1, bdxs_dim2, bdxf_dim1, bdxf_dim2
        integer(c_int), value      :: bdys_dim1, bdys_dim2, bdyf_dim1, bdyf_dim2
        integer(c_int), value      :: bdzs_dim1, bdzs_dim2, bdzf_dim1, bdzf_dim2
        type(c_ptr),    value      :: bdxs_ptr, bdxf_ptr, bdys_ptr, bdyf_ptr, bdzs_ptr, bdzf_ptr, f_ptr
        real(c_double), intent(out) :: pertrb
        integer(c_int), intent(out) :: ierror

        real(c_double), pointer :: bdxs(:, :), bdxf(:, :), bdys(:, :), bdyf(:, :)
        real(c_double), pointer :: bdzs(:, :), bdzf(:, :), f(:, :, :)

        call c_f_pointer(bdxs_ptr, bdxs, [bdxs_dim1, bdxs_dim2])
        call c_f_pointer(bdxf_ptr, bdxf, [bdxf_dim1, bdxf_dim2])
        call c_f_pointer(bdys_ptr, bdys, [bdys_dim1, bdys_dim2])
        call c_f_pointer(bdyf_ptr, bdyf, [bdyf_dim1, bdyf_dim2])
        call c_f_pointer(bdzs_ptr, bdzs, [bdzs_dim1, bdzs_dim2])
        call c_f_pointer(bdzf_ptr, bdzf, [bdzf_dim1, bdzf_dim2])
        call c_f_pointer(f_ptr, f, [ldimf, mdimf, n + 1])

        call hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, bdys, &
            bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, mdimf, f, &
            pertrb, ierror)

    end subroutine pyfp_hw3crt

end module fishpack_c_api
