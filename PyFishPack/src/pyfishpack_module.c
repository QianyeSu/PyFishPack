#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_20_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>

extern void pyfp_genbun(int nperod, int n, int mperod, int m, int idimy,
                        void *a, void *b, void *c, void *y, int *ierror);
extern void pyfp_poistg(int nperod, int n, int mperod, int m, int idimy,
                        void *a, void *b, void *c, void *y, int *ierror);
extern void pyfp_pois3d(int lperod, int l, double c1, int mperod, int m,
                        double c2, int nperod, int n, int ldimf, int mdimf,
                        void *a, void *b, void *c, void *f, int *ierror);
extern void pyfp_hwscrt(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hstcrt(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hwsplr(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hwscyl(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hstplr(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hstcyl(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hwsssp(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hstssp(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hwscsp(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hstcsp(double a_lower, double b_upper, int m, int mbdcnd,
                        void *bda, int len_bda, void *bdb, int len_bdb,
                        double c_lower, double d_upper, int n, int nbdcnd,
                        void *bdc, int len_bdc, void *bdd, int len_bdd,
                        double elmbda, int idimf, void *f, double *pertrb,
                        int *ierror);
extern void pyfp_hw3crt(double xs, double xf, int l, int lbdcnd,
                        void *bdxs, int bdxs_dim1, int bdxs_dim2,
                        void *bdxf, int bdxf_dim1, int bdxf_dim2,
                        double ys, double yf, int m, int mbdcnd,
                        void *bdys, int bdys_dim1, int bdys_dim2,
                        void *bdyf, int bdyf_dim1, int bdyf_dim2,
                        double zs, double zf, int n, int nbdcnd,
                        void *bdzs, int bdzs_dim1, int bdzs_dim2,
                        void *bdzf, int bdzf_dim1, int bdzf_dim2,
                        double elmbda, int ldimf, int mdimf, void *f,
                        double *pertrb, int *ierror);
extern void pyfp_sor_standard1d(int n, void *s, void *a, void *b, void *f,
                                double delx, int bcx, double optarg, double undef,
                                int mxloop, double tolerance, int *overflow,
                                double *relerr, int *loops);
extern void pyfp_sor_standard2d(int ny, int nx, void *s, void *a, void *b,
                                void *c, void *f, double dely, double delx,
                                int bcy, int bcx, double optarg, double undef,
                                int mxloop, double tolerance, int *overflow,
                                double *relerr, int *loops);
extern void pyfp_sor_standard3d(int nz, int ny, int nx, void *s, void *a,
                                void *b, void *c, void *f, double delz,
                                double dely, double delx, int bcz, int bcy,
                                int bcx, double optarg, double undef,
                                int mxloop, double tolerance, int *overflow,
                                double *relerr, int *loops);
extern void pyfp_sor_general2d(int ny, int nx, void *s, void *a, void *b,
                               void *c, void *d, void *e, void *fcoef, void *g,
                               double dely, double delx, int bcy, int bcx,
                               double optarg, double undef, int mxloop,
                               double tolerance, int *overflow, double *relerr,
                               int *loops);
extern void pyfp_sor_general3d(int nz, int ny, int nx, void *s, void *a,
                               void *b, void *c, void *d, void *e,
                               void *fcoef, void *g, void *h, double delz,
                               double dely, double delx, int bcz, int bcy,
                               int bcx, double optarg, double undef,
                               int mxloop, double tolerance, int *overflow,
                               double *relerr, int *loops);
extern void pyfp_sor_biharmonic2d(int ny, int nx, void *s, void *a, void *b,
                                  void *c, void *d, void *e, void *fcoef,
                                  void *g, void *h, void *icoef, void *jcoef,
                                  double dely, double delx, int bcy, int bcx,
                                  double optarg, double undef, int mxloop,
                                  double tolerance, int *overflow,
                                  double *relerr, int *loops);
extern void pyfp_fftpack_real_transform(int n, void *x, int kind, int direction);
extern void pyfp_fftpack_complex_transform(int n, void *c, int direction);

static PyArrayObject *as_double_1d(PyObject *obj, const char *name)
{
    PyArrayObject *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr == NULL) {
        return NULL;
    }
    if (PyArray_NDIM(arr) != 1) {
        PyErr_Format(PyExc_ValueError, "%s must be a one-dimensional float64 array", name);
        Py_DECREF(arr);
        return NULL;
    }
    return arr;
}

static PyArrayObject *as_fortran_double_copy(PyObject *obj, int ndim, const char *name)
{
    PyArray_Descr *descr = PyArray_DescrFromType(NPY_DOUBLE);
    PyArrayObject *arr = (PyArrayObject *)PyArray_FromAny(
        obj,
        descr,
        ndim,
        ndim,
        NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE | NPY_ARRAY_ENSURECOPY,
        NULL);
    if (arr == NULL) {
        return NULL;
    }
    if (!PyArray_ISFARRAY(arr)) {
        PyErr_Format(PyExc_ValueError, "%s could not be converted to a Fortran-contiguous float64 array", name);
        Py_DECREF(arr);
        return NULL;
    }
    return arr;
}

static int require_len(PyArrayObject *arr, npy_intp expected, const char *name)
{
    if (PyArray_DIM(arr, 0) < expected) {
        PyErr_Format(PyExc_ValueError, "%s has length %zd, expected at least %zd",
                     name, (Py_ssize_t)PyArray_DIM(arr, 0), (Py_ssize_t)expected);
        return 0;
    }
    return 1;
}

static int require_2d_shape(PyArrayObject *arr, npy_intp dim0, npy_intp dim1, const char *name)
{
    if (PyArray_DIM(arr, 0) < dim0 || PyArray_DIM(arr, 1) < dim1) {
        PyErr_Format(PyExc_ValueError, "%s has shape (%zd, %zd), expected at least (%zd, %zd)",
                     name,
                     (Py_ssize_t)PyArray_DIM(arr, 0), (Py_ssize_t)PyArray_DIM(arr, 1),
                     (Py_ssize_t)dim0, (Py_ssize_t)dim1);
        return 0;
    }
    return 1;
}

static int require_3d_shape(PyArrayObject *arr, npy_intp dim0, npy_intp dim1, npy_intp dim2, const char *name)
{
    if (PyArray_DIM(arr, 0) < dim0 || PyArray_DIM(arr, 1) < dim1 ||
        PyArray_DIM(arr, 2) < dim2) {
        PyErr_Format(PyExc_ValueError, "%s has shape (%zd, %zd, %zd), expected at least (%zd, %zd, %zd)",
                     name,
                     (Py_ssize_t)PyArray_DIM(arr, 0), (Py_ssize_t)PyArray_DIM(arr, 1),
                     (Py_ssize_t)PyArray_DIM(arr, 2),
                     (Py_ssize_t)dim0, (Py_ssize_t)dim1, (Py_ssize_t)dim2);
        return 0;
    }
    return 1;
}

static int bc_code_from_string(const char *bc, const char *name)
{
    if (strcmp(bc, "fixed") == 0) {
        return 0;
    }
    if (strcmp(bc, "extend") == 0) {
        return 1;
    }
    if (strcmp(bc, "periodic") == 0) {
        return 2;
    }
    PyErr_Format(PyExc_NotImplementedError,
                 "%s must be one of 'fixed', 'extend', or 'periodic'", name);
    return -1;
}

typedef void (*helmholtz_2d_solver)(double, double, int, int, void *, int, void *, int,
                                    double, double, int, int, void *, int, void *, int,
                                    double, int, void *, double *, int *);

static PyObject *fishpack_helmholtz_2d(PyObject *args, helmholtz_2d_solver solver, int centered_grid)
{
    int m, mbdcnd, n, nbdcnd;
    int idimf, ierror = 0;
    double a_lower, b_upper, c_lower, d_upper, elmbda, pertrb = 0.0;
    PyObject *bda_obj, *bdb_obj, *bdc_obj, *bdd_obj, *f_obj;
    PyArrayObject *bda = NULL, *bdb = NULL, *bdc = NULL, *bdd = NULL, *f = NULL;

    if (!PyArg_ParseTuple(args, "ddiiOOddiiOOdO",
                          &a_lower, &b_upper, &m, &mbdcnd, &bda_obj, &bdb_obj,
                          &c_lower, &d_upper, &n, &nbdcnd, &bdc_obj, &bdd_obj,
                          &elmbda, &f_obj)) {
        return NULL;
    }

    bda = as_double_1d(bda_obj, "bda");
    bdb = as_double_1d(bdb_obj, "bdb");
    bdc = as_double_1d(bdc_obj, "bdc");
    bdd = as_double_1d(bdd_obj, "bdd");
    f = as_fortran_double_copy(f_obj, 2, "f");
    if (bda == NULL || bdb == NULL || bdc == NULL || bdd == NULL || f == NULL) {
        goto fail;
    }

    idimf = (int)PyArray_DIM(f, 0);
    const int min_rows = centered_grid ? m + 1 : m;
    const int min_cols = centered_grid ? n + 1 : n;
    if (idimf < min_rows || PyArray_DIM(f, 1) < min_cols) {
        PyErr_Format(PyExc_ValueError,
                     "f must have shape (idimf, n%s) with idimf >= m%s; got (%zd, %zd), m=%d, n=%d",
                     centered_grid ? " + 1" : "",
                     centered_grid ? " + 1" : "",
                     (Py_ssize_t)PyArray_DIM(f, 0), (Py_ssize_t)PyArray_DIM(f, 1), m, n);
        goto fail;
    }
    if (PyArray_DIM(bda, 0) < 1 || PyArray_DIM(bdb, 0) < 1 ||
        PyArray_DIM(bdc, 0) < 1 || PyArray_DIM(bdd, 0) < 1) {
        PyErr_SetString(PyExc_ValueError, "boundary arrays must have at least one element");
        goto fail;
    }

    solver(a_lower, b_upper, m, mbdcnd,
           PyArray_DATA(bda), (int)PyArray_DIM(bda, 0),
           PyArray_DATA(bdb), (int)PyArray_DIM(bdb, 0),
           c_lower, d_upper, n, nbdcnd,
           PyArray_DATA(bdc), (int)PyArray_DIM(bdc, 0),
           PyArray_DATA(bdd), (int)PyArray_DIM(bdd, 0),
           elmbda, idimf, PyArray_DATA(f), &pertrb, &ierror);

    Py_DECREF(bda);
    Py_DECREF(bdb);
    Py_DECREF(bdc);
    Py_DECREF(bdd);
    return Py_BuildValue("Ndi", (PyObject *)f, pertrb, ierror);

fail:
    Py_XDECREF(bda);
    Py_XDECREF(bdb);
    Py_XDECREF(bdc);
    Py_XDECREF(bdd);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_backend_info(PyObject *self, PyObject *Py_UNUSED(args))
{
    return Py_BuildValue(
        "{s:s,s:s,s:s}",
        "backend", "modern-fortran",
        "abi", "fortran-bind-c",
        "wrapper", "numpy-c-api");
}

static PyObject *fishpack_genbun(PyObject *self, PyObject *args)
{
    int nperod, n, mperod, m;
    int idimy, ierror = 0;
    PyObject *a_obj, *b_obj, *c_obj, *y_obj;
    PyArrayObject *a = NULL, *b = NULL, *c = NULL, *y = NULL;

    if (!PyArg_ParseTuple(args, "iiiiOOOO:genbun", &nperod, &n, &mperod, &m,
                          &a_obj, &b_obj, &c_obj, &y_obj)) {
        return NULL;
    }

    a = as_double_1d(a_obj, "a");
    b = as_double_1d(b_obj, "b");
    c = as_double_1d(c_obj, "c");
    y = as_fortran_double_copy(y_obj, 2, "y");
    if (a == NULL || b == NULL || c == NULL || y == NULL) {
        goto fail;
    }

    idimy = (int)PyArray_DIM(y, 0);
    if (!require_len(a, m, "a") || !require_len(b, m, "b") || !require_len(c, m, "c")) {
        goto fail;
    }
    if (PyArray_DIM(y, 1) != n || idimy < m) {
        PyErr_Format(PyExc_ValueError,
                     "y must have shape (idimy, n) with idimy >= m; got (%zd, %zd), m=%d, n=%d",
                     (Py_ssize_t)PyArray_DIM(y, 0), (Py_ssize_t)PyArray_DIM(y, 1), m, n);
        goto fail;
    }

    pyfp_genbun(nperod, n, mperod, m, idimy,
                PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c), PyArray_DATA(y), &ierror);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    return Py_BuildValue("Ni", (PyObject *)y, ierror);

fail:
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(y);
    return NULL;
}

static PyObject *fishpack_poistg(PyObject *self, PyObject *args)
{
    int nperod, n, mperod, m;
    int idimy, ierror = 0;
    PyObject *a_obj, *b_obj, *c_obj, *y_obj;
    PyArrayObject *a = NULL, *b = NULL, *c = NULL, *y = NULL;

    if (!PyArg_ParseTuple(args, "iiiiOOOO:poistg", &nperod, &n, &mperod, &m,
                          &a_obj, &b_obj, &c_obj, &y_obj)) {
        return NULL;
    }

    a = as_double_1d(a_obj, "a");
    b = as_double_1d(b_obj, "b");
    c = as_double_1d(c_obj, "c");
    y = as_fortran_double_copy(y_obj, 2, "y");
    if (a == NULL || b == NULL || c == NULL || y == NULL) {
        goto fail;
    }

    idimy = (int)PyArray_DIM(y, 0);
    if (!require_len(a, m, "a") || !require_len(b, m, "b") || !require_len(c, m, "c")) {
        goto fail;
    }
    if (PyArray_DIM(y, 1) != n || idimy < m) {
        PyErr_Format(PyExc_ValueError,
                     "y must have shape (idimy, n) with idimy >= m; got (%zd, %zd), m=%d, n=%d",
                     (Py_ssize_t)PyArray_DIM(y, 0), (Py_ssize_t)PyArray_DIM(y, 1), m, n);
        goto fail;
    }

    pyfp_poistg(nperod, n, mperod, m, idimy,
                PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c), PyArray_DATA(y), &ierror);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    return Py_BuildValue("Ni", (PyObject *)y, ierror);

fail:
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(y);
    return NULL;
}

static PyObject *fishpack_pois3d(PyObject *self, PyObject *args)
{
    int lperod, l, mperod, m, nperod, n;
    int ldimf, mdimf, ierror = 0;
    double c1, c2;
    PyObject *a_obj, *b_obj, *c_obj, *f_obj;
    PyArrayObject *a = NULL, *b = NULL, *c = NULL, *f = NULL;

    if (!PyArg_ParseTuple(args, "iidiidiiOOOO:pois3d", &lperod, &l, &c1,
                          &mperod, &m, &c2, &nperod, &n,
                          &a_obj, &b_obj, &c_obj, &f_obj)) {
        return NULL;
    }

    a = as_double_1d(a_obj, "a");
    b = as_double_1d(b_obj, "b");
    c = as_double_1d(c_obj, "c");
    f = as_fortran_double_copy(f_obj, 3, "f");
    if (a == NULL || b == NULL || c == NULL || f == NULL) {
        goto fail;
    }

    ldimf = (int)PyArray_DIM(f, 0);
    mdimf = (int)PyArray_DIM(f, 1);
    if (!require_len(a, n, "a") || !require_len(b, n, "b") || !require_len(c, n, "c")) {
        goto fail;
    }
    if (PyArray_DIM(f, 2) != n || ldimf < l || mdimf < m) {
        PyErr_Format(PyExc_ValueError,
                     "f must have shape (ldimf, mdimf, n) with ldimf >= l and mdimf >= m; got (%zd, %zd, %zd)",
                     (Py_ssize_t)PyArray_DIM(f, 0), (Py_ssize_t)PyArray_DIM(f, 1),
                     (Py_ssize_t)PyArray_DIM(f, 2));
        goto fail;
    }

    pyfp_pois3d(lperod, l, c1, mperod, m, c2, nperod, n, ldimf, mdimf,
                PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c), PyArray_DATA(f), &ierror);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    return Py_BuildValue("Ni", (PyObject *)f, ierror);

fail:
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_hwscrt(PyObject *self, PyObject *args)
{
    int m, mbdcnd, n, nbdcnd;
    int idimf, ierror = 0;
    double a_lower, b_upper, c_lower, d_upper, elmbda, pertrb = 0.0;
    PyObject *bda_obj, *bdb_obj, *bdc_obj, *bdd_obj, *f_obj;
    PyArrayObject *bda = NULL, *bdb = NULL, *bdc = NULL, *bdd = NULL, *f = NULL;

    if (!PyArg_ParseTuple(args, "ddiiOOddiiOOdO:hwscrt",
                          &a_lower, &b_upper, &m, &mbdcnd, &bda_obj, &bdb_obj,
                          &c_lower, &d_upper, &n, &nbdcnd, &bdc_obj, &bdd_obj,
                          &elmbda, &f_obj)) {
        return NULL;
    }

    bda = as_double_1d(bda_obj, "bda");
    bdb = as_double_1d(bdb_obj, "bdb");
    bdc = as_double_1d(bdc_obj, "bdc");
    bdd = as_double_1d(bdd_obj, "bdd");
    f = as_fortran_double_copy(f_obj, 2, "f");
    if (bda == NULL || bdb == NULL || bdc == NULL || bdd == NULL || f == NULL) {
        goto fail;
    }

    idimf = (int)PyArray_DIM(f, 0);
    if (idimf < m + 1 || PyArray_DIM(f, 1) < n + 1) {
        PyErr_Format(PyExc_ValueError,
                     "f must have shape (idimf, n + 1) with idimf >= m + 1; got (%zd, %zd), m=%d, n=%d",
                     (Py_ssize_t)PyArray_DIM(f, 0), (Py_ssize_t)PyArray_DIM(f, 1), m, n);
        goto fail;
    }
    if (PyArray_DIM(bda, 0) < 1 || PyArray_DIM(bdb, 0) < 1 ||
        PyArray_DIM(bdc, 0) < 1 || PyArray_DIM(bdd, 0) < 1) {
        PyErr_SetString(PyExc_ValueError, "boundary arrays must have at least one element");
        goto fail;
    }

    pyfp_hwscrt(a_lower, b_upper, m, mbdcnd,
                PyArray_DATA(bda), (int)PyArray_DIM(bda, 0),
                PyArray_DATA(bdb), (int)PyArray_DIM(bdb, 0),
                c_lower, d_upper, n, nbdcnd,
                PyArray_DATA(bdc), (int)PyArray_DIM(bdc, 0),
                PyArray_DATA(bdd), (int)PyArray_DIM(bdd, 0),
                elmbda, idimf, PyArray_DATA(f), &pertrb, &ierror);

    Py_DECREF(bda);
    Py_DECREF(bdb);
    Py_DECREF(bdc);
    Py_DECREF(bdd);
    return Py_BuildValue("Ndi", (PyObject *)f, pertrb, ierror);

fail:
    Py_XDECREF(bda);
    Py_XDECREF(bdb);
    Py_XDECREF(bdc);
    Py_XDECREF(bdd);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_hstcrt(PyObject *self, PyObject *args)
{
    int m, mbdcnd, n, nbdcnd;
    int idimf, ierror = 0;
    double a_lower, b_upper, c_lower, d_upper, elmbda, pertrb = 0.0;
    PyObject *bda_obj, *bdb_obj, *bdc_obj, *bdd_obj, *f_obj;
    PyArrayObject *bda = NULL, *bdb = NULL, *bdc = NULL, *bdd = NULL, *f = NULL;

    if (!PyArg_ParseTuple(args, "ddiiOOddiiOOdO:hstcrt",
                          &a_lower, &b_upper, &m, &mbdcnd, &bda_obj, &bdb_obj,
                          &c_lower, &d_upper, &n, &nbdcnd, &bdc_obj, &bdd_obj,
                          &elmbda, &f_obj)) {
        return NULL;
    }

    bda = as_double_1d(bda_obj, "bda");
    bdb = as_double_1d(bdb_obj, "bdb");
    bdc = as_double_1d(bdc_obj, "bdc");
    bdd = as_double_1d(bdd_obj, "bdd");
    f = as_fortran_double_copy(f_obj, 2, "f");
    if (bda == NULL || bdb == NULL || bdc == NULL || bdd == NULL || f == NULL) {
        goto fail;
    }

    idimf = (int)PyArray_DIM(f, 0);
    if (idimf < m || PyArray_DIM(f, 1) < n) {
        PyErr_Format(PyExc_ValueError,
                     "f must have shape (idimf, n) with idimf >= m; got (%zd, %zd), m=%d, n=%d",
                     (Py_ssize_t)PyArray_DIM(f, 0), (Py_ssize_t)PyArray_DIM(f, 1), m, n);
        goto fail;
    }
    if (PyArray_DIM(bda, 0) < 1 || PyArray_DIM(bdb, 0) < 1 ||
        PyArray_DIM(bdc, 0) < 1 || PyArray_DIM(bdd, 0) < 1) {
        PyErr_SetString(PyExc_ValueError, "boundary arrays must have at least one element");
        goto fail;
    }

    pyfp_hstcrt(a_lower, b_upper, m, mbdcnd,
                PyArray_DATA(bda), (int)PyArray_DIM(bda, 0),
                PyArray_DATA(bdb), (int)PyArray_DIM(bdb, 0),
                c_lower, d_upper, n, nbdcnd,
                PyArray_DATA(bdc), (int)PyArray_DIM(bdc, 0),
                PyArray_DATA(bdd), (int)PyArray_DIM(bdd, 0),
                elmbda, idimf, PyArray_DATA(f), &pertrb, &ierror);

    Py_DECREF(bda);
    Py_DECREF(bdb);
    Py_DECREF(bdc);
    Py_DECREF(bdd);
    return Py_BuildValue("Ndi", (PyObject *)f, pertrb, ierror);

fail:
    Py_XDECREF(bda);
    Py_XDECREF(bdb);
    Py_XDECREF(bdc);
    Py_XDECREF(bdd);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_hw3crt(PyObject *self, PyObject *args)
{
    int l, lbdcnd, m, mbdcnd, n, nbdcnd;
    int ldimf, mdimf, ierror = 0;
    double xs, xf, ys, yf, zs, zf, elmbda, pertrb = 0.0;
    PyObject *bdxs_obj, *bdxf_obj, *bdys_obj, *bdyf_obj, *bdzs_obj, *bdzf_obj, *f_obj;
    PyArrayObject *bdxs = NULL, *bdxf = NULL, *bdys = NULL, *bdyf = NULL;
    PyArrayObject *bdzs = NULL, *bdzf = NULL, *f = NULL;

    if (!PyArg_ParseTuple(args, "ddiiOOddiiOOddiiOOdO:hw3crt",
                          &xs, &xf, &l, &lbdcnd, &bdxs_obj, &bdxf_obj,
                          &ys, &yf, &m, &mbdcnd, &bdys_obj, &bdyf_obj,
                          &zs, &zf, &n, &nbdcnd, &bdzs_obj, &bdzf_obj,
                          &elmbda, &f_obj)) {
        return NULL;
    }

    bdxs = as_fortran_double_copy(bdxs_obj, 2, "bdxs");
    bdxf = as_fortran_double_copy(bdxf_obj, 2, "bdxf");
    bdys = as_fortran_double_copy(bdys_obj, 2, "bdys");
    bdyf = as_fortran_double_copy(bdyf_obj, 2, "bdyf");
    bdzs = as_fortran_double_copy(bdzs_obj, 2, "bdzs");
    bdzf = as_fortran_double_copy(bdzf_obj, 2, "bdzf");
    f = as_fortran_double_copy(f_obj, 3, "f");
    if (bdxs == NULL || bdxf == NULL || bdys == NULL || bdyf == NULL ||
        bdzs == NULL || bdzf == NULL || f == NULL) {
        goto fail;
    }

    ldimf = (int)PyArray_DIM(f, 0);
    mdimf = (int)PyArray_DIM(f, 1);
    if (ldimf < l + 1 || mdimf < m + 1 || PyArray_DIM(f, 2) < n + 1) {
        PyErr_Format(PyExc_ValueError,
                     "f must have shape (ldimf, mdimf, n + 1) with ldimf >= l + 1 and mdimf >= m + 1; got (%zd, %zd, %zd)",
                     (Py_ssize_t)PyArray_DIM(f, 0), (Py_ssize_t)PyArray_DIM(f, 1),
                     (Py_ssize_t)PyArray_DIM(f, 2));
        goto fail;
    }
    if (lbdcnd == 3 || lbdcnd == 4) {
        if (!require_2d_shape(bdxs, m + 1, n + 1, "bdxs")) goto fail;
    }
    if (lbdcnd == 2 || lbdcnd == 3) {
        if (!require_2d_shape(bdxf, m + 1, n + 1, "bdxf")) goto fail;
    }
    if (mbdcnd == 3 || mbdcnd == 4) {
        if (!require_2d_shape(bdys, l + 1, n + 1, "bdys")) goto fail;
    }
    if (mbdcnd == 2 || mbdcnd == 3) {
        if (!require_2d_shape(bdyf, l + 1, n + 1, "bdyf")) goto fail;
    }
    if (nbdcnd == 3 || nbdcnd == 4) {
        if (!require_2d_shape(bdzs, l + 1, m + 1, "bdzs")) goto fail;
    }
    if (nbdcnd == 2 || nbdcnd == 3) {
        if (!require_2d_shape(bdzf, l + 1, m + 1, "bdzf")) goto fail;
    }

    pyfp_hw3crt(xs, xf, l, lbdcnd,
                PyArray_DATA(bdxs), (int)PyArray_DIM(bdxs, 0), (int)PyArray_DIM(bdxs, 1),
                PyArray_DATA(bdxf), (int)PyArray_DIM(bdxf, 0), (int)PyArray_DIM(bdxf, 1),
                ys, yf, m, mbdcnd,
                PyArray_DATA(bdys), (int)PyArray_DIM(bdys, 0), (int)PyArray_DIM(bdys, 1),
                PyArray_DATA(bdyf), (int)PyArray_DIM(bdyf, 0), (int)PyArray_DIM(bdyf, 1),
                zs, zf, n, nbdcnd,
                PyArray_DATA(bdzs), (int)PyArray_DIM(bdzs, 0), (int)PyArray_DIM(bdzs, 1),
                PyArray_DATA(bdzf), (int)PyArray_DIM(bdzf, 0), (int)PyArray_DIM(bdzf, 1),
                elmbda, ldimf, mdimf, PyArray_DATA(f), &pertrb, &ierror);

    Py_DECREF(bdxs);
    Py_DECREF(bdxf);
    Py_DECREF(bdys);
    Py_DECREF(bdyf);
    Py_DECREF(bdzs);
    Py_DECREF(bdzf);
    return Py_BuildValue("Ndi", (PyObject *)f, pertrb, ierror);

fail:
    Py_XDECREF(bdxs);
    Py_XDECREF(bdxf);
    Py_XDECREF(bdys);
    Py_XDECREF(bdyf);
    Py_XDECREF(bdzs);
    Py_XDECREF(bdzf);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_sor_standard1d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *f_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *f = NULL;
    const char *bcx_name;
    int bcx, mxloop, overflow = 0, loops = 0;
    double delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOdsddid:sor_standard1d",
                          &s_obj, &a_obj, &b_obj, &f_obj,
                          &delx, &bcx_name, &optarg, &undef, &mxloop, &tolerance)) {
        return NULL;
    }

    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 1, "s");
    a = as_fortran_double_copy(a_obj, 1, "a");
    b = as_fortran_double_copy(b_obj, 1, "b");
    f = as_fortran_double_copy(f_obj, 1, "f");
    if (s == NULL || a == NULL || b == NULL || f == NULL) {
        goto fail;
    }

    const npy_intp n = PyArray_DIM(s, 0);
    if (!require_len(a, n, "a") || !require_len(b, n, "b") || !require_len(f, n, "f")) {
        goto fail;
    }
    if (n < 3) {
        PyErr_SetString(PyExc_ValueError, "sor_standard1d requires at least 3 grid points");
        goto fail;
    }

    pyfp_sor_standard1d((int)n, PyArray_DATA(s), PyArray_DATA(a), PyArray_DATA(b),
                        PyArray_DATA(f), delx, bcx, optarg, undef, mxloop,
                        tolerance, &overflow, &relerr, &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(f);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_sor_standard2d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *c_obj, *f_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *c = NULL, *f = NULL;
    const char *bcy_name, *bcx_name;
    int bcy, bcx, mxloop, overflow = 0, loops = 0;
    double dely, delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOOddssddid:sor_standard2d",
                          &s_obj, &a_obj, &b_obj, &c_obj, &f_obj,
                          &dely, &delx, &bcy_name, &bcx_name, &optarg,
                          &undef, &mxloop, &tolerance)) {
        return NULL;
    }

    bcy = bc_code_from_string(bcy_name, "bcy");
    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcy < 0 || bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 2, "s");
    a = as_fortran_double_copy(a_obj, 2, "a");
    b = as_fortran_double_copy(b_obj, 2, "b");
    c = as_fortran_double_copy(c_obj, 2, "c");
    f = as_fortran_double_copy(f_obj, 2, "f");
    if (s == NULL || a == NULL || b == NULL || c == NULL || f == NULL) {
        goto fail;
    }

    const npy_intp ny = PyArray_DIM(s, 0);
    const npy_intp nx = PyArray_DIM(s, 1);
    if (!require_2d_shape(a, ny, nx, "a") || !require_2d_shape(b, ny, nx, "b") ||
        !require_2d_shape(c, ny, nx, "c") || !require_2d_shape(f, ny, nx, "f")) {
        goto fail;
    }
    if (ny < 3 || nx < 3) {
        PyErr_SetString(PyExc_ValueError, "sor_standard2d requires at least a 3 by 3 grid");
        goto fail;
    }

    pyfp_sor_standard2d((int)ny, (int)nx, PyArray_DATA(s), PyArray_DATA(a),
                        PyArray_DATA(b), PyArray_DATA(c), PyArray_DATA(f),
                        dely, delx, bcy, bcx, optarg, undef, mxloop, tolerance,
                        &overflow, &relerr, &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    Py_DECREF(f);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_sor_standard3d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *c_obj, *f_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *c = NULL, *f = NULL;
    const char *bcz_name, *bcy_name, *bcx_name;
    int bcz, bcy, bcx, mxloop, overflow = 0, loops = 0;
    double delz, dely, delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOOdddsssddid:sor_standard3d",
                          &s_obj, &a_obj, &b_obj, &c_obj, &f_obj,
                          &delz, &dely, &delx, &bcz_name, &bcy_name, &bcx_name,
                          &optarg, &undef, &mxloop, &tolerance)) {
        return NULL;
    }

    bcz = bc_code_from_string(bcz_name, "bcz");
    bcy = bc_code_from_string(bcy_name, "bcy");
    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcz < 0 || bcy < 0 || bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 3, "s");
    a = as_fortran_double_copy(a_obj, 3, "a");
    b = as_fortran_double_copy(b_obj, 3, "b");
    c = as_fortran_double_copy(c_obj, 3, "c");
    f = as_fortran_double_copy(f_obj, 3, "f");
    if (s == NULL || a == NULL || b == NULL || c == NULL || f == NULL) {
        goto fail;
    }

    const npy_intp nz = PyArray_DIM(s, 0);
    const npy_intp ny = PyArray_DIM(s, 1);
    const npy_intp nx = PyArray_DIM(s, 2);
    if (!require_3d_shape(a, nz, ny, nx, "a") || !require_3d_shape(b, nz, ny, nx, "b") ||
        !require_3d_shape(c, nz, ny, nx, "c") || !require_3d_shape(f, nz, ny, nx, "f")) {
        goto fail;
    }
    if (nz < 3 || ny < 3 || nx < 3) {
        PyErr_SetString(PyExc_ValueError, "sor_standard3d requires at least a 3 by 3 by 3 grid");
        goto fail;
    }

    pyfp_sor_standard3d((int)nz, (int)ny, (int)nx, PyArray_DATA(s),
                        PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c),
                        PyArray_DATA(f), delz, dely, delx, bcz, bcy, bcx,
                        optarg, undef, mxloop, tolerance, &overflow, &relerr,
                        &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    Py_DECREF(f);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(f);
    return NULL;
}

static PyObject *fishpack_sor_general2d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *c_obj, *d_obj, *e_obj, *fcoef_obj, *g_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *c = NULL;
    PyArrayObject *d = NULL, *e = NULL, *fcoef = NULL, *g = NULL;
    const char *bcy_name, *bcx_name;
    int bcy, bcx, mxloop, overflow = 0, loops = 0;
    double dely, delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOOOOOddssddid:sor_general2d",
                          &s_obj, &a_obj, &b_obj, &c_obj, &d_obj, &e_obj,
                          &fcoef_obj, &g_obj, &dely, &delx, &bcy_name, &bcx_name,
                          &optarg, &undef, &mxloop, &tolerance)) {
        return NULL;
    }

    bcy = bc_code_from_string(bcy_name, "bcy");
    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcy < 0 || bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 2, "s");
    a = as_fortran_double_copy(a_obj, 2, "a");
    b = as_fortran_double_copy(b_obj, 2, "b");
    c = as_fortran_double_copy(c_obj, 2, "c");
    d = as_fortran_double_copy(d_obj, 2, "d");
    e = as_fortran_double_copy(e_obj, 2, "e");
    fcoef = as_fortran_double_copy(fcoef_obj, 2, "fcoef");
    g = as_fortran_double_copy(g_obj, 2, "g");
    if (s == NULL || a == NULL || b == NULL || c == NULL ||
        d == NULL || e == NULL || fcoef == NULL || g == NULL) {
        goto fail;
    }

    const npy_intp ny = PyArray_DIM(s, 0);
    const npy_intp nx = PyArray_DIM(s, 1);
    if (!require_2d_shape(a, ny, nx, "a") || !require_2d_shape(b, ny, nx, "b") ||
        !require_2d_shape(c, ny, nx, "c") || !require_2d_shape(d, ny, nx, "d") ||
        !require_2d_shape(e, ny, nx, "e") || !require_2d_shape(fcoef, ny, nx, "fcoef") ||
        !require_2d_shape(g, ny, nx, "g")) {
        goto fail;
    }
    if (ny < 3 || nx < 3) {
        PyErr_SetString(PyExc_ValueError, "sor_general2d requires at least a 3 by 3 grid");
        goto fail;
    }

    pyfp_sor_general2d((int)ny, (int)nx, PyArray_DATA(s), PyArray_DATA(a),
                       PyArray_DATA(b), PyArray_DATA(c), PyArray_DATA(d),
                       PyArray_DATA(e), PyArray_DATA(fcoef), PyArray_DATA(g),
                       dely, delx, bcy, bcx, optarg, undef, mxloop, tolerance,
                       &overflow, &relerr, &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    Py_DECREF(d);
    Py_DECREF(e);
    Py_DECREF(fcoef);
    Py_DECREF(g);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(d);
    Py_XDECREF(e);
    Py_XDECREF(fcoef);
    Py_XDECREF(g);
    return NULL;
}

static PyObject *fishpack_sor_general3d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *c_obj, *d_obj, *e_obj, *fcoef_obj, *g_obj, *h_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *c = NULL;
    PyArrayObject *d = NULL, *e = NULL, *fcoef = NULL, *g = NULL, *h = NULL;
    const char *bcz_name, *bcy_name, *bcx_name;
    int bcz, bcy, bcx, mxloop, overflow = 0, loops = 0;
    double delz, dely, delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOdddsssddid:sor_general3d",
                          &s_obj, &a_obj, &b_obj, &c_obj, &d_obj, &e_obj,
                          &fcoef_obj, &g_obj, &h_obj, &delz, &dely, &delx,
                          &bcz_name, &bcy_name, &bcx_name, &optarg, &undef,
                          &mxloop, &tolerance)) {
        return NULL;
    }

    bcz = bc_code_from_string(bcz_name, "bcz");
    bcy = bc_code_from_string(bcy_name, "bcy");
    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcz < 0 || bcy < 0 || bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 3, "s");
    a = as_fortran_double_copy(a_obj, 3, "a");
    b = as_fortran_double_copy(b_obj, 3, "b");
    c = as_fortran_double_copy(c_obj, 3, "c");
    d = as_fortran_double_copy(d_obj, 3, "d");
    e = as_fortran_double_copy(e_obj, 3, "e");
    fcoef = as_fortran_double_copy(fcoef_obj, 3, "fcoef");
    g = as_fortran_double_copy(g_obj, 3, "g");
    h = as_fortran_double_copy(h_obj, 3, "h");
    if (s == NULL || a == NULL || b == NULL || c == NULL || d == NULL ||
        e == NULL || fcoef == NULL || g == NULL || h == NULL) {
        goto fail;
    }

    const npy_intp nz = PyArray_DIM(s, 0);
    const npy_intp ny = PyArray_DIM(s, 1);
    const npy_intp nx = PyArray_DIM(s, 2);
    if (!require_3d_shape(a, nz, ny, nx, "a") || !require_3d_shape(b, nz, ny, nx, "b") ||
        !require_3d_shape(c, nz, ny, nx, "c") || !require_3d_shape(d, nz, ny, nx, "d") ||
        !require_3d_shape(e, nz, ny, nx, "e") ||
        !require_3d_shape(fcoef, nz, ny, nx, "fcoef") ||
        !require_3d_shape(g, nz, ny, nx, "g") || !require_3d_shape(h, nz, ny, nx, "h")) {
        goto fail;
    }
    if (nz < 3 || ny < 3 || nx < 3) {
        PyErr_SetString(PyExc_ValueError, "sor_general3d requires at least a 3 by 3 by 3 grid");
        goto fail;
    }

    pyfp_sor_general3d((int)nz, (int)ny, (int)nx, PyArray_DATA(s),
                       PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c),
                       PyArray_DATA(d), PyArray_DATA(e), PyArray_DATA(fcoef),
                       PyArray_DATA(g), PyArray_DATA(h), delz, dely, delx,
                       bcz, bcy, bcx, optarg, undef, mxloop, tolerance,
                       &overflow, &relerr, &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    Py_DECREF(d);
    Py_DECREF(e);
    Py_DECREF(fcoef);
    Py_DECREF(g);
    Py_DECREF(h);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(d);
    Py_XDECREF(e);
    Py_XDECREF(fcoef);
    Py_XDECREF(g);
    Py_XDECREF(h);
    return NULL;
}

static PyObject *fishpack_sor_biharmonic2d(PyObject *self, PyObject *args)
{
    (void)self;
    PyObject *s_obj, *a_obj, *b_obj, *c_obj, *d_obj, *e_obj, *fcoef_obj;
    PyObject *g_obj, *h_obj, *icoef_obj, *jcoef_obj;
    PyArrayObject *s = NULL, *a = NULL, *b = NULL, *c = NULL, *d = NULL;
    PyArrayObject *e = NULL, *fcoef = NULL, *g = NULL, *h = NULL;
    PyArrayObject *icoef = NULL, *jcoef = NULL;
    const char *bcy_name, *bcx_name;
    int bcy, bcx, mxloop, overflow = 0, loops = 0;
    double dely, delx, optarg, undef, tolerance, relerr = 1.0;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOddssddid:sor_biharmonic2d",
                          &s_obj, &a_obj, &b_obj, &c_obj, &d_obj, &e_obj,
                          &fcoef_obj, &g_obj, &h_obj, &icoef_obj, &jcoef_obj,
                          &dely, &delx, &bcy_name, &bcx_name, &optarg,
                          &undef, &mxloop, &tolerance)) {
        return NULL;
    }

    bcy = bc_code_from_string(bcy_name, "bcy");
    bcx = bc_code_from_string(bcx_name, "bcx");
    if (bcy < 0 || bcx < 0) {
        return NULL;
    }

    s = as_fortran_double_copy(s_obj, 2, "s");
    a = as_fortran_double_copy(a_obj, 2, "a");
    b = as_fortran_double_copy(b_obj, 2, "b");
    c = as_fortran_double_copy(c_obj, 2, "c");
    d = as_fortran_double_copy(d_obj, 2, "d");
    e = as_fortran_double_copy(e_obj, 2, "e");
    fcoef = as_fortran_double_copy(fcoef_obj, 2, "fcoef");
    g = as_fortran_double_copy(g_obj, 2, "g");
    h = as_fortran_double_copy(h_obj, 2, "h");
    icoef = as_fortran_double_copy(icoef_obj, 2, "icoef");
    jcoef = as_fortran_double_copy(jcoef_obj, 2, "jcoef");
    if (s == NULL || a == NULL || b == NULL || c == NULL || d == NULL ||
        e == NULL || fcoef == NULL || g == NULL || h == NULL ||
        icoef == NULL || jcoef == NULL) {
        goto fail;
    }

    const npy_intp ny = PyArray_DIM(s, 0);
    const npy_intp nx = PyArray_DIM(s, 1);
    if (!require_2d_shape(a, ny, nx, "a") || !require_2d_shape(b, ny, nx, "b") ||
        !require_2d_shape(c, ny, nx, "c") || !require_2d_shape(d, ny, nx, "d") ||
        !require_2d_shape(e, ny, nx, "e") || !require_2d_shape(fcoef, ny, nx, "fcoef") ||
        !require_2d_shape(g, ny, nx, "g") || !require_2d_shape(h, ny, nx, "h") ||
        !require_2d_shape(icoef, ny, nx, "icoef") || !require_2d_shape(jcoef, ny, nx, "jcoef")) {
        goto fail;
    }
    if (ny < 5 || nx < 5) {
        PyErr_SetString(PyExc_ValueError, "sor_biharmonic2d requires at least a 5 by 5 grid");
        goto fail;
    }

    pyfp_sor_biharmonic2d((int)ny, (int)nx, PyArray_DATA(s),
                          PyArray_DATA(a), PyArray_DATA(b), PyArray_DATA(c),
                          PyArray_DATA(d), PyArray_DATA(e), PyArray_DATA(fcoef),
                          PyArray_DATA(g), PyArray_DATA(h), PyArray_DATA(icoef),
                          PyArray_DATA(jcoef), dely, delx, bcy, bcx, optarg,
                          undef, mxloop, tolerance, &overflow, &relerr, &loops);

    Py_DECREF(a);
    Py_DECREF(b);
    Py_DECREF(c);
    Py_DECREF(d);
    Py_DECREF(e);
    Py_DECREF(fcoef);
    Py_DECREF(g);
    Py_DECREF(h);
    Py_DECREF(icoef);
    Py_DECREF(jcoef);
    return Py_BuildValue("Ndii", (PyObject *)s, relerr, overflow, loops);

fail:
    Py_XDECREF(s);
    Py_XDECREF(a);
    Py_XDECREF(b);
    Py_XDECREF(c);
    Py_XDECREF(d);
    Py_XDECREF(e);
    Py_XDECREF(fcoef);
    Py_XDECREF(g);
    Py_XDECREF(h);
    Py_XDECREF(icoef);
    Py_XDECREF(jcoef);
    return NULL;
}

static PyObject *fishpack_fftpack_real(PyObject *args, int kind, int direction, const char *name)
{
    PyObject *x_obj;
    PyArrayObject *x = NULL;

    if (!PyArg_ParseTuple(args, "O", &x_obj)) {
        return NULL;
    }
    x = as_fortran_double_copy(x_obj, 1, name);
    if (x == NULL) {
        return NULL;
    }
    const npy_intp n = PyArray_DIM(x, 0);
    if (n < 1) {
        PyErr_Format(PyExc_ValueError, "%s requires at least one sample", name);
        Py_DECREF(x);
        return NULL;
    }
    pyfp_fftpack_real_transform((int)n, PyArray_DATA(x), kind, direction);
    return (PyObject *)x;
}

static PyObject *fishpack_fftpack_complex(PyObject *args, int direction, const char *name)
{
    PyObject *x_obj;
    PyArrayObject *x = NULL;

    if (!PyArg_ParseTuple(args, "O", &x_obj)) {
        return NULL;
    }
    x = as_fortran_double_copy(x_obj, 1, name);
    if (x == NULL) {
        return NULL;
    }
    const npy_intp n2 = PyArray_DIM(x, 0);
    if (n2 < 2 || n2 % 2 != 0) {
        PyErr_Format(PyExc_ValueError, "%s expects an interleaved real/imag float64 array with even length", name);
        Py_DECREF(x);
        return NULL;
    }
    pyfp_fftpack_complex_transform((int)(n2 / 2), PyArray_DATA(x), direction);
    return (PyObject *)x;
}

static PyObject *fishpack_rfftf(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 0, 0, "rfftf");
}

static PyObject *fishpack_rfftb(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 0, 1, "rfftb");
}

static PyObject *fishpack_sint(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 1, 0, "sint");
}

static PyObject *fishpack_cost(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 2, 0, "cost");
}

static PyObject *fishpack_sinqf(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 3, 0, "sinqf");
}

static PyObject *fishpack_sinqb(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 3, 1, "sinqb");
}

static PyObject *fishpack_cosqf(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 4, 0, "cosqf");
}

static PyObject *fishpack_cosqb(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_real(args, 4, 1, "cosqb");
}

static PyObject *fishpack_cfftf(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_complex(args, 0, "cfftf");
}

static PyObject *fishpack_cfftb(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_fftpack_complex(args, 1, "cfftb");
}

static PyObject *fishpack_hwsplr(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hwsplr, 1);
}

static PyObject *fishpack_hwscyl(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hwscyl, 1);
}

static PyObject *fishpack_hstplr(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hstplr, 0);
}

static PyObject *fishpack_hstcyl(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hstcyl, 0);
}

static PyObject *fishpack_hwsssp(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hwsssp, 1);
}

static PyObject *fishpack_hstssp(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hstssp, 0);
}

static PyObject *fishpack_hwscsp(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hwscsp, 1);
}

static PyObject *fishpack_hstcsp(PyObject *self, PyObject *args)
{
    (void)self;
    return fishpack_helmholtz_2d(args, pyfp_hstcsp, 0);
}

static PyMethodDef fishpack_methods[] = {
    {"backend_info", (PyCFunction)fishpack_backend_info, METH_NOARGS,
     "Return details about the compiled PyFishPack backend."},
    {"genbun", (PyCFunction)fishpack_genbun, METH_VARARGS,
     "Solve a centered-grid separable elliptic system using the modern Fortran GENBUN backend."},
    {"poistg", (PyCFunction)fishpack_poistg, METH_VARARGS,
     "Solve a staggered-grid separable elliptic system using the modern Fortran POISTG backend."},
    {"pois3d", (PyCFunction)fishpack_pois3d, METH_VARARGS,
     "Solve a 3D separable elliptic system using the modern Fortran POIS3D backend."},
    {"hwscrt", (PyCFunction)fishpack_hwscrt, METH_VARARGS,
     "Solve a centered-grid Cartesian Helmholtz equation using the modern Fortran HWSCRT backend."},
    {"hstcrt", (PyCFunction)fishpack_hstcrt, METH_VARARGS,
     "Solve a staggered-grid Cartesian Helmholtz equation using the modern Fortran HSTCRT backend."},
    {"hw3crt", (PyCFunction)fishpack_hw3crt, METH_VARARGS,
     "Solve a 3D centered-grid Cartesian Helmholtz equation using the modern Fortran HW3CRT backend."},
    {"hwsplr", (PyCFunction)fishpack_hwsplr, METH_VARARGS,
     "Solve a centered-grid polar Helmholtz equation using the modern Fortran HWSPLR backend."},
    {"hwscyl", (PyCFunction)fishpack_hwscyl, METH_VARARGS,
     "Solve a centered-grid cylindrical Helmholtz equation using the modern Fortran HWSCYL backend."},
    {"hstplr", (PyCFunction)fishpack_hstplr, METH_VARARGS,
     "Solve a staggered-grid polar Helmholtz equation using the modern Fortran HSTPLR backend."},
    {"hstcyl", (PyCFunction)fishpack_hstcyl, METH_VARARGS,
     "Solve a staggered-grid cylindrical Helmholtz equation using the modern Fortran HSTCYL backend."},
    {"hwsssp", (PyCFunction)fishpack_hwsssp, METH_VARARGS,
     "Solve a centered-grid spherical Helmholtz equation using the modern Fortran HWSSSP backend."},
    {"hstssp", (PyCFunction)fishpack_hstssp, METH_VARARGS,
     "Solve a staggered-grid spherical Helmholtz equation using the modern Fortran HSTSSP backend."},
    {"hwscsp", (PyCFunction)fishpack_hwscsp, METH_VARARGS,
     "Solve a centered-grid axisymmetric spherical Helmholtz equation using the modern Fortran HWSCSP backend."},
    {"hstcsp", (PyCFunction)fishpack_hstcsp, METH_VARARGS,
     "Solve a staggered-grid axisymmetric spherical Helmholtz equation using the modern Fortran HSTCSP backend."},
    {"sor_standard1d", (PyCFunction)fishpack_sor_standard1d, METH_VARARGS,
     "Solve a one-dimensional variable-coefficient elliptic equation using the modern Fortran SOR backend."},
    {"sor_standard2d", (PyCFunction)fishpack_sor_standard2d, METH_VARARGS,
     "Solve a two-dimensional variable-coefficient elliptic equation using the modern Fortran SOR backend."},
    {"sor_standard3d", (PyCFunction)fishpack_sor_standard3d, METH_VARARGS,
     "Solve a three-dimensional variable-coefficient elliptic equation using the modern Fortran SOR backend."},
    {"sor_general2d", (PyCFunction)fishpack_sor_general2d, METH_VARARGS,
     "Solve a two-dimensional general-form elliptic equation using the modern Fortran SOR backend."},
    {"sor_general3d", (PyCFunction)fishpack_sor_general3d, METH_VARARGS,
     "Solve a three-dimensional general-form elliptic equation using the modern Fortran SOR backend."},
    {"sor_biharmonic2d", (PyCFunction)fishpack_sor_biharmonic2d, METH_VARARGS,
     "Solve a two-dimensional general-form biharmonic equation using the modern Fortran SOR backend."},
    {"rfftf", (PyCFunction)fishpack_rfftf, METH_VARARGS,
     "Apply the FFTPACK real periodic forward transform to a one-dimensional float64 array."},
    {"rfftb", (PyCFunction)fishpack_rfftb, METH_VARARGS,
     "Apply the FFTPACK real periodic backward transform to a one-dimensional float64 array."},
    {"sint", (PyCFunction)fishpack_sint, METH_VARARGS,
     "Apply the FFTPACK sine transform to a one-dimensional float64 array."},
    {"cost", (PyCFunction)fishpack_cost, METH_VARARGS,
     "Apply the FFTPACK cosine transform to a one-dimensional float64 array."},
    {"sinqf", (PyCFunction)fishpack_sinqf, METH_VARARGS,
     "Apply the FFTPACK quarter-wave sine forward transform to a one-dimensional float64 array."},
    {"sinqb", (PyCFunction)fishpack_sinqb, METH_VARARGS,
     "Apply the FFTPACK quarter-wave sine backward transform to a one-dimensional float64 array."},
    {"cosqf", (PyCFunction)fishpack_cosqf, METH_VARARGS,
     "Apply the FFTPACK quarter-wave cosine forward transform to a one-dimensional float64 array."},
    {"cosqb", (PyCFunction)fishpack_cosqb, METH_VARARGS,
     "Apply the FFTPACK quarter-wave cosine backward transform to a one-dimensional float64 array."},
    {"cfftf", (PyCFunction)fishpack_cfftf, METH_VARARGS,
     "Apply the FFTPACK complex forward transform to an interleaved real/imag float64 array."},
    {"cfftb", (PyCFunction)fishpack_cfftb, METH_VARARGS,
     "Apply the FFTPACK complex backward transform to an interleaved real/imag float64 array."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fishpack_module = {
    PyModuleDef_HEAD_INIT,
    "fishpack",
    "C extension wrapper for the modern Fortran PyFishPack backend.",
    -1,
    fishpack_methods
};

PyMODINIT_FUNC PyInit_fishpack(void)
{
    import_array();
    return PyModule_Create(&fishpack_module);
}
