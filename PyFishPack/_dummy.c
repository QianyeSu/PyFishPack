/* Dummy C extension for Windows compatibility */
#include <Python.h>

static PyObject* dummy_function(PyObject* self, PyObject* args) {
    return Py_None;
}

static PyMethodDef DummyMethods[] = {
    {"dummy", dummy_function, METH_VARARGS, "Dummy function"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef dummymodule = {
    PyModuleDef_HEAD_INIT,
    "_dummy",
    NULL,
    -1,
    DummyMethods
};

PyMODINIT_FUNC PyInit__dummy(void) {
    return PyModule_Create(&dummymodule);
}