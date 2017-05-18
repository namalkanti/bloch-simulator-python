#define NPY_1_7_API_VERSION 0x00000007
#include <Python.h>
#include <numpy/arrayobject.h>

#include "bloch.c"

//Docstrings
static char module_docstring[] = "Hargreaves Bloch Equation simulator implemented as a C extension for python.";

static char bloch_docstring[] = "Bloch equation simulator.";

static PyObject* bloch(PyObject *self, PyObject *args){

    //Arguement declarations
    //double t1, t2;
    int nf, mode, n_pos;
    PyObject *py_b1_real, *py_b1_imag, *py_grx, *py_gry, *py_grz, *py_t1, *py_t2, *py_tp, *py_df, *py_dx, *py_dy, *py_dz;
    PyObject *py_mx, *py_my, *py_mz;

    //Bloch sim arugments declarations
    PyObject *b1_real_arr, *b1_imag_arr, *grx_arr, *gry_arr, *grz_arr, *tp_arr, *t1_arr, *t2_arr, *df_arr, *dx_arr, *dy_arr, *dz_arr, *mx_arr, *my_arr, *mz_arr;
    int ntime;
    double *b1_real, *b1_imag, *grx, *gry, *grz, *tp, *t1, *t2, *df, *dx, *dy, *dz, *mx, *my, *mz;

    if (!PyArg_ParseTuple(args, "OOOOOOiOOOiOOOiiOOO", &py_b1_real, &py_b1_imag, &py_grx, &py_gry, &py_grz, &py_tp, &ntime, &py_t1, &py_t2,
                &py_df, &nf, &py_dx, &py_dy, &py_dz, &n_pos, &mode, &py_mx, &py_my, &py_mz)){
        return NULL;
    }

    b1_real_arr = PyArray_FROM_OTF(py_b1_real, NPY_DOUBLE, NPY_IN_ARRAY);
    b1_imag_arr = PyArray_FROM_OTF(py_b1_imag, NPY_DOUBLE, NPY_IN_ARRAY);
    grx_arr = PyArray_FROM_OTF(py_grx, NPY_DOUBLE, NPY_IN_ARRAY);
    gry_arr = PyArray_FROM_OTF(py_gry, NPY_DOUBLE, NPY_IN_ARRAY);
    grz_arr = PyArray_FROM_OTF(py_grz, NPY_DOUBLE, NPY_IN_ARRAY);
    tp_arr = PyArray_FROM_OTF(py_tp, NPY_DOUBLE, NPY_IN_ARRAY);
    t1_arr = PyArray_FROM_OTF(py_t1, NPY_DOUBLE, NPY_IN_ARRAY);
    t2_arr = PyArray_FROM_OTF(py_t2, NPY_DOUBLE, NPY_IN_ARRAY);
    df_arr = PyArray_FROM_OTF(py_df, NPY_DOUBLE, NPY_IN_ARRAY);
    dx_arr = PyArray_FROM_OTF(py_dx, NPY_DOUBLE, NPY_IN_ARRAY);
    dy_arr = PyArray_FROM_OTF(py_dy, NPY_DOUBLE, NPY_IN_ARRAY);
    dz_arr = PyArray_FROM_OTF(py_dz, NPY_DOUBLE, NPY_IN_ARRAY);
    mx_arr = PyArray_FROM_OTF(py_mx, NPY_DOUBLE, NPY_INOUT_ARRAY);
    my_arr = PyArray_FROM_OTF(py_my, NPY_DOUBLE, NPY_INOUT_ARRAY);
    mz_arr = PyArray_FROM_OTF(py_mz, NPY_DOUBLE, NPY_INOUT_ARRAY);

    b1_real = (double *) PyArray_DATA(b1_real_arr);
    b1_imag = (double *) PyArray_DATA(b1_imag_arr);
    grx = (double *) PyArray_DATA(grx_arr);
    gry = (double *) PyArray_DATA(gry_arr);
    grz = (double *) PyArray_DATA(grz_arr);
    tp = (double *) PyArray_DATA(tp_arr);
    t1 = (double *) PyArray_DATA(t1_arr);
    t2 = (double *) PyArray_DATA(t2_arr);
    df = (double *) PyArray_DATA(df_arr);
    dx = (double *) PyArray_DATA(dx_arr);
    dy = (double *) PyArray_DATA(dy_arr);
    dz = (double *) PyArray_DATA(dz_arr);
    mx = (double *) PyArray_DATA(mx_arr);
    my = (double *) PyArray_DATA(my_arr);
    mz = (double *) PyArray_DATA(mz_arr);

    blochsimfz(b1_real, b1_imag, grx, gry, grz, tp, ntime, t1, t2, df, nf, dx, dy, dz, n_pos, mx, my, mz, mode);

    Py_DECREF(b1_real_arr);
    Py_DECREF(b1_imag_arr);
    Py_DECREF(grx_arr);
    Py_DECREF(gry_arr);
    Py_DECREF(grz_arr);
    Py_DECREF(tp_arr);
    Py_DECREF(t1_arr);
    Py_DECREF(t2_arr);
    Py_DECREF(df_arr);
    Py_DECREF(dx_arr);
    Py_DECREF(dy_arr);
    Py_DECREF(dz_arr);
    Py_DECREF(mx_arr);
    Py_DECREF(my_arr);
    Py_DECREF(mz_arr);

    return Py_BuildValue("");
}

static PyMethodDef module_methods[] = {
    {"bloch_c", bloch, METH_VARARGS, bloch_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef bloch_module = {
    PyModuleDef_HEAD_INIT,
    "bloch_simulator",
    module_docstring,
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_bloch_simulator(void){
    import_array();
    return PyModule_Create(&bloch_module);
}
