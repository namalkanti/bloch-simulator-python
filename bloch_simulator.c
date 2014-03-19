#include <Python.h>
#include <numpy/arrayobject.h>

#include "bloch.c"

//Docstrings
static char module_docstring[] = "Hargreaves Bloch Equation simulator implemented as a C extension."

static char bloch_docstring[] = "Bloch equation simulator."

static PyObject* bloch(PyObject *self, PyObject *args){

    double t1, t2;
    int mode;
    PyObject *b1, *gr, *tp, *df, *dp, mx, my, mz;

    if (!PyArg_ParseTuple(args, "000dd00i000", &b1, &gr, &tp, &t1, &t2, 
                &df, &dp, &mode, &mx, &my, &mz)){
        return NULL;
    }
    
    ntime =;

    b1_parameters b1 = process_b1_field();

    gradient_parameters gradients = process_gradients();

    time_parameters time_points = process_time_points();

    t1 = ;
    t2 = ;

    df = ;
    nf = ;

    position_parameters position_points = process_position_points();

    mode_parameters mode = process_mode_information();

    magnetization_parameters magnetization = process_magnetization();

    blochsimfz(b1.blr, b1.bli, gradients.gx, gradients.gy, gradients.gz, time_points.tp,
            ntime, t1, t2, df, nf, position_points.dx, position_points.dy, position_points.dz,
            magnetization.mx, magnetization.my, magnetization.mz, mode.md);
}

static PyMethodDef module_methods[] = {
    {"bloch", bloch, METH_VARARGS, bloch_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_bloch(void){
    PyObject* m = PyInitModule("bloch_simulator", modules_methods, module_docstring);
    if (NULL == m){
        return;
    }
    import_array();
}
