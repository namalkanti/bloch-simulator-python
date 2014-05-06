#include <Python.h>
#include <numpy/arrayobject.h>

#include "bloch.c"

//Declarations for argument processing helper functions.
b1_parameters process_b1_field(Py_complex b1, int ntime);
gradient_parameters process_gradients();
time_paremters process_time_points();
position_parameters position_points();
mode_parameters process_mode_information();
magnetization_parameters process_magnetization();

//Docstrings
static char module_docstring[] = "Hargreaves Bloch Equation simulator implemented as a C extension."

static char bloch_docstring[] = "Bloch equation simulator."

static PyObject* bloch(PyObject *self, PyObject *args){

    double t1, t2;
    int mode;
    PyObject *gr, *tp, *df, *dp, mx, my, mz;
    Py_complex b1_input;

    if (!PyArg_ParseTuple(args, "D00dd00i000", &b1_input, &gr, &tp, &t1, &t2, 
                &df, &dp, &mode, &mx, &my, &mz)){
        return NULL;
    }

    int ntime;
    
    ntime =;

    b1_parameters b1 = process_b1_field(b1_input);

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

//Function takes in and returns b1 parameters ready for bloch sim.
//If complex, split up; if real, allocate an imaginary part.
b1_parameters process_b1_field(Py_complex b1){
    b1_parameters parameters = malloc(sizeof(b1_parameters));
    parameters.blr = b1.real;
    parameters.bli = b1.imag; 
    return parameters;
}

