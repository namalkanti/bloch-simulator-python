//Structs for clean input initialization

typedef struct {
    double* blr;
    double* bli;
}b1_parameters;

typedef struct {
    double* gx;
    double* gy;
    double* gz;
    int gyaflag;
    int gzaflag;
}gradient_parameters;

typedef struct {
    double* tp;
    double* ti; 
}time_parameters;

typedef struct {
    double* npos;
    double* nposM;
    double* nposN;
    double* nfnpos;
    double* dx;
    double* dy;
    double* dz;
    int dyaflag;
    int dzaflag;
}position_parameters;

typedef struct {
    int md;
    int ntout;
    double ntnfnpos;
}mode_parameters;

typedef struct {
    double* mx;
    double* my;
    double* mz;
    double* mxout;
    double* myout;
    double* mzout;
    double* mxin;
    double* myin;
    double* mzin;
}magnetization_parameters;
