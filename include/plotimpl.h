#if !defined (PLOTIMPL)
#define PLOTIMPL

#include <mgl2/mgl_cf.h>
#include <mgl2/mpi.h>
#include "plot.h"

struct __Plot{

    HMGL   *gr;
    int    *figsize;
    double *relplot;
    double *rotate;
    double *ranges;
    double *crange;

    int    nppe;
    int    *uid;
    double uscale;

    const char *fname;
};

#endif
