#if !defined (SOLNIMPL)
#define SOLNIMPL

#include "soln.h"

struct __Soln{

    double *U_hist;
    int    nDof;
    int    *nDof_array;
    int    *iDof_displ_array;
    int    nDof_global;
    int    *universalIdx;
};

#endif
