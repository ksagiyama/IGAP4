#if !defined (IOCOMMIMPL)
#define IOCOMMIMPL

#include "iocomm.h"

struct __IOcomm{

    int nDof;
    int *nDof_array;
    int *iDof_displ_array;
    int nDof_global;
    int *universalIdx;

};

#endif
