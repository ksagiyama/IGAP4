#if !defined (ASSMIMPL)
#define ASSMIMPL

#include "petscsnes.h"
#include "assm.h"

typedef PetscErrorCode (*Residual_ptr)(SNES,Vec,Vec,void *);
typedef PetscErrorCode (*Tangent_ptr)(SNES,Vec,Mat,Mat,void *);

struct __Assm{
    Residual_ptr Residual;
    Tangent_ptr  Tangent;
    Vec          RR;
    Mat          TT;
};

#endif
