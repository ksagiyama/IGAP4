#if !defined (SOLVIMPL)
#define SOLVIMPL

#include "petscsnes.h"
#include "solv.h"

struct __Solv{
    int               converged;    
    SNES              snes;
    Vec               U;
};

#endif
