#if !defined (COREIMPL)
#define COREIMPL

#include "core.h"

#include "physimpl.h"
#include "meshimpl.h"
#include "bcimpl.h"
#include "commimpl.h"
#include "solnimpl.h"

struct __Core{

    Mesh     *mesh;
    Phys     *phys;
    BC       *bc;
    Comm     *comm;
    Soln     *soln;
    int      nmesh;
    int      nphys;
    int      nbc;
    int      ncomm;
    int      nsoln;
};

#endif
