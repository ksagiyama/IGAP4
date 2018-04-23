#if !defined (COREIMPL)
#define COREIMPL

#include "core.h"

#include "physimpl.h"
#include "meshimpl.h"
#include "bcimpl.h"
#include "commimpl.h"

struct __Core{

    Phys     *phys;
    Mesh     *mesh;
    BC       *bc;
    Comm     *comm;
    double   *U_hist;

};

#endif
