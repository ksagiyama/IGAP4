#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "meshimpl.h"
#include "igaoshape.h"
#include "meshigaoshape.h"

Mesh MeshIGAOshape(double R, double r, double th, int mref, int porder)
{
    Mesh mesh=(Mesh)malloc(1*sizeof(struct __Mesh));
    int *bc_periodic;
    int *nelem_global;
    double **knotVector_global;
    double *wVector_global;
    double *XVector_global;
    igaoshape(R,r,th,mref,porder,&bc_periodic,&nelem_global,&knotVector_global,&wVector_global,&XVector_global);
    mesh->porder            = porder;
    mesh->bc_periodic       = bc_periodic;
    mesh->nelem_global      = nelem_global;
    mesh->knotVector_global = knotVector_global;
    mesh->wVector_global    = wVector_global;
    mesh->XVector_global    = XVector_global;

    return mesh;
}
