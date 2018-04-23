#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "meshimpl.h"
#include "igabeam.h"
#include "meshigabeam.h"

Mesh MeshIGABeam(double Lx, double Ly, double Lz, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z)
{
    Mesh mesh=(Mesh)malloc(1*sizeof(struct __Mesh));
    int *bc_periodic;
    int *nelem_global;
    double **knotVector_global;
    double *wVector_global;
    double *XVector_global;
    igabeam(Lx,Ly,Lz,mref,porder,bc_periodic_x,bc_periodic_y,bc_periodic_z,&bc_periodic,&nelem_global,&knotVector_global,&wVector_global,&XVector_global);
    mesh->porder            = porder;
    mesh->bc_periodic       = bc_periodic;
    mesh->nelem_global      = nelem_global;
    mesh->knotVector_global = knotVector_global;
    mesh->wVector_global    = wVector_global;
    mesh->XVector_global    = XVector_global;

    return mesh;
}
