#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"
#include "igabeamtest.h"

void igabeamtest(double Lx, double Ly, double Lz, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z, int **bc_periodic_, int **nelem_global_, double ***knotVector_global_, double **wVector_global_, double **XVector_global_)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int ndim=3;

    // ---- allocate memory first
    *bc_periodic_       = (int*)malloc(ndim*sizeof(int));
    *nelem_global_      = (int*)malloc(ndim*sizeof(int));
    *knotVector_global_ = (double **)malloc(ndim*sizeof(double*));

    int    *bc_periodic        = *bc_periodic_;
    int    *nelem_global       = *nelem_global_;
    double **knotVector_global = *knotVector_global_;

    // ---- parametric (xi,eta,zeta)
    bc_periodic[0]=bc_periodic_x;
    bc_periodic[1]=bc_periodic_y;
    bc_periodic[2]=bc_periodic_z;

    for (int idim=0; idim<ndim; idim++) { nelem_global[idim]=int_pow(2,mref); }
    if (rank==0) { for (int idim=0; idim<ndim; idim++) { knotVector_global[idim]=(double*) malloc((nelem_global[idim]+1+2*porder)*sizeof(double)); } }
    else         { for (int idim=0; idim<ndim; idim++) { knotVector_global[idim]=(double*) malloc(0*sizeof(double)); } }

    if (rank==0)
    {
        double l[ndim]; for (int idim=0; idim<ndim; idim++) { l[idim]=1.0/nelem_global[idim]; }
        for (int idim=0; idim<ndim; idim++)
        {
            int iknot=0;
            if (bc_periodic[idim]==0)
            {
                while (iknot < porder+1)                      { knotVector_global[idim][iknot++] = 0.; }
                while (iknot < nelem_global[idim]+porder)     { knotVector_global[idim][iknot]   = (iknot-porder)*l[idim]; iknot++; }
                while (iknot < nelem_global[idim]+1+2*porder) { knotVector_global[idim][iknot++] = 1.; }
            }
            else
            {
                while (iknot < nelem_global[idim]+1+2*porder) { knotVector_global[idim][iknot] = (iknot-porder)*l[idim]; iknot++; }
            }
        }
    }

    // ---- reference (X)
    int nbasis_global[] = {nelem_global[0]+porder,nelem_global[1]+porder,nelem_global[2]+porder};
    *wVector_global_    = (double*)malloc(     nbasis_global[0]*nbasis_global[1]*nbasis_global[2]*sizeof(double));
    *XVector_global_    = (double*)malloc(ndim*nbasis_global[0]*nbasis_global[1]*nbasis_global[2]*sizeof(double));

    double *wVector_global = *wVector_global_;
    double *XVector_global = *XVector_global_;

    assert(porder==2);
    double L[]={Lx,Ly,Lz};
    int ia=0;
    for (int i=0; i<nbasis_global[0]; i++) {
    for (int j=0; j<nbasis_global[1]; j++) {    
    for (int k=0; k<nbasis_global[2]; k++) {
        wVector_global[ia]=1.0;
        XVector_global[ndim*ia+0]=( bc_periodic[0] ? ((double)i-0.5)/((double)(nbasis_global[0]-2)) : (double)i/((double)(nbasis_global[0]-1)) )*L[0];
        XVector_global[ndim*ia+1]=( bc_periodic[1] ? ((double)j-0.5)/((double)(nbasis_global[1]-2)) : (double)j/((double)(nbasis_global[1]-1)) )*L[1];
        XVector_global[ndim*ia+2]=( bc_periodic[2] ? ((double)k-0.5)/((double)(nbasis_global[2]-2)) : (double)k/((double)(nbasis_global[2]-1)) )*L[2];
        XVector_global[ndim*ia+0]+=0.5*Lx/nbasis_global[0]*(double)(i>2&&j>2&&k>2&&i<nbasis_global[0]-2&&j<nbasis_global[1]-2&&k<nbasis_global[2]-2)*sin(3*i*k*j);
        XVector_global[ndim*ia+1]+=0.5*Ly/nbasis_global[1]*(double)(i>2&&j>2&&k>2&&i<nbasis_global[0]-2&&j<nbasis_global[1]-2&&k<nbasis_global[2]-2)*sin(5*i*k*j);
        XVector_global[ndim*ia+2]+=0.5*Lz/nbasis_global[2]*(double)(i>2&&j>2&&k>2&&i<nbasis_global[0]-2&&j<nbasis_global[1]-2&&k<nbasis_global[2]-2)*sin(7*i*k*j);
        ia++;
    }}}
}
