#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"
#include "igaannulus.h"

void igaannulus(double R, double r, double th, int mref, int porder, int **bc_periodic_, int **nelem_global_, double ***knotVector_global_, double **wVector_global_, double **XVector_global_)
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
    bc_periodic[0]=1;
    bc_periodic[1]=0;
    bc_periodic[2]=0;

    nelem_global[0]=8;
    nelem_global[1]=1;
    nelem_global[2]=1;

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
if (rank==0){
    knotVector_global[0][0] =0.;
    knotVector_global[0][1] =1.;
    knotVector_global[0][2] =1.+0*1.e-14;
    knotVector_global[0][3] =2.;
    knotVector_global[0][4] =2.+0*1.e-14;
    knotVector_global[0][5] =3;
    knotVector_global[0][6] =3.+0*1.e-14;
    knotVector_global[0][7] =4;
    knotVector_global[0][8] =4.+0*1.e-14;
    knotVector_global[0][9] =5.;
    knotVector_global[0][10]=5.+0*1.e-14;
    knotVector_global[0][11]=6.;
    knotVector_global[0][12]=6.+0*1.e-14;
}
    // ---- reference (X)
    int nbasis_global[] = {nelem_global[0]+porder,nelem_global[1]+porder,nelem_global[2]+porder};
    *wVector_global_    = (double*)malloc(     nbasis_global[0]*nbasis_global[1]*nbasis_global[2]*sizeof(double));
    *XVector_global_    = (double*)malloc(ndim*nbasis_global[0]*nbasis_global[1]*nbasis_global[2]*sizeof(double));

    double *wVector_global = *wVector_global_;
    double *XVector_global = *XVector_global_;

    assert(porder==2);
    assert(R>r);
    double Xs[]={1,1,0,-1,-1,-1,0,1,1,1};
    double Ys[]={0,1,1,1,0,-1,-1,-1,0,1};
    int ia=0;
    for (int i=0; i<nbasis_global[0]; i++) {
    for (int j=0; j<nbasis_global[1]; j++) {    
    for (int k=0; k<nbasis_global[2]; k++) {
        wVector_global[ia]=( i%2==0 ? 1.0 : 1.0/sqrt(2.0) );
        XVector_global[ndim*ia+0]=Xs[i]*(R*(2-j)/2.0+r*j/2.0);
        XVector_global[ndim*ia+1]=Ys[i]*(R*(2-j)/2.0+r*j/2.0);
        XVector_global[ndim*ia+2]=th*k/2.0;
        ia++;
    }}}
}
