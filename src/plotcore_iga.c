#if defined (__with_IGAPPlot)

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "mathutil.h"

#include <mgl2/mgl_cf.h>
#include <mgl2/mpi.h>

#include "coreimpl.h"
#include "meshimpl.h"
#include "physimpl.h"
#include "bcimpl.h"
#include "commimpl.h"
#include "solnimpl.h"

#include "meshcore.h"
#include "plotcore.h"

void iga_plotfield(HMGL* gr, void *app_, int ifield, int uid[], double uscale, int nppe, const char *sch, const char *opt)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    /* ================ input parameters ================ */
    Core app=(Core)app_;

    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    int nfield           = phys->nfield;
    field_ptr field      = phys->field;
    
    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
    int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    double *wVector      = mesh->wVector;
    double *XVector      = mesh->XVector;
    //int nboun            = mesh->nboun;
    //int *boun            = mesh->boun;
    int nsendPart        = mesh->nsendPart;
    int *sendPart        = mesh->sendPart;
    int nrecvPart        = mesh->nrecvPart;
    int *recvPart        = mesh->recvPart;

    BC bc               = app->bc[0];
    //int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    //double *par_dirichlet= bc->par_dirichlet;
    //double *par_neumann  = bc->par_neumann;
    
    // ---- MPI
    Comm comm            = app->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    //int *globalIdx       = comm->globalIdx;

    Soln soln            = app->soln[0];
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = soln->U_hist;
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],50+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],50+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble  ================ */
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double u[(torder+1)*nddim*ndof];
    double field_[nfield];
    // ---- Mesh
    int ibpe;
    int nbpe = int_pow(porder+1,3);                   // number of active basis per element
    int ibasis;                                                    // local basis id
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    double wV[nbpe];
    double XV[ndim*nbpe];
    // ---- quadrature
    double xq[ndim];
    // ---- IGA
    double N[nddim*nbpe];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim];//,detG;
    int         ia;
    int nbasis_y_active=nknot[1]-porder-1;
    int nbasis_z_active=nknot[2]-porder-1;
    /* ================ set up gr ================ */
    double machine_tol=1.e-15;
    int boff=8*4*0;
    int nknot_p[]={0,0,0};
    for (int iknot_x=porder; iknot_x<nknot[0]-porder-1; iknot_x++) nknot_p[0] = ( knotVector[0][iknot_x+1] < knotVector[0][iknot_x]+1.e-12 ? nknot_p[0] : nknot_p[0]+1 );
    for (int iknot_y=porder; iknot_y<nknot[1]-porder-1; iknot_y++) nknot_p[1] = ( knotVector[1][iknot_y+1] < knotVector[1][iknot_y]+1.e-12 ? nknot_p[1] : nknot_p[1]+1 );
    for (int iknot_z=porder; iknot_z<nknot[2]-porder-1; iknot_z++) nknot_p[2] = ( knotVector[2][iknot_z+1] < knotVector[2][iknot_z]+1.e-12 ? nknot_p[2] : nknot_p[2]+1 );
    HMDT cplot=mgl_create_data_size(nppe*(nknot_p[0]-2*boff)+1,nppe*(nknot_p[1]-2*boff)+1,nppe*(nknot_p[2]-2*boff)+1);
    HMDT aplot=mgl_create_data_size(nppe*(nknot_p[0]-2*boff)+1,nppe*(nknot_p[1]-2*boff)+1,nppe*(nknot_p[2]-2*boff)+1);
    HMDT xplot=mgl_create_data_size(nppe*(nknot_p[0]-2*boff)+1,nppe*(nknot_p[1]-2*boff)+1,nppe*(nknot_p[2]-2*boff)+1);
    HMDT yplot=mgl_create_data_size(nppe*(nknot_p[0]-2*boff)+1,nppe*(nknot_p[1]-2*boff)+1,nppe*(nknot_p[2]-2*boff)+1);
    HMDT zplot=mgl_create_data_size(nppe*(nknot_p[0]-2*boff)+1,nppe*(nknot_p[1]-2*boff)+1,nppe*(nknot_p[2]-2*boff)+1);

    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    // ---- loop over elements
    int ielem_x, ielem_y, ielem_z;
    int iknot_xp,iknot_yp,iknot_zp;
    iknot_xp=0;
    for (int iknot_x=porder; iknot_x<nknot[0]-porder-1; iknot_x++) {
        X0[0]=knotVector[0][iknot_x]; X1[0]=knotVector[0][iknot_x+1];
        if (X1[0]<X0[0]+1.e-12) continue;
    iknot_yp=0;
    for (int iknot_y=porder; iknot_y<nknot[1]-porder-1; iknot_y++) {
        X0[1]=knotVector[1][iknot_y]; X1[1]=knotVector[1][iknot_y+1];
        if (X1[1]<X0[1]+1.e-12) continue;
    iknot_zp=0;
    for (int iknot_z=porder; iknot_z<nknot[2]-porder-1; iknot_z++) {
        X0[2]=knotVector[2][iknot_z]; X1[2]=knotVector[2][iknot_z+1];
        if (X1[2]<X0[2]+1.e-12) continue;
    ielem_x=iknot_x-porder;
    ielem_y=iknot_y-porder;
    ielem_z=iknot_z-porder;
    kV[0]=knotVector[0]+iknot_x-porder;
    kV[1]=knotVector[1]+iknot_y-porder;
    kV[2]=knotVector[2]+iknot_z-porder;
    //for (int ielem_x=boff; ielem_x<(nelem[0]-boff); ielem_x++) {
    //for (int ielem_y=boff; ielem_y<(nelem[1]-boff); ielem_y++) {
    //for (int ielem_z=boff; ielem_z<(nelem[2]-boff); ielem_z++) {
    //X0[0]=knotVector[0][ielem_x+porder]; X1[0]=knotVector[0][ielem_x+porder+1];
    //X0[1]=knotVector[1][ielem_y+porder]; X1[1]=knotVector[1][ielem_y+porder+1];
    //X0[2]=knotVector[2][ielem_z+porder]; X1[2]=knotVector[2][ielem_z+porder+1];
    //kV[0]=knotVector[0]+ielem_x;
    //kV[1]=knotVector[1]+ielem_y;
    //kV[2]=knotVector[2]+ielem_z;

    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
        wV[ia++]=wVector[(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+(nbasis_z_active)*(ielem_y+j)+(ielem_z+k)];
    }}}
    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
    ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+nbasis_z_active*(ielem_y+j)+(ielem_z+k);
    for (int idim=0; idim<ndim; idim++) {
        XV[ia++]=XVector[ndim*ibasis+idim];
    }}}}
    // ---- loop over quadrature points
    for (int ippe_x=0; ippe_x<nppe+(int)(iknot_xp==(nknot_p[0]-boff)-1); ippe_x++) {
    for (int ippe_y=0; ippe_y<nppe+(int)(iknot_yp==(nknot_p[1]-boff)-1); ippe_y++) {
    for (int ippe_z=0; ippe_z<nppe+(int)(iknot_zp==(nknot_p[2]-boff)-1); ippe_z++) {
        // ---- evaluate quadrature coord.
        xq[0]=X0[0]+machine_tol+(X1[0]-X0[0])*ippe_x/nppe-2*machine_tol*(int)(ippe_x==nppe);
        xq[1]=X0[1]+machine_tol+(X1[1]-X0[1])*ippe_y/nppe-2*machine_tol*(int)(ippe_y==nppe);
        xq[2]=X0[2]+machine_tol+(X1[2]-X0[2])*ippe_z/nppe-2*machine_tol*(int)(ippe_z==nppe);
        // ---- evaluate N
        eval_Bspline(N,kV,porder,xq);
        // W
        for (int iddim=0; iddim<nddim; iddim++) W[iddim]=0.0;
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { W[iddim]+=wV[ibpe]*N[nddim*ibpe+iddim]; } }
        // R (parametric)
	for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { R[nddim*ibpe+iddim] =N[nddim*ibpe+iddim]*wV[ibpe]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=1; iddim<nddim; iddim++) { R[nddim*ibpe+iddim]-=R[nddim*ibpe+0]*W[iddim]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            for (int iddim=1; iddim<1+ndim; iddim++) {
            for (int jddim=iddim; jddim<1+ndim; jddim++) {
                R[ia++]-=(R[nddim*ibpe+jddim]*W[iddim]+R[nddim*ibpe+iddim]*W[jddim])/W[0];
            }}
        }
        // X
        for (int i=0; i<ndim*nddim; i++) { X[i]=0.0; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                for (int iddim=0; iddim<nddim; iddim++) { X[nddim*idim+iddim]+=XV[ndim*ibpe+idim]*R[nddim*ibpe+iddim]; }
            }
        }
        // G
        for (int idim=0; idim<ndim; idim++) { for (int jdim=0; jdim<ndim; jdim++) { G[ndim*idim+jdim]=X[nddim*idim+(1+jdim)]; } }
        //detG=G[0]*(G[4]*G[8]-G[5]*G[7])-G[1]*(G[3]*G[8]-G[5]*G[6])+G[2]*(G[3]*G[7]-G[4]*G[6]);
        matinv(G,ndim);
        // R (reference)
        for (int i=0; i<nddim*nbpe; i++) { N[i]=R[i]; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                R[nddim*ibpe+(1+idim)]=0.0;
                for (int jdim=0; jdim<ndim; jdim++)
                {
                    R[nddim*ibpe+(1+idim)]+=N[nddim*ibpe+(1+jdim)]*G[ndim*jdim+idim];
                }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int iddim=1+ndim; iddim<nddim; iddim++)
            {
                for (int idim=0; idim<ndim; idim++) { N[nddim*ibpe+iddim]-=R[nddim*ibpe+(1+idim)]*X[nddim*idim+iddim]; }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            R[ia+0]=G[0]*(N[ia+0]*G[0]+N[ia+1]*G[3]+N[ia+2]*G[6])+G[3]*(N[ia+1]*G[0]+N[ia+3]*G[3]+N[ia+4]*G[6])+G[6]*(N[ia+2]*G[0]+N[ia+4]*G[3]+N[ia+5]*G[6]);
            R[ia+1]=G[0]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[3]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[6]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+2]=G[0]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[3]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[6]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+3]=G[1]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[4]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[7]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+4]=G[1]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[4]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[7]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+5]=G[2]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[5]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[8]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
        }
        // ---- evaluate u
        for (ia=0; ia<(torder+1)*nddim*ndof; ia++) { u[ia]=0.0; }
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
            for (int ihist=0; ihist<torder+1; ihist++)
            {
                for (int idof=0; idof<ndof; idof++)
                {
                    ia=ndof*ibasis+idof;
                    for (int iddim=0; iddim<nddim; iddim++) { u[(nddim*ndof)*ihist+nddim*idof+iddim]+=R[nddim*ibpe+iddim]*U[(nselfIdx+nrecvIdx)*ihist+ia]; }
                }
            }
        }}}
        for (int ihist=0; ihist<torder+1; ihist++)
        {
            for (int idof=0; idof<ndof; idof++) { u[(nddim*ndof)*ihist+nddim*idof+0]+=par_periodic[ndim*idof+0]*xq[0]+par_periodic[ndim*idof+1]*xq[1]+par_periodic[ndim*idof+2]*xq[2]; }
            for (int idof=0; idof<ndof; idof++) { for (int idim=0; idim<ndim; idim++) { u[(nddim*ndof)*ihist+nddim*idof+(1+idim)]+=par_periodic[ndim*idof+idim]; } }
        }
        // ---- evaluate energy
        field(field_,u,par_mat);
        mgl_data_set_value(cplot,field_[0],                  nppe*(iknot_xp-boff)+ippe_x,nppe*(iknot_yp-boff)+ippe_y,nppe*(iknot_zp-boff)+ippe_z); // *cvalPtr 0 <-> 1 : transparent <-> opaque (if used for transparency)
        mgl_data_set_value(aplot,field_[ifield],             nppe*(iknot_xp-boff)+ippe_x,nppe*(iknot_yp-boff)+ippe_y,nppe*(iknot_zp-boff)+ippe_z);
        mgl_data_set_value(xplot,X[nddim*0]+uscale*u[uid[0]],nppe*(iknot_xp-boff)+ippe_x,nppe*(iknot_yp-boff)+ippe_y,nppe*(iknot_zp-boff)+ippe_z);
        mgl_data_set_value(yplot,X[nddim*1]+uscale*u[uid[1]],nppe*(iknot_xp-boff)+ippe_x,nppe*(iknot_yp-boff)+ippe_y,nppe*(iknot_zp-boff)+ippe_z);
        mgl_data_set_value(zplot,X[nddim*2]+uscale*u[uid[2]],nppe*(iknot_xp-boff)+ippe_x,nppe*(iknot_yp-boff)+ippe_y,nppe*(iknot_zp-boff)+ippe_z);
    }}} // iquad

    iknot_zp++;
    }
    iknot_yp++;
    }
    iknot_xp++;
    } // iknot


    const char *str_xyz[]={"x","","z"};
    for (int idim=0; idim<ndim; idim++)
    {
        int    naval;
        double *aval;
        char   sch_[512]; strcpy(sch_,sch);
        int    space=(nknot_p[idim]-2*boff)*nppe;//upto nelem*nppe;
        naval=nppe*(nknot_p[idim]-2*boff)/space+1;
        aval=(double*)malloc(naval*sizeof(double));
        for (int iaval=0; iaval<naval; iaval++) { aval[iaval]=space*iaval; }
        for (int iaval=0; iaval<naval; iaval++) { mgl_dens3_xyz(gr[0],xplot,yplot,zplot,aplot,strcat(sch_,str_xyz[idim]),aval[iaval],opt); }
        free(aval);
    }
    // ---- finalize
    mgl_delete_data(cplot);
    mgl_delete_data(aplot);
    mgl_delete_data(xplot);
    mgl_delete_data(yplot);
    mgl_delete_data(zplot);
    free(U);
}

/*

    faceID's for plotting:

      ______________________
     /|                  3 /|
    / |                   / |
   /____5________________/  |
   |  |               2  |  |
   | 0|_  _  _  _  _  _  |_1|
   |  /                  |  /
   | /                   | /
   |/___4________________|/

*/

void iga_plotface(HMGL* gr, void *app_, int uid[], double uscale, int nppe, const char *sch, const char *opt)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    /* ================ input parameters ================ */
    Core app=(Core)app_;

    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    int nfield           = phys->nfield;
    field_ptr field      = phys->field;
    
    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
    //int *nbasis          = mesh->nbasis;
    //int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    double *wVector      = mesh->wVector;
    double *XVector      = mesh->XVector;
    //int nboun            = mesh->nboun;
    //int *boun            = mesh->boun;
    int nsendPart        = mesh->nsendPart;
    int *sendPart        = mesh->sendPart;
    int nrecvPart        = mesh->nrecvPart;
    int *recvPart        = mesh->recvPart;

    BC bc               = app->bc[0];
    //int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    //double *par_dirichlet= bc->par_dirichlet;
    //double *par_neumann  = bc->par_neumann;
    
    // ---- MPI
    Comm comm            = app->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    //int *globalIdx       = comm->globalIdx;

    Soln soln            = app->soln[0];
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = soln->U_hist;
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],50+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],50+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble  ================ */
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double u[(torder+1)*nddim*ndof];
    double field_[nfield];
    // ---- Mesh
    int ibpe;
    int nbpe = int_pow(porder+1,3);                   // number of active basis per element
    int ibasis;                                                    // local basis id
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    double wV[nbpe];
    double XV[ndim*nbpe];
    // ---- quadrature
    double xq[ndim];
    // ---- IGA
    double N[nddim*nbpe];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim];//,detG;
    int         ia;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    /* ================ set up gr ================ */
    double machine_tol=1.e-15;
    int boff=8*4*0;
    HMDT cplot=mgl_create_data_size(nppe+1,nppe+1,1);
    //HMDT aplot=mgl_create_data_size(nppe+1,1,1);
    HMDT xplot=mgl_create_data_size(nppe+1,nppe+1,1);
    HMDT yplot=mgl_create_data_size(nppe+1,nppe+1,1);
    HMDT zplot=mgl_create_data_size(nppe+1,nppe+1,1);

    /* ================================================================   VOLUME INTEGRAL   ================================================================ */

    // ---- loop over elements
    for (int ielem_x=boff; ielem_x<(nelem[0]-boff); ielem_x++) {
    for (int ielem_y=boff; ielem_y<(nelem[1]-boff); ielem_y++) {
    for (int ielem_z=boff; ielem_z<(nelem[2]-boff); ielem_z++) {
    X0[0]=knotVector[0][ielem_x+porder]; X1[0]=knotVector[0][ielem_x+porder+1];
    X0[1]=knotVector[1][ielem_y+porder]; X1[1]=knotVector[1][ielem_y+porder+1];
    X0[2]=knotVector[2][ielem_z+porder]; X1[2]=knotVector[2][ielem_z+porder+1];
    kV[0]=knotVector[0]+ielem_x;
    kV[1]=knotVector[1]+ielem_y;
    kV[2]=knotVector[2]+ielem_z;
    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
        wV[ia++]=wVector[(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+(nbasis_z_active)*(ielem_y+j)+(ielem_z+k)];
    }}}
    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
    ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+nbasis_z_active*(ielem_y+j)+(ielem_z+k);
    for (int idim=0; idim<ndim; idim++) {
        XV[ia++]=XVector[ndim*ibasis+idim];
    }}}}
    // ---- loop over quadrature points
    int plotface[6]; for (int i=0; i<6; i++) { plotface[i]=1; }
    if (ielem_x!=(nelem[0]-boff)-1) { plotface[1]=0; }
    if (ielem_y!=(nelem[1]-boff)-1) { plotface[3]=0; }
    if (ielem_z!=(nelem[2]-boff)-1) { plotface[5]=0; }
    for (int iplotface=0; iplotface<6; iplotface++) {
    if (plotface[iplotface]==1) {
    for (int ippe_x=0; ippe_x<nppe*((int)(iplotface/2!=0))+1; ippe_x++) {
    for (int ippe_y=0; ippe_y<nppe*((int)(iplotface/2!=1))+1; ippe_y++) {
    for (int ippe_z=0; ippe_z<nppe*((int)(iplotface/2!=2))+1; ippe_z++) {
        int ippe_w0;
        int ippe_w1;
        switch (iplotface/2)
        {
        case 0: ippe_w0=ippe_y; ippe_w1=ippe_z; break;
        case 1: ippe_w0=ippe_z; ippe_w1=ippe_x; break;
        case 2: ippe_w0=ippe_x; ippe_w1=ippe_y; break;
        }
        // ---- evaluate quadrature coord.
        int temp0=iplotface/2;
        xq[(temp0+1)%ndim]=X0[(temp0+1)%ndim]+machine_tol+(X1[(temp0+1)%ndim]-X0[(temp0+1)%ndim])*ippe_w0/nppe-2*machine_tol*(int)(ippe_w0==nppe);
        xq[(temp0+2)%ndim]=X0[(temp0+2)%ndim]+machine_tol+(X1[(temp0+2)%ndim]-X0[(temp0+2)%ndim])*ippe_w1/nppe-2*machine_tol*(int)(ippe_w1==nppe);
        int temp1=iplotface-2*temp0;
        switch (temp1)
        {
        case 0: xq[temp0]=X0[temp0]+machine_tol; break;
        case 1: xq[temp0]=X1[temp0]-machine_tol; break;
        }
        // ---- evaluate N
        eval_Bspline(N,kV,porder,xq);
        // W
        for (int iddim=0; iddim<nddim; iddim++) W[iddim]=0.0;
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { W[iddim]+=wV[ibpe]*N[nddim*ibpe+iddim]; } }
        // R (parametric)
	for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { R[nddim*ibpe+iddim] =N[nddim*ibpe+iddim]*wV[ibpe]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=1; iddim<nddim; iddim++) { R[nddim*ibpe+iddim]-=R[nddim*ibpe+0]*W[iddim]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            for (int iddim=1; iddim<1+ndim; iddim++) {
            for (int jddim=iddim; jddim<1+ndim; jddim++) {
                R[ia++]-=(R[nddim*ibpe+jddim]*W[iddim]+R[nddim*ibpe+iddim]*W[jddim])/W[0];
            }}
        }
        // X
        for (int i=0; i<ndim*nddim; i++) { X[i]=0.0; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                for (int iddim=0; iddim<nddim; iddim++) { X[nddim*idim+iddim]+=XV[ndim*ibpe+idim]*R[nddim*ibpe+iddim]; }
            }
        }
        // G
        for (int idim=0; idim<ndim; idim++) { for (int jdim=0; jdim<ndim; jdim++) { G[ndim*idim+jdim]=X[nddim*idim+(1+jdim)]; } }
        //detG=G[0]*(G[4]*G[8]-G[5]*G[7])-G[1]*(G[3]*G[8]-G[5]*G[6])+G[2]*(G[3]*G[7]-G[4]*G[6]);
        matinv(G,ndim);
        // R (reference)
        for (int i=0; i<nddim*nbpe; i++) { N[i]=R[i]; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                R[nddim*ibpe+(1+idim)]=0.0;
                for (int jdim=0; jdim<ndim; jdim++)
                {
                    R[nddim*ibpe+(1+idim)]+=N[nddim*ibpe+(1+jdim)]*G[ndim*jdim+idim];
                }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int iddim=1+ndim; iddim<nddim; iddim++)
            {
                for (int idim=0; idim<ndim; idim++) { N[nddim*ibpe+iddim]-=R[nddim*ibpe+(1+idim)]*X[nddim*idim+iddim]; }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            R[ia+0]=G[0]*(N[ia+0]*G[0]+N[ia+1]*G[3]+N[ia+2]*G[6])+G[3]*(N[ia+1]*G[0]+N[ia+3]*G[3]+N[ia+4]*G[6])+G[6]*(N[ia+2]*G[0]+N[ia+4]*G[3]+N[ia+5]*G[6]);
            R[ia+1]=G[0]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[3]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[6]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+2]=G[0]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[3]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[6]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+3]=G[1]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[4]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[7]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+4]=G[1]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[4]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[7]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+5]=G[2]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[5]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[8]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
        }
        // ---- evaluate u
        for (ia=0; ia<(torder+1)*nddim*ndof; ia++) { u[ia]=0.0; }
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
            for (int ihist=0; ihist<torder+1; ihist++)
            {
                for (int idof=0; idof<ndof; idof++)
                {
                    ia=ndof*ibasis+idof;
                    for (int iddim=0; iddim<nddim; iddim++) { u[(nddim*ndof)*ihist+nddim*idof+iddim]+=R[nddim*ibpe+iddim]*U[(nselfIdx+nrecvIdx)*ihist+ia]; }
                }
            }
        }}}
        for (int ihist=0; ihist<torder+1; ihist++)
        {
            for (int idof=0; idof<ndof; idof++) { u[(nddim*ndof)*ihist+nddim*idof+0]+=par_periodic[ndim*idof+0]*xq[0]+par_periodic[ndim*idof+1]*xq[1]+par_periodic[ndim*idof+2]*xq[2]; }
            for (int idof=0; idof<ndof; idof++) { for (int idim=0; idim<ndim; idim++) { u[(nddim*ndof)*ihist+nddim*idof+(1+idim)]+=par_periodic[ndim*idof+idim]; } }
        }
        // ---- evaluate energy
        //field(field_,u,par_mat);
        mgl_data_set_value(cplot,1.0,ippe_w0,ippe_w1,0); // *cvalPtr 0 <-> 1 : transparent <-> opaque (if used for transparency)
        //mgl_data_set_value(aplot,field_[32],nppe*(ielem_x-boff)+ippe_x,nppe*(ielem_y-boff)+ippe_y,nppe*(ielem_z-boff)+ippe_z);
        mgl_data_set_value(xplot,X[nddim*0]+uscale*u[uid[0]],ippe_w0,ippe_w1,0);
        mgl_data_set_value(yplot,X[nddim*1]+uscale*u[uid[1]],ippe_w0,ippe_w1,0);
        mgl_data_set_value(zplot,X[nddim*2]+uscale*u[uid[2]],ippe_w0,ippe_w1,0);
    }}} // ippe
    mgl_surfc_xy(gr[0],xplot,yplot,zplot,cplot,sch,"");
    } // if
    } // iplotedge
    }}} // ielem


    // ---- combine and output graphics
    // ---- finalize
    mgl_delete_data(cplot);
    //mgl_delete_data(aplot);
    mgl_delete_data(xplot);
    mgl_delete_data(yplot);
    mgl_delete_data(zplot);
    free(U);
}

/*

    edgeID's for plotting:

      __________3___________
     /|                    /|
    6 |                   7 |
   /__9_______1__________/  11
   |  |                  |  |
   |  |_ _ _ _ _2 _ _ _ _|_ |
   8  /                  10 /
   | 4                   | 5
   |/_________0__________|/

*/

void iga_plotedge(HMGL* gr, void *app_, int uid[], double uscale, int nppe, const char *sch, const char *opt)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    /* ================ input parameters ================ */
    Core app=(Core)app_;

    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    int nfield           = phys->nfield;
    field_ptr field      = phys->field;
    
    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
    //int *nbasis          = mesh->nbasis;
    //int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    double *wVector      = mesh->wVector;
    double *XVector      = mesh->XVector;
    //int nboun            = mesh->nboun;
    //int *boun            = mesh->boun;
    int nsendPart        = mesh->nsendPart;
    int *sendPart        = mesh->sendPart;
    int nrecvPart        = mesh->nrecvPart;
    int *recvPart        = mesh->recvPart;

    BC bc               = app->bc[0];
    //int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    //double *par_dirichlet= bc->par_dirichlet;
    //double *par_neumann  = bc->par_neumann;
    
    // ---- MPI
    Comm comm            = app->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    //int *globalIdx       = comm->globalIdx;

    Soln soln            = app->soln[0];
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = soln->U_hist;
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],50+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],50+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble  ================ */
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double u[(torder+1)*nddim*ndof];
    double field_[nfield];
    // ---- Mesh
    int ibpe;
    int nbpe = int_pow(porder+1,3);                   // number of active basis per element
    int ibasis;                                                    // local basis id
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    double wV[nbpe];
    double XV[ndim*nbpe];
    // ---- quadrature
    double xq[ndim];
    // ---- IGA
    double N[nddim*nbpe];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim];//,detG;
    int         ia;
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    /* ================ set up gr ================ */
    double machine_tol=1.e-15;
    int boff=8*4*0;
    //HMDT cplot=mgl_create_data_size(nppe+1,1,1);
    //HMDT aplot=mgl_create_data_size(nppe+1,1,1);
    HMDT xplot=mgl_create_data_size(nppe+1,1,1);
    HMDT yplot=mgl_create_data_size(nppe+1,1,1);
    HMDT zplot=mgl_create_data_size(nppe+1,1,1);

    /* ================================================================   VOLUME INTEGRAL   ================================================================ */

    // ---- loop over elements
    for (int ielem_x=boff; ielem_x<(nelem[0]-boff); ielem_x++) {
    for (int ielem_y=boff; ielem_y<(nelem[1]-boff); ielem_y++) {
    for (int ielem_z=boff; ielem_z<(nelem[2]-boff); ielem_z++) {
    X0[0]=knotVector[0][ielem_x+porder]; X1[0]=knotVector[0][ielem_x+porder+1];
    X0[1]=knotVector[1][ielem_y+porder]; X1[1]=knotVector[1][ielem_y+porder+1];
    X0[2]=knotVector[2][ielem_z+porder]; X1[2]=knotVector[2][ielem_z+porder+1];
    kV[0]=knotVector[0]+ielem_x;
    kV[1]=knotVector[1]+ielem_y;
    kV[2]=knotVector[2]+ielem_z;
    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
        wV[ia++]=wVector[(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+(nbasis_z_active)*(ielem_y+j)+(ielem_z+k)];
    }}}
    ia=0; 
    for (int i=0; i<porder+1; i++) {
    for (int j=0; j<porder+1; j++) {
    for (int k=0; k<porder+1; k++) {
    ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+i)+nbasis_z_active*(ielem_y+j)+(ielem_z+k);
    for (int idim=0; idim<ndim; idim++) {
        XV[ia++]=XVector[ndim*ibasis+idim];
    }}}}
    // ---- loop over quadrature points
    int plotedge[12]; for (int i=0; i<12; i++) { plotedge[i]=1; }
    if (ielem_x!=(nelem[0]-boff)-1) { plotedge[5]=0; plotedge[7]=0; plotedge[10]=0; plotedge[11]=0; }
    if (ielem_y!=(nelem[1]-boff)-1) { plotedge[2]=0; plotedge[3]=0; plotedge[9]=0; plotedge[11]=0; }
    if (ielem_z!=(nelem[2]-boff)-1) { plotedge[1]=0; plotedge[3]=0; plotedge[6]=0; plotedge[7]=0; }
    for (int iplotedge=0; iplotedge<12; iplotedge++) {
    if (plotedge[iplotedge]==1) {
    for (int ippe_x=0; ippe_x<nppe*((int)(iplotedge/4==0))+1; ippe_x++) {
    for (int ippe_y=0; ippe_y<nppe*((int)(iplotedge/4==1))+1; ippe_y++) {
    for (int ippe_z=0; ippe_z<nppe*((int)(iplotedge/4==2))+1; ippe_z++) {
        int ippe_w=ippe_x+ippe_y+ippe_z;
        // ---- evaluate quadrature coord.
        int temp0=iplotedge/4;
        xq[ temp0        ]=X0[ temp0        ]+machine_tol+(X1[ temp0        ]-X0[ temp0        ])*ippe_w/nppe-2*machine_tol*(int)(ippe_w==nppe);
        int temp1=iplotedge-4*temp0;
        switch (temp1)
        {
        case 0: xq[(temp0+1)%ndim]=X0[(temp0+1)%ndim]+machine_tol; xq[(temp0+2)%ndim]=X0[(temp0+2)%ndim]+machine_tol; break;
        case 1: xq[(temp0+1)%ndim]=X0[(temp0+1)%ndim]+machine_tol; xq[(temp0+2)%ndim]=X1[(temp0+2)%ndim]-machine_tol; break;
        case 2: xq[(temp0+1)%ndim]=X1[(temp0+1)%ndim]-machine_tol; xq[(temp0+2)%ndim]=X0[(temp0+2)%ndim]+machine_tol; break;
        case 3: xq[(temp0+1)%ndim]=X1[(temp0+1)%ndim]-machine_tol; xq[(temp0+2)%ndim]=X1[(temp0+2)%ndim]-machine_tol; break;
        }
        // ---- evaluate N
        eval_Bspline(N,kV,porder,xq);
        // W
        for (int iddim=0; iddim<nddim; iddim++) W[iddim]=0.0;
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { W[iddim]+=wV[ibpe]*N[nddim*ibpe+iddim]; } }
        // R (parametric)
	for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=0; iddim<nddim; iddim++) { R[nddim*ibpe+iddim] =N[nddim*ibpe+iddim]*wV[ibpe]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++) { for (int iddim=1; iddim<nddim; iddim++) { R[nddim*ibpe+iddim]-=R[nddim*ibpe+0]*W[iddim]/W[0]; } }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            for (int iddim=1; iddim<1+ndim; iddim++) {
            for (int jddim=iddim; jddim<1+ndim; jddim++) {
                R[ia++]-=(R[nddim*ibpe+jddim]*W[iddim]+R[nddim*ibpe+iddim]*W[jddim])/W[0];
            }}
        }
        // X
        for (int i=0; i<ndim*nddim; i++) { X[i]=0.0; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                for (int iddim=0; iddim<nddim; iddim++) { X[nddim*idim+iddim]+=XV[ndim*ibpe+idim]*R[nddim*ibpe+iddim]; }
            }
        }
        // G
        for (int idim=0; idim<ndim; idim++) { for (int jdim=0; jdim<ndim; jdim++) { G[ndim*idim+jdim]=X[nddim*idim+(1+jdim)]; } }
        //detG=G[0]*(G[4]*G[8]-G[5]*G[7])-G[1]*(G[3]*G[8]-G[5]*G[6])+G[2]*(G[3]*G[7]-G[4]*G[6]);
        matinv(G,ndim);
        // R (reference)
        for (int i=0; i<nddim*nbpe; i++) { N[i]=R[i]; }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int idim=0; idim<ndim; idim++)
            {
                R[nddim*ibpe+(1+idim)]=0.0;
                for (int jdim=0; jdim<ndim; jdim++)
                {
                    R[nddim*ibpe+(1+idim)]+=N[nddim*ibpe+(1+jdim)]*G[ndim*jdim+idim];
                }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            for (int iddim=1+ndim; iddim<nddim; iddim++)
            {
                for (int idim=0; idim<ndim; idim++) { N[nddim*ibpe+iddim]-=R[nddim*ibpe+(1+idim)]*X[nddim*idim+iddim]; }
            }
        }
        for (int ibpe=0; ibpe<nbpe; ibpe++)
        {
            ia=nddim*ibpe+(1+ndim);
            R[ia+0]=G[0]*(N[ia+0]*G[0]+N[ia+1]*G[3]+N[ia+2]*G[6])+G[3]*(N[ia+1]*G[0]+N[ia+3]*G[3]+N[ia+4]*G[6])+G[6]*(N[ia+2]*G[0]+N[ia+4]*G[3]+N[ia+5]*G[6]);
            R[ia+1]=G[0]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[3]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[6]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+2]=G[0]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[3]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[6]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+3]=G[1]*(N[ia+0]*G[1]+N[ia+1]*G[4]+N[ia+2]*G[7])+G[4]*(N[ia+1]*G[1]+N[ia+3]*G[4]+N[ia+4]*G[7])+G[7]*(N[ia+2]*G[1]+N[ia+4]*G[4]+N[ia+5]*G[7]);
            R[ia+4]=G[1]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[4]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[7]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
            R[ia+5]=G[2]*(N[ia+0]*G[2]+N[ia+1]*G[5]+N[ia+2]*G[8])+G[5]*(N[ia+1]*G[2]+N[ia+3]*G[5]+N[ia+4]*G[8])+G[8]*(N[ia+2]*G[2]+N[ia+4]*G[5]+N[ia+5]*G[8]);
        }
        // ---- evaluate u
        for (ia=0; ia<(torder+1)*nddim*ndof; ia++) { u[ia]=0.0; }
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
            for (int ihist=0; ihist<torder+1; ihist++)
            {
                for (int idof=0; idof<ndof; idof++)
                {
                    ia=ndof*ibasis+idof;
                    for (int iddim=0; iddim<nddim; iddim++) { u[(nddim*ndof)*ihist+nddim*idof+iddim]+=R[nddim*ibpe+iddim]*U[(nselfIdx+nrecvIdx)*ihist+ia]; }
                }
            }
        }}}
        for (int ihist=0; ihist<torder+1; ihist++)
        {
            for (int idof=0; idof<ndof; idof++) { u[(nddim*ndof)*ihist+nddim*idof+0]+=par_periodic[ndim*idof+0]*xq[0]+par_periodic[ndim*idof+1]*xq[1]+par_periodic[ndim*idof+2]*xq[2]; }
            for (int idof=0; idof<ndof; idof++) { for (int idim=0; idim<ndim; idim++) { u[(nddim*ndof)*ihist+nddim*idof+(1+idim)]+=par_periodic[ndim*idof+idim]; } }
        }
        // ---- evaluate energy
        //field(field_,u,par_mat);
        //mgl_data_set_value(cplot,field_[0],nppe*(ielem_x-boff)+ippe_x,nppe*(ielem_y-boff)+ippe_y,nppe*(ielem_z-boff)+ippe_z); // *cvalPtr 0 <-> 1 : transparent <-> opaque (if used for transparency)
        //mgl_data_set_value(aplot,field_[32],nppe*(ielem_x-boff)+ippe_x,nppe*(ielem_y-boff)+ippe_y,nppe*(ielem_z-boff)+ippe_z);
        mgl_data_set_value(xplot,X[nddim*0]+uscale*u[uid[0]],ippe_w,0,0);
        mgl_data_set_value(yplot,X[nddim*1]+uscale*u[uid[1]],ippe_w,0,0);
        mgl_data_set_value(zplot,X[nddim*2]+uscale*u[uid[2]],ippe_w,0,0);
    }}} // ippe
    mgl_plot_xyz(gr[0],xplot,yplot,zplot,sch,opt);
    } // if
    } // iplotedge
    }}} // ielem


    // ---- combine and output graphics
    // ---- finalize
    //mgl_delete_data(cplot);
    //mgl_delete_data(aplot);
    mgl_delete_data(xplot);
    mgl_delete_data(yplot);
    mgl_delete_data(zplot);
    free(U);
}
#endif
