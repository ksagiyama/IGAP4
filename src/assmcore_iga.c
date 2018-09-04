#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "mathutil.h"
#include "quadrature.h"

#include "coreimpl.h"
#include "meshcore.h"
#include "assmcore.h"

#undef __FUNCT__
#define __FUNCT__ "Residual_iga"
PetscErrorCode Residual_iga(SNES snes, Vec U_, Vec RR, void *app_)
{
    /* ================ input parameters ================ */
    Core app=(Core)app_;

    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    residual_ptr residual= phys->residual;
    
    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
    int *nbasis          = mesh->nbasis;
    int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    double *wVector      = mesh->wVector;
    double *XVector      = mesh->XVector;
    int nboun            = mesh->nboun;
    int *boun            = mesh->boun;
    int nsendPart        = mesh->nsendPart;
    int *sendPart        = mesh->sendPart;
    int nrecvPart        = mesh->nrecvPart;
    int *recvPart        = mesh->recvPart;

    BC bc                = app->bc[0];
    int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    double *par_dirichlet= bc->par_dirichlet;
    double *par_neumann  = bc->par_neumann;
    neumann_ptr neumann  = bc->neumann;

    Comm comm            = app->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    int *globalIdx       = comm->globalIdx;

    Soln soln            = app->soln[0];
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = soln->U_hist;
    const double *U_curr; VecGetArrayRead(U_,&U_curr); for (int i=0; i<nselfIdx; i++) { U_hist[i]=U_curr[i]; } VecRestoreArrayRead(U_,&U_curr);
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],20+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],20+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble RR ================ */
    // ---- Mesh
    int ibpe;
    int nbpe = int_pow(porder+1,3);                   // number of active basis per element
    int ibasis;                                                    // local basis id
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    double wV[nbpe];
    double XV[ndim*nbpe];
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double *rr = (double*) malloc(ndof*nddim*sizeof(double));
    double u[(torder+1)*nddim*ndof];
    // ---- quadrature
    double xq[ndim],xi;
    double weight;
    // ---- IGA
    double N[nbpe*nddim];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim],detG;
    VecSet(RR,0.0);
    double val[ndof*nbpe];
    int    idx[ndof*nbpe];
    double temp;
    int    ia;
    int nbasis_y_active=nknot[1]-porder-1;
    int nbasis_z_active=nknot[2]-porder-1;
    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    int    nquad=4;
    double cquad[nquad];
    double wquad[nquad];
    legendre_handle(cquad,wquad,nquad,0.,1.);
    // ---- loop over elements
    int ielem_x, ielem_y, ielem_z;
    for (int iknot_x=porder; iknot_x<nknot[0]-porder-1; iknot_x++) {
        X0[0]=knotVector[0][iknot_x]; X1[0]=knotVector[0][iknot_x+1];
        if (X1[0]<X0[0]+1.e-12) continue;
    for (int iknot_y=porder; iknot_y<nknot[1]-porder-1; iknot_y++) {
        X0[1]=knotVector[1][iknot_y]; X1[1]=knotVector[1][iknot_y+1];
        if (X1[1]<X0[1]+1.e-12) continue;
    for (int iknot_z=porder; iknot_z<nknot[2]-porder-1; iknot_z++) {
        X0[2]=knotVector[2][iknot_z]; X1[2]=knotVector[2][iknot_z+1];
        if (X1[2]<X0[2]+1.e-12) continue;
    ielem_x=iknot_x-porder;
    ielem_y=iknot_y-porder;
    ielem_z=iknot_z-porder;
    kV[0]=knotVector[0]+iknot_x-porder;
    kV[1]=knotVector[1]+iknot_y-porder;
    kV[2]=knotVector[2]+iknot_z-porder;
    for (ia=0; ia<ndof*nbpe; ia++) { val[ia]=0.0; }
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
    for (int iquad_x=0; iquad_x<nquad; iquad_x++) {
    for (int iquad_y=0; iquad_y<nquad; iquad_y++) {
    for (int iquad_z=0; iquad_z<nquad; iquad_z++) {
        // ---- evaluate quadrature coord.
        xi=cquad[iquad_x];xq[0]=(-xi+1.0)*X0[0]+xi*X1[0];
        xi=cquad[iquad_y];xq[1]=(-xi+1.0)*X0[1]+xi*X1[1];
        xi=cquad[iquad_z];xq[2]=(-xi+1.0)*X0[2]+xi*X1[2];
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
        detG=G[0]*(G[4]*G[8]-G[5]*G[7])-G[1]*(G[3]*G[8]-G[5]*G[6])+G[2]*(G[3]*G[7]-G[4]*G[6]);
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
        // ---- evaluate rr
        residual(rr,u,par_mat);
        // ---- add to val
        weight=wquad[iquad_x]*wquad[iquad_y]*wquad[iquad_z]*(X1[0]-X0[0])*(X1[1]-X0[1])*(X1[2]-X0[2])*detG;
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            for (int idof=0; idof<ndof; idof++)
            {
                temp=0.0;
                for (int iddim=0; iddim<nddim; iddim++) { temp+=R[nddim*ibpe+iddim]*rr[nddim*idof+iddim]; }
                val[ndof*ibpe+idof]+=weight*temp;
            }
        }}}
    }}} // iquad
    // ---- set idx
    for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
    for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
    for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
        ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
        ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
        for (int idof=0; idof<ndof; idof++)
        {
            idx[ndof*ibpe+idof]=globalIdx[ndof*ibasis+idof];
        }
    }}}
    // ---- add to RR
    VecSetValues(RR,ndof*nbpe,idx,val,ADD_VALUES);
    }}} // ielem

    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    int ax_,ay_,az_;
    int nbasis_[ndim],ibasis_z_0,ibasis_z_1;
    /* ----------------------------------------------------------------     Neumann BCs     ---------------------------------------------------------------- */
    double neumann_;
    double X0_[ndim-1],X1_[ndim-1];
    int    nelem_[ndim],nknot_[ndim];
    double *knotVector_[ndim];
    double N_[2*int_pow(porder+1,2)];
    double val_[2*int_pow(porder+1,2)];
    int    idx_[2*int_pow(porder+1,2)];
    double *xq_,*yq_,*zq_;
    for (int iboun=0; iboun<nboun; iboun++)
    { 
        // ---- map local -> global
        switch (boun[iboun]/2)
        {
        case 0: //  YZ-surface: x->yglobal, y->zglobal, z->xglobal
            nelem_[0]=nelem[1]; nknot_[0]=nknot[1]; knotVector_[0]=knotVector[1];
            nelem_[1]=nelem[2]; nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2];
            nelem_[2]=nelem[0]; nknot_[2]=nknot[0]; knotVector_[2]=knotVector[0];
            xq_=xq+1;
            yq_=xq+2;
            zq_=xq+0;
            ax_=nbasis_z_active;
            ay_=1;
            az_=nbasis_y_active*nbasis_z_active;
            break;
        case 1: // XZ-surface: x->zglobal, y->xglobal, z->yglobal
            nelem_[0]=nelem[0]; nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0];
            nelem_[1]=nelem[2]; nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2];
            nelem_[2]=nelem[1]; nknot_[2]=nknot[1]; knotVector_[2]=knotVector[1];
            xq_=xq+0;
            yq_=xq+2;
            zq_=xq+1;
            ax_=nbasis_y_active*nbasis_z_active;
            ay_=1;
            az_=nbasis_z_active;
            break;
        case 2: // XY-surface: x->xglobal, y->yglobal, z->zglobal
            nelem_[0]=nelem[0]; nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0];
            nelem_[1]=nelem[1]; nknot_[1]=nknot[1]; knotVector_[1]=knotVector[1];
            nelem_[2]=nelem[2]; nknot_[2]=nknot[2]; knotVector_[2]=knotVector[2];
            xq_=xq+0;
            yq_=xq+1;
            zq_=xq+2;
            ax_=nbasis_y_active*nbasis_z_active;
            ay_=nbasis_z_active;
            az_=1;
            break;
        default:
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            *zq_=knotVector_[2][0];
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            *zq_=knotVector_[2][nknot_[2]-1];
            ibasis_z_0=nelem_[2]+porder-1;
            ibasis_z_1=nelem_[2]+porder-1-1;
            break;
        default: exit(EXIT_FAILURE);
        }
        // ---- standard Neumann
        for (int idof=0; idof<ndof; idof++)
        {
            if (bc_type[(2*ndof)*boun[iboun]+0*ndof+idof]==1)
            {
                for (int ielem_x_=0; ielem_x_<nelem_[0]; ielem_x_++) {
                for (int ielem_y_=0; ielem_y_<nelem_[1]; ielem_y_++) {
                    X0_[0]=knotVector_[0][ielem_x_+porder]; X1_[0]=knotVector_[0][ielem_x_+porder+1];
                    X0_[1]=knotVector_[1][ielem_y_+porder]; X1_[1]=knotVector_[1][ielem_y_+porder+1];
                    // set val_
                    for (int i=0; i<int_pow(porder+1,2); i++) { val_[i]=0.0; }
                    for (int iquad_x_=0; iquad_x_<nquad; iquad_x_++) {
                    for (int iquad_y_=0; iquad_y_<nquad; iquad_y_++) {
                        xi=cquad[iquad_x_];*xq_=(-xi+1.0)*X0_[0]+xi*X1_[0];
                        xi=cquad[iquad_y_];*yq_=(-xi+1.0)*X0_[1]+xi*X1_[1];
                        ia=0;
                        for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++) {
                        for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++) {
                            N_[ia++]=evalN(knotVector_[0],nknot_[0],0,ielem_x_+ibpe_x_,porder,*xq_)
                                    *evalN(knotVector_[1],nknot_[1],0,ielem_y_+ibpe_y_,porder,*yq_);
                        }}
                        neumann(&neumann_,par_neumann,boun[iboun],0,idof,xq);
                        weight=wquad[iquad_x_]*wquad[iquad_y_]*(X1_[0]-X0_[0])*(X1_[1]-X0_[1]);
                        for (int i=0; i<int_pow(porder+1,2); i++) { val_[i]-=weight*N_[i]*neumann_; }
                    }}
                    // set idx_
                    ia=0;
                    for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++) {
                    for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++) {
                        ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_0; idx_[ia++]=globalIdx[ndof*ibasis+idof];
                    }}
                    // add val_ & idx_
                    VecSetValues(RR,int_pow(porder+1,2),idx_,val_,ADD_VALUES);
                }} // ielem
            } // if bc_type
        } // for idof
        // ---- high-order Neumann
        for (int idof=0; idof<ndof; idof++)
        {
            if (bc_type[(2*ndof)*boun[iboun]+1*ndof+idof]==1)
            {
                for (int ielem_x_=0; ielem_x_<nelem_[0]; ielem_x_++) {
                for (int ielem_y_=0; ielem_y_<nelem_[1]; ielem_y_++) {
                    X0_[0]=knotVector_[0][ielem_x_+porder]; X1_[0]=knotVector_[0][ielem_x_+porder+1];
                    X0_[1]=knotVector_[1][ielem_y_+porder]; X1_[1]=knotVector_[1][ielem_y_+porder+1];
                    // set val_
                    for (int i=0; i<2*int_pow(porder+1,2); i++) { val_[i]=0.0; }
                    for (int iquad_x_=0; iquad_x_<nquad; iquad_x_++) {
                    for (int iquad_y_=0; iquad_y_<nquad; iquad_y_++) {
                        xi=cquad[iquad_x_];*xq_=(-xi+1.0)*X0_[0]+xi*X1_[0];
                        xi=cquad[iquad_y_];*yq_=(-xi+1.0)*X0_[1]+xi*X1_[1];
                        ia=0;
                        for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++) {
                        for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++) {
                            N_[ia++]=evalN(knotVector_[0],nknot_[0],0,ielem_x_+ibpe_x_,porder,*xq_)
                                    *evalN(knotVector_[1],nknot_[1],0,ielem_y_+ibpe_y_,porder,*yq_)
                                    *evalN(knotVector_[2],nknot_[2],1,ibasis_z_0      ,porder,*zq_);
                            N_[ia++]=evalN(knotVector_[0],nknot_[0],0,ielem_x_+ibpe_x_,porder,*xq_)
                                    *evalN(knotVector_[1],nknot_[1],0,ielem_y_+ibpe_y_,porder,*yq_)
                                    *evalN(knotVector_[2],nknot_[2],1,ibasis_z_1      ,porder,*zq_);
                        }}
                        neumann(&neumann_,par_neumann,boun[iboun],1,idof,xq);
                        weight=wquad[iquad_x_]*wquad[iquad_y_]*(X1_[0]-X0_[0])*(X1_[1]-X0_[1]);
                        for (int i=0; i<2*int_pow(porder+1,2); i++) { val_[i]-=weight*N_[i]*neumann_; }
                    }}
                    // set idx_
                    ia=0;
                    for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++) {
                    for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++) {
                        ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_0; idx_[ia++]=globalIdx[ndof*ibasis+idof];
                        ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_1; idx_[ia++]=globalIdx[ndof*ibasis+idof];
                    }}
                    // add val_ & idx_
                    VecSetValues(RR,2*int_pow(porder+1,2),idx_,val_,ADD_VALUES);
                }} // ielem
            } // if bc_type
        } // for idof
    } // for iboun

    // ----
    VecAssemblyBegin(RR);
    VecAssemblyEnd(RR);

    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int nrow_dirichlet;
    int *row_dirichlet;
    double *val_dirichlet;
    double *dirichlet_;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nknot_[0]=nknot[1]; knotVector_[0]=knotVector[1]; nbasis_[0]=nbasis[1]; ax_=nbasis_z_active;
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1;
            nknot_[2]=nknot[0]; knotVector_[2]=knotVector[0]; nbasis_[2]=nbasis[0]; az_=nbasis_y_active*nbasis_z_active;
            break;
        case 1:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1; 
            nknot_[2]=nknot[1]; knotVector_[2]=knotVector[1]; nbasis_[2]=nbasis[1]; az_=nbasis_z_active;
            break;
        case 2:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nknot_[1]=nknot[1]; knotVector_[1]=knotVector[1]; nbasis_[1]=nbasis[1]; ay_=nbasis_z_active;
            nknot_[2]=nknot[2]; knotVector_[2]=knotVector[2]; nbasis_[2]=nbasis[2]; az_=1;
            break;
        default:
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_[2]-1;
            ibasis_z_1=nbasis_[2]-2;
            break;
        default: exit(EXIT_FAILURE);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+0*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                row_dirichlet = (int*) malloc(nrow_dirichlet*sizeof(int));
                val_dirichlet = (double*) malloc(nrow_dirichlet*sizeof(double));
                dirichlet_    = (double*) malloc(nrow_dirichlet*sizeof(double));
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    row_dirichlet[ia]=globalIdx[ndof*ibasis+idof];
                    val_dirichlet[ia]=       U[ndof*ibasis+idof];
                    ia++;
                }}
                dirichlet_iga(dirichlet_,par_dirichlet,boun[iboun],0,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                for (int i=0; i<nrow_dirichlet; i++) { val_dirichlet[i]-=dirichlet_[i]; }
                VecSetValues(RR,nrow_dirichlet,row_dirichlet,val_dirichlet,INSERT_VALUES);
                free(row_dirichlet);
                free(val_dirichlet);
                free(dirichlet_);
        }}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+1*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                row_dirichlet = (int*) malloc(nrow_dirichlet*sizeof(int));
                val_dirichlet = (double*) malloc(nrow_dirichlet*sizeof(double));
                dirichlet_    = (double*) malloc(nrow_dirichlet*sizeof(double));
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    row_dirichlet[ia]=globalIdx[ndof*ibasis+idof];
                    val_dirichlet[ia]=       U[ndof*ibasis+idof];
                    ia++;
                }}
                dirichlet_iga(dirichlet_,par_dirichlet,boun[iboun],1,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                for (int i=0; i<nrow_dirichlet; i++) { val_dirichlet[i]-=dirichlet_[i]; }
                VecSetValues(RR,nrow_dirichlet,row_dirichlet,val_dirichlet,INSERT_VALUES);
                free(row_dirichlet);
                free(val_dirichlet);
                free(dirichlet_);
        }}
    } // for iboun

    VecAssemblyBegin(RR);
    VecAssemblyEnd(RR);
    /* ---------------------------------------------------------------- FINAL ---------------------------------------------------------------- */
    free(U);
    free(rr);
    return 0;
}


#undef __FUNCT__
#define __FUNCT__ "Tangent_iga"
PetscErrorCode Tangent_iga(SNES snes, Vec U_, Mat TT, Mat Pmat_, void *app_)
{
    /* ================ input parameters ================ */
    Core app=(Core)app_;
    
    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    tangent_ptr tangent  = phys->tangent;
    
    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
    int *nbasis          = mesh->nbasis;
    int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    double *wVector      = mesh->wVector;
    double *XVector      = mesh->XVector;
    int nboun            = mesh->nboun;
    int *boun            = mesh->boun;
    int nsendPart        = mesh->nsendPart;
    int *sendPart        = mesh->sendPart;
    int nrecvPart        = mesh->nrecvPart;
    int *recvPart        = mesh->recvPart;

    BC bc                = app->bc[0];
    int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    //double *par_dirichlet= bc->par_dirichlet;
    //double *par_neumann  = bc->par_neumann;

    Comm comm            = app->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    int *globalIdx       = comm->globalIdx;

    Soln soln            = app->soln[0];
    /* ================ assemble U ================ */
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    double *sendbuff=(double*) malloc(nsendIdx*sizeof(double));
    double *recvbuff=(double*) malloc(nrecvIdx*sizeof(double));
    double *U = (double*) malloc((torder+1)*(nselfIdx+nrecvIdx)*sizeof(double));
    double *U_hist = soln->U_hist;
    const double *U_curr; VecGetArrayRead(U_,&U_curr); for (int i=0; i<nselfIdx; i++) { U_hist[i]=U_curr[i]; } VecRestoreArrayRead(U_,&U_curr);
    for (int ihist=0; ihist<torder+1; ihist++)
    {
        for (int i=0; i<nselfIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+selfIdx[i]]=U_hist[ihist*nselfIdx+i]; }
        for (int i=0; i<nsendIdx; i++) { sendbuff[i]=U[ihist*(nselfIdx+nrecvIdx)+sendIdx[i]]; }
        for (int isendPart=0; isendPart<nsendPart; isendPart++) 
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],30+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],30+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble TT ================ */
    // ---- Mesh
    int ibpe,lbpe;
    int nbpe = int_pow(porder+1,3);
    int ibasis;
    double X0[ndim],X1[ndim];
    double *kV[ndim];
    double wV[nbpe];
    double XV[ndim*nbpe];
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double *tt = (double*) malloc(int_pow(ndof*nddim,2)*sizeof(double));
    double u[(torder+1)*nddim*ndof];
    // ---- quadrature
    double xq[ndim],xi;
    double weight;
    // ---- IGA
    double N[nbpe*nddim];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim],detG;
    MatZeroEntries(Pmat_);
    double val[int_pow(ndof*nbpe,2)];
    int         idx[ndof*nbpe];
    double      temp;
    int         ia,row,col;
    int nbasis_y_active=nknot[1]-porder-1;
    int nbasis_z_active=nknot[2]-porder-1;
    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    int    nquad=4;
    double cquad[nquad];
    double wquad[nquad];
    legendre_handle(cquad,wquad,nquad,0.,1.);
    // ---- loop over elements
    int ielem_x, ielem_y, ielem_z;
    for (int iknot_x=porder; iknot_x<nknot[0]-porder-1; iknot_x++) {
        X0[0]=knotVector[0][iknot_x]; X1[0]=knotVector[0][iknot_x+1];
        if (X1[0]<X0[0]+1.e-12) continue;
    for (int iknot_y=porder; iknot_y<nknot[1]-porder-1; iknot_y++) {
        X0[1]=knotVector[1][iknot_y]; X1[1]=knotVector[1][iknot_y+1];
        if (X1[1]<X0[1]+1.e-12) continue;
    for (int iknot_z=porder; iknot_z<nknot[2]-porder-1; iknot_z++) {
        X0[2]=knotVector[2][iknot_z]; X1[2]=knotVector[2][iknot_z+1];
        if (X1[2]<X0[2]+1.e-12) continue;
    ielem_x=iknot_x-porder;
    ielem_y=iknot_y-porder;
    ielem_z=iknot_z-porder;
    kV[0]=knotVector[0]+iknot_x-porder;
    kV[1]=knotVector[1]+iknot_y-porder;
    kV[2]=knotVector[2]+iknot_z-porder;

    for (int i=0; i<int_pow(ndof*nbpe,2); i++) { val[i]=0.0; }
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
    for (int iquad_x=0; iquad_x<nquad; iquad_x++) {
    for (int iquad_y=0; iquad_y<nquad; iquad_y++) {
    for (int iquad_z=0; iquad_z<nquad; iquad_z++) {
        // ---- evaluate quadrature coord.
        xi=cquad[iquad_x];xq[0]=(-xi+1.0)*X0[0]+xi*X1[0];
        xi=cquad[iquad_y];xq[1]=(-xi+1.0)*X0[1]+xi*X1[1];
        xi=cquad[iquad_z];xq[2]=(-xi+1.0)*X0[2]+xi*X1[2];
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
        detG=G[0]*(G[4]*G[8]-G[5]*G[7])-G[1]*(G[3]*G[8]-G[5]*G[6])+G[2]*(G[3]*G[7]-G[4]*G[6]);
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
        // ---- evaluate tt
        tangent(tt,u,par_mat);
        // ---- add to val
        weight=wquad[iquad_x]*wquad[iquad_y]*wquad[iquad_z]*(X1[0]-X0[0])*(X1[1]-X0[1])*(X1[2]-X0[2])*detG;
        for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
        for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
        for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
            ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
            for (int idof=0; idof<ndof; idof++)
            {
                row=ndof*ibpe+idof;
                for (int lbpe_x=0; lbpe_x<porder+1; lbpe_x++){
                for (int lbpe_y=0; lbpe_y<porder+1; lbpe_y++){
                for (int lbpe_z=0; lbpe_z<porder+1; lbpe_z++){
                    lbpe=(porder+1)*(porder+1)*lbpe_x+(porder+1)*lbpe_y+lbpe_z;
                    for (int ldof=0; ldof<ndof; ldof++)
                    {
                        col=ndof*lbpe+ldof;
                        temp=0.0;
                        for (int iddim=0; iddim<nddim; iddim++)
                        {
                            for (int lddim=0; lddim<nddim; lddim++)
                            {
                                temp+=R[nddim*ibpe+iddim]*tt[(ndof*nddim)*(nddim*idof+iddim)+(nddim*ldof+lddim)]*R[nddim*lbpe+lddim];
                            } // lddim
                        } // iddim
                        val[(ndof*nbpe)*row+col]+=weight*temp;
                    } // ldof
                }}} // lbpe
            } // idof
        }}} // ibpe
    }}} // iquad
    // ---- set idx
    for (int ibpe_x=0; ibpe_x<porder+1; ibpe_x++){
    for (int ibpe_y=0; ibpe_y<porder+1; ibpe_y++){
    for (int ibpe_z=0; ibpe_z<porder+1; ibpe_z++){
        ibpe=(porder+1)*(porder+1)*ibpe_x+(porder+1)*ibpe_y+ibpe_z;
        ibasis=(nbasis_y_active*nbasis_z_active)*(ielem_x+ibpe_x)+nbasis_z_active*(ielem_y+ibpe_y)+(ielem_z+ibpe_z);
        for (int idof=0; idof<ndof; idof++)
        {
            idx[ndof*ibpe+idof]=globalIdx[ndof*ibasis+idof];
        }
    }}}
    // ---- add to Pmat_
    MatSetValues(Pmat_,ndof*nbpe,idx,ndof*nbpe,idx,val,ADD_VALUES);
    }}} // ielem
    // ---- 
    MatAssemblyBegin(Pmat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Pmat_,MAT_FINAL_ASSEMBLY);

/*Mat Atrans;
MatTranspose(Pmat_,MAT_INITIAL_MATRIX,&Atrans);
MatAXPY(Atrans,-1.0,Pmat_,SAME_NONZERO_PATTERN);
double nm;
MatNorm(Atrans,NORM_INFINITY,&nm);
printf("%e\n",nm);
*/
    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    int nbasis_x_,nbasis_y_,nbasis_z_,ibasis_z_0,ibasis_z_1;
    int ax_,ay_,az_;
    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int nrow_dirichlet=0;
    int itemp;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        if      (boun[iboun]/2==0) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[1]*nbasis[2]; }
        else if (boun[iboun]/2==1) { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[2]*nbasis[0]; }
        else                       { itemp=(2*ndof); for (int i=0; i<(2*ndof); i++) { itemp-=bc_type[boun[iboun]*(2*ndof)+i]; } nrow_dirichlet += itemp*nbasis[0]*nbasis[1]; }
    }
    int *row_dirichlet = (int*) malloc(nrow_dirichlet*sizeof(int));
    ia=0;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nbasis_x_=nbasis[1]; ax_=nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[0]; az_=nbasis_y_active*nbasis_z_active;
            break;
        case 1:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[2]; ay_=1;
            nbasis_z_=nbasis[1]; az_=nbasis_z_active;
            break;
        case 2:
            nbasis_x_=nbasis[0]; ax_=nbasis_y_active*nbasis_z_active;
            nbasis_y_=nbasis[1]; ay_=nbasis_z_active;
            nbasis_z_=nbasis[2]; az_=1;
            break;
        default:
            printf("error: \n");
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_z_-1;
            ibasis_z_1=nbasis_z_-2;
            break;
        default: exit(EXIT_FAILURE);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+0*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[boun[iboun]*(2*ndof)+1*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_x_; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_y_; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    row_dirichlet[ia++]=globalIdx[ndof*ibasis+idof];
        }}}}
    } // for iboun
    assert(ia==nrow_dirichlet);

    MatZeroRowsColumns(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // symmetrized
    //MatZeroRows(Pmat_,nrow_dirichlet,row_dirichlet,1.0,PETSC_NULL,PETSC_NULL); // unsymmetrized
    free(row_dirichlet);
    
    /* ---------------------------------------------------------------- FINAL ASSEMBLY ---------------------------------------------------------------- */
    if (TT != Pmat_)
    {
        MatAssemblyBegin(TT,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(TT,MAT_FINAL_ASSEMBLY);
    }

    free(U);
    free(tt);
    return 0;
}


void dirichlet_iga(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder)
{
    // ---- current parametrization of Dirichlet B.C.:
    // u_x = a0 X + a1 Y + a2 Z              |1    |   |a0 a1 a2|
    // u_y = a3 X + a4 Y + a5 Z    -->   F = |  1  | + |a3 a4 a5|
    // u_z = a6 X + a7 Y + a8 Z              |    1|   |a6 a7 a8|


    double Z_;
    //int ibasis_z_0;
    int ibasis_z_1;
    // ---- set Z_(=const.)
    switch (face_id%2)
    {
        case 0:
            Z_=knotVector_[ndim-1][0];
            //ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            Z_=knotVector_[ndim-1][nknot_[ndim-1]-1]-1.e-15; // need fixing for end values.
            //ibasis_z_0=nbasis_[2]-1;
            ibasis_z_1=nbasis_[2]-2;
            break;
        default: exit(EXIT_FAILURE);
    }
    // ---- 
    switch (bc_order)
    {
        case 0:
        {
            // ---- map par to the face-local coord.
            double par_[ndim];
            switch(face_id/2)
            {
                case 0: par_[0]=par[ndim*idof+1]; par_[1]=par[ndim*idof+2]; par_[2]=par[ndim*idof+0]; break; // YZ
                case 1: par_[0]=par[ndim*idof+0]; par_[1]=par[ndim*idof+2]; par_[2]=par[ndim*idof+1]; break; // XZ
                case 2: par_[0]=par[ndim*idof+0]; par_[1]=par[ndim*idof+1]; par_[2]=par[ndim*idof+2]; break; // XY
                default: exit(0);
            }
            // ---- apply affine deformation on a boundary: par_[0]*X_+par_[1]*Y_+par_[2]*Z_.
            double *dirichlet_1d[ndim-1];
            double vec_p[porder];
            double mat_p[porder*porder];
            double W_;
            for (int bdim=0; bdim<ndim-1; bdim++)
            {
                dirichlet_1d[bdim]=(double*)malloc(nbasis_[bdim]*sizeof(double));
                // 0th - (p-1)th deriv. at knotVector_[bdim][porder].
                W_=knotVector_[bdim][porder];
                for (int i=0; i<porder; i++) { vec_p[i]=0.0; } vec_p[0]=par_[bdim]*W_; vec_p[1]=par_[bdim];
                for (int i=0; i<porder; i++) { for (int j=0; j<porder; j++) { mat_p[porder*i+j]=evalN(knotVector_[bdim],nknot_[bdim],i,j,porder,W_); }}
                matinv_dot_vec(mat_p,vec_p,porder);
                for (int i=0; i<porder; i++) { dirichlet_1d[bdim][i]=vec_p[i]; }
                // 1st deriv.(=par_[bdim]=const.) at knotVector_[bdim][i].
                for (int i=porder; i<nbasis_[bdim]; i++) { dirichlet_1d[bdim][i]=dirichlet_1d[bdim][i-1]+(knotVector_[bdim][i+porder]-knotVector_[bdim][i])/((double)porder)*par_[bdim]; }
            }
            // ----
            for (int i=0; i<nbasis_[0]; i++) { for (int j=0; j<nbasis_[1]; j++) { dirichlet[nbasis_[1]*i+j]=dirichlet_1d[0][i]+dirichlet_1d[1][j]+par_[2]*Z_; }}
            // ----
            for (int bdim=0; bdim<ndim-1; bdim++) { free(dirichlet_1d[bdim]); }
            break;
        }
        case 1:
        {
            // ---- first, set u_{i,n}=0.
            dirichlet_iga(dirichlet,par,face_id,0,idof,knotVector_,nknot_,nbasis_,ndim,porder);
            // ---- then,  set u_{i,n}=some val.
            double m=0.0/evalN(knotVector_[ndim-1],nknot_[ndim-1],1,ibasis_z_1,porder,Z_);// evalN need fixing for end values...
            for (int i=0; i<nbasis_[0]; i++) { for (int j=0; j<nbasis_[1]; j++) { dirichlet[nbasis_[1]*i+j]+=m; }}
            break;
        }
    }
}


void dirichletinit_iga(double *U_, void *app_)
{
    /* ================ input parameters ================ */
    Core app=(Core)app_;

    Phys phys            = app->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    //double *par_mat      = app->par_mat;

    Mesh mesh            = app->mesh[0];
    int porder           = mesh->porder;
    //int *nelem           = mesh->nelem;
    int *nbasis          = mesh->nbasis;
    int *nknot           = mesh->nknot;
    double **knotVector  = mesh->knotVector;
    int nboun            = mesh->nboun;
    int *boun            = mesh->boun;

    BC bc                = app->bc[0];
    int *bc_type         = bc->bc_type;
    //double *par_periodic = bc->par_periodic;
    double *par_dirichlet= bc->par_dirichlet;
    //double *par_neumann  = bc->par_neumann;

    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    /* ----------------------------------------------------------------    Dirichlet BCs    ---------------------------------------------------------------- */
    int ia,ibasis;
    int ax_,ay_,az_;
    int nbasis_[ndim],ibasis_z_0,ibasis_z_1;
    int    nknot_[ndim];
    double *knotVector_[ndim];
    int nrow_dirichlet;
    double *dirichlet;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nknot_[0]=nknot[1]; knotVector_[0]=knotVector[1]; nbasis_[0]=nbasis[1]; ax_=nbasis[2];
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1;
            nknot_[2]=nknot[0]; knotVector_[2]=knotVector[0]; nbasis_[2]=nbasis[0]; az_=nbasis[1]*nbasis[2];
            break;
        case 1:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nknot_[1]=nknot[2]; knotVector_[1]=knotVector[2]; nbasis_[1]=nbasis[2]; ay_=1; 
            nknot_[2]=nknot[1]; knotVector_[2]=knotVector[1]; nbasis_[2]=nbasis[1]; az_=nbasis[2];
            break;
        case 2:
            nknot_[0]=nknot[0]; knotVector_[0]=knotVector[0]; nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nknot_[1]=nknot[1]; knotVector_[1]=knotVector[1]; nbasis_[1]=nbasis[1]; ay_=nbasis[2];
            nknot_[2]=nknot[2]; knotVector_[2]=knotVector[2]; nbasis_[2]=nbasis[2]; az_=1;
            break;
        default:
            exit(0);
        }
        switch (boun[iboun]%2)
        {
        case 0:
            ibasis_z_0=0;
            ibasis_z_1=1;
            break;
        case 1:
            ibasis_z_0=nbasis_[2]-1;
            ibasis_z_1=nbasis_[2]-2;
            break;
        default: exit(EXIT_FAILURE);
        }
        // standard
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+0*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                dirichlet     = (double*) malloc(nrow_dirichlet*sizeof(double));
                dirichlet_iga(dirichlet,par_dirichlet,boun[iboun],0,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    U_[ndof*ibasis+idof]=dirichlet[ia++];
                }}
                free(dirichlet);
        }}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+1*ndof+idof]==0) {
                nrow_dirichlet=nbasis_[0]*nbasis_[1];
                dirichlet     = (double*) malloc(nrow_dirichlet*sizeof(double));
                dirichlet_iga(dirichlet,par_dirichlet,boun[iboun],1,idof,knotVector_,nknot_,nbasis_,ndim,porder);
                ia=0;
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    U_[ndof*ibasis+idof]=dirichlet[ia++];
                }}
                free(dirichlet);
        }}
    } // for iboun
}

void neumann_default_iga(double *neumann_, double par_neumann[], int face_id, int bc_order, int idof, double xq[])
{
    int ndof=3;

    //                                          0                          |                          1                               // bc_order
    //                -------------------------------------------------------------------------------------------------------------------------
    //                        0        ,        1        ,        2        |        0        ,        1        ,        2             // idof
    //                -------------------------------------------------------------------------------------------------------------------------
    double neumann[]={        0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 0
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 1
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 2
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 3
                              0        ,        0        ,        0        ,        0        ,        0        ,        0        ,    // face 4
                              0        ,        0        ,        0        ,        0        ,        0        ,        0         };  // face 5

    *neumann_=neumann[(2*ndof)*face_id+ndof*bc_order+idof];
}

void setnz_iga(int *d_nz, int *o_nz, Core core)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int porder    = core->mesh[0]->porder;
    int *nbasis   = core->mesh[0]->nbasis;
    int nsendPart = core->mesh[0]->nsendPart;
    int *sendPart = core->mesh[0]->sendPart;
    int nrecvPart = core->mesh[0]->nrecvPart;
    //int *recvPart = core->mesh[0]->recvPart;
    int *sendBlock= core->mesh[0]->sendBlock;
    int *recvBlock= core->mesh[0]->recvBlock;
    int ndof      = core->phys[0]->ndof;

    int sd[]={0,0,0}; for (int i=0; i<nsendPart; i++) { if (sendBlock[i]<3) sd[sendBlock[i]]=1; }
    int rv[]={0,0,0}; for (int i=0; i<nrecvPart; i++) { if (recvBlock[i]<3) rv[recvBlock[i]]=1; }

    int ps[]={0,0,0}; //periodic and single_part
    for (int i=0; i<nsendPart; i++) { if ( sendBlock[i]<3 && sendPart[i]==rank) { ps[sendBlock[i]]=1; } }
    int ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
    for (int idof=0; idof<ndof; idof++) {
        d_nz[ia]=ndof*( ( (ps[2]==0 && ibasis_x < porder) ? ibasis_x : porder )+1+( (ps[2]==0 && ibasis_x > (nbasis[0]-1)-porder) ? (nbasis[0]-1)-ibasis_x : porder ) )
                     *( ( (ps[1]==0 && ibasis_y < porder) ? ibasis_y : porder )+1+( (ps[1]==0 && ibasis_y > (nbasis[1]-1)-porder) ? (nbasis[1]-1)-ibasis_y : porder ) )
                     *( ( (ps[0]==0 && ibasis_z < porder) ? ibasis_z : porder )+1+( (ps[0]==0 && ibasis_z > (nbasis[2]-1)-porder) ? (nbasis[2]-1)-ibasis_z : porder ) );
        o_nz[ia]=ndof*( ( (sd[2]==0 && ibasis_x < porder) ? ibasis_x : porder )+1+( (rv[2]==0 && ibasis_x > (nbasis[0]-1)-porder) ? (nbasis[0]-1)-ibasis_x : porder ) )
                     *( ( (sd[1]==0 && ibasis_y < porder) ? ibasis_y : porder )+1+( (rv[1]==0 && ibasis_y > (nbasis[1]-1)-porder) ? (nbasis[1]-1)-ibasis_y : porder ) )
                     *( ( (sd[0]==0 && ibasis_z < porder) ? ibasis_z : porder )+1+( (rv[0]==0 && ibasis_z > (nbasis[2]-1)-porder) ? (nbasis[2]-1)-ibasis_z : porder ) )
                 -d_nz[ia];
        ia++;
    }}}}
}
// if send to/recv from self: d_nz=d_nz+o_nz, o_nz=0.

