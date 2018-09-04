#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"
#include "quadrature.h"

#include "coreimpl.h"
#include "meshcore.h"
#include "postcore.h"

void iga_intdensity(double HH_[], void *core_)
//void iga_intdensity(double HH_[], Core core)
{
    /* ================ input parameters ================ */
    Core core=(Core)core_;

    Phys phys            = core->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    int torder           = phys->torder;
    double *par_mat      = phys->par_mat;
    int ndensity         = phys->ndensity;
    density_ptr  density = phys->density;
    densityb_ptr densityb= phys->densityb;

    Mesh mesh            = core->mesh[0];
    int porder           = mesh->porder;
    int *nelem           = mesh->nelem;
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

    BC bc                = core->bc[0];
    int *bc_type         = bc->bc_type;
    double *par_periodic = bc->par_periodic;
    //double *par_dirichlet= bc->par_dirichlet;
    double *par_neumann  = bc->par_neumann;
    neumann_ptr neumann  = bc->neumann;

    Comm comm            = core->comm[0];
    int nselfIdx         = comm->nselfIdx;
    int *selfIdx         = comm->selfIdx;
    int nsendIdx         = comm->nsendIdx;
    int *sendIdx         = comm->sendIdx;
    int *sendPtr         = comm->sendPtr;
    int nrecvIdx         = comm->nrecvIdx;
    int *recvIdx         = comm->recvIdx;
    int *recvPtr         = comm->recvPtr;
    //int *globalIdx       = comm->globalIdx;

    Soln soln            = core->soln[0];
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
        { assert(MPI_Isend(sendbuff+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_DOUBLE,sendPart[isendPart],40+ihist,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS);}
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) 
        { assert(MPI_Irecv(recvbuff+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_DOUBLE,recvPart[irecvPart],40+ihist,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS);}
        for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
        for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
        for (int i=0; i<nrecvIdx; i++) { U[ihist*(nselfIdx+nrecvIdx)+recvIdx[i]]=recvbuff[i]; }
    }
    free(sendbuff);
    free(recvbuff);
    /* ================ assemble Residual_ ================ */
    // ---- PDE
    int nddim=1+ndim+ndim*(ndim+1)/2;
    double u[(torder+1)*nddim*ndof];
    double HH[ndensity]; for (int i=0; i<ndensity; i++) HH[i]=0.0;
    double hh[ndensity];
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
    double xi;
    double weight;
    // ---- IGA
    double N[nddim*nbpe];
    double R[nbpe*nddim];
    double W[nddim];
    double X[ndim*nddim];
    double G[ndim*ndim],detG;
    int         ia;
    //int nbasis_x_active=nelem[0]+porder;
    int nbasis_y_active=nknot[1]-porder-1;
    int nbasis_z_active=nknot[2]-porder-1;
    /* ================================================================   VOLUME INTEGRAL   ================================================================ */
    int    nquad=5;
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
        // ---- evaluate hh
        density(hh,u,par_mat);
        //energy=pow(u1[nddim*0]-u0[nddim*0],2.)+pow(u1[nddim*1]-u0[nddim*1],2.)+pow(u1[nddim*2]-u0[nddim*2],2.);
        // ---- add to HH
        weight=wquad[iquad_x]*wquad[iquad_y]*wquad[iquad_z]*(X1[0]-X0[0])*(X1[1]-X0[1])*(X1[2]-X0[2])*detG;
        for (int i=0; i<ndensity; i++) HH[i]+=weight*hh[i];
    }}} // iquad
    }}} // ielem

    /* ================================================================ BOUNDARY CONDITIONS ================================================================ */
    int ax_,ay_,az_;
    int ibasis_z_0,ibasis_z_1;
    /* ----------------------------------------------------------------     Neumann BCs     ---------------------------------------------------------------- */
    double      neumann_,ub[torder+1];
    double      X0_[ndim-1],X1_[ndim-1];
    int         nelem_[ndim];
    int         nknot_[ndim];
    double      *knotVector_[ndim];
    double      N_[2*int_pow(porder+1,2)];
    double      *xq_,*yq_,*zq_;
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
        default: exit(EXIT_FAILURE);
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
                    // 
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
                        ia=0; for (int i=0; i<torder+1; i++) ub[i]=0.0;
                        for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++){
                        for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++){
                            ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_0;
                            for (int i=0; i<torder+1; i++) { ub[i]+=N_[ia]*U[(nselfIdx+nrecvIdx)*i+ndof*ibasis+idof]; } ia++;
                        }}
                        neumann(&neumann_,par_neumann,boun[iboun],0,idof,xq);
                        densityb(hh,ub,neumann_,par_mat,0);
                        weight=wquad[iquad_x_]*wquad[iquad_y_]*(X1_[0]-X0_[0])*(X1_[1]-X0_[1]);
                        for (int i=0; i<ndensity; i++)  HH[i]-=weight*hh[i];
                    }}
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
                    // 
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
                        ia=0; for (int i=0; i<torder+1; i++) ub[i]=0.0;
                        for (int ibpe_x_=0; ibpe_x_<porder+1; ibpe_x_++) {
                        for (int ibpe_y_=0; ibpe_y_<porder+1; ibpe_y_++) {
                            ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_0; for (int i=0; i<torder+1; i++) { ub[i]+=N_[ia]*U[(nselfIdx+nrecvIdx)*i+ndof*ibasis+idof]; } ia++;
                            ibasis=ax_*(ielem_x_+ibpe_x_)+ay_*(ielem_y_+ibpe_y_)+az_*ibasis_z_1; for (int i=0; i<torder+1; i++) { ub[i]+=N_[ia]*U[(nselfIdx+nrecvIdx)*i+ndof*ibasis+idof]; } ia++;
                        }}
                        neumann(&neumann_,par_neumann,boun[iboun],1,idof,xq);
                        densityb(hh,ub,neumann_,par_mat,1);
                        weight=wquad[iquad_x_]*wquad[iquad_y_]*(X1_[0]-X0_[0])*(X1_[1]-X0_[1]);
                        for (int i=0; i<ndensity; i++) HH[i]-=weight*hh[i];
                    }}
                }} // ielem
            } // if bc_type
        } // for idof
    } // for iboun

    assert(MPI_Allreduce(HH,HH_,ndensity,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD)==MPI_SUCCESS);
    
    free(U);
}
