#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"

#include "meshimpl.h"
#include "meshcore.h"
#define MTOL 1.e-15

void mpipart_iga(Mesh mesh)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    //int ndim   = mesh->ndim;
    int ia;
    int ndim                   = 3;
    int porder                 = mesh->porder;
    int *bc_periodic           = mesh->bc_periodic;
    int *nelem_global          = mesh->nelem_global;
    double **knotVector_global = mesh->knotVector_global;

    int npart[ndim],ipart[ndim];
    int *nelem=(int*) malloc(ndim*sizeof(int));
    int *nknot=(int*) malloc(ndim*sizeof(int));
    double **knotVector=(double**)malloc(ndim*sizeof(double*));
    int ielem_displ[ndim];

    int nepp_approx;
    nepp_approx=round(pow((double)(nelem_global[0]*nelem_global[1]*nelem_global[2])/(double)nproc,1./3.));
    for (int i=nelem_global[0]/nepp_approx+1; i>0; i--)  { if (nproc%i ==0) { npart[0]=i; break;} } 
    nepp_approx=round(pow((double)(nelem_global[1]*nelem_global[2])/(double)(nproc/npart[0]),1./2.));
    for (int i=nelem_global[1]/nepp_approx+1; i>0; i--)  { if (nproc/npart[0]%i ==0) { npart[1]=i; break;} }
    npart[2]=nproc/npart[0]/npart[1];
    if (rank==0) { printf("npart[]={%d,%d,%d}\n",npart[0],npart[1],npart[2]); }
    if ( nelem_global[0]<npart[0]*porder || nelem_global[1]<npart[1]*porder || nelem_global[2]<npart[2]*porder)
    {
        if (rank==0) { printf("\n\nAt least one of the partitions contains too few nodes (control points) in at least one spatial direction.\n\n"); }
        exit(EXIT_SUCCESS);
    }

    ipart[0]=rank/(npart[1]*npart[2]);
    ipart[1]=(rank%(npart[1]*npart[2]))/npart[2];
    ipart[2]=rank%npart[2];

    int *nelem_w_array;       if (rank==0) { nelem_w_array      =(int*) malloc(nproc*sizeof(int)); } else { nelem_w_array      =(int*) malloc(0); }
    int *nknot_w_array;       if (rank==0) { nknot_w_array      =(int*) malloc(nproc*sizeof(int)); } else { nknot_w_array      =(int*) malloc(0); }
    int *ielem_w_displ_array; if (rank==0) { ielem_w_displ_array=(int*) malloc(nproc*sizeof(int)); } else { ielem_w_displ_array=(int*) malloc(0); }
    for (int idim=0; idim<ndim; idim++)
    {
        if (rank==0)
        {
            int *nelem_w_1d;
            int *ielem_w_displ_1d;
            int ipart_w;
            nelem_w_1d      =(int*) malloc(npart[idim]*sizeof(int));
            ielem_w_displ_1d=(int*) malloc(npart[idim]*sizeof(int));
            for (int i=0; i<npart[idim]; i++) { nelem_w_1d[i]=nelem_global[idim]/npart[idim]; }
            for (int i=0; i<nelem_global[idim]%npart[idim]; i++) { nelem_w_1d[i]++; }
            ielem_w_displ_1d[0]=0;
            for (int i=1; i<npart[idim]; i++) { ielem_w_displ_1d[i]=ielem_w_displ_1d[i-1]+nelem_w_1d[i-1]; }
            for (int iproc=0; iproc<nproc; iproc++)
            {
                switch (idim) {
                case 0: ipart_w=iproc/(npart[1]*npart[2]); break;
                case 1: ipart_w=(iproc%(npart[1]*npart[2]))/npart[2]; break;
                case 2: ipart_w=iproc%npart[2]; break;
                default: exit(EXIT_FAILURE);
                }
                nelem_w_array[iproc]=nelem_w_1d[ipart_w];
                nknot_w_array[iproc]=nelem_w_1d[ipart_w]+1+2*porder;
                ielem_w_displ_array[iproc]=ielem_w_displ_1d[ipart_w];
            }
            free(nelem_w_1d);
            free(ielem_w_displ_1d);
        }
        assert(MPI_Scatter(nelem_w_array,1,MPI_INT,nelem+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        assert(MPI_Scatter(nknot_w_array,1,MPI_INT,nknot+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        assert(MPI_Scatter(ielem_w_displ_array,1,MPI_INT,ielem_displ+idim,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        knotVector[idim]=(double*) malloc(nknot[idim]*sizeof(double));
        assert(MPI_Scatterv(knotVector_global[idim],nknot_w_array,ielem_w_displ_array,MPI_DOUBLE,knotVector[idim],nknot[idim],MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    }
    free(nelem_w_array);
    free(nknot_w_array);
    free(ielem_w_displ_array);

    // ---- reference (X)
    int nbasis_active[ndim]; for (int i=0; i<ndim; i++) { nbasis_active[i]=nelem[i]+porder; }
    int nbasis_global[ndim]; for (int i=0; i<ndim; i++) { nbasis_global[i]=nelem_global[i]+porder; }
    int nbuf,*nbuf_array,*nbuf_displ,nbuf_global=0,*buf_index_self,*buf_index;
    double *buf_global;
    // ---- wVector 
    nbuf=nbasis_active[0]*nbasis_active[1]*nbasis_active[2];
    if (rank==0) { nbuf_array=(int*)malloc(nproc*sizeof(int)); } else { nbuf_array=(int*)malloc(0); }
    if (rank==0) { nbuf_displ=(int*)malloc((nproc+1)*sizeof(int)); } else { nbuf_displ=(int*)malloc(0); }
    assert(MPI_Gather(&nbuf,1,MPI_INT,nbuf_array,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0) { nbuf_displ[0]=0; for (int i=1; i<nproc+1; i++) { nbuf_displ[i]=nbuf_displ[i-1]+nbuf_array[i-1]; } }
    if (rank==0) { nbuf_global=nbuf_displ[nproc]; }
    buf_index_self=(int*)malloc(nbuf*sizeof(int));
    ia=0;
    for (int i=0; i<nbasis_active[0]; i++) {
    for (int j=0; j<nbasis_active[1]; j++) {
    for (int k=0; k<nbasis_active[2]; k++) {
        buf_index_self[ia++]=nbasis_global[1]*nbasis_global[2]*(ielem_displ[0]+i)+nbasis_global[2]*(ielem_displ[1]+j)+(ielem_displ[2]+k);
    }}}
    if (rank==0) { buf_index =(int*)malloc(nbuf_global*sizeof(int)); }    else { buf_index =(int*)malloc(0); }
    if (rank==0) { buf_global=(double*)malloc(nbuf_global*sizeof(double)); } else { buf_global=(double*)malloc(0); }
    assert(MPI_Gatherv(buf_index_self,nbuf,MPI_INT,buf_index,nbuf_array,nbuf_displ,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    double *wVector_global=mesh->wVector_global;
    if (rank==0) { for (int i=0; i<nbuf_global; i++) { buf_global[i]=wVector_global[buf_index[i]]; } }
    double *wVector = (double*)malloc(nbuf*sizeof(double));
    assert(MPI_Scatterv(buf_global,nbuf_array,nbuf_displ,MPI_DOUBLE,wVector,nbuf,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(nbuf_array);
    free(nbuf_displ);
    free(buf_index_self);
    free(buf_index);
    free(buf_global);
    // ---- XVector
    nbuf=ndim*nbasis_active[0]*nbasis_active[1]*nbasis_active[2];
    if (rank==0) { nbuf_array=(int*)malloc(nproc*sizeof(int)); } else { nbuf_array=(int*)malloc(0); }
    if (rank==0) { nbuf_displ=(int*)malloc((nproc+1)*sizeof(int)); } else { nbuf_displ=(int*)malloc(0); }
    assert(MPI_Gather(&nbuf,1,MPI_INT,nbuf_array,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0) { nbuf_displ[0]=0; for (int i=1; i<nproc+1; i++) { nbuf_displ[i]=nbuf_displ[i-1]+nbuf_array[i-1]; } }
    if (rank==0) { nbuf_global=nbuf_displ[nproc]; }
    buf_index_self=(int*)malloc(nbuf*sizeof(int));
    ia=0;
    for (int i=0; i<nbasis_active[0]; i++) {
    for (int j=0; j<nbasis_active[1]; j++) {
    for (int k=0; k<nbasis_active[2]; k++) {
        for (int idim=0; idim<ndim; idim++) { buf_index_self[ia++]=ndim*(nbasis_global[1]*nbasis_global[2]*(ielem_displ[0]+i)+nbasis_global[2]*(ielem_displ[1]+j)+(ielem_displ[2]+k))+idim; }
    }}}
    if (rank==0) { buf_index =(int*)malloc(nbuf_global*sizeof(int)); }    else { buf_index =(int*)malloc(0); }
    if (rank==0) { buf_global=(double*)malloc(nbuf_global*sizeof(double)); } else { buf_global=(double*)malloc(0); }
    assert(MPI_Gatherv(buf_index_self,nbuf,MPI_INT,buf_index,nbuf_array,nbuf_displ,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    double *XVector_global=mesh->XVector_global;
    if (rank==0) { for (int i=0; i<nbuf_global; i++) { buf_global[i]=XVector_global[buf_index[i]]; } }
    double *XVector = (double*)malloc(nbuf*sizeof(double));
    assert(MPI_Scatterv(buf_global,nbuf_array,nbuf_displ,MPI_DOUBLE,XVector,nbuf,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(nbuf_array);
    free(nbuf_displ);
    free(buf_index_self);
    free(buf_index);
    free(buf_global);

    mesh->wVector = wVector;
    mesh->XVector = XVector;

    int cond0,cond1,cond2;
    // ---- send
    cond0=( bc_periodic[0]==1 || ipart[0] != 0 );
    cond1=( bc_periodic[1]==1 || ipart[1] != 0 );
    cond2=( bc_periodic[2]==1 || ipart[2] != 0 );
    int nsendPart=(int)(                   cond2 )
                 +(int)(          cond1          )
                 +(int)( cond0                   )
                 +(int)(          cond1 && cond2 )
                 +(int)( cond0          && cond2 )
                 +(int)( cond0 && cond1          )
                 +(int)( cond0 && cond1 && cond2 );
    int *sendPart =(int*)malloc(nsendPart*sizeof(int));
    int *sendBlock=(int*)malloc(nsendPart*sizeof(int));
    ia=0;
    if (                   cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=0; ia++; }
    if (          cond1          ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=1; ia++; }
    if ( cond0                   ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=2; ia++; }
    if (          cond1 && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]  +npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=3; ia++; }
    if ( cond0          && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]  +npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=4; ia++; }
    if ( cond0 && cond1          ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]  +npart[2])%npart[2]); sendBlock[ia]=5; ia++; }
    if ( cond0 && cond1 && cond2 ) { sendPart[ia]=npart[1]*npart[2]*((ipart[0]-1+npart[0])%npart[0])+npart[2]*((ipart[1]-1+npart[1])%npart[1])+((ipart[2]-1+npart[2])%npart[2]); sendBlock[ia]=6; ia++; }
    // ---- recv
    cond0=( bc_periodic[0]==1 || ipart[0] != npart[0]-1 );
    cond1=( bc_periodic[1]==1 || ipart[1] != npart[1]-1 );
    cond2=( bc_periodic[2]==1 || ipart[2] != npart[2]-1 );
    int nrecvPart=(int)(                   cond2 )
                 +(int)(          cond1          )
                 +(int)( cond0                   )
                 +(int)(          cond1 && cond2 )
                 +(int)( cond0          && cond2 )
                 +(int)( cond0 && cond1          )
                 +(int)( cond0 && cond1 && cond2 );
    int *recvPart =(int*)malloc(nrecvPart*sizeof(int));
    int *recvBlock=(int*)malloc(nrecvPart*sizeof(int));
    ia=0;
    if (                   cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=0; ia++; }
    if (          cond1          ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=1; ia++; }
    if ( cond0                   ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=2; ia++; }
    if (          cond1 && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]           )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=3; ia++; }
    if ( cond0          && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]           )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=4; ia++; }
    if ( cond0 && cond1          ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]           )%npart[2]); recvBlock[ia]=5; ia++; }
    if ( cond0 && cond1 && cond2 ) { recvPart[ia]=npart[1]*npart[2]*((ipart[0]+1         )%npart[0])+npart[2]*((ipart[1]+1         )%npart[1])+((ipart[2]+1         )%npart[2]); recvBlock[ia]=6; ia++; }
    
    // ----
    mesh->nsendPart = nsendPart;
    mesh->sendPart  = sendPart;
    mesh->sendBlock = sendBlock;
    mesh->nrecvPart = nrecvPart;
    mesh->recvPart  = recvPart;
    mesh->recvBlock = recvBlock;
    
    int nboun = (int)(ipart[0]==0 && bc_periodic[0]==0)+(int)(ipart[0]==npart[0]-1 && bc_periodic[0]==0)
              + (int)(ipart[1]==0 && bc_periodic[1]==0)+(int)(ipart[1]==npart[1]-1 && bc_periodic[1]==0)
              + (int)(ipart[2]==0 && bc_periodic[2]==0)+(int)(ipart[2]==npart[2]-1 && bc_periodic[2]==0);
    int *boun=(int*)malloc(nboun*sizeof(int));
    ia=0;
    if (ipart[0]==0 && bc_periodic[0]==0) { boun[ia++]=0; } if (ipart[0]==npart[0]-1 && bc_periodic[0]==0) { boun[ia++]=1; }
    if (ipart[1]==0 && bc_periodic[1]==0) { boun[ia++]=2; } if (ipart[1]==npart[1]-1 && bc_periodic[1]==0) { boun[ia++]=3; }
    if (ipart[2]==0 && bc_periodic[2]==0) { boun[ia++]=4; } if (ipart[2]==npart[2]-1 && bc_periodic[2]==0) { boun[ia++]=5; }

    int *nbasis=(int*) malloc(ndim*sizeof(int));
    for (int idim=0; idim<ndim; idim++) { nbasis[idim] = ( bc_periodic[idim]==0 && ipart[idim]==npart[idim]-1 ? nelem[idim]+porder : nelem[idim]); }

    // set universalIdx_temp
    mesh->universalIdx_temp=(int*)malloc(nbasis[0]*nbasis[1]*nbasis[2]*sizeof(int));
    ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
        mesh->universalIdx_temp[ia++]=(nelem_global[1]+((int)(bc_periodic[1]==0))*porder)*(nelem_global[2]+((int)(bc_periodic[2]==0))*porder)*(ielem_displ[0]+ibasis_x)+(nelem_global[2]+((int)(bc_periodic[2]==0))*porder)*(ielem_displ[1]+ibasis_y)+(ielem_displ[2]+ibasis_z);
    }}}

    mesh->nelem        = nelem;
    mesh->nbasis       = nbasis;
    mesh->nknot        = nknot;
    mesh->knotVector   = knotVector;
    mesh->nboun        = nboun;
    mesh->boun         = boun;

    free(mesh->bc_periodic);
    free(mesh->nelem_global);
    for (int idim=0; idim<ndim; idim++) { free(mesh->knotVector_global[idim]); }
    free(mesh->knotVector_global);
    free(mesh->wVector_global);
    free(mesh->XVector_global);

    // F-periodic B.C.s: x=X+(F0.X+u)
    //double par_periodic[] ={0.0000*((double)app.bc_periodic[0]),0.0000*((double)app.bc_periodic[1]),0.0000*((double)app.bc_periodic[2]),
    //                        0.0000*((double)app.bc_periodic[0]),0.0000*((double)app.bc_periodic[1]),0.0000*((double)app.bc_periodic[2]),
    //                        0.0000*((double)app.bc_periodic[0]),0.0000*((double)app.bc_periodic[1]),0.0000*((double)app.bc_periodic[2]) };
}

void eval_Bspline(double N[], double *kV[], int p, double xq[])
/*

requires p==2.
specificaclly for quadrature point evaluations (no end points).



kV[idim]=knotVector[idim]+ielem

Basic formulas:
                  /  1  if kV[i] =< xi < kV[i+1],
   N_{i,0}(xi) = | 
                  \  0  otherwise.

                  xi   - kV[i]                 kV[i+1+p] -   xi
   N_{i,p} =  ----------------- N_{i,p-1}  +  --------------------- N_{i+1,p-1}
               kV[i+p] - kV[i]                 kV[i+1+p] - kV[i+1]

1st-order derivative:
 dN_{i,p}             p                                 p
---------- = ----------------- N_{i,p-1}  -  --------------------- N_{i+1,p-1}
   dxi        kV[i+p] - kV[i]                 kV[i+1+p] - kV[i+1]

2nd-order derivative:
                                                                             d^2N_{i,p}
                                                                            ------------
                                                                               dxi^2 
                                                                        /                 \
                                                                       /                   \
                                                                      /                     \

                        p         dN_{i,p-1}                                                              p           dN_{i+1,p-1}
               ----------------- ---------                                       -             --------------------- ------------
                kV[i+p] - kV[i]     dxi                                                         kV[i+1+p] - kV[i+1]       dxi
                                   /   \                                                                                 /   \
                                  /     \                                                                               /     \
                                 /       \                                                                             /       \
                 p-1                             p-1                                           p-1                                   p-1
       ------------------- N_{i,p-2}  -  ------------------- N_{i+1,p-2}               ------------------- N_{i+1,p-2}  -  --------------------- N_{i+2,p-2}
        kV[i+p-1] - kV[i]                 kV[i+p] - kV[i+1]                             kV[i+p] - kV[i+1]                   kV[i+1+p] - kV[i+2]


supp(N_{i,p}) = supp(dN_{i,p}) = supp(d^2N_{i,p}) = [ kV[i] , kV[i+p+1] ]. (by induction)
Backward dependency:
N_{i,p}
N_{i,p-1} N_{i+1,p-1}
N_{i,p-2} N_{i+1,p-2} N_{i+2,p-2}
*/
{
    if (p!=2) exit(EXIT_FAILURE);
int ndim=3;
int nddim=10; 
    int ia,ibpe;
    int pplus1=p+1;
    double Np[(3*pplus1)*ndim];
    double N1[pplus1+1];
    double N0[pplus1+2]; for (int i=0; i<pplus1+2; i++) { N0[i]=0; } N0[p]=1;
    for (int idim=0; idim<ndim; idim++)
    {
        for (int i=0; i<pplus1+1; i++) { N1[i]                         =(xq[idim]-kV[idim][i])/(kV[idim][i+1]-kV[idim][i]+MTOL)*N0[i]+(kV[idim][i+2]-xq[idim])/(kV[idim][i+2]-kV[idim][i+1]+MTOL)*N0[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+i]         =(xq[idim]-kV[idim][i])/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]+(kV[idim][i+3]-xq[idim])/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+pplus1+i]  =                   2.0/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]-                     2.0/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
        for (int i=0; i<pplus1+1; i++) { N1[i]                         =                   1.0/(kV[idim][i+1]-kV[idim][i]+MTOL)*N0[i]-                     1.0/(kV[idim][i+2]-kV[idim][i+1]+MTOL)*N0[i+1]; }
        for (int i=0; i<pplus1; i++)   { Np[(3*pplus1)*idim+2*pplus1+i]=                   2.0/(kV[idim][i+2]-kV[idim][i]+MTOL)*N1[i]-                     2.0/(kV[idim][i+3]-kV[idim][i+1]+MTOL)*N1[i+1]; }
    }

        
        for (int i=0; i<pplus1; i++){
        for (int j=0; j<pplus1; j++){
        for (int k=0; k<pplus1; k++){
            ibpe=pplus1*pplus1*i+pplus1*j+k;
            ia=nddim*ibpe;
            N[ia+0]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+1]=Np[pplus1+i]  *Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+2]=Np[i]         *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+k];
            N[ia+3]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+pplus1+k];
            N[ia+4]=Np[2*pplus1+i]*Np[3*pplus1+j]         *Np[6*pplus1+k];
            N[ia+5]=Np[pplus1+i]  *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+k];
            N[ia+6]=Np[pplus1+i]  *Np[3*pplus1+j]         *Np[6*pplus1+pplus1+k];
            N[ia+7]=Np[i]         *Np[3*pplus1+2*pplus1+j]*Np[6*pplus1+k];
            N[ia+8]=Np[i]         *Np[3*pplus1+pplus1+j]  *Np[6*pplus1+pplus1+k];
            N[ia+9]=Np[i]         *Np[3*pplus1+j]         *Np[6*pplus1+2*pplus1+k];
        }}}
}

double evalN(double *kV, int nk, int k, int i, int p, double xi)
/*

evaluates N_{i,p}(\xi)　using recursion.

k  : order of derivative
kV : knotVector
nk : # of knots
p  : porder
i  : ibasis
xi : coordinate of the quadrature point

*/
{

if (k==0)
{

//if (i==(nk-p-1-1) && xi==kV[nk-1]) // due to imperfection of the definition.
if (i==(nk-p-1-1) && double_abs(xi-kV[nk-1])<1.e-15 && double_abs(kV[nk-1]-kV[nk-2])<1.e-15) // temp, due to imperfection of the definition.
{ return 1.0; }
else if (p==0)
{ return ((xi>=kV[i] && kV[i+1]>xi) ? 1.0:0.0); }
else if (kV[i+p]==kV[i] && kV[i+p+1]==kV[i+1])
{ return 0.0; }
else if (kV[i+p]==kV[i])
{ return (kV[i+p+1]-xi)/(kV[i+p+1]-kV[i+1])*evalN(kV,nk,k,i+1,p-1,xi); }
else if (kV[i+p+1]==kV[i+1])
{ return (xi-kV[i])/(kV[i+p]-kV[i])*evalN(kV,nk,k,i,p-1,xi); }
else
{ return (xi-kV[i])/(kV[i+p]-kV[i])*evalN(kV,nk,k,i,p-1,xi)+(kV[i+p+1]-xi)/(kV[i+p+1]-kV[i+1])*evalN(kV,nk,k,i+1,p-1,xi); }

}
else if ( k>0 && k<=p )
{
double alpha[k+1][k+1];
int icol, irow;
alpha[0][0]=1.0;
for (irow=1; irow<k+1; irow++)
{
    alpha[irow][0]=( (kV[i+p-irow+1]==kV[i]) ? 0.0:alpha[irow-1][0]/(kV[i+p-irow+1]-kV[i]) );
    for (icol=1; icol<irow; icol++)
    {
        alpha[irow][icol]=( (kV[i+p+icol-irow+1]==kV[i+icol]) ? 0.0:(alpha[irow-1][icol]-alpha[irow-1][icol-1])/(kV[i+p+icol-irow+1]-kV[i+icol]) );
    }
    alpha[irow][irow]=( (kV[i+p+1]==kV[i+irow]) ? 0.0:-alpha[irow-1][irow-1]/(kV[i+p+1]-kV[i+irow]) );
}

double sum=0.0;
for (int j=0; j<k+1; j++){ sum+=alpha[k][j]*evalN(kV,nk,0,i+j,p-k,xi); }
return sum*factorial(p)/factorial(p-k);
}
else
{ printf("error in evalN.c: 0 =< k <= p has to be satisfied.\n"); return 0.0; }

} // evalN

void nurbs_knotinsert(double *U_universal, const int nelem[], int nref, int porder, int ndim, int ndof, int bc_periodic[])
{
    assert(porder==2); // currently only supports porder==2
    double vx[4],vy[4],vz[4];
    int nelem_old[ndim];
    int nbasis_old[ndim];
    int nbasis_new[ndim];
    int nbasis[ndim];
    double *Ui_new,*Ui_old;

    Ui_new=(double*) malloc((nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder)*sizeof(double));
    Ui_old=(double*) malloc((nelem[0]+porder)*(nelem[1]+porder)*(nelem[2]+porder)*sizeof(double));
    for (int idof=0; idof<ndof; idof++)
    {
        for (int idim=0; idim<ndim; idim++) {
            nelem_old[idim] =nelem[idim]/int_pow(2,nref);
            nbasis_old[idim]=nelem_old[idim]+porder;
            nbasis[idim]    =nelem_old[idim]+(bc_periodic[idim]==0)*porder;
        }
        for (int i=0; i<nbasis_old[0]; i++) { 
        for (int j=0; j<nbasis_old[1]; j++) { 
        for (int k=0; k<nbasis_old[2]; k++) { 
            Ui_old[nbasis_old[1]*nbasis_old[2]*i+nbasis_old[2]*j+k]=U_universal[ndof*(nbasis[1]*nbasis[2]*(i%nbasis[0])+nbasis[2]*(j%nbasis[1])+(k%nbasis[2]))+idof];
        }}}
        for (int iref=0; iref<nref; iref++)
        {
            for (int idim=0; idim<ndim; idim++) { nbasis_new[idim]=2*nelem_old[idim]+porder; }
            for (int i=0; i<nbasis_new[0]*nbasis_new[1]*nbasis_new[2]; i++) { Ui_new[i]=0; }

            for (int ix=0; ix<nbasis_old[0]; ix++) {
            for (int iy=0; iy<nbasis_old[1]; iy++) {
            for (int iz=0; iz<nbasis_old[2]; iz++) {
                switch (bc_periodic[0]) {
                case 0:
                    if (ix==0)                    { vx[0]=0    ; vx[1]=0    ; vx[2]=1    ; vx[3]=1./2.; }
                    else if (ix==1)               { vx[0]=0    ; vx[1]=1./2.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix<nelem_old[0])     { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0])    { vx[0]=1./4.; vx[1]=3./4.; vx[2]=1./2.; vx[3]=0    ; }
                    else if (ix==nelem_old[0]+1)  { vx[0]=1./2.; vx[1]=1    ; vx[2]=0    ; vx[3]=0    ; }
                    break;
                case 1:
                    if (ix==0)                    { vx[0]=0    ; vx[1]=0    ; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==1)               { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix<nelem_old[0])     { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0])    { vx[0]=1./4.; vx[1]=3./4.; vx[2]=3./4.; vx[3]=1./4.; }
                    else if (ix==nelem_old[0]+1)  { vx[0]=1./4.; vx[1]=3./4 ; vx[2]=0    ; vx[3]=0    ; }
                    break;
                }
                switch (bc_periodic[1]) {
                case 0:
                    if (iy==0)                    { vy[0]=0    ; vy[1]=0    ; vy[2]=1    ; vy[3]=1./2.; }
                    else if (iy==1)               { vy[0]=0    ; vy[1]=1./2.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy<nelem_old[1])     { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1])    { vy[0]=1./4.; vy[1]=3./4.; vy[2]=1./2.; vy[3]=0    ; }
                    else if (iy==nelem_old[1]+1)  { vy[0]=1./2.; vy[1]=1    ; vy[2]=0    ; vy[3]=0    ; }
                    break;
                case 1:
                    if (iy==0)                    { vy[0]=0    ; vy[1]=0    ; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==1)               { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy<nelem_old[1])     { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1])    { vy[0]=1./4.; vy[1]=3./4.; vy[2]=3./4.; vy[3]=1./4.; }
                    else if (iy==nelem_old[1]+1)  { vy[0]=1./4.; vy[1]=3./4.; vy[2]=0    ; vy[3]=0    ; }
                    break;
                }
                switch (bc_periodic[2]) {
                case 0:
                    if (iz==0)                    { vz[0]=0    ; vz[1]=0    ; vz[2]=1    ; vz[3]=1./2.; }
                    else if (iz==1)               { vz[0]=0    ; vz[1]=1./2.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz<nelem_old[2])     { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2])    { vz[0]=1./4.; vz[1]=3./4.; vz[2]=1./2.; vz[3]=0    ; }
                    else if (iz==nelem_old[2]+1)  { vz[0]=1./2.; vz[1]=1    ; vz[2]=0    ; vz[3]=0    ; }
                    break;
                case 1:
                    if (iz==0)                    { vz[0]=0    ; vz[1]=0    ; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==1)               { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz<nelem_old[2])     { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2])    { vz[0]=1./4.; vz[1]=3./4.; vz[2]=3./4.; vz[3]=1./4.; }
                    else if (iz==nelem_old[2]+1)  { vz[0]=1./4.; vz[1]=3./4.; vz[2]=0    ; vz[3]=0    ; }
                    break;
                }
                for (int jx=0; jx<4; jx++) {
                for (int jy=0; jy<4; jy++) {
                for (int jz=0; jz<4; jz++) {
                if (vx[jx]*vy[jy]*vz[jz] != 0) 
                { Ui_new[(nbasis_new[1]*nbasis_new[2])*(2*ix-2+jx)+nbasis_new[2]*(2*iy-2+jy)+(2*iz-2+jz)]+=Ui_old[(nbasis_old[1]*nbasis_old[2])*ix+nbasis_old[2]*iy+iz]*vx[jx]*vy[jy]*vz[jz]; }
                }}}
            }}}
            for (int i=0; i<nbasis_new[0]*nbasis_new[1]*nbasis_new[2]; i++) { Ui_old[i]=Ui_new[i]; }
            for (int idim=0; idim<ndim; idim++) { nelem_old[idim]*=2; }
            for (int idim=0; idim<ndim; idim++) { nbasis_old[idim]=nelem_old[idim]+porder; }
        }
        for (int idim=0; idim<ndim; idim++) { nbasis[idim]=nelem_old[idim]+(bc_periodic[idim]==0)*porder; }
        for (int i=0; i<nbasis[0]; i++) { 
        for (int j=0; j<nbasis[1]; j++) { 
        for (int k=0; k<nbasis[2]; k++) { 
            U_universal[ndof*(nbasis[1]*nbasis[2]*i+nbasis[2]*j+k)+idof]=Ui_old[nbasis_old[1]*nbasis_old[2]*i+nbasis_old[2]*j+k];
        }}}
    } // for ndof
    free(Ui_new);
    free(Ui_old);
}
