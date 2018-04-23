#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "igapimpl.h"

#include "meshcore.h"
#include "commcore.h"
#include "assmcore.h"
#include "postcore.h"
#include "plotcore.h"

#include "mathutil.h"
#include "ioutil.h"

#include "interface.h"

// ---- Initialize
void IGAPInitialize(IGAP *igap_ptr)
{
    *igap_ptr  = NULL;

    IGAP igap  = (IGAP)malloc(1*sizeof(struct __IGAP));
    igap->core = (Core)malloc(1*sizeof(struct __Core));

    igap->post = (Post)malloc(1*sizeof(struct __Post));
    igap->plot = (Plot)malloc(1*sizeof(struct __Plot));

    *igap_ptr  = igap;
}

// ---- Mesh
void IGAPMeshAdd(IGAP igap, Mesh mesh)
{
    igap->core->mesh=(Mesh*)malloc(1*sizeof(Mesh));
    igap->core->mesh[0]=mesh;
}

// ---- Phys
void IGAPPhysAdd(IGAP igap, Phys phys)
{
    igap->core->phys=(Phys*)malloc(1*sizeof(Phys));
    igap->core->phys[0]=phys;

    igap->post->HH=(double**)malloc(1*sizeof(double*));
    igap->post->HH[0]=(double*)malloc(igap->core->phys[0]->ndensity*sizeof(double));
}

// ---- BC
void IGAPBCAdd(IGAP igap, BC bc)
{
    igap->core->bc=(BC*)malloc(1*sizeof(BC));
    igap->core->bc[0]=bc;
}

void IGAPBCAddDefault(IGAP igap)
{
    BC bc = (BC)malloc(1*sizeof(struct __BC));

    int ndim = igap->core->phys[0]->ndim;
    int ndof = igap->core->phys[0]->ndof;

    bc->par_periodic  = (double*) malloc(ndim*ndof*sizeof(double)); for (int i=0; i<ndim*ndof; i++) bc->par_periodic[i]  = 0.0;
    bc->bc_type       = (int*)    malloc(2*ndof*6*sizeof(int));     for (int i=0; i<2*ndof*6; i++)  bc->bc_type[i]       = 1;
    bc->par_dirichlet = (double*) malloc(ndim*ndof*sizeof(double)); for (int i=0; i<ndim*ndof; i++) bc->par_dirichlet[i] = 0.0;
    bc->par_neumann   = (double*) malloc(1*sizeof(double));         for (int i=0; i< 1; i++)        bc->par_neumann[i]   = 0.0;
    bc->neumann       = neumann_default_iga;

    IGAPBCAdd(igap,bc);
}

double* IGAPBCGetParPeriodic(IGAP igap)
{
    return igap->core->bc[0]->par_periodic;
}

// ---- Comm
void IGAPCommInit(IGAP igap)
{
    igap->core->comm = (Comm*)malloc(1*sizeof(Comm));
    igap->core->comm[0] = (Comm)malloc(1*sizeof(struct __Comm));
    // ---- mpi partition
    mpipart_iga(igap->core->mesh[0]);
    // ---- mpi communication
    mpicomm_iga(igap->core);
}

// ---- IOcomm
void IGAPIOcommInit(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int *nbasis    = igap->core->mesh[0]->nbasis;
    int ndof       = igap->core->phys[0]->ndof;
    int torder     = igap->core->phys[0]->torder;

    igap->iocomm = (IOcomm)malloc(1*sizeof(struct __IOcomm));
    int nDof=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *nDof_array;          if (rank==0) { nDof_array      =(int*) malloc(nproc*sizeof(int)); }     else { nDof_array      =(int*) malloc(0); }
    int *iDof_displ_array;    if (rank==0) { iDof_displ_array=(int*) malloc((nproc+1)*sizeof(int)); } else { iDof_displ_array=(int*) malloc(0); }
    assert(MPI_Gather(&nDof,1,MPI_INT,nDof_array,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    if (rank==0) { iDof_displ_array[0]=0; for (int iproc=1; iproc<nproc+1; iproc++) { iDof_displ_array[iproc]=iDof_displ_array[iproc-1]+nDof_array[iproc-1]; } }
    int nDof_global;
    if (rank==0) { nDof_global=iDof_displ_array[nproc]; } assert(MPI_Bcast(&nDof_global,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    int *universalIdx_self=(int*)malloc(nDof*sizeof(int));
    int ia=0, ja=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++) {
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++) {
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++) {
        for (int idof=0; idof<ndof; idof++) { universalIdx_self[ia++]=ndof*(igap->core->mesh[0]->universalIdx_temp[ja])+idof; } ja++;
    }}}
    free(igap->core->mesh[0]->universalIdx_temp);
    if (rank==0) { igap->iocomm->universalIdx=(int*)malloc(nDof_global*sizeof(int)); } else { igap->iocomm->universalIdx=(int*)malloc(0); }
    assert(MPI_Gatherv(universalIdx_self,nDof,MPI_INT,igap->iocomm->universalIdx,nDof_array,iDof_displ_array,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(universalIdx_self);

    igap->iocomm->nDof             = nDof;
    igap->iocomm->nDof_array       = nDof_array;
    igap->iocomm->iDof_displ_array = iDof_displ_array;
    igap->iocomm->nDof_global      = nDof_global;

    igap->core->U_hist=(double*)malloc((torder+1)*nDof*sizeof(double));
}

// ---- Assm
void IGAPAssmInit(IGAP igap)
{
    // ---- Comm
    IGAPCommInit(igap);

    // ---- IOcomm
    IGAPIOcommInit(igap);
    
    // ---- globalIdx (maps local idx to global idx) (required by SNES)
    int *nbasis    = igap->core->mesh[0]->nbasis;
    int nsendPart  = igap->core->mesh[0]->nsendPart;
    int *sendPart  = igap->core->mesh[0]->sendPart;
    int nrecvPart  = igap->core->mesh[0]->nrecvPart;
    int *recvPart  = igap->core->mesh[0]->recvPart;
    int ndof       = igap->core->phys[0]->ndof;
    int nselfIdx   = igap->core->comm[0]->nselfIdx;
    int *selfIdx   = igap->core->comm[0]->selfIdx;
    int nsendIdx   = igap->core->comm[0]->nsendIdx;
    int *sendIdx   = igap->core->comm[0]->sendIdx;
    int *sendPtr   = igap->core->comm[0]->sendPtr;
    int nrecvIdx   = igap->core->comm[0]->nrecvIdx;
    int *recvIdx   = igap->core->comm[0]->recvIdx;
    int *recvPtr   = igap->core->comm[0]->recvPtr;
    int *iDof_displ_array = igap->iocomm->iDof_displ_array;
    int *globalIdx = (int*) malloc((nselfIdx+nrecvIdx)*sizeof(int));
    MPI_Request request_send[nsendPart],request_recv[nrecvPart];
    MPI_Status status;
    int ia=0;
    int ibasis; 
    int iDof_displ;
    assert(MPI_Scatter(iDof_displ_array,1,MPI_INT,&iDof_displ,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    int *globalIdx_self = (int*) malloc(nselfIdx*sizeof(int));
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis[1]*nbasis[2])*ibasis_x+nbasis[2]*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { globalIdx_self[ia++]=iDof_displ+ndof*ibasis+idof; }
    }}}
    for (int i=0; i<nselfIdx; i++) { globalIdx[selfIdx[i]]=globalIdx_self[i]; }
    free(globalIdx_self);
    int *globalIdx_send=(int*) malloc(nsendIdx*sizeof(int));
    int *globalIdx_recv=(int*) malloc(nrecvIdx*sizeof(int));
    for (int i=0; i<nsendIdx; i++) { globalIdx_send[i]=globalIdx[sendIdx[i]]; }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) { assert(MPI_Isend(globalIdx_send+sendPtr[isendPart],sendPtr[isendPart+1]-sendPtr[isendPart],MPI_INT,sendPart[isendPart],1,MPI_COMM_WORLD,request_send+isendPart)==MPI_SUCCESS); }
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) { assert(MPI_Irecv(globalIdx_recv+recvPtr[irecvPart],recvPtr[irecvPart+1]-recvPtr[irecvPart],MPI_INT,recvPart[irecvPart],1,MPI_COMM_WORLD,request_recv+irecvPart)==MPI_SUCCESS); }
    for (int isendPart=0; isendPart<nsendPart; isendPart++) assert(MPI_Wait(request_send+isendPart,&status)==MPI_SUCCESS);
    for (int irecvPart=0; irecvPart<nrecvPart; irecvPart++) assert(MPI_Wait(request_recv+irecvPart,&status)==MPI_SUCCESS);
    for (int i=0; i<nrecvIdx; i++) { globalIdx[recvIdx[i]]=globalIdx_recv[i]; }
    free(globalIdx_send);
    free(globalIdx_recv);
    igap->core->comm[0]->globalIdx = globalIdx;

    // ----
    int nDof        = igap->iocomm->nDof;
    int nDof_global = igap->iocomm->nDof_global;
    igap->assm = (Assm)malloc(1*sizeof(struct __Assm));
    VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&igap->assm->RR);
    int *d_nz=(int*) malloc(nDof*sizeof(int));
    int *o_nz=(int*) malloc(nDof*sizeof(int));
    setnz_iga(d_nz,o_nz,igap->core);
    MatCreateAIJ(MPI_COMM_WORLD,nDof,nDof,nDof_global,nDof_global,0,d_nz,0,o_nz,&igap->assm->TT);
    free(d_nz);
    free(o_nz);
    free(igap->core->mesh[0]->sendBlock);
    free(igap->core->mesh[0]->recvBlock);
    MatSetOption(igap->assm->TT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(igap->assm->TT,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
    MatSetOption(igap->assm->TT,MAT_SUBSET_OFF_PROC_ENTRIES,PETSC_TRUE);
    //MatSetOption(TT,MAT_SYMMETRIC,PETSC_TRUE);
    //MatSetOption(TT,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);
    igap->assm->Residual=Residual_iga;
    igap->assm->Tangent =Tangent_iga;
}

// ---- Solv

void IGAPSolvInit(IGAP igap)
{
    igap->solv = (Solv)malloc(1*sizeof(struct __Solv));

    KSP            ksp;
    PC             pc;
    SNESLineSearch ls;
    // ---- U
    VecDuplicate(igap->assm->RR,&igap->solv->U);
    // ---- SNES
    SNESCreate(MPI_COMM_WORLD,&igap->solv->snes);
    SNESSetFunction(igap->solv->snes,igap->assm->RR,igap->assm->Residual,igap->core);
    SNESSetJacobian(igap->solv->snes,igap->assm->TT,igap->assm->TT,igap->assm->Tangent,igap->core);
    SNESSetTolerances(igap->solv->snes,1.e-10,1.e-20,1.e-20,1,1e9);
    SNESSetType(igap->solv->snes,SNESNEWTONLS);
    SNESGetLineSearch(igap->solv->snes,&ls);
    SNESLineSearchSetType(ls,SNESLINESEARCHBT); //SNESLINESEARCHBT, SNESLINESEARCHL2, SNESLINESEARCHCP, SNESLINESEARCHBASIC, SNESLINESEARCHSHELL
    SNESLineSearchSetOrder(ls,SNES_LINESEARCH_ORDER_CUBIC); //LINEAR, QUADRATIC, CUBIC
    // ---- KSP/PC
    SNESGetKSP(igap->solv->snes,&ksp); 
    KSPGetPC(ksp,&pc);
    KSPSetInitialGuessNonzero(ksp,1);
    KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,1e8,1e8);
/*
    KSPSetType(ksp,KSPMINRES); // The operator and the preconditioner must be symmetric and the preconditioner must be positive definite for this method.
    KSPSetPCSide(ksp,PC_LEFT); //
    PCSetType(pc,PCJACOBI);
    PCJacobiSetUseAbs(pc,1);
*/
    KSPSetType(ksp,KSPGMRES);
    KSPGMRESSetRestart(ksp,10000); // eats memory. 2e4 doesn't work on edison; srun error.
    KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);//KSPGMRESClassicalGramSchmidtOrthogonalization
    KSPSetPCSide(ksp,PC_RIGHT);//PC_RIGHT, PC_LEFT, PC_SYMMETRIC
    PCSetType(pc,PCASM);
    PCASMSetType(pc,PC_ASM_RESTRICT);//PC_ASM_BASIC, PC_ASM_INTERPOLATE, PC_ASM_RESTRICT, PC_ASM_NONE
/*
    KSPSetType(ksp,KSPCG);
    KSPCGSetType(ksp,KSP_CG_SYMMETRIC);
    PCSetType(pc,PCJACOBI);
    PCJacobiSetUseAbs(pc,1);
*/
    // ---- overwrite with runtime option
    SNESSetFromOptions(igap->solv->snes);
}

void IGAPSolvGuess(IGAP igap)
{
    double *U_;
    double *U_hist=igap->core->U_hist;
    VecGetArray(igap->solv->U,&U_);
    for (int i=0; i<igap->iocomm->nDof; i++) { U_[i]=U_hist[i]; }
    VecRestoreArray(igap->solv->U,&U_);
}

void IGAPSolvRand(IGAP igap, int seed, double scale)
{
    srand(seed);
    int nDof=igap->iocomm->nDof;
    int torder = igap->core->phys[0]->torder;
    double *U_hist=igap->core->U_hist;
    for (int i=0; i<(torder+1)*nDof; i++) { U_hist[i]=scale*((double)(rand()%999999999))/999999999.;  }
    IGAPSolvGuess(igap);
}

void IGAPSolvLoadSolution(IGAP igap, const char fname[])
{
    int torder = igap->core->phys[0]->torder;
    double *U_hist=igap->core->U_hist;
    int nDof=igap->iocomm->nDof;
    if ( torder == 0 )
    {
        read_solution(fname,U_hist,igap->iocomm,0,1);
    }
    else
    {
        read_solution(fname,U_hist,igap->iocomm,1,torder+1);
        for (int i=0; i<nDof; i++) { U_hist[i]=U_hist[nDof+i]; }
    }
    IGAPSolvGuess(igap);
}

void IGAPSolvSaveSolution(IGAP igap, const char fname[])
{
    int torder = igap->core->phys[0]->torder;
    double *U_hist=igap->core->U_hist;
    if ( torder == 0 )
    {
        write_solution(fname,U_hist,igap->iocomm,0,1);
    }
    else
    {
        write_solution(fname,U_hist,igap->iocomm,1,torder+1);
    }
}

void IGAPSolvSolve(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Solv solv = igap->solv;
    clock_t time0;

    KSP ksp; SNESGetKSP(solv->snes,&ksp);
    int                 ksp_niter;
    SNESConvergedReason snes_reason;

    // ---- subiterate
    double Fnorm_new,Fnorm_old=999.;
    for (int isubiter=0; isubiter<10000; isubiter++)
    {
        time0=clock();
        SNESSolve(solv->snes,NULL,solv->U);
        KSPGetIterationNumber(ksp,&ksp_niter);
        SNESGetConvergedReason(solv->snes,&snes_reason);
        SNESGetFunctionNorm(solv->snes,&Fnorm_new); if ( double_abs(Fnorm_old-Fnorm_new)/Fnorm_old < 1.e-6 ) { snes_reason=SNES_DIVERGED_LINE_SEARCH; } Fnorm_old=Fnorm_new;
        if (rank==0) { printf("%6d: Fnorm = %e, #lin.ite.= %6d, (%08.2f[min])\n",isubiter+1,Fnorm_new,ksp_niter,((float)(clock()-time0))/CLOCKS_PER_SEC/60.0); }
        if (snes_reason != SNES_DIVERGED_MAX_IT) { break; }
    }
    if (snes_reason<=0 && rank==0) { printf("%s.\n",SNESConvergedReasons[snes_reason]); }
    if (snes_reason>0) { solv->converged=1; } else { solv->converged=0; }
}

int IGAPSolvGetConverged(IGAP igap)
{
    return igap->solv->converged;
}

// ---- Post
double **IGAPPostGetHH(IGAP igap)
{
    return igap->post->HH;
}

void IGAPPostComputeHH(IGAP igap)
{
    iga_intdensity(igap->post->HH[0],igap->core);
}


// ---- 
void IGAPPrintf( const char * format, ... )
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        va_list args;
        va_start (args, format);
        vprintf (format, args);
        va_end (args);
    }
}

// ---- Finalize
void IGAPFinalize(IGAP *igap_ptr)
{
    IGAP igap  = *igap_ptr;

    // ---- Core : Mesh
    free(igap->core->mesh[0]->nelem);
    free(igap->core->mesh[0]->nbasis);
    free(igap->core->mesh[0]->nknot);
    for (int idim=0; idim<3; idim++) { free(igap->core->mesh[0]->knotVector[idim]); }
    free(igap->core->mesh[0]->knotVector);
    free(igap->core->mesh[0]->wVector);
    free(igap->core->mesh[0]->XVector);
    free(igap->core->mesh[0]->boun);
    free(igap->core->mesh[0]->sendPart);
    free(igap->core->mesh[0]->recvPart);
    free(igap->core->mesh[0]);
    free(igap->core->mesh);

    // ---- Core : Phys
    free(igap->core->phys[0]->par_mat); 
    free(igap->core->phys[0]);
    free(igap->core->phys);
 
    // ---- Core : BC
    free(igap->core->bc[0]->bc_type);
    free(igap->core->bc[0]->par_periodic);
    free(igap->core->bc[0]->par_dirichlet);
    free(igap->core->bc[0]->par_neumann);
    free(igap->core->bc[0]);
    free(igap->core->bc);

    // ---- Core : Comm
    free(igap->core->comm[0]->selfIdx);
    free(igap->core->comm[0]->sendIdx);
    free(igap->core->comm[0]->sendPtr);
    free(igap->core->comm[0]->recvIdx);
    free(igap->core->comm[0]->recvPtr);
    free(igap->core->comm[0]->globalIdx);
    free(igap->core->comm[0]);
    free(igap->core->comm);

    // ---- Core : U_hist
    free(igap->core->U_hist);
    
    // ---- Core
    free(igap->core);

    // ---- Assm
    VecDestroy(&igap->assm->RR);
    MatDestroy(&igap->assm->TT);
    free(igap->assm);

    // ---- Solv
    SNESDestroy(&igap->solv->snes);
    VecDestroy(&igap->solv->U);
    free(igap->solv);

    // ---- IOcomm
    free(igap->iocomm->nDof_array);
    free(igap->iocomm->iDof_displ_array);
    free(igap->iocomm->universalIdx);
    free(igap->iocomm);

    // ---- Post
    free(igap->post->HH[0]);
    free(igap->post->HH);
    free(igap->post);

    // ---- Plot
    free(igap->plot);

    // ----
    free(igap);
}
