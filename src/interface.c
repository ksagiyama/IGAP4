#include "petscsnes.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "igapimpl.h"
#include "coreimpl.h"
#include "meshimpl.h"
#include "physimpl.h"
#include "bcimpl.h"
#include "commimpl.h"
#include "solnimpl.h"
#include "assmimpl.h"
#include "solvimpl.h"
#include "postimpl.h"

#include "meshcore.h"
#include "commcore.h"
#include "assmcore.h"
#include "postcore.h"

#if defined (__with_IGAPPlot)
#include "plotimpl.h"
#include "plotcore.h"
#endif

#include "mathutil.h"
#include "ioutil.h"

#include "interface.h"


/* ---------------------------------------------------------------------- */
/* -------------------------------- IGAP -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPInitialize(IGAP *igap_ptr, int *argc, char ***argv)
{
    PetscInitialize( argc , argv , (char*)0 , "" );

    *igap_ptr  = NULL;

    IGAP igap         = (IGAP)malloc(1*sizeof(struct __IGAP));
    igap->core        = (Core)malloc(1*sizeof(struct __Core));
    igap->core->mesh  = (Mesh*)malloc(0*sizeof(Mesh));
    igap->core->phys  = (Phys*)malloc(0*sizeof(Phys));
    igap->core->bc    = (BC*)  malloc(0*sizeof(BC));
    igap->core->comm  = (Comm*)malloc(0*sizeof(Comm));
    igap->core->soln  = (Soln*)malloc(0*sizeof(Soln));
    igap->core->nmesh = 0;
    igap->core->nphys = 0;
    igap->core->nbc   = 0;
    igap->core->ncomm = 0;
    igap->core->nsoln = 0;
    igap->assm        = (Assm)malloc(0*sizeof(struct __Assm));
    igap->solv        = (Solv)malloc(0*sizeof(struct __Solv));
    igap->post        = (Post)malloc(0*sizeof(struct __Post));
#if defined (__with_IGAPPlot)
    igap->plot        = (Plot)malloc(0*sizeof(struct __Plot));
#endif
    igap->nassm       = 0;
    igap->nsolv       = 0;
    igap->npost       = 0;
    igap->nplot       = 0;

    *igap_ptr        = igap;
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Mesh -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPMeshAdd(IGAP igap, Mesh mesh)
{
    Mesh* mesh_ptr  = igap->core->mesh;
    igap->core->mesh=(Mesh*)malloc((igap->core->nmesh+1)*sizeof(Mesh));
    for (int i=0; i<igap->core->nmesh; i++) igap->core->mesh[i]=mesh_ptr[i];
    free(mesh_ptr);
    igap->core->mesh[igap->core->nmesh]=mesh;

    int ndim=3;
    mesh->nelem      = (int*)malloc(0);
    mesh->nbasis     = (int*)malloc(0);
    mesh->nknot      = (int*)malloc(0);
    mesh->knotVector = (double**)malloc(ndim*sizeof(double*));
    for (int i=0; i<ndim; i++) mesh->knotVector[i]=(double*)malloc(0);
    mesh->wVector    = (double*)malloc(0);
    mesh->XVector    = (double*)malloc(0);
    mesh->boun       = (int*)malloc(0);
    mesh->sendPart   = (int*)malloc(0);
    mesh->recvPart   = (int*)malloc(0);
    mesh->sendBlock  = (int*)malloc(0);
    mesh->recvBlock  = (int*)malloc(0);
    mesh->universalIdx_temp = (int*)malloc(0);

    igap->core->nmesh++;
}

void IGAPMeshRemove(IGAP igap)
{
    Mesh mesh = igap->core->mesh[igap->core->nmesh-1];
    free(mesh->nelem);
    free(mesh->nbasis);
    free(mesh->nknot);
    for (int idim=0; idim<3; idim++) { free(mesh->knotVector[idim]); }
    free(mesh->knotVector);
    free(mesh->wVector);
    free(mesh->XVector);
    free(mesh->boun);
    free(mesh->sendPart);
    free(mesh->recvPart);
    free(mesh->bc_periodic);
    free(mesh->sendBlock);
    free(mesh->recvBlock);
    free(mesh->nelem_global);
    for (int idim=0; idim<3; idim++) { free(mesh->knotVector_global[idim]); }
    free(mesh->knotVector_global);
    free(mesh->wVector_global);
    free(mesh->XVector_global);
    free(mesh->universalIdx_temp);
    free(mesh);

    igap->core->nmesh--;
}

void IGAPMeshKnotRefine(IGAP igap, int nref_x, int nref_y, int nref_z)
{
    IGAPSolnKnotRefine(igap,nref_x,nref_y,nref_z);
}

void IGAPMeshOrderElev(IGAP igap)
{
    IGAPSolnOrderElev(igap);
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Phys -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPPhysAdd(IGAP igap, Phys phys)
{
    Phys* phys_ptr  = igap->core->phys;
    igap->core->phys=(Phys*)malloc((igap->core->nphys+1)*sizeof(Phys));
    for (int i=0; i<igap->core->nphys; i++) igap->core->phys[i]=phys_ptr[i];
    free(phys_ptr);
    igap->core->phys[igap->core->nphys]=phys;
    igap->core->nphys++;
}

void IGAPPhysRemove(IGAP igap)
{
    Phys phys = igap->core->phys[igap->core->nphys-1];
    free(phys->par_mat); 
    free(phys);

    igap->core->nphys--;
}

void IGAPPhysSetParam(IGAP igap, int ipar, double val)
{
    igap->core->phys[0]->par_mat[ipar]=val;
}

/* ---------------------------------------------------------------------- */
/* --------------------------------  BC  -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPBCAdd(IGAP igap, BC bc)
{
    BC* bc_ptr  = igap->core->bc;
    igap->core->bc=(BC*)malloc((igap->core->nbc+1)*sizeof(BC));
    for (int i=0; i<igap->core->nbc; i++) igap->core->bc[i]=bc_ptr[i];
    free(bc_ptr);
    igap->core->bc[igap->core->nbc]=bc;
    igap->core->nbc++;
}

void IGAPBCRemove(IGAP igap)
{
    BC bc = igap->core->bc[igap->core->nbc-1];
    free(bc->par_periodic);
    free(bc->bc_type);
    free(bc->par_dirichlet);
    free(bc->par_neumann);
    free(bc);

    igap->core->nbc--;
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

void IGAPBCSetBCType(IGAP igap, int iboun_, int idof_, int iorder_, int type_)
{
    BC bc = igap->core->bc[0];
    
    int ndof = igap->core->phys[0]->ndof;

    bc->bc_type[(2*ndof)*iboun_+ndof*iorder_+idof_]=type_;
}

double* IGAPBCGetParPeriodic(IGAP igap)
{
    return igap->core->bc[0]->par_periodic;
}

/* ---------------------------------------------------------------------- */
/* -------------------------------- Comm -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPCommAdd(IGAP igap)
{
    Comm* comm_ptr  = igap->core->comm;
    igap->core->comm=(Comm*)malloc((igap->core->ncomm+1)*sizeof(Comm));
    for (int i=0; i<igap->core->ncomm; i++) igap->core->comm[i]=comm_ptr[i];
    free(comm_ptr);
    igap->core->comm[igap->core->ncomm]=(Comm)malloc(1*sizeof(struct __Comm));

    Comm comm = igap->core->comm[igap->core->ncomm];
    comm->selfIdx   = (int*)malloc(0);
    comm->sendIdx   = (int*)malloc(0);
    comm->sendPtr   = (int*)malloc(0);
    comm->recvIdx   = (int*)malloc(0);
    comm->recvPtr   = (int*)malloc(0);
    comm->globalIdx = (int*)malloc(0);

    igap->core->ncomm++;
}

void IGAPCommRemove(IGAP igap)
{
    Comm comm = igap->core->comm[igap->core->ncomm-1];
    free(comm->selfIdx);
    free(comm->sendIdx);
    free(comm->sendPtr);
    free(comm->recvIdx);
    free(comm->recvPtr);
    free(comm->globalIdx);
    free(comm);

    igap->core->ncomm--;
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Soln -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPSolnAdd(IGAP igap)
{
    Soln* soln_ptr  = igap->core->soln;
    igap->core->soln=(Soln*)malloc((igap->core->nsoln+1)*sizeof(Soln));
    for (int i=0; i<igap->core->nsoln; i++) igap->core->soln[i]=soln_ptr[i];
    free(soln_ptr);
    igap->core->soln[igap->core->nsoln]=(Soln)malloc(1*sizeof(struct __Soln));

    Soln soln = igap->core->soln[igap->core->nsoln];
    soln->U_hist           = (double*)malloc(0);
    soln->nDof_array       = (int*)malloc(0);
    soln->iDof_displ_array = (int*)malloc(0);
    soln->universalIdx     = (int*)malloc(0);

    igap->core->nsoln++;
}

void IGAPSolnRemove(IGAP igap)
{
    Soln soln = igap->core->soln[igap->core->nsoln-1];
    free(soln->U_hist);
    free(soln->nDof_array);
    free(soln->iDof_displ_array);
    free(soln->universalIdx);
    free(soln);

    igap->core->nsoln--;
}

void IGAPSolnMPI(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Soln soln = igap->core->soln[0];
    free(soln->U_hist);
    free(soln->nDof_array);
    free(soln->iDof_displ_array);
    free(soln->universalIdx);

    int *nbasis    = igap->core->mesh[0]->nbasis;
    int ndof       = igap->core->phys[0]->ndof;
    int torder     = igap->core->phys[0]->torder;

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
    //free(igap->core->mesh[0]->universalIdx_temp);
    if (rank==0) { igap->core->soln[0]->universalIdx=(int*)malloc(nDof_global*sizeof(int)); } else { igap->core->soln[0]->universalIdx=(int*)malloc(0); }
    assert(MPI_Gatherv(universalIdx_self,nDof,MPI_INT,igap->core->soln[0]->universalIdx,nDof_array,iDof_displ_array,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    free(universalIdx_self);

    igap->core->soln[0]->nDof             = nDof;
    igap->core->soln[0]->nDof_array       = nDof_array;
    igap->core->soln[0]->iDof_displ_array = iDof_displ_array;
    igap->core->soln[0]->nDof_global      = nDof_global;

    igap->core->soln[0]->U_hist=(double*)malloc((torder+1)*nDof*sizeof(double));
}

void IGAPSolnZero(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if ( igap->core->nsoln == 0 )
    {
        IGAPCommAdd(igap);
        IGAPSolnAdd(igap);
        mpipart_iga(igap->core->mesh[0]);
        mpicomm_iga(igap->core);
        IGAPSolnMPI(igap);
    }

    int nDof=igap->core->soln[0]->nDof;
    int torder = igap->core->phys[0]->torder;
    double *U_hist=igap->core->soln[0]->U_hist;
    for (int i=0; i<(torder+1)*nDof; i++) { U_hist[i]=0.0;  }

    /*if ( torder == 0)
    {
        IGAPSolnApplyDirichletBC(igap,U_hist);
    }
    else
    {
        for (int ihist=0; ihist<torder; ihist++) IGAPSolnApplyDirichletBC(igap,U_hist+nDof*(ihist+1));
        for (int i=0; i<nDof; i++) { U_hist[i]=U_hist[nDof+i]; }
    }*/
}

void IGAPSolnRand(IGAP igap, int seed, double scale)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if ( igap->core->nsoln == 0 )
    {
        IGAPCommAdd(igap);
        IGAPSolnAdd(igap);
        mpipart_iga(igap->core->mesh[0]);
        mpicomm_iga(igap->core);
        IGAPSolnMPI(igap);
    }

    int torder = igap->core->phys[0]->torder;
    int ihist0, ihist1;
    if ( torder == 0 ) { ihist0=0; ihist1=1;        }
    else               { ihist0=1; ihist1=torder+1; }

    Soln soln             = igap->core->soln[0];
    double *U_hist        = soln->U_hist;
    int nDof              = soln->nDof;
    int nDof_global       = soln->nDof_global;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc((ihist1-ihist0)*nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    srand(seed);
    if (rank==0)
    {
        for (int i=0; i<(ihist1-ihist0)*nDof_global; i++) { U_universal[i]=scale*((double)(rand()%999999999))/999999999.;  }
    }
    IGAPSolnScatter(igap,U_universal);
    free(U_universal);

    if ( torder == 0)
    {
        IGAPSolnApplyDirichletBC(igap,U_hist);
    }
    else
    {
        for (int ihist=0; ihist<torder; ihist++) IGAPSolnApplyDirichletBC(igap,U_hist+nDof*(ihist+1));
        for (int i=0; i<nDof; i++) { U_hist[i]=U_hist[nDof+i]; }
    }
}

void IGAPSolnLoad(IGAP igap, const char fname[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if ( igap->core->nsoln == 0 )
    {
        IGAPCommAdd(igap);
        IGAPSolnAdd(igap);
        mpipart_iga(igap->core->mesh[0]);
        mpicomm_iga(igap->core);
        IGAPSolnMPI(igap);
    }

    int torder = igap->core->phys[0]->torder;
    int ihist0, ihist1;
    if ( torder == 0 ) { ihist0=0; ihist1=1;        }
    else               { ihist0=1; ihist1=torder+1; }

    Soln soln             = igap->core->soln[0];
    double *U_hist        = soln->U_hist;
    int nDof              = soln->nDof;
    int nDof_global       = soln->nDof_global;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc((ihist1-ihist0)*nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    if (rank==0)
    {
        FILE *fptr=fopen(fname,"rb");
        assert(fread(U_universal,sizeof(double),(ihist1-ihist0)*nDof_global,fptr)==(ihist1-ihist0)*nDof_global);
        fclose(fptr);
    }
    IGAPSolnScatter(igap,U_universal);
    free(U_universal);

    if ( torder == 0)
    {
        IGAPSolnApplyDirichletBC(igap,U_hist);
    }
    else
    {
        for (int ihist=0; ihist<torder; ihist++) IGAPSolnApplyDirichletBC(igap,U_hist+nDof*(ihist+1));
        for (int i=0; i<nDof; i++) { U_hist[i]=U_hist[nDof+i]; }
    }
}

void IGAPSolnSave(IGAP igap, const char fname[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int torder = igap->core->phys[0]->torder;
    int ihist0, ihist1;
    if ( torder == 0 ) { ihist0=0; ihist1=1;        }
    else               { ihist0=1; ihist1=torder+1; }

    Soln soln             = igap->core->soln[0];
    int nDof_global       = soln->nDof_global;

    double *U_universal     ; if (rank==0) { U_universal     =(double*) malloc((ihist1-ihist0)*nDof_global*sizeof(double)); } else { U_universal     =(double*) malloc(0); }
    IGAPSolnGather(igap,U_universal);
    if (rank==0)
    {
        FILE *fptr=fopen(fname,"wb");
        fwrite(U_universal,sizeof(double),(ihist1-ihist0)*nDof_global,fptr);
        fclose(fptr);
    }
    free(U_universal);
}

void IGAPSolnScatter(IGAP igap, double *U_universal)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int torder = igap->core->phys[0]->torder;
    int ihist0, ihist1;
    if ( torder == 0 ) { ihist0=0; ihist1=1;        }
    else               { ihist0=1; ihist1=torder+1; }

    Soln soln             = igap->core->soln[0];
    double *U_hist        = soln->U_hist;
    int nDof              = soln->nDof;
    int *nDof_array       = soln->nDof_array;
    int *iDof_displ_array = soln->iDof_displ_array;
    int nDof_global       = soln->nDof_global;
    int *universalIdx     = soln->universalIdx;

    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(                nDof_global*sizeof(double)); } else { U_universal_temp=(double*) malloc(0); }
    for (int ihist=ihist0; ihist<ihist1; ihist++)
    {
        if (rank==0)
        {
            for (int i=0; i<nDof_global; i++) { U_universal_temp[i]=U_universal[(ihist-ihist0)*nDof_global+universalIdx[i]]; }
        }
        assert(MPI_Scatterv(U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,U_hist+nDof*ihist,nDof,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    }
    free(U_universal_temp);
}

void IGAPSolnGather(IGAP igap, double *U_universal)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int torder = igap->core->phys[0]->torder;
    int ihist0, ihist1;
    if ( torder == 0 ) { ihist0=0; ihist1=1;        }
    else               { ihist0=1; ihist1=torder+1; }

    Soln soln             = igap->core->soln[0];
    double *U_hist        = soln->U_hist;
    int nDof              = soln->nDof;
    int *nDof_array       = soln->nDof_array;
    int *iDof_displ_array = soln->iDof_displ_array;
    int nDof_global       = soln->nDof_global;
    int *universalIdx     = soln->universalIdx;

    double *U_universal_temp; if (rank==0) { U_universal_temp=(double*) malloc(                nDof_global*sizeof(double)); } else { U_universal_temp=(double*) malloc(0); }
    for (int ihist=ihist0; ihist<ihist1; ihist++)
    {
        assert(MPI_Gatherv(U_hist+nDof*ihist,nDof,MPI_DOUBLE,U_universal_temp,nDof_array,iDof_displ_array,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
        if (rank==0)
        {
            for (int i=0; i<nDof_global; i++) { U_universal[(ihist-ihist0)*nDof_global+universalIdx[i]]=U_universal_temp[i]; }
        }
    }
    free(U_universal_temp);
}

void IGAPSolnUpdate(IGAP igap)
{
    int    torder = igap->core->phys[0]->torder;
    double *U_hist=igap->core->soln[0]->U_hist;
    int    nDof=igap->core->soln[0]->nDof;
    if ( torder == 0 )
    {
    }
    else
    {
        for (int i=nDof*torder-1; i>=0; i--) { U_hist[i+nDof]=U_hist[i]; }
    }
}

void IGAPSolnApplyDirichletBC(IGAP igap, double *U_hist)
{
    Phys phys            = igap->core->phys[0];
    int ndim             = phys->ndim;
    int ndof             = phys->ndof;
    //int torder           = phys->torder;
    
    Mesh mesh            = igap->core->mesh[0];
    int *nbasis          = mesh->nbasis;
    int nboun            = mesh->nboun;
    int *boun            = mesh->boun;

    BC bc                = igap->core->bc[0];
    int *bc_type         = bc->bc_type;
    //double *par_dirichlet= bc->par_dirichlet;

    //Soln soln            = igap->core->soln[0];
    /* ================ assemble U ================ */
    //double *U_hist = soln->U_hist;

    int ax_,ay_,az_;
    int nbasis_[ndim],ibasis_z_0,ibasis_z_1;
    for (int iboun=0; iboun<nboun; iboun++)
    {
        switch (boun[iboun]/2)
        {
        case 0:
            nbasis_[0]=nbasis[1]; ax_=nbasis[2];
            nbasis_[1]=nbasis[2]; ay_=1;
            nbasis_[2]=nbasis[0]; az_=nbasis[1]*nbasis[2];
            break;
        case 1:
            nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nbasis_[1]=nbasis[2]; ay_=1; 
            nbasis_[2]=nbasis[1]; az_=nbasis[2];
            break;
        case 2:
            nbasis_[0]=nbasis[0]; ax_=nbasis[1]*nbasis[2];
            nbasis_[1]=nbasis[1]; ay_=nbasis[2];
            nbasis_[2]=nbasis[2]; az_=1;
            break;
        default:
            exit(EXIT_FAILURE);
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
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    int ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_0;
                    U_hist[ndof*ibasis+idof]=0; // homogeneous only for now
                }}
        }}
        // high-order
        for (int idof=0; idof<ndof; idof++) {
            if (bc_type[(2*ndof)*boun[iboun]+1*ndof+idof]==0) {
                for (int ibasis_x_=0; ibasis_x_<nbasis_[0]; ibasis_x_++) {
                for (int ibasis_y_=0; ibasis_y_<nbasis_[1]; ibasis_y_++) {
                    int ibasis=ax_*ibasis_x_+ay_*ibasis_y_+az_*ibasis_z_1;
                    U_hist[ndof*ibasis+idof]=0; // homogeneous only, for simple geometry.
                }}
        }}
    } // for iboun
}

void IGAPSolnKnotRefine(IGAP igap, int nref_x, int nref_y, int nref_z)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Mesh mesh = igap->core->mesh[0];
    int porder       = mesh->porder;
    int ndim=3;
    int *nelem_global          = mesh->nelem_global;
    int nknot[ndim]; for (int i=0; i<ndim; i++) nknot[i]=nelem_global[i]+1+2*porder;
    double **knotVector = mesh->knotVector_global;
    double *wVector     = mesh->wVector_global;
    double *XVector     = mesh->XVector_global;
    double *UVector;
    int    nknot_[ndim];
    double *knotVector_;
    double *wVector_;
    double *XVector_;
    double *UVector_;

    // uniform refinement
    int nref[]={nref_x,nref_y,nref_z};
    int r[ndim], r0[ndim], r1[ndim];
    for (int idim=0; idim<ndim; idim++) { r[idim]=0; r0[idim]=0; r1[idim]=0; }
    double **X=(double**)malloc(ndim*sizeof(double*));
    if (rank==0)
    {
        for (int idim=0; idim<ndim; idim++)
        {
            for (int i=0; i<nknot[idim]-1; i++)
            {
                if (knotVector[idim][i]+1.e-12<knotVector[idim][i+1])
                {
                    r[idim]+=(int)(pow(2.0,nref[idim]))-1;
                }
            }
            X[idim]=(double*)malloc(r[idim]*sizeof(double));
            r[idim]=0;
            for (int i=0; i<nknot[idim]-1; i++)
            {
                if (knotVector[idim][i]+1.e-12<knotVector[idim][i+1])
                {
                    int np=(int)(pow(2.0,nref[idim]))-1;
                    for (int j=0; j<np; j++)
                    {
                        X[idim][r[idim]+j]=knotVector[idim][i]*(np-j)/(np+1)+knotVector[idim][i+1]*(j+1)/(np+1);
                    }
                    r[idim]+=np;
                }
            }
            if (r[idim]>0)
            {
                while ( X[idim][r0[idim]]             < knotVector[idim][porder] )                 { r0[idim]++; } // only significant if periodic
                while ( X[idim][(r[idim]-1)-r1[idim]] > knotVector[idim][(nknot[idim]-1)-porder] ) { r1[idim]++; } // only significant if periodic
            }
        }
    }
    else
    {
        for (int idim=0; idim<ndim; idim++) X[idim]=(double*)malloc(0);
    }

    int initialized = (int)(igap->core->nsoln>0);
    if (initialized)
    {
        Phys phys  = igap->core->phys[0];
        int torder = phys->torder;
        int ihist0, ihist1;
        if ( torder == 0 ) { ihist0=0; ihist1=1;        }
        else               { ihist0=1; ihist1=torder+1; }
        Soln soln       = igap->core->soln[0];
        int nDof = soln->nDof_global;
        if (rank==0)
        {
            UVector = (double*)malloc((ihist1-ihist0)*nDof*sizeof(double));
        }
        else
        {
            UVector = (double*)malloc(0);
        }
        IGAPSolnGather(igap,UVector);
    }
    else
    {
        UVector = (double*)malloc(0);
    }

    if (rank==0)
    {
        int *bc_periodic = mesh->bc_periodic;

        // nknot
        for (int idim=0; idim<ndim; idim++) nknot_[idim]=nknot[idim];
        for (int idim=0; idim<ndim; idim++) nknot[idim] =nknot[idim]+2*porder*bc_periodic[idim];
        // knotVector
        for (int idim=0; idim<ndim; idim++)
        {
            knotVector_      = knotVector[idim];
            knotVector[idim] = (double*)malloc(nknot[idim]*sizeof(double));
            int l=0;
            for (int i=0; i<porder*bc_periodic[idim]; i++)  knotVector[idim][l++] = knotVector_[0]-1;
            for (int i=0; i<nknot_[idim]; i++)              knotVector[idim][l++] = knotVector_[i];
            for (int i=0; i<porder*bc_periodic[idim]; i++)  knotVector[idim][l++] = knotVector_[nknot_[idim]-1]+1;
            free(knotVector_);
        }
        int nbasis_[ndim]; for (int idim=0; idim<ndim; idim++) nbasis_[idim] = nknot_[idim]-porder-1;
        int nbasis[ndim] ; for (int idim=0; idim<ndim; idim++) nbasis[idim]  = nknot[idim]-porder-1;
        // wVector
        wVector_ = wVector;
        wVector  = (double*)malloc(nbasis[0]*nbasis[1]*nbasis[2]*sizeof(double));
        for (int i=0; i<nbasis[0]*nbasis[1]*nbasis[2]; i++) wVector[i]=0;
        for (int i=0; i<nbasis_[0]; i++) {
        for (int j=0; j<nbasis_[1]; j++) {
        for (int k=0; k<nbasis_[2]; k++) {
            wVector[nbasis[1]*nbasis[2]*(i+porder*bc_periodic[0])+nbasis[2]*(j+porder*bc_periodic[1])+(k+porder*bc_periodic[2])]=wVector_[nbasis_[1]*nbasis_[2]*i+nbasis_[2]*j+k];
        }}}
        free(wVector_);
        // XVector
        XVector_= XVector;
        XVector = (double*)malloc(ndim*nbasis[0]*nbasis[1]*nbasis[2]*sizeof(double));
        for (int i=0; i<ndim*nbasis[0]*nbasis[1]*nbasis[2]; i++) XVector[i]=0;
        for (int i=0; i<nbasis_[0]; i++) {
        for (int j=0; j<nbasis_[1]; j++) {
        for (int k=0; k<nbasis_[2]; k++) {
        int ibasis_ = nbasis_[1]*nbasis_[2]*i+nbasis_[2]*j+k;
        int ibasis  = nbasis[1]*nbasis[2]*(i+porder*bc_periodic[0])+nbasis[2]*(j+porder*bc_periodic[1])+(k+porder*bc_periodic[2]);
        for (int idim=0; idim<ndim; idim++) {
            XVector[ndim*ibasis+idim] = XVector_[ndim*ibasis_+idim];
        }}}}
        free(XVector_);
        // UVector
        if (initialized)
        {
        Phys phys  = igap->core->phys[0];
        int ndof   = phys->ndof;
        int torder = phys->torder;
        int ihist0, ihist1;
        if ( torder == 0 ) { ihist0=0; ihist1=1;        }
        else               { ihist0=1; ihist1=torder+1; }
        int nDof_  = ndof*(nbasis_[0]-porder*bc_periodic[0])*(nbasis_[1]-porder*bc_periodic[1])*(nbasis_[2]-porder*bc_periodic[2]);
        int nDof   = ndof*nbasis[0]*nbasis[1]*nbasis[2];
        UVector_ = UVector;
        UVector  = (double*)malloc((ihist1-ihist0)*nDof*sizeof(double));
        for (int i=0; i<(ihist1-ihist0)*nDof; i++) UVector[i]=0;
        for (int ihist=ihist0; ihist<ihist1; ihist++) {
        for (int i=0; i<nbasis_[0]; i++) {
        for (int j=0; j<nbasis_[1]; j++) {
        for (int k=0; k<nbasis_[2]; k++) {
        int ibasis_ = (nbasis_[1]-porder*bc_periodic[1])*(nbasis_[2]-porder*bc_periodic[2])*(i%(nbasis_[0]-porder*bc_periodic[0]))+(nbasis_[2]-porder*bc_periodic[2])*(j%(nbasis_[1]-porder*bc_periodic[1]))+(k%(nbasis_[2]-porder*bc_periodic[2]));
        int ibasis  = nbasis[1]*nbasis[2]*(i+porder*bc_periodic[0])+nbasis[2]*(j+porder*bc_periodic[1])+(k+porder*bc_periodic[2]);
        for (int idof=0; idof<ndof; idof++) {
            UVector[(ihist-ihist0)*nDof+ndof*ibasis+idof] = UVector_[(ihist-ihist0)*nDof_+ndof*ibasis_+idof];
        }}}}}
        free(UVector_);
        }

        for (int idim=0; idim<ndim; idim++)
        {
            if (r[idim]==0) continue;
            //
            double *T;
            int *colIdx;
            int *rowPtr;
            nknot_[idim]     = nknot[idim];
            nknot[idim]+=r[idim];
            knotVector_      = knotVector[idim];
            knotVector[idim] = (double*)malloc(nknot[idim]*sizeof(double));
            knotinsertmap(knotVector[idim],&T,&colIdx,&rowPtr,knotVector_,X[idim],nknot_[idim]-porder-1,r[idim],porder);
            free(knotVector_);

            int ai,aj,ak,ni,nj,nk,ni_,aj_,ak_;
            switch (idim)
            {
                case 0: ni=nknot[0]-porder-1; nj=nknot[1]-porder-1; nk=nknot[2]-porder-1; ai=nj*nk; aj=nk   ; ak=1 ; ni_=nknot_[0]-porder-1; aj_=nk    ; ak_=1  ; break;
                case 1: ni=nknot[1]-porder-1; nj=nknot[0]-porder-1; nk=nknot[2]-porder-1; ai=nk   ; aj=ni*nk; ak=1 ; ni_=nknot_[1]-porder-1; aj_=ni_*nk; ak_=1  ; break;
                case 2: ni=nknot[2]-porder-1; nj=nknot[0]-porder-1; nk=nknot[1]-porder-1; ai=1    ; aj=ni*nk; ak=ni; ni_=nknot_[2]-porder-1; aj_=ni_*nk; ak_=ni_; break;
            }
            double *Pw         = (double*)malloc(ni_*sizeof(double));
            double *Qw         = (double*)malloc(ni *sizeof(double));
            // wVector
            wVector_ = wVector;
            wVector  = (double*)malloc(ni*nj*nk*sizeof(double));
            for (int j=0; j<nj; j++)
            {
                for (int k=0; k<nk; k++)
                {
                    for (int i=0; i<ni_; i++) { Pw[i]=wVector_[ai*i+aj_*j+ak_*k]; }
                    for (int i=0; i<ni ; i++) { Qw[i]=0.0; }
                    for (int i=0; i<ni ; i++) { for (int l=rowPtr[i]; l<rowPtr[i+1]; l++) { Qw[i]+=T[l]*Pw[colIdx[l]]; } }
                    for (int i=0; i<ni ; i++) { wVector[ai*i+aj*j+ak*k]=Qw[i]; }
                }
            }
            // XVector
            XVector_ = XVector;
            XVector  = (double*)malloc(ndim*ni*nj*nk*sizeof(double));
            for (int jdim=0; jdim<ndim; jdim++)
            {
            for (int j=0; j<nj; j++)
            {
                for (int k=0; k<nk; k++)
                {
                    for (int i=0; i<ni_; i++) { Pw[i]=XVector_[ndim*(ai*i+aj_*j+ak_*k)+jdim]*wVector_[ai*i+aj_*j+ak_*k]; }
                    for (int i=0; i<ni ; i++) { Qw[i]=0.0; }
                    for (int i=0; i<ni ; i++) { for (int l=rowPtr[i]; l<rowPtr[i+1]; l++) { Qw[i]+=T[l]*Pw[colIdx[l]]; } }
                    for (int i=0; i<ni ; i++) { XVector[ndim*(ai*i+aj*j+ak*k)+jdim]=Qw[i]/wVector[ai*i+aj*j+ak*k]; }
                }
            }
            }
            // U_hist
            if (initialized)
            {
            Phys phys  = igap->core->phys[0];
            int ndof   = phys->ndof;
            int torder = phys->torder;
            int ihist0, ihist1;
            if ( torder == 0 ) { ihist0=0; ihist1=1;        }
            else               { ihist0=1; ihist1=torder+1; }
            int nDof_ = ndof*ni_*nj*nk;
            int nDof  = ndof*ni *nj*nk;
            UVector_ = UVector;
            UVector  = (double*)malloc((ihist1-ihist0)*nDof*sizeof(double));
            for (int ihist=ihist0; ihist<ihist1; ihist++)
            {
            for (int jdof=0; jdof<ndof; jdof++)
            {
            for (int j=0; j<nj; j++)
            {
                for (int k=0; k<nk; k++)
                {
                    for (int i=0; i<ni_; i++) { Pw[i]=UVector_[(ihist-ihist0)*nDof_+ndof*(ai*i+aj_*j+ak_*k)+jdof]*wVector_[ai*i+aj_*j+ak_*k]; }
                    for (int i=0; i<ni ; i++) { Qw[i]=0.0; }
                    for (int i=0; i<ni ; i++) { for (int l=rowPtr[i]; l<rowPtr[i+1]; l++) { Qw[i]+=T[l]*Pw[colIdx[l]]; } }
                    for (int i=0; i<ni ; i++) { UVector[(ihist-ihist0)*nDof+ndof*(ai*i+aj*j+ak*k)+jdof]=Qw[i]/wVector[ai*i+aj*j+ak*k]; }
                }
            }
            }
            }
            free(UVector_);
            }
            free(XVector_);
            free(wVector_);
            free(Pw);
            free(Qw);
            free(T);
            free(colIdx);
            free(rowPtr);
        }

        // nknot
        for (int idim=0; idim<ndim; idim++) nknot_[idim]=nknot[idim];
        for (int idim=0; idim<ndim; idim++) nknot[idim]-=(2*porder+r0[idim]+r1[idim])*bc_periodic[idim];
        // knotVector_global
        for (int idim=0; idim<ndim; idim++)
        {
            knotVector_      = knotVector[idim];
            knotVector[idim] = (double*)malloc(nknot[idim]*sizeof(double));
            for (int i=0; i<nknot[idim]; i++) knotVector[idim][i] = knotVector_[(porder+r0[idim])*bc_periodic[idim]+i];
            free(knotVector_);
        }

        for (int idim=0; idim<ndim; idim++) nbasis_[idim] = nknot_[idim]-porder-1;
        for (int idim=0; idim<ndim; idim++) nbasis[idim]  = nknot[idim]-porder-1;
        // wVector_global
        wVector_= wVector;
        wVector = (double*)malloc(nbasis[0]*nbasis[1]*nbasis[2]*sizeof(double));
        for (int i=0; i<nbasis[0]; i++) {
        for (int j=0; j<nbasis[1]; j++) {
        for (int k=0; k<nbasis[2]; k++) {
            wVector[nbasis[1]*nbasis[2]*i+nbasis[2]*j+k]=wVector_[nbasis_[1]*nbasis_[2]*(i+(porder+r0[0])*bc_periodic[0])+nbasis_[2]*(j+(porder+r0[1])*bc_periodic[1])+(k+(porder+r0[2])*bc_periodic[2])];
        }}}
        free(wVector_);
        // XVector_global
        XVector_= XVector;
        XVector = (double*)malloc(ndim*nbasis[0]*nbasis[1]*nbasis[2]*sizeof(double));
        for (int i=0; i<nbasis[0]; i++) {
        for (int j=0; j<nbasis[1]; j++) {
        for (int k=0; k<nbasis[2]; k++) {
        int ibasis_ = nbasis_[1]*nbasis_[2]*(i+(porder+r0[0])*bc_periodic[0])+nbasis_[2]*(j+(porder+r0[1])*bc_periodic[1])+(k+(porder+r0[2])*bc_periodic[2]);
        int ibasis  = nbasis[1]*nbasis[2]*i+nbasis[2]*j+k;
        for (int idim=0; idim<ndim; idim++) {
            XVector[ndim*ibasis+idim] = XVector_[ndim*ibasis_+idim];
        }}}}
        free(XVector_);
        // UVector_global
        if (initialized)
        {
        Phys phys  = igap->core->phys[0];
        int ndof   = phys->ndof;
        int torder = phys->torder;
        int ihist0, ihist1;
        if ( torder == 0 ) { ihist0=0; ihist1=1;        }
        else               { ihist0=1; ihist1=torder+1; }
        int nDof_  = ndof*nbasis_[0]*nbasis_[1]*nbasis_[2];
        int nDof   = ndof*(nbasis[0]-porder*bc_periodic[0])*(nbasis[1]-porder*bc_periodic[1])*(nbasis[2]-porder*bc_periodic[2]);
        UVector_= UVector;
        UVector = (double*)malloc((ihist1-ihist0)*nDof*sizeof(double));
        for (int ihist=ihist0; ihist<ihist1; ihist++) {
        for (int i=0; i<nbasis[0]-porder*bc_periodic[0]; i++) {
        for (int j=0; j<nbasis[1]-porder*bc_periodic[1]; j++) {
        for (int k=0; k<nbasis[2]-porder*bc_periodic[2]; k++) {
        int ibasis_ = nbasis_[1]*nbasis_[2]*(i+(porder+r0[0])*bc_periodic[0])+nbasis_[2]*(j+(porder+r0[1])*bc_periodic[1])+(k+(porder+r0[2])*bc_periodic[2]);
        int ibasis  = (nbasis[1]-porder*bc_periodic[1])*(nbasis[2]-porder*bc_periodic[2])*i+(nbasis[2]-porder*bc_periodic[2])*j+k;
        for (int idof=0; idof<ndof; idof++) {
            UVector[(ihist-ihist0)*nDof+ndof*ibasis+idof] = UVector_[(ihist-ihist0)*nDof_+ndof*ibasis_+idof];
        }}}}}
        free(UVector_);
        }
    }

    for (int idim=0; idim<ndim; idim++) nelem_global[idim] = nknot[idim]-2*porder-1;
    assert(MPI_Bcast(nelem_global,ndim,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
    mesh->wVector_global = wVector;
    mesh->XVector_global = XVector;
    if (initialized)
    {
        mpipart_iga(mesh);
        mpicomm_iga(igap->core);
        IGAPSolnMPI(igap);

        IGAPSolnScatter(igap,UVector);
        Phys phys      = igap->core->phys[0];
        int torder     = phys->torder;
        Soln soln      = igap->core->soln[0];
        double *U_hist = soln->U_hist;
        int nDof       = igap->core->soln[0]->nDof;
        if ( torder != 0)
        {
            for (int i=0; i<nDof; i++) { U_hist[i]=U_hist[nDof+i]; }
        }
    }
    free(UVector);

    for (int idim=0; idim<ndim; idim++) free(X[idim]);
    free(X);
}

void IGAPSolnOrderElev(IGAP igap)
{}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Assm -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPAssmAdd(IGAP igap)
{
    free(igap->assm);
    igap->assm = (Assm)malloc(1*sizeof(struct __Assm));

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
    int *iDof_displ_array = igap->core->soln[0]->iDof_displ_array;
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
    free(igap->core->comm[0]->globalIdx);
    igap->core->comm[0]->globalIdx = globalIdx;

    // ----
    int nDof        = igap->core->soln[0]->nDof;
    int nDof_global = igap->core->soln[0]->nDof_global;
    VecCreateMPI(MPI_COMM_WORLD,nDof,nDof_global,&igap->assm->RR);
    int *d_nz=(int*) malloc(nDof*sizeof(int));
    int *o_nz=(int*) malloc(nDof*sizeof(int));
    setnz_iga(d_nz,o_nz,igap->core);
    MatCreateAIJ(MPI_COMM_WORLD,nDof,nDof,nDof_global,nDof_global,0,d_nz,0,o_nz,&igap->assm->TT);
    free(d_nz);
    free(o_nz);
    //free(igap->core->mesh[0]->sendBlock);
    //free(igap->core->mesh[0]->recvBlock);
    MatSetOption(igap->assm->TT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(igap->assm->TT,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
    MatSetOption(igap->assm->TT,MAT_SUBSET_OFF_PROC_ENTRIES,PETSC_TRUE);
    //MatSetOption(TT,MAT_SYMMETRIC,PETSC_TRUE);
    //MatSetOption(TT,MAT_SYMMETRY_ETERNAL,PETSC_TRUE);
    igap->assm->Residual=Residual_iga;
    igap->assm->Tangent =Tangent_iga;

    igap->nassm++;
}

void IGAPAssmRemove(IGAP igap)
{
    VecDestroy(&igap->assm->RR);
    MatDestroy(&igap->assm->TT);

    igap->nassm--;
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Solv -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPSolvAdd(IGAP igap)
{
    free(igap->solv);
    igap->solv = (Solv)malloc(1*sizeof(struct __Solv));

    if (igap->nassm==0) IGAPAssmAdd(igap);

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
    KSPGMRESSetRestart(ksp,1000); //10000,20000, eat memory.
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

    igap->nsolv++;
}

void IGAPSolvRemove(IGAP igap)
{
    SNESDestroy(&igap->solv->snes);
    VecDestroy(&igap->solv->U);

    igap->nsolv--;
}

void IGAPSolvSolve(IGAP igap)
{
    if ( igap->core->nsoln ==0) exit(EXIT_FAILURE);

    if (igap->nsolv==0) IGAPSolvAdd(igap);

    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Solv solv = igap->solv;
    clock_t time0;

    KSP ksp; SNESGetKSP(solv->snes,&ksp);
    int                 ksp_niter;
    SNESConvergedReason snes_reason;

    // ---- set initial guess
    double *U_;
    double *U_hist=igap->core->soln[0]->U_hist;
    VecGetArray(igap->solv->U,&U_);
    for (int i=0; i<igap->core->soln[0]->nDof; i++) { U_[i]=U_hist[i]; }
    VecRestoreArray(igap->solv->U,&U_);

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
    if (snes_reason<=0) { IGAPPrintf("%s.\n",SNESConvergedReasons[snes_reason]); }
    if (snes_reason>0) { solv->converged=1; } else { solv->converged=0; }
    IGAPSolnUpdate(igap);
}

int IGAPSolvGetConverged(IGAP igap)
{
    if (igap->nsolv==0) IGAPSolvAdd(igap);

    return igap->solv->converged;
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Post -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPPostAdd(IGAP igap)
{
    free(igap->post);
    igap->post        = (Post)malloc(1*sizeof(struct __Post));

    igap->post->HH=(double**)malloc(1*sizeof(double*));
    igap->post->HH[0]=(double*)malloc(igap->core->phys[0]->ndensity*sizeof(double));

    igap->npost++;
}

void IGAPPostRemove(IGAP igap)
{
    free(igap->post->HH[0]);
    free(igap->post->HH);

    igap->npost--;
}

double **IGAPPostGetHH(IGAP igap)
{
    if (igap->npost==0) IGAPPostAdd(igap);

    return igap->post->HH;
}

void IGAPPostComputeHH(IGAP igap)
{
    if ( igap->core->nsoln ==0) exit(EXIT_FAILURE);

    if (igap->npost==0) IGAPPostAdd(igap);

    iga_intdensity(igap->post->HH[0],igap->core);
}


/* ---------------------------------------------------------------------- */
/* -------------------------------- Plot -------------------------------- */
/* ---------------------------------------------------------------------- */

#if !defined (__with_IGAPPlot)
void IGAPPlotAdd(IGAP igap)                                                                         { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotRemove(IGAP igap)                                                                      { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotInitGR(IGAP igap)                                                                      { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotBegin(IGAP igap, const char fname[])                                                   { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotEnd(IGAP igap)                                                                         { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotBeginGIF(IGAP igap, const char fname[], int ms)                                        { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotBeginFrame(IGAP igap)                                                                  { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotEndFrame(IGAP igap)                                                                    { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotEndGIF(IGAP igap)                                                                      { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotMesh(IGAP igap, const char sch[], const char opt[])                                    { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotField(IGAP igap, int ifield, const char sch[], const char opt[])                       { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetFigsize(IGAP igap, int a0, int a1)                                                  { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetRelplot(IGAP igap, double a0, double a1, double a2, double a3)                      { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetRotate(IGAP igap, double a0, double a1, double a2)                                  { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetRanges(IGAP igap, double a0, double a1, double a2, double a3, double a4, double a5) { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetCrange(IGAP igap, double a0, double a1)                                             { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetNppe(IGAP igap, int a0)                                                             { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetUid(IGAP igap, int a0, int a1, int a2)                                              { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
void IGAPPlotSetUscale(IGAP igap, double a0)                                                        { IGAPPrintf("\n\nIGAP is not installed with plot.\n\n"); }
#else
void IGAPPlotAdd(IGAP igap)
{
    free(igap->plot);
    igap->plot        = (Plot)malloc(1*sizeof(struct __Plot));

    Plot plot        = igap->plot;
    int ndim         = igap->core->phys[0]->ndim;
    plot->gr         = (HMGL*)malloc(1*sizeof(HMGL));
    plot->gr[0]      = mgl_create_graph(1200,1200);
    plot->figsize    = (int*)malloc(2*sizeof(int));
    plot->figsize[0] = 1200;
    plot->figsize[1] = 1200;
    plot->nppe       = 4;
    plot->uid        = (int*)malloc(ndim*sizeof(int));
    plot->uid[0]     = 0;
    plot->uid[1]     = 10;
    plot->uid[2]     = 20;
    plot->uscale     = 1.0;
    plot->relplot    = (double*)malloc(4*sizeof(double));
    plot->relplot[0] = -0.5;
    plot->relplot[1] =  1.5;
    plot->relplot[2] = -0.5;
    plot->relplot[3] =  1.5;
    plot->rotate     = (double*)malloc(3*sizeof(double));
    plot->rotate[0]  = 60.;
    plot->rotate[1]  = 30.;
    plot->rotate[2]  =  0.;
    plot->ranges     = (double*)malloc(2*ndim*sizeof(double));
    //plot->ranges[0]  = -0.14;    plot->ranges[1]   =  1.14;
    //plot->ranges[2]  = -0.14;    plot->ranges[3]   =  1.14;
    //plot->ranges[4]  = -0.14;    plot->ranges[5]   =  1.14;
    plot->ranges[0]  = -3;    plot->ranges[1]   =  3;
    plot->ranges[2]  = -3;    plot->ranges[3]   =  3;
    plot->ranges[4]  = -3;    plot->ranges[5]   =  3;
    plot->crange     = (double*)malloc(2*sizeof(double));
    plot->crange[0]  =     0;    plot->crange[1]   =     3;

    igap->nplot++;
}

void IGAPPlotRemove(IGAP igap)
{
    Plot plot        = igap->plot;
    mgl_delete_graph(plot->gr[0]);
    free(plot->gr);
    free(plot->figsize);
    free(plot->uid);
    free(plot->relplot);
    free(plot->rotate);
    free(plot->ranges);
    free(plot->crange);
    
    igap->nplot--;
}

void IGAPPlotInitGR(IGAP igap)
{
    Plot plot = igap->plot;

    double *relplot = plot->relplot;
    double *rotate  = plot->rotate;
    double *ranges  = plot->ranges;
    double *crange  = plot->crange;
    mgl_set_origin(plot->gr[0],0,0,0);
    mgl_relplot(plot->gr[0],relplot[0],relplot[1],relplot[2],relplot[3]); // Warning: relplot changes fontsize, too. -1.4 would results in a diff. size
    mgl_rotate(plot->gr[0],rotate[0],rotate[1],rotate[2]);
    mgl_aspect(plot->gr[0],1,1,1);
    mgl_set_ranges(plot->gr[0],ranges[0],ranges[1],ranges[2],ranges[3],ranges[4],ranges[5]);
    //mgl_set_ticks(plot->gr,'x',0.1,4,0.0);
    //mgl_set_ticks(plot->gr,'y',0.1,4,0.0);
    //mgl_set_ticks(plot->gr,'z',0.1,4,0.0);
    //mgl_axis(plot->gr,"xyzU","","");
    mgl_set_range_val(plot->gr[0],'c',crange[0],crange[1]);
    mgl_set_alpha(plot->gr[0],0);
}

void IGAPPlotBegin(IGAP igap, const char fname[])
{
    if (igap->nplot==0) IGAPPlotAdd(igap);
    Plot plot = igap->plot;
    int    *figsize = plot->figsize;
    mgl_set_size(plot->gr[0],figsize[0],figsize[1]);
    plot->fname = fname;
}

void IGAPPlotEnd(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;

    if (rank!=0)
    {
        mgl_mpi_send(gr[0],0); 
        //mgl_clf(gr); 
    }
    else
    {
        HMGL gr_temp;
        int *figsize = plot->figsize;
        gr_temp = mgl_create_graph(figsize[0],figsize[1]);
        for (int iproc=1; iproc<nproc; iproc++)
        {
            mgl_mpi_recv(gr_temp,iproc);
            mgl_combine_gr(gr[0],gr_temp);
        }
        mgl_delete_graph(gr_temp);
    }
    if (rank==0) { mgl_write_frame(gr[0],plot->fname,""); }

    mgl_clf(gr[0]);
}

void IGAPPlotBeginGIF(IGAP igap, const char fname[], int ms)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);
    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;
    int    *figsize = plot->figsize;
    mgl_set_size(plot->gr[0],figsize[0],figsize[1]);

    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0) mgl_start_gif(gr[0],fname,ms);
}

void IGAPPlotBeginFrame(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;

    if (rank==0) mgl_new_frame(gr[0]);
}

void IGAPPlotEndFrame(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;

    if (rank!=0)
    {
        mgl_mpi_send(gr[0],0); 
    }
    else
    {
        HMGL gr_temp;
        int *figsize = plot->figsize;
        gr_temp = mgl_create_graph(figsize[0],figsize[1]);
        for (int iproc=1; iproc<nproc; iproc++)
        {
            mgl_mpi_recv(gr_temp,iproc);
            mgl_combine_gr(gr[0],gr_temp);
        }
        mgl_delete_graph(gr_temp);
    }

    if (rank==0) mgl_end_frame(gr[0]);
    else mgl_clf(gr[0]); 
}

void IGAPPlotEndGIF(IGAP igap)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;

    if (rank==0) mgl_close_gif(gr[0]);
}

void IGAPPlotTitle(IGAP igap, const char txt[], const char stl[], double size)
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0) mgl_title(igap->plot->gr[0],txt,stl,size);
}

void IGAPPlotMesh(IGAP igap, const char sch[], const char opt[])
{
    if ( igap->core->nsoln ==0) IGAPSolnZero(igap);

    IGAPPlotInitGR(igap);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;
    int    nppe     = plot->nppe;
    int    *uid     = plot->uid;
    double uscale   = plot->uscale*0;

    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    char randcolor[10];
    srand(999);
    sprintf(randcolor,"{x%06X}",((unsigned int)(rand()*(rank+1)))%((unsigned int)int_pow(16,6)));
    iga_plotface(gr,igap->core,uid,uscale,nppe,randcolor,"");
    iga_plotedge(gr,igap->core,uid,uscale,nppe,"k","");
}

void IGAPPlotField(IGAP igap, int ifield, const char sch[], const char opt[])
{
    if ( igap->core->nsoln ==0) exit(EXIT_FAILURE);

    IGAPPlotInitGR(igap);

    Plot plot       = igap->plot;
    HMGL *gr        = plot->gr;
    int    nppe     = plot->nppe;
    int    *uid     = plot->uid;
    double uscale   = plot->uscale;
    iga_plotfield(gr,igap->core,ifield,uid,uscale,nppe,sch,opt);
}

void IGAPPlotSetFigsize(IGAP igap, int figsize0_, int figsize1_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->figsize[0] = figsize0_;
    igap->plot->figsize[1] = figsize1_;
}

void IGAPPlotSetRelplot(IGAP igap, double relplot0_, double relplot1_, double relplot2_, double relplot3_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->relplot[0] = relplot0_;
    igap->plot->relplot[1] = relplot1_;
    igap->plot->relplot[2] = relplot2_;
    igap->plot->relplot[3] = relplot3_;
}

void IGAPPlotSetRotate(IGAP igap, double rotate0_, double rotate1_, double rotate2_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->rotate[0] = rotate0_;
    igap->plot->rotate[1] = rotate1_;
    igap->plot->rotate[2] = rotate2_;
}

void IGAPPlotSetRanges(IGAP igap, double ranges0_, double ranges1_, double ranges2_, double ranges3_, double ranges4_, double ranges5_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->ranges[0] = ranges0_;
    igap->plot->ranges[1] = ranges1_;
    igap->plot->ranges[2] = ranges2_;
    igap->plot->ranges[3] = ranges3_; 
    igap->plot->ranges[4] = ranges4_;
    igap->plot->ranges[5] = ranges5_;
}

void IGAPPlotSetCrange(IGAP igap, double crange0_, double crange1_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->crange[0] = crange0_;
    igap->plot->crange[1] = crange1_;
}

void IGAPPlotSetNppe(IGAP igap, int nppe_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->nppe = nppe_;
}

void IGAPPlotSetUid(IGAP igap, int uid0_, int uid1_, int uid2_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->uid[0] = uid0_;
    igap->plot->uid[1] = uid1_;
    igap->plot->uid[2] = uid2_;
}

void IGAPPlotSetUscale(IGAP igap, double uscale_)
{
    if (igap->nplot==0) IGAPPlotAdd(igap);

    igap->plot->uscale = uscale_;
}

#endif


/* ---------------------------------------------------------------------- */
/* -------------------------------- Util -------------------------------- */
/* ---------------------------------------------------------------------- */

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


/* ---------------------------------------------------------------------- */
/* -------------------------------- IGAP -------------------------------- */
/* ---------------------------------------------------------------------- */

void IGAPFinalize(IGAP *igap_ptr)
{
    IGAP igap  = *igap_ptr;

    // ---- Core : Mesh
    while (igap->core->nmesh > 0) IGAPMeshRemove(igap);
    free(igap->core->mesh);

    // ---- Core : Phys
    while (igap->core->nphys > 0) IGAPPhysRemove(igap);
    free(igap->core->phys);
 
    // ---- Core : BC
    while (igap->core->nbc > 0) IGAPBCRemove(igap);
    free(igap->core->bc);

    // ---- Core : Comm
    while (igap->core->ncomm > 0) IGAPCommRemove(igap);
    free(igap->core->comm);

    // ---- Core : Soln
    while (igap->core->nsoln > 0) IGAPSolnRemove(igap);
    free(igap->core->soln);

    // ---- Core
    free(igap->core);

    // ---- Assm
    while (igap->nassm > 0) IGAPAssmRemove(igap);
    free(igap->assm);

    // ---- Solv
    while (igap->nsolv > 0) IGAPSolvRemove(igap);
    free(igap->solv);

    // ---- Post
    while (igap->npost > 0) IGAPPostRemove(igap);
    free(igap->post);

    // ---- Plot
#if defined (__with_IGAPPlot)
    while (igap->nplot > 0) IGAPPlotRemove(igap);
    free(igap->plot);
#endif
    // ----
    free(igap);

    PetscFinalize();
}
