#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mathutil.h"

#include "coreimpl.h"
#include "commcore.h"

void mpicomm_iga(Core core)
{
    Comm comm = core->comm[0];
    free(comm->selfIdx);
    free(comm->sendIdx);
    free(comm->sendPtr);
    free(comm->recvIdx);
    free(comm->recvPtr);

    // ----
    int porder     = core->mesh[0]->porder;
    int *nbasis    = core->mesh[0]->nbasis;
    int *nelem     = core->mesh[0]->nelem;
    int nsendPart  = core->mesh[0]->nsendPart;
    int nrecvPart  = core->mesh[0]->nrecvPart;
    int *sendBlock = core->mesh[0]->sendBlock;
    int *recvBlock = core->mesh[0]->recvBlock;
    int ndof       = core->phys[0]->ndof;

    // ----
    int nbasis_y_active=nelem[1]+porder;
    int nbasis_z_active=nelem[2]+porder;
    
    int ibasis,ia;
    // ---- self
    int nselfIdx=ndof*(nbasis[0]*nbasis[1]*nbasis[2]);
    int *selfIdx=(int*)malloc(nselfIdx*sizeof(int));
    ia=0;
    for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
    for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
    for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
        ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
        for (int idof=0; idof<ndof; idof++) { selfIdx[ia++]=ndof*ibasis+idof; }
    }}}

    // ---- send
    int nsendIdx=0;
    int *sendPtr=(int*)malloc((nsendPart+1)*sizeof(int)); sendPtr[0]=nsendIdx;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: nsendIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nsendIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nsendIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nsendIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nsendIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nsendIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nsendIdx+=ndof*(porder*porder*porder); break;
        }
        sendPtr[i+1]=nsendIdx;
    }
    int *sendIdx=(int*)malloc(nsendIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nsendPart; i++)
    {
        switch (sendBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=0; ibasis_x<porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { sendIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }

    // ---- recv
    int nrecvIdx=0;
    int *recvPtr=(int*)malloc((nrecvPart+1)*sizeof(int)); recvPtr[0]=nrecvIdx;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: nrecvIdx+=ndof*(nbasis[0]*nbasis[1]*porder); break;
        case 1: nrecvIdx+=ndof*(nbasis[0]*nbasis[2]*porder); break;
        case 2: nrecvIdx+=ndof*(nbasis[1]*nbasis[2]*porder); break;
        case 3: nrecvIdx+=ndof*(nbasis[0]*porder*porder); break;
        case 4: nrecvIdx+=ndof*(nbasis[1]*porder*porder); break;
        case 5: nrecvIdx+=ndof*(nbasis[2]*porder*porder); break;
        case 6: nrecvIdx+=ndof*(porder*porder*porder); break;
        }
        recvPtr[i+1]=nrecvIdx;
    }
    int *recvIdx=(int*)malloc(nrecvIdx*sizeof(int));
    ia=0;
    for (int i=0; i<nrecvPart; i++)
    {
        switch (recvBlock[i]) {
        case 0: // (face0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 1: // (face1)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 2: // (face2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 3: // (edge0)
            for (int ibasis_x=0; ibasis_x<nbasis[0]; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 4: // (edge1)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=0; ibasis_y<nbasis[1]; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 5: // (edge2)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=0; ibasis_z<nbasis[2]; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        case 6: // (corner)
            for (int ibasis_x=nbasis[0]; ibasis_x<nbasis[0]+porder  ; ibasis_x++){
            for (int ibasis_y=nbasis[1]; ibasis_y<nbasis[1]+porder  ; ibasis_y++){
            for (int ibasis_z=nbasis[2]; ibasis_z<nbasis[2]+porder  ; ibasis_z++){
                ibasis=(nbasis_y_active*nbasis_z_active)*ibasis_x+nbasis_z_active*ibasis_y+ibasis_z;
                for (int idof=0; idof<ndof; idof++) { recvIdx[ia++]=ndof*ibasis+idof; }
            }}}
            break;
        }
    }

    // ----
    comm->nselfIdx = nselfIdx;
    comm->selfIdx  = selfIdx;
    comm->nsendIdx = nsendIdx;
    comm->sendIdx  = sendIdx;
    comm->sendPtr  = sendPtr;
    comm->nrecvIdx = nrecvIdx;
    comm->recvIdx  = recvIdx;
    comm->recvPtr  = recvPtr;
}
