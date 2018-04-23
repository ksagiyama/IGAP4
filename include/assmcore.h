#if !defined (ASSMCORE)
#define ASSMCORE

#include "petscsnes.h"
#include "coreimpl.h"
#include "assmimpl.h"

PetscErrorCode Residual_iga(SNES snes, Vec U_, Vec RR, void *app_);
PetscErrorCode Tangent_iga(SNES snes, Vec U_, Mat TT, Mat Pmat_, void *app_);
void dirichlet_iga(double *dirichlet, double par[], int face_id, int bc_order, int idof, double *knotVector_[], int nknot_[], int nbasis_[], int ndim, int porder);
void dirichletinit_iga(double *U_, void *app_);
void neumann_default_iga(double *neumann_, double par_neumann[], int face_id, int bc_order, int idof, double xq[]);
void setnz_iga(int *d_nz, int *o_nz, Core core);

#endif
