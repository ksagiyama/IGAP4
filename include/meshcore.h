#if !defined (MESHCORE)
#define MESHCORE

#include "mesh.h"

void mpipart_iga(Mesh mesh);
void nurbs_knotinsert(double *U_universal, const int nelem[], int nref, int porder, int ndim, int ndof, int bc_periodic[]);
void eval_Bspline(double N[], double *kV[], int p, double xi[]);
double evalN(double *kV, int nk, int k, int i, int p, double xi);

#endif
