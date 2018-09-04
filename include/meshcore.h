#if !defined (MESHCORE)
#define MESHCORE

#include "mesh.h"

void mpipart_iga(Mesh mesh);
void eval_Bspline(double N[], double *kV[], int p, double xi[]);
double evalN(double *kV, int nk, int k, int i, int p, double xi);
void knotinsert(double *Ubar, double *Qw, double *U, double *Pw, double *X, int n, int r, int p);
void knotinsertmap(double *Ubar, double **T_, int **colIdx_, int **rowPtr_, double *U, double *X, int n, int r, int p);

#endif
