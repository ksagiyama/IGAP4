#if !defined (BCIMPL)
#define BCIMPL

#include "bc.h"

typedef void (*neumann_ptr)(double *neumann_, double par_neumann[], int face_id, int bc_order, int idof, double xq[]);

struct __BC{

    int         *bc_type;
    double      *par_periodic;
    double      *par_dirichlet;
    double      *par_neumann;
    neumann_ptr neumann;
};

#endif
