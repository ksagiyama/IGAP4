#if !defined (MESHIMPL)
#define MESHIMPL

#include "mesh.h"

struct __Mesh{

    int porder;
    int *nelem;
    int *nbasis;
    int *nknot;
    double **knotVector;
    double *wVector;
    double *XVector;
    int nboun;
    int *boun;
    int nsendPart;
    int *sendPart;
    int nrecvPart;
    int *recvPart;
    int *sendBlock;              // MESH -> IGA (transient storage)
    int *recvBlock;              // MESH -> IGA (transient storage)
    int *universalIdx_temp;      // MESH -> IGA (transient storage)
    int *nelem_global;           // Global MESH Info (transient storage)
    double **knotVector_global;  // Global MESH Info (transient storage)
    double *wVector_global;      // Global MESH Info (transient storage)
    double *XVector_global;      // Global MESH Info (transient storage)
    int *bc_periodic;            // Global MESH Info (transient storage) -- GEOM Info
};

#endif
