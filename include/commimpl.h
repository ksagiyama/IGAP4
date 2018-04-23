#if !defined (COMMIMPL)
#define COMMIMPL

#include "comm.h"

struct __Comm{

    int nselfIdx;
    int *selfIdx;
    int nsendIdx;
    int *sendIdx;
    int *sendPtr;
    int nrecvIdx;
    int *recvIdx;
    int *recvPtr;
    int *globalIdx;
};

#endif
