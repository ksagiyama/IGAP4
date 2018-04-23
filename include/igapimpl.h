#if !defined (IGAPIMPL__)
#define IGAPIMPL__

#include "igap.h"

#include "coreimpl.h"
#include "assmimpl.h"
#include "solvimpl.h"
#include "iocommimpl.h"
#include "postimpl.h"
#include "plotimpl.h"

struct __IGAP{

Core core;
Assm assm;
Solv solv;
IOcomm iocomm;
Post post;
Plot plot;

};

#endif
