#if !defined (IGAPIMPL__)
#define IGAPIMPL__

#include "igap.h"

#include "core.h"
#include "assm.h"
#include "solv.h"
#include "post.h"
#include "plot.h"

struct __IGAP{

Core core;
Assm assm;
Solv solv;
Post post;
Plot plot;
int  nassm;
int  nsolv;
int  npost;
int  nplot;
};

#endif
