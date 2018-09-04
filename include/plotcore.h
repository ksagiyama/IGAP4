#if !defined (PLOTCORE)
#define PLOTCORE

#include <mgl2/mgl_cf.h>
#include <mgl2/mpi.h>
#include "plot.h"

void iga_plotfield(HMGL *gr, void *app_, int ifield, int uid[], double uscale, int nppe, const char *sch, const char *opt);
void iga_plotface(HMGL *gr, void *app_, int uid[], double uscale, int nppe, const char *sch, const char *opt);
void iga_plotedge(HMGL *gr, void *app_, int uid[], double uscale, int nppe, const char *sch, const char *opt);

#endif
