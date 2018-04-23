#if !defined (PLOTCORE)
#define PLOTCORE

#include "plot.h"

void iga_plotfield(const char *fname, void *app_, int uid[], double uscale, int nppe, int figsize[], double relplot[],  double rotate[], double ranges[], double crange[], const char *sch, const char *opt);

#endif
