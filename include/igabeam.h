#if !defined (IGABEAM)
#define IGABEAM

void igabeam(double Lx, double Ly, double Lz, int mref, int porder, int bc_periodic_x, int bc_periodic_y, int bc_periodic_z, int **bc_periodic_, int **nelem_global_, double ***knotVector_global_, double **wVector_global_, double **XVector_global_);

#endif
