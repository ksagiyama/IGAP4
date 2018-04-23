#if !defined (GRADELAST)
#define GRADELAST

void gradelast_residual(double *residual, double *u, double *par);
void gradelast_tangent(double *tangent, double *u, double *par);
void gradelast_density(double h[], double u[], double *par);
void gradelast_densityb(double h[], double u[], double traction, double *par, int order);
void gradelast_field(double *field, double *u, double *par);
void gradelast_assert_par_mat(double *par);

#endif
