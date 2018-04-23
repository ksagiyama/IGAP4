#if !defined (GRADELASTTIME)
#define GRADELASTTIME

void gradelasttime_residual(double *residual, double *u, double *par);
void gradelasttime_tangent(double *tangent, double *u, double *par);
void gradelasttime_density(double h[], double u[], double *par);
void gradelasttime_densityb(double h[], double u[], double traction, double *par, int order);
void gradelasttime_field(double *field, double *u, double *par);
void gradelasttime_assert_par_mat(double *par);

#endif
