#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "physimpl.h"
#include "physgradelasttime.h"
#include "gradelasttime.h"

Phys PhysGradElastTime()
{
    Phys phys=(Phys)malloc(1*sizeof(struct __Phys));
    
    phys->ndim=3;
    phys->ndof=3;
    phys->torder=2;
    phys->sorder=2;
    phys->ndensity=2;
    phys->nfield=6;

    //                    r ,  curv,   e0 , e345 , le,   rho,  cii,    dt,    t
    double par_mat[]={ 0.25,  0.50,  500.,  250.,  0.025, 1.0,  1.0,  1.e-3, 0.0 };
    phys->par_mat=(double*)malloc(9*sizeof(double));
    for (int i=0; i<9; i++) phys->par_mat[i]=par_mat[i];
    
    phys->residual       = gradelasttime_residual;
    phys->tangent        = gradelasttime_tangent;
    phys->field          = gradelasttime_field;
    phys->density        = gradelasttime_density;
    phys->densityb       = gradelasttime_densityb;
    phys->assert_par_mat = gradelasttime_assert_par_mat;
    
    phys->assert_par_mat(par_mat);

    return phys;
}
