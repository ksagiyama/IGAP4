static char help[] = "gradelasttime on unitcube.\n";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "petscsnes.h"

#include "igap4.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    // ---- Initialize

    PetscInitialize( &argc , &argv , (char*)0 , help );
    IGAP igap;
    IGAPInitialize( &igap );

    // ---- Mesh : add generated beam mesh

    IGAPMeshAdd ( igap , MeshIGABeam(1.0,1.0,1.0,3,2,1,1,1) );

    // ---- Phys : add transient gradient elasticity application

    IGAPPhysAdd ( igap , PhysGradElastTime() );

    // ---- BC   : add default boundary conditions

    IGAPBCAddDefault ( igap );

    // ---- Assm : allocate memory, etc...

    IGAPAssmInit ( igap );

    // ---- Solv : set random I.C. and solve one step

    IGAPSolvInit ( igap );
    IGAPSolvRand ( igap , 123456 , 0.01 );

    IGAPSolvSolve ( igap );
    if ( IGAPSolvGetConverged ( igap ) )
    {
        // ---- Post : compute energy and total volume(=1)

        IGAPPostComputeHH ( igap );
        double **HH = IGAPPostGetHH ( igap );
        IGAPPrintf("PI = %e.\n\n",HH[0][0]);
        IGAPPrintf("Volume = %e.\n\n",HH[0][1]);
    }

    // ---- Finalize

    IGAPFinalize( &igap );
    PetscFinalize();

    return 0;
}
