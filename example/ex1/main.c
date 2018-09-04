#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "igap4.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    // ---- Initialize

    IGAP igap;
    IGAPInitialize( &igap , &argc , &argv );

    // ---- Mesh : add generated beam mesh

    IGAPMeshAdd ( igap , MeshIGABeam (8.0,1.0,1.0,0,2,0,0,0) );
    IGAPSolnKnotRefine ( igap , 5 , 2 , 2 );

    // ---- Phys : add transient gradient elasticity application

    IGAPPhysAdd ( igap , PhysGradElastTime() );
    IGAPPhysSetParam ( igap , 6 , 1. );
    IGAPPhysSetParam ( igap , 7 , 0.005 );

    // ---- BC   : add default boundary conditions

    IGAPBCAddDefault ( igap );
    IGAPBCSetBCType ( igap , 0 , 0 , 0 , 0);
    IGAPBCSetBCType ( igap , 0 , 1 , 0 , 0);
    IGAPBCSetBCType ( igap , 0 , 2 , 0 , 0);
    IGAPBCSetBCType ( igap , 1 , 0 , 0 , 0);
    IGAPBCSetBCType ( igap , 1 , 1 , 0 , 0);
    IGAPBCSetBCType ( igap , 1 , 2 , 0 , 0);

    // ---- Soln

    IGAPSolnRand ( igap , 123456 , 0.01 );

    // ---- Solv : set random I.C. and solve one step

    for (int i=0; i<3; i++) 
    {
        IGAPSolvSolve ( igap );
        if ( IGAPSolvGetConverged ( igap ) )
        {
            // ---- Post : compute energy and total volume(=1)
            char fname[512];
            sprintf(fname,"data_%d.bin",i);
            IGAPSolnSave ( igap , fname );

            IGAPPostComputeHH ( igap );
            double **HH = IGAPPostGetHH ( igap );
            IGAPPrintf("PI = %e.\n\n",HH[0][0]);
            IGAPPrintf("Volume = %e.\n\n",HH[0][1]);
        }
    }

    // ---- Plot

    IGAPPlotBeginGIF ( igap , "field.gif" , 200 );
    IGAPPlotSetNppe ( igap , 2 );
    IGAPPlotSetRanges ( igap, -0.14 , 8.14 , -4.14 , 4.14 , -4.14 , 4.14 );
    IGAPPlotSetRotate( igap , 60 , 30 , 0 );
    for (int i=0; i<3; i++)
    {
        char fname[512];
        sprintf(fname,"data_%d.bin",i);
        IGAPSolnLoad ( igap , fname );
        IGAPPlotBeginFrame ( igap );
        IGAPPlotField ( igap , 5 , "{x202020}{xFFA500}{xFFA500}{x3CB371}{x3CB371}{xD2691E}|#", "meshnum 7");
        IGAPPlotEndFrame ( igap );
    }
    IGAPPlotEndGIF ( igap );

    // ---- Finalize

    IGAPFinalize( &igap );

    return 0;
}
