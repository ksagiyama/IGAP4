#if !defined (PLOTIMPL)
#define PLOTIMPL

#include "plot.h"

struct __Plot{

    int dummy;
/*
    {
        
        int uid[]        = {0,10,20};// index for displacement
        double uscale    = 1.;
        int nppe         = 4;
        int figsize[]    = {1200,1200};
        double relplot[] = {-0.5,1.5,-0.5,1.5};
        double rotate[]  = {60,30,0};
        double ranges[]  = {0.5-0.64,0.5+0.64,0.5-0.64,0.5+0.64,0.5-0.64,0.5+0.64};
        double crange[]  = {0,3};
        for (int i=-12;i<-11;i++)
        {
            par_periodic[2]=0.01*i;
            sprintf(fname1,"%d_%d_999_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f_%08.8f",6,nurbs[0].porder,par_mat[0],par_mat[1],par_mat[2],par_mat[3],par_mat[4],par_periodic[0],par_periodic[1],par_periodic[2],par_periodic[3],par_periodic[4],par_periodic[5],par_periodic[6],par_periodic[7],par_periodic[8]);
            sprintf(fname,"data/U_%s.bin",fname1); read_solution(fname,U_hist,igap->iocomm,0,1);
            sprintf(fname,"plot/phase_%s.png",fname1);
            iga_plotfield(fname,igap->core,uid,uscale,nppe,figsize,relplot,rotate,ranges,crange,"{x202020}{xFFA500}{xFFA500}{x3CB371}{x3CB371}{xD2691E}|#","meshnum 15");
        }
    }
*/
};

#endif
