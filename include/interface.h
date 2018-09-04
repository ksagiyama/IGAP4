#if !defined (INTERFACE)
#define INTERFACE

#include "igap.h"
#include "mesh.h"
#include "phys.h"
#include "bc.h"

// ---- Initialize
void IGAPInitialize(IGAP*,int*,char***);

// ---- Mesh
void IGAPMeshAdd(IGAP,Mesh);
void IGAPMeshRemove(IGAP);
void IGAPMeshKnotRefine(IGAP,int,int,int);
void IGAPMeshOrderElev(IGAP);

// ---- Phys
void IGAPPhysAdd(IGAP,Phys);
void IGAPPhysRemove(IGAP);
void IGAPPhysSetParam(IGAP,int,double);

// ---- BC
void IGAPBCAdd(IGAP,BC);
void IGAPBCRemove(IGAP);
void IGAPBCAddDefault(IGAP);
void IGAPBCSetBCType(IGAP,int,int,int,int);
double* IGAPBCGetParPeriodic(IGAP);

// ---- Comm
void IGAPCommAdd(IGAP);
void IGAPCommRemove(IGAP);

// ---- Soln
void IGAPSolnAdd(IGAP);
void IGAPSolnRemove(IGAP);
void IGAPSolnUpdate(IGAP);
void IGAPSolnMPI(IGAP);
void IGAPSolnZero(IGAP);
void IGAPSolnRand(IGAP,int,double);
void IGAPSolnLoad(IGAP,const char []);
void IGAPSolnSave(IGAP,const char []);
void IGAPSolnScatter(IGAP, double *);
void IGAPSolnGather(IGAP, double *);
void IGAPSolnApplyDirichletBC(IGAP,double*);
void IGAPSolnKnotRefine(IGAP,int,int,int);
void IGAPSolnOrderElev(IGAP);

// ---- Assm
void IGAPAssmAdd(IGAP);
void IGAPAssmRemove(IGAP);

// ---- Solv
void IGAPSolvAdd(IGAP);
void IGAPSolvRemove(IGAP);
void IGAPSolvSolve(IGAP);
int  IGAPSolvGetConverged(IGAP);

// ---- Post
void IGAPPostAdd(IGAP);
void IGAPPostRemove(IGAP);
double **IGAPPostGetHH(IGAP);
void IGAPPostComputeHH(IGAP);

// ---- Plot
void IGAPPlotAdd(IGAP);
void IGAPPlotRemove(IGAP);
void IGAPPlotInitGR(IGAP);
void IGAPPlotBegin(IGAP,const char[]);
void IGAPPlotEnd(IGAP);
void IGAPPlotBeginGIF(IGAP,const char[],int);
void IGAPPlotBeginFrame(IGAP);
void IGAPPlotEndFrame(IGAP);
void IGAPPlotEndGIF(IGAP);
void IGAPPlotTitle(IGAP,const char[],const char[],double);
void IGAPPlotMesh(IGAP,const char[],const char[]);
void IGAPPlotField(IGAP,int,const char[],const char[]);
void IGAPPlotSetFigsize(IGAP,int,int);
void IGAPPlotSetRelplot(IGAP,double,double,double,double);
void IGAPPlotSetRotate(IGAP,double,double,double);
void IGAPPlotSetRanges(IGAP,double,double,double,double,double,double);
void IGAPPlotSetCrange(IGAP,double,double);
void IGAPPlotSetNppe(IGAP,int);
void IGAPPlotSetUid(IGAP,int,int,int);
void IGAPPlotSetUscale(IGAP,double);

// ---- 
void IGAPPrintf( const char * format, ... );

// ---- Finalize
void IGAPFinalize(IGAP*);

#endif
