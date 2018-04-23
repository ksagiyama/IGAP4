#if !defined (INTERFACE)
#define INTERFACE

#include "igap.h"
#include "mesh.h"
#include "phys.h"
#include "bc.h"

// ---- Initialize
void IGAPInitialize(IGAP*);

// ---- Phys
void IGAPPhysAdd(IGAP,Phys);

// ---- Mesh
void IGAPMeshAdd(IGAP,Mesh);

// ---- BC
void IGAPBCAdd(IGAP,BC);
void IGAPBCAddDefault(IGAP);
double* IGAPBCGetParPeriodic(IGAP);

// ---- Comm
void IGAPCommInit(IGAP);

// ---- IOcomm
void IGAPIOcommInit(IGAP);

// ---- Assm
void IGAPAssmInit(IGAP);

// ---- Solv
void IGAPSolvInit(IGAP);
void IGAPSolvGuess(IGAP);
void IGAPSolvRand(IGAP,int,double);
void IGAPSolvSolve(IGAP);
void IGAPSolvLoadSolution(IGAP,const char []);
void IGAPSolvSaveSolution(IGAP,const char []);
int  IGAPSolvGetConverged(IGAP);

// ---- Post
double **IGAPPostGetHH(IGAP);
void IGAPPostComputeHH(IGAP);

// ---- 
void IGAPPrintf( const char * format, ... );

// ---- Finalize
void IGAPFinalize(IGAP*);

#endif
