#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include "ioutil.h"

#include "solnimpl.h"

void mpi_printf( const char * format, ... )
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        va_list args;
        va_start (args, format);
        vprintf (format, args);
        va_end (args);
    }
}

int file_exist(char fname[])
{
    int nproc; assert(MPI_Comm_size(MPI_COMM_WORLD,&nproc)==MPI_SUCCESS);
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    int fexist;
    FILE *fptr;
    if (rank == 0)
    {
        fexist = ( ((fptr=fopen(fname,"rb")) != NULL) ? 1:0 );
        if (fexist == 1) { fclose(fptr); }
    }
    assert(MPI_Bcast(&fexist,1,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);

    return fexist;
}

void write_iarray(const char *fname, const char *mode, int iarray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        fwrite(iarray,sizeof(int),size,fptr);
        fclose(fptr);
    }
}

void read_iarray(const char *fname, const char *mode, int iarray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        assert(fread(iarray,sizeof(int),size,fptr)==size);
        fclose(fptr);
    }
    assert(MPI_Bcast(iarray,size,MPI_INT,0,MPI_COMM_WORLD)==MPI_SUCCESS);
}

void write_darray(const char *fname, const char *mode, double darray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);

    if (rank==0)
    {
        FILE *fptr=fopen(fname,mode);
        fwrite(darray,sizeof(double),size,fptr);
        fclose(fptr);
    }
}

void read_darray(const char *fname, const char *mode, double darray[], int size)
{
    int rank;  assert(MPI_Comm_rank(MPI_COMM_WORLD,&rank) ==MPI_SUCCESS);
    if (rank==0)
    { 
        FILE *fptr=fopen(fname,mode);
        assert(fread(darray,sizeof(double),size,fptr)==size);
        fclose(fptr);
    }
    assert(MPI_Bcast(darray,size,MPI_DOUBLE,0,MPI_COMM_WORLD)==MPI_SUCCESS);
}













