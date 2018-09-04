#if !defined (IOUTIL)
#define IOUTIL

void mpi_printf( const char * format, ... );
int file_exist(char fname[]);
void write_iarray(const char *fname, const char *mode, int iarray[], int size);
void read_iarray(const char *fname, const char *mode, int iarray[], int size);
void write_darray(const char *fname, const char *mode, double darray[], int size);
void read_darray(const char *fname, const char *mode, double darray[], int size);

#endif
