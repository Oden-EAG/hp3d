! fortran 90 format to replace syscom.blk
!
#if C_MODE
#define VTYPE complex(8)
#else
#define VTYPE real(8)
#endif

#if C_MODE
#define MPI_VTYPE MPI_COMPLEX16
#else
#define MPI_VTYPE MPI_REAL8
#endif
