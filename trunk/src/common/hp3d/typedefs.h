
#if HP3D_COMPLEX
#define VTYPE complex(8)
#else
#define VTYPE real(8)
#endif

#if HP3D_COMPLEX
#define MPI_VTYPE MPI_COMPLEX16
#else
#define MPI_VTYPE MPI_REAL8
#endif
