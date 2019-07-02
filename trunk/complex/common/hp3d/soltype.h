#ifdef C_MODE
#define SOL_TYPE COMPLEX(DP)
#warning "USING COMPLEX MODE"
#else
#define SOL_TYPE REAL(DP)
#warning "USING REAL MODE"
#endif


