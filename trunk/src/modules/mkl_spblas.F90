
#if HP3D_USE_INTEL_MKL
include '/opt/intel/mkl/include/mkl_spblas.f90'
#else
module mkl_spblas
end module mkl_spblas
#endif
