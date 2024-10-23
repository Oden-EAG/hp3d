
#if HP3D_USE_INTEL_MKL
include '/opt/intel/oneapi/mkl/2023.1.0/include/mkl_spblas.f90'
#else
module mkl_spblas
end module mkl_spblas
#endif
