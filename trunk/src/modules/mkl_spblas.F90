
#if HP3D_USE_INTEL_MKL
!include '/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_spblas.f90'
include '/opt/intel/oneapi/mkl/2023.1.0/include/mkl_spblas.f90'
#else
module mkl_spblas
end module mkl_spblas
#endif
