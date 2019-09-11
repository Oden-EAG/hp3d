
Currently Loaded Modules:
  1) git/2.9.0       7) libfabric/1.7.0  13) qt5/5.11.2
  2) autotools/1.1   8) hdf5/1.10.4      14) ospray/1.8.0
  3) xalt/2.6.5      9) impi/18.0.2      15) paraview/5.6.0
  4) TACC           10) python2/2.7.15   16) tau/2.27.2
  5) cmake/3.10.2   11) petsc/3.11
  6) intel/18.0.2   12) swr/18.3.3

 

Fatal error in PMPI_Reduce_scatter: Other MPI error, error stack:
PMPI_Reduce_scatter(1662).....: MPI_Reduce_scatter(sbuf=0x2b639b5fd060, rbuf=0x2b62d81fd518, rcnts=0xcff06c0, MPI_INTEGER, MPI_SUM, comm=0x84000001) failed
MPIR_Reduce_scatter_impl(1199): fail failed
MPIR_Reduce_scatter(1156).....: fail failed
MPIR_Reduce_scatter_intra(357): Out of memory
