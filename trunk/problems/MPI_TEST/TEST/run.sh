#!/bin/bash

#
# MPI Ranks
nproc=3
#
# OMP THREADS
nthreads=2

# TACC MPI
# ibrun -n <mpi_procs>

#
# 1)
mpirun -np ${nproc} ./test_mpi -nthreads ${nthreads}
#ibrun -n ${nproc} ./test_mpi -nthreads ${nthreads}
#ibrun -n ${nproc} xterm -hold -e ./test_mpi -nthreads ${nthreads}
