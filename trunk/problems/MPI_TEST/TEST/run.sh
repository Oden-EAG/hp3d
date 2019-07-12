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
#ibrun -n ${nproc} ./test -nthreads ${nthreads}
ibrun -n ${nproc} xterm -hold -e ./test -nthreads ${nthreads}
