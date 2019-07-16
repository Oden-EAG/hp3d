#!/bin/bash

#
# MPI Ranks
nproc=3
#
# OMP THREADS
nthreads=1

# Configure problem params
# - initial polynomial order of approximation
p=3
# - manufactured solution
#   0: polynomial solution
#   1: smooth sin*sinh solution, uniform in z
#   2: smooth sin*sinh solution
isol=0

# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
#ibrun -n ${nproc} ./test -p ${p} -isol ${isol} -nthreads ${nthreads}
ibrun -n ${nproc} xterm -hold -e ./test -p ${p} -isol ${isol} -nthreads ${nthreads}
