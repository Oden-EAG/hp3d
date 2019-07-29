#!/bin/bash

#
# MPI Procs
nproc=2
#
# OMP THREADS
nthreads=24

# Configure problem params
# - initial polynomial order of approximation
p=5
# - manufactured solution
#   0: polynomial solution
#   1: smooth sin*sinh solution, uniform in z
#   2: smooth sin*sinh solution
isol=0
# - number of refinements (if job=1)
imax=6

# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
ibrun -n ${nproc} ./pois -job 1 -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads}
#ibrun -n ${nproc} xterm -hold -e ./pois -job 0 -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads}
