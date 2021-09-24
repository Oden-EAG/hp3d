#!/bin/bash

#
# MPI Procs
nproc=1
#
# OMP THREADS
nthreads=1

# Configure problem params
# - initial polynomial order of approximation
p=4
# - manufactured solution
#   0: polynomial solution
#   1: smooth sin*sinh solution, uniform in z
#   2: smooth sin*sinh solution
isol=0
# - number of refinements (if job=1)
imax=2

# - MAXNODS
maxnods=525000

file='./geometries/tet_hexa_prism_pyram'


export KMP_STACKSIZE=64M # p=5

# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
mpirun -np ${nproc} ./tests -maxnods ${maxnods} -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads} -file-geometry ${file}
#ibrun -n ${nproc} ./tests -job 1 -maxnods ${maxnods} -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads}
#ibrun -n ${nproc} xterm -hold -e ./tests -job 0 -p ${p} -isol ${isol} -nthreads ${nthreads}
