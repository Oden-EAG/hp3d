#!/bin/bash

#
# MPI Procs
nproc=1
#
# OMP THREADS
nthreads=1

# Configure problem params
# - initial polynomial order of approximation
p=2
# - manufactured solution
#   0: polynomial solution
#   1: smooth sin*sinh solution, uniform in z
#   2: smooth sin*sinh solution
isol=0
# Interactive mode   (job=0)
# Pre-configured run (job=1)
job=0
# - number of refinements (if job=1)
imax=3

# - MAXNODS
maxnods=525000

export KMP_STACKSIZE=32M

# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
mpirun -np ${nproc} ./maxw -job ${job} -maxnods ${maxnods} -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads}
#ibrun -n ${nproc} ./maxw -job ${job} -maxnods ${maxnods} -p ${p} -isol ${isol} -imax ${imax} -nthreads ${nthreads}
#ibrun -n ${nproc} xterm -hold -e ./maxw -job ${job} -p ${p} -isol ${isol} -nthreads ${nthreads}
