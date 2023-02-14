#!/bin/bash

#
# MPI Procs
nproc=4
#
# OMP THREADS
nthreads=1

# Configure problem params
# - initial polynomial order of approximation
p=2

# set enriched order (p+dp)
dp=1
# - manufactured solution
#   0: polynomial solution
#   1: smooth sin*sinh solution, uniform in z
#   2: smooth sin*sinh solution
#   3: arc tan solution
#   4: x boundary layer
isol=4
# - user mode
#   interactive run    (job=0)
#   pre-configured run (job=1)
job=0
# - number of refinements (if job=1)
imax=10

# - MAXNODS
maxnods=525000

export KMP_STACKSIZE=64M

args=" -job ${job} -maxnods ${maxnods} -p ${p}"
args+=" -isol ${isol} -imax ${imax} -nthreads ${nthreads}"

mpirun -np ${nproc} ./pois ${args}
