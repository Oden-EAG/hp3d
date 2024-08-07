#!/bin/bash

#
# MPI Procs
nproc=1
#
# OMP THREADS
nthreads=1

# Configure problem parameters:
# - initial polynomial order of approximation
p=1
# - manufactured solution
#   0: polynomial solution
isol=0
# - user mode
#   interactive run    (job=0)
#   pre-configured run (job=1)
job=0
# - number of refinements (if job=1)
imax=3

# MAXNODS
maxnods=525000

# OMP Stacksize
export KMP_STACKSIZE=32M

args=" -job ${job} -maxnods ${maxnods} -p ${p}"
args+=" -isol ${isol} -imax ${imax} -nthreads ${nthreads}"

mpirun -np ${nproc} ./maxw ${args}
#ibrun -n ${nproc} ./maxw ${args}
