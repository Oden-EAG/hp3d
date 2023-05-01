#!/bin/bash

#
# MPI Procs
nproc=6
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
job=1
# - number of refinements (if job=1)
imax=12

# - MAXNODS
maxnods=525000

export KMP_STACKSIZE=64M

args=" -job ${job} -maxnods ${maxnods} -p ${p}"
args+=" -isol ${isol} -imax ${imax} -nthreads ${nthreads}"

mpirun -np ${nproc} ./pois ${args}

# to run fichera corner
#   set NEXACT = 0 in control
#   set geomtype = 1 in set_initial_mesh.F90
#   set csn = 1 in ficheracornerdirichlet.F90 (always double check this as one may forget)
#   set Fval = 0 in getF.F90 (always double check this as may forget)
#   if u want to run exact solution with fichera corner then comment Nexact if statement in exact error to compute exact solution