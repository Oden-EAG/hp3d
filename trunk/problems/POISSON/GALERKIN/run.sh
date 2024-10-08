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
# - user mode
#   interactive run    (job=0)
#   pre-configured run (job=1)
job=0
# - number of refinements (if job=1)
imax=2

# Initial mesh
geom='./geometries/hexa_orient0'

# MAXNODS
maxnods=525000

# OMP stack size
export KMP_STACKSIZE=64M

args=" -job ${job} -maxnods ${maxnods} -p ${p}"
args+=" -file-geometry ${geom}"
args+=" -isol ${isol} -imax ${imax} -nthreads ${nthreads}"

mpirun -np ${nproc} ./pois ${args}
#ibrun -n ${nproc} ./pois ${args}
