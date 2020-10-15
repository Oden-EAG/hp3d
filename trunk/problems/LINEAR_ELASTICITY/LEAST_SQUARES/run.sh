#!/bin/bash

#
# MPI Procs
nproc=1
#
# OMP THREADS
nthreads=1

# Configure problem params
# - initial polynomial order of approximation
p=1
# - manufactured solution

isol=5
# - number of refinements (if job=1)
imax=10

# - MAXNODS
maxnods=100000

# - geometry file
geom='../common/geometries/L_shape'

# - control file
ctrl='../common/control/control_L_shape'

# - boundary condition
bc=8
# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
# mpirun -np ${nproc} ./elast_DPG_pr -file-control ${ctrl} -file-geometry ${geom} -maxnods ${maxnods} -p ${p} -bc ${bc} -exact ${isol} -imax ${imax} -nthreads ${nthreads}
#mpirun -np ${nproc} xterm -hold -e ./elast_BG -job 0 -p ${p} -exact ${isol} -nthreads ${nthreads}
mpirun -np ${nproc} ./elast_FOSLS -file-control ${ctrl} -file-geometry ${geom} -maxnods ${maxnods} -p ${p} -bc ${bc} -exact ${isol} -imax ${imax}
