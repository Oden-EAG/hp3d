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

isol=0
# - number of refinements (if job=1)
imax=10

# - MAXNODS
maxnods=525000

# - geometry file
geom='../common/geometries/hexa'

# - control file
ctrl='../common/control/control'

# - materials file
matr='../common/materials/my_material'

# - boundary condition
bc=1
# TACC MPI
# ibrun -n <mpi_procs>
#
# 1)
mpirun -np ${nproc} ./hyperelast_BG -file-control ${ctrl} -file-geometry ${geom} -p ${p} -bc ${bc} -exact ${isol} -file-materials ${matr} -imax ${imax} -job 0 -maxnods ${maxnods} 
#mpirun -np ${nproc} xterm -hold -e ./elast_BG -job 0 -p ${p} -exact ${isol} -nthreads ${nthreads}
