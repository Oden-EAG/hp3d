#!/bin/bash

#
# MPI Ranks
nproc=2
#
# OMP THREADS
nthreads=1

# TACC MPI
# ibrun -n <mpi_procs>

#
# 1)
ibrun -n ${nproc} ./test -nthreads ${nthreads}

