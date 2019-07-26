#!/bin/bash
#----------------------------------------------------
# Example Slurm job script
# for TACC Stampede2 SKX nodes
#
#   *** Hybrid Job on SKX Normal Queue ***
# 
#       This sample script specifies:
#         10 nodes (capital N)
#         40 total MPI tasks (lower case n); this is 4 tasks/node
#         12 OpenMP threads per MPI task (48 threads per node)
#
# Last revised: 20 Oct 2017
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch skx.slurm" on Stampede2 login node.
#
#   -- Check current queue with: "squeue -u sh43394"
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- In most cases it's best to keep
#      ( MPI ranks per node ) x ( threads per rank )
#      to a number no more than 48 (total cores).
#
#   -- If you're running out of memory, try running
#      fewer tasks and/or threads per node to give each 
#      process access to more memory.
#
#   -- IMPI and MVAPICH2 both do sensible process pinning by default.
#
#----------------------------------------------------

#SBATCH -J hp3d            # Job name
#SBATCH -o hp3d.o%j        # Name of stdout output file
#SBATCH -e hp3d.e%j        # Name of stderr error file
#SBATCH -p skx-normal      # Queue (partition) name
#SBATCH -N 16              # Total # of nodes
#SBATCH -n 16              # Total # of mpi tasks
#SBATCH -t 00:59:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@ices.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-ti2           # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

# export PROFILEDIR=/work/05246/sh43394/stampede2/par_hp3d/trunk/problems/MPI_POISSON/GALERKIN/tau/profiles
# export TRACEDIR=/work/05246/sh43394/stampede2/par_hp3d/trunk/problems/MPI_POISSON/GALERKIN/tau/traces
# export TAU_PROFILE=1
# export TAU_TRACE=1

module list
pwd
date

# Set thread count (default value is 1)...

# export OMP_NUM_THREADS=24

# Launch MPI code...

ibrun ./pois -job 1 -p 5 -isol 0 -imax 15 -nthreads 48

date
# ---------------------------------------------------
