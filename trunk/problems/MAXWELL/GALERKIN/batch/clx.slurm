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
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 00:02:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DMSXXXXXX       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

export PROFILEDIR=/work2/05246/sh43394/frontera/hp3d/trunk/problems/MAXWELL/GALERKIN/tau/profiles
export TRACEDIR=/work2/05246/sh43394/frontera/hp3d/trunk/problems/MAXWELL/GALERKIN/tau/traces
export TAU_PROFILE=1
export TAU_TRACE=0

module list
pwd
date

# Set thread count (default value is 1)...
# export KMP_NUM_THREADS=12
nthreads=1

# Set OpenMP stack size per thread
export KMP_STACKSIZE=64M

# Launch MPI code...
ibrun ./maxw -job 1 -p 5 -isol 0 -imax 2 -maxnods 1025000 -nthreads ${nthreads}

date
# ---------------------------------------------------
