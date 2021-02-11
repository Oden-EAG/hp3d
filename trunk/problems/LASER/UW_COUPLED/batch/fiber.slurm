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
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 4               # Total # of mpi tasks
#SBATCH -t 00:15:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
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

# set OpenMP threads
nthreads=12
#
# set OpenMP stack size
#export KMP_STACKSIZE=24M   # p=3
#export KMP_STACKSIZE=32M   # p=4
export KMP_STACKSIZE=48M   # p=5
#export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8
#
# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE
#
gain=1.0d4

# set fiber length, length of PML, #refs, maxnodes
#  1_2 (16 wavelengths), 5 refs (16*  32 elems), 6250 nodes
#zl=1.2d0; pmlfrac=0.25d0; imax=5; maxnods=6250
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1_2'
#gain=4.0d4

#  2_4 (32 wavelengths), 6 refs (16*  64 elems), 12500 nodes
zl=2.4d0; pmlfrac=0.125d0; imax=6; maxnods=12500
file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_2_4'
#gain=2.0d4

#  4_8 (64 wavelengths), 7 refs (16*  128 elems), 25000 nodes
#zl=4.8d0; pmlfrac=0.0625d0; imax=7; maxnods=25000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_4_8'
#gain=1.0d4

#  9_6 (128 wavelengths), 8 refs (16*  256 elems), 51000 nodes
#zl=9.6d0; pmlfrac=0.03125d0; imax=8; maxnods=51000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_9_6'
#gain=5.0d3

#  19_2 (256 wavelengths), 9 refs (16*  512 elems), 101000 nodes
#zl=19.2d0; pmlfrac=0.015625d0; imax=9; maxnods=101000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_19_2'
#gain=2.5d3

#  38_4 (512 wavelengths), 10 refs (16*  1024 elems), 201000 nodes
#zl=38.4d0; pmlfrac=8.0d-3; imax=10; maxnods=201000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_38_4'
#gain=1.25d3

#  76_8 (1024 wavelengths), 11 refs (16*  2048 elems), 401000 nodes
#zl=76.8d0; pmlfrac=4.0d-3; imax=11; maxnods=401000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_76_8'
#gain=6.25d2

#  153_6 (2048 wavelengths), 12 refs (16*  4096 elems), 801000 nodes
#zl=153.6d0; pmlfrac=2.0d-3; imax=12; maxnods=801000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_153_6'
#gain=3.125d2

#  307_2 (4096 wavelengths), 13 refs (16*  8192 elems), 1601000 nodes
#zl=307.2d0; pmlfrac=1.0d-3; imax=13; maxnods=1601000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_307_2'
#gain=1.5625d2
#
# set job type
job=1
#
# overwrite number of max nodes
#maxnods=1505000
#
# set number of adaptive refs
jmax=0
#
# set polynomial order p in (xy,z)
px=6; py=6; pz=6
#
# set enriched order (p+dp)
dp=1
#
# set refractive index
ref_core=1.4512d0
ref_clad=1.4500d0
#
# set BC flag
ibc=0
#
# set PML flag
usepml=1
#
# set nonlinear flags
nlflag=0
heat=0
#
# set Raman gain for nonlinear case
raman=0.d0
#
# set anisotropic heat operator
aniso_heat=1
# set anisotropic refractive index (CORE_NX,CORE_NY,CORE_NZ)
aniso_ref_index=0
# set artificial refractive index grating
art_grating=0
#
# set number and size of steps
nsteps=10
dt=0.1d0
#
# other params
gamma=1.0d0
ctrl='../COMMON_FILES/control_0'
#
# set output
dir_output='../outputs/'

# Launch MPI code...
ibrun ./uwLaser -geom 5 -isol 13 -comp 1 -ref_core ${ref_core} -ref_clad ${ref_clad} -gamma ${gamma} -file_control ${ctrl} -job ${job} -maxnods ${maxnods} -imax ${imax} -jmax ${jmax} -px ${px} -py ${py} -pz ${pz} -dp ${dp} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag ${nlflag} -gain ${gain} -raman ${raman} -heat ${heat} -aniso_heat ${aniso_heat} -aniso_ref_index ${aniso_ref_index} -art_grating ${art_grating} -nsteps ${nsteps} -dt ${dt} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -file_geometry ${file_geometry}

date
# ---------------------------------------------------
