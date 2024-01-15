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
# Last revised: June 2021
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch fiber.slurm" on Stampede2 login node.
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
#SBATCH -n 64              # Total # of mpi tasks
#SBATCH -t 00:29:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-ti2           # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

module list
pwd
date

# set OpenMP threads
nthreads=12
#
# set OpenMP stack size
export KMP_STACKSIZE=24M   # p=3
#export KMP_STACKSIZE=32M   # p=4
#export KMP_STACKSIZE=48M   # p=5
#export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8
#
# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# weight for l2 term in test norm
alpha=1.0d-4

gain=1.0d4
envelope=true

# config assuming wavenum_signal=85.0d0

# set fiber length, length of PML, #refs, maxnodes
#  76_8 (1024 wavelengths), 5 refs (16*  32 elems), 8000 nodes
#zl=76.8d0; pmlfrac=0.5d0; imax=5; maxnods=8000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_76_8'
#gain=5.0d2

#  153_6 (2048 wavelengths), 6 refs (16*  64 elems), 16000 nodes
#zl=153.6d0; pmlfrac=0.25d0; imax=6; maxnods=16000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_153_6'
#gain=2.5d2

#  307_2 (4096 wavelengths), 7 refs (16*  128 elems), 32000 nodes
#zl=307.2d0; pmlfrac=0.125d0; imax=7; maxnods=32000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_307_2'
#gain=1.25d2
#
#  614_4 (8192 wavelengths), 8 refs (16*  256 elems), 63000 nodes
#zl=614.4d0; pmlfrac=0.0625d0; imax=8; maxnods=63000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_614_4'
#gain=6.25d1
#
#  1228_8 (16384 wavelengths), 9 refs (16*  512 elems), 126000 nodes
zl=1228.8d0; pmlfrac=0.03125d0; imax=9; maxnods=126000
file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1228_8'
#gain=3.125d1
#
# set job type
job=1
#
# overwrite number of max nodes
#maxnods=1505000
#
# set number of adaptive refs
jmax=1
#
# set polynomial order p in (xy,z)
px=5; py=5; pz=6
#
# set enriched order (p+dp)
dp=1
#
# set refractive index
ref_core=1.4512d0
ref_clad=1.4500d0
#ref_core=1.1520d0
#ref_clad=1.1500d0
#
# set BC flag
ibc=0
#
# set PML flag
usepml=true
#
# set nonlinear flags
nlflag=0
heat=0
#
# set pump field
plane_pump=0
plane_pump_power=1.d3
#
# set pump direction
copump=1
#
# set Raman gain
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
mkdir ${SLURM_JOB_ID}
cd ${SLURM_JOB_ID}
mkdir outputs
cd outputs
mkdir paraview
mkdir power
mkdir temp
cd ../..
dir_output="${SLURM_JOB_ID}/outputs/"
vis_level=3

args=" -geom 5 -isol 17 -alpha ${alpha} -gamma ${gamma} -comp 1"
#args+=" -omega ${omega}"
args+=" -ref_core ${ref_core} -ref_clad ${ref_clad}"
args+=" -aniso_ref_index ${aniso_ref_index}"
args+=" -art_grating ${art_grating}"
args+=" -job ${job} -imax ${imax} -jmax ${jmax}"
args+=" -ibc ${ibc}"
args+=" -px ${px} -py ${py} -pz ${pz} -dp ${dp}"
args+=" -copump ${copump} -nlflag ${nlflag} -gain ${gain} -raman ${raman}"
args+=" -plane_pump ${plane_pump} -plane_pump_power ${plane_pump_power}"
args+=" -heat ${heat} -aniso_heat ${aniso_heat} -nsteps ${nsteps} -dt ${dt}"
args+=" -dir_output ${dir_output} -vis_level ${vis_level}"
args+=" -file_geometry ${file_geometry} -zl ${zl}"
args+=" -file_control ${ctrl}"
args+=" -maxnods ${maxnods}"
args+=" -nthreads ${nthreads}"
if [ "$usepml" = true ] ; then
   args+=" -usepml -pmlfrac ${pmlfrac}"
fi
if [ "$envelope" = true ] ; then
   args+=" -envelope"
   args+=" -wavenum_signal 85.0d0 -wavenum_pump 92.7d0"
fi

# Launch MPI code...
ibrun ./uwLaser ${args}

date
# ---------------------------------------------------
