#!/bin/bash
#----------------------------------------------------
# Example Slurm job script
# for TACC Frontera CLX nodes
#
#   *** Hybrid Job on CLX Normal Queue ***
#
#       This sample script specifies:
#         10 nodes (capital N)
#         40 total MPI tasks (lower case n); this is 4 tasks/node
#         14 OpenMP threads per MPI task (56 threads per node)
#
# Last revised: June 2021
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch fiber.slurm" on Frontera login node.
#
#   -- Check current queue with: "squeue -u sh43394"
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- In most cases it's best to keep
#      ( MPI ranks per node ) x ( threads per rank )
#      to a number no more than 56 (total cores).
#
#   -- If you're running out of memory, try running
#      fewer tasks and/or threads per node to give each
#      process access to more memory.
#
#   -- IMPI and MVAPICH2 both do sensible process pinning by default.
#
#   -- QUEUES: 
#        development  1-40  | max  2 hrs
#        small        1-2   | max 48 hrs
#        normal       3-512 | max 48 hrs 
#        large      513-2048| max 48 hrs
#
#----------------------------------------------------

#SBATCH -J hp3d            # Job name
#SBATCH -o hp3d.o%j        # Name of stdout output file
#SBATCH -e hp3d.e%j        # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 16              # Total # of nodes
#SBATCH -n 64              # Total # of mpi tasks
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DMS22025        # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

module list
pwd
date

# set OpenMP threads
nthreads=14
#
# set OpenMP stack size
export KMP_STACKSIZE=24M   # p=6
#export KMP_STACKSIZE=48M   # p=7
#
# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# weight for l2 term in test norm
alpha=1.0d-4

# Set up envelope formulation:
envelope=true
envelope_config=3
# envelope configs:
#   #1
#   assuming wavenum_signal=85.0d0
#   LP01 envelope beat ca. 125 wavelengths
#   LP02 envelope beat ca. 135 wavelengths
#   use PML length of ~0.4 mm (ca. 545 wavelengths)
if [ "$envelope_config" = 1 ] ; then
   wavenum_signal=85.0d0 
fi
#   #2
#   assuming wavenum_signal=85.5119d0
#   LP01 envelope beat ca. 500 wavelengths
#   LP02 envelope beat ca. 710 wavelengths
#   use PML length of ~1.6 mm (ca. 2180 wavelengths)
if [ "$envelope_config" = 2 ] ; then
   wavenum_signal=85.5119d0 
fi
#   #3
#   assuming wavenum_signal=85.5976d0
#   LP01 envelope beat ca. 1000 wavelengths
#   LP02 envelope beat ca. 2475 wavelengths
#   use PML length of ~6.4 mm (ca. 8720 wavelengths)
if [ "$envelope_config" = 3 ] ; then
   wavenum_signal=85.5976d0 
fi

# set fiber length, length of PML, #refs, maxnodes
if [ "$envelope_config" = 1 ] ; then
   # 1 NODE
   #  76_8 (1024 wavelengths), 5 refs (16*  32 elems), 8000 nodes
   #zl=76.8d0; pmlfrac=0.5d0; imax=5; maxnods=8000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_76_8'

   # 2 NODES
   #  153_6 (2048 wavelengths), 6 refs (16*  64 elems), 16000 nodes
   #zl=153.6d0; pmlfrac=0.25d0; imax=6; maxnods=16000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_153_6'

   # 4 NODES
   #  307_2 (4096 wavelengths), 7 refs (16*  128 elems), 32000 nodes
   #zl=307.2d0; pmlfrac=0.125d0; imax=7; maxnods=32000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_307_2'
   #
   # 8 NODES
   #  650_0 (6.0mm = 8160 wavelengths), 8 refs (16*  256 elems), 63000 nodes
   #zl=650.0d0; pmlfrac=0.0625d0; imax=8; maxnods=63000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_650'
   #
   # 16 NODES
   #  1300_0 (1.25cm = 17,000 wavelengths), 9 refs (16*  512 elems), 126000 nodes
   #zl=1300.0d0; pmlfrac=0.03125d0; imax=9; maxnods=126000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1300'
   #
   # 32 NODES
   #  2600_0 (2.5cm = 34,000 wavelengths), 10 refs (16*  1024 elems), 251000 nodes
   zl=2600.0d0; pmlfrac=0.015625d0; imax=10; maxnods=251000
   file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_2600'
   #
   # 64 NODES
   #  5100_0 (5.0cm = 68,000 wavelengths), 11 refs (16*  2048 elems), 501000 nodes
   #zl=5100.0d0; pmlfrac=7.8125d-3; imax=11; maxnods=501000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_5100'
   #
   # 128 NODES
   #  10100_0 (10.0cm = 136000 wavelengths), 12 refs (16*  4096 elems), 1001000 nodes
   #zl=10100.0d0; pmlfrac=3.90625d-3; imax=12; maxnods=1001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_10100'
   #
   # 256 NODES
   #  20100_0 (20.0cm = 272000 wavelengths), 13 refs (16*  8192 elems), 2001000 nodes
   #zl=20100.0d0; pmlfrac=1.953125d-3; imax=13; maxnods=2001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_20100'
   #
   # 512 NODES
   #  40100_0 (40.0cm = 544000 wavelengths), 14 refs (16*  16384 elems), 4001000 nodes
   #zl=40100.0d0; pmlfrac=9.765625d-4; imax=14; maxnods=4001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_40100'
fi

if [ "$envelope_config" = 2 ] ; then
   # 1 NODE
   #  L=400 (2.0 mm = 2720 wavelengths), 5 refs (16*  32 elems), 8000 nodes
   #zl=400.d0; pmlfrac=0.5d0; imax=5; maxnods=8000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_400'
   #
   # 2 NODES
   #  L=800 (6.0 mm = 8160 wavelengths), 6 refs (16*  64 elems), 16000 nodes
   #zl=800.0d0; pmlfrac=0.25d0; imax=6; maxnods=16000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_800'
   #
   # 4 NODES
   #  L=1500 (1.25 cm = 17,000 wavelengths), 7 refs (16*  128 elems), 32000 nodes
   #zl=1500.0d0; pmlfrac=0.125d0; imax=7; maxnods=32000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1500'
   #
   # 8 NODES
   #  L=2750 (2.5 cm = 34,000 wavelengths), 8 refs (16*  256 elems), 63000 nodes
   #zl=2750.0d0; pmlfrac=0.0625d0; imax=8; maxnods=63000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_2750'
   #
   # 16 NODES
   #  L=5250 (5 cm = 68,000 wavelengths), 9 refs (16*  512 elems), 126000 nodes
   #zl=5250.0d0; pmlfrac=0.03125d0; imax=9; maxnods=126000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_5250'
   #
   # 32 NODES
   #  L=10250 (10 cm = 136,000 wavelengths), 10 refs (16*  1024 elems), 251000 nodes
   #zl=10250.0d0; pmlfrac=0.015625d0; imax=10; maxnods=251000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_10250'
   #
   # 64 NODES
   #  L=20250 (20 cm = 272,000 wavelengths), 11 refs (16*  2048 elems), 501000 nodes
   zl=20250.0d0; pmlfrac=7.8125d-3; imax=11; maxnods=501000
   file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_20250'
   #
   # 128 NODES
   #  L=40250 (40 cm = 544,000 wavelengths), 12 refs (16*  4096 elems), 1001000 nodes
   #zl=40250.0d0; pmlfrac=3.90625d-3; imax=12; maxnods=1001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_40250'
   #
   # 256 NODES
   #  L=80250 (80 cm = 1,088,000 wavelengths), 13 refs (16*  8192 elems), 2001000 nodes
   #zl=80250.0d0; pmlfrac=1.953125d-3; imax=13; maxnods=2001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_80250'
   #
   # 512 NODES
   #  L=160250 (160 cm = 2,176,000 wavelengths), 14 refs (16*  16384 elems), 4001000 nodes
   #zl=160250.0d0; pmlfrac=9.765625d-4; imax=14; maxnods=4001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_160250'
fi

# NOTE: p=666 appears too low for longer than 40 cm fibers
#       p=777 needs twice the number nodes
if [ "$envelope_config" = 3 ] ; then
   # 4 NODES
   #  L=3000 (2 cm = 27,000 wavelengths), 7 refs (16*  128 elems), 32000 nodes
   #zl=3000.0d0; pmlfrac=0.25d0; imax=7; maxnods=32000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_3000'
   #
   # 8 NODES
   #  L=6000 (5 cm = 68,000 wavelengths), 8 refs (16*  256 elems), 63000 nodes
   #zl=6000.0d0; pmlfrac=0.125d0; imax=8; maxnods=63000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_6000'
   #
   # 16 NODES
   #  L=11000 (10 cm = 136,000 wavelengths), 9 refs (16*  512 elems), 126000 nodes
   zl=11000.0d0; pmlfrac=0.0625d0; imax=9; maxnods=126000
   file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_11000'
   #
   # 32 NODES p=666 | 64 NODES p=777
   #  L=21000 (20 cm = 272,000 wavelengths), 10 refs (16*  1024 elems), 251000 nodes
   #zl=21000.0d0; pmlfrac=0.03125d0; imax=10; maxnods=251000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_21000'
   #
   # 64 NODES p=666 | 128 NODES p=777
   #  L=41000 (40 cm = 544,000 wavelengths), 11 refs (16*  2048 elems), 501000 nodes
   #zl=41000.0d0; pmlfrac=0.015625d0; imax=11; maxnods=501000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_41000'
   #
   # 128 NODES p=666 | 256 NODES p=777
   #  L=81000 (80 cm = 1,088,000 wavelengths), 12 refs (16*  4096 elems), 1001000 nodes
   #zl=81000.0d0; pmlfrac=7.8125d-3; imax=12; maxnods=1001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_81000'
   #
   # 256 NODES p=666 | 512 NODES p=777
   #  L=161000 (160 cm = 2,176,000 wavelengths), 13 refs (16*  8192 elems), 2001000 nodes
   #zl=161000.0d0; pmlfrac=3.90625d-3; imax=13; maxnods=2001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_161000'
   #
   # 512 NODES p=666 | 1024 NODES p=777
   #  L=321000 (320 cm = 4,352,000 wavelengths), 14 refs (16*  16384 elems), 4001000 nodes
   #zl=321000.0d0; pmlfrac=1.953125d-3; imax=14; maxnods=4001000
   #file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_321000'
fi

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
px=6; py=6; pz=6
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
nlflag=1
heat=1
#
# set pump field
# 0: Maxwell
# 1: Constant
# 2: ODE model
plane_pump=2
plane_pump_power=1.0d3
#
# set pump direction
copump=1
#
# set active gain
gain=1.0d0
#
# set Raman gain
raman=0.d0
#
# set anisotropic heat operator
aniso_heat=0
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
vis_level=2

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
   args+=" -wavenum_signal ${wavenum_signal} -wavenum_pump 92.7d0"
fi

# Launch MPI code...
ibrun ./uwLaserThermal ${args}

date
# ---------------------------------------------------
