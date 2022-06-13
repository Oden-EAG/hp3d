#!/bin/bash
#----------------------------------------------------
# Example Slurm job script
# for TACC Frontera CLX nodes
#
#   *** Hybrid Job on CLX Development Queue ***
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
#----------------------------------------------------

#SBATCH -J hp3d            # Job name
#SBATCH -o hp3d.o%j        # Name of stdout output file
#SBATCH -e hp3d.e%j        # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 2               # Total # of nodes
#SBATCH -n 8               # Total # of mpi tasks
#SBATCH -t 00:05:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A FTA-SUB-Ghattas # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

# export PROFILEDIR=/work/05246/sh43394/frontera/hp3d/trunk/problems/LASER/UW_COUPLED/tau/profiles
# export TRACEDIR=/work/05246/sh43394/frontera/hp3d/trunk/problems/LASER/UW_COUPLED/tau/traces
# export TAU_PROFILE=1
# export TAU_TRACE=1

module list
pwd
date

# set OpenMP threads
nthreads=14
#
# set OpenMP stack size
#export KMP_STACKSIZE=24M    # p=3
#export KMP_STACKSIZE=32M    # p=4
 export KMP_STACKSIZE=48M    # p=5
#export KMP_STACKSIZE=64M    # p=6
#export KMP_STACKSIZE=80M    # p=7
#export KMP_STACKSIZE=96M    # p=8
#
# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# weight for l2 term in test norm
alpha=1.0d-4

#
# Rectangular waveguide for linear Maxwell (e.g., pollution study)
#    / with dirichlet BC (-ibc 0), or impedance BC (-ibc 3), or PML
#    / -omega sqrt(5.d0)/2.d0*PI -gamma sqrt(1.d0-PI*PI/(w*w))
#      --> exactly 1 wavelength per 4 unit lengths in z-direction
#    / can be run with NEXACT=1 or NEXACT=0 (unless PML --use NEXACT=0)

# set waveguide length, #refs, maxnodes
# L=4        (1 wavelength ), 3  refs (1 aniso x,  2 aniso z --> 2 *     4  elems),    250 nodes
#zl=4.0d0    ; imax=3 ; maxnods=250    ; file_geometry='../GEOMETRIES/waveguide/rect_4'

# L=8        (2 wavelengths), 4  refs (1 aniso x,  3 aniso z --> 2 *     8  elems),    550 nodes
#zl=8.0d0    ; imax=4 ; maxnods=550    ; file_geometry='../GEOMETRIES/waveguide/rect_8'

# L=16       (4 wavelengths), 5  refs (1 aniso x,  4 aniso z --> 2 *    16  elems),   1050 nodes
#zl=16.0d0   ; imax=5 ; maxnods=1050   ; file_geometry='../GEOMETRIES/waveguide/rect_16'

# L=32       (8 wavelengths), 6  refs (1 aniso x,  5 aniso z --> 2 *    32  elems),   2050 nodes
#zl=32.0d0   ; imax=6 ; maxnods=2050   ; file_geometry='../GEOMETRIES/waveguide/rect_32'

# L=64      (16 wavelengths), 7  refs (1 aniso x,  6 aniso z --> 2 *    64  elems),   4050 nodes
#zl=64.0d0   ; imax=7 ; maxnods=4050   ; file_geometry='../GEOMETRIES/waveguide/rect_64'

# L=128     (32 wavelengths), 8  refs (1 aniso x,  7 aniso z --> 2 *   128  elems),   8050 nodes
 zl=128.0d0  ; imax=8 ; maxnods=8050   ; file_geometry='../GEOMETRIES/waveguide/rect_128'

# L=256     (64 wavelengths), 9  refs (1 aniso x,  8 aniso z --> 2 *   256  elems),  16050 nodes
#zl=256.0d0  ; imax=9 ; maxnods=16050  ; file_geometry='../GEOMETRIES/waveguide/rect_256'

# L=512    (128 wavelengths), 10 refs (1 aniso x,  9 aniso z --> 2 *   512  elems),  32050 nodes
#zl=512.0d0  ; imax=10; maxnods=32050  ; file_geometry='../GEOMETRIES/waveguide/rect_512'

# L=1024   (256 wavelengths), 11 refs (1 aniso x, 10 aniso z --> 2 *  1024  elems),  64050 nodes
#zl=1024.0d0 ; imax=11; maxnods=64050  ; file_geometry='../GEOMETRIES/waveguide/rect_1024'

# L=2048   (512 wavelengths), 12 refs (1 aniso x, 11 aniso z --> 2 *  2048  elems), 105000 nodes
#zl=2048.0d0 ; imax=12; maxnods=105000 ; file_geometry='../GEOMETRIES/waveguide/rect_2048'

# L=4096  (1024 wavelengths), 13 refs (1 aniso x, 12 aniso z --> 2 *  4096  elems), 205000 nodes
#zl=4096.0d0 ; imax=13; maxnods=205000 ; file_geometry='../GEOMETRIES/waveguide/rect_4096'

# L=8192  (2048 wavelengths), 14 refs (1 aniso x, 13 aniso z --> 2 *  8192  elems), 405000 nodes
#zl=8192.0d0 ; imax=14; maxnods=405000 ; file_geometry='../GEOMETRIES/waveguide/rect_8192'

# L=16384 (4096 wavelengths), 15 refs (1 aniso x, 14 aniso z --> 2 * 16384  elems), 805000 nodes
#zl=16384.0d0; imax=15; maxnods=805000 ; file_geometry='../GEOMETRIES/waveguide/rect_16384'

# L=32768 (8192 wavelengths), 16 refs (1 aniso x, 15 aniso z --> 2 * 32768  elems), 1605000 nodes
#zl=32768.0d0; imax=16; maxnods=1605000; file_geometry='../GEOMETRIES/waveguide/rect_32768'

# set wave number, #refs, maxnodes, 4 elems/wavelength
#zl=4.0d0; file_geometry='../GEOMETRIES/waveguide/rect_4'
#    1 wavelength  (  0.25 waves/unit length)
#imax=3 ; maxnods=250    ; omega=3.512407365520363196578d0
#    2 wavelengths (  0.50 waves/unit length)
#imax=4 ; maxnods=550    ; omega=4.442882938158366247015d0
#    4 wavelengths (  1.00 waves/unit length)
#imax=5 ; maxnods=1050   ; omega=7.024814731040726393156d0
#    8 wavelengths (  2.00 waves/unit length)
#imax=6 ; maxnods=2050   ; omega=12.95311834341519078982d0
#   16 wavelengths (  4.00 waves/unit length)
#imax=7 ; maxnods=4050   ; omega=25.32832971340211194512d0
#   32 wavelengths (  8.00 waves/unit length)
#imax=8 ; maxnods=8050   ; omega=50.36356154085972864896d0
#   64 wavelengths ( 16.00 waves/unit length)
#imax=9 ; maxnods=16050  ; omega=100.5800403217089219730d0
#  128 wavelengths ( 32.00 waves/unit length)
#imax=10; maxnods=32050  ; omega=201.0864720245076612134d0
#  256 wavelengths ( 64.00 waves/unit length)
#imax=11; maxnods=64050  ; omega=402.1361313185488198126d0
#  512 wavelengths (128.00 waves/unit length)
#imax=12; maxnods=105000 ; omega=804.2538552187321009202d0
# 1024 wavelengths (256.00 waves/unit length)
#imax=13; maxnods=205000 ; omega=1608.498506596624078797d0
# 2048 wavelengths (512.00 waves/unit length)
#imax=14; maxnods=405000 ; omega=3216.992411256370432481d0
# 4096 wavelengths (1024.00 waves/unit length)
#imax=15; maxnods=805000 ; omega=6433.982521542244779007d0
# 8192 wavelengths (2048.00 waves/unit length)
#imax=16; maxnods=1525000; omega=12867.96389259898436162d0
# 16384 wavelengths (4096.00 waves/unit length)
#imax=17; maxnods=3025000; omega=25735.92720995518398087d0
# 32768 wavelengths (8192.00 waves/unit length)
#imax=18; maxnods=6025000; omega=51471.85413228897157253d0
# 65536 wavelengths (16384.00 waves/unit length)
#imax=19; maxnods=12025000; omega=102943.7081207672444482d0

#
# set job type
job=1
#
# set polynomial order p in (xy,z)
px=5; py=5; pz=5
#
# set enriched order (p+dp)
dp=1
#
# set BC flag
# 0: Dirichlet
# 2: impedance via penalty method
# 3: impedance via elimination
ibc=2
#
# set PML flag
usepml=false
pmlfrac=0.25d0
#
# other params
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

args=" -geom 1 -isol 5 -comp 2 -alpha ${alpha}"
args+=" -job ${job} -imax ${imax} -jmax 0"
args+=" -px ${px} -py ${py} -pz ${pz} -dp ${dp}"
args+=" -ibc ${ibc}"
args+=" -nlflag 0 -heat 0 -copump 1"
args+=" -dir_output ${dir_output} -vis_level ${vis_level}"
args+=" -file_geometry ${file_geometry} -zl ${zl}"
args+=" -file_control ${ctrl}"
args+=" -maxnods ${maxnods}"
args+=" -nthreads ${nthreads}"
if [ "$usepml" = true ] ; then
   args+=" -usepml -pmlfrac ${pmlfrac}"
fi

# ADD optionally -omega ${omega}

# Launch MPI code...
ibrun ./uwLaser ${args}

date
# ---------------------------------------------------
