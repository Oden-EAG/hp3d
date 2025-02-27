#
# Script to run the Maxwell problem with bending and envelope ansatz
#
# last modified: October 2024
#
# =======================
# I N S T R U C T I O N S
# =======================
#
# For parallel runs with OpenMP, make sure that the OMP Stacksize is set large enough.
#
# In addition to the parameters in this file,
# the following parameters must be set accordingly:
#
# - NEXACT in 'control' file
#
# ===================
# P A R A M E T E R S
# ===================

# paraview parameters
dir_output='../outputs/'
vis_level=2

# MPI Procs
nproc=4

# OMP THREADS
nthreads=1

# Set polynomial order p
p=1

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=0
imax=3

# max NODES
maxnods=123456

# export KMP_STACKSIZE=24M   # p=3
#export KMP_STACKSIZE=32M   # p=4
export KMP_STACKSIZE=48M   # p=5
#export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8

# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# component number for manufactured solution
comp=1
# solution number
isol=15
# DPG test norm scaling
alpha=0.01d0
# vacuum permeability
mu=7.91582400800000E-01
# vacuum permittivity
epsilon=9.06838390340768E-01
# angular frequency
omega=1.77034921739554E+02
# impedance coefficient
gamma=1.07032799231899E+00
# base electric field
e0=9.66588353875592E-04
# base magnetic field
h0=1.03456657220257E-03
# envelope wavenumber
k=1.49393360127023E+02   # exactly k = omega/c
# bending radius
rbend=1300.0
# half width
halfwidth=0.5 # 12.7d0
# spanning angle of bent fiber
spangle=0.000174532925199 # 0.01 DEGREES
# set BC flag -> 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=0

# prefix for paraview
pref='bslab_Dir1_EnvAns1' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir1_EnvAns0' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir2_EnvAns1' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir2_EnvAns0' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'

#
# ==================
# RUN CONFIGURATIONS
# ==================
#
# NEW RUN CONFIGS (MPI+OpenMP CODE)
# ===========================================================================
# file_geometry='./geometries/torus_part2'
# file_geometry='./geometries/bent_fiber_test_30'
# file_geometry='./geometries/partly_bent_slab_d2_cR420_DEG30'
# file_geometry='./geometries/bent_waveguide_sq1_R10_DEG90'
# file_geometry='./geometries/bent_waveguide_sq1_R10_DEG30'
# file_geometry='./geometries/bent_waveguide_sq1_R500_DEG1'
# file_geometry='./geometries/bent_waveguide_d2_cR420_DEG30'
# file_geometry='./geometries/bent_waveguide_d12p7_cR35000_DEG5'
# file_geometry='./geometries/bent_waveguide_d0p5_cR5_DEG5' 
file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'


ctrl='control/control'
#
args=" -file_control ${ctrl}"
args+=" -file_geometry ${file_geometry}"
args+=" -p ${p} -dp ${dp}"
args+=" -comp ${comp} -isol ${isol} -imax ${imax} -job ${job}"
args+=" -maxnods ${maxnods} -alpha ${alpha}"
args+=" -mu ${mu} -epsilon ${epsilon} -e0 ${e0} -h0 ${h0}"
args+=" -omega ${omega} -gamma ${gamma}"
args+=" -k ${k} -rbend ${rbend} -halfwidth ${halfwidth} -spangle ${spangle}"
args+=" -ibc ${ibc}"
args+=" -prefix ${pref} -vis_level ${vis_level} -dir_output ${dir_output}"
args+=" -nthreads ${nthreads}"

mpirun -np ${nproc} ./bending ${args}
