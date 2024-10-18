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
nproc=1

# OMP THREADS
nthreads=1

# Set polynomial order p
p=2

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=0
imax=3

# max NODES
maxnods=123456

export KMP_STACKSIZE=24M   # p=3
#export KMP_STACKSIZE=32M   # p=4
#export KMP_STACKSIZE=48M   # p=5
#export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8

# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# set BC flag
# 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=0

# DPG test norm scaling
alpha=1.d0

# solution number
isol=10

# component number for manufactured solution
comp=1

# angular frequency
omega=6.283185307179586d0

# envelope wavenumber
k=4.d0
# bending radius
rbend=10.d0

#
# ==================
# RUN CONFIGURATIONS
# ==================
#
# NEW RUN CONFIGS (MPI+OpenMP CODE)
# ===========================================================================
file_geometry='./geometries/bent_square_waveguide_new'
ctrl='control/control'
#
args=" -file_control ${ctrl}"
args+=" -file_geometry ${file_geometry}"
args+=" -p ${p} -dp ${dp}"
args+=" -comp ${comp} -isol ${isol} -imax ${imax} -job ${job}"
args+=" -maxnods ${maxnods} -alpha ${alpha} -omega ${omega} -gamma 1.0d0"
args+=" -k ${k} -rbend ${rbend}"
args+=" -ibc ${ibc}"
args+=" -vis_level ${vis_level} -dir_output ${dir_output}"
args+=" -nthreads ${nthreads}"

mpirun -np ${nproc} ./bending ${args}
