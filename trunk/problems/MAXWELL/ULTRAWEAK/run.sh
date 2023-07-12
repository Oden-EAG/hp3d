#
# Script to run the Maxwell problem
#
# last modified: July 2023
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
vis_level=0

# MPI Procs
nproc=4

# OMP THREADS
nthreads=1

# Set polynomial order p
p=2

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=1
imax=0
jmax=1

# max NODES
maxnods=100000

#export KMP_STACKSIZE=24M   # p=3
export KMP_STACKSIZE=32M   # p=4
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
isol=4

# component number for manufactured solution
comp=1

#
# ==================
# RUN CONFIGURATIONS
# ==================
#
# NEW RUN CONFIGS (MPI+OpenMP CODE)
# ===========================================================================
file_geometry='./geometries/cube'
ctrl='control/control_1'
#
args=" -isol ${isol} -omega ${omega} -gamma 1.0d0 -comp ${comp}"
args+=" -job ${job} -imax ${imax} -jmax ${jmax}"
args+=" -ibc ${ibc}"
args+=" -p ${p} -dp ${dp}"
args+=" -dir_output ${dir_output} -vis_level ${vis_level}"
args+=" -file_geometry ${file_geometry}"
args+=" -file_control ${ctrl}"
args+=" -maxnods ${maxnods}"
args+=" -nthreads ${nthreads} -alpha ${alpha}"

mpirun -np ${nproc} ./maxwell_UW ${args}
