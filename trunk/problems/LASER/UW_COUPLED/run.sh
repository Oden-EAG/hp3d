#
# Script to run the laser code
#
# last modified: Jul 2018
#
# =======================
# I N S T R U C T I O N S
# =======================
#
#
# For parallel runs with OpenMP, make sure that the OMP Stacksize is set large enough.
# export KMP_STACKSIZE=512M
# export OMP_NUM_THREADS=48 (#cores on skx)
#
# In addition to the parameters in this file,
# the following parameters must be set accordingly:
#
# - NEXACT in 'control' file
# - geometry in set_environment.F90
# - for GEOM=4,5; set r_core and r_clad in LaserParam.F90
#
# TACC MPI
# ibrun -n <mpi_procs>
#
#
# ===============
# E X A M P L E S
# ===============
# 1) testing convergence rates on a cube
ibrun -n 1 ./uwLaser -geom 1 -isol 1 -p 3 -dp 1 -px 3 -py 3 -pz 3 -usepml 0 -nlflag 0 -lasermode 0 -copump 1 
#-nthreads 1
#
# 2) testing convergence rates on the fiber (core only)
#./uwLaser -geom 4 -isol 1 -p 5 -dp 1 -px 5 -py 5 -pz 5 -usepml 0 -nlflag 0 -lasermode 0 -copump 1 
#
# 3) nonlinear solve of Maxwell (signal+pump fields)
#./uwLaser -geom 5 -p 5 -dp 1 -usepml 1 -nlflag 1 -lasermode 0 -copump 1 -raman 0.0001d0
#
# 4) nonlinear solve of Maxwell weakly coupled with the heat equation
#
# ===================
# P A R A M E T E R S
# ===================
#
#
#eof
