#
# Script to run the laser code
#
# last modified: Apr 2019
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
# - geometry in set_environment.F90
# - for GEOM=4,5; set r_core and r_clad in laserParam.F90
# - OMEGA/GAMMA; set frequency and BCs in set_environment.F90
#
# ===================
# P A R A M E T E R S
# ===================

# paraview parameters
dir_output='../outputs/'
vis_level=3
#
# MPI Procs
nproc=1
#
# OMP THREADS
nthreads=12
#
# set polynomial order p in (xy,z)
px=5; py=5; pz=5
# set enriched order (p+dp)
dp=1
# set length of waveguide
zl=1.2d0
# set length of PML
usepml=1
pmlfrac=0.25d0         #  1_2,   16 wavelengths,   6250 nodes,  5 refs (16*  32 elems)
#pmlfrac=0.125d0        #  2_4,   32 wavelengths,  12500 nodes,  6 refs (16*  64 elems)
#pmlfrac=0.0625d0        #  4_8,   64 wavelengths,  25000 nodes,  7 refs (16* 128 elems)
#pmlfrac=0.03125d0      #  9_6,  128 wavelengths,  51000 nodes,  8 refs (16* 256 elems)
#pmlfrac=0.015625d0     # 19_2,  256 wavelengths, 101000 nodes,  9 refs (16* 512 elems)
#pmlfrac=0.0078125d0    # 38_4,  512 wavelengths, 201000 nodes, 10 refs (16*1024 elems)
#pmlfrac= 0.00390625d0  # 76_8, 1024 wavelengths, 401000 nodes, 11 refs (16*2048 elems)
#
job=1
imax=5
maxnods=25000
#
#export KMP_STACKSIZE=128M    # p=4
export KMP_STACKSIZE=192M     # p=5
#
# set BC flag
ibc=0
#
# set nonlinear flags
nlflag=0
heat=0
#
# set anisotropic heat operator
aniso_heat=1
aniso_ref_index=0
#
# set gain for nonlinear case
gain=1.0d4
raman=0.d0
#
# set number and size of steps
nsteps=200
dt=0.1d0
#
# ==================
# RUN CONFIGURATIONS
# ==================
#
# NEW RUN CONFIGS (MPI+OpenMP CODE)

# TESTING (CONVERGENCE RATES ON A CUBE) -- use NEXACT=1
file_geometry='../GEOMETRIES/cubes/cube'
#ibrun -n ${nproc} ./uwLaser -geom 1 -isol 2 -omega 1.0d0 -comp 1 -job ${job} -maxnods ${maxnods} -imax ${imax} -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -npx 3 -npy 3 -npz 3 -usepml 0 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -file_geometry ${file_geometry}
#ibrun -n ${nproc} xterm -hold -e ./uwLaser -geom 1 -isol 2 -omega 1.0d0 -comp 1 -job ${job} -imax ${imax} -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -npx 3 -npy 3 -npz 3 -usepml 0 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -file_geometry ${file_geometry}

# FULL FIBER (non-)linear solve of Maxwell (with heat/gain) (signal+pump fields) -- use NEXACT=0
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1_2'
file_geometry='../GEOMETRIES/fiber/fiber_prism_1_2'
ibrun -n ${nproc} ./uwLaser -geom 5 -isol 13 -comp 1 -job ${job} -maxnods ${maxnods} -imax ${imax} -px ${px} -py ${py} -pz ${pz} -dp ${dp} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag ${nlflag} -gain ${gain} -raman ${raman} -heat ${heat} -aniso_heat ${aniso_heat} -aniso_ref_index ${aniso_ref_index} -nsteps ${nsteps} -dt ${dt} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -file_geometry ${file_geometry}

# OLD RUN CONFIGS (OpenMP CODE)
#
# 1) testing convergence rates on a cube -- use NEXACT=1
#file_geometry='../GEOMETRIES/cubes/cube'
#./uwLaser -geom 1 -isol 2 -omega 1.0d0 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -npx 3 -npy 3 -npz 3 -usepml 0 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 2) testing convergence rates on the fiber (core only) -- use NEXACT=1
file_geometry='../GEOMETRIES/fiber/fhcor_curv'
#./uwLaser -geom 4 -isol 2 -omega 0.1d0 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc 0 -usepml 0 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 3) (non-)linear solve of Maxwell (with heat/gain) (signal+pump fields)
#file_geometry='../GEOMETRIES/fiber/fiber_hexa9'
#file_geometry='../GEOMETRIES/fiber/fiber_hexa13'
file_geometry='../GEOMETRIES/fiber/fiber_prism'
#ibrun -n ${nproc} ./uwLaser -geom 5 -isol 13 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag ${nlflag} -gain ${gain} -raman ${raman} -heat ${heat} -aniso_heat ${aniso_heat} -aniso_ref_index ${aniso_ref_index} -nsteps ${nsteps} -dt ${dt} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 4) rectangular waveguide for linear Maxwell
file_geometry='../GEOMETRIES/waveguide/rect'
#    / with dirichlet BC (-ibc 0), or impedance BC (-ibc 3), or PML
#    / -omega sqrt(5.d0)/2.d0*PI -gamma sqrt(1.d0-PI*PI/(w*w))
#      --> 1 wavelength per 4 unit lengths in z-direction
#    / can be run with NEXACT=1 or NEXACT=0 (unless PML --use NEXACT=0)
#./uwLaser -geom 1 -isol 3 -comp 2 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 5) circular hollow waveguide for linear Maxwell
file_geometry='../GEOMETRIES/fiber/fhcor_curv'
#./uwLaser -geom 4 -isol 12 -comp 1 -bessord 1.d0 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 6) dielectric fiber waveguide for linear Maxwell
file_geometry='../GEOMETRIES/fiber/fiber_prism'
#./uwLaser -geom 5 -isol 13 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml ${usepml} -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
#eof

