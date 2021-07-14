#
# Script to run the laser code
#
# last modified: June 2021
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
nthreads=1
#
# set polynomial order p in (xy,z)
#px=1; py=1; pz=1
#px=2; py=2; pz=2
#px=3; py=3; pz=3
px=4; py=4; pz=4
#px=5; py=5; pz=5
#px=6; py=6; pz=6
#px=7; py=7; pz=7
#px=8; py=8; pz=8
#
# set enriched order (p+dp)
dp=1
#
# set length of waveguide
zl=1.0d0
#
# set if PML should be used (true/false)
usepml=false

# EXEC JOB options
job=0
imax=0
jmax=0

# max NODES
maxnods=105000

export KMP_STACKSIZE=24M   # p=3
#export KMP_STACKSIZE=32M   # p=4
#export KMP_STACKSIZE=48M   # p=5
#export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8

# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

#
# set BC flag
# 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=0
#
# set refractive index of the fiber
ref_core=1.4512d0
ref_clad=1.4500d0
#ref_core=1.1520d0
#ref_clad=1.1500d0
#
# set nonlinear flags
nlflag=0
heat=0
#
# set anisotropic heat operator
aniso_heat=1
#
# set anisotropic refractive index (CORE_NX,CORE_NY,CORE_NZ)
aniso_ref_index=0
#
# set artificial refractive index grating
art_grating=0
#
# set gain for nonlinear case
gain=1.0d4
raman=0.d0
#
# set number and size of time steps
nsteps=10
dt=0.1d0
#
# envelope formulation
envelope=false
#
# ==================
# RUN CONFIGURATIONS
# ==================
#
# NEW RUN CONFIGS (MPI+OpenMP CODE)

# ============================================================================================
# A. TESTING (CONVERGENCE RATES ON A CUBE/PRISM) -- use NEXACT=1
file_geometry='../GEOMETRIES/cubes/cube'
#file_geometry='../GEOMETRIES/prisms/prism_curv1'
ctrl='../COMMON_FILES/control_1'
#
args=" -geom 1 -isol 1 -omega 1.0d0 -gamma 1.0d0 -comp 1"
args+=" -job ${job} -imax ${imax} -jmax ${jmax}"
args+=" -ibc ${ibc}"
args+=" -px ${px} -py ${py} -pz ${pz} -dp ${dp} -npx 4 -npy 4 -npz 4"
args+=" -nlflag 0 -heat ${heat} -aniso_heat ${aniso_heat}"
args+=" -dir_output ${dir_output} -vis_level ${vis_level}"
args+=" -file_geometry ${file_geometry} -zl ${zl}"
args+=" -file_control ${ctrl}"
args+=" -maxnods ${maxnods}"
args+=" -nthreads ${nthreads}"
if [ "$envelope" = true ] ; then
   args+=" -envelope"
   args+=" -wavenum_signal 1.0d0 -wavenum_pump 1.0d0"
fi

mpirun -np ${nproc} ./uwLaser ${args}
#ibrun -n ${nproc} ./uwLaser ${args}
#ibrun -n ${nproc} xterm -hold -e ./uwLaser ${args}

# ============================================================================================
# B. FULL FIBER (non-)linear solve of Maxwell (with heat/gain) (signal+pump fields) -- use NEXACT=0
# set fiber length, length of PML, #refs, maxnodes
#  1_2 (16 wavelengths), 5 refs (16*  32 elems), 6250 nodes
zl=1.2d0; pmlfrac=0.25d0; imax=5; maxnods=8000 #6250
file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_1_2'
#gain=3.2d4
#dir_output='../gain/outputs16/'

#  2_4 (32 wavelengths), 6 refs (16*  64 elems), 12500 nodes
#zl=2.4d0; pmlfrac=0.125d0; imax=6; maxnods=16000 #12500
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_2_4'
#gain=1.6d4
#dir_output='../gain/outputs32/'

#  4_8 (64 wavelengths), 7 refs (16*  128 elems), 25000 nodes
#zl=4.8d0; pmlfrac=0.0625d0; imax=7; maxnods=32000 #25000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_4_8'
#gain=8.0d3
#dir_output='../gain/outputs64/'

#  9_6 (128 wavelengths), 8 refs (16*  256 elems), 51000 nodes
#zl=9.6d0; pmlfrac=0.03125d0; imax=8; maxnods=63000 #51000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_9_6'
#gain=4.0d3
#dir_output='../gain/outputs128/'

#  19_2 (256 wavelengths), 9 refs (16*  512 elems), 101000 nodes
#zl=19.2d0; pmlfrac=0.015625d0; imax=9; maxnods=126000 #101000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_19_2'
#gain=2.0d3
#dir_output='../gain/outputs256/'

#  38_4 (512 wavelengths), 10 refs (16*  1024 elems), 201000 nodes
#zl=38.4d0; pmlfrac=8.0d-3; imax=10; maxnods=251000 #201000
#file_geometry='../GEOMETRIES/fiber/fiber_prism/fiber_prism_38_4'
#gain=1.0d3
#dir_output='../gain/outputs512/'

#omega=59.05249348852994809140307d0
gamma=1.0d0
ctrl='../COMMON_FILES/control_0'
usepml=true
#
args=" -geom 5 -isol 13 -gamma ${gamma} -comp 1"
#args+=" -omega ${omega}"
args+=" -ref_core ${ref_core} -ref_clad ${ref_clad}"
args+=" -aniso_ref_index ${aniso_ref_index}"
args+=" -art_grating ${art_grating}"
args+=" -job ${job} -imax ${imax} -jmax ${jmax}"
args+=" -ibc ${ibc}"
args+=" -px ${px} -py ${py} -pz ${pz} -dp ${dp} -npx 4 -npy 4 -npz 4"
args+=" -copump 1 -nlflag ${nlflag} -gain ${gain} -raman ${raman}"
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
   args+=" -wavenum_signal 1.0d0 -wavenum_pump 1.0d0"
fi

#mpirun -np ${nproc} ./uwLaser ${args}
#ibrun -n ${nproc} xterm -hold -e ./uwLaser ${args}
#ibrun -n ${nproc} ./uwLaser ${args}
#valgrind --tool=memcheck --leak-check=full --tool=memcheck ibrun -n ${nproc} ./uwLaser ${args}

# ============================================================================================
# C. Rectangular waveguide for linear Maxwell (e.g., pollution study)
#    / with dirichlet BC (-ibc 0), or impedance BC (-ibc 2 or -ibc 3), or PML
#    / -omega sqrt(5.d0)/2.d0*PI -gamma sqrt(1.d0-PI*PI/(w*w))
#      --> 1 wavelength per 4 unit lengths in z-direction
#    / can be run with NEXACT=1 or NEXACT=0 (with PML: --use NEXACT=0)
usepml=false
ibc=0

# VARYING LENGTH OF WAVEGUIDE
# set waveguide length, #refs, maxnodes, 4 elems/wavelength
# L=4       (1 wavelength ), 3  refs (1 aniso x,  2 aniso z --> 2 *    4  elems),    250 nodes
zl=4.0d0   ; imax=3 ; maxnods=250   ; file_geometry='../GEOMETRIES/waveguide/rect_4'
# L=8       (2 wavelengths), 4  refs (1 aniso x,  3 aniso z --> 2 *    8  elems),    550 nodes
#zl=8.0d0   ; imax=4 ; maxnods=550   ; file_geometry='../GEOMETRIES/waveguide/rect_8'
# L=16      (4 wavelengths), 5  refs (1 aniso x,  4 aniso z --> 2 *   16  elems),   1050 nodes
#zl=16.0d0  ; imax=5 ; maxnods=1050  ; file_geometry='../GEOMETRIES/waveguide/rect_16'
# L=32      (8 wavelengths), 6  refs (1 aniso x,  5 aniso z --> 2 *   32  elems),   2050 nodes
#zl=32.0d0  ; imax=6 ; maxnods=2050  ; file_geometry='../GEOMETRIES/waveguide/rect_32'
# L=64     (16 wavelengths), 7  refs (1 aniso x,  6 aniso z --> 2 *   64  elems),   4050 nodes
#zl=64.0d0  ; imax=7 ; maxnods=4050  ; file_geometry='../GEOMETRIES/waveguide/rect_64'
# L=128    (32 wavelengths), 8  refs (1 aniso x,  7 aniso z --> 2 *  128  elems),   8050 nodes
#zl=128.0d0 ; imax=8 ; maxnods=8050  ; file_geometry='../GEOMETRIES/waveguide/rect_128'
# L=256    (64 wavelengths), 9  refs (1 aniso x,  8 aniso z --> 2 *  256  elems),  16050 nodes
#zl=256.0d0 ; imax=9 ; maxnods=16050 ; file_geometry='../GEOMETRIES/waveguide/rect_256'
# L=512   (128 wavelengths), 10 refs (1 aniso x,  9 aniso z --> 2 *  512  elems),  32050 nodes
#zl=512.0d0 ; imax=10; maxnods=32050 ; file_geometry='../GEOMETRIES/waveguide/rect_512'
# L=1024  (256 wavelengths), 11 refs (1 aniso x, 10 aniso z --> 2 * 1024  elems),  64050 nodes
#zl=1024.0d0; imax=11; maxnods=64050 ; file_geometry='../GEOMETRIES/waveguide/rect_1024'

# VARYING FREQUENCY (WAVE NUMBER) OF WAVEGUIDE
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
#imax=12; maxnods=128050 ; omega=804.2538552187321009202d0
# 1024 wavelengths (256.00 waves/unit length)
#imax=13; maxnods=256050 ; omega=1608.498506596624078797d0

ctrl='../COMMON_FILES/control_1'
maxnods=10000

args=" -geom 1 -isol 5 -comp 2"
#args+=" -gamma ${gamma}"
#args+=" -omega ${omega}"
args+=" -ref_core ${ref_core} -ref_clad ${ref_clad}"
args+=" -aniso_ref_index ${aniso_ref_index}"
args+=" -art_grating ${art_grating}"
args+=" -job ${job} -imax ${imax} -jmax ${jmax}"
args+=" -ibc ${ibc}"
args+=" -px ${px} -py ${py} -pz ${pz} -dp ${dp} -npx 4 -npy 4 -npz 4"
args+=" -copump 1 -nlflag ${nlflag} -gain ${gain} -raman ${raman}"
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
   args+=" -wavenum_signal 1.0d0 -wavenum_pump 1.0d0"
fi

#mpirun -np ${nproc} ./uwLaser ${args}
#ibrun -n ${nproc} ./uwLaser ${args}

# ============================================================================================
# ============================================================================================
# ============================================================================================
# =============================
# OLD RUN CONFIGS (OpenMP CODE)
#
# 1) testing convergence rates on a cube -- use NEXACT=1
#file_geometry='../GEOMETRIES/cubes/cube'
#./uwLaser -geom 1 -isol 2 -omega 1.0d0 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -npx 3 -npy 3 -npz 3 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 2) testing convergence rates on the fiber (core only) -- use NEXACT=1
#file_geometry='../GEOMETRIES/fiber/fhcor_curv'
#./uwLaser -geom 4 -isol 2 -omega 0.1d0 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc 0 -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 3) (non-)linear solve of Maxwell (with heat/gain) (signal+pump fields)
#file_geometry='../GEOMETRIES/fiber/fiber_hexa9'
#file_geometry='../GEOMETRIES/fiber/fiber_hexa13'
#file_geometry='../GEOMETRIES/fiber/fiber_prism'
#ibrun -n ${nproc} ./uwLaser -geom 5 -isol 13 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -usepml -pmlfrac ${pmlfrac} -nlflag ${nlflag} -gain ${gain} -raman ${raman} -heat ${heat} -aniso_heat ${aniso_heat} -aniso_ref_index ${aniso_ref_index} -nsteps ${nsteps} -dt ${dt} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 4) rectangular waveguide for linear Maxwell (e.g., pollution study)
#file_geometry='../GEOMETRIES/waveguide/rect_4'
#    / with dirichlet BC (-ibc 0), or impedance BC (-ibc 2 or -ibc 3), or PML
#    / -omega sqrt(5.d0)/2.d0*PI -gamma sqrt(1.d0-PI*PI/(w*w))
#      --> 1 wavelength per 4 unit lengths in z-direction
#    / can be run with NEXACT=1 or NEXACT=0 (unless PML --use NEXACT=0)
#./uwLaser -geom 1 -isol 5 -comp 2 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 5) circular hollow waveguide for linear Maxwell
#file_geometry='../GEOMETRIES/fiber/fhcor_curv'
#./uwLaser -geom 4 -isol 12 -comp 1 -bessord 1.d0 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
# 6) dielectric fiber waveguide for linear Maxwell
#file_geometry='../GEOMETRIES/fiber/fiber_prism'
#./uwLaser -geom 5 -isol 13 -comp 1 -px ${px} -py ${py} -pz ${pz} -dp ${dp} -ibc ${ibc} -usepml -pmlfrac ${pmlfrac} -nlflag 0 -heat 0 -aniso_heat ${aniso_heat} -copump 1 -zl ${zl} -nthreads ${nthreads} -dir_output ${dir_output} -vis_level ${vis_level} -file_geometry ${file_geometry}
#
#eof

