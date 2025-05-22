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
p=3

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=0
imax=3

# max NODES
maxnods=123456

export KMP_STACKSIZE=24M   # p=3
# export KMP_STACKSIZE=32M   # p=4
# export KMP_STACKSIZE=48M   # p=5
# export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8

# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# component number for manufactured solution
comp=1
# solution number
isol=200
# DPG test norm scaling
alpha=0.01
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
k=217.455
# bending radius
rbend=1300.0
# coordinate bounds
thup=0.1553118594095
thlo=0.0
rup=1305.0
rlo=1295.0
xup=0.5
xlo=-0.5
rhoup=0.0
rholo=0.0
# pml proportions
pmlthup=0.5
pmlthlo=0.0
pmlrup=0.225
pmlrlo=0.225
pmlxup=0.0
pmlxlo=0.0
pmlrhoup=0.0
pmlrholo=0.0
# flags for toroidal pml and slab guide geometry
torpml=0
slab=1
# fiber radii (or slab halfwidths)
rcore=0.5
rclad=5.0
rcoat=10.0
# refractive indices
ncore=1.4512
nclad=1.45
ncoat=1.38
# coating attenuation coefficient (ratio/[unit length] NOT in dB...)
attncoat=1.75456984086146E-07
# set BC flag -> 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=0

# prefix for paraview
# pref='bslab_Dir1_EnvAns1_short' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir1_EnvAns0_short' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir2_EnvAns1' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='bslab_Dir2_EnvAns0' # with file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# pref='fiber_15deg_lp01'
# pref='fiber_2deg_lp01'
pref='stepslab_12wl_m00'
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
# file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# file_geometry='./geometries/bent_fiber_test_15'
# file_geometry='./geometries/partly_bent_fiber_test_4wl'
file_geometry='./geometries/bent_stepslab_0p5_5_1300_12wl'

ctrl='control/control'
#
args="  -file_control ${ctrl} -file_geometry ${file_geometry}"
args+=" -p ${p} -dp ${dp}"
args+=" -comp ${comp} -isol ${isol} -imax ${imax} -job ${job}"
args+=" -maxnods ${maxnods} -alpha ${alpha}"
args+=" -mu ${mu} -epsilon ${epsilon} -e0 ${e0} -h0 ${h0}"
args+=" -omega ${omega} -gamma ${gamma}"
args+=" -k ${k} -rbend ${rbend}"
args+=" -thup ${thup} -thlo ${thlo}"
args+=" -rup ${rup} -rlo ${rlo}"
args+=" -xup ${xup} -xlo ${xlo}"
args+=" -rhoup ${rhoup} -rholo ${rholo}"
args+=" -pmlthup ${pmlthup} -pmlthlo ${pmlthlo}"
args+=" -pmlrup ${pmlrup} -pmlrlo ${pmlrlo}"
args+=" -pmlxup ${pmlxup} -pmlthlo ${pmlxlo}"
args+=" -pmlrhoup ${pmlrhoup} -pmlrholo ${pmlrholo}"
args+=" -torpml ${torpml} -slab ${slab}" 
args+=" -rcore ${rcore} -rclad ${rclad} -rcoat ${rcoat}"
args+=" -ncore ${ncore} -nclad ${nclad} -ncoat ${ncoat} -attncoat ${attncoat}"
args+=" -ibc ${ibc}"
args+=" -prefix ${pref} -vis_level ${vis_level} -dir_output ${dir_output}"
args+=" -nthreads ${nthreads}"

mpirun -np ${nproc} ./bending ${args}
