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

# Set polynomial order p
p=6

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=0
imax=3

# max NODES
maxnods=10000000

# export KMP_STACKSIZE=24M   # p=3
# export KMP_STACKSIZE=32M   # p=4
# export KMP_STACKSIZE=48M   # p=5
export KMP_STACKSIZE=64M   # p=6
#export KMP_STACKSIZE=80M   # p=7
#export KMP_STACKSIZE=96M   # p=8

# Fix Intel18 OMP issue
export KMP_INIT_AT_FORK=FALSE

# component number for manufactured solution
comp=1
# solution number
isol=102
# DPG test norm scaling
alpha=100.0
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
k=217.48
# bending radius
rbend=1300.0
# coordinate bounds
thup=0.1553118594095
thlo=0.0
rup=1306.25
rlo=1293.75
xup=-6.25
xlo=6.25
rhoup=6.25
rholo=0.0
# pml proportions
pmlthup=0.125
pmlthlo=0.0
pmlrup=0.0
pmlrlo=0.0
pmlxup=0.0
pmlxlo=0.0
pmlrhoup=0.2
pmlrholo=0.0
# flags for toroidal pml and slab guide geometry
torpml=1
slab=0
# fiber radii (or slab halfwidths)
rcore=0.5
rclad=5.0
rcoat=6.25
# refractive indices
ncore=1.4512
nclad=1.45
ncoat=1.38
# coating attenuation coefficient (ratio/[unit length] NOT in dB...)
attncoat=1.75456984086146E-07
# set BC flag -> 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=0

# mode, bending radius, max theta, prefix for paraview and mesh file
isol=101
rbend=1300.0
thup=0.10713716101510534
pref='fiber_R1300_4wl_lp01'
file_geometry='./geometries/fiber_R1300_4wl_lp01'

# isol=102
# rbend=1300.0
# thup=0.3835459596931083
# pref='fiber_R1300_4wl_lp02'
# file_geometry='./geometries/fiber_R1300_4wl_lp02'

# isol=101
# rbend=2600.0
# rup=2606.25
# rlo=2593.75
# thup=0.05356858050755267
# pref='fiber_R2600_4wl_lp01'
# file_geometry='./geometries/fiber_R2600_4wl_lp01'

# isol=102
# rbend=2600.0
# rup=2606.25
# rlo=2593.75
# thup=0.3804840429151485
# pref='fiber_R2600_4wl_lp02'
# file_geometry='./geometries/fiber_R2600_4wl_lp02'

# pref='newtest'

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

### RUNNING ON PERSONAL LAPTOP
# Set MPI Procs and OpenMP threads
nproc=4
nthreads=1
mpirun -np ${nproc} ./bending ${args} -nthreads ${nthreads}
###

# ### RUNNING ON TACC'S FRONTERA'S CLX NODE (TACC)
# # Set MPI Procs and OpenMP threads
# nproc=32
# nthreads=56
# ibrun -n ${nproc} ./bending ${args} -nthreads ${nthreads}
# ###