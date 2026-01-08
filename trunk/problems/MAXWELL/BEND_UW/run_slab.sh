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
p=5

# Set enriched order (p+dp)
dp=1

# EXEC JOB options
job=0
imax=3

# max NODES
maxnods=1000000

# export KMP_STACKSIZE=24M   # p=3
export KMP_STACKSIZE=32M   # p=4
# export KMP_STACKSIZE=48M   # p=5
# export KMP_STACKSIZE=64M   # p=6
# export KMP_STACKSIZE=80M   # p=7
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
gamma=1.55197558886253
# base electric field
e0=9.66588353875592E-04
# base magnetic field
h0=1.03456657220257E-03
# envelope wavenumber
k=217.43
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
pmlthup=0.25
pmlthlo=0.0
pmlrup=0.225
pmlrlo=0.225
pmlxup=0.0
pmlxlo=0.0
pmlrhoup=0.0
pmlrholo=0.0
# flags for toroidal pml and slab guide geometry
tordmn=0
slab=1
# fiber radii (or slab halfwidths)
rcore=0.5
rclad=5.0
rcoat=5.0
# refractive indices
ncore=1.4512
nclad=1.45
ncoat=1.45
nair=1.0 #   actual nair=1.00026897
# coating attenuation coefficient (ratio/[unit length] NOT in dB...)
attncoat=0.0
# elasto-optical effect flag -> 0: off, 1: on
elast=0
# set BC flag -> 0: dirichlet, 2: impedance via penalty term, 3: impedance via elimination
ibc=3


# mode, bending radius, max theta, prefix for paraview and mesh file
# isol=200
# rbend=1300.0
# thup=0.0965092161439455
# pmlthup=0.0
# file_geometry='./geometries/bentstepslab_R1300_4wl_m00'
# pref='ref_stepslab_R1300_4wl_m00'
# ibc=0

# isol=202
# rbend=1300.0
# thup=0.17043375398290803
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R1300'   # ~8 wavelengths of this straight waveguide mode for this mesh and k_env=217.43
# ibc=0
# pref='verif5_stepslab_R1300_m02_straight'


# isol=200
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# thup=0.08521687699145401
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R2600'   # ~8 wavelengths of this straight waveguide mode for this mesh and k_env=217.43
# ibc=0
# pref='verif4_stepslab_R2600_m00_straight'


# isol=200
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# thup=0.04825460807197276
# file_geometry='./geometries/bentstepslab_R2600_4wl_m00'
# pref='stepslab_R2600_4wl_m00'

# isol=202
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# pmlrlo=0.0
# thup=0.15238724408876186
# file_geometry='./geometries/bentstepslab_R2600_4wl_m02'
# pref='stepslab_R2600_4wl_m02'

# isol=11
# rbend=1300.0
# file_geometry='./geometries/bent_waveguide_d0p5_cR1300_DEG3'
# k=149.76815150765 # with k0 = 149.993333460866, mode 2: nu = 149.86045924807476 * 1300
# pmlthup=0.25
# pmlrlo=0.0
# pmlrup=0.0
# rlo=1299.5
# rup=1300.5
# ibc=0
# thup=0.05235987755983 # 3 degrees             0.008726646259971648 # 0.5 degrees
# pref='pec_r1300_p4_mode3'
# ncore=1.0
# nclad=1.0
# ncoat=1.0
# rclad=0.5
# rcoat=0.5
# gamma=1.07032799231899E+00


# # cases bent slab waveguide with semianalytical solutions obtained with Bessel-Frobenius PML approach
isol=300
rbend=1300.0
thup=0.17043375398290803
pmlrlo=0.0
file_geometry='./geometries/verif_bentstepslab_R1300'   # 7.95 wavelengths of this mode for this mesh and k_env=217.43
# the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
pref='verif5_stepslab_R1300_m00'
ibc=0

# isol=301
# rbend=1300.0
# thup=0.17043375398290803
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R1300'   # 6.40 wavelengths of this mode for this mesh and k_env=217.43
# # the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
# pref='verif5_stepslab_R1300_m01'
# ibc=0

# isol=302
# rbend=1300.0
# thup=0.17043375398290803
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R1300'   # 4 wavelengths of this mode for this mesh and k_env=217.43
# # the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
# pref='verif4_stepslab_R1300_m02'
# ibc=0

# isol=300
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# thup=0.08521687699145401
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R2600'   # 7.95 wavelengths of this mode for this mesh and k_env=217.43
# # the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
# pref='verif4_stepslab_R2600_m00'
# ibc=0

# isol=301
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# thup=0.08521687699145401
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R2600'     # 6.40 wavelengths of this mode for this mesh and k_env=217.43
# # the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
# pref='verif4_stepslab_R2600_m01'
# ibc=0

# isol=302
# rbend=2600.0
# rup=2605.0
# rlo=2595.0
# thup=0.08521687699145401
# pmlthup=0.25
# pmlrlo=0.0
# file_geometry='./geometries/verif_bentstepslab_R2600'      # 4 wavelengths of this mode for this mesh and k_env=217.43
# # the difference in this mesh file is only the bc flags: 1 on the inner face, 0 elsewhere
# pref='verif3_stepslab_R2600_m02'
# ibc=0


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
args+=" -tordmn ${tordmn} -slab ${slab}" 
args+=" -rcore ${rcore} -rclad ${rclad} -rcoat ${rcoat}"
args+=" -ncore ${ncore} -nclad ${nclad} -ncoat ${ncoat} -attncoat ${attncoat}"
args+=" -elast ${elast}"
args+=" -ibc ${ibc}"
args+=" -prefix ${pref} -vis_level ${vis_level} -dir_output ${dir_output}"

### RUNNING ON PERSONAL LAPTOP
# # Set MPI Procs and OpenMP threads
nproc=10
nthreads=1
mpirun -np ${nproc} ./bending ${args} -nthreads ${nthreads}
# ###

# ### RUNNING ON TACC'S FRONTERA'S CLX NODE (TACC)
# # Set MPI Procs and OpenMP threads
# nproc=24
# nthreads=56
# ibrun -n ${nproc} ./bending ${args} -nthreads ${nthreads}
###