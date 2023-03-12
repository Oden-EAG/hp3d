#!/bin/bash

formulation=ultraweak

### HP3D PARAMETERS
#
#
#control file
control=./control/control
#
#physics file
physics=./input/physics
#
#history file
history=./input/history
#
#error file
err=./output/errorlogs/log.tx
#
#visualization level (0-3)
level=2
#
#PML flag
# ipml=1
# paraview directory
paraview=../output/gaussbeam/paraview/
#matlab=../output/MATLAB/
#
# problem: 0 = freespace, 1 = Cavity,      2 = sphere, 
#                         3 = scat_cavity, 4 = scat_sphere
#                         5 = scat_cavity_pml   
# prob=1
prob=1
#geometry
if [ $prob == 0 ]; then
    prob_verb="Free space"
    # geometry=./geometries/hexa_orient0
    geometry=./geometries/hexa
    # geometry=./geometries/MATLAB/hexa342_orient0_unit_cube
    # geometry=./geometries/hexa27
    # geometry=./geometries/hemisphere_orient0
    bc=3
elif
    [ $prob == 1 ]; then
    prob_verb="Cavity"
    # geometry=./geometries/hexa26    
    # geometry=./geometries/hexa26_orient0    
    # geometry=./geometries/MATLAB/hexa124_orient0_alongated    
    # geometry=./geometries/MATLAB/hexa124_orient0_unit_cube
    # geometry=./geometries/MATLAB/hexa124_orient0_alongated_unit_cube
    geometry=./geometries/MATLAB/hexa342_orient0_unit_cube
    bc=4
elif
    [ $prob == 2 ]; then
    prob_verb="sphere"
    geometry=./geometries/MATLAB/sphere_matlab
    bc=6
elif
    [ $prob == 3 ]; then
    prob_verb="Cavityscattering"
    # geometry=./geometries/MATLAB/hexa124_orient0    
    # geometry=./geometries/MATLAB/hexa124_orient0_alongated    
    # geometry=./geometries/MATLAB/hexa124_orient0_unit_cube
    geometry=./geometries/MATLAB/hexa342_orient0_unit_cube
    # geometry=./geometries/hexa26_orient0    
    # geometry=./geometries/hexa26
    bc=5
    # bc=1

    # geometry=./geometries/hexa26
elif
    [ $prob == 4 ]; then
    prob_verb="Sphere scattering"
    # geometry=./geometries/MATLAB/hemisphere
    # geometry=./geometries/MATLAB/sphere_matlab
    geometry=./geometries/MATLAB/sphere_matlab_orient0
    bc=6
    # geometry=./geometries/hexa26
elif [ $prob == 5 ]; then
    prob_verb="Cubic scattering with PML"
    # geometry=./geometries/hexa_orient0
    # geometry=./geometries/hexa
    geometry=./geometries/MATLAB/hexa342_orient0_unit_cube
    # geometry=./geometries/hexa27
    # geometry=./geometries/hemisphere_orient0
    bc=1
    # bc=2
elif
    [ $prob == 6 ]; then
    prob_verb="Sphere scattering WITH PML"
    # geometry=./geometries/MATLAB/hemisphere
    # geometry=./geometries/MATLAB/sphere_matlab
    geometry=./geometries/MATLAB/sphere_matlab_orient0
    bc=1
    # geometry=./geometries/hexa26    
elif [ $prob == 7 ]; then
    prob_verb="Thin plane scattering with PML"
    geometry=./geometries/MATLAB/hexa64_orient0_unit_cube
    bc=3
fi

#
#boundary condition (1=Dirichlet, 2=Neumann, 3=Impedance)
# bc=5
if [ $bc == 1 ]; then
    bc_verb="Dirichlet"
elif [ $bc == 2 ]; then
    bc_verb="Neumann"
elif [ $bc == 3 ]; then
    bc_verb="Impedance"
elif [ $bc == 4 ]; then
    bc_verb="bc_cavity"
elif [ $bc == 5 ]; then
    bc_verb="bc_cavity_scattering"    
elif [ $bc == 6 ]; then
    bc_verb="bc_sphere_scattering"    
fi
#
# Polynomial order
p=1
#
# Enrichment order
dp=1
#
# Test norm (1 = Adj graph, 2 = Mathematician's)
testNorm=1
if [ $testNorm == 1 ]; then
    testnorm_verb="Adjoint Graph norm"
else
    testnorm_verb="Mathematician's norm"
fi
#
# Exact solution
exact=5
if [ $exact == 0 ]; then
    exact_verb="Polynomial in x dimension"
elif [ $exact == 1 ]; then
    exact_verb="Sine solution"
elif [ $exact == 2 ]; then
    exact_verb="Polynomial solution vanishing on the boundary"
elif [ $exact == 3 ]; then
    exact_verb="Plane wave"
elif [ $exact == 4 ]; then
    exact_verb="Point source wave"
elif [ $exact == 5 ]; then
    exact_verb="Gaussian beam"
elif [ $exact == 6 ]; then
    exact_verb="Gaussian beam simplified"       
elif [ $exact == 7 ]; then
    exact_verb="Plane wave in z vanishing on x y boundary"          
fi
#
#number of wavelengths
rnum=5.0d0
#
#Test norm scaling parameter
alpha=1.0d0
#
#solver (1=Multigrid, 2=Smoother)
solver=1
if [ $solver == 1 ]; then
    solver_verb="Multigrid"
else
    solver_verb="Block Jacobi"
fi
#
#OMP THREADS
    
if [ $HOSTNAME == typhoon.ices.utexas.edu ]; then
    nthreads=28
elif [ $HOSTNAME == ritz.ices.utexas.edu ]; then
    nthreads=24
fi    
    nthreads=1
    
# Executable
echo
if [ $formulation == ultraweak ]; then
    exec=DPG_uweak_acoustics
    echo "Running ultraweak formulation..."
    echo "Problem          = " $prob_verb
    echo "geometry file    = " $geometry
    echo "Order            = " $p
    echo "Enriched order   = " $dp
    echo "Test norm        = " $testnorm_verb
    echo "bc               = " $bc_verb
    echo "Exact solution   = " $exact_verb
    echo "# of wavelengths = " $rnum
    echo "solver           = " $solver_verb
    echo "OMP threads      = " $nthreads
fi
echo
# read -p "Press [Enter] key to continue..."

mpirun -np 2 ./$exec \
    -file-geometry $geometry\
    -file-control $control \
    -file-phys $physics \
    -file-history $history \
    -file-err $err \
    -vis-level $level \
    -dir-paraview $paraview \
    -dir-matlab $matlab \
    -prob $prob \
    -p $p \
    -dp $dp \
    -norm-test $testNorm \
    -exact $exact \
    -bc $bc \
    -rnum $rnum   \
    -alpha $alpha \
    -solver $solver \
    -nthreads $nthreads

