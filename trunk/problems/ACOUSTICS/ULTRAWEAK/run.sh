#!/bin/bash

formulation=ultraweak

### HP3D PARAMETERS

#
#geometry
geometry=./geometries/hexa
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
# paraview directory
paraview=./output/figures/
#
# problem: 0 = freespace, 1=Cavity
prob=0
if [ $prob == 0 ]; then
    prob_verb="Free space"
else
    prob_verb="Cavity"
fi
#
# Polynomial order
p=2
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
fi
#
#boundary condition (1=Dirichlet, 2=Neumann, 3=Impedance)
bc=3
if [ $bc == 1 ]; then
    bc_verb="Dirichlet"
elif [ $bc == 2 ]; then
    bc_verb="Neumann"
elif [ $bc == 3 ]; then
    bc_verb="Impedance"
fi
#
#number of wavelengths
rnum=20.3d0
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
nthreads=48

# # Paraview file prefix
# if [ $# = 0 ]; then
#     prefix=Sheathe_
# else
#     prefix=$1
# fi

# Executable
echo
if [ $formulation == ultraweak ]; then
    exec=DPG_uweak_acoustics
    echo "Running ultraweak formulation..."
    echo "Problem          = " $prob_verb
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
read -p "Press [Enter] key to continue..."

# L-shape domain problem
# gdb --args \
./$exec \
    -file-geometry $geometry\
    -file-control $control \
    -file-phys $physics \
    -file-history $history \
    -file-err $err \
    -vis-level $level \
    -dir-paraview $paraview \
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