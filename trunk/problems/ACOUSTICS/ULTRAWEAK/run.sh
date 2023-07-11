#!/bin/bash

### HP3D PARAMETERS
# Control file
control=./control/control

# Physics file
physics=./input/physics

# History file
history=./input/history

# Error file
err=./output/errorlogs/log.tx

# Paraview visualization level (0-3)
level=1

# paraview directory
paraview=../output/paraview/

# Polynomial order
p=2

# Enrichment order
dp=1

# Number of OpenMP threads
nthreads=1

# Problem: 0 = freespace, 1 = Cavity
prob=1

# Geometry
if [ $prob == 0 ]; then
    prob_verb="Free space"
    geometry=./geometries/hexa27
    bc=3
elif
    [ $prob == 1 ]; then
    prob_verb="Cavity"
    geometry=./geometries/hexa342_cavity
    bc=4
fi

bc=1
#boundary condition (1=Dirichlet, 2=Neumann, 3=Impedance)
if [ $bc == 1 ]; then
    bc_verb="Dirichlet"
elif [ $bc == 2 ]; then
    bc_verb="Neumann"
elif [ $bc == 3 ]; then
    bc_verb="Impedance"
elif [ $bc == 4 ]; then
    bc_verb="bc_cavity"
fi

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
fi

# Number of wavelengths
rnum=5.0d0

# DPG test norm scaling parameter
alpha=1.0d0
    
# Executable
echo

# Print setup info
exec=DPG_uweak_acoustics
echo "Running ultraweak formulation..."
echo "Problem          = " $prob_verb
echo "geometry file    = " $geometry
echo "Order            = " $p
echo "Enriched order   = " $dp
echo "bc               = " $bc_verb
echo "Exact solution   = " $exact_verb
echo "# of wavelengths = " $rnum
echo "OMP threads      = " $nthreads
echo

mpirun -np 4 ./$exec \
    -file-geometry $geometry \
    -file-control $control \
    -file-phys $physics \
    -file-history $history \
    -file-err $err \
    -vis-level $level \
    -dir-paraview $paraview \
    -prob $prob \
    -p $p \
    -dp $dp \
    -exact $exact \
    -bc $bc \
    -rnum $rnum   \
    -alpha $alpha \
    -nthreads $nthreads


