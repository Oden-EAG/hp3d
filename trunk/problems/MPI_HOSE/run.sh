#!/bin/bash

# formulation=Projection
# formulation=Primal
formulation=Hybrid

### HP3D PARAMETERS

# Polynomial order
p=2

# Enrichment order
dp=1

# Test norm (1 = Adj graph, 2 = Mathematician's)
testNorm=2

# Trial norm (1 = L2, 2 = Mathematician's)
trialNorm=1

# Exact solution
exact=5
control=./common/control_sheathe
# control=./common/control_sheathe_exact
# control=./common/control_sheathe_test

# Boundary condition
bc=3

### UHM PARAMETERS
# dataType=
NThreads=44
blockSize=128
numRHS=1
randomMatrixType=1
treeCuttingLevel=0
overAllocationLevel=100


# Paraview file prefix
if [ $# = 0 ]; then
    prefix=Sheathe_
else
    prefix=$1
fi

# Executable
echo
if [ $formulation == Primal ]; then
# if [ -z $formulation ]; then
    elast=DPG_incomp_primal
    physics=INCOMP_PRIMAL/input/physics
    history=INCOMP_PRIMAL/input/history
    err=INCOMP_PRIMAL/output/errorlogs/log.txt
    paraview=INCOMP_PRIMAL/output/figures/
    echo "Running constrained primal formulation..."
elif [ $formulation == PrimalImproved ]; then
    elast=DPG_impr_primal
    physics=IMPROVED_PRIMAL/input/physics
    history=IMPROVED_PRIMAL/input/history
    err=IMPROVED_PRIMAL/output/errorlogs/log.txt
    paraview=IMPROVED_PRIMAL/output/figures/
    echo "Running improved constrained primal formulation..."
elif [ $formulation == Hybrid ]; then
    elast=DPG_hybrid
    physics=HYBRID/input/physics
    history=HYBRID/input/history
    err=HYBRID/output/errorlogs/log.txt
    paraview=HYBRID/output/figures/
    echo "Running hybrid formulation..."
elif [ $formulation == Weak_Symm ]; then
    elast=DPG_uweak
    physics=WEAK_SYMMETRY/input/physics
    history=WEAK_SYMMETRY/input/history
    err=WEAK_SYMMETRY/output/errorlogs/log.txt
    paraview=WEAK_SYMMETRY/output/figures/
    echo "Running ultraweak with weak symmetry..."
elif [ $formulation == Strong_Symm ]; then
    elast=DPG_uweak_symm
    physics=STRONG_SYMMETRY/input/physics
    history=STRONG_SYMMETRY/input/history
    err=STRONG_SYMMETRY/output/errorlogs/log.txt
    paraview=STRONG_SYMMETRY/output/figures/
    echo "Running ultraweak with strong symmetry..."
elif [ $formulation == Projection ]; then
    elast=proj
    physics=PROJECTION/input/physics
    history=PROJECTION/input/history
    err=PROJECTION/output/errorlogs/log.txt
    paraview=PROJECTION/output/figures/
    echo "Running projection..."
fi
echo

# L-shape domain problem
# gdb --args \
./$elast \
    -file-geometry ../common/geometries/other/sheathed_tube_cyl \
    -file-control $control \
    -file-phys $physics \
    -file-history $history \
    -file-err $err \
    -dir-paraview $paraview \
    -p $p \
    -dp $dp \
    -norm-test $testNorm \
    -norm-trial $trialNorm \
    -exact $exact \
    -bc $bc \
    -prefix $prefix \
    -blk $blockSize \
    -rhs $numRHS \
    -spd $randomMatrixType \
    -level $treeCuttingLevel \
    -allc $overAllocationLevel \
    -nthreads $NThreads

