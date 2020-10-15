#!/bin/bash

# polynomial order
p=3

# error calculation (1= ...)
error=3

# exact solution (1 = ...)
exact=5

if [ $# = 0 ]; then
    elast=BG_primal
    echo "Running Bubnov Galerkin primal problem..."
    echo
else
    elast=$1
    echo "Running $elast"
    echo
fi

# L-shape domain problem
./$elast \
    -p $p \
    -error $error \
    -exact $exact \
    -bc 8 \
    -file-geometry ./files/hexa/L_shape \
    -file-control ./files/control/control_L_shape
