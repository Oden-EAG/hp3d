#!/bin/bash

# set order
P=3
echo "Testing prism orientations (p = $P) ..."

./proj $* \
    -file-control   ./files/pris/control \
    -file-geometry  ./files/pris/geom \
    -file-phys      ./files/pris/physics \
    -file-err       ./files/dump_error02 \
    -file-tags      ./files/dump_tags02 \
    -vis-level      '4' \
    -paraview-attr  \
    -prefix         'pris_' \
    -problem        'pris' \
    -solver         1 \
    -testing        2 \
    -order-pris     $P \
    -quiet-mode

# check error dump
./check_err.sh dump_error02
