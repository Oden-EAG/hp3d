#!/bin/bash

# set order
P=3
echo "Testing hexa orientations (p = $P) ..."

./proj $* \
    -file-control   ./files/hexa/control \
    -file-geometry  ./files/hexa/geom \
    -file-phys      ./files/hexa/physics \
    -file-err       ./files/dump_error04 \
    -file-tags      ./files/dump_tags04 \
    -vis-level      '3' \
    -paraview-attr  \
    -prefix         'hexa_' \
    -problem        'hexa' \
    -solver         1 \
    -testing        2 \
    -order-bric     $P \
    -quiet-mode

# check error dump
./check_err.sh dump_error04
