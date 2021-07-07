#!/bin/bash

# set order
P=4
echo "Testing orientations on hybrid mesh (p = $P) ..."

./proj $* \
    -file-control   ./files/hybr/control \
    -file-geometry  ./files/hybr/geom \
    -file-phys      ./files/hybr/physics \
    -file-err       ./files/dump_error08 \
    -file-tags      ./files/dump_tags08 \
    -vis-level      '4' \
    -paraview-attr  \
    -prefix         'hybr_' \
    -problem        'hybr' \
    -solver         1 \
    -testing        2 \
    -order-tetr     $P \
    -order-bric     $P \
    -order-pris     $P \
    -order-pyra     $P \
    -quiet-mode

# check error dump
./check_err.sh dump_error08
