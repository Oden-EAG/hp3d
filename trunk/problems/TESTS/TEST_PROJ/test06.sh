#!/bin/bash

# set order
P=3
echo "Testing tet orientations (p = $P) ..."

./proj $* \
    -file-control   ./files/tetr/control \
    -file-geometry  ./files/tetr/geom \
    -file-phys      ./files/tetr/physics \
    -file-err       ./files/dump_error06 \
    -file-tags      ./files/dump_tags06 \
    -vis-level      '4' \
    -paraview-attr  \
    -prefix         'tetr_' \
    -problem        'tetr' \
    -solver         1 \
    -testing        2 \
    -order-tetr     $P \
    -quiet-mode

# check error dump
./check_err.sh dump_error06
