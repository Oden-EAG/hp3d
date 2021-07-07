#!/bin/bash

echo "Testing tet shape functions ..."

./proj $* \
    -file-control   ./files/tetr/control \
    -file-geometry  ./files/tetr/geom \
    -file-phys      ./files/tetr/physics \
    -file-err       ./files/dump_error05 \
    -file-tags      ./files/dump_tags05 \
    -vis-level      '4' \
    -paraview-attr  \
    -prefix         'tetr_' \
    -problem        'tetr' \
    -solver         1 \
    -testing        1 \
    -order-tetr     6 \
    -quiet-mode

# check error dump
./check_err.sh dump_error05
