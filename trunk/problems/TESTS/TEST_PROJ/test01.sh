#!/bin/bash

# set 
echo "Testing prism shape functions ..."

./proj $* \
    -file-control   ./files/pris/control \
    -file-geometry  ./files/pris/geom \
    -file-phys      ./files/pris/physics \
    -file-err       ./files/dump_error01 \
    -file-tags      ./files/dump_tags01 \
    -vis-level      '4' \
    -paraview-attr  \
    -prefix         'pris_' \
    -problem        'pris' \
    -solver         1 \
    -testing        1 \
    -order-pris     6 \
    -quiet-mode

# check error dump
./check_err.sh dump_error01
