#!/bin/bash

echo "Testing hexa shape functions ..."

./proj $* \
    -file-control   ./files/hexa/control \
    -file-geometry  ./files/hexa/geom \
    -file-phys      ./files/hexa/physics \
    -file-err       ./files/dump_error03 \
    -file-tags      ./files/dump_tags03 \
    -vis-level      '3' \
    -paraview-attr  \
    -prefix         'hexa_' \
    -problem        'hexa' \
    -solver         1 \
    -testing        1 \
    -order-bric     6 \
    -quiet-mode

# check error dump
./check_err.sh dump_error03
