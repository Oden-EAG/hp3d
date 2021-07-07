#!/bin/bash

echo "Testing random refinements ..."

./proj $* \
    -file-control   ./files/ref/control \
    -file-geometry  ./files/ref/geom \
    -file-phys      ./files/ref/physics \
    -vis-level      '0' \
    -paraview-geom  \
    -prefix         'ref_' \
    -problem        'hybr' \
    -testing        4 \
    -order-tetr     1 \
    -order-bric     1 \
    -order-pris     1 \
    -order-pyra     1 \
    -percentage     0.05 \
    -ref-levels     20 \
    -quiet-mode
