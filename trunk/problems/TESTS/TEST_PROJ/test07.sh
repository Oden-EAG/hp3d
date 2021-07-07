#!/bin/bash

# set mesh order
P=4

echo "Testing constrained approximation on hybrid mesh (p = $P) ..."

./proj $* \
    -file-control   ./files/hybr/control \
    -file-geometry  ./files/hybr/geom \
    -file-phys      ./files/hybr/physics \
    -file-err       ./files/dump_error07 \
    -file-tags      ./files/dump_tags07 \
    -problem        'hybr' \
    -solver         1 \
    -testing        3 \
    -order-tetr     $P \
    -order-bric     $P \
    -order-pris     $P \
    -order-pyra     $P \
    -quiet-mode

# check error dump
./check_err.sh dump_error07
