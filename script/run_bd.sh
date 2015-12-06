#!/bin/bash

prog=../build/bin/test_bounded_dist
source_mesh=../dat/kitty/orig.tet.vtk
initial_mesh=../dat/kitty/polycube.tet.vtk
pos_constraint=../dat/kitty/constraint.fv
output_folder=../build/bin/BD/kitty
K=2

time $prog -s $source_mesh -i $initial_mesh -c $pos_constraint -o $output_folder -k $K