#!/bin/bash

prog=../build/bin/test_bounded_dist
source_mesh=../dat/kitty/orig.tet.vtk
initial_mesh=../dat/kitty/polycube.tet.vtk
pos_constraint=../dat/kitty/constraint.fv
output_folder=../build/bin/large_scale_bounded_distortion
K=1.2
spectral_radius=0.973

if [ ! -d "$output_folder" ]; then
  mkdir -p $output_folder
fi
time $prog -s $source_mesh -i $initial_mesh -c $pos_constraint -o $output_folder --method $1 -k $K -e 1e-9 --spectral_radius $spectral_radius