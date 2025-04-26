#!/bin/bash

exe=build/examples/test_spin_trans

# data 1
mesh=dat/spin_sphere.obj
curvature_change=dat/sphere_curv_change.txt
output_folder=result/spin_transform/sphere/

# data 2
mesh=dat/capsule.obj
curvature_change=dat/capsule_curv_change.txt
output_folder=result/spin_transform/capsule/

${exe} -i ${mesh} -c ${curvature_change} -o ${output_folder}
