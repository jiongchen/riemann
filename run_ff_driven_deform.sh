#!/bin/bash


exe=build/examples/test_ff_deform
mesh=dat/surf_sphere.obj
ff_file=dat/surf_sphere_cons.ff

${exe} ${mesh} ${ff_file}

