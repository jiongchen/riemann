#!/bin/bash

exe=build/examples/test_wave_construction

mesh=dat/plane_re.obj
frame=dat/plane_re.ef
feature=dat/plane_re.fl

mesh=dat/wave_special_mesh.obj
frame=dat/wave_special_mesh.ef
feature=dat/wave_special_mesh.fl

${exe} ${mesh} ${frame} ${feature}

mv standing_wave result/
