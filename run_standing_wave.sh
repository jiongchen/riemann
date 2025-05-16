#!/bin/bash

exe=build/examples/test_wave_construction

mesh=dat/plane_re_sub2.obj
frame=dat/plane_re_sub2.ef
feature=dat/plane_re_sub2.fl

${exe} ${mesh} ${frame} ${feature}
