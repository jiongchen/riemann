#!/bin/bash

exe=build/examples/test_deform_transfer

source_mesh=dat/dt_horse_ref.obj
target_mesh=dat/dt_camel_ref.obj
marker=dat/dt_horse_camel.cons
source_defo=dat/dt_horse_01.obj

${exe} ${source_mesh} ${target_mesh} ${marker} ${source_defo}
