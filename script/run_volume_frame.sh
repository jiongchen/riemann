#!/bin/bash

EXE=../build/bin/test_vol_frame
MESH=../dat/cube.tet
OUT_DIR=../result/volume_frame/cube

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

# $EXE -i $MESH -o $OUT_DIR


# Visualize sigularity of frame field
cd ../bin
FF_EXE=./frame_field 
TET_FILE=../dat/sphere265k_tet.vtk
ZYZ_FILE=../dat/sphere265k_f1e3_s1e5_n1000.zyz.txt
OUT_FILE=$OUT_DIR/sing.vtk

$FF_EXE prog=draw_3d_frame_sing tet=$TET_FILE zyz=$ZYZ_FILE out=$OUT_FILE
