#!/bin/bash

EXE=../build/bin/test_vol_frame
MESH=../dat/cube.tet
OUT_DIR=../result/volume_frame/cube

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

$EXE -i $MESH -o $OUT_DIR | tee LOG.txt


# # === Visualize sigularity of frame field === #
# cd ../bin

# FF_EXE=./frame_field
# MESH_FILE=$OUT_DIR/tet.vtk
# ZYZ_FILE=$OUT_DIR/zyz.txt
# SING_FILE=$OUT_DIR/sing.vtk

# $FF_EXE prog=draw_3d_frame_sing tet=$MESH_FILE zyz=$ZYZ_FILE out=$SING_FILE
