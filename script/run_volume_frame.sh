#!/bin/bash

echo -e "# ================================== Run cross frame fild ========================================= #"
EXE=../build/bin/test_vol_frame
MESH=../dat/tets/fandisk.c58k.tet.vtk
OUT_DIR=../result/volume_frame/fandisk
WS=1e0
WA=1e3
EPS=1e-5
MAXITS=1000

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

$EXE -i $MESH -t vtk -o $OUT_DIR --ws=$WS --wa=$WA --epsf=$EPS --maxits=$MAXITS | tee LOG.txt


echo -e "# ============================== Visualize sigularity of frame field ============================== #"
cd ../bin

FF_EXE=./frame_field
MESH_FILE=$OUT_DIR/tet.vtk
ZYZ_FILE=$OUT_DIR/zyz.txt
SING_FILE=$OUT_DIR/sing.vtk

$FF_EXE prog=draw_3d_frame_sing tet=$MESH_FILE zyz=$ZYZ_FILE out=$SING_FILE

echo -e "# ================================================================================================= #"
