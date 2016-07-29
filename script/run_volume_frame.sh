#!/bin/bash

echo -e "# ================================== Run cross frame field ======================================== #"

EXE=../build/bin/test_vol_frame
MESH=../dat/tets/eight.vtk
WS=1e0
WA=$1
EPS=1e-8
MAXITS=3000
OUT_DIR=../result/volume_frame/eight-ws$WS-wa$WA

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

$EXE mesh=$MESH out_dir=$OUT_DIR weight.smooth=$WS weight.align=$WA lbfgs.epsf=$EPS lbfgs.maxits=$MAXITS | tee $OUT_DIR/log.txt


echo -e "# ============================== Visualize sigularity of frame field ============================== #"

cd ../bin
export LD_LIBRARY_PATH=./3rd:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./sys:$LD_LIBRARY_PATH

FF_EXE=./frame_field
MESH_FILE=$OUT_DIR/tet.vtk
ZYZ_FILE=$OUT_DIR/zyz.txt
SING_FILE=$OUT_DIR/sing.vtk

$FF_EXE prog=draw_3d_frame_sing tet=$MESH_FILE zyz=$ZYZ_FILE out=$SING_FILE

echo -e "# ================================================================================================= #"
