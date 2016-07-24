#!/bin/bash

echo -e "# ================================== Run cross frame field ======================================== #"

EXE=../build/bin/compare_cubic_smooth
MESH=../dat/tets/sphere.c13k.vtk

TYPE=SH

WS=1e0
WP=1e3

ABS=0.1
WO=1e3

ESPF=1e-8
MAXITS=1000

OUT_DIR=../result/compare_smooth/sphere-abs$ABS-ws$WS-wo$WO-wp$WP-tp$TYPE

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

if [ "$TYPE" == "SH" ]; then
    $EXE mesh=$MESH out_dir=$OUT_DIR sm_type=$TYPE weight.smooth=$WS weight.boundary=$WP lbfgs.maxits=$MAXITS lbfgs.epsf=$ESPF | tee $OUT_DIR/log.txt
elif [ "$TYPE" == "L1" ]; then
    $EXE mesh=$MESH out_dir=$OUT_DIR sm_type=$TYPE weight.smooth=$WS weight.orth=$WO weight.boundary=$WP abs_eps=$ABS lbfgs.maxits=$MAXITS lbfgs.epsf=$ESPF | tee $OUT_DIR/log.txt
else
    echo -e "{error return}"
fi


echo -e "# ============================== Visualize sigularity of frame field ============================== #"

cd ../bin
export LD_LIBRARY_PATH=./3rd:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./sys:$LD_LIBRARY_PATH

FF_EXE=./frame_field
MESH_FILE=$OUT_DIR/tet.vtk
FRAME_FILE=$OUT_DIR/frames.txt
SING_FILE=$OUT_DIR/sing.vtk

if [ "$TYPE" == "SH" ]; then
    $FF_EXE prog=draw_3d_frame_sing tet=$MESH_FILE zyz=$FRAME_FILE out=$SING_FILE
elif [ "$TYPE" == "L1" ]; then
    $FF_EXE prog=draw_3d_frame_sing tet=$MESH_FILE ff=$FRAME_FILE out=$SING_FILE
else
    echo -e "{error} return"
fi

echo -e "# ================================================================================================= #"
