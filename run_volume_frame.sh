#!/bin/bash

echo -e "# ================================== Run cross frame field ======================================== #"

EXE=build/examples/test_vol_frame
MESH=dat/tets/eight.vtk
NAME=${MESH##*/}
WS=1e0
WA=1e3
EPS=1e-8
MAXITS=1000
OUT_DIR=result/volume_frame/${NAME}-ws$WS-wa$WA

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

$EXE mesh=$MESH out_dir=$OUT_DIR weight.smooth=$WS weight.align=$WA lbfgs.epsf=$EPS lbfgs.maxits=$MAXITS | tee $OUT_DIR/log.txt
