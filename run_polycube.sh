#!/bin/bash

echo -e "# ================================== Run L1 Polycube ======================================== #"

EXE=build/examples/test_polycube
MESH=dat/tets/eight.vtk
NAME=${MESH##*/}
EPS=1
W1=1
Wd=0.1
MAXITS=500

OUT_DIR=result/polycube/${NAME}-w1-${W1}-wd-${Wd}-eps-${EPS}

mkdir -p $OUT_DIR
$EXE mesh=$MESH outdir=$OUT_DIR epsilon=${EPS} weight.onenorm=${W1} weight.distortion=${Wd} maxits=${MAXITS} | tee $OUT_DIR/log.txt

