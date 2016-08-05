#!/bin/bash

EXE=../build/bin/test_polycube
MESH=../dat/tets/cube.c3k.vtk
MAXITS=300
OUTDIR=../result/polycube/cube/

if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

$EXE mesh=$MESH outdir=$OUTDIR weight.distortion=10 weight.onenorm=0.1 epsilon=1e-5 maxits=$MAXITS
