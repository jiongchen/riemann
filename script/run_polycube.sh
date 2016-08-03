#!/bin/bash

EXE=../build/bin/test_polycube
MESH=../dat/tets/sphere.c1k.vtk
W1=0.1
MAXITS=10
OUTDIR=../result/polycube/sphere/

if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

$EXE mesh=$MESH outdir=$OUTDIR weight.onenorm=$W1 weight.distortion=1 epsilon=1 maxits=$MAXITS
