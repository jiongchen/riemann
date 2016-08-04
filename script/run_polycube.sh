#!/bin/bash

EXE=../build/bin/test_polycube
MESH=../dat/tets/sphere.c1k.vtk
MAXITS=100
OUTDIR=../result/polycube/sphere/

if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

$EXE mesh=$MESH outdir=$OUTDIR weight.distortion=1 weight.onenorm=0.1 epsilon=1 maxits=$MAXITS
