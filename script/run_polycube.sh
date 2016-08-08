#!/bin/bash

EXE=../build/bin/test_polycube
MESH=../dat/tets/torus.c2k.vtk
MAXITS=200
OUTDIR=../result/polycube/torus/

if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

$EXE mesh=$MESH outdir=$OUTDIR weight.distortion=1 weight.onenorm=0.1 epsilon=1 maxits=$MAXITS
