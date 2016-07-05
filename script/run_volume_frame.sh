#!/bin/bash

EXE=../build/bin/test_vol_frame
MESH=../dat/cube.tet
OUT_DIR=../result/volume_frame/cube

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

$EXE -i $MESH -o $OUT_DIR
