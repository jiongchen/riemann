#!/bin/bash

EXE=build/examples/test_arap_param
MODEL=dat/beetle.obj

${EXE} ${MODEL}
mv arap_param result/
