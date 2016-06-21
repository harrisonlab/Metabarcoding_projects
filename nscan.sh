#!/bin/bash

DATA=$1
OUT=$2
EVAL=${3:-0.0015}
HMM=${4:-58endjames23.HMM}


SCRIPT_DIR=$(readlink -f ${0%/*})

#echo "Data:" $DATA
#echo "Outdir:" $OUT
#echo "Eval:" $EVAL
#echo "Hmm:" $HMM

qsub $SCRIPT_DIR/submit_nscan.sh $DATA $OUT $EVAL $HMM
