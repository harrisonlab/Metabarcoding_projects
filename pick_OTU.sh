#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

DATA=$1
OUT=$2
PARAMS=$3
REF=$4
WITHTABLE=$5

qsub $SCRIPT_DIR/submit_pick_OTU.sh $DATA $OUT $PARAMS $REF $WITHTABLE
