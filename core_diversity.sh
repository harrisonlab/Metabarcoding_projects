#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

OTU=$1
OUT=$2
MAP=$3
TREE=$4
SEQ=$5


qsub $SCRIPT_DIR/submit_core_diversity.sh $OTU $OUT $MAP $TREE $SEQ
