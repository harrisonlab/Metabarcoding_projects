#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
QUAL=$5
MINL=$6


qsub $SCRIPT_DIR/submit_ujoin.sh $F $R $OUTFILE $OUTDIR $QUAL $MINL
