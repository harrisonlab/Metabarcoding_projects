#!/bin/bash

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
PRIMERS=$5
MINL=$6
MAXR2=$7
QUAL=$8


SCRIPT_DIR=$(readlink -f ${0%/*})

qsub $SCRIPT_DIR/submit_ITSpre.sh $F $R $OUTFILE $OUTDIR $PRIMERS $MINL $MAXR2 $QUAL $SCRIPT_DIR