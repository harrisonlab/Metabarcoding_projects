#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

FASTQ=$1
OUTFILE=$2
OUTDIR=$3
QUAL=$4
MINL=$5
LABEL=$6

qsub $SCRIPT_DIR/submit_utrim.sh $FASTQ $OUTFILE $OUTDIR $QUAL $MINL $LABEL
