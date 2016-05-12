#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

FASTQ=$1
DB=$2
OUTFILE=$3
OUTDIR=$4

qsub $SCRIPT_DIR/submit_oligo.sh $FASTQ $DB $OUTFILE $OUTDIR 
