#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

FASTA=$1
DB=$2
OUTFILE=$3
OUTDIR=$4

qsub $SCRIPT_DIR/submit_chimeras.sh $FASTA $DB $OUTFILE $OUTDIR
