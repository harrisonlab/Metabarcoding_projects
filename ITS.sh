#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

RFILE=$1
DIR=$2
REG1=$3
REG2=$4
FASTA=$5
ID=$6
LOWQUAL=$7

qsub $SCRIPT_DIR/submit_ITS.sh $RFILE $DIR $REG1 $REG2 $FASTA $ID $LOWQUAL