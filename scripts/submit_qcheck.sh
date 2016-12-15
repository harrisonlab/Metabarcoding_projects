#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FILE=$1
OUTDIR=$2

for SCRIPT_DIR; do true; done

fastqc $FILE -o $OUTDIR

