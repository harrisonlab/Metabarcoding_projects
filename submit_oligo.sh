#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FASTQ=$1
DB=$2
OUTFILE=$3
OUTDIR=$4

mkdir -p $OUTDIR
cd $OUTDIR

usearch8.1 -search_oligodb $FASTQ -db $DB -strand both -userout $OUTFILE -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
