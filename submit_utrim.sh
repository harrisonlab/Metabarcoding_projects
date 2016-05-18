#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FASTQ=$1
OUTFILE=$2
OUTDIR=$3
QUAL=$4
MINL=$5
LABEL=$6


mkdir -p $OUTDIR
cd $OUTDIR

usearch8 -fastq_filter $FASTQ -fastq_minlen $MINL -fastq_maxee_rate $QUAL -relabel $LABEL -fastaout $OUTFILE
