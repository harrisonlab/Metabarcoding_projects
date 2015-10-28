#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
QUAL=$5
MINL=$6


mkdir -p $OUTDIR
cd $OUTDIR

usearch8 -fastq_mergepairs $F -reverse $R -fastqout $OUTFILE -fastq_merge_maxee $QUAL -fastq_minlen $MINL

