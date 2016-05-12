#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4


mkdir -p $OUTDIR
cd $OUTDIR

usearch8 -fastq_mergepairs $F -reverse $R -fastqout $OUTFILE

