#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

WORK_DIR=$TMPDIR
FORWARD=$1
REVERSE=$2
OUTDIR=$3
TRIMLOC=$4


cd $WORK_DIR

java -jar $TRIMLOC/trimmomatic-0.33.jar PE -phred33 $1 $2 $1.trimmed.fq /dev/null $2.trimmed.fq /dev/null ILLUMINACLIP:$TRIMLOC/illumina_full_adapters.fa:2:30:10 SLIDINGWINDOW:8:15 MINLEN:200
