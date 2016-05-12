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

usearch8 -uchime_ref $FASTQ -db $DB -nonchimerasq $OUTFILE -strand plus
