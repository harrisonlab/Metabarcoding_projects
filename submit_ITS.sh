#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

RFILE=$1
DIR=$2
REG1=$3
REG2=$4
FASTA=$5
ID=$6

cd $DIR

Rscript $RFILE $DIR $REG1 $REG2 $FASTA $ID