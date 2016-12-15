#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=12G

RFILE=$1
HF=$2
HR=$3
OTUOUT=$4

Rscript $RFILE $HF $HR $OTUOUT