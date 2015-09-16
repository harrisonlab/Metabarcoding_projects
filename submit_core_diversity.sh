#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

WORK_DIR=$HOME/tmp

OTU=$1
OUT=$2
MAP=$3
TREE=$4
SEQ=$5

cd $WORK_DIR

#core_diversity_analyses.py -a -o $OUT -i $OTU -m $MAP -t $TREE -e $SEQ 
core_diversity_analyses.py -a -o $OUT -i $OTU -m $MAP -e $SEQ --nonphylogenetic_diversity
