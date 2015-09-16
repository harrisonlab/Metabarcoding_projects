#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

WORK_DIR=$HOME/tmp

DATA=$1
OUT=$2
PARAMS=$3
REF=$4
WITHTABLE=$5

cd $WORK_DIR

if [ -n $WITHTABLE ] ;
then
  echo "running with table"
  pick_open_reference_otus.py -a -i $DATA -o $OUT -p $PARAMS -r $REF 
else
  echo "running without table"
  pick_open_reference_otus.py -a -i $DATA -o $OUT -p $PARAMS -r $REF --suppress_align_and_tree
fi
