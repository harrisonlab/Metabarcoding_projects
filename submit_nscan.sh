#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

DATA=$(sed -n -e "$SGE_TASK_ID p" split_files.txt)
OUT=$DATA.$1
EVAL=$2
HMM=$3

nhmmscan --noali --cpu 8 --incT 6 --tblout $OUT -E $EVAL $HMM $DATA >/dev/null
