#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

DATA=$(sed -n -e "$SGE_TASK_ID p" split_files.txt)
HMM1=$1
HMM2=$2
EVAL=$5

nhmmscan --noali --cpu 8 --incT 6 --tblout $DATA.$3 -E $EVAL $HMM1 $DATA >/dev/null
nhmmscan --noali --cpu 8 --incT 6 --tblout $DATA.$4 -E $EVAL $HMM2 $DATA >/dev/null