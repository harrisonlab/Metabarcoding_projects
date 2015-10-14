#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

DATA=$1
OUT=$2
EVAL=$3
HMM=$4

nhmmscan --noali --cpu 8 --incT 6 --tblout $OUT -E $EVAL $HMM $DATA >/dev/null
