#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

DATA=$1
OUT=$2
#DATA=$(sed -n -e "$SGE_TASK_ID p" split_files.txt)
#OUT=$DATA.$1

usearch9 -unoise $DATA -tabbedout ${OUT}.txt -fastaout ${OUT}.fa


