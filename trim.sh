#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

#R1=$1
#R2=$2
OUT=$1
TRIMLOC=$2
counter=0

for f in $OUT/*
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
		qsub $SCRIPT_DIR/submit_trim.sh $R1 $R2 $OUT $TRIMLOC 
	fi
	R1=$f
done
