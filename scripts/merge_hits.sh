#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

RFILE=$1
HF=$2
HR=$3
OTUOUT=$4

qsub $SCRIPT_DIR/submit_merge_hits.sh $RFILE $HF $HR $OTUOUT $SCRIPT_DIR