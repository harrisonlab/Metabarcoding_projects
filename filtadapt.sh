#!/bin/bash

READ1=$1
ADAPTERS=$2

SCRIPT_DIR=$(readlink -f ${0%/*})

qsub $SCRIPT_DIR/submit_filtadapt.sh $READ1 $ADAPTERS