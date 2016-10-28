#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})
S=$SCRIPT_DIR

qsub $SCRIPT_DIR/submit_demulti.sh $@ $SCRIPT_DIR
