#!/bin/bash

OPTIND=1 

SCRIPT_DIR=$(readlink -f ${0%/*})
S=$SCRIPT_DIR


while getopts ":hc:" options; do
	case "$options" in
	c)  
 	    program=$OPTARG
	    ;;
	h)  
	    exit 1
 	    ;;
	\?) 
	    echo "Invalid option: -$OPTARG" >&2
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1
      	    ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift


case $program in
demulti|demultiplex) 
	echo qsub $SCRIPT_DIR/submit_demulti.sh $@ $SCRIPT_DIR
	exit 1
	;;
16Spre)
	qsub $SCRIPT_DIR/submit_16Spre.sh $@ $SCRIPT_DIR
	exit 1
	;;
ITSpre)
	qsub $SCRIPT_DIR/submit_ITSpre.sh $@ $SCRIPT_DIR
	exit 1
	;;
ITS)
	qsub $SCRIPT_DIR/submit_ITS.sh $@ $SCRIPT_DIR
	exit 1
	;;
merge_hits)
	qsub $SCRIPT_DIR/submit_merge_hits.sh $@ $SCRIPT_DIR
	exit 1
	;;
UPARSE|uparse)
	qsub $SCRIPT_DIR/submit_uparse.sh $@ $SCRIPT_DIR
	exit 1
	;;
TEST)
	echo $program $@ $SCRIPT_DIR
	exit 1
	;;
*)
	echo "Invalid program: $program" >&2
	exit 1
esac
	