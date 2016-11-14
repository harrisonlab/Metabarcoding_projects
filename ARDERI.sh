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
	LOC=$1
	shift
 	for f in loc
	do
		R1=$f     
		R2=$(echo $R1|sed 's/_R1_/_R2_/')    
		qsub $SCRIPT_DIR/submit_demulti.sh $R1 $R2 $@ $SCRIPT_DIR
	done
	exit 1
	;;
16Spre)
	LOC=$1
	shift
	for f in $LOC
	do
		R1=$f
   		R2=$(echo $R1|sed 's/_R1_/_R2_/')
		S=$(echo $f|awk -F"_" -v D=$RUN '{print $2"D"D}')	
		qsub $SCRIPT_DIR/submit_16Spre.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
ITSpre)
	LOC=$1
	for f in $LOC
	do     
		R1=$f;     
		R2=$(echo $R1|sed 's/_R1_/_R2_/');     
		S=$(echo $f|awk -F"_" -v D=$RUN '{print $2"D"D}');
		qsub $SCRIPT_DIR/submit_ITSpre.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
procends)
	LOC=$1
	shift	

	cd $LOC
	for f in *.fa
	do
		d=$(echo $f|awk -F"." '{print $1}')
		mkdir $d
		split -l 2000 $f -a 3 -d ${d}/$f.
		cd $d
		find $PWD -name '*.fa.*' >split_files.txt
		TASKS=$(wc -l split_files.txt|awk -F" " '{print $1}')
        		qsub -t 1-$TASKS:1 $SCRIPT_DIR/submit_nscan.sh $@
		cd ..
	done
	exit 1
	;;
ITS)
	LOC=$1
	shift
	for d in $LOC
	do
		S=$(echo $d|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}')
		qsub $SCRIPT_DIR/submit_ITS.sh $d $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
merge_hits)
	qsub $SCRIPT_DIR/submit_merge_hits.sh $@ $SCRIPT_DIR
	exit 1
	;;
UPARSE|uparse)
	qsub $SCRIPT_DIR/submit_uparse_v2.sh $@ $SCRIPT_DIR
	exit 1
	;;
OTU|otu)
	qsub $SCRIPT_DIR/submit_otu.sh $@ $SCRIPT_DIR
	exit 1
	;;
tax_assign)
	qsub $SCRIPT_DIR/submit_taxonomy.sh $@ $SCRIPT_DIR
	exit 1
	;;
denoise)
	qsub $SCRIPT_DIR/submit_denoise.sh $@ $SCRIPT_DIR
	exit 1
	;;
TEST)
	X=$1
	shift
	echo $program $@ $X $SCRIPT_DIR
	exit 1
	;;
*)
	echo "Invalid program: $program" >&2
	exit 1
esac
	