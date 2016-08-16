# apple_replant
Metagenomic study of apple replant disease  
 1. HMM Preperation for ITS analysis  
 2. Common workflow 
 6. 16S workflow
 18. ITS workflow
 35. Oomycetes workflow
 36. Combine samples
 37. Old (Qiime method)


## HMM Preperation for ITS analysis
Using HHMMER v 3.1b2 (http://hmmer.janelia.org/)

Used HMM files from ITSx (http://microbiology.se/software/itsx/)

```shell
perl $METAGENOMICS/scripts/cut_hmm v.3.1 $METAGENOMICS/hmm/chopped_hmm fungi
cd $METAGENOMICS/hmm/chopped_hmm
cat *SSU*> t1
cat *58S_start* > t2
cat *58S_end* > t3
cat *LSU* > t4
hmmconvert t1 > ssu_end.hmm
hmmconvert t2 > 58s_end.hmm
hmmconvert t3 > 58s_start.hmm
hmmconvert t4 > lsu_start.hmm

rm t1 t2 t3 t4

for f in *.hmm
do
	sed -i -e'/^LENG/a MAXL  90' $f
done

hmmpress ssu_end.hmm
hmmpress 58s_end.hmm
hmmpress 58s_start.hmm
hmmpress lsu_start.hmm
```
Ouptut files were copied to $METAGENOMICS/hmm. Hacked the HMM files to include a MAXL satement (required) and manually split out SSU,58S and LSU into seperate files (only fungal hmms are implemented in this pipeline)

## Utax reference databases
Reference databases were downloaded from:
http://drive5.com/usearch/manual/utax_downloads.html
(Unite V7 and RDP trainset 15)
```shell
usearch8.1 -makeudb_utax refdb.fa -output 16s_ref.udb -report 16s_report.txt
usearch8.1 -makeudb_utax refdb.fa -utax_trainlevels kpcofgs ‑utax_splitlevels NVpcofgs -output ITS_ref.udb -report ITS_report.txt
```

## Set directories
The following directories should be created prior to starting the workflow:
```shell
#run once
mkdir $METAGENOMICS
mkdir $METAGENOMICS/analysis
mkdir $METAGENOMICS/data
mkdir $METAGENOMICS/hmm
mkdir $METAGENOMICS/scripts
mkdir $METAGENOMICS/taxonomies
```

___
https://github.com/harrisonlab/apple_replant/blob/master/Common%20workflow.md
##Common workflow
[Common workflow](../Common%20workflow.md)

##[16S workflow](https://github.com/harrisonlab/apple_replant/blob/master/16S%20%20workflow.md)
##ITS workflow
https://github.com/harrisonlab/apple_replant/blob/master/ITS%20%20workflow.md


##oomycetes
```shell
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/$RUN/ITS/final/ITS.all.fa $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/Silva119/97_18S_only/Silva_119_rep_set97_aligned_18S_only.fna FALSE
```

##Combine samples
Biom table for samples from multiple NGS runs are required.

This will mean the names of each fasta will need to be made unique and the sequence lengths will need to be set to the same.

Something like the below will copy samples with the wholename string to a new location. Uses original fastq file name and reconsructing the sample ID for each sample used in the workflow.
```shell
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/16S/filtered/"$3"D"$1".filtered.fa 16S/filtered/."}' > runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/16S/unfiltered/"$3"D"$1".unfiltered.fastq 16S/unfiltered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/filtered/"$3"D"$1".fa ITS/filtered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/unfiltered/"$3"D"$1".r1.unfiltered.fastq ITS/unfiltered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R2*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/unfiltered/"$3"D"$1".r2.unfiltered.fastq ITS/unfiltered/."}' >> runme.sh
./runme.sh
```

## old stuff

##### Split fasta into chunks 
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
for f in *.fa;
do
  S=$(echo $f|awk -F"." '{print $1}')
    mkdir $S
    split -l 2000 $f -a 3 -d ${S}/$f.
done
```
##### Create lists of file paths to chunks

```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
for d in */
do
	cd $d
	find $PWD -name '*.fa.*' >split_files.txt
	cd ..
done
```
##### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
counter=0
for d in */
do counter=$((counter+1))
	cd $d
	TASKS=$(wc -l split_files.txt|awk -F" " '{print $1}')
	if (( $counter % 2 == 0 ))
	then
		qsub -t 1-$TASKS:1 $METAGENOMICS/scripts/submit_nscan.sh lsu 20 $METAGENOMICS/hmm/lsu_start.hmm
		qsub -t 1-$TASKS:1 $METAGENOMICS/scripts/submit_nscan.sh 58se 20 $METAGENOMICS/hmm/58s_end.hmm
	else
		qsub -t 1-$TASKS:1 $METAGENOMICS/scripts/submit_nscan.sh ssu 20 $METAGENOMICS/hmm/ssu_end.hmm
		qsub -t 1-$TASKS:1 $METAGENOMICS/scripts/submit_nscan.sh 58ss 20 $METAGENOMICS/hmm/58s_start.hmm
	fi
	cd ..
done
```
