# apple_replant
Metagenomic study of apple replant disease  
 1. Installing Qiime to a local directory  
  2. Parellel Qiime  
 3. Common workflow
  4. QC  
  5. Trimming
 6. 16S workflow
  7. Join PE reads 
  8. Rename files  
  9. Convert joined fastq to fasta  
  10. Remove chimeras  
  11. Concatenate files
  12. OTU Picking and descriptive statistics
  13. Statistical analysis  
 14. ITS workflow  
  15. Convert to unpaired fasta files  
  16. Rename files  
  17. SSU/58S/LSU removal  
    18. Split fasta into chunks for SSU/58S/LSu removal
    19. Remove SSU/LSU
    20. Merge output
  21. Remove chimeras
  22. Return merged common ITS1 and ITS2, unique ITS1 and unique ITS2
  23. OTU Picking and descriptive statistics
    24. Common and unique (ITS1 and ITS2)  
    25. Common ITS
    26. Unique ITS1 only
    27. Unique ITS2 only
  28. Statistical analysis  

## Installing Qiime to a local directory
Downloaded Python 2.7.9 tar ball and unzipped.  
From root of Python 2.7.9 directory ran :

	./configure --prefix=/home/deakig/usr/local --exec-prefix=/home/deakig/usr/local --enable-unicode=ucs4
	make
	make install

Downloaded pip tarball amd unzipped to pip directory then ran:

	~/usr/local/bin/python ~/pip/getpip.py


Set Qiime path with below (not permanent)

	export PYTHONUSERBASE=/home/deakig/usr/local/
	
	
	
To install Qiime and dependencies

	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall numpy
	~/usr/local/bin/python -m pip install --user --upgrade --force-reinstall qiime
	
(the upgrade and force-reinstall flags may not be necessary)

To test qiime, ensure ~/usr/local/bin (the qiime script directory) is in path

	export PATH=$PATH:/home/deakig/usr/local/bin

then

	 ~/usr/local/bin/python ~/usr/local/bin/print_qiime_config.py -t

should retun something like

	$> Ran 9 test in 0.05s
	$> OK

### Parallel qiime
for single machine throw in -a -O (no. processes) to the workflow script

using HPC... 
create qimme_config in home root

	cd ~
	touch .qiime_config

added to qimme_config:  
jobs_to_start 8  
temp_dir $HOME/tmp  
cluster_jobs_fp start_parallel_jobs_sc.py	

Hacked start_parallel_jobs_sc.py for use in our environment. Changed the qsub template settings as bellow:    
\# qsub template  
QSUB_TEXT = """#!/bin/bash  
\#$ -S %s  
\#$ -l %s  
\#$ -cwd  

___
## Common workflow
The following directories should be created prior to starting the workflow:

	$METAGENOMICS
	$METAGENOMICS/analysis
	$METAGENOMICS/analysis/16S
	$METAGENOMICS/analysis/ITS	
	$METAGENOMICS/data
	$METAGENOMICS/data/fastq
	$METAGENOMICS/data/trimmed
	$METAGENOMICS/data/joined
	$METAGENOMICS/data/fasta
	$METAGENOMICS/data/16S
	$METAGENOMICS/data/16S/de_chimeraed
	$METAGENOMICS/data/ITS
	$METAGENOMICS/data/ITS/de_chimeraed
	$METAGENOMICS/data/ITS/final
	$METAGENOMICS/hmm
	$METAGENOMICS/scripts
	$METAGENOMICS/taxonomies
	
The $METAGENOMICS directory should be set to something appropriate (e.g. /home/bob/metagenomics). The realtive path is used in the scripts below - depending on your config you may have to specify full paths.	

### QC
Qualtiy checking was performed with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

From same folder containing fastq files ran:

	fastqc *

#### Demulitplexing
We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Takes as input paired data and will output two files for each. Sequence which doesn't match either index is written to both fungal and bacterial fastq files.

```shell
counter=0
for f in $METAGENOMICS/data/fastq/*
do counter=$((counter+1))
    if (( $counter % 2 == 0 ))
    then
  	R2=$f
	echo $f
	# replace index_1 and 2 with a regular expression for each index
	$METAGENOMICS/scripts/demulti.pl $R1 $R2 "^index_1" "^index_2"	
    fi
    R1=$f
done
```


### Trimming
####THIS NEEDS TO BE CHANGED - I've scrapped Trimmomatic
Paired end trimming was preformed with Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic).

The following settings were used:
- minimum length 200
- sliding window
 + 8 bases
 + quality 15
- illumina adapter clipping
 + 2 mismatches 
 + palindrome quality 30
 + clip threshold quality 10

Shell script trim.sh used to submit trimming jobs to cluster. 

The first argument specifies the input/output directory of paired end fastq files (all file in folder will be processed - paired files must sort adjacently)

Second argument specifies location of illumina adapter file

UPDATE - trim now contains parameters to pass trimming phred quality and minumum length

	
	$METAGENOMICS/scripts/trim.sh %METAGENOMICS/data/fastq $METAGENOMICS/scripts quality minlength


Then

	mv $METAGENOMICS/data/fastq/*trimmed* $METAGENOMICS/data/trimmed/.


The following 16S and ITS workflows are dependent on specific naming conventions for the samples - most of the bash scripts have been written to work with sequential sample naming  (S86 - S96 for the data presented below). One of the R scripts also uses the same sample name format for fast sorting via a regex. 

## 16s workflow

### Join PE reads
####USE ALTERNATE METHOD ONLY
The 'join PE script' was run form the trimmed directrory  
The counter used in the next couple of commands was set to match the names of the samples, i.e. S85, S86 and etc.
```shell

cd $METAGENOMICS/data/trimmed

###join PE script
counter=0;
X=85;
for f in *trimmed*; 
do counter=$((counter+1)); 
	if (( $counter % 2 == 0 )); 
		then R2=$f;
		echo join_paired_ends.py -f $R1 -r $R2 -o S$X;
		join_paired_ends.py -f $R1 -r $R2 -o S$X; 
		X=$((X+1));
	fi; 
	R1=$f; 
done
```
##### Alternative method using usearch (for untrimmed data)
usearch trims based on the expected error for the entire joined sequence.
Expected error set to 1 in below and min length set to 200
```shell
counter=0
for f in $METAGENOMICS/data/1910/fastq/16S/*
do counter=$((counter+1))
	if (( $counter % 2 == 0 ))
	then
		R2=$f
		S=$(echo $f|awk -F"_" '{print $2}')
		$METAGENOMICS/scripts/ujoin.sh $R1 $R2 ${S}.joined.fq $METAGENOMICS/data/1910/joined 1 200
	fi
	R1=$f
done
```

### Rename files - THIS IS NO LONGER NECESSARY 
Moved joined directories/files to the $METAGENOMICS/data/joined directory (it is important to ensure there are no files in the root of the joined directory or you risk renaming all files in lower level directories)
	
	mv $METAGENOMICS/data/trimmed/S* $METAGENOMICS/data/joined/.

Then ran the below:
```shell
cd $METAGENOMICS/data/joined
counter=85
for d in * 
do 
	cd S$counter
	for f in *
	do
		echo mv -i "${f}" "S${f/fastqjoin/$counter}"
	done
	cd ..
	counter=$((counter+1));
done
```
### Convert joined fastq to fasta
must be run from root of joined directory 

cd  $METAGENOMICS/data/$RUN/16S/joined/	
	
```shell
for f in  *join*
do
	S=$(echo $f|awk -F"." '{print $1}')
	$METAGENOMICS/scripts/fq2fa.pl $f $f.fa $S
	mv $f.fa $METAGENOMICS/data/$RUN/16S/fasta/.
done

```
### Remove chimeras
Downloaded usearch 8.0 and RDP gold reference database from http://drive5.com/usearch/manual/cmd_uchime_ref.html

Ran the 'remove chimeras script'

```shell
#remove chimeras script 	
counter=85
for f in $METAGENOMICS/data/fasta/16S/*
do 
	$METAGENOMICS/scripts/chimeras.sh $f $METAGENOMICS/taxonomies/RDP_gold.fasta S${counter}.cfree.fa
	$METAGENOMICS/data/fasta/16S/de_chimeraed/
	counter=$((counter+1));
done
```
#### Concatenate files
Concatenated all the de-chimeraed files and copied the output to the $METAGENOMICS/data/fasta/16S directory

	cat $METAGENOMICS/data/fasta/16S/de_chimeraed/*cfree* > $METAGENOMICS/data/fasta/16S/16S.joined.fa

### OTU Picking and descriptive statistics
Run the 2nd and 3rd commands below only after the cluster jobs created by the 1st command have finished
```shell
$METAGENOMICS/scripts/pick_OTU.sh   $METAGENOMICS/data/fasta/16S/16S.joined.fa  $METAGENOMICS/analysis/16S/16S_otus $METAGENOMICS/scripts/parameters.txt $PYTHONUSERBASE/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta TRUE
 X=`biom summarize-table -i METAGENOMICS/analysis/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"`
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom $METAGENOMICS/analysis/16S/16s_cdout/ $METAGENOMICS/data/map.tsv $METAGENOMICS/analysis/16S/16S_otus/rep_set.tre $X
```
### Statistical analysis
analysis.R biom_table colData median/geomean outfile  

Requires a file (colData) which describes condition (e.g. infected or uninfected) for each sample 
```shell
cd $METAGENOMICS/analysis/16S
Rscript $METAGENOMICS/scripts/analysis.R "analysis/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom" colData median res.sig.csv
```	
## ITS workflow

### Alternative trimming method with usearch
utrim is using the expected error per base. The settings below (which also set minimum length to 200) will discard sequences of 200 bases if expected error is > 1 - this is for the forward read only, the reverse read is not as stringent due to fairly poor quality of data in this example.
```shell
counter=0;
for f in $METAGENOMICS/data/1910/fastq/ITS/*
do counter=$((counter+1))
	S=$(echo $f|awk -F"_" '{print $2}')
	if (( $counter % 2 == 0 ))
	then
		$METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.2.fq $METAGENOMICS/data/1910/trimmed 0.02 200
	else
		echo $METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.1.fq $METAGENOMICS/data/1910/trimmed 0.005 200
	fi
done
```

### Convert to unpaired fasta files

```shell
X=91
counter=0
for f in $METAGENOMICS/data/trimmed/*trimmed*;
do counter=$((counter+1));
  if [ "$counter" -gt 12 ]
  then
    if (( $counter % 2 == 0 ))
    then
      $METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/fasta/ITS/${f}.fa S$X ;
      X=$((X+1))
    else
      $METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/fasta/ITS/${f}.fa S$X ;
    fi
  fi
done
```
### Rename files 
Moved to fasta/ITS directory then ran: 
```shell
cd $METAGENOMICS/fasta/ITS

counter=0
for f in *trimmed*;
do counter=$((counter+1));
	S=$(echo $f|awk -F"." '{print $1}')
	if (( $counter % 2 == 0 ))
	then
		mv -i "${f}" "${S}_R2.fa"
	else
		mv -i "${f}" "${S}_R1.fa"
	fi
done

```
### SSU/58S/LSU removal
Using HHMMER v 3.1b2 (http://hmmer.janelia.org/)

Used HMM files from ITSx (http://microbiology.se/software/itsx/)

```shell
perl $METAGENOMICS/scripts/cut_hmm v.3.1 $METAGENOMICS/hmm/chopped_hmm fungi
cd $METAGENOMICS/hmm/chopped_hmm
cat *SSU*> ssu_end.hmm
cat *58S_start* > 58s_start.hmm
cat *58S_end* > 58s_end.hmm
cat *LSU* > lsu_start.hmm
hmmconvert ssu_end.hmm > ssu_end.hmm
hmmconvert 58s_end.hmm > 58s_end.hmm
hmmconvert 58s_start.hmm > 58s_start.hmm
hmmconvert lsu_start.hmm > lsu_start.hmm

for f in *.hmm
do
	sed -e'/^LENG/a MAXL  90' $f >> $f.temp
done

hmmpress ssu_end.hmm
hmmpress 58s_end.hmm
hmmpress 58s_start.hmm
hmmpress lsu_start.hmm
```
Ouptut files were copied to $METAGENOMICS/hmm. Hacked the HMM files to include a MAXL satement (required) and manually split out SSU,58S and LSU into seperate files (using fungal only)

##### Split fasta into chunks for SSU/58S/LSU removal
```shell
cd $METAGENOMICS/fasta/ITS
X=91 
counter=0
for f in *.fa;
do counter=$((counter+1));
  if (( $counter % 2 == 0 ))
  then
    mkdir S${X}_R2
    split -l 2000 $f -a 3 -d S${X}_R2/$f.
    X=$((X+1))
  else
    mkdir S${X}_R1
    split -l 2000 $f -a 3 -d S${X}_R1/$f.
  fi
done
```
##### Remove SSU/LSU
Note - creates a file with the paths to all of the split files in each sample directory then submits cluster array job
(nscan.sh is no longer used)
```shell
cd $METAGENOMICS/fasta/ITS

for d in */
do
cd $d
find $PWD -name '*.fa.*' >split_files.txt
cd ..
done

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
##### Merge output
(bug fixed)
```shell
for d in $METAGENOMICS/data/fasta/ITS/*R1
do
	 $METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_SSU_58Ss.R $d "*.\\.ssu" "*.\\.58" $d.fa
done

for d in $METAGENOMICS/data/trimmed_q20/fasta/ITS/*R2
do
	 $METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_58Se_LSU.R $d "*.\\.58" "*.\\.lsu" $d.fa
done
```
### Remove empty fastas
```shell
cd $METAGENOMICS/data/fasta/ITS
counter=0
for d in */;
do counter=$((counter+1));
	cd $d
	if (( $counter % 2 == 0 ))
	then
		awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ITS2.fa > ITS2.t.fa
	else
		awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ITS1.fa > ITS1.t.fa
	fi
	cd ..
done
```
### Remove chimeras
Using UNITE v 7.0 ITS database for chimeras (UCHIME reference dataset) https://unite.ut.ee/repository.php#uchime

```shell
counter=91
counter2=1
for d in $METAGENOMICS/data/fasta/ITS/*R[0-9]
do 
  if (( $counter2==1 ))
  then
    $METAGENOMICS/scripts/chimeras.sh $d/ITS1.t.fa $METAGENOMICS/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta S${counter}.${counter2}.cfree.fa $METAGENOMICS/data/fasta/ITS/de_chimerad/
    counter2=2
  else
    $METAGENOMICS/scripts/chimeras.sh $d/ITS2.t.fa $METAGENOMICS/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta S${counter}.${counter2}.cfree.fa $METAGENOMICS/data/ITS/fasta/de_chimerad/
    counter2=1
    counter=$((counter+1));	
  fi
done
```
### Return merged common ITS1 and ITS2, unique ITS1 and unique ITS2
```shell	
cd $METAGENOMICS/data/fasta/ITS/final

counter=0;
X=91
for f in $METAGENOMICS/data/fasta/ITS/de_chimeraed/*
do counter=$((counter+1)); 
	if (( $counter % 2 == 0 ))
	then
		R2=$f;
		$METAGENOMICS/scripts/catfiles.pl $R1 $R2 S$X;
		X=$((X+1));
	fi
	R1=$f
done
```
### OTU Picking and descriptive statistics
Multiple analyses were perfomed on:

1. Common and unique
2. Common ITS1 and ITS2
3. Unique ITS1
4. Unique ITS2

##### Common and unique (ITS1 and ITS2)
```shell
cd $METAGENOMICS/data/fasta/ITS/final
cat *.fa > ITS.all.fa
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/fasta/ITS/final/ITS.all.fa $METAGENOMICS/analysis/ITS/ITS_all_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom
X=`biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom $METAGENOMICS/analysis/ITS/ITS_all_cdout/ $METAGENOMICS/data/map.tsv . $X
```
##### Common ITS
```shell
cat S91.fa S92.fa S93.fa S94.fa S95.fa S96.fa > ITS.common.fa
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/fasta/ITS/final/ITS.common.fa $METAGENOMICSs/analysis/ITS/ITS_common_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS_common_otus/otu_table_mc2_w_tax.biom
X=`biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS_common_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/ITS/ITS_common_otus/otu_table_mc2_w_tax.biom $METAGENOMICS/analysis/ITS/ITS_common_cdout/ $METAGENOMICS/data/map.tsv . $X
```
##### Unique ITS1 only
```shell
cat *r1* >ITS1.only.fa
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/fasta/ITS/final/ITS1.only.fa $METAGENOMICS/analysis/ITS/ITS1_only_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS1_only_otus/otu_table_mc2_w_tax.biom
X=`biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS1_only_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/ITS/ITS1_only_otus/otu_table_mc2_w_tax.biom $METAGENOMICS/analysis/ITS/ITS1_only_cdout/ $METAGENOMICS/data/map.tsv . $X
```
##### Unique ITS2 only
```shell
cat *r2* >ITS2.only.fa
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/fasta/ITS/final/ITS2.only.fa /$METAGENOMICS/analysis/ITS/ITS2_only_otus $METAGENOMICSs/scripts/params.txt $METAGENOMICS/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS2_only_otus/otu_table_mc2_w_tax.biom
X=`biom summarize-table -i $METAGENOMICS/analysis/ITS/ITS2_only_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/ITS/ITS2_only_otus/otu_table_mc2_w_tax.biom $METAGENOMICS/analysis/ITS/ITS2_only_cdout/ $METAGENOMICS/data/map.tsv . $X
```
### Statistical analysis
analysis.R biom_table colData median/geomean outfile  

Requires a file (colData) which describes condition (e.g. infected or uninfected) for each sample

As there were multiple OTU picking steps, multiple statistical analyses are necessary.
```shell
Rscript $METAGENOMICS/scripts/analysis.R $METAGENOMICS/analysis/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom colData median ITS.all.median.csv
Rscript $METAGENOMICS/scripts/analysis.R $METAGENOMICS/analysis/ITS/ITS_common_otus/otu_table_mc2_w_tax.biom colData median ITS.median.csv
Rscript $METAGENOMICS/scripts/analysis.R $METAGENOMICS/analysis/ITS/ITS1_only_otus/otu_table_mc2_w_tax.biom colData median ITS1.median.csv
Rscript $METAGENOMICS/scripts/analysis.R $METAGENOMICS/analysis/ITS/ITS2_only_otus/otu_table_mc2_w_tax.biom colData median ITS2.median.csv
```
