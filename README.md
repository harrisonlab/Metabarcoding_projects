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

	./configure --prefix=$HOME/usr/local --exec-prefix=$HOME/usr/local --enable-unicode=ucs4
	make
	make install

Downloaded pip tarball amd unzipped to pip directory then ran:

	~/usr/local/bin/python ~/pip/getpip.py


Set Qiime path with below (not permanent)

	export PYTHONUSERBASE=$HOME/usr/local/
	
	
	
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
```shell
#run once
mkdir $METAGENOMICS
mkdir $METAGENOMICS/analysis
mkdir $METAGENOMICS/data
mkdir $METAGENOMICS/hmm
mkdir $METAGENOMICS/scripts
mkdir $METAGENOMICS/taxonomies
```
```shell
#run for each analysis
mkdir $METAGENOMICS/analysis/$RUN/16S
mkdir $METAGENOMICS/analysis/$RUN/ITS	
mkdir $METAGENOMICS/data/$RUN
mkdir $METAGENOMICS/data/$RUN/fastq
mkdir $METAGENOMICS/data/$RUN/PhiX
mkdir $METAGENOMICS/data/$RUN/16S
mkdir $METAGENOMICS/data/$RUN/16S/fastq
mkdir $METAGENOMICS/data/$RUN/16S/fasta
mkdir $METAGENOMICS/data/$RUN/16S/joined
mkdir $METAGENOMICS/data/$RUN/16S/de_chimeraed
mkdir $METAGENOMICS/data/$RUN/ITS
mkdir $METAGENOMICS/data/$RUN/ITS/fastq
mkdir $METAGENOMICS/data/$RUN/ITS/fasta
mkdir $METAGENOMICS/data/$RUN/ITS/trimmed
mkdir $METAGENOMICS/data/$RUN/ITS/de_chimeraed
mkdir $METAGENOMICS/data/$RUN/ITS/final

```	
The $METAGENOMICS directory should be set to something appropriate (e.g. /home/bob/metagenomics) and $RUN to the name of the NGS run. The realtive path is used in the scripts below - depending on your config you may have to specify full paths.	

### QC
Qualtiy checking was performed with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

From same folder containing fastq files ran:

	fastqc *

### PhiX filtering
For the particular sequencing protocal we don't get much (or any) PhiX contamination. Removal of any contaminants is simple via aligning to the Illumina PhiX genome <ln>http://support.illumina.com/sequencing/sequencing_software/igenome.html </ln> Bowtie2 method implemented here
```shell
 counter=0
 for f in $METAGENOMICS/data/$RUN/fastq/*.fastq
 do counter=$((counter+1))
 if (( $counter % 2 == 0 ))
     then
         R2=$f
         S=$(echo $f|awk -F"_" '{print $2}')
         $METAGENOMICS/scripts/bowtie.sh $R1 $R2 $HOME/Data/PhiX/Illumina/RTA/Sequence/Bowtie2Index/genome $METAGENOMICS/data/$RUN/PhiX ${S}.phix.fq 250 500
     fi
     R1=$f
done
```

#### Demulitplexing
We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Takes as input paired data and will output two files for each. Sequence which doesn't match either index is written to both fungal and bacterial fastq files.

Running something like the below should give a good indication of what index_1 and index_2 should be. 
```shell
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-8|sort|uniq > expressions.txt
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-8|sort|uniq|xargs -I r grep -c ^r $(ls|head -n1) >counts.txt
```

```shell
counter=0
for f in $METAGENOMICS/data/$RUN/PhiX/*
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

mv *bacterial* ../16S/fastq/.
mv *fungal* ../ITS/fastq/.
```


## 16s workflow

### Join PE reads
usearch trims based on the expected error for the entire joined sequence.
Expected error set to 1 in below and min length set to 200
```shell
counter=0
for f in $METAGENOMICS/data/$RUN/16S/fastq/*
do counter=$((counter+1))
	if (( $counter % 2 == 0 ))
	then
		R2=$f
		S=$(echo $f|awk -F"_" '{print $2}')
		$METAGENOMICS/scripts/ujoin.sh $R1 $R2 ${S}.joined.fq $METAGENOMICS/data/$RUN/16S/joined 1 200
	fi
	R1=$f
done
```
### Convert joined fastq to fasta
must be run from root of joined directory 

```shell
cd  $METAGENOMICS/data/$RUN/16S/joined/	

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
for f in $METAGENOMICS/data/$RUN/16S/fasta/*
do
	S=$(echo $f|awk -F"." '{print $1}')
	$METAGENOMICS/scripts/chimeras.sh $f $METAGENOMICS/taxonomies/RDP_gold.fasta ${S}.cfree.fa $METAGENOMICS/data/$RUN/16S/de_chimeraed/
done
```
#### Concatenate files
Concatenated all the de-chimeraed files and copied the output to the $METAGENOMICS/data/$RUN/16S directory

	cat $METAGENOMICS/data/$RUN/16S/de_chimeraed/*cfree* > $METAGENOMICS/data/$RUN/16S/16S.fa

### OTU Picking and descriptive statistics
Run the 2nd and 3rd commands below only after the cluster jobs created by the 1st command have finished
```shell
$METAGENOMICS/scripts/pick_OTU.sh   $METAGENOMICS/data/$RUN/16S/16S.fa  $METAGENOMICS/analysis/$RUN/16S/16S_otus $METAGENOMICS/scripts/parameters.txt $PYTHONUSERBASE/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta TRUE
 X=`biom summarize-table -i METAGENOMICS/analysis/$RUN/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"`
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/$RUN/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom $METAGENOMICS/analysis/$RUN/16S/16s_cdout/ $METAGENOMICS/data/map.tsv $METAGENOMICS/analysis/$RUN/16S/16S_otus/rep_set.tre $X
```
### Statistical analysis
analysis.R biom_table colData median/geomean outfile  

Requires a file (colData) which describes condition (e.g. infected or uninfected) for each sample 
```shell
cd $METAGENOMICS/analysis/$RUN/16S/16S_otus
Rscript $METAGENOMICS/scripts/analysis.R "otu_table_mc2_w_tax_no_pynast_failures.biom" colData median res.sig.csv
```	
## ITS workflow

### Alternative trimming method with usearch
utrim is using the expected error per base. The settings below (which also set minimum length to 200) will discard sequences of 200 bases if expected error is > 1 - this is for the forward read only, the reverse read is not as stringent due to fairly poor quality of data in this example.
```shell
counter=0;
for f in $METAGENOMICS/data/$RUN/ITS/fastq/*
do counter=$((counter+1))
	S=$(echo $f|awk -F"_" '{print $2}')
	if (( $counter % 2 == 0 ))
	then
		$METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.2.fq $METAGENOMICS/data/$RUN/ITS/trimmed 0.02 200
	else
		$METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.1.fq $METAGENOMICS/data/$RUN/ITS/trimmed 0.005 200
	fi
done
```

### Convert to unpaired fasta files
This might have some use, but can't remeber what - alternative simpler method should work...
```shell
X=91
counter=0
for f in $METAGENOMICS/data/$RUN/ITS/trimmed/*trimmed*;
do counter=$((counter+1));
  if [ "$counter" -gt 12 ]
  then
    if (( $counter % 2 == 0 ))
    then
      $METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/ITS/fasta/${f}.fa S$X ;
      X=$((X+1))
    else
      $METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/ITS/fasta/${f}.fa S$X ;
    fi
  fi
done
```
Alternative method
```shell
cd $METAGENOMICS/data/$RUN/ITS/trimmed

for f in *trimmed*;
do
	S=$(echo $f|awk -F"." '{print $1}');
	$METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/ITS/fasta/${f}.fa $S;
done
```

### Rename files
(should edit fq2fa.pl to name the files correctly...)
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
rename 's/\.trimmed\.1\.fq.fa/_R1.fa/' *1.fq.fa
rename 's/\.trimmed\.2\.fq.fa/_R2.fa/' *2.fq.fa
```
### SSU/58S/LSU removal 

#### Preperation - this is run once only
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


##### Split fasta into chunks for SSU/58S/LSU removal
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
for f in *.fa;
do
  S=$(echo $f|awk -F"." '{print $1}')
    mkdir $S
    split -l 2000 $f -a 3 -d ${S}/$f.
done
```
##### Remove SSU/LSU
Note - creates a file with the paths to all of the split files in each sample directory then submits cluster array job

This will create a large number of array jobs on the cluster

```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta

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
for d in $METAGENOMICS/data/$RUN/ITS/fasta/*R1
do
	 $METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_SSU_58Ss.R $d "*.\\.ssu" "*.\\.58" $d.fa
done

for d in $METAGENOMICS/data/$RUN/ITS/fasta/*R2
do
	 $METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_58Se_LSU.R $d "*.\\.58" "*.\\.lsu" $d.fa
done
```
### Remove empty fastas
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
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
counter=1
for d in $METAGENOMICS/data/$RUN/ITS/fasta/*R[0-9]
do 
S=$(echo $d|awk -F"_" '{print $1}'|awk -F"/" '{print $NF}')
  if (( $counter==1 ))
  then
    $METAGENOMICS/scripts/chimeras.sh $d/ITS1.t.fa $METAGENOMICS/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS1.fasta ${S}.${counter}.cfree.fa $METAGENOMICS/data/$RUN/ITS/de_chimerad/
    counter=2
  else
    $METAGENOMICS/scripts/chimeras.sh $d/ITS2.t.fa $METAGENOMICS/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta ${S}.${counter}.cfree.fa $METAGENOMICS/data/$RUN/ITS/de_chimerad/
    counter=1
  fi
done
```
### Return merged common ITS1 and ITS2, unique ITS1 and unique ITS2
```shell	
cd $METAGENOMICS/data/$RUN/ITS/final

counter=0;
for f in $METAGENOMICS/data/$RUN/ITS/de_chimeraed/*cfree*
do counter=$((counter+1)); 
S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	if (( $counter % 2 == 0 ))
	then
		R2=$f;
		$METAGENOMICS/scripts/catfiles.pl $R1 $R2 $S "nojoin";
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
cd $METAGENOMICS/data/$RUN/ITS/final
cat *.fa > ITS.all.fa
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/$RUN/ITS/final/ITS.all.fa $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/its/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta FALSE
biom summarize-table -i $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom
X=`biom summarize-table -i $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"` 
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus/otu_table_mc2_w_tax.biom $METAGENOMICS/analysis/$RUN/ITS/ITS_all_cdout/ $METAGENOMICS/data/$RUN/map.tsv . $X
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

```shell
cd $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus
Rscript $METAGENOMICS/scripts/analysis.R otu_table_mc2_w_tax.biom colData median ITS.median.csv
```
