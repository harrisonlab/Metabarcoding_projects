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
mkdir $METAGENOMICS/analysis/$RUN
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

### Utax reference databases
Reference databases were downloaded from:
http://drive5.com/usearch/manual/utax_downloads.html
(Unite V7 and RDP trainset 15)
```shell
usearch8.1 -makeudb_utax refdb.fa -output 16s_ref.udb -report 16s_report.txt
usearch8.1 -makeudb_utax refdb.fa -utax_trainlevels kpcofgs ‑utax_splitlevels NVpcofgs -output ITS_ref.udb -report ITS_report.txt
```

### QC
Qualtiy checking was performed with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

From same folder containing fastq files ran:

	fastqc *

### Demultiplexing

#### OLD (fast) method
We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Ambiguous reads are written to two (f & r) seperate files.

Running something like the below should give a good indication of what index_1 and index_2 should be - this is useful if you don't knwo what the primer sequences are and to get a feel of how many mismatches to use. 
```shell
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq > zzexpressions.txt
grep -x "[ATCG]\+" $(ls|head -n1)| cut -c-16|sort|uniq|xargs -I r grep -c ^r $(ls|head -n1) >zzcounts.txt
```

I typically use the first 8 nucleotides of the primer and allow 2 mismatches (the final parameter)
```shell
for f in $METAGENOMICS/data/$RUN/fastq/*_R1_*
do     
	R1=$f     
	R2=$(echo $R1|sed 's/_R1_/_R2_/')    
	S=$(echo $f|awk -F"_" '{print $2}')     
	echo $f    
	$METAGENOMICS/scripts/demulti.sh $R1 $R2 "CCTACGGG" "GACTACHV" "CTTGGTCA" "ATATGCTT" 2  
done   
```
```shell
mv *bacterial* ../16S/fastq/.
mv *fungal* ../ITS/fastq/.
```

## UPARSE 16s workflow

### Pre-processing
Script will join PE reads (and save joined files to unfiltered folder), remove adapter contamination and filter on minimum size and quality threshold.
```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*
do
    R1=$f
    R2=$(echo $R1|sed 's/_R1_/_R2_/')
    S=$(echo $f|awk -F"_" '{print $2}')
    $METAGENOMICS/scripts/16Spre.sh $R1 $R2 ${S}.filtered.fa  $METAGENOMICS/data/$RUN/16S/filtered $METAGENOMICS/primers/adapters.db 300 0.005 ${S}.
done   

```

### OTU fasta creation

##### Concatenate files

```shell
cat $METAGENOMICS/data/$RUN/16S/filtered/*filtered* > $METAGENOMICS/data/$RUN/16S/16S.t.fa
```

##### Truncate and pad
Remove multiplex primers and pad reads to same length.

```shell
X=`cat 16S.t.fa|awk '{if ($1~/>/) {print $0} else {print length($0)};}'|awk '{if ($1~/>/) {y=0} else{y+=$0}};y>x{x=y};END{print x}'`
usearch8.1 -fastx_truncate 16S.t.fa -stripleft 17 -stripright 21 -trunclen $X -padlen $X -fastaout 16S.fa
rm 16S.t.fa
```
Problem with (free version) usearch running out of memory for this and subsequent steps. Cutting and recombining data during dereplication phase gives a fairly unsatisfactory, but working method. 

##### Dereplication 
Required for usearch 8.x otu clustering

```shell
usearch8.1 -derep_fulllength 16S.fa -fastaout 16S.uniques.fasta -sizeout 
usearch8.1 -sortbysize 16S.uniques.fasta -fastaout 16S.sorted.fasta -minsize 2
rm 16S.fa 16S.uniques.fasta
```
get_uniq.pl will give output comparable to derep_fulllength for larger sequence collections
combine_uniq.pl will combine several sets of dereplicated sequences, maintaining the counts.
The sorting algorithm may run out of memory as well - it shouldn't be too difficult to adjust combine_uniq.pl to sort and filter on size (though the cluster algorithm will also filter on min size)

##### Clustering
Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras)
```shell
usearch8.1 -cluster_otus 16S.sorted.fasta -otus 16S.otus.fa -uparseout 16S.out.up -relabel OTU -minsize 2
```
 
##### Assign Taxonomy
The taxa file output by utax is difficult to manipulate in R. Therefore the script mod_taxa.pl should be used to produce an R friendly taxa file.

```shell
usearch8.1 -utax 16S.otus.fa -db $METAGENOMICS/taxonomies/utax/16s_ref.udb -strand both -utaxout 16S.reads.utax -rdpout 16S.rdp -alnout 16S.aln.txt
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl > 16S.taxa
```

### OTU table 

##### Concatenate unfiltered reads
Unfiltered fastq will need to be converted to fasta first 
```shell
for f in $METAGENOMICS/data/$RUN/16S/joined/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/16S/16S.unfiltered.fa $S
done
```	
##### Make table
Creates an OTU table of read counts per OTU per sample
```shell
usearch8.1 -usearch_global 16S.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.97 -biomout 16S.otu_table.biom -otutabout 16S.otu_table.txt
```

### Statistical analysis
Requires analysis2.R and deseq.r

ubiom makes a S3 biom object from the OTU table (16S.otu_table.txt), OTU taxonomy (16S.taxa) and sample description file (colData)
analysis2.R/deseq.r contain scripts to produce deseq objects and run differential analysis + a few graphing options.
	
## ITS workflow

### Pre-processing

#### Remove (and save) reads contain both f & r primers
```shell
for f in $METAGENOMICS/data/$RUN/ITS/fastq/*R1*.fastq; 
do 
    usearch8.1 -search_oligodb $f -db $METAGENOMICS/primers/f_primers.db -strand both -userout ${f}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
done

for f in $METAGENOMICS/data/$RUN/ITS/fastq/*R2*.fastq; 
do 
    usearch8.1 -search_oligodb $f -db $METAGENOMICS/primers/r_primers.db -strand both -userout ${f}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
done
```

```shell
cd $METAGENOMICS/data/$RUN/ITS/fastq
counter=0
for  f in *.fastq.txt
do counter=$((counter+1))
    if (( $counter % 2 == 0 ))
    then
    	R2=$f
    	S1=$(echo $R1|sed 's/.txt//')
    	S2=$(echo $R2|sed 's/.txt//')
	sed 's|^|/|;s|$|/,+3 d|' <(grep p13 $R1|awk -F"\t" '{print $1}') > temp.sed
	sed -f temp.sed $S1 > ${S1}.cleaned.fastq
	sed 's|^|/|;s|$|/,+3 d|' <(grep p14 $R2|awk -F"\t" '{print $1}') > temp.sed
	sed -f temp.sed $S2 > ${S2}.cleaned.fastq	
    fi
    R1=$f
done
#grep -A 3 -F -f <(grep p13 $R1|awk -F"\t" '{print $1}') $S1|grep "\-\-" -v > ${S1}.short.fastq
#grep -A 3 -F -f <(grep p14 $R2|awk -F"\t" '{print $1}') $S2|grep "\-\-" -v > ${S2}.short.fastq 
```

#### Trimming with usearch
utrim is using the expected error per base. The settings below (which also set minimum length to 200) will discard sequences of 200 bases if expected error is > 1 - this is for the forward read only, the reverse read is not as stringent due to (current) fairly poor quality of data. 
Will also save as renamed fasta.
```shell
counter=0;
for f in $METAGENOMICS/data/$RUN/ITS/fastq/*
do counter=$((counter+1))
	S=$(echo $f|awk -F"_" '{print $2}')
	if (( $counter % 2 == 0 ))
	then
		$METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.2.fq $METAGENOMICS/data/$RUN/ITS/trimmed 0.02 200 ${S}.
	else
		$METAGENOMICS/scripts/utrim.sh $f ${S}.trimmed.1.fq $METAGENOMICS/data/$RUN/ITS/trimmed 0.005 200 ${S}
	fi
done
```
but the fasta is multilined - needs to be single lined...
```
for f in *.fa
do
	S=$(echo $f|awk -F"." '{print $1}')
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <$f > ${S}.1.fa
	sed -i -e '1d' ${S}.1.fa
done
```

#### SSU/58S/LSU removal 

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
##### Identify SSU, 5.8S  and LSU regions
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
##### Remove SSU, 5.8S  and LSU regions and merge output
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

##### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2
(script needs editing...)

```shell	
cd $METAGENOMICS/data/$RUN/ITS/....

counter=0;
for f in `ls $METAGENOMICS/data/$RUN/ITS/de_chimeraed/*cfree*| sort -V`
do counter=$((counter+1)); 
S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	if (( $counter % 2 == 0 ))
	then
		R2=$f;
		$METAGENOMICS/scripts/catfiles.pl $R1 $R2 $S;
	fi
	R1=$f
done
```

### UPARSE workflow

##### Pad file (probably not necessary as done in previous step...
```shell
X=`cat ITS1.fa|awk '{if ($1~/>/) {print $0} else {print length($0)};}'|awk '{if ($1~/>/) {y=0} else{y+=$0}};y>x{x=y};END{print x}'`
usearch8.1 -fastx_truncate ITS1.fa -trunclen $X -padlen $X -fastaout ITS1.t.fa
```
##### Dereplication

```shell
usearch8.1 -derep_fulllength ITS1.t.fa -fastaout ITS.uniques.fasta -sizeout
usearch8.1 -sortbysize ITS.uniques.fasta -fastaout ITS.sorted.fasta -minsize 2
rm ITS.uniques.fasta
```
##### Clustering
```shell
usearch8.1 -cluster_otus ITS.sorted.fasta -otus ITS.otus.fa -uparseout ITS.out.up -relabel OTU -minsize 2 
```
##### Assign taxonomy
```shell
usearch8.1 -utax ITS.otus.fa -db $METAGENOMICS/taxonomies/utax/ITS_ref.udb -strand both -utaxout ITS.reads.utax -rdpout ITS.rdp -alnout ITS.aln.txt
cat ITS.rdp|$METAGENOMICS/scripts/mod_taxa.pl > ITS.taxa

```

#### OTU table creation

##### Concatenate unfiltered reads
Unfiltered fastq will need to be converted to fasta first 
```shell
for f in $METAGENOMICS/data/$RUN/ITS/joined/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/ITS/ITS.unfiltered.fa $S
done
```	
##### Make table
Creates an OTU table of read counts per OTU per sample
```shell
usearch8.1 -usearch_global ITS.unfiltered.fa -db ITS.otus.fa -strand plus -id 0.97 -biomout ITS.otu_table.biom -otutabout ITS.otu_table.txt
```

### Statistical analysis
Requires analysis2.R and deseq.r

ubiom makes a S3 biom object from the OTU table (ITS.otu_table.txt), OTU taxonomy (ITS.taxa) and sample description file (colData)
analysis2.R/deseq.r contain scripts to produce deseq objects and run differential analysis + a few graphing options.

The OTU table header is contains a #. To import into R the set comment.char="" in the read.table parameters


##oomycetes
```shell
$METAGENOMICS/scripts/pick_OTU.sh  $METAGENOMICS/data/$RUN/ITS/final/ITS.all.fa $METAGENOMICS/analysis/$RUN/ITS/ITS_all_otus $METAGENOMICS/scripts/params.txt $METAGENOMICS/taxonomies/Silva119/97_18S_only/Silva_119_rep_set97_aligned_18S_only.fna FALSE
```

##Combine samples
Biom table for samples from multiple NGS runs are required.

This will mean the names of each fasta will need to be made unique and the sequence lengths will need to be set to the same.

Concatanate required samples per run. All fastas have common naming format so should be able to change with sed:
```shell
sed -e -i 's/_/_runID_/g' < input file
```

##OLD

### Quiime pipeline 16S Remove chimeras
Downloaded usearch 8.0 and RDP gold reference database from http://drive5.com/usearch/manual/cmd_uchime_ref.html

Ran the 'remove chimeras script'

```shell
#remove chimeras script 	
for f in $METAGENOMICS/data/$RUN/16S/joined/*
do
	S=$(echo $f|awk -F"." '{print $1}')
	$METAGENOMICS/scripts/chimeras.sh $f $METAGENOMICS/taxonomies/RDP_gold.fasta ${S}.cfree.fastq $METAGENOMICS/data/$RUN/16S/de_chimeraed/
done
```

### 16S Convert filtered fastq to fasta

configured utrim.sh to do this

```shell
cd  $METAGENOMICS/data/$RUN/16S/trimmed/	

for f in  *trimmed*
do
 S=$(echo $f|awk -F"." '{print $1}')
 $METAGENOMICS/scripts/fq2fa.pl $f $f.fa $S
 mv $f.fa $METAGENOMICS/data/$RUN/16S/filtered/.
done

cd  $METAGENOMICS/data/$RUN/16S/de_chimeread/	

for f in  *cfree*
do
 S=$(echo $f|awk -F"." '{print $1}')
 $METAGENOMICS/scripts/fq2fa.pl $f $f.fa $S
 mv $f.fa $METAGENOMICS/data/$RUN/16S/unfiltered/.
done

```
Quiime method (usearch 6.x)
```shell
$METAGENOMICS/scripts/pick_OTU.sh   $METAGENOMICS/data/$RUN/16S/16S.fa  $METAGENOMICS/analysis/$RUN/16S/16S_otus $METAGENOMICS/scripts/parameters.txt $PYTHONUSERBASE/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta TRUE
 X=`biom summarize-table -i METAGENOMICS/analysis/$RUN/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom|grep  Min|sed -n "/ Min: */s/ Min: *//p"|sed -n "/\..*/s/\..*//p"`
$METAGENOMICS/scripts/core_diversity.sh $METAGENOMICS/analysis/$RUN/16S/16S_otus/otu_table_mc2_w_tax_no_pynast_failures.biom $METAGENOMICS/analysis/$RUN/16S/16s_cdout/ $METAGENOMICS/data/map.tsv $METAGENOMICS/analysis/$RUN/16S/16S_otus/rep_set.tre $X
```

##### OLD statistical analysis
analysis.R biom_table colData median/geomean outfile  
Requires a file (colData) which describes condition (e.g. infected or uninfected) for each sample 
```shell
cd $METAGENOMICS/analysis/$RUN/16S/16S_otus
Rscript $METAGENOMICS/scripts/analysis.R "otu_table_mc2_w_tax_no_pynast_failures.biom" colData median res.sig.csv
```

### PhiX filtering
Not implemented... For the particular sequencing protocol we don't get much (or any) PhiX contamination. Removal of any contaminants is simple via aligning to the Illumina PhiX genome <ln>http://support.illumina.com/sequencing/sequencing_software/igenome.html </ln> Bowtie2 method implemented here.

NOTE - the below scipts that implement something like 'for f in *' are dependent on the naming convention of the samples. For instance something like s1.1.fq - s20.2.fq will loop through the files in the order  s1.1.fq, s12.1.fq, s12.2.fq, s1.2.fq, which is clearly not what is wanted.
'for f in `ls *| sort -V`' will do a natural sort of the files which should fix any problems - or use a different sample naming convention (e.g. s001. - sxxx.)

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

```shell
 counter=0
 for f in $METAGENOMICS/data/$RUN/fastq/*R1*.fastq
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

####adapter remove
```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*.fastq; 
do 
	usearch8.1 -search_oligodb $f -db $METAGENOMICS/primers/adapters.db -strand both -userout ${f}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
done

for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*.fastq
do
	R1=$f
	R2=$(echo $R1|sed 's/_R1_/_R2_/')
	S=$(echo $f|awk -F"_" '{print $2}')
    	cat ${R1}.txt ${R2}.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' $R1
    	cat ${R1}.txt ${R2}.txt|awk -F"\t" '{print $1}'|sort|uniq|xargs -I ¬ sed -i -ne:t -e"/*\@¬.*/D" -e'$!N;//D;/'"\@¬/{" -e"s/\n/&/3;t" -e'$q;bt' -e\} -e's/\n/&/'"1;tP" -e'$!bt' -e:P  -e'P;D' $R2
done
    
    
```
### Demultiplexing

##### New (slow) method
This is going to be edited to use usearch8.1 search_oligodb - the algorithm used accepts mismatches at multiple positions.

```
for f in *.fastq; 
do 
	usearch8.1 -search_oligodb $f -db ../../../scripts/primers.db -strand both -userout ${f}.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
done
```
The bit below is a bit rubbish - working on a speed improvement
The search_oligodb part will also identify adapter contamination .
```shell
#this is slow as a slow thing (about 2 minutes per sample! - the old method was roughly 100 times faster)
counter=0
for f in *.txt
 do counter=$((counter+1))
 if (( $counter % 2 == 0 ))
      then
        R2=$f
        S1=$(echo $R1|sed 's/\.txt//')
        S2=$(echo $R2|sed 's/\.txt//')
	grep -F -f <(awk -F"\t" '{print $1}'<$R1) $R2 > output.txt
	grep -A 3 -F -f <(awk -F"\t" '{if ($2=="p13") print $1}' <output.txt)  $S1 > ${S1}.bacterial.fq
	grep -A 3 -F -f <(awk -F"\t" '{if ($2=="p13") print $1}' <output.txt)  $S2 > ${S2}.bacterial.fq
	grep -A 3 -F -f <(awk -F"\t" '{if ($2=="p11") print $1}' <output.txt)  $S1 > ${S1}.fungal.fq
	grep -A 3 -F -f <(awk -F"\t" '{if ($2=="p11") print $1}' <output.txt)  $S2 > ${S2}.fungal.fq
 fi
 R1=$f
done
```
### Join PE reads
(do not filter at this stage - unfiltered joined reads are required for later stage)

```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*
do
	R1=$f
	R2=$(echo $R1|sed 's/_R1_/_R2_/')
	S=$(echo $f|awk -F"_" '{print $2}')
	$METAGENOMICS/scripts/ujoin.sh $R1 $R2 ${S}.joined.fastq $METAGENOMICS/data/$RUN/16S/joined
done
```

#### Filter fastq files

##### Remove adapters
```shell
 for f in $METAGENOMICS/data/$RUN/16S/filtered/*.filtered; 
 do 
 	$METAGENOMICS/scripts/filtadapt.sh $f $METAGENOMICS/primers/adapters.db
 done
```

##### Quality trimming
```shell
for f in $METAGENOMICS/data/$RUN/16S/joined/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/utrim.sh $f ${S}.filtered.fa $METAGENOMICS/data/$RUN/16S/filtered 0.005 300 ${S}.
done
```

#### Rename sequences
The sequence renaming of utrim is not working correctly (not unique). The below will produce unique sequence names per sample
```shell
for f in $METAGENOMICS/data/$RUN/16S/filtered/*.filtered
do	
	S=$(echo $f|awk -F"/" '{print $NF}')
	cat $f|awk '{if ($0~/>/){y+=1;gsub(/_.*/,"."y,$0)};print $0}' > ${S}.fa
	rm $f
done
```
### Convert to unpaired fasta files

Alternative method
```shell
cd $METAGENOMICS/data/$RUN/ITS/trimmed

for f in *trimmed*;
do
	S=$(echo $f|awk -F"." '{print $1}');
	$METAGENOMICS/scripts/fq2fa.pl $f $METAGENOMICS/data/$RUN/ITS/fasta/${f}.fa $S;
done
```

This might have some use, but can't remeber what - alternative simpler method works...
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
### Rename files
(should edit fq2fa.pl to name the files correctly...)
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
rename 's/\.trimmed\.1\.fq.fa/_R1.fa/' *1.fq.fa
rename 's/\.trimmed\.2\.fq.fa/_R2.fa/' *2.fq.fa
```
### Remove chimeras
NOT NECESSARY
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
    $METAGENOMICS/scripts/chimeras.sh $d/ITS2.fa $METAGENOMICS/taxonomies/uchime_sh_refs_dynamic_develop_985_11.03.2015.ITS2.fasta ${S}.${counter}.cfree.fa $METAGENOMICS/data/$RUN/ITS/de_chimerad/
    counter=1
  fi
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
###### Remove empty fastas - now incorporated into above step
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
###### Pad files - in above step (though ITS1 and ITS2 have different pad lengths...)
uclust performs better if FASTAs are same length.

Example (of padding):
```shell
X=`grep ">" -v S13_R1.fa|awk '{ print length($0); }'|awk '$0>x{x=$0};END{print x}'`
cat S13_R1.fa| sed -e :a -e "s/^[^>].\{1,`expr $X - 1`\}$/&N/;ta"
```

This could also be done in R as well (in the merge bit) - this would also remove the empty fasta files...
```Rscript
ITS <- ITS[ITS@ranges@width>0]
ITS <- stackStrings(ITS,0,max(ITS@ranges@width),Lpadding.letter="N",Rpadding.letter="N")
ITS <- subseq(ITS,start=2,width = (max(ITS@ranges@width)-1))
writeXStringSet(ITS,"ITS.t.fa")
```

cat ITS.x.fa
