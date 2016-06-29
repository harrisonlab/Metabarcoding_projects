# apple_replant
Metagenomic study of apple replant disease  
 1. HMM Preperation for ITS analysis  
 2. Common workflow
  3. Utax reference databases  
  4. QC
  5. Demultiplexing
 6. 16S workflow
  7. Pre-processing 
  8. UPARSE
    9. Concatenate files  
    10. Truncate and pad  
    11. Dereplication
    12. Clustering
    13. Assign Taxonomy
   14. OTU table
    15. Concatenate unfiltered reads
    16. Make table
  17. Statistical analysis
 18. ITS workflow
  19. Pre-processing
   20. SSU/58S/LSU removal
    21. Split fasta into chunks
    22. Create lists of file paths to chunks
    23. Identify SSU, 5.8S and LSU regions
    24. Remove SSU, 5.8S and LSU regions and merge output
    25. Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2
  26. UPARSE
    27. Pad files
    28. Dereplication
    29. Clustering
    30. Assign taxonomy
   31. OTU table creation
    32. Concatenate unfiltered reads
    33. Make table
  34. Statistical analysis
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
## Common workflow

```shell
#run for each analysis
mkdir -p $METAGENOMICS/analysis/$RUN/16S
mkdir $METAGENOMICS/analysis/$RUN/ITS	
mkdir -p $METAGENOMICS/data/$RUN/fastq
mkdir -p $METAGENOMICS/data/$RUN/16S/fastq
mkdir $METAGENOMICS/data/$RUN/16S/filtered
mkdir $METAGENOMICS/data/$RUN/16S/unfiltered
mkdir -p $METAGENOMICS/data/$RUN/ITS/fastq
mkdir $METAGENOMICS/data/$RUN/ITS/filtered
mkdir $METAGENOMICS/data/$RUN/ITS/unfilterd

```	
The $METAGENOMICS directory should be set to something appropriate (e.g. /home/bob/metagenomics) and $RUN to the name of the NGS run. The realtive path is used in the scripts below - depending on your config you may have to specify full paths.	


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
mkdir -p $METAGENOMICS/data/$RUN/16S/fastq
mkdir -p $METAGENOMICS/data/$RUN/ITS/fastq
mkdir -p $METAGENOMICS/data/$RUN/ambiguous

mv $METAGENOMICS/data/$RUN/fastq/*bacterial* $METAGENOMICS/data/$RUN/16S/fastq/.
mv $METAGENOMICS/data/$RUN/fastq/*fungal* $METAGENOMICS/data/$RUN/ITS/fastq/.
mv $METAGENOMICS/data/$RUN/fastq/*ambig* $METAGENOMICS/data/$RUN/ambiguous/.
```

## 16s workflow

### Pre-processing
Script will join PE reads (and save joined files to unfiltered folder), remove adapter contamination and filter on minimum size and quality threshold.
```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*
do
    R1=$f
    R2=$(echo $R1|sed 's/_R1_/_R2_/')
    S=$(echo $f|awk -F"_" -v D=$RUN '{print $2"D"D}')
    $METAGENOMICS/scripts/16Spre.sh $R1 $R2 $S  $METAGENOMICS/data/$RUN/16S/filtered $METAGENOMICS/primers/adapters.db 300 1 
done   

```

### UPARSE

#### Cluster and assign taxonomy
Problem with (free version) usearch running out of memory for dereplication and subsequent steps. Cutting and recombining data during dereplication phase gives a fairly unsatisfactory, but working method. combine_uniq.pl will combine several sets of dereplicated sequences, maintaining the counts.
The sorting algorithm may run out of memory as well - it shouldn't be too difficult to adjust combine_uniq.pl to sort and filter on size (though the cluster algorithm will also filter on min size)

get_uniq.pl will give output comparable to derep_fulllength and sortbysize for larger sequence collections. get_uniq.pl requires unformatted fasta (as in sequence not split every 80 nucleotides). It also sorts and removes singletons. It is remmed out in the script as low mem version.
... Need to do some testing on speed as I don't think there is much difference between usearch and get_uniq.pl.

The taxa file output by utax is difficult to manipulate in R. Therefore the script mod_taxa.pl should be used to produce an R friendly taxa file.

```shell
##### Concatenate files
cat $METAGENOMICS/data/$RUN/16S/filtered/*filtered* > $METAGENOMICS/data/$RUN/16S.t.fa
##### Truncate and pad (Remove multiplex primers and pad reads to same length.)
X=`cat 16S.t.fa|awk '{if ($1!~/>/) {print length($0)};}'|awk '$0>x{x=$0};END{print x}'`
usearch8.1 -fastx_truncate 16S.t.fa -stripleft 17 -stripright 21 -trunclen $X -padlen $X -fastaout 16S.fa
rm 16S.t.fa
##### Dereplication
usearch8.1 -derep_fulllength 16S.fa -fastaout 16S.uniques.fasta -sizeout 
usearch8.1 -sortbysize 16S.uniques.fasta -fastaout 16S.sorted.fasta -minsize 2
rm 16S.fa 16S.uniques.fasta
#low mem dereplication 
#cat 16S.fa|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|$METAGENOMICS/scripts/get_uniq.pl > 16S.sorted.fasta 
#rm 16S.fa
##### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))
usearch8.1 -cluster_otus 16S.sorted.fasta -otus 16S.otus.fa -uparseout 16S.out.up -relabel OTU -minsize 2
##### Assign Taxonomy
usearch8.1 -utax 16S.otus.fa -db $METAGENOMICS/taxonomies/utax/16s_ref.udb -strand both -utaxout 16S.reads.utax -rdpout 16S.rdp -alnout 16S.aln.txt
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl > 16S.taxa
```
#### OTU table 

fq2fa_v2.pl will convert fastq to fasta and trim left and right ends of reads

```shell
##### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first )
for f in $METAGENOMICS/data/$RUN/16S/unfiltered/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $f $METAGENOMICS/data/$RUN/16S/16S.unfiltered.fa $S 17 21
done
##### Make table (Creates an OTU table of read counts per OTU per sample)
usearch8.1 -usearch_global 16S.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.97 -biomout 16S.otu_table.biom -otutabout 16S.otu_table.txt
```

### Statistical analysis
Requires analysis2.R and deseq.r

ubiom makes a S3 biom object from the OTU table (16S.otu_table.txt), OTU taxonomy (16S.taxa) and sample description file (colData)
analysis2.R/deseq.r contain scripts to produce deseq objects and run differential analysis + a few graphing options.
	
## ITS workflow

### Pre-processing
Script will:
1. Remove reads with both forward and reverse primers
2. Remove reads with adapter contamination
3. Filter for quality and minimum length (with UTRIM)
4. Convert FASTQ to single line FASTA

```shell
for f in $METAGENOMICS/data/$RUN/ITS/fastq/*R1*; 
do     
	R1=$f;     
	R2=$(echo $R1|sed 's/_R1_/_R2_/');     
	S=$(echo $f|awk -F"_" -v D=$RUN '{print $2"D"D}');
	$METAGENOMICS/scripts/ITSpre.sh $R1 $R2 $S  $METAGENOMICS/data/$RUN/ITS/fasta $METAGENOMICS/primers/primers.db 200 200 1; 
done
```

#### SSU/58S/LSU removal 

##### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster
```shell
cd $METAGENOMICS/data/$RUN/ITS/fasta
counter=0
for f in *.fa;
do counter=$((counter+1))
    d=$(echo $f|awk -F"." '{print $1}')
    mkdir $d
    split -l 2000 $f -a 3 -d ${d}/$f.
    cd $d
    find $PWD -name '*.fa.*' >split_files.txt
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

```shell
for d in $METAGENOMICS/data/$RUN/ITS/fasta/*R1
do
	S=$(echo $d|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}');
	$METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_SSU_58Ss.R $d "*.\\.ssu" "*.\\.58" $d.fa $S
done
```
If reverse read quality was poor and it was necessary to truncate reads to get more than a couple of reads past set LOWQUAL to TRUE

LOWQUAL keeps reads which lack 5.8S homology - this is necessary as trimming will in most instances have removed the homologous region. 

```shell
LOWQUAL=TRUE
for d in $METAGENOMICS/data/$RUN/ITS/fasta/*R2
do
	S=$(echo $d|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}');
	$METAGENOMICS/scripts/ITS.sh $METAGENOMICS/scripts/rm_58Se_LSU_v2.R $d "*.\\.58" "*.\\.lsu" $d.fa $S $LOWQUAL
done
```

##### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
mkdir -p $METAGENOMICS/data/$RUN/ITS/filtered
find $METAGENOMICS/data/$RUN/ITS/fasta -type f -name *.r*|xargs -I myfile mv myfile $METAGENOMICS/data/$RUN/ITS/filtered/.

cd $METAGENOMICS/data/$RUN/ITS/filtered
for f in $METAGENOMICS/data/$RUN/ITS/filtered/*r1.fa
do
    R1=$f
    R2=$(echo $R1|sed 's/\.r1\.fa/\.r2\.fa/')
    S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
    $METAGENOMICS/scripts/catfiles_v2.pl $R1 $R2 $S;
done

mkdir R1
mkdir R2
mv *r1* R1/.
mv *r2* R2/.
```

### UPARSE

#### Cluster and assign taxonomy
```shell
##### Concatenate
cat $METAGENOMICS/data/$RUN/ITS/filtered/*.fa > $METAGENOMICS/data/$RUN/ITS.t.fa
##### Pad
X=`cat ITS.t.fa|awk '{if ($1!~/>/) {print length($0)};}'|awk '$0>x{x=$0};END{print x}'`
usearch8.1 -fastx_truncate ITS.t.fa -trunclen $X -padlen $X -fastaout ITS.fa
rm ITS.t.fa
##### Dereplicate
usearch8.1 -derep_fulllength ITS.fa -fastaout ITS.uniques.fasta -sizeout
usearch8.1 -sortbysize ITS.uniques.fasta -fastaout ITS.sorted.fasta -minsize 2
rm ITS.fa ITS.uniques.fasta
##### Cluster
usearch8.1 -cluster_otus ITS.sorted.fasta -otus ITS.otus.fa -uparseout ITS.out.up -relabel OTU -minsize 2 
##### Taxonomy
usearch8.1 -utax ITS.otus.fa -db $METAGENOMICS/taxonomies/utax/ITS_ref.udb -strand both -utaxout ITS.reads.utax -rdpout ITS.rdp -alnout ITS.aln.txt
cat ITS.rdp|$METAGENOMICS/scripts/mod_taxa.pl > ITS.taxa
```

#### OTU table creation

```shell
##### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first)
for f in $METAGENOMICS/data/$RUN/ITS/unfiltered/*.r1.*
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $f $METAGENOMICS/data/$RUN/ITS/ITS.unfiltered.fa $S 22 0
done

for f in $METAGENOMICS/data/$RUN/ITS/unfiltered/*.r2.*
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $f $METAGENOMICS/data/$RUN/ITS/ITS.unfiltered.fa $S 20 0
done
	
##### Make table (creates an OTU table of read counts per OTU per sample)
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
