# 16s workflow

## Pre-processing
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
## UPARSE

### Cluster and assign taxonomy
Problem with (free version) usearch running out of memory for dereplication and subsequent steps. Cutting and recombining data during dereplication phase gives a fairly unsatisfactory, but working method. combine_uniq.pl will combine several sets of dereplicated sequences, maintaining the counts.
The sorting algorithm may run out of memory as well - it shouldn't be too difficult to adjust combine_uniq.pl to sort and filter on size (though the cluster algorithm will also filter on min size)

get_uniq.pl will give output comparable to derep_fulllength and sortbysize for larger sequence collections. get_uniq.pl requires unformatted fasta (as in sequence not split every 80 nucleotides). It also sorts and removes singletons. 

... Need to do some testing on speed as I don't think there is much difference between usearch and get_uniq.pl.

For 1,000,000 reads

|	|derep|sortbysize|get_uniq.pl|
|---|---|---|---|
|real|0m38.037s|0m15.242s|0m16.738s|
|user|0m18.361s|0m5.236s|0m9.585s|
|sys|0m0.744s|0m0.304s|0m1.124s|

I've updated the below to use get_unique rather than usearch.    


The taxa file output by utax is difficult to manipulate in R. Therefore the script mod_taxa.pl should be used to produce an R friendly taxa file.

```shell
#### Concatenate files
cat $METAGENOMICS/data/$RUN/16S/filtered/*filtered* > $METAGENOMICS/data/$RUN/16S.t.fa
#### Truncate and pad (Remove multiplex primers and pad reads to same length.)
X=`cat 16S.t.fa|awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}'|awk '$0>x{x=$0};END{print x}'`
usearch8.1 -fastx_truncate 16S.t.fa -stripleft 17 -stripright 21 -trunclen $X -padlen $X -fastaout 16S.fa
rm 16S.t.fa
#### Dereplication
#usearch8.1 -derep_fulllength 16S.fa -fastaout 16S.uniques.fasta -sizeout 
#usearch8.1 -sortbysize 16S.uniques.fasta -fastaout 16S.sorted.fasta -minsize 2
#rm 16S.fa 16S.uniques.fasta
cat 16S.fa|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|$METAGENOMICS/scripts/get_uniq.pl > 16S.sorted.fasta 
rm 16S.fa
#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))
usearch8.1 -cluster_otus 16S.sorted.fasta -otus 16S.otus.fa -uparseout 16S.out.up -relabel OTU -minsize 2
#### Assign Taxonomy
usearch8.1 -utax 16S.otus.fa -db $METAGENOMICS/taxonomies/utax/16s_ref.udb -strand both -utaxout 16S.reads.utax -rdpout 16S.rdp -alnout 16S.aln.txt
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl > 16S.taxa
```
### OTU table 

fq2fa_v2.pl will convert fastq to fasta and trim left and right ends of reads

usearch_global may run out of memory as well. 16S.unfiltered.fa can be split into chunks (reads from a single sample must not be split across chunks) and seperate OTU tables created for each chunk. It is trivial to recombine these in R.

fq2fa_v2.pl could be replaced with a similar awk script as per ITS - will save a couple of minutes.
```shell
#### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first )
for f in $METAGENOMICS/data/$RUN/16S/unfiltered/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $f $METAGENOMICS/data/$RUN/16S.unfiltered.fa $S 17 21
done
#### Make table (Creates an OTU table of read counts per OTU per sample)
usearch8.1 -usearch_global 16S.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.97 -biomout 16S.otu_table.biom -otutabout 16S.otu_table.txt
```

## Statistical analysis


Requires analysis2.R and deseq.r

ubiom makes a S3 biom object from the OTU table (16S.otu_table.txt), OTU taxonomy (16S.taxa) and sample description file (colData)
analysis2.R/deseq.r contain scripts to produce deseq objects and run differential analysis + a few graphing options.

The OTU table header starts with a hash. To import into R set comment.char="" in the read.table parameters
