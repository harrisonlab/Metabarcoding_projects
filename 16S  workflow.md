# 16s workflow

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size percent_diff max_errrors 

```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*.fastq
do
    R1=$f
    R2=$(echo $R1|sed 's/_R1_/_R2_/')
    S=$(echo $f|awk -F"_" -v D=$RUN '{print $2"D"D}')
    $METAGENOMICS/scripts/16Spre.sh $R1 $R2 $S  $METAGENOMICS/data/$RUN/16S/filtered $METAGENOMICS/primers/adapters.db 300 15 1 
done   

```
## UPARSE

### Cluster and assign taxonomy
Problem with (free version) usearch running out of memory for dereplication and subsequent steps. combine_uniq.pl will combine several sets of dereplicated sequences, maintaining the counts.
The sorting algorithm may run out of memory as well - it shouldn't be too difficult to adjust combine_uniq.pl to sort and filter on size (though the cluster algorithm will also filter on min size)

get_uniq.pl will give output comparable to derep_fulllength and sortbysize for larger sequence collections. get_uniq.pl requires unformatted fasta (as in sequence not split every 80 nucleotides). It also sorts and removes singletons. Performance wise this is slightly faster than the usearch method (though both only take minutes for several gig of data) 

The taxa file output by utax is difficult to manipulate in R. Therefore the script mod_taxa.pl should be used to produce an R friendly taxa file.

```shell
#### Concatenate files
cat $METAGENOMICS/data/$RUN/16S/filtered/*filtered* > $METAGENOMICS/data/$RUN/16S.t.fa
#### Truncate and pad (Remove multiplex primers and pad reads to same length.)
X=`cat 16S.t.fa|awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}'|awk '$0>x{x=$0};END{print x}'`
usearch8.1 -fastx_truncate 16S.t.fa -stripleft 17 -stripright 21 -trunclen $X -padlen $X -fastaout 16S.fa
rm 16S.t.fa
#### Dereplication
cat 16S.fa|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|$METAGENOMICS/scripts/get_uniq.pl > 16S.sorted.fasta 
rm 16S.fa
#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))
usearch8.1 -cluster_otus 16S.sorted.fasta -otus 16S.otus.fa -uparseout 16S.out.up -relabel OTU -minsize 2
#### Assign Taxonomy
usearch8.1 -utax 16S.otus.fa -db $METAGENOMICS/taxonomies/utax/16s_ref.udb -strand both -utaxout 16S.reads.utax -rdpout 16S.rdp -alnout 16S.aln.txt
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl > 16S.taxa


### Not implemented (but keep incase something breaks)
#usearch8.1 -derep_fulllength 16S.fa -fastaout 16S.uniques.fasta -sizeout 
#usearch8.1 -sortbysize 16S.uniques.fasta -fastaout 16S.sorted.fasta -minsize 2
#rm 16S.fa 16S.uniques.fasta
```

There's a new version of usearch (v9) which has a different clustring step (denoising). However, it is a bit slower than cluster_otus (takes about 1.5hrs for 50meg derelicated file rather than about 5 minutes on one of our servers). There isn't a paper for this, so I don't know how it works internally - if each entry is independent there's no problem in splitting and running multiple instances as an array job (well with the free version of usearch anyway).

 - doesn't work

``` shell
#mkdir -p temp 
#split -l 2000 16S.sorted.fasta -a 4 -d temp/xx.
#TASKS=`ls temp|wc -l`
#cd temp 
#find . $PWD -name 'xx*' >split_files.txt
#qsub -t 1-$TASKS:1 $METAGENOMICS/scripts/submit_denoise.sh 16S.denoised  
usearch9 -unoise 16S.sorted.fasta -tabbedout out.txt -fastaout 16S.denoised.fa
```



### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
```shell
usearch8 -calc_distmx 16S.otus.fa -distmxout 16S.phy -distmo fractdiff -format phylip_square
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

Occasionally, due to v.poor reverse read quality, joining of f+r reads fails for the vast majority. The following will cluster f+r reads separately and then merge read counts which align to the same OTU. I've dropped the clustering down to 0.95 similarity - both reads aligning to the same OTU at this similarity, I'd suggest is pretty good evidence they're the same. 
I've also added a rev compliment routine to fq2fa_v2.pl, means the reverse reads can be called as plus strand by usearch_global.

```shell
for f in $METAGENOMICS/data/$RUN/16S/fastq/*R1*
do
	R1=$f
	R2=$(echo $R1|sed 's/\_R1_/\_R2_/')
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $R1 $METAGENOMICS/data/$RUN/16S.r1.unfiltered.fa $S 17 0
	$METAGENOMICS/scripts/fq2fa_v2.pl $R2 $METAGENOMICS/data/$RUN/16S.r2.unfiltered.fa $S 21 30 rev
done
usearch8.1 -usearch_global 16S.r1.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r1.txt -userfields query+target+id
usearch8.1 -usearch_global 16S.r2.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r2.txt -userfields query+target+id
$METAGENOMICS/scripts/merge_hits.sh $METAGENOMICS/scripts/merge_hits.R hits.r1.txt hits.r2.txt 16S.otu_table.txt
$METAGENOMICS/scripts/otu_to_biom.pl row_biom col_biom data_biom >16S.otu_table.biom
rm row_biom col_biom data_biom
```

###[ITS workflow](../master//ITS%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)
