# calculate distance matrix from amplicon fasta
usearch -calc_distmx 16S_otus.fa -distmxout 16S.phy -format phylip_square 
usearch -calc_distmx ITS_otus.fa -distmxout ITS.phy -format phylip_square 

# convert usearch taxonomy utax format to taxa format
~/pipelines/metabarcoding/scripts/mod_taxa_sintax.pl 16S_reads.utax > 16S.taxa
~/pipelines/metabarcoding/scripts/mod_taxa_sintax.pl ITS_reads.utax > ITS.taxa

# species level taxonomy have no values (16S doesn't have any species)
sed -i -e 's/,$/,0/' ITS.taxa
