# Create dostance matirix from amplicon fatsa
usearch -calc_distmx 16S_Endo_otus.fa -distmxout 16S_Endo.phy -format phylip_square
usearch -calc_distmx 16S_Rhiz_otus.fa -distmxout 16S_Rhiz.phy -format phylip_square
usearch -calc_distmx ITS_otus.fa -distmxout ITS.phy -format phylip_square
