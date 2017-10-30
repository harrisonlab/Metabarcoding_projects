# cluster
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU $FPL $RPL

# taxonomy
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU

# OTU distance matrix
usearch8.1 -calc_distmx 16S.otus.fa -distmxout 16S.phy -distmo fractdiff -format phylip_square

# Create OTU table (-c OTUS - 32bit large data version of -c OTU)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN $SSU $FPL $RPL
