# setup links and variables
PROJECT_FOLDER=~/projects/ARD/analysis
ln -s ~/pipelines/metabarcoding/ $PROJECT_FOLDER/metabarcoding_pipeline
RUN="HG"
mkdir ~/projects/ARD/analysis/data # I've hardcoded some stuff in the pipeline to use a folder called data (need to edit this...) 


# ITS
# cluster
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN ITS 0 0
# taxonomy
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN ITS
# OTU distance matrix
usearch8.1 -calc_distmx ITS.zotus.fa -distmxout ITS.phy -distmo fractdiff -format phylip_square
# Create OTU table (-c OTUS  32bit USEARCH large data version of -c OTU)
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN ITS 23 21

# 16S
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 17 21
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
usearch8.1 -calc_distmx 16S.zotus.fa -distmxout 16S.phy -distmo fractdiff -format phylip_square
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC 17 21


