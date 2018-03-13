# setup links and variables
PROJECT_FOLDER=~/projects/ARD/metabarcoding/analysis
ln -s ~/pipelines/metabarcoding/ $PROJECT_FOLDER/metabarcoding_pipeline
RUN="../WP3"
mkdir ~/projects/ARD/analysis/data # I've hardcoded some stuff in the pipeline to use a folder called data (need to edit this...) 

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa  # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN FUN 23 21

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC 17 21

# OO
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN OO 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/OO.zotus.fa 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN NEM
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN OO sintax
# data size too large for 32 bit usearch - unable to calculate denoised otu table
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN OO 21 20 

# NEM
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN NEM 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/NEM.zotus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN NEM
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN NEM sintax
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN NEM 23 18
