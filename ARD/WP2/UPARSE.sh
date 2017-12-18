PROJECT_FOLDER=~/projects/ARD/analysis
RUN="../WP2"

# Fungi
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' FUN.zotus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN FUN $FPL $RPL
