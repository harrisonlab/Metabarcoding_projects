# folder containing project files
PROJECT_FOLDER=~/projects/Endophytes

# sequencer run folder 
RUN=leone_combined

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa  # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN FUN 22 21 FUN.otus.fa

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC 17 21 BAC.otus.fa
