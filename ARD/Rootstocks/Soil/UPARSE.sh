# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN FUN 23 21

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN BAC 0 0

# OO
cp $PROJECT_FOLDER/FM/Stoolbeds/OO*  $PROJECT_FOLDER/data/$RUN/.
cp $PROJECT_FOLDER/FM/Stoolbeds/zOO* $PROJECT_FOLDER/data/$RUN/.

# NEM
cp $PROJECT_FOLDER/FM/Stoolbeds/NEM* $PROJECT_FOLDER/data/$RUN/.
cp $PROJECT_FOLDER/FM/Stoolbeds/zNEM* $PROJECT_FOLDER/data/$RUN/.
