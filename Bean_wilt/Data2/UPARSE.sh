PROJECT_FOLDER=~/projects/Kenyan_beans/analysis
ln -s ~/pipelines/metabarcoding/ $PROJECT_FOLDER/metabarcoding_pipeline
RUN="../data2"
mkdir $PROJECT_FOLDER/data

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa  # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN FUN 23 21

# OO
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN OO 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/OO.zotus.fa 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN OO sintax
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN OO 21 20 
