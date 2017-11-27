# setup links and variables
PROJECT_FOLDER=~/projects/ARD/analysis
ln -s ~/pipelines/metabarcoding/ $PROJECT_FOLDER/metabarcoding_pipeline
RUN="Koekoek"
mkdir ~/projects/ARD/analysis/data # I've hardcoded some stuff in the pipeline to use a folder called data (need to edit this...) 

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER/ ../$RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/$RUN/FUN.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER ../$RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER ../$RUN FUN 23 21

for F in  $PROJECT_FOLDER/$RUN/FUN*otus.fa; do
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS_regions \
  $F \
  $PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
  $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
  20
done

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER ../$RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/$RUN/BAC.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER ../$RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER ../$RUN BAC 17 21
