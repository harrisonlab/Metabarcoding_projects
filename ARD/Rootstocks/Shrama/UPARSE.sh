PROJECT_FOLDER=~/projects/ARD/analysis
RUN="../Rootstock/Shrama"

ln -s $PROJECT_FOLDER/Rootstocks/Soil/FUN/filtered/Sch-* $PROJECT_FOLDER/RUN/FUN/filtered/.
ln -s $PROJECT_FOLDER/Rootstocks/Soil/FUN/unfiltered/Sch-* $PROJECT_FOLDER/RUN/FUN/unfiltered/.
ln -s $PROJECT_FOLDERs/Rootstocks/Soil/BAC/filtered/Sch-* $PROJECT_FOLDER/RUN/BAC/filtered/.
ln -s $PROJECT_FOLDER/Rootstocks/Soil/BAC/unfiltered/Sch-* $PROJECT_FOLDER/RUN/BAC/unfiltered/.
ln -s $PROJECT_FOLDER/Rootstocks/Rhizosphere/BAC/filtered/Sch-* $PROJECT_FOLDER/RUN/BAC/filtered/.
ln -s $PROJECT_FOLDER/Rootstocks/Rhizosphere/BAC/filtered/Sch-* $PROJECT_FOLDER/RUN/BAC/unfiltered/.
ln -s $PROJECT_FOLDER/Rootstocks/Rhizosphere/FUN/filtered/Sch-* $PROJECT_FOLDER/RUN/FUN/filtered/.
ln -s $PROJECT_FOLDER/Rootstocks/Rhizosphere/FUN/unfiltered/Sch-* $PROJECT_FOLDER/RUN/FUN/unfiltered/.


# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN FUN 23 21

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
#$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN BAC 0 0
# file size too big for 32 bit version - using alternative method
PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC 0 0

