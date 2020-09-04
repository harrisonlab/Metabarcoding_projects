# folder containing project files
PROJECT_FOLDER=~/projects/Eggplant

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa  # workaround for uparse bug
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN FUN FUN.otus.fa 22 21

# BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0
sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC BAC.otus.fa 17 21 

#### ADDITIONAL CODE ####


ln -s ~/projects/Eggplant/data/uganda_root/BAC/unfiltered/* .
rename 's/^/RT/' *
mv * ../tanzania_uganda_all/BAC/unfiltered/.

ln -s ~/projects/Eggplant/data/uganda_root/BAC/filtered/* .
rename 's/^/RT/' *
mv * ../tanzania_uganda_all/BAC/filtered/.

ln -s ~/projects/Eggplant/data/uganda_root/FUN/unfiltered/* .
rename 's/^/RT/' *
  mv * ../tanzania_uganda_all/FUN/unfiltered/.

ln -s ~/projects/Eggplant/data/uganda_root/FUN/filtered/* .
rename 's/^/RT/' *
  mv * ../tanzania_uganda_all/FUN/filtered/.


ln -s ~/projects/Eggplant/data/uganda_soil_stem/BAC/unfiltered/* .
mv * ../tanzania_uganda_all/BAC/unfiltered/.

ln -s ~/projects/Eggplant/data/uganda_soil_stem/BAC/filtered/* .
mv * ../tanzania_uganda_all/BAC/filtered/.

ln -s ~/projects/Eggplant/data/uganda_soil_stem/FUN/unfiltered/* .
mv * ../tanzania_uganda_all/FUN/unfiltered/.

ln -s ~/projects/Eggplant/data/uganda_soil_stem/FUN/filtered/* .
mv * ../tanzania_uganda_all/FUN/filtered/.


ln -s ~/projects/Eggplant/data/tanzania/BAC/unfiltered/* .
mv * ../tanzania_uganda_all/BAC/unfiltered/.

ln -s ~/projects/Eggplant/data/tanzania/BAC/filtered/* .
mv * ../tanzania_uganda_all/BAC/filtered/.

ln -s ~/projects/Eggplant/data/tanzania/FUN/unfiltered/* .
mv * ../tanzania_uganda_all/FUN/unfiltered/.

ln -s ~/projects/Eggplant/data/tanzania/FUN/filtered/* .
mv * ../tanzania_uganda_all/FUN/filtered/.

ls ../tanzania_uganda_all/FUN/filtered/*
ls ../tanzania_uganda_all/FUN/unfiltered/*
ls ../tanzania_uganda_all/BAC/filtered/*
ls ../tanzania_uganda_all/BAC/unfiltered/*
mkdir -p  root_stem/BAC/filtered
mkdir -p  root_stem/BAC/unfiltered
mkdir -p  root_stem/FUN/filtered
mkdir -p  root_stem/FUN/unfiltered

 cp -av tanzania_uganda_all/FUN/filtered/* root_stem/FUN/filtered/.
 cp -av tanzania_uganda_all/FUN/unfiltered/* root_stem/FUN/unfiltered/.
 cp -av tanzania_uganda_all/BAC/filtered/* root_stem/BAC/filtered/.
 cp -av tanzania_uganda_all/BAC/unfiltered/* root_stem/BAC/unfiltered/.


# folder containing project files
PROJECT_FOLDER=~/projects/Eggplant

RUN=root_stem

# FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN FUN 0 0
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN BAC 0 0

sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/FUN.zotus.fa  # workaround for uparse bug
#$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN FUN
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN FUN 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN FUN FUN.otus.fa 22 21

# BAC

sed -i -e 's/Zotu/OTU/' $PROJECT_FOLDER/data/$RUN/BAC.zotus.fa
#$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c dist $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN BAC
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN BAC BAC.otus.fa 17 21 


cp ../../data/root_stem/BAC.otu_table.txt .
cp ../../data/root_stem/BAC.utax.taxa .
cp ../../data/root_stem/BAC.sintax.taxa .
cp ../../data/root_stem/BAC.otus.fa .

cp ../../data/root_stem/FUN.otu_table.txt .
cp ../../data/root_stem/FUN.utax.taxa .
cp ../../data/root_stem/FUN.sintax.taxa .
cp ../../data/root_stem/FUN.otus.fa .



