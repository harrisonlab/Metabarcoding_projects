# set variables
PROJECT_FOLDER=~/projects/Endophytes
RUN=COMBINED

### FUNGAL OTUS ###
SSU=FUN

# forward primer length
FPL=23 

# reverse primer length
RPL=21

# find OTUs using denoising and 97% clustering
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU 0 0

# there's a bug in usearch 10 - below is the workaround
cd $PROJECT_FOLDER/data/$RUN
sed -i -e 's/Zotu/OTU/' FUN.zotus.fa

# assign taxonomy to OTUs
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU 

# create OTU tables
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU $FPL $RPL

### BACTERIAL OTUS ###
SSU=BAC

# forward primer length
FPL=17 

# reverse primer length
RPL=21

# find OTUs using denoising and 97% clustering
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU $FPL $RPL

# there's a bug in usearch 10 - below is the workaround
cd $PROJECT_FOLDER/data/$RUN
sed -ie 's/Zotu/OTU/' BAC.zotus.fa

# assign taxonomy to OTUs
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU 

# create OTU tables
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU $FPL $RPL
