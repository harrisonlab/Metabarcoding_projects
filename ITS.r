#################################################
#                                               #  
# No longer implemeted                          #
# replaced with rm_SSU_58Ss.R and rm_58Se_LSU.R #
#                                               #
#################################################


#### Load libraries
library(Biostrings)

args <- commandArgs(TRUE) #get command line variables

fasta_headers <- as.data.frame(read.table("fasta_headers",header=FALSE,,stringsAsFactors=F)) ###Read fasta headers

#### Read hmm tables, combine and remove duplicates
ITS1Starts_1 <- read.table("outdata/ITS1Starts_1", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS1Starts_2 <- read.table("outdata/ITS1Starts_2", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS1Starts_1 <- ITS1Starts_1[,c(3,8)]
ITS1Starts_2 <- ITS1Starts_2[,c(3,8)]
ITS1Start <- rbind(ITS1Starts_1,ITS1Starts_2)
ITS1Start <- ITS1Start[!duplicated(ITS1Start[,1]),]
rm(ITS1Starts_1,ITS1Starts_2)

ITS1Ends_1 <- read.table("outdata/ITS1Ends_1", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS1Ends_2 <- read.table("outdata/ITS1Ends_2", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS1Ends_1 <- ITS1Ends_1[,c(3,7)]
ITS1Ends_2 <- ITS1Ends_2[,c(3,7)]
ITS1End <- rbind(ITS1Ends_1,ITS1Ends_2)
ITS1End <- ITS1End[!duplicated(ITS1End[,1]),]
rm(ITS1Ends_1,ITS1Ends_2)

ITS2Starts_1 <- read.table("outdata/ITS2Starts_1", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS2Starts_2 <- read.table("outdata/ITS2Starts_2", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS2Starts_1 <- ITS2Starts_1[,c(3,7)]
ITS2Starts_2 <- ITS2Starts_2[,c(3,7)]
ITS2Start <- rbind(ITS2Starts_1,ITS2Starts_2)
ITS2Start <- ITS2Start[!duplicated(ITS2Start[,1]),]
rm(ITS2Starts_1,ITS2Starts_2)

ITS2Ends_1 <- read.table("outdata/ITS2Ends_1", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS2Ends_2 <- read.table("outdata/ITS2Ends_2", skip=2,header=FALSE,fill=TRUE,stringsAsFactors=F)
ITS2Ends_1 <- ITS2Ends_1[,c(3,8)]
ITS2Ends_2 <- ITS2Ends_2[,c(3,8)]
ITS2End <- rbind(ITS2Ends_1,ITS2Ends_2)
ITS2End <- ITS2End[!duplicated(ITS2End[,1]),]
rm(ITS2Ends_1,ITS2Ends_2)

#### Add meaningful colnames to data tables
colnames(fasta_headers) <- "ID"
colnames(ITS1Start) <- c("ID","ITS1_Start")
colnames(ITS1End) <- c("ID","ITS1_End")
colnames(ITS2Start) <- c("ID","ITS2_Start")
colnames(ITS2End) <- c("ID","ITS2_End")

#### Concatenate hmm output tables, reorder into same order as fasta file (maybe change the merge to not reorder?)
mytable <- as.data.frame(Reduce(function(...) merge(...,all=T), list(fasta_headers,ITS1Start,ITS1End,ITS2Start,ITS2End)))
mytable[,c(3,5)] <- mytable[,c(3,5)]-1 

#------# specific to fasta labels in this experiment, regex matches will need to be changed for any other labeling system
f1 <- gsub('^[A-Za-z=]','',mytable$ID)
m <- regexpr("[0-9]+$",f1)
f2 <- as.numeric(regmatches(f1,m))
mytable <- mytable[order(f2),]
f1 <- gsub('^S','',mytable$ID)
m <- regexpr("^[0-9]+",f1)
f2 <- as.numeric(regmatches(f1,m))
mytable <- mytable[order(f2),]
#------#

mytable[is.na(mytable)] <- 0 # convert NA values to 1

##### Check ITS regions for consistency such as end not before start, no end or no start. Set bad values to zero
mytable[,2] <-  ifelse(mytable[,2]>=mytable[,3],0,mytable[,2])
mytable[,4] <-  ifelse(mytable[,4]>=mytable[,5],0,mytable[,4])
mytable[,2] <-  ifelse(mytable[,2]*mytable[,3]==0,0,mytable[,2])
mytable[,3] <-  ifelse(mytable[,2]*mytable[,3]==0,0,mytable[,3])
mytable[,4] <-  ifelse(mytable[,4]*mytable[,5]==0,0,mytable[,4])
mytable[,5] <-  ifelse(mytable[,4]*mytable[,5]==0,0,mytable[,5])
mytable[(mytable==0)]<-1

##### Load fasta, convert mytable to IRange objects and extract ITS regions from fasta file
myfasta <- readDNAStringSet(args[1])
ITS1 <- IRanges(c(mytable$ITS1_Start),c(mytable$ITS1_End), names=mytable$ID)
ITS2 <- IRanges(c(mytable$ITS2_Start),c(mytable$ITS2_End), names=mytable$ID)
ITS1 <- DNAStringSet(myfasta,start=ITS1@start+1,width=ITS1@width-1)
ITS2 <- DNAStringSet(myfasta,start=ITS2@start+1,width=ITS2@width-1)
ITS <- xscat(ITS1,ITS2)
names(ITS) <-names(myfasta)
writeXStringSet(ITS,"ITS.all.fa") # write output to file
