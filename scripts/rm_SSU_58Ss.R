#### Load libraries
library(Biostrings)

args <- commandArgs(TRUE) #get command line variables
#args <- c(".","*.\\.ssu","*.\\.58ss","../S77_R1.fa", "ITS1")
# 1: folder
# 2: regex start file names
# 3: regex end file names
# 4: sample fasta
# 5: sample ID
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

load_tables <- function (files,func,End=F) {
  tables <- lapply(files, load_table,func=func,End=End)
  do.call(rbind, tables)
}

load_table <- function(file,func,End) {
 tryCatch({
    X <- read.table(file,skip=2,header=F)
    X <- X[with(X, ave(V14,V3, FUN=max)==V14),]
    if (isTRUE(End)) {
      X <- X[!duplicated(X[,c(3,8)]),c(3,8)]
      X <- X[with(X, ave(V8,V3,FUN=func)==V8),]
      colnames(X) <- c("seq","end")
    } else {
      X <- X[!duplicated(X[,c(3,7)]),c(3,7)]
      X <- X[with(X, ave(V7,V3,FUN=func)==V7),]
      colnames(X) <- c("seq","start")
    }
      return(X) 
    }, error = function(err) {
  } )
} 

print("loading ssu tables")
r_start <- load_tables(list.files(args[1],args[2],full.names=T),max,T)
print("loading 58S tables")
r_end <- load_tables(list.files(args[1],args[3],full.names=T),min)

#SSU_test <- load_tables(head(list.files(".","*.\\.ssu"),2),max,T)
#ss58_test <- load_tables(head(list.files(".","*.\\.58"),2),min)
print("merging tables")
mytable <- merge(r_start,r_end,by.all=T)
print("removing NAs")
mytable <- na.omit(mytable)
mytable <- mytable[((mytable$start-mytable$end)>40),]

print("reading DNA table")
myfasta <- readDNAStringSet(args[4])
#myfasta <- readDNAStringSet("../S91_R1.fa")
#mytable2 <- as.data.frame(myfasta@ranges@NAMES)
#colnames(mytable2) <- "seq"
#print("merging ITS and DNA tables")
#mytable <- merge(mytable2,mytable,all=T)

myfasta<- myfasta[mytable$seq]

#print("convert NAs to 1 or 0")
#mytable[c("start")][is.na(mytable[c("start")])] <- 1
#mytable[c("end")][is.na(mytable[c("end")])] <- 0

#print("Reorder table")
#------# Reorder : specific to fasta labels in this experiment, regex matches will need to be changed for any other labeling system
#f1 <- gsub('^[A-Za-z=]','',mytable$seq)
#m <- regexpr("[0-9]+$",f1)
#f2 <- as.numeric(regmatches(f1,m))
#mytable <- mytable[order(f2),]
#f1 <- gsub('^S','',mytable$seq)
#m <- regexpr("^[0-9]+",f1)
#f2 <- as.numeric(regmatches(f1,m))
#mytable <- mytable[order(f2),]
#------#

print("set ITS_IR")
ITS_IR <- IRanges(start=mytable$end+1,end=mytable$start-1,names=mytable$seq)

print("set ITS")
ITS <- DNAStringSet(myfasta,start=ITS_IR@start,width=ITS_IR@width)

ITS <- ITS[ITS@ranges@width>=140]
#ITS <- ITS[ITS@ranges@width>0]
#ITS <- stackStrings(ITS,0,max(ITS@ranges@width),Lpadding.letter="N",Rpadding.letter="N")
#ITS <- subseq(ITS ,start=2,width = (max(ITS@ranges@width)-1))

print("write ITS")
writeXStringSet(ITS,paste(args[5],".r1.fa",sep=""))

#ITS <- sapply(ITS_IR@NAMES,function(x) DNAStringSet(myfasta[(myfasta@ranges@NAMES==x)],start=ITS_IR@start[(ITS_IR@NAMES==x)],width=ITS_IR@width[(ITS_IR@NAMES==x)],use.names=T))

#lapply(ITS,function(x) writeXStringSet(x,paste(args[5],".fa",sep=""),append=T))