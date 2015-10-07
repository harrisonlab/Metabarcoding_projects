#### Load libraries
library(Biostrings)

args <- commandArgs(TRUE) #get command line variables
#args <- c(".","*.\\.58","*.\\.lsu","../S91_R2.fa", "ITS2")
# 1: folder
# 2: regex start file names
# 3: regex end file names
# 4: sample fasta
print(args[1])
print(args[2])
print(args[3])
print(args[4])

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

r_start <- load_tables(list.files(args[1],args[2],full.names=T),max,T)
r_end <- load_tables(list.files(args[1],args[3],full.names=T),min)


mytable <- merge(r_start,r_end,by.all=T)
mytable <- na.omit(mytable)
mytable <- mytable[((mytable$end-mytable$start)>40),]

myfasta <- readDNAStringSet(args[4])
mytable2 <- as.data.frame(myfasta@ranges@NAMES)
colnames(mytable2) <- "seq"
mytable <- merge(mytable2,mytable,all=T)

mytable[c("start")][is.na(mytable[c("start")])] <- 0
mytable[c("end")][is.na(mytable[c("end")])] <- 1


#------# Reorder : specific to fasta labels in this experiment, regex matches will need to be changed for any other labeling system
f1 <- gsub('^[A-Za-z=]','',mytable$seq)
m <- regexpr("[0-9]+$",f1)
f2 <- as.numeric(regmatches(f1,m))
mytable <- mytable[order(f2),]
f1 <- gsub('^S','',mytable$seq)
m <- regexpr("^[0-9]+",f1)
f2 <- as.numeric(regmatches(f1,m))
mytable <- mytable[order(f2),]
#------#

ITS_IR <- IRanges(start=mytable$start+1,end=mytable$end-1,names=mytable$seq)

ITS <- DNAStringSet(myfasta,start=ITS_IR@start,width=ITS_IR@width)
ITS <- complement(reverse(ITS))
writeXStringSet(ITS,"ITS2.fa")
