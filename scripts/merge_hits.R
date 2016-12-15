#### Load libraries
library(data.table)
args <- commandArgs(TRUE)
r1 <- fread(args[1])
r2 <- fread(args[2])
r2$V1 <- sub("R2","R1",r2$V1)
m <- merge(r1,r2,by=c("V1","V2"))
rm(r1,r2)
m$V1 <- sub("_.*","",m$V1)
m <- m[(V3.x>=95|V3.y>=95)]
m <- m[,1:2,with=F]
m[,count:=.N,by=.(V1,V2)]
m <- unique(m)
m2<-m
m <- dcast(m, V2 ~ V1, value.var = "count")
write.table(m,args[3],sep="\t",row.names=F,quote=F,na="0") 


#biom_data

m2<-m2[order(V1)]
m2<-m2[,V1:=as.factor(V1)]
m2<-m2[,sam_pos:=(as.numeric(V1)-1)]
m2<-m2[order(V2)]
m2<-m2[,V2:=as.factor(V2)]
m2<-m2[,otu_pos:=(as.numeric(V2)-1)]
write.table(m2[,5:3,with=F],"data_biom",sep=",",row.names=F,col.names=F,quote=F,na="0")
write.table(levels(m2$V2),"row_biom",sep=",",row.names=F,col.names=F,quote=F)
write.table(levels(m2$V1),"col_biom",sep=",",row.names=F,col.names=F,quote=F)



# m2<-m2[order(otu_pos,sam_pos)]
#m2 <- m2[,5:3,with=F]