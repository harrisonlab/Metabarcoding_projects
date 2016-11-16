
myfiltbiom1 <- prune_samples(sample_data(myfiltbiom)$block!=3,myfiltbiom)
myfiltbiom2 <- prune_samples(sample_data(myfiltbiom)$block!=1,myfiltbiom)

cond<-"Y"
pc.x <- scores(mypca)[rownames(scores(mypca))%in%rownames(sample_data(myfiltbiom1)[sample_data(myfiltbiom1)$condition==cond]),]
col.x <- col.x <- sample_data(myfiltbiom1)[sample_data(myfiltbiom1)$condition==cond,]
pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
tree12<-pc.reshape


pc.x <- scores(mypca)[rownames(scores(mypca))%in%rownames(sample_data(myfiltbiom2)[sample_data(myfiltbiom2)$condition==cond]),]
col.x <- col.x <- sample_data(myfiltbiom2)[sample_data(myfiltbiom2)$condition==cond,]
pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
tree23<-pc.reshape

cond <- "N"

pc.x <- scores(mypca)[rownames(scores(mypca))%in%rownames(sample_data(myfiltbiom1)[sample_data(myfiltbiom1)$condition==cond]),]
col.x <- col.x <- sample_data(myfiltbiom1)[sample_data(myfiltbiom1)$condition==cond,]
pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
aisle12<-pc.reshape


pc.x <- scores(mypca)[rownames(scores(mypca))%in%rownames(sample_data(myfiltbiom2)[sample_data(myfiltbiom2)$condition==cond]),]
col.x <- col.x <- sample_data(myfiltbiom2)[sample_data(myfiltbiom2)$condition==cond,]
pc.dt<- data.table(merge(pc.x,col.x,by="row.names"))
pc.reshape <- dcast(pc.dt,distance~.,fun.aggregate=function(x) mean(x,na.rm=T),value.var=c(names(pc.dt)[grep("PC",names(pc.dt))]))
names(pc.reshape)[grep("PC",names(pc.reshape))] <- sub("_.*","",names(pc.reshape)[grep("PC",names(pc.reshape))])
aisle23<-pc.reshape


cutoff <- 20
pc<-"PC1"

test1 <- correr2_dat(tree12[[pc]])
test2 <- correr2_dat(tree23[[pc]])
test3 <- sapply(seq(1,length(test1)),function(x) rbind(test1[[x]],test2[[x]]))

ct1 <- sapply(seq(1,length(test3)),function(x) cor(test3[[x]][,1],test3[[x]][,2]))

test1 <- correr2_dat(aisle12[[pc]])
test2 <- correr2_dat(aisle23[[pc]])
test3 <- sapply(seq(1,length(test1)),function(x) rbind(test1[[x]],test2[[x]]))

ca1 <- sapply(seq(1,length(test3)),function(x) cor(test3[[x]][,1],test3[[x]][,2]))

d<-as.data.frame(cbind(ct1,ca1,pc.reshape$distance[1:(length(pc.reshape$distance)-2)]))
d<- d[1:(length(d$V3[d$V3<=cutoff])),]
names(d) <- c(paste0("Tree_",pc),paste0("Aisle_",pc),"Distance")
d2 <- melt(d,id="Distance")
colnames(d2)[3] <- c("Correlation")

g <- ggplot(d2)
g <- g + coord_fixed(ratio = 10, xlim = NULL, ylim = NULL, expand = TRUE)
g <- g + theme_bw()
g <- g + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g <- g + theme(axis.line.x = element_line(size=0.5,colour = "black"),axis.line.y = element_line(size=0.5,colour = "black"),axis.text = element_text(colour = "black"))
g <- g + geom_line(na.rm=T,aes(x=Distance, y=Correlation, colour=variable))
g <- g + scale_colour_manual(values=c("red","green"))
g



