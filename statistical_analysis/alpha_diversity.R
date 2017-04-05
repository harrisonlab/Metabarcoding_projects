# filter out duplicates and control samples
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
colnames(sample_data(myfiltbiom))[c(1,6,11)] <- c("Sample","Distance","Orchard")
levels(sample_data(myfiltbiom)[[1]]) <- c("C","Aisle","Tree")
myfiltbiom<-prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)

# Get alpha data
all_alpha <- plot_richness(myfiltbiom,returnData=T)

# ANOVA of alpha data
summary(aov(Chao1~location+(Sample*Orchard),all_alpha))[[1]][2]
summary(aov(Shannon~location+(Sample*Orchard)),all_alpha))[[1]][2]
summary(aov(Simpson~location+(Sample*Orchard),all_alpha))[[1]][2]
#~condition*orchard+location

sample_data(myfiltbiom)$Class <- paste(sample_data(myfiltbiom)$Orchard,sample_data(myfiltbiom)$Sample,sep=" ")

pdf("Alpha.pdf", height=8,width=8)
plot_richness(myfiltbiom,x="Class",color="Distance",measures=c("Chao1", "Shannon", "Simpson"))
dev.off()
