
##################################################################################
#################################################################################

#===============================================================================
#       Alpha diversity analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

pdf("Alpha_diversity.pdf")
lapply(myfiltbioms ,function(obj) plot_richness(obj,x="condition",color="condition",measures=c("Chao1", "Shannon", "Simpson")))
dev.off()

#===============================================================================
#       network analysis
#===============================================================================

myfiltbioms <- lapply(mybioms,function(obj) prune_samples(sample_data(obj)$condition!="C",obj))

cotables <- lapply(myfiltbioms,function(obj) as.data.frame(as.matrix(otu_table(obj))))
cotables_h <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Healthy",myfiltbioms[[i]])))]) 
cotables_s <- lapply(seq(1,length(myfiltbioms)),function(i)  cotables[[i]][,row.names(sample_data(prune_samples(sample_data(myfiltbioms[[i]])$condition=="Symptom",myfiltbioms[[i]])))])

cotables_h <- lapply(cotables_h, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])
cotables_s <- lapply(cotables_s, function(obj) obj[rowSums(obj)>5,colSums(obj)>5])

lapply(seq(1,length(cotables_h)), function(i) cotables_h[[i]][cotables_h[[i]]>0] <<- 1)
lapply(seq(1,length(cotables_s)), function(i) cotables_s[[i]][cotables_s[[i]]>0] <<- 1)

CFcoHmodels <- mclapply(cotables_h, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)
CFcoSmodels <- mclapply(cotables_s, function(obj) cooccur2(obj,type = "spp_site",spp_names = T,thresh = T),mc.cores=4)

lapply(seq(1,length(CFcoHmodels)), function(i) {
	CFcoHmodels[[i]]$results$padj <<- p.adjust(apply(CFcoHmodels[[i]]$results[,8:9],1, min),"BH")
	CFcoSmodels[[i]]$results$padj <<- p.adjust(apply(CFcoSmodels[[i]]$results[,8:9],1, min),"BH")
})

lapply(CFcoHmodels, function(obj) {nrow(obj$results[obj$results$padj<=0.1,]})

lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoHmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))

lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_gt<=0.05,]))
lapply(CFcoSmodels, function(obj) nrow(obj$results[obj$results$padj<=0.05&obj$results$p_lt<=0.05,]))


X <- rbind.fill(lapply(CFcoHmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))
Y <- rbind.fill(lapply(CFcoSmodels, function(obj) head(obj$results[order(obj$results$padj),],6)))


HcoHmodel16$results$p_lt

write.table(X,"cooc.healthy.txt",sep="\t",quote=F,row.names=F)

head(CHcoHmodel$results[order(CHcoHmodel$results$p_lt,CHcoHmodel$results$padj),])
X <- head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_lt,CHcoHmodel16$results$padj),])
X <- rbind(X,head(CHcoHmodel16$results[order(CHcoHmodel16$results$p_gt,CHcoHmodel16$results$padj),]))
