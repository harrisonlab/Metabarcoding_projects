library(DESeq2)
library("BiocParallel")
register(MulticoreParam(8))

# filter samples and set location to factor levels
myfiltbiom <- prune_samples(sample_data(mybiom)[[10]]!="duplicate",mybiom)
myfiltbiom <- prune_samples(sample_data(myfiltbiom)[[1]]!="C",myfiltbiom)
myfiltbiom@sam_data$location <- as.factor(myfiltbiom@sam_data$meters)
#colnames(sample_data(myfiltbiom))[c(1,6,11)] <- c("Sample","Distance","Orchard")
#levels(sample_data(myfiltbiom)[[1]]) <- c("C","Aisle","Tree")

design=~orchard + condition + orchard:location + orchard:condition
dds <- phylo_to_des(myfiltbiom,design=design,parallel=T,fit=T)

contrast=list(c("conditiongrass", "conditiontree" )) # main effect
contrast=list( "orchardCider.conditiongrass","orchardCider.conditiontree") # Grass:Tree effect (cider)
contrast=list( "orchardDessert.conditiongrass","orchardDessert.conditiontree") # Grass:Tree (Dessert)

alpha <- 0.05
res <- results(dds,contrast=contrast,alpha=alpha,parallel=T)
res.merge <- merge(as.data.frame(res),tax_table(myfiltbiom),by="row.names",all.x=TRUE)
rownames(res.merge) <- res.merge$Row.names
res.merge <- res.merge[-1]
sig.res <- subset(res.merge,padj<=alpha)
