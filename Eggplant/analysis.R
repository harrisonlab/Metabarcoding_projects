#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(data.table)
library(tidyverse)
library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(metacoder)
library(ape)
library(metafuncs)
library(BiocParallel)

#environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)
#assignInNamespace("plot_ordination",value=plot_ordination,ns="phyloseq")

#===============================================================================
#       Functions
#===============================================================================


em_splitter <- function(test,col="contrast",f=" - ") {
  DF <- do.call(rbind,
                strsplit(gsub(f,"\t",gsub(",","\t",test[[col]],perl=T)),"\t")
  )
  DT <- as.data.table(cbind(DF,test[,-1]))
  setnames(DT,names(DT)[1:ncol(DF)],paste0("V",names(DT)[1:ncol(DF)]))
}


# adapted version of phyloTaxaTidy to mash names by conf
phyloTidy <- function(taxData,conf,setCount=T) {
  X <- taxData[,1:7]
  Y <- taxData[,9:15]
  Y[,1] <-1.00
  X[apply(Y,2,function(y)y<conf)] <- "unknown"
  n <- c("kingdom","phylum","class","order","family","genus","species")
  X <-apply(X,2,function(x) {
    x <- sub("\\(.*","",x)
    x <- sub("_SH[0-9]+\\..*","",x)
    x <- gsub("_"," ",x)
    #x <- sub("_"," ",x)
    x
  })
  
  td <- as.data.table(X,keep.rownames = "OTU")
  cols <- names(td)[-1]
  num <- td[,.N,by=cols]
  num2 <- invisible(apply(num,1,function(x) as.data.table(t(matrix(x[-8],nrow=7,ncol=as.numeric(x[8]))))))
  invisible(lapply(num2,setnames,cols))
  invisible(lapply(num2,function(x)x[,counts:=1:nrow(x)]))
  td2 <- do.call("rbind",num2)
  td[,(cols):=lapply(.SD,as.character),.SDcols=cols]
  setorderv(td, cols[-1])
  setorderv(td2,cols[-1])
  
  td$counts <- td2$counts
  td[,(cols):=lapply(.SD,as.factor),.SDcols=cols]
  taxData <- as.data.frame(td)
  rownames(taxData) <- taxData$OTU
  taxData <- taxData[,-1]
  if(setCount) {
    taxData$species <- as.factor(paste(taxData$species,taxData$counts,sep=" "))
    taxData <- taxData[,-8]
  }
  taxData
}

funq <- function(lfc,p){lfc[c(which(p>0.1),which(is.na(p)))] <- 0;lfc}

qf <- function(x)sub("Fungi;unknown.*","REMOVE",x)

# get the full taxonomy for each taxon_id
# get_full_taxa <- function(taxon_id,obj) {
#   t_var   <- c(obj$taxon_indexes()[[taxon_id]],obj$supertaxa()[[taxon_id]])
#   t_names <- sapply(t_var,function(x) {obj$taxa[[x]]$get_name()})
#   paste(rev(t_names),collapse = ";")
# }

get_full_taxa <- function(taxon_id,t1=t1,t2=t2) {
  t_names <- t2[c(taxon_id,unlist(t1[taxon_id]))]
  paste(rev(t_names),collapse = ";")
}

MCDESeq <- function (OBJ,formula,Rank,colData) {
  
  o <- OBJ$clone(deep=T)
  formula <- formula(formula)
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == Rank,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData$sample_id)
  
  # metacoder differential analysis - this is bobbins  
  #FUN_ENDO$data$diff_table <- compare_groups(FUN_ENDO, dataset = "otu_table",
  #                                           cols = colData$sample_id, # What columns of sample data to use
  #                                           groups = colData$status) # What category each sample is assigned to
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  
  
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,c(-1)]
  
  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,c(-1)]
  
  # set character columns to factors
  colData <- as.data.frame(unclass(as.data.frame(colData)))
  
  # create new dds object from combined tables
  dds <- DESeqDataSetFromMatrix(rbind(otu_table,tax_abund),colData,~1)
  #cols <- colnames(colData(dds))[sapply(colData(dds),is.factor)]
  #dds$status <- as.factor(dds$status)
  design(dds) <- formula#~status
  dds <- DESeq(dds)
  return(dds)
  # res <- results(dds,alpha=alpha,contrast=contrast)
  # 
  # # make a data table from the resuults
  # res_merge <- as.data.table(as.data.frame(res),keep.rownames="taxon_id")
  # 
  # # add second log fold change column with fold changes for non sig taxa set to 0 (for colouring purposes on graphs)
  # res_merge[,log2FoldChange2:=funq(log2FoldChange,padj)]
  # 
  # # add full taxonomy to results (this is slow - good reason not to use S4 object model!)
  # 
  # # order the results
  # setkey(res_merge,taxon_id)
  # 
  # # these are required due to the restrictive interace of taxa/metacoder 
  # # reduces the time required to get the full taxonomy of each taxon id by about x10,000
  # t1 <- o$supertaxa()
  # t2 <- o$taxon_names()
  # res_merge[,taxonomy:=sapply(1:nrow(res_merge),get_full_taxa,t1,t2)]
  # 
  # # add results as diff_table (could be called anything - but for consistency with metacoder keep same names)
  # #
  # as_tibble(res_merge)
}

MCDESres <- function (OBJ,dds,Rank,colData,contrast,alpha=0.1) {
  
  o <- OBJ$clone(deep=T)
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == Rank,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData$sample_id)
  
  # metacoder differential analysis - this is bobbins  
  #FUN_ENDO$data$diff_table <- compare_groups(FUN_ENDO, dataset = "otu_table",
  #                                           cols = colData$sample_id, # What columns of sample data to use
  #                                           groups = colData$status) # What category each sample is assigned to
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,c(-1,-2)]
  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,c(-1)]
  
  # set character columns to factors
  colData <- as.data.frame(unclass(as.data.frame(colData)))
  
  # create new dds object from combined tables
  #dds <- DESeqDataSetFromMatrix(rbind(otu_table,tax_abund),colData,~1)
  #cols <- colnames(colData(dds))[sapply(colData(dds),is.factor)]
  #dds$status <- as.factor(dds$status)
  #design(dds) <- formula#~status
  #dds <- DESeq(dds)
  res <- results(dds,alpha=alpha,contrast=contrast)
  
  # make a data table from the resuults
  res_merge <- as.data.table(as.data.frame(res),keep.rownames="taxon_id")
  
  # add second log fold change column with fold changes for non sig taxa set to 0 (for colouring purposes on graphs)
  res_merge[,log2FoldChange2:=funq(log2FoldChange,padj)]
  
  # add full taxonomy to results (this is slow - good reason not to use S4 object model!)
  
  # order the results
  setkey(res_merge,taxon_id)
  
  # these are required due to the restrictive interace of taxa/metacoder 
  # reduces the time required to get the full taxonomy of each taxon id by about x10,000
  t1 <- o$supertaxa()
  t2 <- o$taxon_names()
  res_merge[,taxonomy:=sapply(1:nrow(res_merge),get_full_taxa,t1,t2)]
  
  # add results as diff_table (could be called anything - but for consistency with metacoder keep same names)
  #
  as_tibble(res_merge)
}

#===============================================================================
#       Load data
#===============================================================================

ubiome_BAC <- loadData("BAC.otu_table.txt","colData","BAC.utax.taxa",RHB="BAC")
names(ubiome_BAC$countData) <- gsub("_.*","",names(ubiome_BAC$countData))
rownames(ubiome_BAC$colData) <- names(ubiome_BAC$countData)

ubiome_FUN <- loadData("FUN.otu_table.txt","colData","FUN.utax.taxa",RHB="FUN")
names(ubiome_FUN$countData) <- gsub("_.*","",names(ubiome_FUN$countData))
rownames(ubiome_FUN$colData) <- names(ubiome_FUN$countData)

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
# invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# combinedTaxa <- combineTaxa("FUN.utax.taxa")
# countData <- combCounts(combinedTaxa,countData)
# taxData <- combTaxa(combinedTaxa,taxData)
# ubiome_FUN$countData <- countData
# ubiome_FUN$taxData <- taxData

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiome_FUN$dds <- ubiom_to_des(ubiome_FUN,filter=expression(colSums(countData)>=500))
ubiome_BAC$dds <- ubiom_to_des(ubiome_BAC,filter=expression(colSums(countData)>=1000))


#===============================================================================
#      ****FUNGI/BACTERIA****
#===============================================================================
# Fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Sample rarefaction plots
#===============================================================================        

library(grid)
library(gridExtra)
library(viridis)

gfunc <- function(countData,coldata,title) {        
  colData <- colData[names(countData),]
  
  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))
  
  # get cummulative sum of each sample
  DT <- cumsum(DT)    
  
  # log the count values                            
  DT <- log10(DT)
  
  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {x[(which.max(x)+1):length(x)]<- NA;x}))
  
  # remove rows with all NA
  DT <- DT[rowSums(is.na(DT)) != ncol(DT), ]
  
  # add a count column to the data table
  DT$x <- seq(1,nrow(DT))
  
  # melt the data table for easy plotting 
  MDT <- melt(DT,id.vars="x")
  
  # create an empty ggplot object from the data table
  g <- ggplot(data=MDT,aes(x=x,y=value,colour=variable))
  
  # remove plot background and etc.
  g <- g + theme_classic_thin() %+replace% theme(legend.position="none",axis.title=element_blank())
  
  # plot cumulative reads
  g <- g + geom_line(size=1.5) + scale_colour_viridis(discrete=T)
  
  # add axis lables
  g <- g + ggtitle(title)
  #g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")
  
  # print the plot
  g
}

invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))
g1 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Bacteria")
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
g2 <-gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Fungi")

glegend <- get_legend(g) # this won't work
ggsave("rarefaction_all.pdf",grid.arrange(g1,g2,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=2))                              


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Recreate dds object and don't filter for low counts before running Alpha diversity

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
plot_alpha(counts(dds,normalize=T),colData(dds),design="Status",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed"),type="box")		   

ggsave(paste(RHB,"Alpha.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed")))
ggsave(paste(RHB,"Alpha_Chao1.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Chao1"))) # ,limits=c(0,xxx,"Chao1")
ggsave(paste(RHB,"Alpha_Shannon.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Shannon")))
ggsave(paste(RHB,"Alpha_Simpson.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Simpson")))
ggsave(paste(RHB,"Alpha_Observed.pdf",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Observed")))

### permutation based anova on diversity index ranks ###
# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",returnData=T)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,colData,by=c("Samples"="Sample_FB"))) # or sample_on

# perform anova for each index
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
setkey(all_alpha_ord,S.chao1)
print("Chao1")
summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
setkey(all_alpha_ord,shannon)
print("Shannon")
summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
setkey(all_alpha_ord,simpson)
print("simpson")
summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Block + Treatment + Genotype + Treatment * Genotype,all_alpha_ord))
sink()




#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
# ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

# Filter chloroplast/mitochondria (bacteria only)
#dds <- dds[c(-1,-12),]

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# add a new column (not strictly necessary)
dds$field_pair <- as.factor(paste(dds$Field,dds$Pair,sep="_"))

# plot the PCA
plotOrd(d,colData(dds),shape="Institution",design=c("Type","Status"),alpha=0.75,cbPalette=T,axes=c(1,2))
#plotOrd(d,colData,design="parent",shape="sex",axes=c(2,3),alpha=0.75,cbPalette=T)

#colData(dds)$Pair <- as.factor(colData(dds)$Pair)

lapply(1:4,function(i){
  summary(aov(mypca$x[,i]~Institution+Field+ field_pair+Type*Status,data=colData(dds)))
})

# fm1 <- lm(mypca$x[,1]~Institution+Field+ field_pair+Type*Status,data=colData(dds))
# (res <- emmeans(fm1,~Type*Status))
# test <- em_splitter(summary(pairs(res)))
# 
# # within Type
# test[V1==V3&p.value <=0.1,-3]
# 
# # between
# test[V2==V4&p.value <=0.1,-4]
# 
# # plot
# g <- ggplot(summary(res),aes(x=Type:Status,y=emmean)) 
# g + geom_bar(stat="identity") + 
#   xlab("Type Status") +
#   ylab("PC score") +
#   scale_color_manual(values=cbPalette) +
#   geom_errorbar(aes(ymin=emmean -SE, ymax=emmean +SE), width=.05,position=position_dodge(.25)) +
#   theme_facet_blank(angle=0,hjust=0.5) %+replace% theme(legend.position = "bottom",legend.justification = "left")
# 
# dds_stem <-dds[,dds$Type=="stem"]
# dds_root <-dds[,dds$Type=="root"]
# dds_soil <-dds[,dds$Type=="soil"]
# 
# mypca <- des_to_pca(dds_stem)
# d <-t(data.frame(t(mypca$x)*mypca$percentVar))
# plotOrd(d,colData(dds_stem),shape="Institution",design="Status",alpha=0.75,cbPalette=T,axes=c(1,2),ylims = c(-25,25))
# 
# 

sum_squares <- apply(mypca$x,2,function(x) 
  summary(aov(x~Institution+Field+ field_pair+Type*Status,data=colData(dds)))[[1]][2]
)
sum_squares <- do.call(cbind,sum_squares)
x<-t(apply(sum_squares,2,prop.table))
perVar <- x * mypca$percentVar
sink(paste(RHB,"PCA_sum_squares.txt",sep="_"))
  colSums(perVar)
  colSums(perVar)/sum(colSums(perVar))*100
sink()
### ADONIS ###
vg <- vegdist(t(counts(dds,normalize=T)),method="bray")
sink(paste(RHB,"ADONIS.txt",sep="_"))
 set.seed(sum(utf8ToInt("Xiangming Xu")))
 fm1 <- adonis(vg~Institution+Field+ field_pair+Type*Status,colData(dds),permutations = 1000)
sink()

# nmds ordination
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),as.data.frame(colData(dds)),taxData))
ord_rda <- phyloseq::ordinate(myphylo,method="NMDS",distance="bray",formula= ~Field+field_pair + Type*Status)		

otus <- scores(ord_rda,"species")
nmds <- scores(ord_rda)

g <- plotOrd(nmds,colData(dds),design=c("Type","Status"),shape="Institution",alpha=0.75,cbPalette=T)
g#+geom_point(data=as.data.frame(otus),inherit.aes = F,aes(x=NMDS1,y=NMDS2))

taxmerge <-data.table(inner_join(data.table(OTU=rownames(otus),as.data.frame(otus)),data.table(OTU=rownames(taxData),taxData)))
taxmerge$phy <- taxaConfVec(taxData[,-8],conf=0.9,level=which(colnames(taxData)=="phylum"))
taxmerge$cls <- taxaConfVec(taxData[,-8],conf=0.9,level=which(colnames(taxData)=="class"))

phylum <- taxmerge[,lapply(.SD,mean),by=phy,.SDcols=c("NMDS1","NMDS2")]
cls <- taxmerge[,lapply(.SD,mean),by=cls,.SDcols=c("NMDS1","NMDS2")]

g + geom_segment(inherit.aes = F,data=cls,aes(xend=NMDS1,yend=NMDS2,x=0,y=0),size=1.5,arrow=arrow()) +
  geom_text(inherit.aes = F,data=cls,aes(x=NMDS1,y=(NMDS2+sign(NMDS2)*0.05),label=cls))  



#===============================================================================
#       differential analysis
#===============================================================================
dds$FTP <- as.factor(paste0(dds$Field,dds$Type,dds$Pair))
# p value for FDR cutoff
alpha <- 0.1

# the model
design <- ~Pair+Institution + Type*Status

# add design to dds object
design(dds) <- design

# run model
dds <- DESeq(dds,parallel=F)

# difference between sexes
res <- results(dds,alpha=alpha,parallel=F,contrast = c("Status","Diseased","Healthy"))
res2 <- results(dds,alpha=alpha,parallel=F,contrast = c("Type","stem","soil"))
res3 <- results(dds,alpha=alpha,parallel=F,contrast = c("Institution","WorldVeg","MARI"))



res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,"Fungi_res.txt",sep="\t",na="",quote=F)
