#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(data.table)
library(tidyverse)
# library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(metacoder)
library(ape)
library(metafuncs)
library(BiocParallel)
register(MulticoreParam(2))

library(grid)
library(gridExtra)
library(viridis)
library(cowplot)
library(VennDiagram)
#environment(plot_ordination) <- environment(ordinate) <- environment(plot_richness) <- environment(phyloseq::ordinate)
#assignInNamespace("plot_ordination",value=plot_ordination,ns="phyloseq")

#===============================================================================
#       Functions
#===============================================================================

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

em_splitter <- function(test,col="contrast",f=" - ") {
  DF <- do.call(rbind,
                strsplit(gsub(f,"\t",gsub(",","\t",test[[col]],perl=T)),"\t")
  )
  DT <- as.data.table(cbind(DF,test[,-1]))
  setnames(DT,names(DT)[1:ncol(DF)],paste0("V",names(DT)[1:ncol(DF)]))
}

triple.venn <- function(A,B,C,...) {
  ab <- sum(duplicated(c(A,B)))
  ac <- sum(duplicated(c(A,C)))
  bc <- sum(duplicated(c(B,C)))
  abc <- sum(duplicated(c(c(A,B)[duplicated(c(A,B))],C)))
  draw.triple.venn(length(A),length(B),length(C),ab,bc,ac,abc,...)
}

venn.out <- function(A,B,C) {
  ab <- c(A,B)[duplicated(c(A,B))]
  ac <- c(A,C)[duplicated(c(A,C))]
  bc <- c(B,C)[duplicated(c(B,C))]
  abc <- c(ab,C)[duplicated(c(ab,C))]
  list(AuB=ab,AuC=ac,BuC=bc,AuBuC=abc)
}



alpha_limit <- function(g,limits=c(0,4000),p=2,l=17) {
    g2 <- g + coord_cartesian(y=limits)
    g <- ggplotGrob(g)
    g2 <- ggplotGrob(g2)
    g[["grobs"]][[p]] <- g2[["grobs"]][[p]]
    g[["grobs"]][[l]] <- g2[["grobs"]][[l]]
    g
  }

combineByTaxa <- 
  function (taxData,countData, rank = "species", confidence = 0.95, column_order = -8) {
    require(plyr)
    require(data.table)
    
    # reorder taxonomy table (remove rank column by default)
    taxData <- taxData[, column_order]
    
    # get rank level
    level <- which(c("kingdom", "phylum", "class", "order", 
                     "family", "genus", "species")==rank)
    
    # select taxonomy at given confidence
    taxData <- phyloTaxaTidy(taxData, confidence,level=level)
    
    # convert to data table
    TD <- data.table(taxData, keep.rownames = "OTU")
    
    # get list of OTUs in each entry at the given rank
    combinedOTUs <- ddply(TD, ~rank, summarize, OTUS = list(as.character(OTU)))
    
    # combine the OTUs into list format
    combinedOTUs <- combinedOTUs[lapply(TD[, 2], function(x) length(unlist(x))) > 
                                   1, ]
    
    list(countData = combCounts(combinedOTUs,ubiome_FUN$countData),
         taxData   = combTaxa(combinedOTUs,taxData))
   
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
  design(dds) <- formula
  dds <- DESeq(dds)
  dds@NAMES <- Rank
  return(dds)
}

MCDESres <- function (OBJ,contrast,alpha=0.1) {
  
  o <- OBJ$clone(deep=T)
  

  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == o$data$dds@NAMES,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData(o$data$dds)$sample_id)
  
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
  #colData <- as.data.frame(unclass(as.data.frame(colData))) # why this???
  
  # run results
  res <- results(o$data$dds,alpha=alpha,contrast=contrast)
  
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
names(ubiome_BAC$countData) <- gsub("([LZ])(M[0-9]+_)([DH])(.*)","\\1\\3-\\4",names(ubiome_BAC$countData))
ubiome_BAC$colData <-ubiome_BAC$colData[names(ubiome_BAC$countData),] 
ubiome_BAC$colData$Location <- as.factor(ubiome_BAC$colData$Location)

ubiome_FUN <- loadData("FUN.otu_table.txt","colData","FUN.utax.taxa",RHB="FUN")
names(ubiome_FUN$countData) <- gsub("([LZ])(M[0-9]+_)([DH])(.*)","\\1\\3-\\4",names(ubiome_FUN$countData))
ubiome_FUN$colData <-ubiome_FUN$colData[names(ubiome_FUN$countData),] 
ubiome_FUN$colData$Location <- as.factor(ubiome_FUN$colData$Location)

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
#invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
#combinedTaxa <- combineTaxa("FUN.utax.taxa")
#countData <- combCounts(combinedTaxa,countData)
#taxData <- combTaxa(combinedTaxa,taxData)
#ubiome_FUN$countData <- countData
#ubiome_FUN$taxData <- taxData

#===============================================================================
#       Create DEseq objects
#===============================================================================

ubiome_FUN$dds <- ubiom_to_des(ubiome_FUN)
#sizeFactors(ubiome_FUN$dds) <- sizeFactors(ubiome_FUN$dds)/(ubiome_FUN$colData$ITS_qPCR/min(ubiome_FUN$colData$ITS_qPCR))
sizeFactors(ubiome_FUN$dds) <- median(ubiome_FUN$colData$ITS_qPCR)/ubiome_FUN$colData$ITS_qPCR

ubiome_BAC$dds <- ubiom_to_des(ubiome_BAC)
#sizeFactors(ubiome_BAC$dds) <- sizeFactors(ubiome_BAC$dds)/(ubiome_BAC$colData$X16S_qPCR/min(ubiome_BAC$colData$X16S_qPCR))
sizeFactors(ubiome_BAC$dds) <- median(ubiome_BAC$colData$X16S_qPCR)/ubiome_BAC$colData$X16S_qPCR

#===============================================================================
#      ****FUNGI/BACTERIA****
#===============================================================================
# Fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

dds$Location <- as.factor(dds$Location)

#===============================================================================
#       Sample rarefaction plots
#===============================================================================        

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

#invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))
gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Bacteria")
#invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
#g2 <-gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Fungi")

glegend <- get_legend(g) # this won't work
ggsave("rarefaction_all.jpg",grid.arrange(g1,g2,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=2),dpi=600)                              


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

### permutation based anova on diversity index ranks ###
# get the diversity index data
all_alpha_ord <- plot_alpha(counts(dds,normalize=F),colData(dds),design="Treatment",returnData=T)

# join diversity indices and metadata
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,colData%>%mutate(Samples=rownames(colData)),by="Samples")) # or sample_on

# perform anova for each index
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
  setkey(all_alpha_ord,S.chao1)
  print("Chao1")
  summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~Cultivar:Location+Cultivar*Status,all_alpha_ord))
  setkey(all_alpha_ord,shannon)
  print("Shannon")
  summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~Cultivar:Location+Cultivar*Status,all_alpha_ord))
  setkey(all_alpha_ord,simpson)
  print("simpson")
  summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~Cultivar:Location+Cultivar*Status,all_alpha_ord))
sink()


# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
#plot_alpha(counts(dds,normalize=T),colData(dds),design="Cultivar:Status",colour="Cultivar",measures=c("Chao1", "Shannon", "Simpson","Observed"),type="box",cbPalette = T,limits=c(0,2000,"Chao1"))

g <- plot_alpha(counts(dds,normalize=T),colData(dds),design="Cultivar:Status",colour="Cultivar",measures=c("Shannon", "Simpson","Observed"),type="box",cbPalette = T)

ggsave(paste(RHB,"Alpha.jpg",sep="_"),g)#plot_grid(alpha_limit(g,limits = c(0,7500))))
#ggsave(paste(RHB,"Alpha_Chao1.jpg",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Chao1"))) # ,limits=c(0,xxx,"Chao1")
#ggsave(paste(RHB,"Alpha_Shannon.jpg",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Shannon")))
#ggsave(paste(RHB,"Alpha_Simpson.jpg",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Simpson")))
#ggsave(paste(RHB,"Alpha_Observed.jpg",sep="_"),plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour="Genotype",measures=c("Observed")))



#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
# ggsave(paste(RHB,"OTU_counts.jpg",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

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

# anova
sink(paste(RHB,"PCA_ANOVA.txt",sep="_"))
 # percentVar
 mypca$percentVar[1:4]
 lapply(1:4,function(i){
   summary(aov(mypca$x[,i]~Cultivar:Location+Cultivar*Status,data=colData(dds)))
 })
sink()

#summary(aov(mypca$x[,1]~Institution + Field + field_pair + Type*Status,data=colData(dds)))
#summary(aov(mypca$x[,1]~Error(Institution/Field/field_pair) + Type*Status,data=colData(dds)))
#summary(aov(mypca$x[,1]~Field+ field_pair+Type*Status,data=colData(dds)))


# plot the PCA
ggsave(paste0(RHB,"_PCA.jpg"),plotOrd(d,colData(dds),shape="Cultivar",design="Status",alpha=0.75,cbPalette=T,axes=c(1,2),ylims = c(-5,+10)))
ggsave(paste0(RHB,"_PCA_2.jpg"),plotOrd(d,colData(dds),shape="Cultivar",design="Status",alpha=0.75,cbPalette=T,axes=c(2,3),ylims = c(-5,10)))

plotOrd(d,colData(dds),shape="Cultivar",design="Status",alpha=0.75,cbPalette=T,axes=c(2,3))

sum_squares <- apply(mypca$x,2,function(x) 
  summary(aov(x~Cultivar:Location+Cultivar*Status,data=colData(dds)))[[1]][2]
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
 (fm1 <- adonis(vg~ Cultivar:Location + Cultivar*Status,colData(dds),permutations = 1000))
 (fm2 <- mrpp(dat=vg,grouping=colData(dds)$Cultivar:colData(dds)$Status,permutations = 1000))
sink()

# nmds ordination
myphylo <- ubiom_to_phylo(list(counts(dds,normalize=T),as.data.frame(colData(dds)),taxData))
set.seed(sum(utf8ToInt("Xiangming Xu")))
ord_rda <- phyloseq::ordinate(myphylo,method="NMDS",distance="bray",formula= ~Cultivar:Location+Cultivar*Status)		

otus <- scores(ord_rda,"species")
nmds <- scores(ord_rda)

g <- plotOrd(nmds,colData(dds),shape="Status",design="Cultivar",alpha=0.75,cbPalette=T)
ggsave(paste0(RHB,"_NMDS.jpg"),g)#+geom_point(data=as.data.frame(otus),inherit.aes = F,aes(x=NMDS1,y=NMDS2))

taxmerge <-data.table(inner_join(data.table(OTU=rownames(otus),as.data.frame(otus)),data.table(OTU=rownames(taxData),taxData)))
taxmerge$phy <- taxaConfVec(taxmerge[,c(-1,-2,-3,-11)],conf=0.9,level=which(colnames(taxmerge)=="phylum"))
taxmerge$cls <- taxaConfVec(taxmerge[,c(-1,-2,-3,-11)],conf=0.9,level=which(colnames(taxData)=="class"))

phylum <- taxmerge[,lapply(.SD,mean),by=phylum,.SDcols=c("NMDS1","NMDS2")]
cls <- taxmerge[,lapply(.SD,mean),by=cls,.SDcols=c("NMDS1","NMDS2")]

g + geom_segment(inherit.aes = F,data=phylum,aes(xend=NMDS1,yend=NMDS2,x=0,y=0),size=1.5,arrow=arrow()) +
  geom_text(inherit.aes = F,data=phylum,aes(x=NMDS1,y=(NMDS2+sign(NMDS2)*0.05),label=phylum))  



#===============================================================================
#       FunGuild
#===============================================================================

# make a funguild compliant table
td <-  as.data.table(taxData,keep.rownames = "OTU ID")
td[,sample1:=100]
td[,taxonomy:=paste(kingdom,phylum,class,order,family,genus,species,sep=";")]
fwrite(td[,c("OTU ID","sample1","taxonomy")],"funtax.txt",sep="\t",na="",quote=F)

# run on funguild server; http://www.stbates.org/guilds/app.php

funguild <- fread("funutax.guilds.txt")
setnames(funguild,"OTU ID", "OTU")
#===============================================================================
#       differential analysis
#===============================================================================
# p value for FDR cutoff
alpha <- 0.05

# # the model
design <- ~Cultivar:Location+Cultivar*Status

dds$Status <- factor(dds$Status,levels(dds$Status)[2:1])

# # add design to dds object
design(dds) <- design

# # run model
dds <- DESeq(dds,parallel=T)

 # results  
res_LM28 <- results(dds,alpha=alpha,parallel=F,contrast = c("Status","Diseased","Healthy"))
res_ZHM2 <- results(dds,alpha=alpha,parallel=F,contrast=list(c("Status_Diseased_vs_Healthy","CultivarZHM2.StatusDiseased")))
res_int <- results(dds,alpha=alpha,name="CultivarZHM2.StatusDiseased")

sink(paste0(RHB,"_results_summary.txt"))
 print("LM28")
 summary(res_LM28)
 print("ZHM2")
 summary(res_ZHM2)
 print("interaction effect")
 summary(res_int)
sink()

res.merge <- data.table(inner_join(data.table(OTU=rownames(res_LM28),as.data.frame(res_LM28)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,paste0(RHB,"_res_LM28.txt"),sep="\t",na="",quote=F)

res.merge <- data.table(inner_join(data.table(OTU=rownames(res_ZHM2),as.data.frame(res_ZHM2)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,paste0(RHB,"_res_ZHM2.txt"),sep="\t",na="",quote=F)

res.merge <- data.table(inner_join(data.table(OTU=rownames(res_int),as.data.frame(res_int)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,paste0(RHB,"_res_int.txt"),sep="\t",na="",quote=F)

# combined status
dds$pos <- as.factor(paste0(dds$Cultivar,dds$Location))
design <- ~pos+Status
design(dds) <- design
dds <- DESeq(dds,parallel=T)

res <- results(dds,alpha=alpha,parallel=F,contrast = c("Status","Diseased","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,paste0(RHB,"_res.txt"),sep="\t",na="",quote=F)

res.merge <- data.table(inner_join(data.table(OTU=rownames(fun_des$res_LM28),as.data.frame(fun_des$res_LM28)),data.table(OTU=rownames(fun_des$taxData),fun_des$taxData)))[funguild[,-"taxonomy"],on="OTU"]
fwrite(res.merge,paste0(RHB,"_res_LM28_guilded.txt"),sep="\t",na="",quote=F)

res.merge <- data.table(inner_join(data.table(OTU=rownames(fun_des$res_ZHM2),as.data.frame(fun_des$res_ZHM2)),data.table(OTU=rownames(fun_des$taxData),fun_des$taxData)))[funguild[,-"taxonomy"],on="OTU"]
fwrite(res.merge,paste0(RHB,"_res_ZHM2_guilded.txt"),sep="\t",na="",quote=F)

res.merge <- data.table(inner_join(data.table(OTU=rownames(fun_des$res),as.data.frame(fun_des$res)),data.table(OTU=rownames(fun_des$taxData),fun_des$taxData)))[funguild[,-"taxonomy"],on="OTU"]
fwrite(res.merge,paste0(RHB,"_res_guilded.txt"),sep="\t",na="",quote=F)

#===============================================================================
#       Venn diagrams
#===============================================================================

#### GET DDS BACK RESULTS FIRST ####

LM28LM28 <- res_LM28$OTU[which(res_LM28$padj<=0.05&res_LM28$log2FoldChange<0)]
ZHM2 <- res_ZHM2$OTU[which(res_ZHM2$padj<=0.05&res_ZHM2$log2FoldChange<0)]
Combined <- res$OTU[which(res$padj<=0.05&res$log2FoldChange<0)]

#pdf(paste0(RHB,"_Euler_down.jpg"),width=9,height=9)

inp <- list(LM28=LM28,ZHM2=ZHM2,Both=Combined)

venn_bd <-venn.diagram(inp,filename=NULL,euler.d=F,scaled=F,
             col = "transparent",
             fill = cbPalette[2:4],
             alpha = 0.5,
             label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
             cex = 2.5,
             fontfamily = "serif",
             fontface = "bold",
             cat.default.pos = "text",
             cat.col = c("darkred", "darkblue", "darkgreen"),
             cat.cex = 2.5,
             cat.fontfamily = "serif",
             cat.dist = c(0.06, 0.06, 0.03),
             cat.pos = 0)

LM28 <- res_LM28$OTU[which(res_LM28$padj<=0.05&res_LM28$log2FoldChange>0)]
ZHM2 <-  res_ZHM2$OTU[which(res_ZHM2$padj<=0.05&res_ZHM2$log2FoldChange>0)]
Combined <- res$OTU[which(res$padj<=0.05&res$log2FoldChange>0)]
inp <- list(LM28=LM28,ZHM2=ZHM2,Both=Combined)

venn_bu <-venn.diagram(inp,filename=NULL,euler.d=F,scaled=F,
                       col = "transparent",
                       fill = cbPalette[2:4],
                       alpha = 0.5,
                       label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
                       cex = 2.5,
                       fontfamily = "serif",
                       fontface = "bold",
                       cat.default.pos = "text",
                       cat.col = c("darkred", "darkblue", "darkgreen"),
                       cat.cex = 2.5,
                       cat.fontfamily = "serif",
                       cat.dist = c(0.06, 0.06, 0.03),
                       cat.pos = 0)



res_LM28 <- fread("FUN_res_LM28.txt")
res_ZHM2 <- fread("FUN_res_ZHM2.txt")
res <- fread("FUN_res.txt")

LM28 <- rownames(fun_des[[2]])[which(fun_des[[2]]$padj<=0.05&fun_des[[2]]$log2FoldChange>0)]
ZHM2 <- rownames(fun_des[[3]])[which(fun_des[[3]]$padj<=0.05&fun_des[[3]]$log2FoldChange>0)]
Combined <- rownames(fun_des[[4]])[which(fun_des[[4]]$padj<=0.05&fun_des[[4]]$log2FoldChange>0)]
inp <- list(LM28=LM28,ZHM2=ZHM2,Both=Combined)

venn_fu <-venn.diagram(inp,filename=NULL,euler.d=F,scaled=F,
                       col = "transparent",
                       fill = cbPalette[2:4],
                       alpha = 0.5,
                       label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
                       cex = 2.5,
                       fontfamily = "serif",
                       fontface = "bold",
                       cat.default.pos = "text",
                       cat.col = c("darkred", "darkblue", "darkgreen"),
                       cat.cex = 2.5,
                       cat.fontfamily = "serif",
                       cat.dist = c(0.06, 0.06, 0.03),
                       cat.pos = 0)



  ggsave("venn_all.jpg",plot_grid(venn_bu,venn_bd,venn_fu,nrow=1,labels = "AUTO",label_size = 25),width=9,height=9,dpi=600)
  ggsave("venn_fun_up.jpg",plot_grid(venn_fu),width=9,height=9)
  ggsave("venn_bac_up.jpg",plot_grid(venn_bu),width=9,height=9)
  ggsave("venn_bac_down.jpg",plot_grid(venn_bd),width=9,height=9)


res_LM28 <- fread("FUN_res_LM28.txt")
res_ZHM2 <- fread("FUN_res_ZHM2.txt")
res <- fread("FUN_res_guilded.txt")

LM28 <- res_LM28$OTU[which(res_LM28$padj<=0.05&res_LM28$log2FoldChange>0)]
ZHM2 <- res_ZHM2$OTU[which(res_ZHM2$padj<=0.05&res_ZHM2$log2FoldChange>0)]
Combined <- res$OTU[which(res$padj<=0.05&res$log2FoldChange>0)]

venn <- venn.out(LM28,ZHM2,Combined)
venn <- lapply(venn,as.data.table)
lapply(venn,function(DT) DT[,OTU:=V1])
venn <- lapply(seq_along(venn),function(i)setnames(venn[[i]],"V1",names(venn)[[i]]))
#lapply(venn,setnames,"OTU")
venn$OTUs <- data.table(OTU=res$OTU)
v2 <- Reduce(function(...) {merge(..., all = T)}, venn)
fwrite(v2,"shared_OTUs_FUN",sep="\t",na="",quote=F)
i <- c(597,736,855,556,729,595,556)
venn.fun <- function(i,Title) {
  plot.new()
 # title(Title)+
  draw.triple.venn(i[1],i[2],i[3],i[4],i[5],i[6],i[7],
                 category = c("LM82","ZHM2","Both"),
                 col = "transparent",
                 euler.d =F, scaled = F,
                 fill = cbPalette[2:4],
                 alpha = 0.5,
                 label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
                 cex = 1.5,
                 fontfamily = "serif",
                 fontface = "bold",
                 cat.default.pos = "text",
                 cat.col = c("darkred", "darkblue", "darkgreen"),
                 cat.cex = 1.5,
                 cat.fontfamily = "serif",
                 cat.dist = c(0.06, 0.06, 0.03),
                 cat.pos = 0,margin=0.1)
}

res <- fread("FUN_res.txt")
res.merged <- fread("FUN_res_guilded_2.txt",check.names = T)
res.merged <- res.merged[OTU%in%res$OTU,]

res.merged[,.N,by=.(Trophic.Mode)]
res.merged[,c("LM28","ZHM2","BOTH"):=0]
res.merged$LM28[(res.merged$LM28.padj<=0.05&res.merged$LM28.FC>0)] <- 1
res.merged$ZHM2[(res.merged$ZHM2.padj<=0.05&res.merged$ZHM2.FC>0)] <- 1
res.merged$BOTH[(res.merged$res.padj<=0.05&res.merged$res.log2>0)] <- 1

cols <- c("LM28","ZHM2","BOTH")
col <- "Trophic.Mode"
col <- "Guild"
output <- res.merged[,.N,by=c(col,cols)]

v.i <- list(
  LM28= output[LM28==1,sum(N),by=col],
  ZHM2= output[ZHM2==1,sum(N),by=col],
  BOTH= output[BOTH==1,sum(N),by=col],
  AuB=  output[LM28==1&ZHM2==1,sum(N),by=col],
  BuC=  output[ZHM2==1&BOTH==1,sum(N),by=col],
  AuC=  output[LM28==1&BOTH==1,sum(N),by=col],
  AuBuC=output[LM28==1&ZHM2==1&BOTH==1,sum(N),by=col] 
)

v.i <- lapply(seq_along(v.i),function(i)setnames(v.i[[i]],"V1",names(v.i)[[i]]))
v.i <- Reduce(function(...) {merge(..., all = T)}, v.i)
v.i[[col]][v.i[[col]]=="-"] <- "Unknown"
v.i[,(col):=lapply(.SD,as.factor),.SDcols=col]

#v.i[,Trophic.Mode:=factor(Trophic.Mode,levels=levels(as.factor(Trophic.Mode))[c(8,1,5,7,2,4,6,3)])]
#v.i <- v.i[c(1,2,6,8,3,5,7,4),]

# 
# g <- lapply(seq_along(v.i),function(i)venn.fun(unlist(v.i[i,-1]),unlist(v.i[i,1])))
# 
# ggplot("guild.jpg",plot_grid(plotlist=g,labels = unlist(v.i[,1]),label_x = 0.5,hjust = 0.5))
# ggplot()
# 
# 
# plot_grid(plotlist=lapply(seq_along(v.i),function(i)venn.fun(unlist(v.i[i,-1]),unlist(v.i[i,1]))))
# 
# i<-2
# venn.fun(i=as.numeric(unlist(v.i[i,-1])),Title=unlist(v.i[i,1]))
#plot_grid(plotlist=g,labels = unlist(v.i[,1]),


##### Fisher exact tests #####

vv <- output[,sum(N),by=col]
vv[[col]][vv[[col]]=="-"] <- "Unknown"
v.i <- v.i[vv,on=col] 
tots <- colSums(v.i[,-1],na.rm=T)

v.i[is.na(v.i)] <- 0
v.i[,rem:=LM28+ZHM2+BOTH]
v.i <- v.i[rem!=0,]
#i=2;j=3

test <-as.data.table(t(sapply(1:nrow(v.i),function(i){
  unlist(lapply(1:3,function(j){
    as.data.table(fisher.test(rbind(unlist(v.i[i,c(9,(j+1)),with=F]),tots[c(8,j)]))[c(1,3)])
  }))
})))

# Guilds
test[,Guild:=v.i$Guild]
v.i <- v.i[test,on="Guild"]
fwrite(v.i[,c(1,9,2,13,12,3,15,14,4,17,16)],"guilds.tsv",sep="\t")
  
# Trophic.Mode
test[,Trophic.Mode:=v.i$Trophic.Mode]
v.i <- v.i[test,on="Trophic.Mode"]
fwrite(v.i[,c(1,9,2,13,12,3,15,14,4,17,16)],"Trophic.Mode.tsv",sep="\t")


# fwrite(test,"fisher_test.tsv",sep="\t",row.names = T)
# fwrite(v.i,"fun_group_numbers.tsv",sep="\t")
# 
# 
# lapply(seq_along(test),function(i)test[[i]][,guild:=v.i[i,][[col]]])
# 
# names(test) <- v.i[[col]]#[1:8]
# 
# test2 <- Reduce(function(...) {rbind(...)},test)
# 
# test2 <- Reduce(function(...) {merge(..., all = T)},test)
# rownames(test) <- v.i[[col]]#[1:8]


#sink("fisher_test.tsv",sep="\t")
#test
#sink()

#===============================================================================
#       Metacoder
#===============================================================================
##### Create Metacoder obj #####

# colData needs to be modified for (easy) use with metacoder   
  colData <- colData(dds)
  colData$Sample_id <- rownames(colData)
  
  # set taxData columns with confidencer less than conf to unkown - and add an extra column for unique names for identical species 
  taxData <- phyloTidy(taxData,conf=0.65,setCount = T)
  
  #create metacoder object
  OBJ <- parse_phyloseq(phyloseq(otu_table(counts(dds), taxa_are_rows = T), tax_table(as.matrix(taxData[rownames(dds),]))))
  
  # calculate abundances of taxons above species level
  OBJ$data$tax_abund <- calc_taxon_abund(OBJ, "otu_table",cols = colData$sample_id)
  
  
  #####   Metacoder DESeq analysis #####
  
  # Run deseq with metacoder at given rank
  OBJ$data$dds <- MCDESeq(OBJ,~Cultivar:Location+Cultivar*Status,"genus",colData)
  
  # Get results
  OBJ$data$diff_table_lm28 <- MCDESres(OBJ,contrast = c("Status","Diseased","Healthy"),alpha=0.05)
  OBJ$data$diff_table_zhm2 <- MCDESres(OBJ,contrast=list(c("Status_Healthy_vs_Diseased","CultivarZHM2.StatusHealthy")),alpha=0.05)
  
  OBJ$data$diff_table_zhm2$log2FoldChange <- OBJ$data$diff_table_zhm2$log2FoldChange*-1
  OBJ$data$diff_table_zhm2$log2FoldChange2 <- OBJ$data$diff_table_zhm2$log2FoldChange2*-1
  
  #res_merge <- as.data.table(OBJ$data$diff_table)


#####  Filter Metacoder #####

# create filter of T/F values to keep - must be same length as 
OBJ$data$diff_table <- OBJ$data$diff_table_zhm2

my_filter_1 <- unlist(OBJ$taxon_ids()) %in%  OBJ$data$diff_table$taxon_id  #OBJ$taxon_ranks()=="species"
my_filter_2 <- (qf(OBJ$data$diff_table$taxonomy)!="REMOVE") &
  (OBJ$data$diff_table$padj<=0.05) & (complete.cases(OBJ$data$diff_table)) #&
#(OBJ$data$diff_table$log2FoldChange<0)

obj <- OBJ %>% 
  filter_taxa(my_filter_1,supertaxa=F) %>%
  filter_taxa(my_filter_2,supertaxa = T) #%>%

# Rename tips - watch out these are R6 objects need to clone before copying
obj$taxa <- lapply(obj$taxa,function(x) {
  x <- x$clone(deep=T)
  x$name$name<-sub("unknown.*","",x$name$name)
  x$name$name<-sub("[0-9].*","",x$name$name)
  x$name$name<-sub("cls ","",x$name$name)
  x
})

# Rename BACTERIA tips - watch out these are R6 objects need to clone before copying
obj$taxa <- lapply(obj$taxa,function(x) {
  x <- x$clone(deep=T)
  x$name$name<-sub("unknown.*","",x$name$name)
  x$name$name<-sub("candi","Candi",x$name$name)
  x$name$name<-sub("Candidatus","",x$name$name)
  x
})


# Remove tip names below "rank" 
obj$taxa <- lapply(obj$taxa,function(x) {
  x <- x$clone(deep=T)
  x$name$name[which(names(taxData)==x$rank$name)>3] <- ""
  x
})


##### Heat trees #####

my_range <- c(low="#D9944D",mid="#DDDDDD", high="#34B8A4")

set.seed(sum(utf8ToInt("Xiangming Xu")))
g <-heat_tree(obj,
            node_size = res_merge[taxon_id%in%obj$taxon_ids(),]$baseMean, 
            #log(t2[t2>100],2), # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_label = taxon_names,
            node_label_size_range = c(0.014,0.028),
            node_color = log2FoldChange2, # A column from `obj$data$diff_table`
            node_color_range = my_range, #diverging_palette(), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            node_color_interval = c(0, 5), # The range of `log2_median_ratio` to display
            edge_color_interval = c(0, 5), # The range of `log2_median_ratio` to display
            node_size_axis_label = "Abundance",
            node_color_axis_label = "Log(2) fold change",#bquote(~Log[2]~' fold change'),
            #initial_layout = "re", layout = "da",
            make_node_legend = F, make_edge_legend = F,
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford"#, # The layout algorithm that initializes node locations
            #output_file = "/home/greg/heat_tree_test.jpg" # Saves the plot as a pdf file
  ) 

g

#g_lm28 <- g
g_zhm2 <- g

plot_grid(g_lm28,g_zhm2,nrow=1)


legend <- get_legend(ggplot(data.frame(x=c(-6:6)),aes(x=x,y=x,colour=x))+ geom_point() + 
                       scale_color_gradient2(low=my_range[1],mid=my_range[2],high=my_range[3],name=bquote(~Log[2]~' fold change'))+
                       theme(legend.position="bottom",legend.justification = "left")+
                       guides(colour = guide_colourbar(barwidth = 10, barheight = 1,ticks = FALSE,title.position = "top")))

plot_grid(g_lm28,g_zhm2,legend,nrow=2,rel_heights = c(1,0.1))
  ggsave("heat_tree_bac.jpg",plot_grid(g_lm28,g_zhm2,legend,nrow=2,rel_heights = c(1,0.1)))
  ggsave("heat_tree_bac_lm28.jpg",plot_grid(g_lm28,legend,nrow=2,rel_heights = c(1,0.1)))
  ggsave("heat_tree_bac_zhm2.jpg",plot_grid(g_zhm2,legend,nrow=2,rel_heights = c(1,0.1)))





