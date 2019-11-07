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

ubiom_BAC <- loadData("BAC.otus_table.txt","colData","BAC.taxa",RHB="BAC")
# ubiom_BAC$countData <- ubiom_BAC$countData[,colnames(ubiom_BAC$countData)%in%rownames(ubiom_BAC$colData)]
ubiom_FUN <- loadData("FUN.otus_table.txt","colData","FUN.taxa",RHB="FUN")
# ubiom_FUN$countData <- ubiom_FUN$countData[,colnames(ubiom_FUN$countData)%in%rownames(ubiom_FUN$colData)]

#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- combineTaxa("FUN.taxa")
countData <- combCounts(combinedTaxa,countData)
taxData <- combTaxa(combinedTaxa,taxData)
ubiom_FUN$countData <- countData
ubiom_FUN$taxData <- taxData

#===============================================================================
#       Create DEseq objects
#===============================================================================
		   
ubiom_FUN$dds <- ubiom_to_des(ubiom_FUN,filter=expression(colSums(countData)>=1000&colData$Block!="R"))
ubiom_BAC$dds <- ubiom_to_des(ubiom_BAC,filter=expression(colSums(countData)>=1000&colData$Block!="R"))
# ubiom_OO$dds <- ubiom_to_des(ubiom_OO,filter=expression(colSums(countData)>=1000&colData$Block!="R"),calcFactors=geoMeans)
# ubiom_NEM$dds <- ubiom_to_des(ubiom_NEM,filter=expression(colSums(countData)>=500&colData$Block!="R"))

#===============================================================================
#      ****FUNGI/BACTERIA****
#===============================================================================
# Fungi
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Sample rarefaction plots
#===============================================================================        

library(grid)
library(gridExtra)
library(viridis)
				    
gfunc <- function(countData,coldata,title) {        
  colData <- colData[names(countData),]

  # remove low count and control samples
  myfilter <- colData$Treatment!="Control"

  # apply filter
  colData <- droplevels(colData[myfilter,])
  countData <- countData[,myfilter]

  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))

  # get cummulative sum of each sample
  DT <- cumsum(DT)    

  # log the count values                            
  DT <- log10(DT)

  # relabel columns
  colnames(DT) <- sub("(X)([0-9]+[HS])(.*)","\\2",colnames(DT))

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

invisible(mapply(assign, names(ubiom_BAC), ubiom_BAC, MoreArgs=list(envir = globalenv())))
g1 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Bacteria")
invisible(mapply(assign, names(ubiom_FUN), ubiom_FUN, MoreArgs=list(envir = globalenv())))
g2 <- gfunc(as.data.frame(counts(dds)),as.data.frame(colData(dds)),"Fungi")

glegend <- get_legend(g) # this won't work
ggsave("rarefaction_all.pdf",grid.arrange(g1,g2,left=textGrob(label=expression("Log"[10] * " aligned sequenecs"),rot=90),bottom="OTU count",nrow=2))                              


#===============================================================================
#       Alpha diversity analysis
#===============================================================================

# Recreate dds object and don't filter for low counts before running Alpha diversity

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
plot_alpha(counts(dds,normalize=T),colData(dds),design="Treatment",colour=NULL,measures=c("Chao1", "Shannon", "Simpson","Observed"),type="box")		   
		   
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

dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
plotOrd(d,colData,design="parent",shape="sex",alpha=0.75,cbPalette=T)
plotOrd(d,colData,design="parent",shape="sex",axes=c(2,3),alpha=0.75,cbPalette=T)














#===============================================================================
#       differential analysis
#===============================================================================

# p value for FDR cutoff
alpha <- 0.1

# the model
design <- ~sex*parent

# add design to dds object
design(dds) <- design

# run model
dds <- DESeq(dds,parallel=F)

# difference between sexes
res <- results(dds,alpha=alpha,parallel=F,contrast=c("sex","M","F"))
		   
# difference between 

# results COD vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm1 <- res.merge
       
# results AOD vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm2 <- res.merge

# results Remission vs Healthy
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","Remission","Healthy"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"Rem_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm3 <- res.merge

# results COD vs AOD
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","AOD"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_AOD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm4 <- res.merge

# results COD vs Remission
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","COD","Remission"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_REM_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm5 <- res.merge

# results AOD vs Remission
res <- lapply(list_dds,results,alpha=alpha,parallel=T,contrast=c("condition","AOD","Remission"))
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"AOD_REM_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
rm6 <- res.merge
     
       
all.tabs <- lapply(seq(1,4),function(i) data.table(list(
	rm1[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_COD = log2FoldChange)  %>% rename(padj_COD = padj),
	rm2[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_AOD = log2FoldChange)  %>% rename(padj_AOD = padj),
	rm3[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_REM = log2FoldChange)  %>% rename(padj_REM = padj),
	rm4[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_CA  = log2FoldChange)  %>% rename(padj_CA  = padj),
	rm5[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_CR  = log2FoldChange)  %>% rename(padj_CR  = padj),
	rm6[[i]] %>% select(-lfcSE,-stat,-pvalue) %>% rename(FC_AR  = log2FoldChange)  %>% rename(padj_AR  = padj)
        ) %>% Reduce(function(dt1,dt2) full_join(dt1,dt2), .))
)	
sapply(seq(1,4),function(i) write.table(all.tabs[[i]],paste(RHB,names(res.merge)[i],"diffs.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# chestnuts and bigwood
res <- lapply(list_dds[5:6],results,alpha=alpha)
res.merge <- lapply(res,function(res)data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData))))
sapply(names(res.merge),function(x) write.table(res.merge[[x]],paste(RHB,x,"COD_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F))
       
# output sig fasta
writeXStringSet(readDNAStringSet(paste0(RHB,".otus.fa"))[ res.merge[padj<=0.05]$OTU],paste0(RHB,".sig.fa")) 

# Merge

rm <- lapply(rm2,function(rm) rm[,c(1,2,3,7),drop=F])      
suffixes <- paste0("_",names(rm))
test <- merge(merge(rm[[1]],rm[[2]],by="OTU",all=T,suffixes=suffixes[1:2]),
	     merge(rm[[3]],rm[[4]],by="OTU",all=T,suffixes=suffixes[3:4]),
	     all=T )
res.merge <- data.table(inner_join(test,data.table(OTU=rownames(taxData),taxData)))
fwrite(res.merge,"AOD.diff.txt",sep="\t",quote=F,na="") 
