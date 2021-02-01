#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library(BiocParallel)
library(data.table)
library(plyr)
library(tidyverse)
library(Biostrings)
library(vegan)
library(lmPerm)
library(phyloseq)
library(ape)


library(devtools)
install_github("eastmallingresearch/Metabarcoding_pipeline/scripts")
library(metafuncs)

# set if using multiprocessing - Linux only
register(MulticoreParam(2))

# colour blind palette
cbPalette <-c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#===============================================================================
#       Load data
#===============================================================================

# load otu count table
countData <- read.table("BAC.otus_table.txt",header=T,sep="\t",row.names=1, comment.char = "")

# load sample metadata
colData <- read.table("colData",header=T,sep="\t",row.names=1)

# load taxonomy data
taxData <- read.table("BAC.taxa",header=F,sep=",",row.names=1)

# reorder columns
taxData<-taxData[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]

# add best "rank" at 0.65 confidence and tidy-up the table
taxData<-phyloTaxaTidy(taxData,0.75)

# get unifrac dist
#phylipData <- fread.phylip("BAC.phy")

#njtree <- nj(as.dist(phylipData))

# save data into a list
ubiome_BAC <- list(
  countData=countData,
  colData=colData,
  taxData=taxData,
  RHB="BAC"
)

# change names to match between metadata and sequence data
rownames(ubiome_BAC$colData) <- sub(".*_L1_","",rownames(ubiome_BAC$colData))
names(ubiome_BAC$countData) <- sub(".*_L1_","",names(ubiome_BAC$countData))

# Fungi all in one call
ubiome_FUN <- list(
  countData=read.table("FUN.otus_table.txt",header=T,sep="\t",row.names=1,comment.char = ""),
  colData=read.table("colData",header=T,sep="\t",row.names=1),
  taxData=phyloTaxaTidy(read.table("FUN.taxa",header=F,sep=",",row.names=1)[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)],0.65),
  RHB="FUN"
)
rownames(ubiome_FUN$colData) <- sub(".*_L1_","",rownames(ubiome_FUN$colData))
names(ubiome_FUN$countData) <- sub(".*_L1_","",names(ubiome_FUN$countData))



#===============================================================================
#       Combine species
#===============================================================================

#### combine species at 0.95 (default) confidence (if they are species)

# Fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
combinedTaxa <- metafuncs::combineTaxa("FUN.taxa")
countData <- metafuncs::combCounts(combinedTaxa,countData)
taxData <- metafuncs::combTaxa(combinedTaxa,taxData)
ubiome_FUN$countData <- countData
ubiome_FUN$taxData <- taxData

#===============================================================================
#      ****FUNGI/BACTERIA****
#===============================================================================

#### RUN ONLY ONE OF THESE #####
# copies objects within a list (ubiom_FUN/BAC) into the glogal environment

# Fungi
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

#===============================================================================
#       Create DEseq objects
#===============================================================================

# ensure colData rows and countData columns have the same order
colData <- colData[names(countData),]

# set a minimal design
design<-~1

#create DES object
dds<-DESeqDataSetFromMatrix(countData,colData,design)

# set the size factors - this will change when qPCR is done
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Alpha diversity analysis - RUN BEFORE FILTERING OUT ANY LOW COUNT OTUS
#===============================================================================

# plot alpha diversity - plot_alpha will convert normalised abundances to integer values
g <- metafuncs::plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",measures=c("Chao1", "Shannon", "Simpson","Observed"),limits=c(500,10000,"Chao1"))
g

# save plot
ggsave(paste(RHB,"Alpha.pdf",sep="_"),g)


### permutation based anova on diversity index ranks ###

# get the diversity index data
all_alpha_ord <- metafuncs::plot_alpha(counts(dds,normalize=T),colData(dds),design="site",colour="condition",returnData=T)

# join diversity indices and metadata
dds$Samples <- rownames(colData(dds)) # or could use tibble/rownames construct in join by syntax)
all_alpha_ord <- as.data.table(left_join(all_alpha_ord,as.data.frame(colData(dds))))

# remove specultation site (not necessary as I now have 
# all_alpha_ord <- all_alpha_ord[site!="Speculation",]

# perform anova for each index (this may need editing as the design will be unbalanced)
sink(paste(RHB,"ALPHA_stats.txt",sep="_"))
setkey(all_alpha_ord,S.chao1)
print("Chao1")
summary(aovp(as.numeric(as.factor(all_alpha_ord$S.chao1))~condition*site,all_alpha_ord))
setkey(all_alpha_ord,shannon)
print("Shannon")
summary(aovp(as.numeric(as.factor(all_alpha_ord$shannon))~condition*site,all_alpha_ord))
setkey(all_alpha_ord,simpson)
print("simpson")
summary(aovp(as.numeric(as.factor(all_alpha_ord$simpson))~condition*site,all_alpha_ord))
sink()

#===============================================================================
#       Filter data
#============================================================================

# plot cummulative reads (will also produce a data table "dtt" in the global environment)
# ggsave(paste(RHB,"OTU_counts.pdf",sep="_"),plotCummulativeReads(counts(dds,normalize=T)))

# remove OTUs with < 5 total read (could be a bit more agressive here)
dds <- dds[rowSums(counts(dds, normalize=T))>4,]

#===============================================================================
#       Plot/table frequencies
#===============================================================================

# bacteria
#invisible(mapply(assign, names(ubiome_BAC_Rhiz), ubiome_BAC_Rhiz, MoreArgs=list(envir = globalenv())))
rhiz <- metafuncs::sumTaxa(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),
                conf=0.8,
                proportional=T,
                design="site",
                taxon = "order"
          )

rhiz <- metafuncs::sumTaxaAdvanced(list(as.data.frame(counts(dds,normalize=T)),taxData,colData),
                conf=0.8,
                proportional=T,
                design="site",
                taxon = "phylum",
                others=T
)

# fungi
levels(rhiz$order) <- gsub("_.*\\(","\\(",levels(rhiz$order))
#levels(rhiz$order) <- gsub("_.*","",levels(rhiz$order))


md1 <- melt(rhiz,id=colnames(rhiz)[1])

md1[[1]]  <- factor(md1[[1]], levels=unique(md1[[1]][order(md1$value,decreasing=T)]),ordered=T)
md1$value <- as.numeric(md1$value)
names(md1)[2] <- "Site"
levels(md1$Site)[4] <- "Gt Monk"
#1015 881

sites <- c("Attingham","Gt Monk","Langdale","Winding")
scaleFUN<-function(x) sprintf("%.0f", x)

g <- ggplot(md1[md1$Site%in%sites,],aes(x=order,y=value,fill=Site)) + theme_classic_thin()
g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
g <- g + scale_fill_manual(values=cbPalette)
g <- g  + xlab("")+ ylab("Frequency (%)")
g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=NULL)
g <- g + guides(fill=guide_legend(title =""))
g <- g + theme_blank()
g <- g + theme(legend.position = "bottom",legend.direction = "horizontal",legend.justification = "left",
  axis.text.x = element_text(angle = 45, hjust = 1,size=14),
  plot.margin=unit(c(0.2,0,0.2,1.5),"cm"),
  axis.line.y = element_line(colour = "black",size=1),
  axis.ticks.x=element_blank(),
  text=element_text(size=14),
  axis.title.y=element_text(size=(14-2)))
g
#ggsave("BAC_frequency.pdf",g,width=8,height=7)

#### identifiable to each level props

invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs=list(envir = globalenv())))
# Bacteria
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs=list(envir = globalenv())))

taxData <- as.data.table(taxData,keep.rownames = "OTU")

cols <- names(taxData)[grep("factor",sapply(taxData,class))]

taxData[,(cols):=lapply(.SD,as.number),.SDcols=cols]

nrow(taxData[p_conf>=0.8,])/nrow(taxData)

qf <- function(taxon,taxData,conf) sum(taxData[[taxon]]>=as.numeric(conf))/nrow(taxData)

cols_rank <- names(taxData)[grep("character",sapply(taxData,class))][c(-1,-9)]
conf <- c(0.5,0.6,0.7,0.8)

sites <- c(Attingham="Attingham",GT_Monk="Gt_Monk",Langdale="Langdale",Winding="Winding")

lapply(sites,function(s) {
  cd <- countData[,colData[colData$site==s,]$id]
  cd <- cd[rowSums(cd)>0,]
  otus <- rownames(cd)
  td <- copy(taxData)
  td <- td[OTU%in%otus,]
  X <- t(sapply(conf,function(i) unlist(lapply(cols,qf,td,i))))
  colnames(X) <- cols_rank
  rownames(X) <- conf
  X
})

#X <- t(sapply(conf,function(i) unlist(lapply(cols,qf,taxData,i))))
#colnames(X) <- cols_rank
#rownames(X) <- conf


#===============================================================================
#       Taxonomy plots
#===============================================================================

td <- taxData

#test <- taxonomyTidy(test[1:2,])
# diversity
#obj <- list(data.frame(cbind(td[,1,drop=F],1)[,2,drop=F]), td,as.data.frame(colData(dds)))

# abundance - with variance stabilizing transformation (probably best)
obj <-list(as.data.frame(counts(dds,normalize=T)),td,as.data.frame(colData(dds)))

# sum at the phylum level (will produce warnings for diversity- can be ignored)
dd <- metafuncs::sumTaxa( 
  obj,
  design = "site",
  conf=0.8,
  proportional=T,
  taxon="class"
)

# Combine and filter low abundance groups (for designed exp)
all <- rowSums(dd[,-1])/(ncol(dd)-1)

#all <- dd$all

dc <- rbind(dd[all>=1.1,],c("others",colSums(dd[all<1.1,-1])))
md <- melt(dc,id=colnames(dd)[1])

#levels(md$variable)[1] <- "OTUs"
md$value <- as.numeric(md$value)
#md$Normalisation <- md$variable

#fungi
#md[[1]] <- sub(".*_"," ",md[[1]])
# bacteria
#md[[1]] <- sub("_Na5|_FOP1","",sub("\\/.*","",sub("_.*\\._"," ",md[[1]])))
#tef
#md[[1]] <- sub("\\/.*","",sub("_.*","",sub("_fsp_"," ",md[[1]])))


md[[1]]  <- factor(md[[1]], levels=c(unique(md[[1]][md[[1]]!="others"][order(md$value[md[[1]]!="others"],decreasing=T)]),"others"),ordered=T)
colnames(md)[1]<-"taxon"

levels(md$variable)[levels(md$variable )=="Gt_Monk"] <- "Gt.Monk"

md$variable <- factor(md$variable,levels=levels(md$variable)[c(4,1,3,2,5:7)])
# split variable column
#md <- separate(md,"variable",into=c("sample","dilution"))

g <- ggplot(md,aes(x=taxon,y=value))
g <- g + geom_bar(stat="identity",colour="black")
g <- g  + xlab("")+ ylab("")
scaleFUN<-function(x) sprintf("%.0f", x)
g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=NULL)
g <- g + guides(fill=guide_legend(ncol=1))
#g <- g+ facet_grid(dilution~sample) + theme_blank() 
g <- g+ facet_wrap(~variable,nrow=6) + theme_blank() 
g <- g + theme(
  strip.background = element_rect(colour="white", fill="white"),
  strip.text = element_text(size=16),
  axis.text.x = element_text(angle = 45, hjust = 1,size=14),
  plot.margin=unit(c(0.2,0,0.2,1.5),"cm"),
  #  panel.spacing.y = unit(0.5, "lines"),
  panel.spacing.x = unit(2, "lines"),
  axis.line.y = element_line(colour = "black",size=1),
  axis.ticks.x=element_blank(),
  text=element_text(size=12),
  plot.title = element_text(hjust = -0.11),
  axis.title.y=element_text(size=(14-2)))
g

ggsave(paste0(RHB,".phylum.pdf"),g)
ggsave(paste0(RHB,".class.pdf"),g)


#===============================================================================
#       Beta diversity PCA/NMDS
#===============================================================================

### PCA ###

# perform PC decomposition of DES object
mypca <- metafuncs::des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# plot the PCA
metafuncs::plotOrd(d,colData,design="Condition",xlabel="PC1",ylabel="PC2")
metafuncs::plotOrd(d,colData,shape="Condition",design="Location",continuous=T,xlabel="PC1",ylabel="PC2")
dev.off()

#===============================================================================
#       differential analysis
#===============================================================================

# remove extra rows from colData
dds <- dds[,dds$condition!="Sandra"]

# p value for FDR cutoff
alpha <- 0.05

##### COMPARISON - ALL SITES######
dds2 <- dds
levels(dds$condition)[1] <- "Symptom"
levels(dds$condition)[2] <- "Symptom"
design <- ~site + condition
design(dds) <- design
dds <- DESeq(dds,parallel=T)
res <- results(dds,alpha=alpha,parallel=F,contrast=c("condition","Symptom","Healthy"))
res_merge <-  data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
summary(res)

##### COMPARISON WITH METAGENOMICS - ALL SITES######


# the model
design <- ~site + condition
#or
design <- ~block_pair + site + condition

dds$condition <- droplevels(dds$condition)
dds$block_pair <- as.factor(dds$block_pair)
dds$site  <- droplevels(dds$site)

# add full model to dds object
design(dds) <- design

# calculate fit - parallel only useful for subbins
dds_res <- DESeq(dds,parallel=T)

healthy_symptom <- results(dds_res,alpha=alpha,parallel=F,contrast=c("condition","Symptom","Healthy"))
healthy_AOD <- results(dds_res,alpha=alpha,parallel=F,contrast=c("condition","AOD","Healthy"))
healthy_COD <- results(dds_res,alpha=alpha,parallel=F,contrast=c("condition","COD","Healthy"))
healthy_COD_all <-  results(dds_res,alpha=alpha,parallel=F,
                            contrast=list(c("condition_COD_vs_AOD","condition_Symptom_vs_AOD"),c("condition_Healthy_vs_AOD","condition_Healthy_vs_AOD")),
                            listValues=c(.5,-.5))



sink(paste0(RHB,"metabarcoding.noblock.dds_res.txt"))
print("COD_Sites")
summary(healthy_symptom)
print("AOD")
summary(healthy_AOD)
print("COD")
summary(healthy_COD)
print("COD_all")
summary(healthy_COD_all)
sink()

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_symptom),as.data.frame(healthy_symptom)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"healthy_symptom.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_AOD),as.data.frame(healthy_AOD)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"healthy_AOD.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_COD),as.data.frame(healthy_COD)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"healthy_COD.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_COD_all),as.data.frame(healthy_COD_all)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"healthy_COD_all.txt",sep="."),sep="\t",quote=F,na="")


sink(paste0(RHB,"Blocks.dds_res.txt"))
print("COD_Sites")
summary(healthy_symptom)
print("AOD")
summary(healthy_AOD)
print("COD")
summary(healthy_COD)
print("COD_all")
summary(healthy_COD_all)
sink()

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_symptom),as.data.frame(healthy_symptom)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"blocks.healthy_symptom.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_AOD),as.data.frame(healthy_AOD)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"blocks.healthy_AOD.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_COD),as.data.frame(healthy_COD)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"blocks.healthy_COD.txt",sep="."),sep="\t",quote=F,na="")

res_merge <-  data.table(inner_join(data.table(OTU=rownames(healthy_COD_all),as.data.frame(healthy_COD_all)),data.table(OTU=rownames(taxData),taxData)))
fwrite(res_merge,paste(RHB,"blocks.healthy_COD_all.txt",sep="."),sep="\t",quote=F,na="")


### END ###


# split dds object into per wood
# first get rid of bigwood as it only has 2 samples - and remove it from the levels of site
list_dds <-list(Attingham   = dds[,dds$site=="Attingham"],
		            Gt_Monk     = dds[,dds$site=="Gt_Monk"],
		            Langdale    = dds[,dds$site=="Langdale"],
		            Winding     = dds[,dds$site=="Winding"],
		            Chestnuts   = dds[,dds$site=="Chestnuts"],
		            Bigwood     = dds[,dds$site=="Bigwood"]
)		
# Filter low count OTUs
list_dds <- lapply(list_dds,function(dds)  dds[rowSums(counts(dds, normalize=T))>0,])

# add full model to dds object
list_dds <- lapply(list_dds,function(dds) {
	design(dds) <- design;
	dds <- dds[,dds$condition!="Sandra"];
	colData(dds) <- droplevels(colData(dds));
	dds}
)

# calculate fit
list_dds2 <- lapply(list_dds,DESeq,parallel=T)
list_dds >- list_dds2[1:4]
		   
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

## across sites analysis

# AOD/COD across sites 
design <- ~site + block_pair + condition
dds2 <- dds[,dds$site=="Attingham"|dds$site=="Langdale"|dds$site=="Winding"]
design(dds2) <- design;
dds2 <- dds2[,dds2$condition!="Sandra"];
colData(dds2) <- droplevels(colData(dds2));

# calculate fit
dds2 <- DESeq(dds2,parallel=T)
res <- results(dds2,alpha=alpha,parallel=T,contrast=c("condition","AOD","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
res.merge[padj<=0.1,]
fwrite(res.merge,paste(RHB,"AOD_sites_paired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

res <- results(dds2,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
res.merge[padj<=0.1,]
fwrite(res.merge,paste(RHB,"COD_sites_paired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# COD - including chestnuts wood
dds2 <- dds[,dds$site=="Attingham"|dds$site=="Langdale"|dds$site=="Winding"|dds$site=="Chestnuts"]
design(dds2) <- design;
dds2 <- dds2[,dds2$condition!="Sandra"];
colData(dds2) <- droplevels(colData(dds2))

res <- results(dds2,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
res.merge[padj<=0.1,]
fwrite(res.merge,paste(RHB,"COD_sites_chest_paired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# unpiared AOD
design <- ~site + condition
dds2 <- dds[,dds$site!="Bigwood"&dds$site!="Chestnuts"&dds$site!="Speculation"]
design(dds2) <- design;
dds2 <- dds2[,dds2$condition!="Sandra"];
colData(dds2) <- droplevels(colData(dds2));

# calculate fit
dds2 <- DESeq(dds2,parallel=T)
res <- results(dds2,alpha=alpha,parallel=T,contrast=c("condition","AOD","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
res.merge[padj<=0.1,]
fwrite(res.merge,paste(RHB,"AOD_sites_unpaired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)

# unpaired COD
design <- ~site + condition
dds2 <- dds
design(dds2) <- design;
dds2 <- dds2[,dds2$condition!="Sandra"];
levels(dds2$condition)[levels(dds2$condition)=="Symptom"] <- "COD"
colData(dds2) <- droplevels(colData(dds2));

# calculate fit
dds2 <- DESeq(dds2,parallel=T)
res <- results(dds2,alpha=alpha,parallel=T,contrast=c("condition","COD","Healthy"))
res.merge <- data.table(inner_join(data.table(OTU=rownames(res),as.data.frame(res)),data.table(OTU=rownames(taxData),taxData)))
res.merge[padj<=0.1,]
fwrite(res.merge,paste(RHB,"COD_sites_unpaired_diff.txt",sep="_"),quote=F,sep="\t",na="",row.names=F)
