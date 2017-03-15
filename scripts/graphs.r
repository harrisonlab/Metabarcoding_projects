library(ggplot2)
library(phyloseq)


test1 <- prune_taxa(row.names(gb.res[abs(gb.res$log2FoldChange)>=1&(gb.res$padj<=0.05),]),mybiom_gb)
test2 <- prune_taxa(row.names(hb.res[abs(hb.res$log2FoldChange)>=1&(hb.res$padj<=0.05),]),mybiom_hb)
test3<- prune_taxa(row.names(gf.res[abs(gf.res$log2FoldChange)>=1&(gf.res$padj<=0.05),]),mybiom_gf)
test4<- prune_taxa(row.names(hf.res[abs(hf.res$log2FoldChange)>=1&(hf.res$padj<=0.05),]),mybiom_hf)

sample_data(test1)$condition <- as.character(sample_data(test1)$condition)
sample_data(test2)$condition <- as.character(sample_data(test2)$condition)
sample_data(test1)$condition[sample_data(test1)$condition=="Y"] <- "Goatham tree"
sample_data(test1)$condition[sample_data(test1)$condition=="N"] <- "Goatham aisle"
sample_data(test2)$condition[sample_data(test2)$condition=="Y"] <- "Heineken tree"
sample_data(test2)$condition[sample_data(test2)$condition=="N"] <- "Heineken aisle"
sample_data(test3)$condition <- as.character(sample_data(test3)$condition)
sample_data(test4)$condition <- as.character(sample_data(test4)$condition)
sample_data(test3)$condition[sample_data(test3)$condition=="Y"] <- "Goatham tree"
sample_data(test3)$condition[sample_data(test3)$condition=="N"] <- "Goatham aisle"
sample_data(test4)$condition[sample_data(test4)$condition=="Y"] <- "Heineken tree"
sample_data(test4)$condition[sample_data(test4)$condition=="N"] <- "Heineken aisle"

mytransbiom  <- merge_phyloseq(test1,test2)
mytransbiom2 <- merge_phyloseq(test3,test4)


A <- plotTaxa(mytransbiom,"class","condition",type=1, others=T,trans=F,ordered=F,ncol=2,fixed=T)+theme(text=element_text(size=16))
B <- plotTaxa(mytransbiom2,"class","condition",type=1, others=T,trans=F,ordered=F,ncol=2,fixed=T)+theme(text=element_text(size=16))
mylist <- list(A,B)
ml <- mylist
for (i in 1:length(ml)) {
ml[[i]] <- ml[[i]] + geom_text(aes(label = LETTERS[i], x = 8, y = 11), hjust = -1, size=7)+theme(text = element_text(size=14))
ml[[i]] <- ggplot_gtable(ggplot_build(ml[[i]]))
ml[[i]]$layout$clip[ml[[i]]$layout$name == "panel"] <- "off"
}
pdf("test.pdf",height=10,width=9)
grid.arrange(grobs=ml,ncol=1,nrow=2)
dev.off()

mylegendA<-g_legend(A)
mylegendB<-g_legend(B)
A <- A + theme(legend.position="none")
B <- B + theme(legend.position="none")
ml<-list(A,B)
for (i in 1:length(ml)) {
ml[[i]] <- ml[[i]] + geom_text(aes(label = LETTERS[i], x = 5, y = 100), hjust = -1, size=5)
ml[[i]] <- ggplot_gtable(ggplot_build(ml[[i]]))
ml[[i]]$layout$clip[ml[[i]]$layout$name == "panel"] <- "off"
}
grid.arrange(grobs=ml)
dev.off()
ml <- list(mylegendA,mylegendB)
for (i in 1:length(ml)) {
ml[[i]] <- ml[[i]] + geom_text(aes(label = LETTERS[i], x = 5, y = 100), hjust = -1, size=5)
ml[[i]] <- ggplot_gtable(ggplot_build(ml[[i]]))
ml[[i]]$layout$clip[ml[[i]]$layout$name == "panel"] <- "off"
}
grid.arrange(grobs=ml)
dev.off()

mx <- list(g4,g3,g2,g1)
for (i in 1:length(mx)) {
  mx[[i]] <- mx[[i]] + geom_text(aes(label = LETTERS[i], x = 8, y = 11), hjust = -1, size=5)+theme(text = element_text(size=14))
  mx[[i]] <- ggplot_gtable(ggplot_build(mx[[i]]))
  mx[[i]]$layout$clip[mx[[i]]$layout$name == "panel"] <- "off"
}
mylegend<-g_legend(g5)
mx[[5]]<-mylegend
pdf("dge_reordered.pdf",height=10,width=9)
# grid.arrange(g1,g2,g3,g4,mylegend,ncol=2,nrow=3,layout_matrix = rbind(c(1,2),c(3,4), c(5,5)),widths = c(2.7, 2.7), heights = c(2.5,2.5, 0.2))
grid.arrange(grobs=mx,ncol=2,nrow=3,layout_matrix = rbind(c(1,2),c(3,4), c(5,5)),widths = c(2.7, 2.7), heights = c(2.5,2.5, 0.2))
dev.off()
dev.off()





get_legend <- 
function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g_legend
function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
}
