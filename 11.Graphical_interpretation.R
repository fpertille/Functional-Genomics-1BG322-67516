###very raw

########################################################
#####################MAIN LIBRARIES####################
#######################################################
library(GenomicRanges)
library(qtl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(gghighlight)
library(tidyverse)
library(ggbio)
library(clusterProfiler)
library(data.table)
library(stringr)
library(DESeq2)
library("limma")
library(readr)
library( "gplots" )
library( "RColorBrewer" )
library (gtools)
library(org.Gg.eg.db)
#############################################################
####################LOADING NECESSARY FILES##################
#############################################################
#setwd ("/where/the/metadata/is/located")
#load(file="PstIgal7_metadata.rda")
#setwd("/where/the/statistics/are/located")
#load (file="stats.rda")
#metadata
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/branchfile/PstIgal7_metadata.rda")
#MWstats
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/branchfile/master.table.mann.whitney.rda")
#TtestStats
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/branchfile/master.table.moderated.t.test.rda")
#DMR inf.
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/t.test.limma.rbc.rda")
#Count matrix
load(file="C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/fit.rbc.rda")

########################################################
#####################PREPARING THE DATA####################
#######################################################
names(master.table.moderated.t.test)[1] <- "Location"
names(master.table.moderated.t.test)[3] <- "Pval_T"
names(master.table.moderated.t.test)[4] <- "FDR_T"
names(master.table.mann.whitney)[1] <- "Location"
names(master.table.mann.whitney)[3] <- "Pval_MW"
names(master.table.mann.whitney)[4] <- "FDR_MW"
row.names <- row.names(top.p.values.t.test.limma) 
top.p.values.t.test.limma$Location <- row.names
top.p.values.t.test.limma <- top.p.values.t.test.limma[, c("logFC", "AveExpr", "B", "Location")]
#Header
#logFC: Logarithm of the fold change in gene expression between two conditions. It indicates the magnitude and direction of the change, with positive values indicating up-regulation and negative values indicating down-regulation.
#AveExpr: Average expression level of a gene across two conditions. It represents the mean expression of a gene in two different biological states.
#B: The B-statistic is a measure of the strength of evidence for differential expression. It takes into account the variability of the gene expression measurements and the sample sizes, and is used to rank genes based on their differential expression significance.

# merge the genomic range objects based on column_name
metadata <- merge(master.table.mann.whitney, master.table.moderated.t.test, all=TRUE)
metadata2 <- merge(metadata, top.p.values.t.test.limma)
metadata <- merge(metadata2, PstIgal7_metadata, by= "Location", all.x = TRUE)

setwd(C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/)
save(metadata, file="metadata_RBCstats.rda")
save(metadata2, file="metadata_onlyStats_RBCstats.rda")

#MW_metadata <- metadata[complete.cases(metadata2$Pval_MW), ]
#T_metadata <- metadata[complete.cases(metadata2$Pval_MW), ]
#same
#metastats <- MW_metadata
#########################################################
###############LOADING THE METADATA WITH STATS############
#########################################################
##########################################################
setwd("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/")
load(file="metadata_RBCstats.rda")
#metadata$start <- as.numeric(metadata$start)
#metadata$end <- as.numeric(metadata$end)
metadata$Pval_MW <- as.numeric(metadata$Pval_MW)

#metadata2 <- metadata[complete.cases(metadata), ]

metadata_GR <- makeGRangesFromDataFrame(metadata, start.field = "start", end.field = "end", 
                                        seqnames.field = "seqnames", strand.field = "strand", 
                                        keep.extra.columns = TRUE, na.rm=TRUE)

metadata_GR <- sortSeqlevels(metadata_GR)

#metadata_GR$chrom_factor <- factor(seqnames(metadata_GR), levels = chrom_order)
#metadata_GR <- metadata_GR[order(metadata_GR$chrom_factor), ]

seqlevels(metadata_GR)

#########################################################
###############defining the output directory############
setwd("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/graphs")
getwd()
#Pval_MW    FDR_MW    Pval_T     FDR_T
#logFC    AveExpr

tiff("Vulcano_MW.tiff", units="cm", width=17.35, height=17.35, res=600, compression = "lzw")
# Make a basic volcano plot
with(metadata_GR, plot(logFC, -log10(Pval_MW), pch=20, col="white", main="MWTest", xlim=c(-3,3)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(metadata_GR, Pval_MW<=1), points(logFC, -log10(Pval_MW), pch="o", col="black"))
with(subset(metadata_GR, logFC <= -1 | logFC >= 1), points(logFC, -log10(Pval_MW), pch="o", col="red"))
#with(subset(metadata, edgeR.p.value<=1), points(edgeR.logFC, -log10(edgeR.p.value), pch="o", col="black"))
#abline(0.69897, 0, col = "red") #-log10(0.5)
#abline(0.5228787, 0, col = "yellow") #-log10(0.3)
#abline(1.30103, 0, col = "green")
abline(1.30103, 0, col = "red")#-log10(0.05)
abline(2.30103, 0, col = "yellow")#-log10(0.005)
abline(3.30103, 0, col = "green")#-log10(0.0005)
#with(subset(CE, edgeR.adj.p.value<.5), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.4, offset = 1))
#legend("bottomright", legend=levels(resbindC$sample.groups), inset=.02, horiz=F, cex=0.9, pch="o", col=c("cyan","yellow2","orangered", "green"))
dev.off() 

tiff("Vulcano_MW_FDR.tiff", units="cm", width=17.35, height=17.35, res=600, compression = "lzw")
with(metadata_GR, plot(logFC, -log10(FDR_MW), pch=20, col="white", main="MWTest", xlim=c(-3,3)))
with(subset(metadata_GR, FDR_MW<=1), points(logFC, -log10(FDR_MW), pch="o", col="black"))
with(subset(metadata_GR, logFC <= -1 | logFC >= 1), points(logFC, -log10(FDR_MW), pch="o", col="red"))
abline(0.69897, 0, col = "red") #-log10(0.5)
abline(0.5228787, 0, col = "yellow") #-log10(0.3)
abline(1.30103, 0, col = "green")
dev.off() 
#################
tiff("Vulcano_T.tiff", units="cm", width=17.35, height=17.35, res=600, compression = "lzw")
# Make a basic volcano plot
with(metadata_GR, plot(logFC, -log10(Pval_T), pch=20, col="white", main="TTest", xlim=c(-3,3)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(metadata_GR, Pval_T<=1), points(logFC, -log10(Pval_T), pch="o", col="black"))
with(subset(metadata_GR, logFC <= -1 | logFC >= 1), points(logFC, -log10(Pval_T), pch="o", col="red"))
#with(subset(metadata, edgeR.p.value<=1), points(edgeR.logFC, -log10(edgeR.p.value), pch="o", col="black"))
#abline(0.69897, 0, col = "red") #-log10(0.5)
#abline(0.5228787, 0, col = "yellow") #-log10(0.3)
#abline(1.30103, 0, col = "green")
abline(1.30103, 0, col = "red")#-log10(0.05)
abline(2.30103, 0, col = "yellow")#-log10(0.005)
abline(3.30103, 0, col = "green")#-log10(0.0005)
#with(subset(CE, edgeR.adj.p.value<.5), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.4, offset = 1))
#legend("bottomright", legend=levels(resbindC$sample.groups), inset=.02, horiz=F, cex=0.9, pch="o", col=c("cyan","yellow2","orangered", "green"))
dev.off() 

tiff("Vulcano_T_FDR.tiff", units="cm", width=17.35, height=17.35, res=600, compression = "lzw")
with(metadata_GR, plot(logFC, -log10(FDR_T), pch=20, col="white", main="MWTest", xlim=c(-3,3)))
with(subset(metadata_GR, FDR_T<=1), points(logFC, -log10(FDR_T), pch="o", col="black"))
with(subset(metadata_GR, logFC <= -1 | logFC >= 1), points(logFC, -log10(FDR_T), pch="o", col="red"))
abline(0.69897, 0, col = "red") #-log10(0.5)
abline(0.5228787, 0, col = "yellow") #-log10(0.3)
abline(1.30103, 0, col = "green")
dev.off() 

#chrorder
#library(Gviz)
#chrom_order=paste("", c(1:39, "Z", "W"), sep="")
#seqlevels(metadata_GR)
#############Manhatan
tiff("Manhattan_MW.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
plotGrandLinear(metadata_GR, coord = "genome", geom = "point", shape=21, aes(y = -log10(Pval_MW), color= seqnames), space.skip = 0.01, spaceline = T, legend = F) +
 # scale_color_manual(values=rainbow(length(levels(factor(seqnames(metadata_GR, levels=seqlevels(metadata_GR))))), start=0, end=.7)) +
  geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 
dev.off()
tiff("Manhattan_FC_MW.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
plotGrandLinear(metadata_GR, coord = "genome", geom = "point", aes(y = logFC, color = seqnames), space.skip = 0.01, spaceline = TRUE, legend = FALSE) +
  geom_point(color = ifelse(metadata_GR$Pval_MW<.005, "black", NA), shape=23, size = 3) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6))
dev.off()



#############Manhatan
tiff("Manhattan_T.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
plotGrandLinear(metadata_GR, coord = "genome", geom = "point", shape=21, aes(y = -log10(Pval_T), color= seqnames), space.skip = 0.01, spaceline = TRUE, scale=list(y=1.5), legend = F) +
  geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 
dev.off()
tiff("Manhattan_FC_T.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
plotGrandLinear(metadata_GR, coord = "genome", geom = "point", aes(y = logFC, color = seqnames), space.skip = 0.01, spaceline = TRUE, legend = FALSE) +
  geom_point(color = ifelse(metadata_GR$Pval_T<.005, "black", NA), shape=23, size = 3) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6))
dev.off()

######
#Pval_MW    FDR_MW    Pval_T     FDR_T
#logFC    AveExpr

#subseting by pvalue 0.05
metadata_GR_MW<- subset(metadata_GR, Pval_MW<=0.05)
metadata_GR_T<- subset(metadata_GR, Pval_T<=0.05)
mcols(metadata_GR_MW) <- data.frame(type = "MW", row.names = names(metadata_GR_MW))
mcols(metadata_GR_T) <- data.frame(type = "Tt", row.names = names(metadata_GR_T))
resbind <- c(metadata_GR_T, metadata_GR_MW)

#subseting by pvalue 0.005
metadata_GR_MW<- subset(metadata_GR, Pval_MW<=0.005)
metadata_GR_T<- subset(metadata_GR, Pval_T<=0.005)
mcols(metadata_GR_MW) <- data.frame(type = "MW", row.names = names(metadata_GR_MW))
mcols(metadata_GR_T) <- data.frame(type = "Tt", row.names = names(metadata_GR_T))
resbind <- c(metadata_GR_T, metadata_GR_MW)


#########Venn diagrams#########
library(ChIPpeakAnno)
library(Cairo)
CairoPDF("venn_MWvsTt_0.005.pdf")
grl <- splitAsList(resbind, resbind$type)
#res <- makeVennDiagram(Peaks=grl, NameOfPeaks=c("MW", "TT"), main="All", col = "black",fill = c("blue3","red3"), alpha = 0.50, cex= 1.5, connectedPeaks = "keepAll")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=c("MW", "TT"), main="All", col = "black",fill = c("blue3","red3"), alpha = 0.50, cex= 1.5)
dev.off()

CairoPDF("venn_MWvsTt_0.05.pdf")
grl <- splitAsList(resbind, resbind$type)
#res <- makeVennDiagram(Peaks=grl, NameOfPeaks=c("MW", "TT"), main="All", col = "black",fill = c("blue3","red3"), alpha = 0.50, cex= 1.5, connectedPeaks = "keepAll")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=c("MW", "TT"), main="All", col = "black",fill = c("blue3","red3"), alpha = 0.50, cex= 1.5)
dev.off()


OVERLAP_all <- findOverlapsOfPeaks(unique(metadata_GR_MW), unique(metadata_GR_T))

OVERLAP_all$venn_cnt
OVERLAP_all$peaklist

Unique_MW <- OVERLAP_all[["all.peaks"]]$unique.metadata_GR_MW.
Unique_T <- OVERLAP_all[["all.peaks"]]$unique.metadata_GR_T.
Overlap <- OVERLAP_all[["overlappingPeaks"]]$`unique.metadata_GR_MW.///unique.metadata_GR_T.`

####STOPED HERE#########
for_path_DMR_MW <- c()
#for_path_DMR <- na.omit(unique(metadata$EntrezID[metadata$edgeR.p.value<=0.05 & metadata$Tissue=='A']))
for_path_DMR_MW$all <- na.omit(unique(metadata$EntrezID[metadata$Pval_MW<=0.05]))
for_path_DMR_MW$gain <- na.omit(unique(metadata$EntrezID[metadata$Pval_MW<=0.05 & metadata$logFC>0]))
for_path_DMR_MW$loss <- na.omit(unique(metadata$EntrezID[metadata$Pval_MW<=0.05 & metadata$logFC<0]))

for_path_DMR_T <- c()
#for_path_DMR <- na.omit(unique(metadata$EntrezID[metadata$edgeR.p.value<=0.05 & metadata$Tissue=='A']))
for_path_DMR_T$all <- na.omit(unique(metadata$EntrezID[metadata$Pval_T<=0.05]))
for_path_DMR_T$gain <- na.omit(unique(metadata$EntrezID[metadata$Pval_T<=0.05 & metadata$logFC>0]))
for_path_DMR_T$loss <- na.omit(unique(metadata$EntrezID[metadata$Pval_T<=0.05 & metadata$logFC<0]))


tiff("GO_ANALYSIS_MW_05.tiff", units="cm", width=25.00, height=30.00, res=600, compression = "lzw")
#jpeg("GO_ANALYSIS_ALL.jpeg" , width = 1300 , height = 700)
GO_analysis <- compareCluster(geneClusters = for_path_DMR_MW , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db, pvalueCutoff = 1, pAdjustMethod = "fdr", readable=T)
dotplot(GO_analysis)
dev.off()
write.table(as.data.frame(GO_analysis@compareClusterResult) , "PathwayGO_Enrich_MW.05.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

tiff("GO_ANALYSIS_T_05.tiff", units="cm", width=25.00, height=30.00, res=600, compression = "lzw")
#jpeg("GO_ANALYSIS_ALL.jpeg" , width = 1300 , height = 700)
GO_analysis <- compareCluster(geneClusters = for_path_DMR_T , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db, pvalueCutoff = 1, pAdjustMethod = "fdr", readable=T)
dotplot(GO_analysis)
dev.off()
write.table(as.data.frame(GO_analysis@compareClusterResult) , "PathwayGO_Enrich_T.05.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

############kegg is not working
tiff("KEGG_ANALYSIS_overlap.tiff", units="cm", width=25.00, height=30.00, res=600, compression = "lzw")
#jpeg("KEGG_ANALYSIS_ALL.jpeg" , width = 1200 , height = 700)
kegg_analysis <- compareCluster(geneClusters = for_path_DMR_MW , fun = "enrichKEGG", organism="gga", pvalueCutoff = 1, qvalueCutoff=1, pAdjustMethod = "fdr")
dotplot(kegg_analysis)
dev.off()

#################################################
#################################################
#################################################
###################selecting significant windows from MEDIPS###################
setwd("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/")
load(file="metadata_RBCstats.rda")
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/count_matrix/counts.final.cpm.tmm.normalized.RBC.rda")
load("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/branchfile/design.list")

megamatrix <- merge(metadata, counts.final.cpm.tmm.normalized.RBC, by.x= c("Location"), by.y=c("location"), all.x = TRUE)

metadata_GR_MW <- subset(megamatrix, Pval_MW<=0.005)
metadata_GR_T <- subset(megamatrix, Pval_T<=0.005)

library(dplyr)
match_df <- inner_join(metadata_GR_MW, metadata_GR_T, by = "Location")
unique_df1 <- anti_join(metadata_GR_MW, metadata_GR_T, by = "Location")
unique_df2 <- anti_join(metadata_GR_T, metadata_GR_MW, by= "Location")

match_df <- filter(megamatrix, Location %in% match_df$Location)
unique_df1 <- filter(megamatrix, Location %in% unique_df1$Location)
unique_df2 <- filter(megamatrix, Location %in% unique_df2$Location)


metadata_sig<- match_df
#metadata_sig <- metadata_sig[complete.cases(metadata_sig), ]
#PREPARING DATA TO HEAD MAP - MW
rownames(metadata_sig) <- metadata_sig$Location
library(gtools)
library(DESeq2)
#ordering rows
metadata_sig<- metadata_sig[ mixedsort(row.names(metadata_sig)),]
metadata_sig <- metadata_sig [, c(31:188)] #subset check it
colnames(metadata_sig) <- design.list$ind
metadata_sig[is.na(metadata_sig)] <- 0
metadata_sig[] <- lapply(metadata_sig, as.integer)
#metadata_sig[] <- lapply(metadata_sig, as.integer)
metadata_sig <- metadata_sig[complete.cases(metadata_sig), ]
samplesheat_metadata_sig <- data.frame(row.names=colnames(metadata_sig), condition= ifelse(design.list$blindtreatment == "elephant", "A", "B"))
ddsheat <- DESeqDataSetFromMatrix(countData=metadata_sig,colData=samplesheat_metadata_sig,design=~condition)
rld <- ddsheat
sampleDist<-dist(t(assay(rld)))
sampleDistMatrix <-as.matrix(sampleDist)
rownames(sampleDistMatrix)<-rld$condition
colnames(sampleDistMatrix)<-NULL
na_cols <- which(colSums(is.na(metadata_sig)) > 0)

setwd("C:/Users/fabpe757/Box/GEroNIMOproject/GEroNIMO workshop on GBS-MeDIP data analysis/GBSMeDIP_workpipeline/outputs/graphs/")

tiff("MeDIP_HeatMap_common.tiff", units="in", width=11, height=8.5, res=600, compress="lzw")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)
graphics.off()

assay(rld)[rownames(assay(rld))%in%rownames(metadata_sig),] -> sig2heatmap

tiff("heatmap.2_common.tiff", units="in", width=11, height=8.0, res=600, compress="lzw")
heatmap.2( sig2heatmap, scale="row", Rowv = FALSE, Colv = FALSE,
           trace="none", dendrogram="none", margins = c(5, 15),    #needs 22 to ML
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
graphics.off()

tiff("heatmap.2_dendr_common.tiff", units="in", width=11, height=8.0, res=600, compress="lzw")
heatmap.2( sig2heatmap, scale="row", Rowv = T, Colv = FALSE,
           trace="none", dendrogram="both", margins = c(5, 15),    #needs 22 to ML
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
graphics.off()
#dendrogram = c("both","row","column","none")

#############################
############################
#######LAPS OF VENN DIAGRAM#
############################
for_path_DMR <- c()
for_path_DMR$A <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='A']))
for_path_DMR$C <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='C']))
for_path_DMR$H <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='H']))
for_path_DMR$AC <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='AC']))
for_path_DMR$AE <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='AE']))
for_path_DMR$AN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='AN']))
for_path_DMR$AEN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='AEN']))
for_path_DMR$CC <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='CC']))
for_path_DMR$CE <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='CE']))
for_path_DMR$CN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='CN']))
for_path_DMR$CEN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='CEN']))
for_path_DMR$HC <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='HC']))
for_path_DMR$HE <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='HE']))
for_path_DMR$HN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='HN']))
for_path_DMR$HEN <- na.omit(unique(resbind$external_gene_name[resbind$edgeR.p.value<=0.05 & resbind$Test=='HEN']))

`%notin%` <- Negate(`%in%`)

pathway_list_overlaps <- c()
pathway_list_overlaps$A <- na.omit(unique(for_path_DMR$A[for_path_DMR$A %notin% for_path_DMR$C & for_path_DMR$A %notin% for_path_DMR$H]))
pathway_list_overlaps$C <- na.omit(unique(for_path_DMR$C[for_path_DMR$C %notin% for_path_DMR$A & for_path_DMR$C %notin% for_path_DMR$H]))
pathway_list_overlaps$H <- na.omit(unique(for_path_DMR$H[for_path_DMR$H %notin% for_path_DMR$C & for_path_DMR$H %notin% for_path_DMR$A]))
pathway_list_overlaps$AiC <- na.omit(unique(for_path_DMR$A[for_path_DMR$A %in% for_path_DMR$C]))
pathway_list_overlaps$AiH <- na.omit(unique(for_path_DMR$A[for_path_DMR$A %in% for_path_DMR$H]))
pathway_list_overlaps$CiH <- na.omit(unique(for_path_DMR$C[for_path_DMR$C %in% for_path_DMR$H]))

pathway_list_overlaps$A <- na.omit(unique(for_path_DMR$A[for_path_DMR$A %notin% for_path_DMR$C & for_path_DMR$A %notin% for_path_DMR$H]))

pathway_list_overlaps <- c()
pathway_list_overlaps$AC <- na.omit(unique(for_path_DMR$AC[for_path_DMR$AC %notin% for_path_DMR$AE & for_path_DMR$AC %notin% for_path_DMR$AN & for_path_DMR$AC %notin% for_path_DMR$AEN]))
pathway_list_overlaps$AE <- na.omit(unique(for_path_DMR$AE[for_path_DMR$AE %notin% for_path_DMR$AC & for_path_DMR$AE %notin% for_path_DMR$AN & for_path_DMR$AE %notin% for_path_DMR$AEN]))
pathway_list_overlaps$AN <- na.omit(unique(for_path_DMR$AN[for_path_DMR$AN %notin% for_path_DMR$AE & for_path_DMR$AN %notin% for_path_DMR$AC & for_path_DMR$AN %notin% for_path_DMR$AEN]))
pathway_list_overlaps$AEN <- na.omit(unique(for_path_DMR$AEN[for_path_DMR$AEN %notin% for_path_DMR$AE & for_path_DMR$AEN %notin% for_path_DMR$AN & for_path_DMR$AEN %notin% for_path_DMR$AC]))

pathway_list_overlaps <- c()
pathway_list_overlaps$CC <- na.omit(unique(for_path_DMR$CC[for_path_DMR$CC %notin% for_path_DMR$CE & for_path_DMR$CC %notin% for_path_DMR$CN & for_path_DMR$CC %notin% for_path_DMR$CEN]))
pathway_list_overlaps$CE <- na.omit(unique(for_path_DMR$CE[for_path_DMR$CE %notin% for_path_DMR$CC & for_path_DMR$CE %notin% for_path_DMR$CN & for_path_DMR$CE %notin% for_path_DMR$CEN]))
pathway_list_overlaps$CN <- na.omit(unique(for_path_DMR$CN[for_path_DMR$CN %notin% for_path_DMR$CE & for_path_DMR$CN %notin% for_path_DMR$CC & for_path_DMR$CN %notin% for_path_DMR$CEN]))
pathway_list_overlaps$CEN <- na.omit(unique(for_path_DMR$CEN[for_path_DMR$CEN %notin% for_path_DMR$CE & for_path_DMR$CEN %notin% for_path_DMR$CN & for_path_DMR$CEN %notin% for_path_DMR$CC]))

pathway_list_overlaps <- c()
pathway_list_overlaps$HC <- na.omit(unique(for_path_DMR$HC[for_path_DMR$HC %notin% for_path_DMR$HE & for_path_DMR$HC %notin% for_path_DMR$HN & for_path_DMR$HC %notin% for_path_DMR$HEN]))
pathway_list_overlaps$HE <- na.omit(unique(for_path_DMR$HE[for_path_DMR$HE %notin% for_path_DMR$HC & for_path_DMR$HE %notin% for_path_DMR$HN & for_path_DMR$HE %notin% for_path_DMR$HEN]))
pathway_list_overlaps$HN <- na.omit(unique(for_path_DMR$HN[for_path_DMR$HN %notin% for_path_DMR$HE & for_path_DMR$HN %notin% for_path_DMR$HC & for_path_DMR$HN %notin% for_path_DMR$HEN]))
pathway_list_overlaps$HEN <- na.omit(unique(for_path_DMR$HEN[for_path_DMR$HEN %notin% for_path_DMR$HE & for_path_DMR$HEN %notin% for_path_DMR$HN & for_path_DMR$HEN %notin% for_path_DMR$HC]))

####overlaps_A
pathway_list_overlaps <- c()
pathway_list_overlaps$ACiAN <- na.omit(unique(for_path_DMR$AC[for_path_DMR$AC %in% for_path_DMR$AN]))
pathway_list_overlaps$ACiAE <- na.omit(unique(for_path_DMR$AC[for_path_DMR$AC %in% for_path_DMR$AE]))
pathway_list_overlaps$ACiAEN <- na.omit(unique(for_path_DMR$AC[for_path_DMR$AC %in% for_path_DMR$AEN]))
pathway_list_overlaps$ANiAE <- na.omit(unique(for_path_DMR$AN[for_path_DMR$AN %in% for_path_DMR$AE]))
pathway_list_overlaps$ANiAEN <- na.omit(unique(for_path_DMR$AN[for_path_DMR$AN %in% for_path_DMR$AEN]))
pathway_list_overlaps$AENiAE <- na.omit(unique(for_path_DMR$AEN[for_path_DMR$AEN %in% for_path_DMR$AE]))
####overlaps_C
pathway_list_overlaps <- c()
pathway_list_overlaps$CCiCN <- na.omit(unique(for_path_DMR$CC[for_path_DMR$CC %in% for_path_DMR$CN]))
pathway_list_overlaps$CCiCE <- na.omit(unique(for_path_DMR$CC[for_path_DMR$CC %in% for_path_DMR$CE]))
pathway_list_overlaps$CCiCEN <- na.omit(unique(for_path_DMR$CC[for_path_DMR$CC %in% for_path_DMR$CEN]))
pathway_list_overlaps$CNiCE <- na.omit(unique(for_path_DMR$CN[for_path_DMR$CN %in% for_path_DMR$CE]))
pathway_list_overlaps$CNiCEN <- na.omit(unique(for_path_DMR$CN[for_path_DMR$CN %in% for_path_DMR$CEN]))
pathway_list_overlaps$CENiCE <- na.omit(unique(for_path_DMR$CEN[for_path_DMR$CEN %in% for_path_DMR$CE]))
####overlaps_H
pathway_list_overlaps <- c()
pathway_list_overlaps$HCiHN <- na.omit(unique(for_path_DMR$HC[for_path_DMR$HC %in% for_path_DMR$HN]))
pathway_list_overlaps$HCiHE <- na.omit(unique(for_path_DMR$HC[for_path_DMR$HC %in% for_path_DMR$HE]))
pathway_list_overlaps$HCiHEN <- na.omit(unique(for_path_DMR$HC[for_path_DMR$HC %in% for_path_DMR$HEN]))
pathway_list_overlaps$HNiHE <- na.omit(unique(for_path_DMR$HN[for_path_DMR$HN %in% for_path_DMR$HE]))
pathway_list_overlaps$HNiHEN <- na.omit(unique(for_path_DMR$HN[for_path_DMR$HN %in% for_path_DMR$HEN]))
pathway_list_overlaps$HENiHE <- na.omit(unique(for_path_DMR$HEN[for_path_DMR$HEN %in% for_path_DMR$HE]))

library(UpSetR)
library(RColorBrewer)
tiff("VennBar_Amigdala_geneID.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
upset(data = fromList(for_path_DMR) , keep.order = T, sets= c("AC","AE", "AN", "AEN"),
      matrix.color = "blue", point.size = 5, text.scale= c(1, 1, 1, 1),
      sets.bar.color = brewer.pal(4, "Set3"))
dev.off() 

##make a comparison list
library(RVenn)
library(purrr)
library(ggplot2)
library (ggvenn)
for_path <- c()
for_path$A <- na.omit(unique(resbind$EntrezID[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='A']))
for_path$C <- na.omit(unique(resbind$EntrezID[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='C']))
for_path$H <- na.omit(unique(resbind$EntrezID[resbind$edgeR.p.value<=0.05 & resbind$Tissue=='H']))

library(VennDiagram)
set.seed(1) # For reproducibility of results
venn.diagram(for_path, filename ="Venn_Tissues_genes.tiff", height = 1000, width = 1000)
dev.off()

