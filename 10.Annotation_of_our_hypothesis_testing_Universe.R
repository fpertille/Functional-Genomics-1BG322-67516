#############################################################################################
#### Annotation of our hypothesis testing Universe####
#############################################################################################


#install_local("//proj/geronimo_2023/private/ForgeBSGenome/BSgenome.Ggallus.NCBI.galGal7_1.7.tar.gz")
#this is a forged genome chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
library(devtools)
library(readr)
library(RMariaDB)
library(GenomicFeatures)
library(ChIPseeker)
library(ChIPpeakAnno)
library(rtracklayer)
library("GenomicRanges")
library('seqinr')
library("rGADEM")
library("Biostrings")
library("org.Gg.eg.db")
library("GenomicRanges")
library("BSgenome.Ggallus.NCBI.galGal7")
library(seqinr)
library("rGADEM")
library("Biostrings")
library(IRanges)
library(seqinr)
library(calibrate)
require(maptools)
library(calibrate)
require(maptools)
library ("ChIPseeker")
library ("org.Gg.eg.db")
library(AnnotationHub)
library("AnnotationDbi")
require(ensembldb)
library ("ChIPseeker")
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
#this is the pathway where you will find the branchfiles
#setwd("/proj/geronimo_2023/private/branchfile")

#Extracting the Universe
genome= BSgenome.Ggallus.NCBI.galGal7
mdf.Ggallus=data.frame();
for (i in seq_along(Ggallus)){
  m<-matchPattern(c("CTGCAG"), Ggallus[[i]])    #PstI
  starts<-start(gaps(m))
  ends<-end(gaps(m))
  temp_df<-data.frame(start=starts-5,end=ends+1,chr=seqnames(Ggallus)[i]) #actually end = ends
  temp_df$start<-replace(temp_df$start, temp_df$start <= 0, 1)
  # Subtract 1 from the end value for the last row of each chromosome
  if (i == length(Ggallus) || seqnames(Ggallus)[i] != seqnames(Ggallus)[i+1]) {
    last_row <- nrow(temp_df)
    temp_df[last_row, "end"] <- temp_df[last_row, "end"] - 1
  }
  mdf.Ggallus<-rbind(mdf.Ggallus,temp_df)
}

setwd("C:/your/folder")
#export(mdf.Ggallus, "in_silico_windows_galgal7.gff", format = "gtf")
#save(mdf.Ggallus, file="mdf.Ggallus.rda")

#library(data.table)
#count.table.RBC = fread(file= paste0(getwd(),"","/count_matrix/counts.matrix.RBC.galgal7.txt"))
#windowsToAnnotate <-count.table.RBC
windowsToAnnotate <-mdf.Ggallus

metadata <- makeGRangesFromDataFrame(windowsToAnnotate, start.field = "Start", end.field = "End", 
                                    seqnames.field = "Chr", strand.field = "strand", 
                                    keep.extra.columns = TRUE)
#Adding Location
mcols(metadata)$Location <- paste0(mdf.Ggallus$chr, ":", mdf.Ggallus$start, "-", mdf.Ggallus$end)
#sorting levels
metadata <- sortSeqlevels(metadata)

#In case you want to filter per library size:
#metadata <- metadata[width(metadata) > 10]
#metadata <- metadata[width(metadata) < 1000]

#Getting the sequencing of the ranges
CG  <- with(metadata, GRanges(seqnames, IRanges(start, end)))
seq_CG = BSgenome::getSeq(BSgenome.Ggallus.NCBI.galGal7, CG)
seq<-seq_CG
#Counting the number of Patterns CG per fragment
seq_CG <- vcountPattern("CG", seq_CG)
library(data.table)
####saving the seq information as character
seq <- as.data.frame(seq)
seq <- as.character(seq$x)

#Adding these columns to the metadata
mcols(metadata)$seq <- seq
mcols(metadata)$seq_CG <- seq_CG

###############################################################################
#####################GENE ANNOTATION#####################################
#downdolad from https://ftp.ensembl.org/pub/release-108/gff3/gallus_gallus/
txdb <- makeTxDbFromGFF(file="?/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.gff3.gz")
columns(txdb)
head(keys(txdb, "GENEID"))

#seqlevelsStyle(txdb) <- "NCBI"
#check that we are working with the same sequence levels
GenomeInfoDb::seqlevels(metadata)
GenomeInfoDb::seqlevels(txdb)
#convert to the wright seqlevel
#check again
#seqlevelsStyle(metadata) <- "NCBI"
#manually convert the remain chromossomes: 39, Z,  W and MT
metadata <- renameSeqlevels(metadata, c("chrW" = "W"))
metadata <- renameSeqlevels(metadata, c("chrZ" = "Z"))
metadata <- renameSeqlevels(metadata, c("chr39" = "39"))
metadata <- renameSeqlevels(metadata, c("chrMT" = "MT"))

peakAnno <- annotatePeak(metadata, tssRegion=c(-3000, 3000),
                         TxDb=txdb)
plotAnnoBar(peakAnno)
#tiff("peakAnnoPie.tiff", units="cm", width=30.00, height=30.00, res=600, compression = "lzw")
plotAnnoPie(peakAnno)
#dev.off()
info_gal <- peakAnno@anno

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
#head(listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")),"chicken")
mart <- useDataset("ggallus_gene_ensembl", mart)
ens <- info_gal$geneId
ensLookup <- gsub("\\.[0-9]*$", "", ens)
info_gal$geneId <- ensLookup 

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
  filter="ensembl_gene_id",
  values=unique(ens),
  uniqueRows=TRUE)
annotLookup <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ensLookup)],
  annotLookup)
colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_gene_id", "external_gene_name", "EntrezID", "description"))

annotLookup
info_gal <- merge(info_gal, annotLookup, by.x="geneId", by.y="ensembl_gene_id", all.x=T)
info_gal <- info_gal[!duplicated(info_gal[,c('Location')]),] #Entrez Id are duplicated

metadata_GR <- makeGRangesFromDataFrame(info_gal, start.field = "start", end.field = "end", 
                                    seqnames.field = "seqnames", strand.field = "strand", 
                                    keep.extra.columns = TRUE)
PstIgal7_metadata <- sortSeqlevels(metadata_GR)

save(PstIgal7_metadata, file="PstIgal7_metadata.rda")
write.table(PstIgal7_metadata, "PstIgal7_metadata.txt" , sep="\t" , row.names = F , quote = F)
