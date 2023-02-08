### Load the needed libraries
# Biostrings is needed for pattern identification
# BSgenome.Ggallus.UCSC.galGal4 is the chicken genome
# plyr and reshape2 are needed for manipulating the data format

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install_local("/local/path/BSgenome.Ggallus.NCBI.galGal7_1.7.tar.gz") #i guess we need to give them this genome?
library("BSgenome.Ggallus.NCBI.galGal7")
library(Biostrings)
library(Rsubread)
library(stringr)
library(genomation)
library(plyr)

setwd("/local/path/")

########### Constructing PstI digestion in silico in Gallus gallus 7 #################
# Identify the PstI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of PstI digested fragments
genome= BSgenome.Ggallus.NCBI.galGal7
mdf.Ggallus=data.frame();
for (i in seq_along(Ggallus)){
  m<-matchPattern(c("CTGCAG"), Ggallus[[i]])    #PstI
  starts<-start(gaps(m))
  ends<-end(gaps(m))
  temp_df<-data.frame(start=starts-5,end=ends+1,chr=seqnames(Ggallus)[i]) #actually end = ends
  temp_df$start<-replace(temp_df$start, temp_df$start <= 0, 1)
  mdf.Ggallus<-rbind(mdf.Ggallus,temp_df)
}
# Now export directly the bed file as gff2/gtf format
export(mdf.Ggallus, "in_silico_windows_galgal7.gff", format = "gtf")
