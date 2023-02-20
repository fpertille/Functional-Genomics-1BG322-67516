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

#The first line of the code defines the variable genome as the genome of Gallus gallus in the BSgenome format, using the galGal7 version from the NCBI.
#The second line of the code initializes an empty data frame called mdf.Ggallus using the data.frame() function.
#The for loop then iterates through each element of the Ggallus object (which is not defined in this code block), using the seq_along() function to generate a sequence of indices to iterate over.
#For each element Ggallus[[i]], the matchPattern() function is used to identify occurrences of the PstI restriction enzyme cut site (CTGCAG). 
#The gaps() function is used to identify the gaps between the matches, and the start() and end() functions are used to obtain the start and end positions of each gap.
#A temporary data frame (temp_df) is then created with columns for the start and end positions of each gap, as well as the chromosome on which the gap occurs (seqnames(Ggallus)[i]). 
#The start position is shifted by 5 bases to the left (i.e., temp_df$start<-replace(temp_df$start, temp_df$start <= 0, 1)), and the end position is shifted by 1 base to the right (i.e., end=ends+1).
#Finally, the temporary data frame temp_df is appended to the mdf.Ggallus data frame using the rbind() function.
#Overall, this code block is identifying PstI restriction enzyme cut sites in the Gallus gallus genome and creating a data frame (mdf.Ggallus) to store the start and end positions of each cut site on each chromosome.

# Now export directly the bed file as gff2/gtf format
export(mdf.Ggallus, "in_silico_windows_galgal7.gff", format = "gtf")

#This code exports the data frame mdf.Ggallus to a GFF (Generic Feature Format) file called "in_silico_windows_galgal7.gff" using the export() function from the rtracklayer package in R.
#The first argument to export() is the data frame to be exported (mdf.Ggallus), and the second argument is the name of the output file ("in_silico_windows_galgal7.gff"). The format argument is specified as "gtf", which is a format that is similar to GFF but with some minor differences.
#The resulting GTF file will contain information on each in silico window created by the previous code block, including the start and end positions of the window and the chromosome on which it occurs. 
#This file can be used in downstream analyses, such as to identify regions of the genome that are enriched for certain features or to design primers for PCR amplification of specific regions.
