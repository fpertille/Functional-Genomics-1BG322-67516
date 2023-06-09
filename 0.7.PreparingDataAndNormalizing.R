#In first place there will be some formatting of the data as FeatureCounts does include information that is not relevant for the whole analysis. 
#From the metadata that we extracted in the SAM files we got the template length, this is useful to filter in the count matrix that range. 
#In this project we are going to compare between two groups, stressed and non-stressed. To avoid bias the groups are coded as crocodile and elephant.
#RBC_154 has been excluded as it does not have a grouping factor.
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
library(edgeR)

setwd("C:/Users/viode560/Documents/GBSMeDIP_workpipeline")
count.table.RBC = fread(file= paste0(getwd(),"","/count_matrix/counts.matrix.RBC.galgal7.txt"))

# We are going to filter for the maximum and minimum template length of the library, script [05.2.FragmentLengthRange.R]
load(paste0(getwd(),"","/stats.rda"))
#ideally we want to load(stats) from script 05.2 and add a command to the next two line to automatically take these two values
#we take the maximum and minum from stats, stats is a summary object, is always going to be a vector where the minimum is the first string of it and the maximum is the sixth string.
count.table.RBC=count.table.RBC%>%filter(Length<as.vector(stats.RBC[6]))
count.table.RBC=count.table.RBC%>%filter(Length>as.vector(stats.RBC[1]))
# Do some formatting of the matrix
count.table.RBC <- count.table.RBC[,2:ncol(count.table.RBC)]
colnames(count.table.RBC) <- sapply(strsplit(colnames(count.table.RBC), "/\\s*"), tail, 1)
colnames(count.table.RBC)<- gsub(pattern = ".sorted.bam", replacement = "", x = colnames(count.table.RBC))
rownames(count.table.RBC) <- paste0(count.table.RBC$Chr, ":", count.table.RBC$Start, "-", count.table.RBC$End)
#save as an R object for the sake of memory performance
count.table.RBC=count.table.RBC%>%mutate(location= paste0(count.table.RBC$Chr, ":", count.table.RBC$Start, "-", count.table.RBC$End))
count.table.RBC=count.table.RBC[,c(193,17:155,157:175)]
save(count.table.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")

#Selecting individual identification and the variable columns
ind.chicken=read.csv(file= paste0(getwd(),"","/EpiAllInfo_comma.csv"),sep = ",")
ind.chicken=ind.chicken[,c(1,41)]
ind.RBC=colnames(count.table.RBC)
ind.RBC=gsub(pattern = "RBC_", replacement = "", x = ind.RBC)
ind.RBC=gsub(pattern = "^0", replacement = "", x = ind.RBC)
ind.RBC=gsub(pattern = "^0", replacement = "", x = ind.RBC)
ind.RBC=as.data.frame(ind.RBC)
colnames(ind.chicken)=c("ind","blindtreatment")
ind.chicken$ind=as.character(ind.chicken$ind)
colnames(ind.RBC)=c("ind")

design.list=left_join(ind.RBC, ind.chicken, by = "ind")
design.list[design.list==""]<-NA
design.list[design.list=="Eggs"]<-NA
design.list <- na.omit(design.list, cols=c("blindtreatment"))
#save as an R object for the sake of memory performance
save(design.list,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/design.list")

# The normalization is done with EdgeR, it uses TMM normalization factor to calculate the effective library sizes (TMM*library size), then this effective library size is used to divide by the raw count.
count.table.RBC.DGElist <- DGEList(counts=count.table.RBC[,-1], genes=count.table.RBC[,1])
rownames(count.table.RBC.DGElist$counts) <- rownames(count.table.RBC.DGElist$genes) <- rownames(count.table.RBC)
count.table.RBC.DGElist <- calcNormFactors(count.table.RBC.DGElist)
effective.library.size.RBC= count.table.RBC.DGElist$samples$lib.size*count.table.RBC.DGElist$samples$norm.factors
counts.final.cpm.tmm.normalized.RBC=data.frame(count.table.RBC[,1],mapply(`*`,count.table.RBC[,c(2:159)],effective.library.size.RBC))
save(counts.final.cpm.tmm.normalized.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
