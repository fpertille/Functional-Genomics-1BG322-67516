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

setwd("C:/Users/your_user/")
count.table.RBC = fread(file= paste0(getwd(),"","/count_matrix/counts.matrix.RBC.galgal7.txt"))

# We are going to filter for the maximum and minimum template length of the library, script [05.2.FragmentLengthRange.R]
count.table.RBC=count.table.RBC%>%filter(Length<<maximum_length>)
count.table.RBC=count.table.RBC%>%filter(Length><minimum_length>)
# Do some formatting of the matrix
count.table.RBC <- count.table.RBC[,2:ncol(count.table.RBC)]
colnames(count.table.RBC) <- sapply(strsplit(colnames(count.table.RBC), "/\\s*"), tail, 1)
colnames(count.table.RBC)<- gsub(pattern = ".sorted.bam", replacement = "", x = colnames(count.table.RBC))
rownames(count.table.RBC) <- paste0(count.table.RBC$Chr, ":", count.table.RBC$Start, "-", count.table.RBC$End)
count.table.RBC$Geneid=lapply(count.table.RBC$Geneid, function(x) paste('window', x))
count.table.RBC$Geneid= gsub(" ","_", count.table.RBC$Geneid)
#save as an R object for the sake of memory performance
save(count.table.RBC,file = "C:/Users/your_user/your_path/count.table.final.counts.RBC.rda")

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
save(design.list,file = "C:/Users/your_user/your_path/design.list")

# The normalization is done with EdgeR, it uses TMM normalization factor to calculate the effective library sizes (TMM*library size), then this effective library size is used to divide by the raw count.
count.table.RBC.DGElist <- DGEList(counts=count.table.RBC[,c(18:156,158:175)], genes=rownames(count.table.RBC))
rownames(count.table.RBC.DGElist$counts) <- rownames(count.table.RBC.DGElist$genes) <- rownames(count.table.RBC)
count.table.RBC.DGElist <- calcNormFactors(count.table.RBC.DGElist)
effective.library.size.RBC= count.table.RBC.DGElist$samples$lib.size*count.table.RBC.DGElist$samples$norm.factors
counts.final.cpm.tmm.normalized.RBC=data.frame(mapply(`*`,count.table.RBC[,c(18:156,158:175)],effective.library.size.RBC))
row.names(counts.final.cpm.tmm.normalized.RBC)=rownames(count.table.RBC)
save(counts.final.cpm.tmm.normalized.RBC,file = "C:/Users/your_user/your_path/counts.final.cpm.tmm.normalized.RBC.rda")

