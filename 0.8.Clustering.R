#clusterization using the k-means method, and the elbow method to know the optimal number of clusters ####
setwd("C:/Users/your_user/")
#install.packages("reshape")
#install.packages("kohonen")
#install.packages("RColorBrewer")
#install.packages("factoextra")
library(writexl)
library(plyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(ggplot2)
library(extrafont)
library(ggpubr)
library(data.table)
library(stringr)
library(dplyr)
library("edgeR")
library(reshape)
library(kohonen)
library(RColorBrewer)
library(factoextra)

load("C:/Users/your_user/counts.final.cpm.tmm.normalized.RBC.rda")
row.names(counts.final.cpm.tmm.normalized.RBC)=NULL
count.scaled.table.RBC=counts.final.cpm.tmm.normalized.RBC
count.scaled.table.RBC=scale(counts.final.cpm.tmm.normalized.RBC)

# Heatmap visualizing all windows distance against all windows
distance.rbc=get_dist(count.scaled.table.RBC)
fviz_dist(distance.rbc, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Perform the "elbow method" or sum of squares error (SSE) to have the number of clusters
wss <- (nrow(count.scaled.table.RBC)-1)*sum(apply(count.scaled.table.RBC,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(count.scaled.table.RBC,centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Perform k-means clusterization method on our set of data
kcluster.3=kmeans(count.scaled.table.RBC,centers = 3, nstart = 100)
k.means.cluster.metadata.RBC=kcluster.3
save(k.means.cluster.metadata.RBC,file = "C:/Users/your_user/k.means.cluster.metadata.RBC.rda")
graph.cluster=fviz_cluster(kcluster.3, data = count.scaled.table.RBC)

# Now add the cluster to the original raw count matrix #####
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")
count.table.RBC=cbind(as.data.frame(row.names(count.table.RBC)),count.table.RBC)
names(count.table.RBC)[names(count.table.RBC) == "row.names(count.table.RBC)"] <- 'gene'
count.table.raw.RBC=left_join(cluster.rbc,count.table.RBC,by="gene")
save(count.table.raw.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")

# Adding a cluster columm to the normalized count matrix #####
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
gene=as.data.frame(row.names(counts.final.cpm.tmm.normalized.RBC))
counts.final.cpm.tmm.normalized.RBC=cbind(as.data.frame(row.names(counts.final.cpm.tmm.normalized.RBC)),counts.final.cpm.tmm.normalized.RBC)
names(counts.final.cpm.tmm.normalized.RBC)[names(counts.final.cpm.tmm.normalized.RBC) == "row.names(counts.final.cpm.tmm.normalized.RBC))"] <- 'gene'
counts.final.cpm.tmm.normalized.RBC=left_join(cluster.rbc,counts.final.cpm.tmm.normalized.RBC,by="gene")
save(counts.final.cpm.tmm.normalized.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
