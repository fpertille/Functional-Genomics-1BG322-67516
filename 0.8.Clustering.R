#clusterization using the k-means method, and the elbow method to know the optimal number of clusters ####
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
setwd("C:/Users/viode560/Documents/GBSMeDIP_workpipeline")
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
count.scaled.table.RBC=scale(counts.final.cpm.tmm.normalized.RBC[,2:159])
row.names(count.scaled.table.RBC)=counts.final.cpm.tmm.normalized.RBC[,1]
# Heatmap visualizing all windows distance against all windows, every computing demanding, will not run
#distance.rbc=get_dist(count.scaled.table.RBC)
#fviz_dist(distance.rbc, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Perform the "elbow method" or sum of squares error (SSE) to have the number of clusters
wss <- (nrow(count.scaled.table.RBC)-1)*sum(apply(count.scaled.table.RBC,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(count.scaled.table.RBC,centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Perform k-means clusterization method on our set of data
k.means.cluster.metadata.RBC=kmeans(count.scaled.table.RBC,centers = 3, nstart = 100)
save(k.means.cluster.metadata.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/clusterization/k.means.cluster.metadata.RBC.rda")
graph.cluster=fviz_cluster(k.means.cluster.metadata.RBC, data = count.scaled.table.RBC)

# Now add the cluster to the original raw count matrix #####
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")
k.means.cluster.metadata.RBC=as.data.frame(k.means.cluster.metadata.RBC$cluster)
k.means.cluster.metadata.RBC=cbind(as.data.frame(rownames(k.means.cluster.metadata.RBC)),k.means.cluster.metadata.RBC)
colnames(k.means.cluster.metadata.RBC)=c("location","cluster")
count.table.raw.RBC=left_join(k.means.cluster.metadata.RBC,count.table.RBC,by="location")
save(count.table.raw.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")

# Adding a cluster columm to the normalized count matrix #####
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
counts.final.cpm.tmm.normalized.RBC=left_join(k.means.cluster.metadata.RBC,counts.final.cpm.tmm.normalized.RBC,by="location")
save(counts.final.cpm.tmm.normalized.RBC,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
