## Benchmarking of different models and approaches to analyze GBS-MeDIP data, so far the two models that are more promissing for cluster 2 and 3 are a non-parametric method: Mann-whitney and a parametric: the moderated t test
#Let us start by setting the working directory and loading the packages:
library(statmod)
library("edgeR")
library(tidyverse)
library(tibble)
library(stringr)
library(palmerpenguins)
library(data.table)
library(dplyr)

# Mann-whitney ####
load("C:/Users/your_path/counts.final.cpm.tmm.normalized.RBC.rda")
load("C:/Users/your_path/design.list")
group.list.RBC=as.data.frame(design.list$blindtreatment)
colnames(group.list.RBC)=c("group")
list.mann.whitney.RBC=list()
for (i in 1:nrow(count.cpm.tmm.norm.rbc)) {
  counts.window.RBC=as.data.frame(t(count.cpm.tmm.norm.rbc[i,-(1:2)]))
  colnames(counts.window.RBC)=c("V1")
  out = wilcox.test(counts.window.RBC$V1~group.list.RBC$group)
  out[8]=count.cpm.tmm.norm.rbc$gene[i]
  list.mann.whitney.RBC[[i]]=out
}

# Now that you have the list extract the p-values retaining the window #####
df.RBC=as.data.frame(do.call(cbind, list.mann.whitney.RBC))
df.RBC=t(df.RBC)
df.RBC=as.data.frame(df.RBC)
p.values.RBC=df.RBC[,c(3,8)]
colnames(p.values.RBC)=c("p_value","window")
p.values.RBC[p.values.RBC=="NaN"]<-NA
p.values.mann.whitney=na.omit(p.values.RBC, cols=c("p_value"))
save(p.values.mann.whitney,file = "C:/Users/your_path/p.values.mann.whitney.rda")

# Perform the Benjamini-Hochberg correction for multiple testing #####
p.values.mann.whitney$corrected_p_value = p.adjust((unlist(p.values.mann.whitney.cluster.2$p_value)), method = "BH")
save(corrected_p.values.mann.whitney,file = "C:/Users/your_path/corrected_p.values.mann.whitney.rda")

# Moderated t test #####
load("C:/your_path/count.table.final.counts.RBC.rda")
load("C:/Users/your_path/design.list")
counts.raw.rbc=count.table.RBC
dge <- DGEList(counts=counts.raw.rbc[,-c(1:2)], genes=counts.raw.rbc[,1])
rownames(dge$counts) <- rownames(dge$genes) <- counts.raw.rbc[,1]
dge <- calcNormFactors(dge)
design=model.matrix(~design.list$blindtreatment)
logCPM <- cpm(dge, log=TRUE)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
top.p.values.t.test.limma=topTable(fit, coef=ncol(design),number = nrow(counts.raw.rbc))
save(top.p.values.t.test.limma,file = "C:/Users/your_path/t.test.limma.rbc.rda")

# Perform the Benjamini-Hochberg correction for multiple testing #####
top.p.values.t.test.limma$corrected_p_value = p.adjust((unlist(top.p.values.t.test.limma$p_value)), method = "BH")
save(corrected_top.p.values.t.test.limma,file = "C:/Users/your_path/corrected_top.p.values.t.test.limma.rda")

