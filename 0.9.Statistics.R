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
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/design.list")
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/counts.final.cpm.tmm.normalized.RBC.rda")
group.list.RBC=as.data.frame(design.list$blindtreatment)
colnames(group.list.RBC)=c("group")
list.mann.whitney.RBC=list()
for (i in 1:nrow(counts.final.cpm.tmm.normalized.RBC)) {
  counts.window.RBC=as.data.frame(t(counts.final.cpm.tmm.normalized.RBC[i,-(1:2)]))
  colnames(counts.window.RBC)=c("V1")
  out = wilcox.test(counts.window.RBC$V1~group.list.RBC$group)
  out[8]=counts.final.cpm.tmm.normalized.RBC$location[i]
  list.mann.whitney.RBC[[i]]=out
}

# Now that you have the list extract the p-values retaining the window #####
df.RBC=as.data.frame(do.call(cbind, list.mann.whitney.RBC))
df.RBC=t(df.RBC)
df.RBC=as.data.frame(df.RBC)
p.values.RBC=df.RBC[,c(3,8)]
colnames(p.values.RBC)=c("p_value","location")
p.values.mann.whitney=as.data.frame(t(rbindlist(lapply(p.values.RBC, as.data.table))))
colnames(p.values.mann.whitney)=c("p_value","location")
save(p.values.mann.whitney,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/p.values.mann.whitney.rda")

# Perform the Benjamini-Hochberg correction for multiple testing #####
p.values.mann.whitney$corrected_p_value = p.adjust(((p.values.mann.whitney$p_value)), method = "BH")
save(p.values.mann.whitney,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/p.values.mann.whitney.rda")

# Moderated t test #####
c=load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/design.list")
counts.raw.rbc=count.table.raw.RBC
dge <- DGEList(counts=counts.raw.rbc[,-c(1:2)], genes=counts.raw.rbc[,1])
rownames(dge$counts) <- rownames(dge$genes) <- counts.raw.rbc[,1]
dge <- calcNormFactors(dge)
design=model.matrix(~design.list$blindtreatment)
logCPM <- cpm(dge, log=TRUE)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
top.p.values.t.test.limma=topTable(fit, coef=ncol(design),number = nrow(counts.raw.rbc))
save(top.p.values.t.test.limma,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/t.test.limma.rbc.rda")

# Perform the Benjamini-Hochberg correction for multiple testing #####
top.p.values.t.test.limma$corrected_p_value = p.adjust(((top.p.values.t.test.limma$P.Value)), method = "BH")
save(top.p.values.t.test.limma,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/top.p.values.t.test.limma.rda")

load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/top.p.values.t.test.limma.rda")
load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/p.values.mann.whitney.rda")
# Create master tables ####
top.p.values.t.test.limma$location=row.names(top.p.values.t.test.limma)
top.p.values.t.test.limma=left_join(top.p.values.t.test.limma,count.table.raw.RBC, by="location")
master.table.moderated.t.test=top.p.values.t.test.limma[,c(8,9,4,7,1)]
save(master.table.moderated.t.test,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/master.table.moderated.t.test.rda")

load("C:/Users/viode560/Documents/GBSMeDIP_workpipeline/count.table.final.counts.RBC.rda")
p.values.mann.whitney=left_join(p.values.mann.whitney,count.table.raw.RBC, by="location")
master.table.mann.whitney=p.values.mann.whitney[,c(2,4,1,3)]
save(master.table.mann.whitney,file = "C:/Users/viode560/Documents/GBSMeDIP_workpipeline/master.table.mann.whitney.rda")
