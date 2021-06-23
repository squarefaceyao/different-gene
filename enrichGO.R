rm(list=ls())
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(org.At.tair.db)
library(topGO)
# BiocManager::install("topGO")  #画GO图用的
exp1 <- read.csv("genes.tsv",sep="\t",header = TRUE)

DEG.gene_symbol = as.character(exp1$Gene.symbol) #获得基因 symbol ID

DEG.entrez_id = mapIds(x = org.At.tair.db,
                       keys = DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
DEG.entrez_id = na.omit(DEG.entrez_id)

erich.go.BP = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.At.tair.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

##分析完成后，作图
dotplot(erich.go.BP)

erich.go.CC = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.At.tair.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
## 画图
barplot(erich.go.CC)

plotGOgraph(erich.go.BP)
pdf(file="./enrich.go.bp.tree.pdf",width = 10,height = 15)
plotGOgraph(erich.go.BP)
dev.off()

