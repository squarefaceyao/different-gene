#【KEGG富集分析】
# 1.下载/导入包并设置当前工作路径
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
BiocManager::install(c("DOSE","topGO","clusterProfiler","pathview"))

library(DOSE)
library(org.At.tair.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
rm(list=ls())
# library(KEGG.db) # 包很久没更新，没有显著富集的通路

# 2.选择差异基因
keytypes(org.At.tair.db)
MyGeneSet_table<-read.table(file.choose()) # chose GEO官方差异基因.txt
MyGeneSet<-as.character(MyGeneSet_table$V1)
typeof(MyGeneSet)

# 3.编号转换SYMBOL->ENTREZID
MyGeneIDSet=bitr(MyGeneSet,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.At.tair.db")
MyGeneIDSet

# 4.代谢通路的富集分析

ego_kegg <- enrichKEGG(MyGeneIDSet$ENTREZID, organism='ath',keyType="ncbi-geneid",
                       pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                       minGSSize=10,maxGSSize=500,use_internal_data=F)

kegg<-as.data.frame(ego_kegg@result)
write.csv(kegg,"KEGG-enrich_ath.csv",row.names = FALSE)
#（1）点图
dotplot(ego_kegg,showCategory=20)
#（2）柱状图
barplot(ego_kegg,showCategory=20,title="ego_KEGG")
#（3）差异基因关联图
emapplot(ego_kegg,showCategory = 10)
#（4）通路图
browseKEGG(ego_kegg,"ath00941")

