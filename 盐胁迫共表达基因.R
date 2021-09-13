load('GSE26983.gset.RData')

exprSet<-exprSet[,c("GSM664647","GSM664648","GSM664649","GSM664650","GSM664655",
           "GSM664656","GSM664657","GSM664658","GSM664663","GSM664664",
           "GSM664665","GSM664666")]

library(WGCNA)
if(T){
  
  exprSet[1:4,1:4]
  RNAseq_voom <- exprSet 
  ## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置
  WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:8000],])
  datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
  datExpr <- datExpr0 
}
datExpr[1:4,1:4]
if(T){
  powers = c(c(1:10), seq(from = 12, to=40, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  png("step2-beta-value.png",width = 800,height = 600)
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate, # 软yuzhi
    maxBlockSize = 30000, # 一次构建的block包含多少个基因 30000 20000 和分析的基因有关系
    TOMType = "unsigned", minModuleSize = 30, 
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE, saveTOMFileBase = 'FPKM-TOM', # 把每次计算的TOM矩阵保存起来，文件名为‘FPKM-TOM’
    loadTOMS = FALSE, verbose = 3     # 如果需要实用计算好的TOM矩阵，把loadTOMS=TRUE
  )
  table(net$colors) 
}

## step 4 ： 模块可视化
if(T){
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, mergedColors )$eigengenes
  MEs = orderMEs(MEs0)
  
  # 输出每个基因所在的模块，以及与该模块的KME值
  #file.remove('All_Gene_KME.txt')
  for(module in substring(colnames(MEs),3)){
    if(module == "grey") next
    ME=as.data.frame(MEs[,paste("ME",module,sep="")])
    colnames(ME)=module
    datModExpr=datExpr[,mergedColors==module]
    datKME = signedKME(datModExpr, ME)
    datKME=cbind(datKME,rep(module,length(datKME)))
    write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
  }
  
  #将表达矩阵转换为一个颜色矩阵，使用log10（FPKM+1）
  expColor=t(numbers2colors(log10(datExpr+1),colors=blueWhiteRed(100),naColor="grey"))
  colnames(expColor)=rownames(datExpr)
  #绘制基因的树形图，模块图，以及每个样品的表达图 ,height = 700,width = 900)
  png("wgcna.dendroColors.png")
  plotDendroAndColors(net$dendrograms[[1]], 
                      colors=cbind(mergedColors[net$blockGenes[[1]]],expColor),
                      c("Module",colnames(expColor)),
                      dendroLabels = F, hang = 0.03,
                      addGuide = T, guideHang = 0.05,
                      cex.rowText=0.5)
  dev.off()
  
  table(mergedColors)
  moduleColors=mergedColors
  # Plot the dendrogram and the module colors underneath
  png("step4-genes-modules.png",width = 800,height = 600)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.
}

