#差异基因分析 主要使用的R包 包括limma+glimma+edgeR
library(limma)
library(Glimma)
library(edgeR)
library(RColorBrewer)
#read in data
data <- read.table('Expression_Profiling/mRNAExpressionProfiling.txt',sep="\t",encoding = 'UTF-8',header=1,fill = 0)

#数据预处理构建完整的表达矩阵
data_count<- data[,c(13:21)]
rownames(data_count) <- data[,2]
colnames(data_count) <- c('PLGA1','PLGA2','PLGA3','Blank1','Blank2','Blank3','PRIP1','PRIP2','PRIP3')
#组织样品信息
data_matrix <- data_count[,c(4,5,6,7,8,9)]
groups <- as.factor(c(rep('Blank',3),rep('PRIP',3)))
x <- DGEList(counts=data_matrix,group = groups)
#过滤低表达的基因，进行标准化操作
keep <- rowSums(cpm(x) > 1 ) >= 2 
x <- x[keep, ,keep.lib.sizes = FALSE]
x_norm <- calcNormFactors(x, method = 'TMM')  #TMM标准化
design <- model.matrix(~groups)    #构建分组矩阵
dge <- estimateDisp(x_norm, design, robust = TRUE) #估算离散值
#差异分析
fit <- glmFit(dge, design, robust = TRUE)     #拟合模型
lrt <- glmLRT(fit)   #统计检验
topTags(lrt)
#write.csv(topTags(lrt, n = nrow(x$counts)), 'plgavsblank差异基因.csv', quote = FALSE) #输出主要结果
dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr',p.value = 0.05)  #查看默认方法获得的差异基因
summary(dge_de)
#计算基因的上调和下调

PRIP_dw <- c()
for (i in 1:length(lrt$table$logFC))
{
  if (lrt$table$logFC[i] <= -2){PRIP_dw = c(PRIP_dw,rownames(lrt$table[i,]))}
}

PRIP_up <- c()
for (i in 1:length(lrt$table$logFC))
{
  if (lrt$table$logFC[i] >= 2){PRIP_up = c(PRIP_up,rownames(lrt$table[i,]))}
}
library(gplots)
feature <- c(PRIP_up,PRIP_dw)
mycol <- colorpanel(1000,"blue","white","red")
degplot <- heatmap.2(x_norm$counts[feature,], scale="row",
                     labRow=FALSE, labCol=c('Blank1','Blank2','Blank3','PRIP1','PRIP2','PRIP3'), 
                     col=mycol, trace="none", density.info="none", 
                     margin=c(6,8), lhei=c(2,10), dendrogram="none",Rowv = FALSE,Colv = FALSE)
#Colv=FALSE
write.csv(lrt$table[c(PRIP_up,PRIP_dw),], 'PLGAvsblank差异表达基因.csv', quote = FALSE)