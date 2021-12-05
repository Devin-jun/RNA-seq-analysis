library(clusterProfiler)
library(org.Rn.eg.db) 
library(tidyverse)
setwd('/Share2/home/lanxun4/ZhuJun/RNA-seq')
data <- read.table('注释结果/pripvsblank_up差异基因名.csv',sep=",",encoding = 'UTF-8',header=1,fill = 0)
genes <- as.vector(data$SYMBOL)
target_gene <- bitr(genes, 'SYMBOL', 'ENTREZID', "org.Rn.eg.db", drop = TRUE)
ego <- enrichKEGG(
  gene = target_gene$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.2,
  pAdjustMethod  = "BH",
)
ego_result <- as.data.frame(ego)
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
barplot(ego,showCategory=20)
write.csv(ego_result, 'KEGG/上调差异kegg富集.csv', quote = FALSE)
