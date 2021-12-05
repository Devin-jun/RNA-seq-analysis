library(clusterProfiler)
library(org.Rn.eg.db) 
library(tidyverse)

data <- read.table('感兴趣的差异基因/下调差异基因.csv',sep=",",encoding = 'UTF-8',header=0,fill = 0)
#colnames(data) <- c('gene','logFC','logCPM','LR','PValue')

ego_MF <- enrichGO(OrgDb="org.Rn.eg.db",
                   gene = target_gene$ENTREZID,
                   pvalueCutoff = 0.05,
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)
ego_CC <- enrichGO(OrgDb="org.Rn.eg.db",
                   gene = target_gene$ENTREZID,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)
ego_BP <- enrichGO(OrgDb="org.Rn.eg.db",
                   gene = target_gene$ENTREZID,
                   pvalueCutoff = 0.01,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP))
display_number = c(length(ego_result_MF$ID), length(ego_result_CC$ID), length(ego_result_BP$ID))
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", display_number[3]), rep("cellular component", display_number[2]),
                                         rep("molecular function", display_number[1])), levels=c("molecular function", "cellular component", "biological process")))            

## numbers as data on x axis
go_enrich_df$number <- factor(1:nrow(go_enrich_df))
#names(labels) = rev(1:nrow(go_enrich_df)
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=go_enrich_df$Description) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
p
go_enrich_df2 <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                            Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                            GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                            pvalue=c(ego_result_BP$pvalue, ego_result_CC$pvalue, ego_result_MF$pvalue),
                            p.adjust=c(ego_result_BP$p.adjust, ego_result_CC$p.adjust, ego_result_MF$p.adjust),
                            geneID=c(ego_result_BP$geneID, ego_result_CC$geneID, ego_result_MF$geneID),
                            type=factor(c(rep("biological process", display_number[3]), rep("cellular component", display_number[2]),
                                          rep("molecular function", display_number[1])), levels=c("molecular function", "cellular component", "biological process")))

