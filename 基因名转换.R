library(clusterProfiler)
library(org.Rn.eg.db) 
library(tidyverse)
#基因名转换
#keytypes(org.Mm.eg.db)
#keys(org.Mm.eg.db, keytype="ENSEMBL") %>% head()
target_gene <- bitr(target_gene_id, 'SYMBOL', 'ENTREZID', "org.Rn.eg.db", drop = TRUE)