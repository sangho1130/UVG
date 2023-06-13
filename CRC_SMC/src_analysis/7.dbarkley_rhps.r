
library(Seurat)
library(plyr)

setwd('../res_various/res_0.3/')

rhps <- read.delim('../../../../../../__geneset__/DBarkley_RHPs.txt', check.names = F)
head(rhps)

subdataset <- readRDS('../../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data); nrow(subdataset@meta.data) # 17334
levels(Idents(subdataset))

label <- subdataset@meta.data[, c(3,4,7)]
colnames(label)[3] <- "res_0.3"
head(label, n=3); nrow(label) # 

rhp_scores <- data.frame(matrix(ncol = length(levels(label$res_0.3)), nrow = ncol(rhps)))
colnames(rhp_scores) <- levels(label$res_0.3)
rownames(rhp_scores) <- colnames(rhps)
rhp_scores

for (clst in colnames(rhp_scores)) {
  tmpObj <- subset(subdataset, idents = clst)
  tmpExpr <- as.matrix(GetAssayData(tmpObj, slot = 'data'))
  
  tmpExprMean <- data.frame(row.names = rownames(tmpExpr),
                            Symbol = unlist(lapply(rownames(tmpExpr), function (x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] )), 
                            cluster = rowMeans(tmpExpr))
  head(tmpExprMean)
  
  for (rhp in colnames(rhps)) {
    rhp_genes <- as.character(rhps[, rhp][rhps[, rhp] != ""])
    tmpExprMean_rhp_genes <- subset(tmpExprMean, Symbol %in% rhp_genes)
    head(tmpExprMean_rhp_genes)
    rhp_scores[rhp, clst] <- mean(tmpExprMean_rhp_genes$cluster)
    
  }
}

head(rhp_scores, n=3)

library(pheatmap)
library(RColorBrewer)

pheatmap(rhp_scores, scale = 'row', filename = 'dbarkley_rhps.pdf',
         treeheight_row = 10, treeheight_col = 10,
         cellwidth = 7, cellheight = 7, fontsize_row = 9, fontsize_col = 9, border_color = NA, 
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete',
         breaks = c(-25:25)/10,
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(n = 51))
