
library(GSEABase)
library(GSVA)
library(pheatmap)

setwd("../res_various/res_0.3/GSVA_flt_res0.3/")

exprSet <- readRDS('tmp/exprSet.Rds'); nrow(exprs(exprSet)) # 

### Get gene sets ###
geneSet <- read.delim("../../../../../../../__geneset__/h.all.v7.2.symbols.gmt", header = F, row.names = 1); nrow(geneSet)
rownames(geneSet) <- unlist(lapply(rownames(geneSet), function(x) unlist(strsplit(as.character(x), split = 'HALLMARK_'))[2] ))
geneSet[1:4, 1:4]
head(geneSet$V2)

gsInfo <- data.frame(row.names = rownames(geneSet), Source = geneSet$V2); head(gsInfo)
geneSet$V2 <- NULL
geneSet <- data.frame(t(geneSet), check.names = F, check.rows = F)
geneSet[1:4, 1:4]
head(geneSet)

gsInfo <- AnnotatedDataFrame(data=gsInfo)
identical(rownames(pData(gsInfo)), colnames(geneSet))

signatureSet <- ExpressionSet(as.matrix(geneSet), phenoData=gsInfo)
exprs(signatureSet)

gs_list <- list()
for (i in c(1:ncol(exprs(signatureSet)))) {
  gs_list[[i]] <- GeneSet(exprs(signatureSet)[, i][exprs(signatureSet)[, i] != ''], setName = as.character(colnames(exprs(signatureSet))[i]))
}

mySignatures <- GeneSetCollection(object = gs_list); mySignatures #
remove(gs_list)

### GSVA ###
head(exprs(exprSet))
head(pData(exprSet))

es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures, method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz = 2)
saveRDS(es.gsva, 'tmp/hallmarks_gsva.Rds')

gsvaScores <- as.data.frame(t(exprs(es.gsva)))
colnames(gsvaScores) <- unlist(lapply(colnames(gsvaScores), function(x) gsub('_', ' ', paste0(c(substring(as.character(x), 1, 1), tolower(substring(as.character(x), 2))), collapse = '')) ))
gsvaScores[1:4, 1:4]
write.table(gsvaScores, 'gsva/hallmarks_gsvaScores.txt', quote = F, col.names = NA, sep = '\t')

min(gsvaScores)
max(gsvaScores)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(t(gsvaScores), 
         filename = 'gsva/hallmarks_gsvaScores.pdf',
         cluster_rows = TRUE, cluster_cols = TRUE, scale = 'none', treeheight_row = 10, treeheight_col = 10,
         cellwidth = 8, cellheight = 8, border_color = NA, fontsize_row = 7, fontsize_col = 7, 
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'complete', clustering_callback = callback,
         breaks = c(-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(n = 16))
