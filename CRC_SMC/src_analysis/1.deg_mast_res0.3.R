library(Seurat)
library(ggplot2)
library(plyr)

setwd('../')
dir.create('res_various/res_0.3/')
dir.create('res_various/res_0.3/degs_mast/')

subdataset <- readRDS('tmp/mini_subdataset_newvari.Rds')
head(subdataset@meta.data); nrow(subdataset@meta.data) # 17334

subdataset@meta.data <- subdataset@meta.data[, c(1:6,10)]
head(subdataset@meta.data, n=3)

subdataset$RNA_snn_res.0.3 <- factor(subdataset$RNA_snn_res.0.3, levels = c(0:7))
Idents(subdataset) <- 'RNA_snn_res.0.3'

markers_mast <- FindAllMarkers(object = subdataset, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST'); nrow(markers_mast) # 1314

markers_mast <- subset(markers_mast, p_val_adj <= 0.05); nrow(markers_mast) # 1286
markers_mast$Symbol <- unlist(lapply(markers_mast$gene, function(x) unlist(strsplit(x, split = '-ENSG'))[1]))
head(markers_mast, n=3)

write.table(markers_mast, 'res_various/res_0.3/degs_mast/markers.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(subdataset, 'tmp/crc_smc.malignantcells.Rds')
