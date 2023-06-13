
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('../res_various/res_0.3/')
dir.create('../res_various/res_0.3/degs_mast/')

subdataset <- readRDS('../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data); nrow(subdataset@meta.data) # 17334

Idents(subdataset) <- 'RNA_snn_res.0.3'
summary(Idents(subdataset))

markers_mast <- FindAllMarkers(object = subdataset, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
nrow(markers_mast) # 2272
markers_mast <- subset(markers_mast, p_val_adj <= 0.05)
nrow(markers_mast) # 2209

markers_mast$Symbol <- unlist(lapply(markers_mast$gene, function(x) unlist(strsplit(x, split = '-ENSG'))[1]))
head(markers_mast, n=3)
write.table(markers_mast, '../res_various/res_0.3/degs_mast/markers.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 

Idents(subdataset) <- 'Patient'
summary(Idents(subdataset))
dir.create('../res_various/res_0.3/ptmarkers/')

markers_mast <- FindAllMarkers(subdataset, min.pct = 0.25, logfc.threshold = 0.25, only.pos = FALSE, test.use = 'MAST')
nrow(markers_mast) # 25567
markers_mast <- subset(markers_mast, p_val_adj <= 0.05)
nrow(markers_mast) # 25044

markers_mast$Symbol <- unlist(lapply(markers_mast$gene, function(x) unlist(strsplit(x, split = '-ENSG'))[1]))
head(markers_mast, n=3)

write.table(markers_mast, '../res_various/res_0.3/ptmarkers/markers.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 

markers_mast <- subset(markers_mast, pct.2 < 0.25)
head(markers_mast, n=3); nrow(markers_mast) # 2114
length(unique(markers_mast$gene)) # 1290
