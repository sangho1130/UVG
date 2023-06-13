
library(Seurat)
library(ggplot2)
source('../../../dataset_gencode34/dotplot_expression.R')

subdataset <- readRDS('../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data); nrow(subdataset@meta.data)
summary(Idents(subdataset))

epcam <- grep(rownames(subdataset), value = T, pattern = "^EPCAM-ENSG")[1]
krt19 <- grep(rownames(subdataset), value = T, pattern = "^KRT19-ENSG")[1]
krt8 <- grep(rownames(subdataset), value = T, pattern = "^KRT8-ENSG")[1]
cd45 <- grep(rownames(subdataset), value = T, pattern = "^PTPRC-ENSG")[1]

nkg7 <- grep(rownames(subdataset), value = T, pattern = "^NKG7-ENSG")[1]
cd8a <- grep(rownames(subdataset), value = T, pattern = "^CD8A-ENSG")[1]
cd3e <- grep(rownames(subdataset), value = T, pattern = "^CD3E-ENSG")[1]

ighg1 <- grep(rownames(subdataset), value = T, pattern = "^IGHG1-ENSG")[1]
cd79b <- grep(rownames(subdataset), value = T, pattern = "^CD79B-ENSG")[1]
cd79a <- grep(rownames(subdataset), value = T, pattern = "^CD79A-ENSG")[1]

dotplot_custom(subdataset, features = c(epcam, krt19, krt8, cd45, # lineage markers
                                        nkg7, cd8a, cd3e, # Cytotoxic cells
                                        ighg1, cd79a, cd79b), # B/Plasma cells
               clstlv = rev(levels(subdataset@meta.data$RNA_snn_res.0.3)))

ggsave('../res_various/res_0.3/doublets.pdf', units = 'cm', width = 8, height = 8)
