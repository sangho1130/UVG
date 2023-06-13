
library(Seurat)
library(GSEABase)

subdataset <- readRDS('../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data); nrow(subdataset@meta.data) # 17334

subdataset <- subset(subdataset, cells = rownames(subset(subdataset@meta.data, !RNA_snn_res.0.3 %in% c('7'))))
subdataset@meta.data <- droplevels(subdataset@meta.data)

subdataset@meta.data$RNA_snn_res.0.3 <- factor(subdataset@meta.data$RNA_snn_res.0.3, 
                                               levels = c('0', '1', '2', '3', '4', '5', '6'))
head(subdataset@meta.data); nrow(subdataset@meta.data) # 17282
summary(subdataset@meta.data$RNA_snn_res.0.3)

label <- subdataset@meta.data
data <- as.matrix(GetAssayData(subdataset, slot = 'data', assay = 'RNA'))
data <- data[names(rowSums(data)[rowSums(data) != 0]), ]
data[1:4, 1:4]; dim(data) # 35817 17282

dir.create('../res_various/res_0.3/GSVA_flt_res0.3')
setwd('../res_various/res_0.3/GSVA_flt_res0.3')

dir.create('tmp')
dir.create('gsva')

### Pseudo-bulk transformation ###
pseudo_data <- data.frame(row.names = rownames(data), matrix(nrow = length(rownames(data)), ncol = length(levels(label$RNA_snn_res.0.3))))
colnames(pseudo_data) <- levels(label$RNA_snn_res.0.3)
head(pseudo_data)

for (clst in levels(label$RNA_snn_res.0.3)) {
  use_bc <- rownames(subset(label, RNA_snn_res.0.3 == clst))
  pseudo_data[, clst] <- rowMeans(data[, use_bc])
}
head(pseudo_data)
remove(data)

symbols <- unlist(lapply(rownames(pseudo_data), function(x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] ))
pseudo_data <- data.frame(Symbol = symbols, pseudo_data, check.rows = F, check.names = F)
pseudo_data[1:4, 1:4]

### Isoform selection ###
duplicated <- data.frame(table(symbols)); remove(symbols)
duplicated <- duplicated[duplicated$Freq > 1, ]
head(duplicated)

selected_iso <- c()

for (dup in duplicated$symbols) {
  tmp <- pseudo_data[pseudo_data$Symbol == dup, ]

  summ <- sort(rowSums(tmp[, -1]), decreasing = T)
  selected <- names(summ[summ == max(summ)])
  if (length(selected) == 1) {
    selected_iso <- append(selected_iso, selected)
  }
  remove(tmp)
}
length(selected_iso) # 

iso_pseudo_data <- pseudo_data[selected_iso, ]
iso_pseudo_data[1:4, 1:4]; nrow(iso_pseudo_data)
pseudo_data <- subset(pseudo_data, Symbol %in% setdiff(pseudo_data$Symbol, duplicated$symbols) )
pseudo_data <- rbind(pseudo_data, iso_pseudo_data); remove(iso_pseudo_data)
head(pseudo_data); nrow(pseudo_data)
remove(subdataset)

rownames(pseudo_data) <- pseudo_data$Symbol
pseudo_data$Symbol <- NULL

summary(pseudo_data)
nrow(pseudo_data) # 

pData <- AnnotatedDataFrame(data.frame(row.names = colnames(pseudo_data), 
                                       clustering = c('0', '1', '2', '3', '4', '5', '6')))
exprSet <- ExpressionSet(as.matrix(pseudo_data), phenoData = pData)
saveRDS(exprSet, 'tmp/exprSet.full.Rds')

pseudo_data$Sums <- rowSums(pseudo_data)
use_gene <- rownames(subset(pseudo_data, Sums > 0.01)); length(use_gene)
pseudo_data$Sums <- NULL

for (sample in colnames(pseudo_data)) {
  add_gene <- setdiff(rownames(pseudo_data[pseudo_data[, sample] > 0.01, ]), use_gene)
  use_gene <- c(use_gene, add_gene)
}
pseudo_data <- pseudo_data[use_gene, ]
colnames(pseudo_data) <- unlist(lapply(colnames(pseudo_data), function(x) paste0(c('res0.3_', as.character(x)), collapse = '') ))
head(pseudo_data); nrow(pseudo_data) # 

pData <- AnnotatedDataFrame(data.frame(row.names = colnames(pseudo_data), 
                                       clustering = c('0', '1', '2', '3', '4', '5', '6')))
exprSet <- ExpressionSet(as.matrix(pseudo_data), phenoData = pData)
saveRDS(exprSet, 'tmp/exprSet.Rds')
