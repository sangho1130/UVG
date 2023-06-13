
library(Seurat)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(pheatmap)

subdataset <- readRDS('../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data, n=3); nrow(subdataset@meta.data)

subdataset@meta.data <- subdataset@meta.data[, c(1:6, 10)]
head(subdataset@meta.data, n=3); nrow(subdataset@meta.data)

Idents(subdataset) <- "RNA_snn_res.0.3"
summary(Idents(subdataset))

saveRDS(subdataset, '../tmp/subdataset_newvari.Rds')

summary(subdataset@meta.data$RNA_snn_res.0.3)
label <- subdataset@meta.data[, c('Patient', 'RNA_snn_res.0.3')]
head(label)

prop <- data.frame(matrix(ncol = length(levels(label$RNA_snn_res.0.3)), nrow = length(levels(label$Patient))))
colnames(prop) <- as.character(levels(label$RNA_snn_res.0.3))
rownames(prop) <- levels(label$Patient)

for (cluster in levels(label$RNA_snn_res.0.3)) {
  tmp <- subset(label, RNA_snn_res.0.3 == cluster)
  tmp_prop <- summary(tmp$Patient)/sum(summary(tmp$Patient)) *100
  prop[, cluster] <- tmp_prop
}

prop_t <- melt(data.frame(Patient = rownames(prop), prop, check.rows = F, check.names = F))
head(prop_t)

getPalette <- colorRampPalette(brewer.pal(11, 'Spectral'))
ggplot(prop_t, aes(variable, value, fill = Patient)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = getPalette(22)) +
  labs(x = 'Cluster', y = 'Proportion (%)') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
ggsave('../res_various/res_0.3/prop_res_0.3.pdf', units = 'cm', width = 12, height = 9)

cellcount <- data.frame(Size = summary(label$RNA_snn_res.0.3))
cellcount <- data.frame(res = rownames(cellcount), cellcount, check.rows = F)
cellcount$res <- factor(cellcount$res, levels = levels(label$RNA_snn_res.0.3))
ggplot(cellcount, aes(res, Size)) +
  geom_bar(stat = 'identity', position = 'stack') +
  labs(x = 'Cluster', y = 'Cell count') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
ggsave('../res_various/res_0.3/count_res_0.3.pdf', units = 'cm', width = 8, height = 9)

### correlation ###
pseudo_cluster <- data.frame(matrix(ncol = length(levels(subdataset@meta.data$RNA_snn_res.0.3)), 
                                    nrow = length(rownames(subdataset))))
colnames(pseudo_cluster) <- levels(subdataset@meta.data$RNA_snn_res.0.3)
rownames(pseudo_cluster) <- rownames(subdataset)
head(pseudo_cluster, n=3)

exprs <- data.frame(as.matrix(GetAssayData(subdataset, slot = 'data')), check.names = F, check.rows = F)
for (clst in levels(subdataset@meta.data$RNA_snn_res.0.3)) {
  bcs <- rownames(subset(subdataset@meta.data, RNA_snn_res.0.3 == clst))
  clst_expr <- rowMeans(exprs[,bcs])
  pseudo_cluster[, clst] <- clst_expr
}
head(pseudo_cluster, n=3)

cormat <- cor(pseudo_cluster, method = 'spearman')
min(cormat)
max(cormat)

getPalette <- colorRampPalette(brewer.pal(9, 'RdYlBu'))
pheatmap(cormat, color = rev(getPalette(8)), 
         filename = '../res_various/res_0.3/clustering_similarity_spearman.pdf', height = 3, width = 3.5,
         breaks = c(0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0),
         legend_breaks =  c(0.8, 0.85, 0.9, 0.95, 1.0), legend_labels = c(0.8, 0.85, 0.9, 0.95, 1.0),
         border_color = NA, treeheight_row = 15, treeheight_col = 15)

### Variable genes
lib_data <- data.frame(matrix(nrow = length(VariableFeatures(subdataset)),
                              ncol = length(levels(subdataset@meta.data$Patient))))
colnames(lib_data) <- levels(subdataset@meta.data$Patient)
rownames(lib_data) <- VariableFeatures(subdataset)

for (lib in levels(subdataset@meta.data$Patient)) {
  tmp_bc <- subdataset@meta.data[subdataset@meta.data$Patient == lib, ]
  tmp_obj <- subset(subdataset, cells = rownames(tmp_bc))
  
  tmp_expr_norm <- as.matrix(GetAssayData(tmp_obj, slot = 'data'))
  lib_data[, lib] <- rowMeans(tmp_expr_norm[VariableFeatures(subdataset), ])
}
head(lib_data, n=3)

pheatmap(lib_data, cluster_rows = T, cluster_cols = T, scale = 'row', show_rownames = F,
         treeheight_row = 15, treeheight_col = 15, cellwidth = 10, 
         filename = '../res_various/res_0.3/4.variableFeatures.gmean.pdf')
