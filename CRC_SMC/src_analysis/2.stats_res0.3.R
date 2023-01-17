library(Seurat)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(pheatmap)


setwd('../')
subdataset <- readRDS('tmp/crc_smc.malignantcells.Rds')
head(subdataset@meta.data, n=3); nrow(subdataset@meta.data)

### Proportion of patients in each cluster
prop <- data.frame(matrix(ncol = length(levels(subdataset@meta.data$RNA_snn_res.0.3)), nrow = length(levels(subdataset@meta.data$Patient))))
colnames(prop) <- as.character(levels(subdataset@meta.data$RNA_snn_res.0.3))
rownames(prop) <- levels(subdataset@meta.data$Patient)

for (cluster in levels(subdataset@meta.data$RNA_snn_res.0.3)) {
  tmp <- subset(subdataset@meta.data, RNA_snn_res.0.3 == cluster)
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
#ggsave('res_various/res_0.3/prop_res_0.3.pdf', units = 'cm', width = 12, height = 9)


### The number of cells in each cluster
cellcount <- data.frame(Size = summary(subdataset@meta.data$RNA_snn_res.0.3))
cellcount <- data.frame(res = rownames(cellcount), cellcount, check.rows = F)
cellcount$res <- factor(cellcount$res, levels = levels(subdataset@meta.data$RNA_snn_res.0.3))

ggplot(cellcount, aes(res, Size)) +
  geom_bar(stat = 'identity', position = 'stack') +
  labs(x = 'Cluster', y = 'Cell count') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('res_various/res_0.3/count_res_0.3.pdf', units = 'cm', width = 8, height = 9)

