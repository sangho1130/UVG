
library(Seurat)
library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)

subdataset <- readRDS('../tmp/subdataset_newvari.Rds')
head(subdataset@meta.data, n=3); nrow(subdataset@meta.data)

pvals <- data.frame(subdataset@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, 'PC'][1]-1
pcs_use # 29

#subdataset <- RunTSNE(subdataset, dims=1:pcs_use, reduction = "pca", reduction.key='tSNE', dim.embed=2, seed.use = 23031401)
#subdataset <- RunUMAP(subdataset, dims=1:pcs_use, reduction = "pca", reduction.key='UMAP', n.components=2, min.dist=0.2, seed.use = 23031451)

getPalette <- colorRampPalette(brewer.pal(9, 'Set1'))
DimPlot(subdataset, reduction = 'tsne', pt.size = .5, label = T, label.size = 7, cols = getPalette(8))
DimPlot(subdataset, reduction = 'umap', pt.size = .5, label = T, label.size = 7, cols = getPalette(8))

saveRDS(subdataset, '../tmp/subdataset_newvari.Rds')

umap_coord <- data.frame(as.matrix(Embeddings(subdataset, reduction = 'umap')), check.names = F, check.rows = F)
umap_coord$Cluster <- subdataset@meta.data$RNA_snn_res.0.3
umap_coord$Cluster <- factor(umap_coord$Cluster, levels = as.character(c(0:7)))
head(umap_coord)

ggplot(umap_coord, aes(UMAP_1, UMAP_2, col = Cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = getPalette(8)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
ggsave('../res_various/res_0.3/umap.res_0.3.uvg.pdf', units = 'cm', width = 10, height = 8)

ggplot(umap_coord, aes(UMAP_1, UMAP_2, col = Cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = getPalette(8)) +
  theme_void() +
  theme(legend.position = 'none')
ggsave('../res_various/res_0.3/umap.res_0.3.uvg.png', units = 'cm', width = 6, height = 6)

umap_coord$Patient <- subdataset@meta.data$Patient
head(umap_coord)

ggplot(umap_coord, aes(UMAP_1, UMAP_2, col = Patient)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = getPalette(22)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
ggsave('../res_various/res_0.3/umap.Patient.uvg.pdf', units = 'cm', width = 12, height = 8)

ggplot(umap_coord, aes(UMAP_1, UMAP_2, col = Patient)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = getPalette(22)) +
  theme_void() +
  theme(legend.position = 'none')
ggsave('../res_various/res_0.3/umap.Patient.uvg.png', units = 'cm', width = 6, height = 6)

saveRDS(umap_coord, "../tmp/umap_coord.Rds")
