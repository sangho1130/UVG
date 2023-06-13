
library(Seurat)
library(reshape2)
library(ggplot2)
library(lisi)

uvgObj <- readRDS("../tmp/subdataset_newvari.Rds")
head(uvgObj@meta.data, n=3)

uvg_umap <- data.frame(as.matrix(Embeddings(uvgObj, reduction = 'umap')), check.names = F, check.rows = F)
head(uvg_umap, n=3)

uvg_lisi <- compute_lisi(uvg_umap, uvgObj@meta.data, c('Patient'))
head(uvg_lisi, n=3)
summary(uvg_lisi)

hvgObj <- readRDS("../../original_flt/tmp/mini_subdataset.Rds")
head(hvgObj@meta.data, n=3)

hvg_umap <- data.frame(as.matrix(Embeddings(hvgObj, reduction = 'umap')), check.names = F, check.rows = F)
head(hvg_umap, n=3)

hvg_lisi <- compute_lisi(hvg_umap, hvgObj@meta.data, c('Patient'))
head(hvg_lisi, n=3)
summary(hvg_lisi)

lisi_mat <- cbind(uvg_lisi, hvg_lisi)
colnames(lisi_mat) <- c("UVG", "HVG")
t.test(lisi_mat$UVG, lisi_mat$HVG) # t = 187.76, df = 33598, p-value < 2.2e-16

lisi_mat_m <- melt(lisi_mat)
head(lisi_mat_m, n=3)

ggplot(lisi_mat_m, aes(variable, value, fill = variable)) +
    geom_boxplot() +
    scale_fill_manual(values = c("HVG"="grey80", "UVG"="blue2")) +
    labs(x = "", y = "LISI") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"))
ggsave("../res_various/res_0.3/lisi.pdf", units = "cm", width = 6, height = 5)
