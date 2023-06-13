
library(reshape2)
library(ggplot2)
library(aricode)

uvg <- read.delim('../res_various/labels.txt', row.names = 1)
head(uvg, n = 3); nrow(uvg) # 17334
uvg$seurat_clusters <- NULL

hvg <- read.delim("../../original_flt/res_various/labels.txt", row.names = 1)
head(hvg, n = 3); nrow(hvg) # 17334
hvg$seurat_clusters <- NULL

arimat <- data.frame(matrix(nrow = 20, ncol = 3))
colnames(arimat) <- c("res", "hvg", "uvg")
arimat$res <- paste("res", c(1:20)/10, sep = "")

for (i in c(1:20)) {
  res <- i/10
  arimat[i, "hvg"] <- ARI(hvg[, paste(c("RNA_snn_res.", res), collapse = "")], hvg$Patient)
}
arimat

for (i in c(1:20)) {
  res <- i/10
  arimat[i, "uvg"] <- ARI(uvg[, paste(c("RNA_snn_res.", res), collapse = "")], uvg$Patient)
}
arimat

arimat_m <- melt(arimat)
head(arimat_m)

ggplot(arimat_m, aes(res, value, group = variable, col = variable)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("hvg"="grey80", "uvg"="blue2")) +
  labs(x = "", y = "ARI") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("../res_various/res_0.3/ari_v2.pdf", units = "cm", width = 12, height = 5)

t.test(arimat$uvg, arimat$hvg) # t = -15.105, df = 19.766, p-value = 2.568e-12

ggplot(arimat_m, aes(variable, value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = c("hvg"="grey80", "uvg"="blue2")) +
  labs(x = "", y = "ARI") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("../res_various/res_0.3/ari_boxplot.pdf", units = "cm", width = 6, height = 5)
