library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

source("1.newVariGenes.R")
source("2.various_res.R")

malObj <- readRDS("../tmp/crc_smc.malignantcells.Rds")
head(malObj@meta.data); nrow(malObj@meta.data) # 17347
Idents(malObj) <- "Patient"

# filter minor libraries (less than 0.1% of total cells)
trim <- names(which(summary(malObj@meta.data$Patient) < sum(nrow(malObj@meta.data))*0.001))
print (paste(c(trim, " (n=", summary(malObj@meta.data$Patient)[trim], ") ", "will be removed"), collapse = ""))

malObj <- subset(malObj, idents = setdiff(levels(Idents(malObj)), trim))
malObj@meta.data <- droplevels(malObj@meta.data)
head(malObj@meta.data); nrow(malObj@meta.data) # 17334


### 1. Variance normalization ###
malObj <- newVariGenes(malObj, "./", "Patient", c("Patient", "nCount_RNA"))
pvals <- data.frame(malObj@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, "PC"][1]-1


### 2. Clustering ###
malObj <- various_res(malObj, pcNum = pcs_use, scaleFactor = 10, rangeMax = 2)
saveRDS(malObj, "./tmp/crc_smc.malignantcells.Rds")