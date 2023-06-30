library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

source("../1.new_variGenes.R")
source("../2.various_res.R")

malObj <- readRDS("tmp/10x.crc.Rds")
head(malObj@meta.data); nrow(malObj@meta.data) # 17347
Idents(malObj) <- "Patient"

# filter minor libraries (less than 0.1% of total cells)
trim <- names(which(summary(malObj@meta.data$Patient) < sum(nrow(malObj@meta.data))*0.001))
if (length(trim) != 0) {
  print (paste(c(trim, " (n=", summary(malObj@meta.data$Patient)[trim], ") ", "will be removed"), collapse = ""))
  
  malObj <- subset(malObj, idents = setdiff(levels(Idents(malObj)), trim))
  malObj@meta.data <- droplevels(malObj@meta.data)
  head(malObj@meta.data); nrow(malObj@meta.data) # 17334
} else {
  print ("No cells filtered")
}


### 1. Variance normalization ###
malObj <- newVariGenes(malObj, "./", "Patient", c("Patient", "nCount_RNA"), pcNum = 40, pltf = "3p")
pvals <- data.frame(malObj@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, "PC"][1]-1; pcs_use


### 2. Clustering ###
malObj <- various_res(malObj, pcNum = pcs_use, scaleFactor = 10, rangeMax = 2)
saveRDS(malObj, "./tmp/10x.crc.UVGs.Rds")
