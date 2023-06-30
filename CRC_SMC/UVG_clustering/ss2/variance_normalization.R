library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

source("../1.new_variGenes.R")
source("../2.various_res.R")

malObj <- readRDS("tmp/ss2.hnscc.Rds")
head(malObj@meta.data); nrow(malObj@meta.data) # 1403
Idents(malObj) <- "Library"

# filter minor libraries (less than 0.1% of total cells)
trim <- names(which(summary(malObj@meta.data$Library) < sum(nrow(malObj@meta.data))*0.001))
if (length(trim) != 0) {
  print (paste(c(trim, " (n=", summary(malObj@meta.data$Library)[trim], ") ", "will be removed"), collapse = ""))
  
  malObj <- subset(malObj, idents = setdiff(levels(Idents(malObj)), trim))
  malObj@meta.data <- droplevels(malObj@meta.data)
  head(malObj@meta.data); nrow(malObj@meta.data)
} else {
  print ("No cells filtered")
}


### 1. Variance normalization ###
# Smart-seq2 data expects log normalized expression; e.g., log2(TPM+1) or log1p(TPM)
malObj <- newVariGenes(malObj, "./", "Library", c("Library", "nCount_RNA"), pcNum = 20, 
                       pltf = "ss2", # pltf = c("3p", "ss2") 
                       base = "log") # base = c("log", "log2")

pvals <- data.frame(malObj@reductions$pca@jackstraw$overall.p.values)
pcs_use <- pvals[pvals$Score > 0.001, "PC"][1]-1; pcs_use


### 2. Clustering ###
malObj <- various_res(malObj, pcNum = pcs_use, scaleFactor = 10, rangeMax = 2)
saveRDS(malObj, "./tmp/ss2.hnscc.UVGs.Rds")
