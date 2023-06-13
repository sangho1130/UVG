
# variance normalization
newVariGenes <- function(inputObj, outputPath, libName, regress, pcNum = 60) {
  library(Seurat)
  library(ggplot2)
  library(plyr)
  
  dir.create(paste0(c(outputPath, "stats"), collapse = "/"))
  
  Idents(inputObj) <- libName 
  libIndex <- match(libName, colnames(inputObj@meta.data))
  
  
  ### 1. individual variance
  vst_vars <- data.frame(matrix(ncol = length(levels(inputObj@meta.data[,libIndex])), 
                                nrow = length(rownames(inputObj))))
  colnames(vst_vars) <- levels(inputObj@meta.data[,libIndex])
  rownames(vst_vars) <- rownames(inputObj)
  
  inputObj_spl <- SplitObject(inputObj)
  for (pt in levels(Idents(inputObj))) {
    tmpObj <- inputObj_spl[[pt]]
    tmpObj@meta.data <- droplevels(tmpObj@meta.data)
    tmpObj <- FindVariableFeatures(tmpObj, selection.method = "vst", nfeatures = 2000, clip.max = Inf)
    vst_vars[, pt] <- tmpObj@assays$RNA@meta.features$vst.variance.standardized
  }
  is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
  vst_vars[is.nan(vst_vars)] <- 0
  
  
  ### 2. pseudo-bulk variance
  # pseudo-bulk count for each sequencing library
  lib_count <- data.frame(matrix(nrow = length(rownames(inputObj)),
                                 ncol = length(levels(inputObj@meta.data[,libIndex]))))
  colnames(lib_count) <- levels(inputObj@meta.data[,libIndex])
  rownames(lib_count) <- rownames(inputObj)
  
  for (lib in levels(inputObj@meta.data[,libIndex])) {
    tmp_bc <- inputObj@meta.data[inputObj@meta.data[,libIndex] == lib, ]
    tmp_obj <- subset(inputObj, cells = rownames(tmp_bc))
    tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = "counts"))
    lib_count[, lib] <- rowSums(tmp_expr)
  }
  
  pseudo_obj <- CreateSeuratObject(counts = lib_count, project = "pseudo")
  pseudo_obj <- NormalizeData(object = pseudo_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  pseudo_obj <- FindVariableFeatures(object = pseudo_obj, selection.method = "vst", nfeatures = 2000, clip.max = Inf)
  
  
  ### 3. Set new variable genes (geometric mean)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  vst_vars$vstGmean <- unlist(lapply(c(1:nrow(vst_vars)), function(x) gm_mean(unlist(vst_vars[x, c(1:length(levels(inputObj@meta.data[,libIndex]))) ])) ))
  
  vst_vars$pseudoB <- pseudo_obj@assays$RNA@meta.features$vst.variance.standardized
  vst_vars$vstGmeanDivByPseudo <- vst_vars$vstGmean/vst_vars$pseudoB
  vst_vars <- subset(vst_vars, !is.infinite(vstGmeanDivByPseudo))
  vst_vars <- vst_vars[order(vst_vars$vstGmeanDivByPseudo, decreasing = T), ]
  new_vari_genes <- rownames(vst_vars)[1:2000]
  
  
  inputObj <- NormalizeData(object = inputObj, normalization.method = "LogNormalize", scale.factor = 10000)
  inputObj <- FindVariableFeatures(object = inputObj, selection.method = "vst", nfeatures = 2000)
  VariableFeaturePlot(inputObj)
  ggsave(paste0(c(outputPath, 'stats/1.variableFeatures.default.pdf'), collapse = '/'), units = 'cm', width = 17, height = 10)
  
  inputObj@assays$RNA@meta.features$vst.variable <- F
  inputObj@assays$RNA@meta.features[new_vari_genes, 'vst.variable'] <- T
  inputObj@assays$RNA@var.features <- new_vari_genes
  
  VariableFeaturePlot(inputObj)
  ggsave(paste0(c(outputPath, 'stats/2.variableFeatures.newHVGs.pdf'), collapse = '/'), units = 'cm', width = 17, height = 10)
  
  
  ### 4. Calculate PCs
  if ( !is.null(regress) ) {
    print (paste(c("regress: ", paste(regress, collapse = ", ")), collapse = ""))
    inputObj <- ScaleData(object = inputObj, vars.to.regress = regress, verbose = FALSE)
  } else {
    inputObj <- ScaleData(object = inputObj, verbose = FALSE)
  }
  inputObj <- RunPCA(object = inputObj, npcs = pcNum)
  inputObj <- JackStraw(object = inputObj, num.replicate = 100, dims = pcNum, verbose = FALSE)
  inputObj <- ScoreJackStraw(object = inputObj, dims = 1:pcNum)
  
  JackStrawPlot(inputObj, dims = 1:pcNum) #
  ggsave(paste0(c(outputPath, 'stats/3.jackStrawPlot.pdf'), collapse = '/'), units = 'cm', width = 25, height = 12)
  
  return(inputObj)
}
