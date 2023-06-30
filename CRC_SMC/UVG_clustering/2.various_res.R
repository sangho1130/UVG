
# 2. clustering
various_res <- function(inputObj, outputPath, pcNum, scaleFactor, rangeMax = 2) {
  library(Seurat)
  library(ggplot2)
  library(plyr)
  
  dir.create(paste0(c(outputPath, "res_various"), collapse = "/"))
  dir.create(paste0(c(outputPath, "res_various/markers"), collapse = "/"))
  
  inputObj <- FindNeighbors(inputObj, dims = 1:pcNum)
  
  endpoint <- rangeMax*scaleFactor
  for (res in c(1:endpoint)){
    res <- res/scaleFactor
    
    inputObj <- FindClusters(inputObj, resolution = res)
    if (length(levels(Idents(inputObj))) != 1) {
      markers <- FindAllMarkers(object = inputObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
      markers$Symbol <- unlist(lapply(markers$gene, function(x) unlist(strsplit(x, split = '-'))[1]))
      write.table(markers, paste(c(outputPath, "/res_various/markers/markers_", pcNum, "pcs_.__res_", res, "__.txt"), collapse = ""), sep = "\t", quote = F, col.names = T, row.names = F) 
    }
  }
  
  inputObj@meta.data$orig.ident <- NULL
  inputObj@meta.data$seurat_cluster <- NULL
  label_w <- data.frame(Barcode = rownames(inputObj@meta.data), inputObj@meta.data, check.rows = F, check.names = F)
  write.table(label_w, paste(c(outputPath, "/res_various/labels.txt"), collapse = ""), sep = "\t", quote = F, col.names = T, row.names = F)
  return(inputObj)
}
