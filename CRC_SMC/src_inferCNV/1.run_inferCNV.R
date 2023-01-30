library(infercnv)

counts <- "merged.counts.txt"
label <- "merged.label.txt"
geneorder <- "~/annotation/gencodev34/gencode.v34lift37.annotation.geneorder.txt"

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = label, 
                                     delim = "\t", gene_order_file = geneorder, ref_group_names = c("CD8", "NK"))

infercnv_obj <- infercnv::run(infercnv_obj, cutoff = 0.1, out_dir = "./results/", 
                              cluster_by_groups = TRUE, denoise = TRUE, HMM = TRUE)
