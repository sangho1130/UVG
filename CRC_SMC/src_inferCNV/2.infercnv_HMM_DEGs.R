library(infercnv)
library(plyr)
library(reshape2)

label <- read.delim("../merged.label.txt", header = F, row.names = 1)
label_mal <- subset(label, !V2 %in% c("NK", "CD8"))

ptlabel <- read.delim("../../../res_various/labels.txt", row.names = 1) # an output from the clustering analysis 
ptlabel <- ptlabel[, c(1:6, 10)]; head(ptlabel)
label_mal$Patient <- ptlabel$Patient


### HMM states
# State 1: 0x: complete loss
# State 2: 0.5x: loss of one copy
# State 3: 1x: neutral
# State 4: 1.5x: addition of one copy
# State 5: 2x: addition of two copies
# State 6: 3x: essentially a placeholder for >2x copies but modeled as 3x
hmmres <- readRDS("17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj")
hmmstates <- hmmres@expr.data[, rownames(label_mal)]
remove(hmmres)
usegenes <- rownames(hmmstates)
hmmstates[1:4, 1:4]


### PT markers w/ CNV profile
ptmarkers <- read.delim("../../../res_various/res_0.6/ptmarkers/markers.MAST.txt")
head(ptmarkers, n = 3) # 3287 unique genes

inferredGenes <- intersect(usegenes, as.character(ptmarkers$gene))
length(inferredGenes) # 3027

hmmstates_markers <- as.data.frame(hmmstates[inferredGenes, ])
dim(hmmstates_markers) # 3027 17334


ord <- read.delim('../gencode.v34lift37.annotation.geneorder.txt', header = F, row.names = 1)
ord <- ord[usegenes, ]; head(ord)
hmmstates_markers <- cbind(Chr = as.character(ord[inferredGenes, "V2"]), hmmstates_markers)
hmmstates_markers[1:4, 1:4]
#saveRDS(hmmstates_markers, "hmmstates.markers.Rds")


### CNV scores 
cnv_degs <- data.frame( matrix(nrow = length(inferredGenes), ncol = length(levels(label_mal$Patient))) )
rownames(cnv_degs) <- inferredGenes
colnames(cnv_degs) <- levels(label_mal$Patient)
dim(cnv_degs)

for (pt in colnames(cnv_degs)) {
  pt_bcs <- rownames(subset(label_mal, Patient == pt))
  
  pt_degs <- as.character(subset(ptmarkers, cluster == pt)$gene)
  pt_degs <- intersect(inferredGenes, pt_degs)
  
  for (deg in pt_degs) {
    summ <- summary( as.factor(t(hmmstates_markers[deg, pt_bcs]) ))/length(pt_bcs)*100
    cnv_degs[deg, pt] <- mean(as.numeric(names(summ[summ == max(summ)])))
  }
}

cnv_degs <- cbind(Chr = hmmstates_markers$Chr, cnv_degs)
cnv_degs$Chr <- factor(cnv_degs$Chr, levels = unique(cnv_degs$Chr))
head(cnv_degs)
#saveRDS(cnv_degs, "hmmstates.markers.CNVstates.Rds")


### plot DEGs per chromosome
library(ggplot2)
dir.create("degs_hmmStates")

for (pt in colnames(cnv_degs)[-1]) {
  tmp_cnv_degs <- na.omit(cnv_degs[, c("Chr", pt)])
  tmp_cnv_degs[, pt] <- mapvalues(tmp_cnv_degs[, pt], from = c(1,2,3,4,5,6), to = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"))
  tmp_cnv_degs[, pt] <- factor(tmp_cnv_degs[, pt], levels = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"))
  
  ggplot(tmp_cnv_degs, aes(tmp_cnv_degs[, pt])) + 
    geom_histogram(stat = "count") + 
    geom_vline(xintercept = "1x") + 
    labs(x = pt, y = "DEG counts") +
    facet_wrap(~Chr, ncol = 5) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = 'black'),
          axis.text = element_text(color = 'black'))
  ggsave(paste(c('degs_hmmStates/', pt, '.chr.pdf'), collapse = ''), units = "cm", width = 14, height = 14)
  
  ggplot(tmp_cnv_degs, aes(tmp_cnv_degs[, pt])) + 
    geom_histogram(stat = "count") + 
    geom_vline(xintercept = "1x") + 
    labs(x = pt, y = "DEG counts") +
    theme_bw() + 
    theme(axis.ticks = element_line(color = 'black'),
          axis.text = element_text(color = 'black'))
  ggsave(paste(c('degs_hmmStates/', pt, '.pdf'), collapse = ''), units = "cm", width = 5, height = 5)
}


### Proportions of CNV associated DEGs
cnv_degs_summary <- data.frame( matrix(nrow = 3, ncol = length(levels(label_mal$Patient))) )
rownames(cnv_degs_summary) <- c("Gain", "Loss", "Neutral")
colnames(cnv_degs_summary) <- levels(label_mal$Patient)
cnv_degs_summary

for (pt in colnames(cnv_degs)[-1]) {
  tmp_cnv_degs <- na.omit(cnv_degs[, c("Chr", pt)])
  
  tmp_cnv_degs$cnv <- "Neutral"
  tmp_cnv_degs[rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] > 3) ), "cnv"] <- "Gain"
  tmp_cnv_degs[rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] < 3) ), "cnv"] <- "Loss"
  tmp_cnv_degs$cnv <- factor(tmp_cnv_degs$cnv, levels = c("Gain", "Loss", "Neutral"))
  
  print (pt)
  print (summary(as.factor(tmp_cnv_degs$cnv))/nrow(tmp_cnv_degs) * 100)
  cnv_degs_summary[, pt] <- summary(as.factor(tmp_cnv_degs$cnv))/nrow(tmp_cnv_degs) * 100
}
cnv_degs_summary <- cbind(CNV = rownames(cnv_degs_summary), cnv_degs_summary)
cnv_degs_summary_m <- melt(cnv_degs_summary)

ggplot(cnv_degs_summary_m, aes(variable, value, fill = CNV)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red2", "blue2", "grey80")) +
  labs(x = "", y = "Proportion (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("degs_hmmStates/Proportion.pdf", units = "cm", width = 12, height = 5)


