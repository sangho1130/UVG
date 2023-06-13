
setwd("./results")

label <- read.delim("../merged.label.txt", header = F, row.names = 1)
summary(as.factor(label$V2))
label_mal <- subset(label, !V2 %in% c("NK", "CD8"))

ptlabel <- read.delim("../../../labels.txt", row.names = 1)
ptlabel <- ptlabel[, c(1:6, 10)]
head(ptlabel)
label_mal$Patient <- ptlabel$Patient

ord <- read.delim('/home/sangho/2020_newbuild/annotation/gencodev34/gencode.v34lift37.annotation.geneorder.txt', header = F, row.names = 1)

# or HMM states
# State 1: 0x: complete loss
# State 2: 0.5x: loss of one copy
# State 3: 1x: neutral
# State 4: 1.5x: addition of one copy
# State 5: 2x: addition of two copies
# State 6: 3x: essentially a placeholder for >2x copies but modeled as 3x

hmmres <- readRDS("17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj")
hmmstates <- hmmres@expr.data[, rownames(label_mal)]
hmmstates[1:4, 1:4]; dim(hmmstates) # 8372 17334

usegenes <- rownames(hmmstates)
ord <- ord[usegenes, ]; head(ord)

### PT markers w/ CNV profile
ptmarkers <- read.delim("../../../res_0.3_github/ptmarkers/markers.MAST.txt")
head(ptmarkers, n = 3); nrow(ptmarkers) # 25044

inferredGenes <- intersect(usegenes, as.character(ptmarkers$gene))
length(inferredGenes) # 4519 

hmmstates_markers <- as.data.frame(hmmstates[inferredGenes, ])
dim(hmmstates_markers) # 4519 17334
hmmstates_markers <- cbind(Chr = as.character(ord[inferredGenes, "V2"]), hmmstates_markers)
saveRDS(hmmstates_markers, "hmmstates.markers.v2.Rds")

### CNV scores 
cnv_degs <- data.frame( matrix(nrow = length(inferredGenes), ncol = length(unique(label_mal$Patient))) )
rownames(cnv_degs) <- inferredGenes
colnames(cnv_degs) <- unique(label_mal$Patient)
dim(cnv_degs)

for (pt in colnames(cnv_degs)) {
  pt_bcs <- rownames(subset(label_mal, Patient == pt))
  
  pt_degs <- subset(ptmarkers, cluster == pt); nrow(pt_degs)
  pt_degs <- subset(pt_degs, gene %in% inferredGenes); nrow(pt_degs)

  for (deg in pt_degs$gene) {
    summ <- summary( as.factor(t(hmmstates_markers[deg, pt_bcs]) ))/length(pt_bcs)*100
    cnv_degs[deg, pt] <- mean(as.numeric(names(summ[summ == max(summ)])))
  }
}

cnv_degs <- cbind(Chr = hmmstates_markers$Chr, cnv_degs)
cnv_degs$Chr <- factor(cnv_degs$Chr, levels = unique(cnv_degs$Chr))
head(cnv_degs)
saveRDS(cnv_degs, "hmmstates.markers.CNVstates.v2.Rds")

library(ggplot2)
library(plyr)

dir.create("degs_hmmStates_v2")

for (pt in colnames(cnv_degs)[-1]) {
  tmp_cnv_degs <- na.omit(cnv_degs[, c("Chr", pt)])
  
  pt_degs <- subset(ptmarkers, cluster == pt)
  pt_degs <- subset(pt_degs, gene %in% rownames(tmp_cnv_degs)); nrow(pt_degs)
  rownames(pt_degs) <- pt_degs$gene
  
  tmp_cnv_degs$cnv <- "Neutral"
  
  gains <- rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] > 3) )
  gains_up <- intersect(gains, subset(pt_degs, avg_log2FC > 0)$gene)
  tmp_cnv_degs[gains_up, "cnv"] <- "Gain (up)"
  gains_dn <- intersect(gains, subset(pt_degs, avg_log2FC < 0)$gene)
  tmp_cnv_degs[gains_dn, "cnv"] <- "Gain (dn)"
  
  losses <- rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] < 3) )
  losses_dn <- intersect(losses, subset(pt_degs, avg_log2FC < 0)$gene)
  tmp_cnv_degs[losses_dn, "cnv"] <- "Loss (dn)"
  losses_up <- intersect(losses, subset(pt_degs, avg_log2FC > 0)$gene)
  tmp_cnv_degs[losses_up, "cnv"] <- "Loss (up)"
  
  
  print (pt)
  print (summary(as.factor(tmp_cnv_degs$cnv))/nrow(tmp_cnv_degs) * 100)
  
  tmp_cnv_degs[, pt] <- round(tmp_cnv_degs[, pt])
  tmp_cnv_degs[, pt] <- mapvalues(tmp_cnv_degs[, pt], from = c(1,2,3,4,5,6), to = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"))
  tmp_cnv_degs[, pt] <- factor(tmp_cnv_degs[, pt], levels = c("0x", "0.5x", "1x", "1.5x", "2x", "3x"))

  tmp_cnv_degs$cnv <- factor(tmp_cnv_degs$cnv, 
                             levels = c("Gain (up)", "Loss (dn)", "Gain (dn)", "Loss (up)", "Neutral"))

  ggplot(tmp_cnv_degs, aes(tmp_cnv_degs[, pt], fill = cnv)) + 
    geom_histogram(stat = "count") + 
    scale_fill_manual(values = c("Gain (up)"="red2", "Loss (dn)"="blue2", "Gain (dn)"="pink", 
                                 "Loss (up)"="steelblue1", "Neutral"="grey80")) +
    geom_vline(xintercept = "1x") + 
    labs(x = pt, y = "DEG counts") +
    facet_wrap(~Chr, ncol = 5) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = 'black'),
          axis.text = element_text(color = 'black'))
  ggsave(paste(c('degs_hmmStates_v2/', pt, '.chr.pdf'), collapse = ''), units = "cm", width = 16, height = 14)
  
  ggplot(tmp_cnv_degs, aes(tmp_cnv_degs[, pt], fill = cnv)) + 
    geom_histogram(stat = "count") + 
    scale_fill_manual(values = c("Gain (up)"="red2", "Loss (dn)"="blue2", "Gain (dn)"="pink", 
                                 "Loss (up)"="steelblue1", "Neutral"="grey80")) +
    geom_vline(xintercept = "1x") + 
    labs(x = pt, y = "DEG counts") +
    theme_bw() + 
    theme(axis.ticks = element_line(color = 'black'),
          axis.text = element_text(color = 'black'))
  ggsave(paste(c('degs_hmmStates_v2/', pt, '.pdf'), collapse = ''), units = "cm", width = 8, height = 5)
}

# Proportions
library(reshape2)

cnv_degs_summary <- data.frame( matrix(nrow = 5, ncol = length(unique(label_mal$Patient))) )
rownames(cnv_degs_summary) <- c("Gain (up)", "Loss (dn)", "Gain (dn)", "Loss (up)", "Neutral")
colnames(cnv_degs_summary) <- unique(label_mal$Patient)
cnv_degs_summary[, 1:4]

for (pt in colnames(cnv_degs)[-1]) {
  tmp_cnv_degs <- na.omit(cnv_degs[, c("Chr", pt)])
  
  pt_degs <- subset(ptmarkers, cluster == pt)
  pt_degs <- subset(pt_degs, gene %in% rownames(tmp_cnv_degs)); nrow(pt_degs)
  rownames(pt_degs) <- pt_degs$gene
  
  tmp_cnv_degs$cnv <- "Neutral"
  
  gains <- rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] > 3) )
  gains_up <- intersect(gains, subset(pt_degs, avg_log2FC > 0)$gene)
  tmp_cnv_degs[gains_up, "cnv"] <- "Gain (up)"
  gains_dn <- intersect(gains, subset(pt_degs, avg_log2FC < 0)$gene)
  tmp_cnv_degs[gains_dn, "cnv"] <- "Gain (dn)"
  
  losses <- rownames( subset(tmp_cnv_degs, tmp_cnv_degs[, pt] < 3) )
  losses_dn <- intersect(losses, subset(pt_degs, avg_log2FC < 0)$gene)
  tmp_cnv_degs[losses_dn, "cnv"] <- "Loss (dn)"
  losses_up <- intersect(losses, subset(pt_degs, avg_log2FC > 0)$gene)
  tmp_cnv_degs[losses_up, "cnv"] <- "Loss (up)"
  
  print (pt)
  print (summary(as.factor(tmp_cnv_degs$cnv))/nrow(tmp_cnv_degs) * 100)
  cnv_degs_summary[names(summary(as.factor(tmp_cnv_degs$cnv))), pt] <- summary(as.factor(tmp_cnv_degs$cnv))/nrow(tmp_cnv_degs) * 100
}

cnv_degs_summary <- cbind(CNV = rownames(cnv_degs_summary), cnv_degs_summary)
cnv_degs_summary

cnv_degs_summary_w <- t(cnv_degs_summary[, c(2:ncol(cnv_degs_summary))])
cnv_degs_summary_w <- cbind(Patient = rownames(cnv_degs_summary_w), cnv_degs_summary_w)
write.table(cnv_degs_summary_w, "degs_hmmStates_v2/Proportion.txt", quote = F, sep = "\t", row.names = F, col.names = T)

cnv_degs_summary_m <- melt(cnv_degs_summary)
cnv_degs_summary_m$CNV <- factor(cnv_degs_summary_m$CNV, 
                                 levels = c("Gain (up)", "Loss (dn)", "Gain (dn)", "Loss (up)", "Neutral"))
head(cnv_degs_summary_m)

ggplot(cnv_degs_summary_m, aes(variable, value, fill = CNV)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red2", "blue2", "pink", "steelblue1", "grey80")) +
  labs(x = "", y = "Proportion (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("degs_hmmStates_v2/Proportion.pdf", units = "cm", width = 12, height = 4)


