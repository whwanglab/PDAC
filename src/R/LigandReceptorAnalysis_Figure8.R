library(pheatmap)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(stringr)

### LIGAND-RECEPTOR ANALYSIS
new_dfs <- readRDS('dfs_detrendApproach21-6-2.RDS')
rownames(new_dfs[[5]]) <- new_dfs[[5]]$TargetName
dfs <- readRDS('dfs_Q321-6.RDS')
ss_match2 <- readRDS(file = 'ROI_level_dataFrame_7-20.RDS')            # saved for convienence from previous analysis scripts

# L/R analysis
segs <- unique(dfs[[3]]$segment)
seg_pairs <- unique(combn(c(segs, segs), 2, paste0, collapse = '_'))

lr_pairs <- read.delim('human_lr_pair.txt', sep ='\t', header = TRUE, as.is = TRUE)
head(lr_pairs)

lr_pairs$found <- lr_pairs$ligand_gene_symbol %in% new_dfs[[2]]$Gene &
  lr_pairs$receptor_gene_symbol %in% new_dfs[[2]]$Gene
table(lr_pairs$found)

lr_pairs <- subset(lr_pairs, found)
lr_pairs[, seg_pairs] <- NA
rownames(lr_pairs) <- lr_pairs$lr_pair
# calculate correlations (spearman) for each gene in a L/R Pair across all ROIs (functionalize & shift to apply later)
for(pair in seg_pairs) {
  brk_pair <- str_split(pair, '_')[[1]]
  brk_pair <- paste0(substr(brk_pair, 1, 3), '_ROIID')
  # identify the pairs
  L_ROIs <- ss_match2[ss_match2$Treatment == 'Untreated' & ss_match2$Complete, brk_pair[1]]
  R_ROIs <- ss_match2[ss_match2$Treatment == 'Untreated' & ss_match2$Complete, brk_pair[2]]
  for(i in 1:nrow(lr_pairs)) {
    # ID genes
    L_gene <- lr_pairs$ligand_gene_symbol[i]
    R_gene <- lr_pairs$receptor_gene_symbol[i]
    # calculate correlation
    cor_val <- cor(t(new_dfs[[2]][L_gene, L_ROIs]),
                   t(new_dfs[[2]][R_gene, R_ROIs]),
                   use = 'complete.obs', method = 'spearman')
    lr_pairs[i, pair] <- cor_val
  }
  
  cat(paste0('Done with ', pair,'\n'))
}
# compare against previously written files to ensure integrity:
lr_pairsOLD <- read.csv('Detrended_LR_Analysis_Untreated.csv', row.names = 1)
all(as.matrix(round(lr_pairs[, seg_pairs], 7)) == 
      as.matrix(round(lr_pairsOLD[, seg_pairs], 7))) # allows for f.p. errors when written to CSV

## CRT Analysis
lr_pairsCRT <- lr_pairs
# calculate correlations (spearman) for each gene in a L/R Pair across all ROIs (functionalize & shift to apply later)
for(pair in seg_pairs) {
  brk_pair <- str_split(pair, '_')[[1]]
  brk_pair <- paste0(substr(brk_pair, 1, 3), '_ROIID')
  # identify the pairs
  L_ROIs <- ss_match2[ss_match2$Treatment == 'CRT' & ss_match2$Complete, brk_pair[1]]
  R_ROIs <- ss_match2[ss_match2$Treatment == 'CRT' & ss_match2$Complete, brk_pair[2]]
  for(i in 1:nrow(lr_pairsCRT)) {
    # ID genes
    L_gene <- lr_pairsCRT$ligand_gene_symbol[i]
    R_gene <- lr_pairsCRT$receptor_gene_symbol[i]
    # calculate correlation
    cor_val <- cor(t(new_dfs[[2]][L_gene, L_ROIs]),
                   t(new_dfs[[2]][R_gene, R_ROIs]),
                   use = 'complete.obs', method = 'spearman')
    lr_pairsCRT[i, pair] <- cor_val
  }
  cat(paste0('Done with ', pair,'\n'))
}
# compare against previously written files to ensure integrity:
lr_pairsCRTOLD <- read.csv('Detrended_LR_Analysis_CRT.csv', row.names = 1)
all(as.matrix(round(lr_pairsCRT[, seg_pairs], 7)) == 
      as.matrix(round(lr_pairsCRTOLD[, seg_pairs], 7))) # allows for f.p. errors when written to CSV

# explore non-self interactions & correlation values
non_self_pair <- seg_pairs[!seg_pairs %in% c('Epithelial_Epithelial', 'Immune_Immune', 'CAF_CAF')]

# Cor > .4; more values reaching .4 in treated
table(Untreated = lr_pairs[, non_self_pair] > 0.4,
      Treated = lr_pairsCRT[, non_self_pair] > 0.4)

# list non-self pairs where CRT > .4 *or* difference between CRT & Untreated > .45
lr_pairsCRT[rowSums(lr_pairsCRT[, non_self_pair] > 0.4) > 0 & 
              rowSums(lr_pairsCRT[, non_self_pair] - lr_pairs[, non_self_pair] > 0.45) > 0, ]

# Select pairs for pheatmap
inc_CRT <- rowSums(lr_pairsCRT[, non_self_pair] > 0.3) > 0
inc_names <- rownames(lr_pairs)[inc_CRT]
inc_names2 <- c()
for(i in inc_names) {
  inc_CRT <- which(lr_pairsCRT[i, non_self_pair] > 0.3)
  test_pr <- lr_pairsCRT[i, non_self_pair] - lr_pairs[i, non_self_pair]
  if(any(test_pr[inc_CRT] > 0.4)) {
    inc_names2 <- c(inc_names2, i)
  }
}

pheatmap(t(lr_pairsCRT[inc_names2, seg_pairs]),
         border_color = NA,
         labels_col = gsub('_', ' )- ', inc_names2),
         labels_row = gsub('_', ' )- ', seg_pairs),
         fontsize = 8, clustering_method = 'average',
         clustering_distance_cols = 'manhattan',
         color = colorRampPalette(c('blue3','white','white','red3'))(60),
         breaks = seq(-0.75,0.75,0.025))

# Save LR analysis results tables
write.csv(lr_pairs, 'Detrended_LR_Analysis_Untreated.csv')
write.csv(lr_pairsCRT, 'Detrended_LR_Analysis_CRT.csv')

## Graph results:
LR_plots <- list()
thr <- 0.4
dir.create('LR_Images') # in case this doesn't exist
for(i in seg_pairs[!seg_pairs %in% non_self_pair]) {
  ind <- (lr_pairsCRT[, i] > thr |
            lr_pairs[, i] > thr) &
    !(lr_pairs[, i] > (thr-.1) & lr_pairsCRT[, i] > (thr-.1))
  plt <- ggplot(subset(lr_pairsCRT, ind),
                aes_string(y = i, x = paste0('lr_pairs[ind, "',i,'"]'),
                           label = "gsub('_',')-',lr_pair)")) + 
    geom_point(data = lr_pairsCRT, aes(x = lr_pairs[, i]), color = 'gray') +
    geom_point(color = 'dodgerblue2') +
    geom_text_repel(box.padding = .25, point.padding = .15, min.segment.length = .1, size = 3, segment.color = 'dodgerblue3', segment.alpha = .5) +
    theme_bw() +
    labs(x = 'Untreated, Correlation (Rho)', y = 'CRT, Correlation (Rho)', title = gsub('_',' )- ',i)) + 
    geom_abline(slope = 1, intercept = 0)
  print(plt)
  LR_plots[[i]] <- plt
  #ggsave(paste0('LR_Images/',i,'.png'), plt, width = 10, height = 8, dpi = 300)
}
