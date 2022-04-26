library(pheatmap)
library(scales)
library(RColorBrewer)

### UNSUPERVISED HIERARCHICAL CLUSTERING: EXTENDED DATA FIGURE 9B
new_dfs <- readRDS('dfs_detrendApproach21-6-2.RDS')
rownames(new_dfs[[5]]) <- new_dfs[[5]]$TargetName
dfs <- readRDS('dfs_Q321-6.RDS')

# Load in SnucSeq DE genes:
goi_list <- read.csv('nanostring_snuqseq_DE_genes_snuseq-Unpivot.csv')
dups <- duplicated(goi_list$Gene)
goi_list <- goi_list[!dups, ]
row.names(goi_list) <- unlist(goi_list[, 1])
goi_list$hold <- ''
goi_list <- subset(goi_list, Type != 'Epithelial (non-malignant)')

goi_list$CorNeg <- new_dfs[[5]][goi_list$Gene, 'CorrelationToNegatives']
table(goi_list$Type, goi_list$CorNeg < 0.9)

# set up cols object
cols <- list(Treatment = c(Untreated = 'white',
                           CRT = 'orange3',
                           CRTL = 'lightblue',
                           CRTLN = 'darkblue'),
             Segment = c(CAF = 'magenta2',
                         Immune = 'cyan2',
                         Epithelial = 'green3'))
cols$Type <- unique(goi_list$Type)
names(cols$Type) <- cols$Type
cols$Type['CAF'] <- 'magenta2'
cols$Type['Immune'] <- 'cyan2'
cols$Type['Epithelial (malignant)'] <- 'green3'
cols$Type[4:length(cols$Type)] <- colorRampPalette(c('darkblue','dodgerblue3','lightblue'))(12)

# Detrended pheatmap for graphing:
pheatmap(new_dfs[[2]][goi_list$Gene[goi_list$Gene %in% new_dfs[[2]]$Gene &
                                      goi_list$CorNeg < 0.90], -1],
         annotation_col = new_dfs[[3]][, c('Segment','Treatment')],
         annotation_row = goi_list[goi_list$Gene %in% new_dfs[[2]]$Gene, c('hold','Type')],
         clustering_method = 'average',
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale = 'row', show_rownames = FALSE, show_colnames = FALSE,
         annotation_colors = list(Treatment = cols$Treatment,
                                  Segment = cols$Segment,
                                  Type = cols$Type),
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         breaks = seq(-3,3,.05))

# Original Q3 pheatmap for graphing:
pheatmap(log2(dfs[[2]][goi_list$Gene[goi_list$Gene %in% dfs[[2]]$Gene &
                                      goi_list$CorNeg < 0.90], -1]),
         annotation_col = dfs[[3]][, c('Segment','Treatment')],
         annotation_row = goi_list[goi_list$Gene %in% dfs[[2]]$Gene, c('hold','Type')],
         clustering_method = 'average',
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale = 'row', show_rownames = FALSE, show_colnames = FALSE,
         annotation_colors = list(Treatment = cols$Treatment,
                                  Segment = cols$Segment,
                                  Type = cols$Type),
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         breaks = seq(-3,3,.05))
