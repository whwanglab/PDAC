setwd('~/Projects/DSP - Hwang PDAC Round 2/new_DCCs/')
load('ssGSEA_UpdatedWorkspace_5-23-21.RData')

## Load libraries
library(ggplot2)
library(ggrepel)
library(ggsci)
library(pheatmap)


## Section 1: Compare CAF / Immune / Epithelial ROIs
## Goal: ID systematic signal-bleed between matched ROIs and diagnose issue

# Example scatter of paired immune & CAF ROIs
smoothScatter(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[1]]]),
              log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$CAF_ROI[1]]]))

mod <- lm(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[1]]]) ~
           log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$CAF_ROI[1]]]))
res <- residuals(mod)

# Determine if IQRs vary by AOI Area & segment type
q25 <- apply(dfs[[1]][, -1], 2, quantile, .25, na.rm = TRUE)
q75 <- apply(dfs[[1]][, -1], 2, quantile, .90, na.rm = TRUE)
delq <- q75 - q25
hist(delq)
dfs[[3]]$IQR <- delq

# IQR histogram
ggplot(dfs[[3]], aes(x = log2(IQR), fill = segment)) +
 geom_histogram() +
 facet_wrap(~segment, nrow = 3)
# Scattered vs Area: note that area is biggest driver
ggplot(dfs[[3]], aes(x = AOI_area, y = IQR, color = segment)) +
 geom_point() +
 scale_x_continuous(trans = 'log10') +
 scale_y_continuous(trans = 'log2')

## Section 2: Exploration of residual analysis - ROI2
mod <- lm(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[2]]]) ~
           log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$CAF_ROI[2]]]))
res <- resid(mod)
smoothScatter(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[2]]]),
              res)
smoothScatter(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$CAF_ROI[2]]]),
              res)

tmp_dat <- cbind(dfs[[2]][, c('Gene', 'DSP.1001660003872.E.C04')],
                 data.frame(residuals = as.numeric(res)))
ggplot(data = tmp_dat, aes(x = log2(`DSP.1001660003872.E.C04`),
                           y = residuals,
                           label = Gene)) +
 geom_point() +
 geom_abline(intercept = -2, slope = 0.5) +
 geom_smooth(method = 'loess', color = 'blue', se = FALSE) +
 geom_text_repel(data = subset(tmp_dat, residuals > 3.5))

ggplot(data = tmp_dat, aes(x = log2(`DSP.1001660003872.E.C04`),
                           y = residuals,
                           label = Gene)) +
 geom_point() +
 geom_abline(intercept = -2, slope = 0.5) +
 geom_smooth(method = 'loess', color = 'blue', se = FALSE) +
 geom_text_repel(data = subset(tmp_dat, log2(`DSP.1001660003872.E.C04`) > 9))

ggplot(data = tmp_dat, aes(x = log2(`DSP.1001660003872.E.C01`),
                           y = residuals,
                           label = Gene)) +
 geom_point() +
 geom_abline(intercept = -2, slope = 0.5) +
 geom_smooth(method = 'gam', color = 'blue', se = FALSE) +
 geom_text_repel(data = subset(tmp_dat, log2(`DSP.1001660003872.E.C01`) < 6 &
                                residuals > 2))

# Compare residuals against each compartment
mod <- lm(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$CAF_ROI[2]]]) ~
           log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[2]]]))
res <- resid(mod)
mod2 <- lm(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Epi_ROI[2]]]) ~
            log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[2]]]))
res2 <- resid(mod2)
mod3 <- lm(log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Epi_ROI[2]]]) ~
            log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match$Imm_ROI[2]]]))
res3 <- resid(mod2)

tmp_dat <- cbind(dfs[[2]][, c('Gene', 'DSP.1001660003872.E.C04')],
                 data.frame(CAF = as.numeric(res),
                            EPI = as.numeric(res2)))
ggplot(tmp_dat, aes(x = CAF, y = EPI, label = Gene)) +
 geom_point(color = 'lightblue') +
 geom_text_repel(data = subset(tmp_dat, CAF - EPI > 3)) +
 geom_abline(slope = 1, intercept = 0)
ggplot(tmp_dat, aes(x = CAF, y = EPI, label = Gene)) +
 geom_point(color = 'lightblue') +
 geom_text_repel(data = subset(tmp_dat, EPI - CAF > 3 | EPI > 4)) +
 geom_abline(slope = 1, intercept = 0)
ggplot(tmp_dat, aes(x = CAF, y = EPI, label = Gene)) +
 geom_point(color = 'lightblue') +
 geom_text_repel(data = subset(tmp_dat, CAF < -1 & EPI < -1 & abs(EPI - CAF) < 0.5)) +
 geom_abline(slope = 1, intercept = 0)

# distribution of residual subtraction
hist(res2-res)
dfs[[2]]$Gene[as.numeric(names(head(sort(res2-res, decreasing = TRUE), 25)))]

hist(pmin(res,res2), 30)
smoothScatter(pmax(res,res2), pmin(res,res2))

### Section 3: Develop model for background detrending ###

# Generate log2 expression data for matched ROIs
# ss_match contains the segment ID for each compartment
#   we are using here the 2nd ROI in the dataset as an example,
#     pulling the segment information for each segment
Immu <- log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match[2, 'Imm_ROI']]])
Cafs <- log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match[2, 'CAF_ROI']]])
Epit <- log2(dfs[[2]][, dfs[[3]]$Sample_ID[ss_match[2, 'Epi_ROI']]])

plot(Immu ~ Cafs)
plot(Immu ~ Epit)

# Estimate differences between the immune compartment and the max of other
# compartments

# simply subtraction of the max expression from another compartment (not final)
ggplot(tmp_dat, aes(y = Immu - pmax(Cafs, Epit), x = Immu, label = Gene)) +
 geom_point(color = 'lightblue') +
 theme_bw() +
 geom_text_repel(data = subset(tmp_dat, Immu > 6 & Immu - pmax(Cafs, Epit) > 1.5))
ggplot(tmp_dat, aes(y = Immu - pmax(Cafs, Epit), x = Epit - pmax(Cafs, Immu), label = Gene)) +
 geom_point(color = 'lightblue') +
 theme_bw() +
 geom_text_repel(data = subset(tmp_dat, Immu > 6 & Immu - pmax(Cafs, Epit) > 1.5))

# Incorporate original on-target segment expression
ggplot(tmp_dat, aes(y = (Immu - pmax(Cafs, Epit)) * Immu, x = (Cafs - pmax(Epit, Immu)) * Cafs, label = Gene)) +
 geom_point(color = 'lightblue') +
 theme_bw() +
 geom_text_repel(data = subset(tmp_dat, Immu > 6 & Immu - pmax(Cafs, Epit) > 1.5))

# Transformation calculation:
#   On target segment - max other segments) * On target segment
tmp_dat$TransImm <- (Immu - pmax(Cafs, Epit)) * Immu
tmp_dat$TransCaf <- (Cafs - pmax(Immu, Epit)) * Cafs
tmp_dat$TransEpi <- (Epit - pmax(Cafs, Immu)) * Epit

mod1 <- lm(Immu ~ Cafs)
mod2 <- lm(Immu ~ Epit)
tmp_dat$ModImm <- Immu * pmax(resid(mod1), resid(mod2))
mod1 <- lm(Cafs ~ Immu)
mod2 <- lm(Cafs ~ Epit)
tmp_dat$ModCaf <- Cafs * pmax(resid(mod1), resid(mod2))
mod1 <- lm(Epit ~ Immu)
mod2 <- lm(Epit ~ Cafs)
tmp_dat$ModEpi <- Epit * pmax(resid(mod1), resid(mod2))

cols_gg <- c(None = 'gray',
             CAF = 'magenta2',
             Immune = 'cyan3',
             Epithelial = 'green3')

tmp_dat$ImmFlg <- (Immu - pmax(Cafs, Epit)) * Immu > 10
tmp_dat$CafFlg <- (Cafs - pmax(Immu, Epit)) * Cafs > 10
tmp_dat$EpitFlg <- (Epit - pmax(Immu, Cafs)) * Epit > 10
tmp_dat$Flg <- 'None'
tmp_dat$Flg[tmp_dat$ImmFlg] <- 'Immune'
tmp_dat$Flg[tmp_dat$CafFlg] <- 'CAF'
tmp_dat$Flg[tmp_dat$EpitFlg] <- 'Epithelial'

library(GGally)
ggpairs(tmp_dat, aes(color = Flg),
        columns = c('Immu','Cafs','ModImm','ModCaf','TransImm','TransCaf')) +
 scale_color_manual(values = cols_gg) +
 scale_fill_manual(values = cols_gg)

pheatmap(tmp_dat[tmp_dat$ImmFlg | tmp_dat$CafFlg | tmp_dat$EpitFlg,
                 c('Immu','Cafs','Epit')],
         show_rownames = FALSE,
         main = "original expression, highly enriched in one ROI Type")
pheatmap(tmp_dat[tmp_dat$ImmFlg | tmp_dat$CafFlg | tmp_dat$EpitFlg,
                 c('TransImm','TransCaf','TransEpi')],
         show_rownames = FALSE,
         main = "detrended expression, highly enriched in one ROI Type")

# @RESTART HERE
ggpairs(tmp_dat, aes(color = Flg),
        columns = c('Immu','Epit','Cafs','TransImm','TransEpi','TransCaf')) +
 scale_color_manual(values = cols_gg) +
 scale_fill_manual(values = cols_gg)

ggpairs(tmp_dat, aes(color = Flg),
        columns = c('Cafs','Immu')) +
 scale_color_manual(values = cols_gg) +
 scale_fill_manual(values = cols_gg)

# Plot difference between immune expression and transformed expression
ggplot(tmp_dat, aes(x = Immu,
                    y = Immu - (log2(pmax(1, pmax(TransCaf, TransEpi))) *
                                 (Immu/max(Immu)) * 2), color = Flg, label = Gene)) +
 geom_abline(slope = 1, intercept = 0) +
 geom_point(size = 3) + theme_bw() +
 scale_color_manual(values = cols_gg) +
    labs(x = 'Original Expression (Immune Segment)',
         y = 'Detrended Expression',
         color = 'Highly enriched in:')

# Plot highly enriched genes to show that only non-immune enriched targets are
# impacted
ggplot(subset(tmp_dat, Flg != 'None'),
       aes(x = Immu,
           y = Immu - (log2(pmax(1, pmax(TransCaf, TransEpi))) *
                        (Immu/max(Immu)) * 2), color = Flg, label = Gene)) +
 geom_abline(slope = 1, intercept = 0) +
 geom_point(size = 3) + theme_bw() +
 scale_color_manual(values = cols_gg) +
 facet_wrap(~Flg)

### CALCULATE Detrending ###

# Function trsf_x - set up detrending:
# x = target ROI
# y = off-target ROI 1
# z = off-target ROI 2
# df = data frame of interest
trsf_x <- function(x, y, z, df) {
 log2(df[, x]) * (log2(df[, x]) - pmax(log2(df[, y]), log2(df[, z])))
}

# Create ROI ID information
ss_match$Epi_ROIID <- dfs[[3]]$Sample_ID[ss_match$Epi_ROI]
ss_match$CAF_ROIID <- dfs[[3]]$Sample_ID[ss_match$CAF_ROI]
ss_match$Imm_ROIID <- dfs[[3]]$Sample_ID[ss_match$Imm_ROI]

# Calculate the transformed data on a rowwise bases for each ROI
df_trans <- data.frame(row.names = dfs[[2]]$Gene)
df_clean <- data.frame(row.names = dfs[[2]]$Gene)
ss_match$Complete <- FALSE
for(i in 1:nrow(ss_match)) {
 if(all(!is.na(ss_match[i, c('Epi_ROI','CAF_ROI','Imm_ROI')]))) {
  iImm <- ss_match[i, 'Imm_ROIID']
  iCAF <- ss_match[i, 'CAF_ROIID']
  iEpi <- ss_match[i, 'Epi_ROIID']
  df_trans[, iImm] <- unlist(trsf_x(iImm, iCAF, iEpi, dfs[[2]]))
  df_trans[, iCAF] <- trsf_x(iCAF, iImm, iEpi, dfs[[2]])
  df_trans[, iEpi] <- trsf_x(iEpi, iImm, iCAF, dfs[[2]])
  ss_match[i, 'Complete'] <- TRUE
  for(seg in c('CAF','Immune','Epithelial')) {
   if(seg == 'Immune') {
    scl <- log2(dfs[[2]][, iImm]) / max(log2(dfs[[2]][, iImm]))
    df_clean[, iImm] <- log2(dfs[[2]][, iImm]) -
     (scl * 2 * log2(pmax(1, pmax(df_trans[, iCAF], df_trans[, iEpi]))))
   } else if(seg == 'CAF') {
    scl <- log2(dfs[[2]][, iCAF]) / max(log2(dfs[[2]][, iCAF]))
    df_clean[, iCAF] <- log2(dfs[[2]][, iCAF]) -
     (scl * 2 * log2(pmax(1, pmax(df_trans[, iImm], df_trans[, iEpi]))))
   } else {
    scl <- log2(dfs[[2]][, iEpi]) / max(log2(dfs[[2]][, iEpi]))
    df_clean[, iEpi] <- log2(dfs[[2]][, iEpi]) -
     (scl * 2 * log2(pmax(1, pmax(df_trans[, iCAF], df_trans[, iImm]))))
   }
  }
 }
}

# floor the data as negative transforms are un-informative
df_trans2 <- apply(df_trans, 2, pmax, 0)

### Section 4: Exploration of Detrending ###
# Calculate mean expression for each compartment
mean_Imm <- rowMeans(df_trans2[, colnames(df_trans2) %in% ss_match$Imm_ROIID])
mean_CAF <- rowMeans(df_trans2[, colnames(df_trans2) %in% ss_match$CAF_ROIID])
mean_Epi <- rowMeans(df_trans2[, colnames(df_trans2) %in% ss_match$Epi_ROIID])

# plot heatmaps before and after detrending of genes which are in the top
# expressors after detrending

# first we identify the top genes of interest:
test_plt <- mean_CAF > quantile(mean_CAF, 0.95) |
 mean_Epi > quantile(mean_Epi, 0.95) |
 mean_Imm > quantile(mean_Imm, 0.95)

# Key heatmaps:
# Before detrending:
pheatmap(log2(dfs[[2]][test_plt, -1]),
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         clustering_method = 'average',
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'original, top 5%')
# data appears to be enriched but shows caveats of need for detrending, lots of
# messy overlaps in these highly enriched genes

# After detrending
pheatmap(df_clean[test_plt, -1],
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         clustering_method = 'average',
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'detrend, top 5%')

# Try using a given threshold (slightly smaller sets)
test_plt <- mean_CAF > 2 | mean_Epi > 2 | mean_Imm > 2
pheatmap(df_clean[test_plt, -1],
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         clustering_method = 'average',
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'detrended, genes with frequent detrending (n = 1696)')

pheatmap(log2(dfs[[2]][test_plt, -1]),
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         clustering_method = 'average',
         # clustering_distance_rows =
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'original, genes with frequent detrending (n = 1696)')

# Calculate CV:
all_dat_cv <- apply(df_clean, 1, function(x) {sd(x)/mean(x)})
plot(rowMeans(df_clean), all_dat_cv, main = 'Detrended Data')

df2_cv <- apply(dfs[[2]][, -1], 1, function(x) {sd(log2(x))/mean(log2(x))})
plot(rowMeans(log2(dfs[[2]][, -1])), df2_cv, main = 'Original Data')

plot(df2_cv, all_dat_cv, xlab = 'Original CV', ylab = 'Detrended CV')
abline(a = 0, b = 1)
# Note: detrending does not decrease CV of a gene, when considered across all
# types of ROIs


table(detrend = all_dat_cv > quantile(all_dat_cv, .75),
      original = df2_cv > quantile(df2_cv, 0.75))
# note: not exact overlap but still strong bias towards remaining high CV after
# detrending

pheatmap(df_clean[(all_dat_cv > quantile(all_dat_cv, .55)) &
                   (rowMeans(df_clean) > quantile(rowMeans(df_clean), 0.4)),],
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         clustering_method = 'average',
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'detrended, high expressors (top 60%)')

pheatmap::pheatmap(dfs[[2]][all_dat_cv > quantile(all_dat_cv, .55) &
                             rowMeans(df_clean) > quantile(rowMeans(df_clean), 0.4), -1],
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = 'row',
                   breaks = seq(-3, 3, 0.05),
                   clustering_method = 'average',
                   clustering_distance_rows = 'correlation',
                   clustering_distance_cols = 'correlation',
                   color = colorRampPalette(c('purple3','black','yellow2'))(120),
                   annotation_col = dfs[[3]][, c('Treatment','Segment')],
                   main = 'original, top CV genes (45%), excluding low expressors')

# Determine the number of times a given gene was expressed in a compartment specific manner
freq_Imm <- rowSums(df_trans2[, colnames(df_trans2) %in% ss_match$Imm_ROIID] > 3)
freq_CAF <- rowSums(df_trans2[, colnames(df_trans2) %in% ss_match$CAF_ROIID] > 3)
freq_Epi <- rowSums(df_trans2[, colnames(df_trans2) %in% ss_match$Epi_ROIID] > 3)

# Immune vs CAF: nice distinction in enrichment
plot(freq_Imm, freq_CAF)

# frequency and mean expression are related:
plot(freq_Imm, mean_Imm)

# Pre- vs Post- detrended expression, Epithelial ROIs:
smoothScatter(as.matrix(log2(dfs[[2]][test_plt, colnames(df_clean)[colnames(df_clean) %in% ss_match$Epi_ROIID]])),
              as.matrix(df_clean[test_plt, colnames(df_clean) %in% ss_match$Epi_ROIID]),
              xlab = 'Original Counts, Epithelial ROIs',
              ylab = 'Detrended Counts, Epithelial ROIs')

#Exploratory heatmaps
pheatmap(df_trans2[test_CAF, -1],
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = seq(2, 20, 0.2),
         clustering_method = 'average',
         clustering_distance_cols = 'correlation',
         color = colorRampPalette(c('black','purple3','yellow2','yellow2','yellow2'))(90),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'CAF marker genes, detrended')

pheatmap(log2(dfs[[2]][test_CAF, -1]),
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = seq(-5, 30, 0.2),
         clustering_method = 'average',
         clustering_distance_cols = 'correlation',
         color = colorRampPalette(c('black','purple3','yellow2','yellow2','yellow2'))(175),
         annotation_col = dfs[[3]][, c('Treatment','Segment')],
         main = 'CAF marker genes, original')

# Heatmap of highest expressors
pheatmap(df_trans2[mean_Imm > quantile(mean_Imm, 0.9) |
                    mean_CAF > quantile(mean_CAF, 0.9) |
                    mean_Epi > quantile(mean_Epi, 0.9), -1],
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = seq(0, 40, 0.2),
         color = colorRampPalette(c('black','purple3','yellow2','yellow2','yellow2'))(200),
         annotation_col = dfs[[3]][, c('Treatment','Segment')])

### Section: Extended analysis ##
# Further exploratory analysis - qualitative
# Raising the minimum floor of expression
df_trans2 <- apply(df_trans, 2, pmax, 3)
pheatmap(df_trans2,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = seq(0, 50, 0.2),
         #clustering_method = 'average',
         # clustering_distance_rows =
         color = colorRampPalette(c('black','purple3','yellow2','yellow2','yellow2'))(200),
         annotation_col = dfs[[3]][, c('Treatment','Segment')])

pheatmap(log2(dfs[[2]][mean_Imm < quantile(mean_Imm, 0.9) |
                        mean_CAF < quantile(mean_CAF, 0.9) |
                        mean_Epi < quantile(mean_Epi, 0.9), -1]),
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         #clustering_method = 'average',
         # clustering_distance_rows =
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')])

test_qual <- apply(df_trans2[, colnames(df_trans2) %in% ss_match$Imm_ROIID], 1,
                   function(x) {sd(x) - log2(mean(x+1))})

pheatmap(log2(dfs[[2]][test_qual > 1.75,
                       -1]),#colnames(dfs[[2]]) %in% ss_match$Imm_ROIID]),
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = 'row',
         breaks = seq(-3, 3, 0.05),
         #clustering_method = 'average',
         # clustering_distance_rows =
         color = colorRampPalette(c('purple3','black','yellow2'))(120),
         annotation_col = dfs[[3]][, c('Treatment','Segment')])


pheatmap(df_trans2[test_qual > 1.75, -1],
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = seq(0, 30, 0.2),
         #clustering_method = 'average',
         # clustering_distance_rows =
         color = colorRampPalette(c('black','purple3','yellow2','yellow2','yellow2'))(150),
         annotation_col = dfs[[3]][, c('Treatment','Segment')])

tmp_dat$mean_Imm <- mean_Imm > quantile(mean_Imm, 0.9)
tmp_dat$mean_CAF <- mean_CAF > quantile(mean_CAF, 0.9)
tmp_dat$mean_Epi <- mean_Epi > quantile(mean_Epi, 0.9)

write.csv(tmp_dat, 'Example_Data_SegmentEnrichment.csv')

#### Section 6: ssGSEA analysis & comparison ####
# rerun ssGSEA and see if it holds the same or if it changes dramatically
#

new_dfs <- dfs
new_dfs[[1]] <- dfs[[1]][, c('TargetName', colnames(df_clean))]
new_dfs[[2]] <- cbind(data.frame(Gene = rownames(df_clean)), df_clean)
new_dfs[[3]] <- dfs[[3]][colnames(df_clean), ]

# Convert to entrez id and set as rownames
norm_data <- new_dfs[[2]]
norm_data <- symbol_to_entrez(path_df=norm_data, gene_col="Gene", species = species)
rownames(norm_data) <- norm_data$entrez
norm_data <- as.matrix(norm_data[, -c(1, dim(norm_data)[2])])

# # Convert de results to entrez id
# de <- de_results
# de <- symbol_to_entrez(path_df=de, gene_col="gene", species = species)
gene_ids <- as.character(rownames(norm_data))
ssGSEA_dir <- paste(outdir, "ssGSEA_Detrend", sep = "/")
dir.create(ssGSEA_dir)
ssGSEA_detrend <- gsva(norm_data,
                       pathways,
                       method="ssgsea", min.sz = 5, max.sz=500,
                       verbose=FALSE, parallel.sz=1)

mal_IDs <- colnames(ssGSEA_detrend)[colnames(ssGSEA_detrend) %in% ss_match$Epi_ROIID]
caf_IDs <- colnames(ssGSEA_detrend)[colnames(ssGSEA_detrend) %in% ss_match$CAF_ROIID]

smoothScatter(ssGSEA_results[mal_paths, mal_IDs],
              ssGSEA_detrend[mal_paths, mal_IDs])

as.data.frame(diag(cor(t(ssGSEA_results[mal_paths, mal_IDs]),
                       t(ssGSEA_detrend[mal_paths, mal_IDs]))))

smoothScatter(ssGSEA_results[caf_paths, caf_IDs],
              ssGSEA_detrend[caf_paths, caf_IDs])

diag(cor(t(ssGSEA_results[caf_paths, caf_IDs]),
         t(ssGSEA_detrend[caf_paths, caf_IDs])))
plot(ssGSEA_results['CAF_Immunomodulatory', caf_IDs],
     ssGSEA_detrend['CAF_Immunomodulatory', caf_IDs])

nepi_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'Epithelial']
ncaf_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']
use_ROIs <- nepi_ROIs
for(topic in c(mal_paths, caf_paths)) {
 if(topic %in% caf_paths) {
  use_ROIs <- ncaf_ROIs
 }
 mod <- lmer(ssGSEA_results[topic, use_ROIs] ~ treatment +
              (1|Patient), data = dfs[[3]][use_ROIs, ])
 mod_sum <- ls_means(mod, pairwise = TRUE)
 add_p <- FALSE
 if(any(mod_sum[,7] < 0.1)) {
  add_p <- TRUE
  mod_sum$plt <- mod_sum[,7] < 0.1
  mod_sum$sig <- cut(mod_sum[,7], breaks = c(0,0.01,0.05,1),
                     labels = c("**", "*", ""))
  mod_sum$x1 <- c(1,1,1,2,2,3)
  mod_sum$x2 <- c(2,3,4,3,4,4)
  mod_sum <- as.data.frame(mod_sum)
  y_range <- diff(range(scale(ssGSEA_results[topic, use_ROIs])))
  y_max <- max(scale(ssGSEA_results[topic, use_ROIs]))
  mod_sum <- subset(mod_sum, plt)
  for(i in 1:nrow(mod_sum)) {
   mod_sum[i, 'y'] <- y_range*(i * 0.065) + y_max
  }
 }

 plt <- ggplot(dfs[[3]][use_ROIs, ],
               aes(x = treatment, fill = treatment)) +
  geom_boxplot(aes(y = scale(ssGSEA_results[topic, use_ROIs]))) +
  theme_bw() +
  scale_fill_jama() +
  labs(x = 'Treatment Group',
       y = paste0(gsub('Malignant_', '', topic), ' Enrichment'),
       title = paste0(gsub('Malignant_', '', topic), ' NMF Score'))
 if(add_p) {
  plt <- plt +
   annotate('segment',
            x = mod_sum$x1, xend = mod_sum$x2,
            y = mod_sum$y, yend = mod_sum$y) +
   annotate('text', x = rowMeans(mod_sum[, c('x1','x2')]),
            y = mod_sum$y + y_range*.01,
            label = paste0("P = ", signif(mod_sum[, 9], 3), mod_sum[, 'sig']),
            hjust = 0.5, vjust = 0, size = 3)
 }
 print(plt)
 ggsave(paste0('ssGSEA_Detrend/', gsub('\\/','-',topic), '_TreatmentAnalysis.png'), plt, 'png',
        width = 5, height = 6)
}

saveRDS(new_dfs, 'dfs_detrendApproach21-6-2.RDS')
saveRDS(ssGSEA_detrend, 'ssGSEA_detrendApproach21-6-2.RDS')


#re-exponentiate
new_dfs_exp <- new_dfs
new_dfs_exp[[2]] <- cbind(new_dfs_exp[[2]][,1], 2^new_dfs_exp[[2]][,-1])

saveRDS(new_dfs_exp, 'dfs_detrendApproach21-6-4_exp.RDS')


test_tsne <- Rtsne(X = t(new_dfs[[2]][,-1]), perplexity = 40, partial_pca = TRUE, PCA = TRUE)
new_dfs[[3]]$Tsne2 <- test_tsne$Y[,2]
new_dfs[[3]]$Tsne1 <- test_tsne$Y[,1]
ggplot(new_dfs[[3]], aes(x = Tsne1, y = Tsne2, color = Segment)) + geom_point() + theme_bw()
ggplot(tsne$samples, aes(x = tsne$X1, y = tsne$X2, color = Segment)) + geom_point() + theme_bw()

test_uamp <- umap::umap(t(new_dfs[[2]][,-1]), preserve.seed = TRUE)
new_dfs[[3]]$Umap1 <- test_uamp$layout[,1]
new_dfs[[3]]$Umap2 <- test_uamp$layout[,2]
ggplot(new_dfs[[3]], aes(x = Umap1, y = Umap2, color = Segment)) + geom_point() + theme_bw() + labs(title = 'Detrended')
ggplot(umap$samples, aes(x = umap$X1, y = umap$X2, color = Segment)) + geom_point() + theme_bw() + labs(title = 'Q3 Norm only')

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = Umap1, y = Umap2,
           size = ssGSEA_detrend['CAF_Immunomodulatory',
                                 new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']],
           color = ssGSEA_detrend['CAF_Immunomodulatory',
                                  new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']])) +
 geom_point() + theme_bw() + guides(color = FALSE, size = FALSE)

# iCAF vs myCAF: detrend
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           y = scale(ssGSEA_detrend['CAF_Myofibroblastic',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'iCAF', y = 'myCAF', title = 'Detrended') +
 scale_color_jama() +
 theme(aspect.ratio = 1)

# iCAF vs myCAF: Q3
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Immunomodulatory',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           y = scale(ssGSEA_results['CAF_Myofibroblastic',
                                    new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'iCAF', y = 'myCAF', title = 'Q3 Norm') +
 scale_color_jama() +
 theme(aspect.ratio = 1)


table(ss_match$Complete)
# iCAF vs Mesenchymal
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_detrend['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'iCAF, CAF', y = 'Mesenchymal, Epi', title = 'Detrended') +
 scale_color_jama() +
 theme(aspect.ratio = 1)
cor(scale(ssGSEA_detrend['CAF_Immunomodulatory',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_detrend['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_results['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'iCAF, CAF', y = 'Mesenchymal, Epi', title = 'Q3 Norm') +
 scale_color_jama() +
 theme(aspect.ratio = 1)
cor(scale(ssGSEA_results['CAF_Immunomodulatory',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_results['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

# myCAF vs Mesenchymal: detrend
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Myofibroblastic',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_detrend['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'myCAF, CAF', y = 'Mesenchymal, Epi', title = 'Detrended') +
 scale_color_jama() +
 theme(aspect.ratio = 1)
cor(scale(ssGSEA_detrend['CAF_Myofibroblastic',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_detrend['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_results['CAF_Myofibroblastic',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = scale(ssGSEA_results['Malignant_Mesenchymal',
                                    subset(ss_match, Complete)$Epi_ROIID]),
           color = treatment)) +
 geom_point(size = 2.5) + theme_bw() + #guides(color = FALSE, size = FALSE) +
 labs(x = 'myCAF, CAF', y = 'Mesenchymal, Epi', title = 'Q3 Norm') +
 scale_color_jama() +
 theme(aspect.ratio = 1)
cor(scale(ssGSEA_results['CAF_Myofibroblastic',
                         subset(ss_match, Complete)$CAF_ROIID]),
    scale(ssGSEA_results['Malignant_Mesenchymal',
                         subset(ss_match, Complete)$Epi_ROIID]))

# CD45 vs iCAF
ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(x = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                    subset(ss_match, Complete)$CAF_ROIID]),
           y = unlist(new_dfs[[2]]['CD68',
                                   subset(ss_match, Complete)$CAF_ROIID]),
           color = Treatment)) + geom_point() + theme_bw()


ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(fill = scale(ssGSEA_detrend['CAF_Immunomodulatory',
                                       subset(ss_match, Complete)$CAF_ROIID]) > 0,
           y = unlist(new_dfs[[2]]['PTPRC',
                                   subset(ss_match, Complete)$CAF_ROIID]),
           x = Treatment)) +
 geom_boxplot() + theme_bw() +
 labs(x = 'iCaf', y = 'PTPRC (CD45) in CAF ROIs', title = 'Detrended', fill = 'High iCAF') +
 ylim(c(0,7))

ggplot(subset(new_dfs[[3]], Segment == 'CAF'),
       aes(fill = scale(ssGSEA_results['CAF_Immunomodulatory',
                                       subset(new_dfs[[3]], Segment == 'CAF')$Sample_ID]) > 0,
           y = unlist(log2(dfs[[2]]['PTPRC',
                                    subset(new_dfs[[3]], Segment == 'CAF')$Sample_ID])),
           x = Treatment)) + geom_boxplot() + theme_bw() +
 labs(x = 'iCaf', y = 'PTPRC (CD45) in CAF ROIs',  title = 'Q3 Norm', fill = 'High iCAF')+
 ylim(c(0,7))

