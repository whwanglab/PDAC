library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggsci)

### INTRA-PATIENT VERSUS INTER-PATIENT HETEROGENEITY: FIGURE 6C
new_dfs <- readRDS('dfs_detrendApproach21-6-2.RDS')
ssGSEA_detrend <- readRDS('ssGSEA_detrendApproach21-6-2.RDS')

# capture relevant pathway names if not in environment
all_paths <- rownames(ssGSEA_detrend)[!grepl("_[0-9]", rownames(ssGSEA_detrend))]

epi_lineage <- c('Malignant_Classical-like', 'Malignant_Mesenchymal',
                 'Malignant_Squamoid', 'Malignant_Basaloid', 'Malignant_Acinar-like', 'Malignant_Neuroendocrine-like',
                 'Malignant_Neuronal-like')
caf_paths <- all_paths[grepl('CAF', all_paths)]

# Scale the ssGSEA scores within each compartment as scores from different
# compartments should not be considered during scaling
# 
# these will be re-aligned with the original data
epi_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'Epithelial']
scaled_detrendEpi <- t(apply(ssGSEA_detrend[all_paths, epi_ROIs], 1, scale))
colnames(scaled_detrendEpi) <- epi_ROIs

CAF_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'CAF']
scaled_detrendCAF <- t(apply(ssGSEA_detrend[all_paths, CAF_ROIs], 1, scale))
colnames(scaled_detrendCAF) <- CAF_ROIs

immune_ROIs <- new_dfs[[3]]$Sample_ID[new_dfs[[3]]$Segment == 'Immune']
scaled_detrendImm <- t(apply(ssGSEA_detrend[all_paths, immune_ROIs], 1, scale))
colnames(scaled_detrendImm) <- immune_ROIs

# merge back to gether
scaled_detrend <- cbind(scaled_detrendEpi, scaled_detrendCAF, scaled_detrendImm)
scaled_detrend <- scaled_detrend[, new_dfs[[3]]$Sample_ID]

# add pathways to data frame
new_dfs[[3]][, all_paths] <- t(scaled_detrend)

# melt the data frame for plotting
ann_melt <- melt(new_dfs[[3]], measure.vars = all_paths,
                 id.vars = colnames(new_dfs[[3]])[!colnames(new_dfs[[3]]) %in% all_paths],
                 variable.name = 'Topic', value.name = 'NES')

# Calculating summary statistics from molten data frame
# * want Patient, Segment, Treatment information collapsed; across ROIS
# * calculate IQR & Range for use
Topic_IQR <- dcast(ann_melt, Patient + Segment + Treatment ~ Topic, fun.aggregate = IQR)
Topic_Range <- dcast(ann_melt, Patient + Segment + Treatment ~ Topic, fun.aggregate = function(x) {diff(range(x))})
Topic_mean <- dcast(ann_melt, Patient + Segment + Treatment ~ Topic, fun.aggregate = mean)

# melt the Topic IQR & Range datasets for collapse:
Topic_IQR_m <- melt(Topic_IQR, measure.vars = all_paths,
                    id.vars = colnames(Topic_IQR)[!colnames(Topic_IQR) %in% all_paths],
                    variable.name = 'Topic',
                    value.name = 'IQR')
Topic_Range_m <- melt(Topic_Range, measure.vars = all_paths,
                      id.vars = colnames(Topic_Range)[!colnames(Topic_Range) %in% all_paths],
                      variable.name = 'Topic',
                      value.name = 'Range')

# Create a container dataframe with mean IQR values
comp_II <- data.frame(Interpatient = c(apply(Topic_mean[Topic_mean$Segment == 'Epithelial', epi_lineage],2,IQR),
                                       apply(Topic_mean[Topic_mean$Segment == 'CAF', caf_paths],2,IQR)),
                      Intrapatient = NA,
                      Topic = c(epi_lineage, caf_paths),
                      Type = c(rep('Epithalial', length(epi_lineage)),
                               rep('CAF', length(caf_paths))))
rownames(comp_II) <- c(epi_lineage, caf_paths)

# Fill in ROI level means
for(i in 1:nrow(comp_II)) {
  topic <- rownames(comp_II)[i]
  if(topic %in% epi_lineage) {
    comp_II[i, 'Intrapatient'] <- mean(Topic_IQR_m$IQR[Topic_IQR_m$Segment == 'Epithelial' &
                                                         Topic_IQR_m$Topic == topic])
  } else {
    comp_II[i, 'Intrapatient'] <- mean(Topic_IQR_m$IQR[Topic_IQR_m$Segment == 'CAF' &
                                                         Topic_IQR_m$Topic == topic])
  }
}

# Strip out prefixes
comp_II$Topic <- gsub('Malignant_', '', comp_II$Topic)
comp_II$Topic <- gsub('CAF_', '', comp_II$Topic)

# Plot Figure (scatter plot of intra vs interpatient IQR)
ggplot(comp_II, aes(x = Intrapatient, y = Interpatient, color = Type, label = Topic)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 'dashed') +
  theme_bw() +
  labs(x = 'Mean IQR within a Patient', y = 'IQR of mean scores across Patients') +
  theme(aspect.ratio = 1) +
  geom_text_repel(point.padding = .2, box.padding = .25, min.segment.length = .1, color = 'black')

# Another option for analysis would be to use the global IQR for all ROIs, this
# is even more heavily favored towards inter-patient variability
dcast(ann_melt, Segment ~ Topic, fun.aggregate = IQR)

