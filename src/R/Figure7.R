library(reshape2)
library(ggplot2)
library(ggsci)
library(lmerTest)

### FIGURE 7A AND 7B
new_dfs <- readRDS('dfs_detrendApproach21-6-2.RDS')
ssGSEA_detrend <- readRDS('ssGSEA_detrendApproach21-6-2.RDS')

# capture relevant pathway names if not in environment
all_paths <- rownames(ssGSEA_detrend)[!grepl("_[0-9]", rownames(ssGSEA_detrend))]

epi_lineage <- c('Malignant_Classical-like', 'Malignant_Mesenchymal',
                 'Malignant_Squamoid','Malignant_Neuroendocrine-like',
                 'Malignant_Neuronal-like','Malignant_Basaloid',
                 'Malignant_Acinar-like')
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

# Figure 3C panels:
# Epithelial
plt <- ggplot(subset(ann_melt, Topic %in% epi_lineage &
                Segment == 'Epithelial' & 
                treatment %in% c('Untreated', 'CRT')),
       aes(x = gsub('Malignant_', '', Topic), fill = treatment, y = NES)) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0)) +
  labs(x = 'Malignant Topics', y = 'Enrichment', fill = '') +
  scale_fill_npg()
print(plt)
ggsave('Treatment vs Epithelial Topics (UvsCRT).png', plt, 'png', width = 7, height = 5)

# CAF
plt <- ggplot(subset(ann_melt, Topic %in% caf_paths &
                Segment == 'CAF' & 
                treatment %in% c('Untreated', 'CRT')),
       aes(x = gsub('CAF_', '', Topic), fill = treatment, y = NES)) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = -45, hjust = 0)) +
  labs(x = 'CAF Topics', y = 'Enrichment', fill = '') +
  scale_fill_npg()

print(plt)
ggsave('Treatment vs CAF Topics (UvsCRT).png', plt, 'png', width = 4.5, height = 5)

# Calculate P-values for association with treatment (note ls_means needed if
# extending to CRTL/N as well)
CRT_Path <- data.frame(row.names = c(epi_lineage, caf_paths),
                       Topic = c(epi_lineage, caf_paths),
                       FC = NA,
                       P = NA)
use_ROIs <- epi_ROIs
for(topic in c(epi_lineage, caf_paths)) {
  if(topic %in% caf_paths) {
    use_seg <- 'CAF'
  } else {
    use_seg <- 'Epithelial'
  }
  fx <- paste0('`', topic, '` ~ treatment + (1|Patient)')
  mod <- lmer(fx, data = subset(new_dfs[[3]],
                                Segment == use_seg &
                                  Treatment %in% c('Untreated', 'CRT')))
  CRT_Path[topic, 'FC'] <- coef(summary(mod))[2,1]
  CRT_Path[topic, 'P'] <- coef(summary(mod))[2,5]
}

write.table(CRT_Path, 'TopicVsTreatment.csv')
