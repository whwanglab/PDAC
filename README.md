# Whole transcriptome digital spatial profiling of pancreatic cancer

Pancreatic ductal adenocarcinoma (PDAC) is a highly lethal and treatment-refractory cancer. Molecular stratification in pancreatic cancer remains rudimentary and does not yet inform clinical management or therapeutic development. We construct a high-resolution molecular landscape of the multicellular subtypes and spatial communities that compose PDAC using single-nucleus RNA-seq and whole-transcriptome digital spatial profiling (DSP) of 43 primary PDAC tumor specimens that either received neoadjuvant therapy or were treatment-naïve. We uncovered recurrent expression programs across malignant cells and fibroblasts, including a newly-identified neural-like progenitor malignant cell program that was enriched after chemotherapy and radiotherapy and associated with poor prognosis in independent cohorts. Integrating spatial and cellular profiles revealed three multicellular communities with distinct contributions from malignant, fibroblast, and immune subtypes: classical, squamoid-basaloid, and treatment-enriched. Our refined molecular and cellular taxonomy can provide a framework for stratification in clinical trials and serve as a roadmap for therapeutic targeting of specific cellular phenotypes and multicellular interactions.

In this repository, we present the analysis conducted for the whole transcriptome DSP experiments. 

Our Nature Genetics manuscript (in press) will be available soon. The preprint can be found here. Data is uploaded to GEO and can be found here. Code for the single-nucleus RNA-seq analysis can be found here. 


## Data analysis

### Data preprocessing

FASTQ files (uploaded to GEO) for DSP were aggregated into count matrices using the azorius and hydra pipeline. Normalized expression was detrended to model cell-type specific expression. 


### Program scoring and correlation analysis

Programs were scored for each DSP sample within each ROI using single-sample gene set enrichment analysis (ssGSEA), which were transformed using the z-score. Unsupervised hierarchical clustering was performed on all features (malignant programs, CAF programs, deconvolved immune cell type proportions, compartment areas within ROI) using the Pearson correlation distance and average linkage. Cell deconvolution analysis was performed using the SpatialDecon v0.99.1 package (https://github.com/Nanostring-Biostats/SpatialDecon/). Analysis of expression or program scores used linear mixed effects models to control for multiple sampling within a slide, using Satterthwaite's approximation for degrees of freedom for p-value calculation. Correlation coefficients were calculated using the Spearman rank correlation.


### Receptor ligand analysis

Known receptor-ligand pairs were obtained from CellPhoneDB v2.0 with potential receptor-ligand pairs quantified using the Spearman rank correlation between paired segments within the same ROI across all ROIs with said pairs. Interactions were calculated for non-self (juxtacrine) and self (autocrine) occurring within the same segment. Receptor-ligand interactions were calculated separately for untreated and CRT specimens to determine interactions that are differential between conditions. All analyses were two-sided and used a significant level of p-value ≤ 0.05 and were adjusted for multiple testing where appropriate using the false discovery rate.

