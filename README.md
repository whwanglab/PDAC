# Whole transcriptome digital spatial profiling of pancreatic cancer

Pancreatic ductal adenocarcinoma (PDAC) is a highly lethal and treatment-refractory cancer. Molecular stratification in pancreatic cancer remains rudimentary and does not yet inform clinical management or therapeutic development. We construct a high-resolution molecular landscape of the multicellular subtypes and spatial communities that compose PDAC using single-nucleus RNA-seq and whole-transcriptome digital spatial profiling (DSP) of 43 primary PDAC tumor specimens that either received neoadjuvant therapy or were treatment-naïve. We uncovered recurrent expression programs across malignant cells and fibroblasts, including a newly-identified neural-like progenitor malignant cell program that was enriched after chemotherapy and radiotherapy and associated with poor prognosis in independent cohorts. Integrating spatial and cellular profiles revealed three multicellular communities with distinct contributions from malignant, fibroblast, and immune subtypes: classical, squamoid-basaloid, and treatment-enriched. Our refined molecular and cellular taxonomy can provide a framework for stratification in clinical trials and serve as a roadmap for therapeutic targeting of specific cellular phenotypes and multicellular interactions.

In this repository, we present the analysis conducted for the whole transcriptome DSP experiments. 

Our Nature Genetics manuscript (in press) will be available soon. The preprint can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.25.267336v1.full). Raw and processed data can be found at GEO under accession number GSE199102. Code for the single-nucleus RNA-seq analysis can be found [here](https://github.com/karthikj89/humanpdac). 


## Data analysis

### Data preprocessing

FASTQ files (uploaded to GEO) for DSP were aggregated into count matrices using the [azorius](https://github.com/whwanglab/PDAC/tree/main/src/R/azorius) and [hydra](https://github.com/whwanglab/PDAC/tree/main/src/R/hydra) pipeline. Normalized expression was [detrended](https://github.com/whwanglab/PDAC/blob/main/src/R/Detrending_and_ssGSEA.R) to model cell-type specific expression. 

Normalized data can be found [here](https://github.com/whwanglab/PDAC/blob/main/src/R/dfs_Q321-6.RDS). Detrended data can be found at GEO under accession number GSE199102. In this repository, normalized and detrended data are referenced in analyses described below.

### Cell type deconvolution and program scoring

Programs were scored for each DSP sample within each ROI using [ssGSEA](https://github.com/whwanglab/PDAC/blob/main/src/R/Detrending_and_ssGSEA.R), which were transformed using the z-score. [Unsupervised hierarchical clustering](https://github.com/whwanglab/PDAC/blob/main/src/R/FigureED9.R) was performed on all features (malignant programs, CAF programs, deconvolved immune cell type proportions, compartment areas within ROI) using the Pearson correlation distance and average linkage. Code for the immune cell type deconvolution analysis can be found [here](https://github.com/whwanglab/PDAC/blob/main/src/R/CellTypeDeconvolution.Rmd).

ssGSEA program scores can be found [here](https://github.com/whwanglab/PDAC/blob/main/src/R/ssGSEA_detrendApproach21-6-2.RDS).

### Receptor ligand analysis

Known receptor-ligand pairs were obtained from [CellPhoneDB v2.0](https://github.com/Teichlab/cellphonedb) with potential receptor-ligand pairs quantified using the Spearman rank correlation between paired segments within the same ROI across all ROIs with said pairs. Interactions were calculated for non-self (juxtacrine) and self (autocrine) occurring within the same segment. Receptor-ligand interactions were calculated separately for untreated and CRT specimens to determine interactions that are differential between conditions. All analyses were two-sided and used a significant level of p-value &gt; 0.05 and were adjusted for multiple testing where appropriate using the false discovery rate.

Code for this analysis can be found [here](https://github.com/whwanglab/PDAC/blob/main/src/R/LigandReceptorAnalysis_Figure8.R).

