## R scripts for downstream DSP data analysis

### Version 0.4

### Latest Release
Latest release can be found on Box - ask Jason (jreeves@nanostring), Nicole (nortogero@nanostring.com), or Maddy (mgriswold@nanostring) for access

### Bugs and Feature Requests 
Please submit bugs and feature requests here: 

Microsoft Teams > BioInfo Help > Questions

### If you need assistance
Please ask Jason or Maddy. 


#### File Structure
	/OutputDir

		/log

			log file
			
			copy of config file

		/RData

			RData file

		/Results
		
		  Normalized expression matrix used in analysis
		
			/cellTyping
				
				cell abundance barplot
				
				cell abundance heatmap
				
				cell abundance in xy space
				
				cell decon abundance estimates
				
				cell decon results.RData
				
			/clustering

				/genes

					Eigenvalues

					Coordinates.csv

					Gene vs components.csv

					PCA

					tSNE

					UMAP

				/samples

					Eigenvalues

					Coordinates.csv

					Samples vs components.csv

					PCA

					tSNE

					UMAP

					clustering dendrogram

				DE heatmaps

			/de

				/modelFit

					qq plots

					residual plots

				ROC curves

				de results.csv

				volcano plots

				DE gene set summary

			/pathway
				
				/individual DE tests 
				
					GSEA results.csv

					Main Independent Pathways (figure & .csv)

					Multiple Annot

					ORA results.csv

					Overall coverage

					Overall coverage map

					Overall coverage with ratios

					ssGSEA results.csv

					Top pathway down

					Top pathway up

					Top pathway table

			/qc

				Normalized distributions

				Normalization heatmaps

				Normalization comparison heatmaps

				HKgeomean

				LOD across AOIs

				Norm Factors

				Pairs plots
				
				QC_linePlot
				
				QC_values
				
				QC_genesPerAOI
				
			/spatialPlots
			
				ROI_ExpressionPlots
