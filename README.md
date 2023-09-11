# GastricTME

This repository contains matadata and codes necessary for analysis of the scRNA-seq data of human STAD presented in Wang et al. Cancer Cell (2023). 

Cells with low complexity libraries or likely cell debris have been removed. Batch effect was corrected by Harmony and verified by k-BET. 

Please contact LWang22@mdanderson.org and RWang12@mdanderson.org with any questions.



## Downloading the data

Cell level metadata is available in the provided /Input_data/meta.Rdata, which contains sample, tissue type, major cell types and detailed cell types. 

Specific input data for each figure is also inclued in /Input_data.

The processed data are uploading to GEO: GSE234129.

## Data visualization

### Requirements

Tested on macOS Big Sur

1. R version: 4.1.2
2. R packages
   - ggplot2
   - data.table
   - Seurat
   - dplyr
   - tidyr
   - ggpubr
   - RColorBrewer
   - pheatmap
   - ggsignif
