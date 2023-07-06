# GastricTME

This repository contains matadata and codes necessary for analysis of the scRNA-seq data of human STAD presented in Wang et al. Cancer Cell (2023). 

Cells with low complexity libraries or likely cell debris have been removed. Batch effect was corrected by Harmony and verified by k-BET. 

Please contact LWang22@mdanderson.org and RWang12@mdanderson.org with any questions.



## Downloading the data

cell2subtype have been organized into R objects that can be downloaded from /Input_data/.

Cell level metadata is available in the provided /Input_data/Cell metaData.rda, which contains QC statistics of cells, clinical information of associated samples and cell types. 

Cell_ID to cell type association for all the TME cells are available in the provided /input_data/cellType_All.rda.

All the processed data are uploading to GEO: GSE234129.

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
   - monocle
   - smoother
   - pheatmap
   - Hmisc
   - monocle3
   - ggrepel
   - CytoTRACE
   - shazam

3. igblast_1.17.1
4. Change-O toolkit
