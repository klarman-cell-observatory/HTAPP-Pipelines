Description of scripts
```
source(paste0(user.path, "/code/orr_plotutils_04182019.R"))
source(paste0(user.path, "/code/orr_seuratutils_04182019.R"))
source(paste0(user.path, "/code/orr_color_04182019.R"))
source(paste0(user.path, "/code/plot.umap.feature.R"))
```

1. htapp_methods_lung_toPostJuly2020.Rmd
This script processes the sample out of cumulus and sets files up for annotation, inferCNV, and figure generation

2. htapp_methods_lung_annotate_toPostJuly2020.Rmd
This script is used to annotate the samples, using figures from the first script as guidance

3. htapp_methods_lung_figures_toPostJuly2020.Rmd
This script generates figures and tables of QC stats, annotations, etc. to facilitate protocol comparison and data interpretation

4. htapp_methods_lung_inferCNV_toPostJuly2020.Rmd
Script to run InferCNV; this script was later converted to a .R file and run on our compute cluster
