data_management_packages_to_install <- c("data.table", "tidyverse" ,"readxl","purrr","readr","scales","msa","jsonlite", "broom")
data_visualisation_packages_to_install <- c("ggplot2","cowplot","pheatmap","viridisLite","ggseqlogo","circlize","ggpubr", "ggnewscale", "legendry", "ggrepel", "ggbeeswarm", "ggmosaic")
structural_packages_to_install <- c("bio3d","NGLVieweR", "htmlwidgets","htmltools")

packages_to_install <- c(data_management_packages_to_install, 
                         data_visualisation_packages_to_install,
                         structural_packages_to_install)

# Check and install packages using a for loop
for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

BiocManager_packages_to_install <- c("Biostrings","pwalign","ComplexHeatmap", "DECIPHER")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (package in BiocManager_packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
  library(package, character.only = TRUE)
}