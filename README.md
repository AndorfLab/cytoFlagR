# ${{\color{#236CA6}{\Large{\textsf{cytoFlagR}}}}}$
A comprehensive framework to objectively assess high-dimensional cytometry data for batch effects

## About
With the advancement of cytometry techniques, researchers have the ability to create longitudinal immunological studies that generate high-dimensional cytometry datasets. However, acquiring this data over a long period of time and analyzing this data over multiple experimental runs leads to technical variations between these experimental batches, commonly known as ‘batch effects’. Several different approaches to minimize/correct for these batch effects exist. However, there are currently no computational approaches to objectively evaluate data and identify these potential batch effects. cytoFlagR is a unique, automated approach that assesses high-parameter cytometry datasets and provides the user with a comprehensive evaluation of problematic batches and markers present in their data. Important caveat: while this tool is designed to be primarily applied on control samples, it is also able to assess biological samples. However, users should take into consideration the inherent biological variability between their biological samples while interpreting the outcome of cytoFlagR on their data. 

## Download cytoFlagR
```
install.packages("devtools") ## ignore this line if you have installed devtools previously
library(devtools)
devtools::install_github('bioinfoSE/cytoFlagR')
```

## Dependencies
This tool was developed using R version 4.4.1, the version required to run this tool.
Download and install R from [here](https://cran.r-project.org/)

cytoFlagR requires several R and BioConductor packages to run. Install R [BioConductor](https://www.bioconductor.org/install/) package manager with the instructions provided (click on BioConductor to be re-directed).
Use the package_installer.R function to install the missing required packages from the list below
```
source("package_installer.R")

### required CRAN packages
requiredPackages<-c("dplyr","scales","tidyr","reshape2","readr","matrixStats","readxl",
                    "ggplot2","ggpubr","ggridges","MASS","RColorBrewer","cowplot",
                    "randomcoloR","ggrepel","emdist","circlize","gridExtra","stats",
                    "LaplacesDemon","pheatmap","umap","progress","crayon","patchwork",
                    "ggpmisc","viridis","tidyverse","shiny","shinyjs","DT","bslib","shinyjs",
                    "cluster","rstudioapi")
### install
package_installer(requiredPackages)

### required BioConductor packages
required_BioconductorPackages<-c("flowCore","FlowSOM","ComplexHeatmap","limma","ConsensusClusterPlus")
### install
BioC_package_installer(required_BioconductorPackages)

```
### Required file formats and inputs for cytoFlagR can be found in the [wiki page](https://github.com/bioinfoSE/cytoFlagR/wiki/cytoFlagR-wiki/).
