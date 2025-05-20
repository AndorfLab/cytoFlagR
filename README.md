# ${{\color{#236CA6}{\Large{\textnormal{\textsf{cytoFlagR}}}}}}$
#### A comprehensive framework to objectively assess high-dimensional cytometry data for batch effects

## About
With the advancement of cytometry techniques, researchers have the ability to create longitudinal immunological studies that generate high-dimensional cytometry datasets. However, acquiring this data over a long period of time and analyzing this data over multiple experimental runs leads to technical variations between these experimental batches, commonly known as ‘batch effects’. Several different approaches to minimize/correct for these batch effects exist. However, there are currently no computational approaches to objectively evaluate data and identify these potential batch effects. cytoFlagR is a unique, automated approach that assesses high-parameter cytometry datasets and provides the user with a comprehensive evaluation of problematic batches and markers present in their data. Important caveat: while this tool is designed to be primarily applied on control samples, it is also able to assess biological samples. However, users should take into consideration the inherent biological variability between their biological samples while interpreting the outcome of cytoFlagR on their data. 

<dl>
<h2>This tool consists of five main steps</h2>

1. Transformation of FCS files of control samples containing live, pre-gated single cells and preliminary visual inspection of data
  
 2. An Inter Quartile Range (IQR) based assessment to check for batch effects for each marker in each control sample provided for the negative and positive populations as well as the percent of positive cells. <br>
 
 3. An Earth Mover’s Distance based assessment to identify batch issues within every marker and control sample <br>
 
 4. A comprehensive summary of the assessment metrics described indicating potentially problematic batches and markers present in the data <br>
 
 5. An unsupervised clustering-based assessment to highlight batch issues present within the unique cell populations of the data <br>
 
</dl>

## Download cytoFlagR
```
Download zip from the `<>Code` button above

Or use command line to download cytoFlagR using -

git clone https://github.com/AndorfLab/cytoFlagR.git
```

## Dependencies
This tool was developed using R version 4.4.1, the version required to run this tool.
Download and install R from [here](https://cran.r-project.org/)

cytoFlagR requires several R and BioConductor packages to run.

Use the package_installer.R function to install the missing required packages from the list below
```
source("package_installer.R")

### required CRAN packages
requiredPackages<-c("dplyr","scales","tidyr","reshape2","readr","matrixStats","readxl",
                    "ggplot2","ggpubr","ggridges","MASS","RColorBrewer","cowplot",
                    "randomcoloR","ggrepel","emdist","circlize","gridExtra","stats",
                    "LaplacesDemon","pheatmap","umap","progress","crayon","patchwork",
                    "ggpmisc","viridis","tidyverse","shiny","shinyjs","DT","bslib","shinyjs",
                    "cluster","rstudioapi","grid","shinycssloaders","shinyWidgets")
### install
package_installer(requiredPackages)
### check if the packages can be loaded
lapply(requiredPackages, require, character.only = TRUE)

### install Bioconductor package installer
checkBiocManager_install()

### required BioConductor packages
required_BioconductorPackages<-c("flowCore","FlowSOM","ComplexHeatmap","limma","ConsensusClusterPlus")
### install
BioC_package_installer(required_BioconductorPackages)
### check if the packages can be loaded
lapply(required_BioconductorPackages, require, character.only = TRUE)
```
### The code to run cytoFlagR is found in the `R/` folder when you download cytoFlagR and detailed instructions for using cytoFlagR are available on the [wiki page](https://github.com/AndorfLab/cytoFlagR/wiki/cytoFlagR-wiki/).
