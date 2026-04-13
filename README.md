# ${{\color{#236CA6}{\Large{\textnormal{\textsf{CytoBatchFlagR}}}}}}$ 

## Overview
The CytoBatchFlagR tool intakes cytometry data and outputs information about potential technical batch issues in that data. The tool is automated and provides multiple metrics and visualizations that allow the user to access the quality of their data. 

While this tool is designed to be primarily applied on control samples, it is also able to assess biological samples. However, users should take into consideration the inherent biological variability between their biological samples while interpreting the outcome of CytoBatchFlagR.

This tool consists of 6 main parts (each link leads to the relevant Wiki section):

1.	[Pre-processing and Visual Assessment of the Data](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-1:-Pre%E2%80%90processing-and-Visual-Assessment-of-the-Data)
    * Transforms control samples (FCS files) and creates various plots (a MDS, UMAP, and bar plot) for initial data analysis 
2.	[Interquartile Range (IQR)-Based Assessment](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-2:-Interquartile-Range-(IQR)%E2%80%90Based-Assessment)
    * Analyzes each marker in each control sample for the negative and positive populations, as well as the percent of positive cells
3.	[Earth Mover’s Distance (EMD)-Based Assessment](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-3:-Earth-Mover%E2%80%99s-Distance-(EMD)%E2%80%90Based-Assessment)
    * Uses the EMD equation for pairwise comparisons between every marker in the control samples
4.	[Summary of Results (for Steps 2 and 3)](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-4:-Summary-of-Results-(for-Steps-2-and-3))
    * Creates a plot that summarizes the flags from the IQR-based and EMD-based assessments
5.	[Unsupervised Clustering-Based Assessment](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-5:-Unsupervised-Clustering%E2%80%90Based-Assessment)
    *  Clusters the data to highlight batch issues present within the cell populations
6. [Summary of Results (for Steps 2, 3, and 5)](https://github.com/AndorfLab/CytoBatchFlagR/wiki/Step-6:-Summary-of-Results-(for-Steps-2,-3,-and-5))
   * Creates a plot that summarizes the flags from the IQR-based, EMD-based, and cluster-based assessments
      
Aditional information about the required input files and the example data are also available in the [Wiki](https://github.com/AndorfLab/CytoBatchFlagR/wiki). 

## Download
The ZIP file containing all code can be downloaded by clicking on the *<>Code* button above.

Alternatively, the command line can be used to download the tool. Just copy and paste:

```
git clone https://github.com/AndorfLab/CytoBatchFlagR.git
```

Once you download CytoBatchFlagR, the code to run the various functions can be found in the `R/` folder. 

## Dependencies
This tool was developed using R version 4.4.1. Other versions may not be compatable with running the tool. Download and install R [here](https://cran.r-project.org/).

CytoBatchFlagR requires several R and BioConductor packages to run.

First, make sure you set your directory as the `R/` folder of CytoBatchFlagR:

```
setwd("C:/you/directory/here/CytoBatchFlagR-main/R")
```

Once the directory is correctly set, you can use the package_installer.R function to install the required packages:
```
source("package_installer.R")

# required CRAN packages
requiredPackages<-c("dplyr","scales","tidyr","reshape2","readr","matrixStats","readxl",
                    "ggplot2","ggpubr","ggridges","MASS","RColorBrewer","cowplot",
                    "randomcoloR","ggrepel","emdist","circlize","gridExtra","stats",
                    "LaplacesDemon","pheatmap","umap","progress","crayon","patchwork",
                    "ggpmisc","viridis","tidyverse","shiny","shinyjs","DT","bslib","shinyjs",
                    "cluster","rstudioapi","grid","shinycssloaders","shinyWidgets")
# install
package_installer(requiredPackages)

# check if the packages can be loaded
lapply(requiredPackages, require, character.only = TRUE)

# install Bioconductor package installer
checkBiocManager_install()

# required BioConductor packages
required_BioconductorPackages<-c("flowCore","FlowSOM","ComplexHeatmap","limma","ConsensusClusterPlus")

# install
BioC_package_installer(required_BioconductorPackages)

# check if the packages can be loaded
lapply(required_BioconductorPackages, require, character.only = TRUE)
```
## Citations
S.Eswar, Z. T.Koenig, A. R.Tursi, J.Cobeña-Reyes, T.Tilburgs, and S.Andorf, “CytoBatchFlagR: A Comprehensive Framework to Objectively Assess High-Parameter Cytometry Data for Batch Effects,” Cytometry Part A (2026): 1–16, https://doi.org/10.1002/cyto.a.70024.
