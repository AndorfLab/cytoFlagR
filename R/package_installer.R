#!/usr/bin/env R

###############################################################################################
############################~~~~ cytoFlagR R package installer ~~~~############################
###############################################################################################

### this is a function to check for missing R packages and installing them
package_installer<-function(pckg) {
  for (i in pckg) {
    if(!require(i, character.only = TRUE, quietly = TRUE)) {
      install.packages(i, dependencies = TRUE)
      require(i, character.only = TRUE, quietly = TRUE)
    }
  }
}

### check if Bioconductor is installed
checkBiocManager_install<-function() {
  if(!require("BiocManager", character.only = TRUE, quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install()
}

### this is a function to check for missing Bioconductor packages and installing them
BioC_package_installer<-function(bc) {
  for (x in bc) {
    if(!require(x, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      require(x, character.only = TRUE, quietly = TRUE)
    }
  }
}