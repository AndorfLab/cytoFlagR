
########### ~~~~ Package installer ~~~~ ########### 

### Checks whether each package is installed, if not installs it and loads it
# pckg: the R package
package_installer <- function(pckg) {
  for (i in pckg) {
    if(!require(i, character.only = TRUE, quietly = TRUE)) {
      install.packages(i, dependencies = TRUE)
      require(i, character.only = TRUE, quietly = TRUE)
    }
  }
}

### Tries to load the BiocManager package. If not available, installs it
checkBiocManager_install <- function() {
  if(!require("BiocManager", character.only = TRUE, quietly = TRUE))
    install.packages("BiocManager")
 # BiocManager::install()
}

### Tries to load Bioconductor packages. If not available, installs it
# bc: the Bioconductor package
BioC_package_installer <- function(bc) {
  for (x in bc) {
    if(!require(x, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      require(x, character.only = TRUE, quietly = TRUE)
    }
  }
}
