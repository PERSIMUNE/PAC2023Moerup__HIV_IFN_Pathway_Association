standard.packages = c("SKAT", "OneR", "dplyr", "MASS", "magrittr", "devtools", "R.utils", "purrr", "broom", "tibble")
# bioclite.packages = c("mygene")
bioconductor.packages = c("mygene")
github.packages   = c()#c("Displayr/flipMultivariates")#c("marchtaylor/sinkr")


# NOTE: Assumes 

checkDependencies = function() {
  # Handle biocLite packages
  # source("https://bioconductor.org/biocLite.R")
  # for (bioclite in bioclite.packages) {
  #   if (!require(bioclite, character.only = TRUE, quietly = TRUE))
  #     biocLite(bioclite)
  # }  
  # # bioclite doesn't work for R 4.10. Need to change to biocmanager
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  for (packages in bioconductor.packages) {
    if (!suppressWarnings(suppressPackageStartupMessages({require(packages, character.only = TRUE, quietly = TRUE)})))
      BiocManager::install("mygene")
  }
  
  
  # Handle standard cran packages
  for (standard in standard.packages) {
    if (!suppressWarnings(suppressPackageStartupMessages({require(standard, character.only = TRUE, quietly = TRUE)})))
      install.packages(standard)
  }
  # Handle github packages
  require(devtools)
  for (github in github.packages) {
    if (!suppressWarnings(suppressPackageStartupMessages({require(github, character.only = TRUE, quietly = TRUE)})))
      install_github(github)
  }
}




checkDependencies()