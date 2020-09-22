#Universal Bioconductor package installation function
install.bioc <- function(pkg){
  vers <- getRversion()
  if (vers >= "3.6"){
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }else{
    if (!requireNamespace("BiocInstaller", quietly = TRUE)){
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg, suppressUpdates=TRUE)
    }else{
      BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
    }
  }
}

#Install Bioconductor dependencies
bioc_pkgs <- c("ComplexHeatmap")# ,"org.Mm.eg.db", "org.Rn.eg.db", "KEGG.db", "reactome.db", "GOSim")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.bioc(pkg)
    print("Installed!!!")
  }
}

#Install CRAN dependencies
#cran_pkgs <- c("ggplotify", "RColorBrewer", "reshape", "ggplot2", "shiny", "shinyjs", "tibble",
#               "gProfileR", "DT", "randomcoloR", "readxl", "cellranger", "devtools", "scales",
#               "gtools", "shinycssloaders", "shinyBS", "tidyverse",
#               "gridExtra", "gtable", "grid", "xlsx")
cran_pkgs <-c("vegan", "Hmisc", "data.table", "FunChisq", "psych", "PerformanceAnalytics", "gtools", "ggpubr", "gridExtra", "grid", "rcompanion", "plyr", "tidyr",
              "dplyr", "entropy", "tidyverse", "plotrix", "igraph", "network", "foreach", "doParallel", "stringr", "stringi", "ggcorrplot", "corrplot",
              "scales","ggplot2","reshape2","RelValAnalysis","RFLPtools","dendextend","ape","stringdist","reshape","rgl","factoextra","stats",
              "mclust","gplots","caret","png","ggplot2","circlize","cluster","pastecs","cowplot") #"ComplexHeatmap","doSNOW","doMC",

cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}

## INSTALL "ROWR" Package
## Package ‘rowr’ was removed from the CRAN repository.
## Formerly available versions can be obtained from the archive by clicking on this link:https://cran.r-project.org/src/contrib/Archive/rowr/
