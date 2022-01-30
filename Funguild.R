###### FUNGuild ###########

### set working directories 

setwd("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/") 
setwd("C:/Users/Lya/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/")

## Import packages - bioinformatics
#rem: phyloseq and metagenome seq: manually

packages = c("Rmisc", "gfcanalysis", "dplyr", "ggplot2", "gtable", "ggthemes", "data.table", "multcomp", "vegan", "biomformat", "ape", "lme4", "nlme", "car", "emmeans", "dendextend", "NMF", "psych", "randomForest", "ROCR", "caret", "tidyverse", "viridis", "reshape","reshape2", "ade4", "PMA", "genefilter", "ggrepel", "lattice", "corrplot", "gridExtra", "stats", "DESeq2", "phytools", "rfUtilities", "e1071", "imputeTS", "ggtree", 'multcopView', 'phyloseq')#use this function to check if each package is on the local machine. if a package is installed, it will be loaded. if any are not, the missing package(s) will be installed and loaded

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})





summary(funguild.guilds)
funguild.guilds$Trophic.Mode

 