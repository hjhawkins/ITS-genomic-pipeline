#########################################################################
########################  MICROBIOME DATA ANALYSIS   ####################
#########################################################################

### installation of bioconductor (Biocmanager) and install packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
# last time didn't work => i had to install manually :
#BiocManager::install("S4Vectors")
#BiocManager::install("IRanges")
#BiocManager::install("Rhdf5lib")
BiocManager::install("GenomicRanges")
BiocManager::valid("phyloseq")     
BiocManager::install("biomformat")
BiocManager::install("metagenomeSeq", ref = "RELEASE_3_9") # need a specific version
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("qiime2R")
BiocManager::install("ggtree")
BiocManager::install("PMA")

#rem: a specific version of metagenomeseq needed:you will need the devtools package for that. didn't work 
library(devtools)
install_github("HCBravoLab/metagenomeSeq", ref = "RELEASE_3_9") 


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



load("phyf2.RData")


#verify they are loaded. this function gives the attached packages.
search()
installed.packages()[, c("Package", "LibPath")]
packageVersion("vegan")
install.packages("vegan")
install.packages("devtools")
remove.packages(c("boot", "foreign", "nlme"))
install.packages(c("boot", "foreign", "nlme"))
Sys.getenv("R_LIBS_USER")

#if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
#p_load(Rmisc, gfcanalysis, dplyr, ggplot2, gtable, ggthemes, data.table, multcomp, phyloseq, vegan, biomformat, ape, metagenomeSeq, lme4, nlme, car, emmeans, dendextend, NMF, psych, randomForest, ROCR, caret, tidyverse, viridis, reshape, lattice, corrplot, gridExtra, stats, imputeTS)

library(phyloseq)

library(vegan)
library(biomformat)
library(qiime2R)
library("Biostrings") #for fast manipulation of large biological sequences or sets of sequences.
library(DESeq2) #:to produce a read count table, analyze count tables for differentially expressed genes, visualize the results, add extra gene annotations, and cluster samples and genes using transformed counts
library(ape) #for reading, writing, plotting, and manipulating phylogenetic trees, analyses of comparative data in a phylogenetic framework
library(metagenomeSeq)#differential abundance testing. Statistical analysis for sparse high-throughput sequencing

## Import packages - statistics

library(lme4)#Linear Mixed-Effects Models using 'Eigen' and S4 + GLMs
library(nlme) # for linear mixed effect models
#library(multcomp) #Simultaneous tests and confidence intervals for general linear hypotheses in parametric models, including linear, generalized linear, linear mixed effects, and survival models. 
library(car) #Companion to Applied Regression
library(emmeans) #Obtain estimated marginal means (EMMs) for many linear, generalized linear, and mixed models. 
library(dendextend) #Extending 'dendrogram'. letting you visualize and compare trees of 'hierarchical clusterings'. 
library(NMF) #Algorithms and Framework for Nonnegative Matrix Factorization (NMF)
library(psych)#corr.test. Functions are primarily for multivariate analysis and scale construction using factor analysis, principal component analysis, cluster analysis and reliability analysis, 
#library(matrixStats)#rowSds
library(randomForest) #Breiman and Cutler's Random Forests for Classification and Regression
library(ROCR) #Visualizing the Performance of Scoring Classifiers
library(caret) #Classification and Regression Training
library(e1071)#: R machine learning library.
#install.packages("remotes")
#remotes::install_github("gmteunisse/Fantaxtic")
library(stats)
library(phytools)
library(multcompView) #for the posthoc test on graphs
library(FactoMineR) # for PCA

## Import packages - data wrangling and plotting

#rlang package : install it then tydiverse works
library(tidyverse)
library(ggplot2)
library(viridis) #Default Color Maps from 'matplotlib'
library(reshape) #Flexibly restructure and aggregate data using just two functions: melt and cast
library(gridExtra)#for grid.arrange() to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.
library(lattice)#to improve on base R graphics by providing better defaults and the ability to easily display multivariate relationships.
library(corrplot)#for cor.mtest. A graphical display of a correlation matrix or general matrix.
library(plyr)
library(dplyr)
library("NMF")

library(imputeTS) #to replace na by something else
library(colortools) #to find nice combinations

## create an output directory, and parameters for export

outDir <- getwd() # Specify output directory
nmf.options(grid.patch=TRUE)#set to avoid blank first pdf page being created

## Load custom functions, from Katie lennard
source("C:/Users/VRMMAR024/Dropbox/R code/microbiome_custom_functions_Katie_v2.R")
#source("C:/Users/Lya/Dropbox/R code/microbiome_custom_functions_Katie.R")

##data sets in package phyloseq:

data(package= "phyloseq")
data(package="vegan")
data(varechem)


####################################################################
##########  import the outputs from QIIME2 pipeline  ################
####################################################################

### tax_table()   Taxonomy Table (need to be a matrix):  

taxon<-read.delim("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/metadata_taxon.tsv", comment.char="#")
tax2 <- taxon %>% separate(Taxon, sep=";", into = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"))
head(tax2)
row.names(tax2) <- tax2$Feature.ID
tax <- tax2 %>% dplyr::select(-Feature.ID, -Confidence)
row.names(tax) #show the row names => DNA sequences.
colnames(tax) # gives column names =>"Kingdom,Phylum,Class,Order,Family,Genus,species""Confidence"
tax$Kingdom <- gsub("k__", "", tax$Kingdom)
tax$Phylum <- gsub("p__", "", tax$Phylum)
tax$Class <- gsub("c__", "", tax$Class)
tax$Order <- gsub("o__", "", tax$Order)
tax$Family <- gsub("f__", "", tax$Family)
tax$Genus <- gsub("g__", "", tax$Genus)
tax$Species <- gsub("s__", "", tax$Species)
head(tax)
class(tax)
str(tax)
write.csv(tax, "tax.csv")
nrow(tax) # 95847 => 10 extra !!!


#test creating file for funguild
funguilds <- dplyr::left_join(Seqtab, taxon, by=c("OTU.ID" = "Feature.ID"))
funguild <- dplyr::rename(funguilds, taxonomy= Taxon ) %>% dplyr::select(-Confidence)
write.table(funguild, file = "funguild.txt", sep = "\t", row.names = F)
funguildtest <- dplyr::select(funguild, -(B11:U8)) 
Funguild_test <- funguildtest[1:100,]
write.table(Funguild_test, file = "Funguild_test.txt", sep = "\t", row.names = F)

#import results from funguild:

summary(funguild.guilds)
tax_fun <- dplyr::select(funguild.guilds, OTU.ID, Taxon, Taxon.Level, Trophic.Mode, Guild, Growth.Morphology, Trait, Confidence.Ranking, Notes)
head(tax_fun)
row.names(tax_fun) <- tax_fun$OTU.ID
tax_fun <- dplyr::select(tax_fun, -OTU.ID)
str(tax_fun)
Tax_total <- merge(tax, tax_fun, by = "row.names", all.x= TRUE) #means I want to keep all the columns of the 1st dataframe (x)
str(Tax_tot) #2972 obs. of  15 variables: merge(tax, tax_fun, by = "row.names")
str(Tax_total) #data.frame':	95847 obs. of  16 variables:
rownames(Tax_total) <- Tax_total$Row.names
Tax_total <- dplyr::select(Tax_total, -Row.names)


### otu_table()   OTU/ASV Table (need to be a matrix)

Seqtab<-read.delim("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/feature-tableforR.txt", header = TRUE)
head(Seqtab)
rownames(Seqtab) <- Seqtab$OTU.ID
Seqtab <- dplyr::select(Seqtab, -OTU.ID)
tseqtab <- t(Seqtab)
head(tseqtab)
str(tseqtab)
tseqtable <- as.data.frame(tseqtab)


#for funguild
Seqtab_funguild[1:5, 1:5]
Seqtab_funguild <- rename(Seqtab, OTU.ID = "taxonomy")
rownames(Seqtab_funguild)
tseqtab[1:5, 1:5]
tseqtab <- tseqtab[-1,] #now correct matrix, with sample names as row names, and sequences as column names
#

row.names(tseqtab) #show the row names => samples names.
colnames(tseqtab) # gives column names =>DNA sequence
tseqtab[5:5]
dim(tseqtab) #210 95837
nrow(tseqtab) # 210 samples
ncol(tseqtab) #95837 sequences
class(tseqtab) <- "numeric"
tseqtab[1:5, 1:5]
# make it crash : OTU <-read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/OTU.csv")

write.table(Seqtab_funguild, file = "Seqtab_funguild.txt")

### sample_data() Sample Data 

#ITS_metadata <- read.delim("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/ITS_metadata_tot.txt", header = T)
Env <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_tot.csv", header = T)
row.names(Env) <- Env$Sample.ID#define the sample names as row names
str(Env)
Env$Av_Na_mg_kg <- as.numeric(Env$Av_Na_mg_kg)
Env$Av_P_mg_kg <- as.numeric(Env$Av_P_mg_kg)
Env$Av_K_mg_kg <- as.numeric(Env$Av_K_mg_kg)
Env$V_C_. <- as.numeric(Env$V_C_.)
Env$V_dC13_12 <- as.numeric(Env$V_dC13_12)
Env$V_C.N <- as.numeric(Env$V_C.N)
Env$L_C_. <- as.numeric(Env$L_C_.)
Env$L_dC13_12 <- as.numeric(Env$L_dC13_12)
Env$L_C.N <- as.numeric(Env$L_C.N)
Env$Total_SR <- as.numeric(Env$Total_SR)
Env$spAbundance <- as.numeric(Env$spAbundance)
Env$Grass_SR <- as.numeric(Env$Grass_SR)
Env$Forb_SR <- as.numeric(Env$Forb_SR)
Env$Tree_SR <- as.numeric(Env$Tree_SR)
Env$Tree_10m_SR <- as.numeric(Env$Tree_10m_SR)
Env$X.bareC <- as.numeric(Env$X.bareC)
Env$X.litterC <- as.numeric(Env$X.litterC)
Env$X.FoliarC <- as.numeric(Env$X.FoliarC)

EnvB <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_B.csv", header = T)
row.names(EnvB) <- EnvB$Sample.ID#define the sample names as row names
str(EnvB)
EnvB$Av_Na_mg_kg <- as.numeric(EnvB$Av_Na_mg_kg)
EnvB$Av_P_mg_kg <- as.numeric(EnvB$Av_P_mg_kg)
EnvB$Av_K_mg_kg <- as.numeric(EnvB$Av_K_mg_kg)


EnvU <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_U.csv", header = T)
row.names(EnvU) <- EnvU$Sample.ID#define the sample names as row names
str(EnvU)
head(EnvU)
EnvU$Av_Na_mg_kg <- as.numeric(EnvU$Av_Na_mg_kg)
EnvU$Av_P_mg_kg <- as.numeric(EnvU$Av_P_mg_kg)
EnvU$Av_K_mg_kg <- as.numeric(EnvU$Av_K_mg_kg)
EnvU$Total_SR <- as.numeric(EnvU$Total_SR)
EnvU$spAbundance <- as.numeric(EnvU$spAbundance)
EnvU$X.bareC <- as.numeric(EnvU$X.bareC)
EnvU$X.litterC <- as.numeric(EnvU$X.litterC)
EnvU$X.BasalC <- as.numeric(EnvU$X.BasalC)

EnvN <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_N.csv", header = T)
row.names(EnvN) <- EnvN$Sample.ID#define the sample names as row names
str(EnvN)
EnvN$Av_Na_mg_kg <- as.numeric(EnvN$Av_Na_mg_kg)
EnvN$Av_P_mg_kg <- as.numeric(EnvN$Av_P_mg_kg)
EnvN$Av_K_mg_kg <- as.numeric(EnvN$Av_K_mg_kg)
EnvN$V_C_. <- as.numeric(EnvN$V_C_.)
EnvN$V_dC13_12 <- as.numeric(EnvN$V_dC13_12)
EnvN$V_C.N <- as.numeric(EnvN$V_C.N)
EnvN$L_C_. <- as.numeric(EnvN$L_C_.)
EnvN$L_dC13_12 <- as.numeric(EnvN$L_dC13_12)
EnvN$L_C.N <- as.numeric(EnvN$L_C.N)
EnvN$Total_SR <- as.numeric(EnvN$Total_SR)
EnvN$Grass_SR <- as.numeric(EnvN$Grass_SR)
EnvN$Forb_SR <- as.numeric(EnvN$Forb_SR)
EnvN$Tree_SR <- as.numeric(EnvN$Tree_SR)
EnvN$Tree_10m_SR <- as.numeric(EnvN$Tree_10m_SR)


EnvL <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_L.csv", header = T)
row.names(EnvL) <- EnvL$Sample.ID#define the sample names as row names
str(EnvL)
EnvL$Av_Na_mg_kg <- as.numeric(EnvL$Av_Na_mg_kg)
EnvL$Av_P_mg_kg <- as.numeric(EnvL$Av_P_mg_kg)
EnvL$Av_K_mg_kg <- as.numeric(EnvL$Av_K_mg_kg)
EnvL$Total_SR <- as.numeric(EnvL$Total_SR)
EnvL$Grass_SR <- as.numeric(EnvL$Grass_SR)
EnvL$Forb_SR <- as.numeric(EnvL$Forb_SR)
EnvL$Tree_SR <- as.numeric(EnvL$Tree_SR)
EnvL$Tree_10m_SR <- as.numeric(EnvL$Tree_10m_SR)
EnvL$X.BasalC  <- as.numeric(EnvL$X.BasalC )
EnvL$X.FoliarC <- as.numeric(EnvL$X.FoliarC)

EnvLC <- filter(EnvL, Zone =="C") %>% dplyr::select(-"site_zone")
row.names(EnvLC) <- EnvLC$Sample.ID#define the sample names as row names
EnvLR <- filter(EnvL, Zone =="R") %>% dplyr::select(-"site_zone")
row.names(EnvLR) <- EnvLR$Sample.ID#define the sample names as row names

# DO NOT DO ! sample reference in the table is always helpful. (remove the column sample names as it is now row name:)
#Env2 <- Env %>% select(-DNA_fullref)
# interesting : to add some columns that will be used afterwards. for example, divide a continuous variable sinto categories.
#sample_data(phy)$pH_binned <- cut(sample_data(phy)$pH, breaks = c(3, 4, 5, 6, 7, 8))


### phy_tree()    Phylogenetic Tree

rtree <- read.newick("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/rooted_tree.nwk") #package phytools
rtree #Phylogenetic tree with 95847 tips and 93983 internal nodes. rooted
str(rtree)
class(rtree) #"phylo"
head(rtree$tip.label)


ASV_table <- import_biom("C:/Users/Lya/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/feature-table.biom")


####################################################################
########### PHYLOSEQ: creation of a phyloseq object ################
####################################################################


# check that the sample names and features match in all cases
gplots::venn(list(metadata=rownames(Env), featuretable=colnames(tseqtab)))
gplots::venn(list(metadata=rownames(Env), featuretable=rownames(tseqtab))) # OK 210 in commom
length(intersect(rownames(Env),rownames(tseqtab)))#210 
gplots::venn(list(seqtab=colnames(tseqtab), tax=rownames(tax))) # be careful : 95837 on common, 10 extra in the tax table
gplots::venn(list(seqtab=colnames(tseqtable), tax=rownames(Tax_total)))  # be careful : 95837 on common, 10 extra in the tax table
colnames(tseqtable)
rownames(Tax_total)


length(intersect(colnames(tseqtab), rownames(Tax_tot)))
length(intersect(rownames(Env),sample_names(phy)))#210 (check that the sample names match in all cases)

setdiff(rownames(tax),colnames(tseqtab))
#[1] "41b7630d36f724d01f2b26e59a4ac197" "fae02eabe5967682bea5355ef976f79b" "a7c70b0b5b068473977d86ec50f6539f" "b3abc83dd70477740b0a7286cda17875"
#[5] "19f6abb1164acf08fa1531c58883f1f6" "ee88cdd1e1dcd592bf7c2072f7fe78b3" "539e7df1cacf6d82adedcc1ebf73711a" "c7cbe8863dc468e265d1b2c2517134b4"
#[9] "c94c8380f3d7973ad802381845e51610" "ba99628fbf767f2736722d3b62792182"

identical(rownames(Env), rownames(tseqtab)) #true
all(rownames(tseqtab) %in% rownames(Env)) #true
rownames(Env)
rownames(tseqtab)

tseqtab[1:2, 1:2]

# keep only overlapping features
common.ids <- intersect(colnames(tseqtab), rownames(tax))
tseqtabCID <- tseqtab[,common.ids]
taxCID <- tax[common.ids,]

#same but with funguild tax
common.ids <- intersect(colnames(tseqtable), rownames(Tax_total))
tseqtabCIDfun <- tseqtable[,common.ids]
taxCIDfun <- tax[common.ids,]


# rename ASVs from sequences to IDs:

TAX = tax_table(as.data.frame(taxCID) %>% as.matrix()) #moving the taxonomy to the way phyloseq wants it
TAXfun = tax_table(as.data.frame(taxCIDfun) %>% as.matrix())

#flip table so that taxa are rows
ASVs <- t(tseqtabCID)
identical(rownames(ASVs),rownames(taxCID)) #true
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
#Named vector:
names(seqs) <- ASV.IDs
head(seqs)
seq_lens <- nchar(seqs) #takes a character vector as an argument and returns a vector whose elements contain the sizes of the corresponding elements of x
seq_lens
plot(density(seq_lens))

#funguild
ASVs <- t(tseqtabCID)
identical(rownames(ASVs),rownames(taxCID)) #true
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
#Named vector:
names(seqs) <- ASV.IDs
head(seqs)
seq_lens <- nchar(seqs) #takes a character vector as an argument and returns a vector whose elements contain the sizes of the corresponding elements of x
seq_lens
plot(density(seq_lens))


OTU = otu_table(ASVs, taxa_are_rows = TRUE)

samples = sample_data(Env)

phy <- phyloseq(OTU,TAX, samples)
identical(taxa_names(phy),rownames(ASVs)) #TRUE
taxa_names(phy) <- names(seqs)



save(phy, file = paste0(outDir,"/phy.RData")) 


############################## To change to names of the tip of the tree:

str(phy)
tax_df <- data.frame(tax_table(phy))
tax_df$seq_len <- seq_lens
head(tax_df)
dim(tax_df) #95837     8
length(which(is.na(tax_df[,"Species"]) | tax_df[,"Species"] == "s__")) # 69628

length(intersect(rtree$tip.label,taxa_names(phy))) #=0 ! => I need to rename the tips of the tree
length(intersect(rtree$tip.label, seqs)) # 95837. we changed the names to user friendly=> rename tip labels accordingly

length(taxa_names(phy)) #95837
length(seqs) #95837
length(rtree$tip.label) #95847 (10 extra !!!)

rtreepruned <- drop.tip(rtree, c("41b7630d36f724d01f2b26e59a4ac197", "fae02eabe5967682bea5355ef976f79b", "a7c70b0b5b068473977d86ec50f6539f", "b3abc83dd70477740b0a7286cda17875", "19f6abb1164acf08fa1531c58883f1f6", "ee88cdd1e1dcd592bf7c2072f7fe78b3", "539e7df1cacf6d82adedcc1ebf73711a", "c7cbe8863dc468e265d1b2c2517134b4", "c94c8380f3d7973ad802381845e51610", "ba99628fbf767f2736722d3b62792182"))
length(rtreepruned$tip.label) #95837 

#Rename GTR_Tree tip labels with ASV.IDs
for_tree.IDs <- seqs[match(rtreepruned$tip.label,seqs)]
head(for_tree.IDs)
#sanity check:
head(match(rtreepruned$tip.label,seqs)) # [1] 69272 94282 93457  9596 53568 22101
identical(as.character(unlist(seqs[69272])),rtreepruned$tip.label[1])#TRUE => means that the first tip label correspond to the 69272 sequence.
rtreepruned$tip.label <- names(for_tree.IDs)

rtreephy <- merge_phyloseq(phy,rtreepruned)#Here we add the tree to the phyloseq object's tree slot
rtreephy
save(rtreephy, file = paste0(outDir,"/rtreephy.RData")) # Save phyloseq object: annotated object as a .RData object



######################################################################
###################### PRE-PROCESSING ################################
############ = to detect taxa/samples with very few features #########
######################################################################

rank_names(phy) # show the different PCOFGS
colnames(tax_table(phy)) #gives PCOFGS
sample_names(phy) # list samples
length(sample_names(phy))# 210
sample_variables(phy) # list the variables in the environmental table
taxa_names(phy)[1:10] # show the 10 first DNA sequences
otu_table(phy)[1:5, 1:5] # show the rows 1 to 5 (counts of sequences, for the first 5 samples) and column 1 to 5 (5 dna sequences)
tax_table(phy)[1:5, 1:4]
sample_data(phy)

######################################################################
########################## FILTER SAMPLES ############################

## (1) check the number of sequences/reads per sample => SEQUENCING DEPTH

nsamples(phy) #210
sum(sample_sums(phy) ==0) # no sample with no count
seqps <- data.frame(sample_sums(phy)) #gives the number of sequences (reads) per sample. 

#other option:
reads <- sample_sums(phy) #number of reads per sample
reads
length(which(reads<5000)) # 0 saples have less than 5000 reads
min(reads) # = L25 : 14845 reads. 
sort(reads)
summary(reads)
#Min.   1st Qu.  Median  Mean   3rd Qu.    Max. 
#14845   80746  113670  147727  178719  695575 
 


plot(seqps$sample_sums.phy.)
boxplot(seqps$sample_sums.phy.) # 11 outliers
pSeqDepth <- ggplot(seqps, aes(sample_sums.phy.)) + geom_histogram() + ggtitle("Sequencing Depth") + theme_bw() + xlab("Total nb reads") #plot the sequencing depth
pSeqDepth

# remove the samples that were merged.
phyf <- subset_samples(phy, Sample.ID != "B37" & Sample.ID != "B38" & Sample.ID != "B39" & Sample.ID != "B40" & Sample.ID != "U41" & Sample.ID != "U42" & Sample.ID != "N91" & Sample.ID != "N92" & Sample.ID != "N93" & Sample.ID != "N94" & Sample.ID != "N95" & Sample.ID != "L97" & Sample.ID != "L98" & Sample.ID != "L99" & Sample.ID != "L100" & Sample.ID != "L101" & Sample.ID != "L102" & Sample.ID != "L103")
seqpsf <- data.frame(sample_sums(phyf))
pSeqDepthf <- ggplot(seqpsf, aes(sample_sums.phyf.)) + geom_histogram() + ggtitle("Sequencing Depth") + theme_bw() + xlab("Total nb reads") #plot the sequencing depth
pSeqDepthf
plot(seqpsf$sample_sums.phyf.)
boxplot(seqpsf$sample_sums.phyf.)

# add sequencing depth to sample table
Enviseqdepth = data.table(as(sample_data(phyf), "data.frame"), TotalReads = sample_sums(phyf), keep.rownames = TRUE)

save(phyf, file = paste0(outDir,"/phyf.RData")) 

rtreephyf <- subset_samples(rtreephy, Sample.ID != "B37" & Sample.ID != "B38" & Sample.ID != "B39" & Sample.ID != "B40" & Sample.ID != "U41" & Sample.ID != "U42" & Sample.ID != "N91" & Sample.ID != "N92" & Sample.ID != "N93" & Sample.ID != "N94" & Sample.ID != "N95" & Sample.ID != "L97" & Sample.ID != "L98" & Sample.ID != "L99" & Sample.ID != "L100" & Sample.ID != "L101" & Sample.ID != "L102" & Sample.ID != "L103")

save(rtreephyf, file = paste0(outDir,"/rtreephyf.RData")) 


#######################################################################
######################## FILTER TAXA: #################################

#### (1) check the TOTAL COUNT per TAXA

#calculate total number of sequences/reads per taxa and deal with 0, singletons (sequence that appear only once in one sample) and doubleton (sequence that appear twice).

ntaxa(phyf) # 95837 tax
nsamples(phyf) #192
sum(taxa_sums(phyf) ==0) # 8839
sum(taxa_sums(phyf) ==1) # 0
sum(taxa_sums(phyf) ==2) # 8028
sum(taxa_sums(phyf) ==3) # 6683

TS <- taxa_sums(phyf)
plot(TS)
TSdf <- data.frame(TS)
TSdf2 <- cbind(rownames(TSdf), data.frame(TSdf, row.names=NULL)) # = to create a data frame with only 2 columns
TSdf3 <- arrange(TSdf2, desc(TS))
head(TSdf3)#max read per taxa = 476525
taxtotreads = data.table(tax_table(phyf), TotalCounts = taxa_sums(phyf), OTU = taxa_names(phyf)) #= to add sum per taxa in the taxtable and the nb of times 
head(taxtotreads)
# idem : tax2 <- data.frame(tax_table(phy), TotalCounts = taxa_sums(phy), ASV = taxa_names(phy)) # add the total count + the sequence in the tax table
ggplot(taxtotreads, aes(TotalCounts)) + geom_histogram() + ggtitle("Histogram of Total Counts")

taxtotreads[(TotalCounts <= 0),] # = the ASV with no observation
taxtotreads[(TotalCounts <= 1),] # = the ASV with singletons or less

# taxa cumulative sum
taxcumsum = taxtotreads[, .N, by = TotalCounts] # = to show the nb of times each sum per taxa happened (8000 times 0, ... times 2, etc.) . ".N" is a special built-in variable that holds the number of observations in the current group (nb de lignes, by group). command means they look at the column totalCounts, and group the same nb together, counting how many times it appears.

setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum
pCumSum + xlim(0, 100)# Zoom-in

#------------------------------------
phyf2 <- filter_taxa(phyf, function (x) {sum(x) > 0}, prune=TRUE)
# other way : GP <- prune_taxa(taxa_sums(phy) > 0, phy)

ntaxa(phyf2) # 86998 taxa
sum(taxa_sums(phyf2) ==0) # 0
save(phyf2, file = paste0(outDir,"/phyf2.RData")) 

rtreephyf2 <- filter_taxa(rtreephyf, function (x) {sum(x) > 0}, prune=TRUE)
#--------------------------------------


#### (2) check the TOTAL COUNT per PHYLUM and Remove OTUs that don't even have Phylum-level annotation (<NA> and "unidentified") and phyla for which only one feature was observed may also be worth filtering


t= which(is.na(tax_table(phyf2)[,"Phylum"]))
t
class(t) # integer
length(t) # 48111


#-------------------------
# create a table with the nb of feature for each phyla. only 1 : Zoopagomycota. Unidentified : 8567. <NA> : 48111
table(tax_table(phyf2)[, "Phylum"], exclude = NULL)  
# remove the NA and features with ambiguous phylum annotation are also removed. REM : check the database, sometimes "unidentified" sometimes "uncharacterized")
phyf3 <- subset_taxa(phyf2, !is.na(Phylum) & !Phylum %in% c("", "unidentified"))
table(tax_table(phyf3)[, "Phylum"], exclude = NULL)  
ntaxa(phyf3) # 30320 taxa
length(sample_names(phyf3)) #192 samples
save(phyf3, file = paste0(outDir,"/phyf3.RData")) 

rtreephyf3 <- subset_taxa(rtreephyf2, !is.na(Phylum) & !Phylum %in% c("", "unidentified"))
save(rtreephyf3, file = paste0(outDir,"/rtreephyf3.RData")) 
#--------------------------



### (3) Check the PREVALENCE


# (1)Lets generate a prevelance table (number of samples where the taxa appeared) for each taxa (2) Add taxonomy and total read counts to this data.frame

prevdf = apply(X = otu_table(phyf3),
               MARGIN = ifelse(taxa_are_rows(phyf3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phyf3),
                    tax_table(phyf3))
prevdf[1:10,] # a lot of taxa observed only twice in only one sample


#(2) SUPERVISED FILTERING : Compute the total and average prevalences of the features in each phylum and check if there are phyla that are comprised of mostly low-prevalence features.ddply = split-apply-combine strategy (a function from the package plyr). prevds= input."phylum" = how to split the input ( can be one or more factor). function = what to do with it.

plyr::ddply(prevdf, "Phylum", function(x){ data.frame(mean_prevalence=mean(x$Prevalence),sum_prevalence=sum(x$Prevalence), total_abundance=sum(x$TotalAbundance,na.rm = T),stringsAsFactors = F)})

#                   Phylum mean_prevalence sum_prevalence total_abundance
#1         Aphelidiomycota        4.000000             28             311
#2              Ascomycota        4.657726          85786        13033999
#3       Basidiobolomycota        1.142857              8             295
#4           Basidiomycota        2.914311          29385         4372529
#5      Blastocladiomycota        1.125000              9             406
#6  Calcarisporiellomycota        5.333333            608           25998
#7         Chytridiomycota        3.111913            862           42838
#8     Entomophthoromycota        1.250000             40             367
#9        Entorrhizomycota        2.387097             74            3321
#10          Glomeromycota        2.090308            949           22852
#11        Kickxellomycota        2.609375            501           12157
#12      Mortierellomycota        6.291765           2674          589761
#13           Mucoromycota        3.940828            666           40434
#14  Neocallimastigomycota        1.000000              3              15
#15          Rozellomycota        2.878788            285           17652
#16          Zoopagomycota        1.000000              1               4

## Define phyla to filter. I start with the one observed in only one sample.

filterPhyla = c("Zoopagomycota")
phyf4 = subset_taxa(phyf3, !Phylum %in% filterPhyla)
phyf4
save(phyf4, file = paste0(outDir,"/phyf4.RData")) 
rtreephyf4 = subset_taxa(rtreephyf3, !Phylum %in% filterPhyla)
save(rtreephyf4, file = paste0(outDir,"/rtreephyf4.RData")) 


###################################################
#(3) UNSUPERVISED PREVALENCE FILTERING: (without using the taxonomy). First, explore the relationship of prevalence and total read count for each feature. Subset to the remaining phyla

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phyf4, "Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phyf4),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme_bw()

# Define prevalence threshold as 5% of total samples, then Execute prevalence filter, using `prune_taxa()` function

prevalenceThreshold = 0.05 * nsamples(phyf4)
prevalenceThreshold # 9.6
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phyf5a = prune_taxa(keepTaxa, phyf4)
ntaxa(phyf5a) #  2718
save(phyf5a, file = paste0(outDir,"/phyf5a.RData")) 

# Method 2: prevalence filtering KATIE. filter by number of samples in which a taxa appears at least once (number of samples in which each OTU was non-zero). filter out taxa that are only present at very low numbers in a small minority of samples. The filter below retains only ASVs that are present at at least 10 counts at least 10% of samples OR that have a total relative abundance of > 0.001 of the total number of reads/sample

total = median(sample_sums(phyf4))#find median sample read count
phyf5b = filter_taxa(phyf4,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE) #from katie: M.std (here phy.filtered) = standardized to the median
ntaxa(phyf5b) # (katie's name : M.f)  we have retained 10274 from the 30319 taxa
save(phyf5b, file = paste0(outDir,"/phyf5b.RData")) 

rtreephyf5b = filter_taxa(rtreephyf4,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE) #from katie: M.std (here phy.filtered) = standardized to the median
ntaxa(rtreephyf5b) # (katie's name : M.f)  we have retained 10274 from the 30319 taxa
save(rtreephyf5b, file = paste0(outDir,"/rtreephyf5b.RData")) 

#-------------
prevdf5b1 = apply(X = otu_table(phyf5b),
               MARGIN = ifelse(taxa_are_rows(phyf5b), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf5b <- data.frame(Prevalence = prevdf5b1,
                     TotalAbundance = taxa_sums(phyf5b),
                     tax_table(phyf5b))


subprevdf5b = subset(prevdf5b, Phylum %in% get_taxa_unique(phyf5b, "Phylum"))

ggplot(subprevdf5b, aes(TotalAbundance, Prevalence / nsamples(phyf5b),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme_bw()

#gpsf = filter_taxa(M.std, function(x) sd(x)/mean(x) > 3.0, TRUE)#Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation



########################################################################
################# ABUNDANCE VALUE TRANSFORMATION #######################
########################################################################


## standardize sample read count so we can compare between sample (Normalize number of reads in each sample using median sequencing depth).

total = median(sample_sums(phyf5b))#find median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count(Standardize abundances to the median sequencing depth)
phyMRA = transform_sample_counts(phyf5b, standf)#apply to phyloseq object
sample_sums(phyMRA)
save(phyMRA, file = paste0(outDir,"/phyMRA.RData")) 

ttotal = median(sample_sums(rtreephyf5b))#find median sample read count
tstandf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count(Standardize abundances to the median sequencing depth)
rtreephyMRA = transform_sample_counts(rtreephyf5b, tstandf)#apply to phyloseq object
sample_sums(rtreephyMRA)
save(rtreephyMRA, file = paste0(outDir,"/rtreephyMRA.RData")) 


## transform to relative abundance, creating the new object, which is then filtered such that only OTUs with a mean greater than 10^-5 are kept.

phyRA  = transform_sample_counts(phyf5b, function(x) x / sum(x) )
#rel_abund_filtered = filter_taxa(rel_abund, function(x) mean(x) > 1e-5, TRUE)
save(phyRA, file = paste0(outDir,"/phyRA.RData")) 
rtreephyRA  = transform_sample_counts(rtreephyf5b, function(x) x / sum(x) )
#rel_abund_filtered = filter_taxa(rel_abund, function(x) mean(x) > 1e-5, TRUE)
save(rtreephyRA, file = paste0(outDir,"/rtreephyRA.RData")) 


## log transform sample counts

phylog  = transform_sample_counts(phyf5b, function(x) log(1 + x) )
save(phylog, file = paste0(outDir,"/phylog.RData")) 
rtreephylog = transform_sample_counts(rtreephyf5b, function(x) log(1 + x) )
save(rtreephylog, file = paste0(outDir,"/rtreephylog.RData")) 

#make an hist to check if it is sufficient for normalizing the abundance data
ss_df_log <- data.frame(sum = sample_sums(phyf5log)) # data frame with a column for the read counts of each sample
# Histogram of sample read counts
ggplot(ss_df_log, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 30) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
taxtotreadslog = data.table(tax_table(phyf5log), TotalCounts = taxa_sums(phyf5log), OTU = taxa_names(phyf5log)) #= to add sum per taxa in the taxtable and the nb of times 
ggplot(taxtotreadslog, aes(TotalCounts)) + geom_histogram() + ggtitle("Histogram of Total Counts")



# to do a hellinger transformation :

ps2 <- transform_sample_counts(ps1, function(x) sqrt(x / sum(x)))




#########################################################################
####################### SUBSETTING/agglomerate DATA FRAME ###############
#########################################################################

######### Per SAMPLE ######

# select only one part of samples

#nkuhlu <- subset_samples(M.std, site =="N") # select only samples from Nkuhlu
#nkuhlu


# regroup together samples from the same treatment / variable. by default, calculates the mean.

#nkuhlu_treatment <- merge_samples(nkuhlu, "Fhzone") 


######### Per TAXA #######

# Agglomerate taxa of the same type: 

phyglom <- tax_glom(phy, taxrank=rank_names(phy)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

phyphylum <- tax_glom(phyf5b, "Phylum", NArm = TRUE)
phyphylum
otu_table(phyphylum)
tax_table(phyphylum)


# create a table of abundance for each taxonomic rank (per phylum, per genus, etc). 

phy_phylum = factor(tax_table(phy)[, "Phylum"]) #Create a factor corresponding to the Phylums#
phy_phylum 
phy_phylum_table = t(apply(otu_table(phy), MARGIN=1, function(x) {tapply(x, INDEX = phy_phylum, FUN = sum, na.rm = TRUE, simplify = TRUE)}))

dim(phy_phylum_table)
?tapply

# create a table of abundance for each genus

phy_genus = factor(tax_table(phy)[, "Genus"]) #Create a factor corresponding to the genu column in the tax matrix
phy_genus_table = t(apply(otu_table(phy), MARGIN=1, function(x) {tapply(x, INDEX = phy_genus, FUN = sum, na.rm = TRUE, simplify = TRUE)}))
dim(phy_genus_table)

  
#keep only certain phylums/taxonomic ranks. method 1 = specify the ones to keep, method 2 the ones to remove

phy-filtered2 <- subset_taxa(phyl_filtered, Phylum %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", "Haptophyta", "Ochrophyta", "Cercozoa"))
phy-filtered2 <- subset_taxa(phyl_filtered, !(Class %in% c("Syndiniales", "Sarcomonadea")))

order1 = subset_taxa(phyl_filtered, Order == "Lactobacillales")
plot_abundance(order1, Facet = "Genus", Color = NULL)


# Katie: It's useful to collapse OTUs down to their lowest available taxonomic annotations:
#For this we use the tax_glom.kv() function from the 'microbiome_custom_functions.R' script loaded at the beginning of this script

taxglom_phy <- tax_glom.kv(phyf5b) #[1] "There are now 1416 merged taxa"
ntaxa(taxglom_phy) #1416
save(taxglom_phy, file = paste0(outDir,"/taxglom_phy.RData")) 

taxglom_rtreephy <- tax_glom.kv(rtreephyf5b) # rem: "Removing phylogenetic tree"

























#########################################################################
######################### PLOT THE RESULTS #############################
#########################################################################



#########################################################################
################### barplots tot abundance###############################


### basic plot = total abundance, for all samples. The dataset is plotted with every sample mapped individually to the horizontal (x) axis, and abundance values mapped to the veritcal (y) axis. At each sample's horizontal position, the abundance values for each OTU are stacked in order from greatest to least, separate by a thin horizontal line. 

plot_bar(phy)


### Add fill color to represent the Genus to which each OTU belongs.

plot_bar(phy, fill="Genus") #doesn't work (too big..)

plot_bar(Phy_f2, fill = "Phylum") 


### Group the samples together by the SampleType/treatment variable; abundance values for the same OTU from the same SampleType will be stacked as separate bar segments, and so the segment lines may not accurately portray the observed richness (because the same OTU might be shown more than once for the same horizontal axis grouping). The total stacked bar height at each horizontal position indicate the sum of all reads for that sample(s). There is not attempt by plot_bar to normalize or standardize your data, which is your job to do before attempting to interpret/compare these values between samples.

plot_bar(phy, x="Site.ttm.zone", fill="Genus")

plot_bar(phyglom, x = "FH", fill = "Phylum")


### More Sophisticated Organization using Facets. subset samples, and plot per site. You can also add extra layers, using ggplot2: p <- plot_bar() ; then do p + geom_point() for example

plot_bar(phy, x = "FH", fill = "Phylum", facet_grid = ~site)

plot_bar(phy, x= "FH", fill = "Phylum", facet_grid = Zone ~ site) #start 2:37 stop 3:08 was all black !

nkuhlu <- subset_samples(phy, site =="N")
plot_bar(nkuhlu, x= "FH", fill = "Phylum", facet_grid = Zone ~ site) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + labs(title="Total abundance per site and treatment", x="Treatment") 

plot_bar(phyglom, x= "FH", fill = "Phylum", facet_grid = Zone ~ site) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + labs(title="Total abundance per site and treatment", x="Treatment") # MUCH FASTER !!!!


#http://joey711.github.io/phyloseq-demo/Restroom-Biogeography
# 1. merge by category . 
#phy2 <- merge_samples(phy, "Site.ttm.zone")
#Repair the merged values associated with each surface after merge.
#sample_data(phy2)$Site.ttm.zone <- levels(sample_data(phy2)$Site.ttm.zone)

#2. transform to percentage of total available
#phy2 = transform_sample_counts(phy2, function(x) 100 * x/sum(x))
#if you want to stack horizontally: add + coord_flip()
title = "Figure 1 Part A (remake), attempt 2"
plot_bar(phy2, "Site.ttm.zone", fill = "family", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")



##########################################################################
################## barplot relative abundance#############################

# normaly, to make a stacked percent barplot, instead of "stacked",  position="fill"
## option 1 : bar stacked vertically. (1) melt to long format (for ggploting) prune out phyla below 2% in each sample

sample_data(phyf5b)$site_ttm_zone <- do.call(paste0, sample_data(phyf5b)[c("site", "Treatm", "Zone")], na.rm = TRUE)

sample_data(phyf5b) %>% unite("site_ttm_zone", c("site", "Treatm", "Zone"), remove = FALSE, na.rm = FALSE)
sample_data(phyf5b) <- sample_data(phyf5b)[,!"site_ttm_zone"]

head(sample_data(phyf5b))

phy_phylum <- phyf5b %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# (2) Set colors for plotting and plot
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

RA <- ggplot(phy_phylum, aes(x = FH, y = Abundance, fill = Phylum)) + 
  facet_grid(Zone~site) +
  geom_bar(stat = "identity", position="fill") + #I had to add "fill, because the fact that the amount of sample per treatment is not the same impacted the stack
  scale_fill_manual(values = phylum_colors) +
  #scale_x_discrete(
   # breaks = c("7/8", "8/4", "9/2", "10/6"),
    #labels = c("Jul", "Aug", "Sep", "Oct"), 
 #   drop = FALSE
 # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  xlab("Treatment") +
  ggtitle("Phylum Composition of fungal Communities by Sampling Site") +
  theme_bw()

RA



##option 2: katie: % of sequences in x axix, treatment in y axix, then barplot of relative abundance, sorted by order of importance, with the bar.plots function from the microbiome_custom_functions.R. M_phy = phyloseq object with collapsed OTUs down to their lowest available taxonomic annotations

barplot <- bar.plots(physeq = phyf5b,cat = "Site_ttm_zone",level = "Phylum", count = 500, perc = 0.25, outDir=outDir, filen = 'Barplots_by_treatment') + ggtitle("Phylum Composition of fungal Communities by Sampling Site") +xlab("")
print(barplot)

barplot2 <- bar.plots(physeq = phyf5b,cat = "Treatm",level = "Phylum", count = 500, perc = 0.25, outDir=outDir, filen = 'Barplots_by_treatment') + facet_grid(Zone ~ site)### doesn t work!
barplot2



## from callahan: Violin plot. defining a custom plot function, that uses phyloseq's psmelt()function to define a relative abundance graphic. We will use this to compare differences in scale and distribution of the abundance values in our phyloseq object before and after transformation.

plot_abundance = function(physeq,title = "", Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Site.ttm.zone",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

plotBefore = plot_abundance(phyl_filtered,"")
plotAfter = plot_abundance(ps3ra,"")
grid.arrange(nrow = 2, plotBefore, plotAfter)# Combine each plot into one graphic.



##########################################################################
###################### Plot the tree #####################################


#You can plto the whole tree, or froup the tips. (before hand, using tip_glom or tax_glom ,  or within the plot tree function)
#plot_tree(closedps, color="Treatment", size="abundance", sizebase=2, label.tips="taxa_names")

#x1 = merge_taxa(closedps, taxa_names(closedps)[3:27], 2)
#plot_tree(x1, color="Treatment", size="abundance", sizebase=2, label.tips="taxa_names")

#plot_tree(physeq2, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
#plot_tree(physeq2, color="Sample", ladderize="left") + coord_polar(theta="y")

phy.test <- merge_phyloseq(phy,GTR_Tree) #merge tree into phyloseq object
phy.test
plot_tree(phy.test)

plot_tree(phyglom) #by default, plot at each tip one dot for each sample this OTU was observed
plot_tree(phyglom, "treeonly") # just the branches
plot_tree(phyglom, "treeonly", ladderize = "left") # makes the way branches are organised better
plot_tree(phyglom, label.tips = "Phylum", ladderize = "left")
plot_tree(phyglom, color="site", ladderize="left", label.tips = "Phylum") + coord_polar(theta="y")
plot_tree(phyglom, label.tips = "Phylum") + coord_polar(theta="y")

plot_tree(phyglom, label.tips = "Phylum", color = "site", size = "abundance", plot.margin = 0.5, ladderize = "left", title = "By phylum")
plot_tree(phyglom, label.tips = "Phylum", color = "site", ladderize = TRUE, title = "By phylum")

plot_tree(phyglom, size = "Abundance", color = "site", justify = "yes please", ladderize = "left", label.tips = "Phylum") +
  scale_size_continuous(range = c(1, 3))


##########################################################################
################################## HEATMAP ###############################


# base plot: all samples in X axis and all OTU's in Y axis

plot_heatmap(rtreephy)

#The following two lines subset the dataset to just the top 300 most abundant Bacteria taxa across all samples (in this case, with no prior preprocessing. Not recommended, but quick). Then group samples per sample type.

Phy_bact <- subset_taxa(phy, Kingdom=="Bacteria") #not needed in our case: only bacteria
Phy_bact <- prune_taxa(names(sort(taxa_sums(Phy_bact),TRUE)[1:300]), Phy_bact)

plot_heatmap(phy, "NMDS", "bray", "FH", "Phylum") + facet_grid(Zone ~ site) #too big, takes a long time and doesn't work
plot_heatmap(phyglom, "NMDS", "bray", "Site.ttm.zone", "Phylum") 

N <- subset_samples(phyf2, site =="N")
Nphylum <- tax_glom(N, "Phylum", NArm = TRUE)
nkuhlu_treatment <- merge_samples(Nphylum, "FH") 
plot_heatmap(nkuhlu_treatment, "NMDS", "bray")

tot_treatm <- merge_samples(phyf2, "site_ttm_zone") 
TP <- tax_glom(tot_treatm, "Phylum", NArm = TRUE)
plot_heatmap(TP, "NMDS", "bray")


# you can also add a grouping per taxonomic rank

p <- plot_heatmap(phy, method = "NMDS", distance = "bray", "Site.ttm.zonee", "Phylum")
p$scales$scales[[1]]$name <- "My X-Axis" #if you wanted to change the axis labels, but not the labels on individual features?
p$scales$scales[[2]]$name <- "My Y-Axis"
print(p)

# If, for whatever reason, you need to change the default color scheme, it is possible through the low, high, and na.value arguments. 
# By default, the plot_heatmap color scale is a log transformation with base 4, using log_trans(4) from the scales package. 

plot_heatmap(phy, "NMDS", "bray", "SampleType", "Family", low="#000033", high="#CCFF66") # dark-blue to green
plot_heatmap(phy, "NMDS", "bray", "SampleType", "Family", low="#000033", high="#FF3300") # dark-blue to red scheme.


#you can also ry different ordination methods, distances
#plot_heatmap(phy, "PCoA", distance="bray", sample.label="...", taxa.label="Genus", low="#FFFFCC", high="#000033", na.value="white")


#It is better to only consider the most abundant OTUs for heatmaps. For example one can only take OTUs that represent at least 20% of reads in at least one sample. Remember we normalized all the sampples to median number of reads (total). We are left with only 33 OTUS which makes the reading much more easy.

total = median(sample_sums(nkuhlu))
standf = function(x, t=total) round(t * (x / sum(x)))
nkuhlu_norm = transform_sample_counts(nkuhlu, standf)
nkuhlu_abund <- filter_taxa(nkuhlu_norm, function(x) sum(x > total*0.20) > 0, TRUE)
plot_heatmap(nkuhlu_abund, method = "NMDS", distance = "bray")



#############################################
#katie # didn't work... Create sample pairwise distance matrix (Bray-Curtis distance) and cluster with hclust():

filename <- c("heatmap_merged_taxa")
main <- paste("Merged taxa, Bray-Curtis distance")
f = paste0(outDir,"/",filename,".pdf")

# Color specification for column annotations above heatmap:
D.cols = c("FH" = "#FFFF99", "NN" = "#FFFFCC", "NH" = "#FFFF00", "FN" = "#FFCC33")
colours = list(FH=D.cols)

# Create sample pairwise distance matrix (Bray-Curtis distance) and cluster with hclust():
# 1. agglomerate lowest avail taxa, 2.
M.f = filter_taxa(taxglom_phy,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)

set.seed(2)

diss <- distance(filterphyf2glom,method = "bray", type = "samples")
clust.res<-hclust(diss)
sample.order = clust.res$order

# Heatmap is output to file (the heatmap.k function can be found in the 'microbiome_custom_functions.R' script)
hm = heatmap.k(physeq= filterphyf2glom, annot.cols = 9, main = main,filename = f,colours=colours,cexRow = 0.5, Colv = sample.order,labrow = TRUE)  

print(hm)

sample_variables(filterphyf2glom)



######################## HEAT TREES ###################################


#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#library(taxa)
#obj <- parse_tax_data(otu_data,class_cols = "taxonomy", class_sep = ";") # The column in the input table# What each taxon is seperated by
#print(obj)





























##########################################################################
################### ALPHA DIVERSITY (within sample) ######################
##########################################################################

#you can prune OTUs that are not present in any of the samples - BUT DON'T TRIM MORE THAN THAT! many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.

### Plot with the plot_richness function from the phyloseq package

p <- plot_richness(phyf2,x = "FH",measures=c("Shannon"), title = "Alpha Diveristy - Shannon") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw() 
p

p <- plot_richness(phyf2,x = "FH",measures=c("Observed"), title = "Alpha Diveristy - Observed") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("Chao1"), title = "Alpha Diveristy - Chao1") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("ACE"), title = "Alpha Diveristy - ACE") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("Simpson"), title = "Alpha Diveristy - Simpson") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("InvSimpson"), title = "Alpha Diveristy - InvSimpson") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p <- plot_richness(phyf2,x = "FH",measures=c("Fisher"), title = "Alpha Diveristy - Fisher") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()


pdf(paste0(outDir,"/alpha_diversity_by_treatment.pdf"))
p
dev.off()


# Compute the Shannon diversity associated with each sample and join it with sample annotation.

alpha_Shannon <- estimate_richness(phyf2, split = TRUE, measure = "Shannon")
alpha_Shannon$Sample.ID <- rownames(alpha_Shannon) %>% as.factor()
ps_samp <- Env %>% left_join(alpha_Shannon, by = "Sample.ID") 
#Rem : probl = some samples were removed from sample table

#other way:
est <- estimate_richness(phyf2, split = TRUE, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) # gives "shannon" for all samples.
shann <- cbind(est,sample_data(phyf2)[,c("Sample.ID", "Sample_name", "site_zone", "fire", "Fire.freq", "herbiv", "FH", "Treatm", "rep", "site", "Zone", "site_FH_zone", "site_ttm_zone", "Block", "plot")]) #add a column "fire " and "herbiv"

write.csv(shann, file = "alpha_div.csv")
shann <- read.csv("shannon.csv")

alpha_div_B <- dplyr::filter(shann, site == "B") 
alpha_div_U <- dplyr::filter(shann, site == "U") 
alpha_div_LC <- dplyr::filter(shann, site_zone == "LC") 
alpha_div_LR <- dplyr::filter(shann, site_zone == "LR") 
alpha_div_NC <- dplyr::filter(shann, site_zone == "NC") 
alpha_div_NR <- dplyr::filter(shann, site_zone == "NR") 

alpha_B <- alpha_div_B %>% dplyr::group_by(Treatm) %>% dplyr::summarize(mean(Observed), SD_obs = sd(Observed), mean_shann = mean(Shannon), SD_shan = sd(Shannon))
alpha_U <- dplyr::filter(shann, site == "U") %>% dplyr::group_by(Treatm) %>% dplyr::summarize(mean(Observed), SD_obs = sd(Observed), mean_shann = mean(Shannon), SD_shan = sd(Shannon))
alpha_L <- dplyr::filter(shann, site == "L") %>% dplyr::group_by(site_FH_zone) %>% dplyr::summarize(mean(Observed), SD_obs = sd(Observed), mean_shann = mean(Shannon), SD_shan = sd(Shannon))
alpha_N <- dplyr::filter(shann, site == "N") %>% dplyr::group_by(site_FH_zone) %>% dplyr::summarize(mean(Observed), SD_obs = sd(Observed), mean_shann = mean(Shannon), SD_shan = sd(Shannon))


summaryalpha <- bind_rows(alpha_B, alpha_U, alpha_N, alpha_L)
write.csv(summaryalpha, file = "summary_alpha_div.csv")

aovShan <- aov(Shannon ~ site_zone, shann)
aovO <- aov(Observed ~ site_zone, shann)

aovB <- aov(Shannon ~ Treatm, alpha_div_B)
aovU <- aov(Shannon ~ Treatm, alpha_div_U)
aovLC <- aov(Shannon ~ Treatm, alpha_div_LC) 
aovLR <- aov(Shannon ~ Treatm, alpha_div_LR) 
aovNC <- aov(Shannon ~ Treatm, alpha_div_NC)  
aovNR <- aov(Shannon ~ Treatm, alpha_div_NR) 

summary(aovShan)
summary(aovO)
summary(aovB)
summary(aovU) 
summary(aovLC)  
summary(aovLR)
summary(aovNC)  
summary(aovNR)

plot(aovShan)
plot(aovB)
plot(aovU) 
plot(aovLC)  
plot(aovLR)
plot(aovNC)  
plot(aovNR)

aovBobs <- aov(Observed ~ Treatm, alpha_div_B)
aovUobs <- aov(Observed ~ Treatm, alpha_div_U)
aovLCobs <- aov(Observed ~ Treatm, alpha_div_LC) 
aovLRobs <- aov(Observed ~ Treatm, alpha_div_LR) 
aovNCobs <- aov(Observed ~ Treatm, alpha_div_NC)  
aovNRobs <- aov(Observed ~ Treatm, alpha_div_NR) 

summary(aovBobs)
summary(aovUobs) 
summary(aovLCobs)  
summary(aovLRobs)
summary(aovNCobs)  
summary(aovNRobs)

plot(aovBobs)
plot(aovUobs) 
plot(aovLCobs)  
plot(aovLRobs)
plot(aovNCobs)  
plot(aovNRobs)


##### Post HOc test: different options 

# from stat course:
phtB <- TukeyHSD(aovB, ordered=T) 
phtB
plot(TukeyHSD(aovB))
library(multcomp)
tukB <- glht(aovB, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
tukB
summary(tukB)          # standard display #Get p values for comparisons
tuk.cldB <- cld(tukB)   # letter-based display
#Mike: tuk.cld <- cld(summary(tuk), level = 0.05, decreasing = FALSE)#Assign letters for signficant comparisons
tuk.cldB
tuk_DFB <- data.frame(rownames(as.data.frame(tuk.cldB$mcletters$Letters)), as.data.frame(tuk.cldB$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFB) <- c("Treatm", "shann_Letter") #Give names to data frame

##################
op <- par(mai=c(1,1,1.5,1))# set top margin bigger
plot(tuk.cldB) #amazing boxplot with letters !!!!
# main="Multiple comparison results", xlab="", ylab="")
par <- par(op)
##################

#OVERALL

tukshan <- glht(aovShan, linfct = mcp(site_zone = "Tukey"))#Estimate multi-comp means
summary(tukshan)
tuk.cldshan <- cld(tukshan)   # letter-based display
tuk_DFshan <- data.frame(rownames(as.data.frame(tuk.cldshan$mcletters$Letters)), as.data.frame(tuk.cldshan$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFshan) <- c("site zone", "shann_Letter") #Give names to data frame
tukO <- glht(aovO, linfct = mcp(site_zone = "Tukey"))#Estimate multi-comp means
summary(tukO)
tuk.cldO<- cld(tukO)   # letter-based display
tuk_DFO <- data.frame(rownames(as.data.frame(tuk.cldO$mcletters$Letters)), as.data.frame(tuk.cldO$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFO) <- c("site zone", "observed") #Give names to data frame


#UKULINGA

tukU <- glht(aovUobs, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
summary(tukU)
tuk.cldU <- cld(tukU)   # letter-based display
tuk_DFU <- data.frame(rownames(as.data.frame(tuk.cldU$mcletters$Letters)), as.data.frame(tuk.cldU$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFU) <- c("Treatm", "shann_Letter") #Give names to data frame

#LETABA

tukLC <- glht(aovLC, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
tuk.cldLC <- cld(tukLC)   # letter-based display
tuk_DFLC <- data.frame(rownames(as.data.frame(tuk.cldLC$mcletters$Letters)), as.data.frame(tuk.cldLC$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFLC) <- c("Treatm", "shann_Letter") #Give names to data frame

tukLR <- glht(aovLRobs, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
tuk.cldLR <- cld(tukLR)   # letter-based display
tuk_DFLR <- data.frame(rownames(as.data.frame(tuk.cldLR$mcletters$Letters)), as.data.frame(tuk.cldLR$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFLR) <- c("Treatm", "shann_Letter") #Give names to data frame

#NKUHLU

tukNC <- glht(aovNC, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
tuk.cldNC <- cld(tukNC)   # letter-based display
tuk_DFNC <- data.frame(rownames(as.data.frame(tuk.cldNC$mcletters$Letters)), as.data.frame(tuk.cldNC$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFNC) <- c("Treatm", "shann_Letter") #Give names to data frame

tukNR <- glht(aovNR, linfct = mcp(Treatm = "Tukey"))#Estimate multi-comp means
tuk.cldNR <- cld(tukNR)   # letter-based display
tuk_DFNR <- data.frame(rownames(as.data.frame(tuk.cldNR$mcletters$Letters)), as.data.frame(tuk.cldNR$mcletters$Letters))#Copy the letters into a data.frame
names(tuk_DFNR) <- c("Treatm", "shann_Letter") #Give names to data frame




#library(agricolae)
#HSD.test(aovB, "tx", group=TRUE)
#also: function multcompLetters 
# my first attempt. also gives the letters
posthocB <- lsmeans(aovB, pairwise~Treatm, adjust="Tukey")
CLD(posthocB, Letters=letters, alpha = .05) #but they say it is depreciated




# plot the calculated shannon

alpha_boxplot <- ggplot(ps_samp, aes(FH, Shannon)) + geom_boxplot()+facet_grid(site ~ Zone)
alpha_boxplot

alpha_boxplot_savanna <- ggplot(subset(ps_samp, site == "N" | site == "L"), aes(FH, Shannon)) +geom_boxplot()+facet_grid(site ~ Zone)
alpha_boxplot_savanna
alpha_boxplot_grassland <- ggplot(subset(ps_samp, site == "U" | site == "B"), aes(FH, Shannon)) +geom_boxplot()+facet_grid(. ~ site)
alpha_boxplot_grassland


# make a table with the average alpha diversity per site, zone and treatment

diversity_means <- ps_samp %>%
  group_by(site, Zone, Site.ttm.zone) %>%
  summarise(mean_div = mean(Shannon)) %>%
  arrange(mean_div)


### test the normality : do an histogram of shannon measures + Kolmogorov-Smirnov test + Shapiro-Wilk test. rem: need to it per group we test, not for the whole dataset.

h<- ggplot(shann, aes(x = Shannon)) + geom_histogram() + theme_bw()
h

ks.test(shann$Shannon, "pnorm", mean = mean(shann$Shannon), sd=sd(shann$Shannon)) #The null hypothesis of the K-S test is that the distribution is normal.Therefore, if p-value of the test is >0.05, we do not reject the null hypothesis and conclude that the distribution in question is not statistically different from a normal distribution.
#D = 0.088895, p-value = 0.09619
#alternative hypothesis: two-sided

shapiro.test(shann$Shannon) # tests the null hypothesis is that the population is normally distributed.have greater power when compared to the K-S test.
#W = 0.92997, p-value = 5.6e-08 => a lot smaller than 0.05 => significantly differ from a normal distribution.


###	Kruskal-Wallis rank sum test : tot test the impact of fire and herbivores.

kruskal.test(alpha_div_B$Shannon~alpha_div_B$fire)

t_fire <- kruskal.test(shann$Shannon~shann$fire)
t_fire # Kruskal-Wallis chi-squared = 0.10579, df = 1, p-value = 0.745

t_herbiv <- kruskal.test(shann$Shannon~shann$herbiv)
t_herbiv # Kruskal-Wallis chi-squared = 0.076529, df = 1, p-value = 0.7821

### ANOVA
d <- data.frame(lev=shann$site_FH_zone, y=shann$Shannon)

aov <- aov(shann$Shannon~shann$site_FH_zone, data=shann)
summary(aov)
#                       Df Sum Sq Mean Sq F value     Pr(>F)    
#shann$site_FH_zone     22  19.62  0.8920   5.773 6.45e-12 ***
#Residuals             169  26.11  0.1545                     
#That tells us that there is a difference, but does not tell us which means are different. A Tukey's Honest Significant Difference (HSD) test can do pairwise comparisons of the means to find this out. We will use the HSD.test function from the agricolae package since it provides grouping codes that are useful for graphing.https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#rem:try this for the letters!

#other way: to try do a post hoc
TUKEY <- TukeyHSD(aov, ordered = FALSE, conf.level = 0.95) # add "shann$site_FH_zone" ?
print(TUKEY)
generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "shann$site_FH_zone")
# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)
# Draw the basic boxplot
a <- boxplot(shann$Shannon~shann$site_FH_zone , ylim=c(min(data$value) , 1.1*max(data$value)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="value" , main="")
# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )
#Add the labels
text( c(1:nlevels(data$treatment)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )



### linear modeling on alpha

alpha_div_model <- lm(Shannon ~ fire + herbiv, data = shann)
summary(alpha_div_model)
plot(alpha_div_model)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  4.37681    0.13945  31.387   <2e-16 ***
#fireN       -0.03652    0.17241  -0.212    0.832    
#herbivN      0.14346    0.17001   0.844    0.400  
#Residual standard error: 0.4915 on 189 degrees of freedom
#Multiple R-squared:  0.001791,	Adjusted R-squared:  -0.008772 
#F-statistic: 0.1695 on 2 and 189 DF,  p-value: 0.8442

Pval <- anova(alpha_div_model) [] #creates a nice table

qqnorm(rstudent(alpha_div_model)); abline(0,1) 

#We can test whether the residuals are normally distributed (then we don't need the data to be normally distributed): Kolmogorov-Smirnov test. If p value >0.05 => we are OK

ks.test(rstudent(alpha_div_model), pnorm, mean = mean(rstudent(alpha_div_model)), sd=sd((rstudent(alpha_div_model))))
#data:  rstudent(alpha_div_model)
#D = 0.084009, p-value = 0.133
#alternative hypothesis: two-sided    

alpha_div_model2 <- lm(Shannon ~ fire * herbiv, data = temp)
summary(alpha_div_model2)
plot(alpha_div_model2)

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    5.39133    0.06515  82.755   <2e-16 ***
#fireN          0.02101    0.10301   0.204    0.839    
#herbivN        0.04210    0.09213   0.457    0.648    
#fireN:herbivN -0.12259    0.14458  -0.848    0.398    
#Residual standard error: 0.4919 on 188 degrees of freedom
#Multiple R-squared:  0.005593,	Adjusted R-squared:  -0.01027 
#F-statistic: 0.3525 on 3 and 188 DF,  p-value: 0.7874

#REM : U1, U22 and U2 are outliers (but within the cook's distance...)

Pval <- anova(alpha_div_model2) [] #creates a nice table

qqnorm(rstudent(alpha_div_model2)); abline(0,1) 

#We can test whether the residuals are normally distributed (then we don't need the data to be normally distributed): Kolmogorov-Smirnov test. If p value >0.05 => we are OK

ks.test(rstudent(alpha_div_model2), pnorm, mean = mean(rstudent(alpha_div_model2)), sd=sd((rstudent(alpha_div_model2))))
#data:  rstudent(alpha_div_model2)
#D = 0.093381, p-value = 0.07027


shannon_U <- read.csv("shannon_U.csv")
head(shannon_U)
alpha_div_U <- lm(Shannon ~ fire.freq * herbiv, data = shannon_U)
summary(alpha_div_U)
plot(alpha_div_U)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         5.278607   0.227896  23.162   <2e-16 ***
#fire.freqC          0.076155   0.322293   0.236   0.8149    
#fire.freqT         -0.001613   0.338024  -0.005   0.9962    
#herbivN            -0.027993   0.322293  -0.087   0.9314    
#fire.freqC:herbivN -1.179046   0.455791  -2.587   0.0152 *  
#fire.freqT:herbivN -0.010332   0.478038  -0.022   0.9829    
#Residual standard error: 0.5582 on 28 degrees of freedom
#Multiple R-squared:  0.4231,	Adjusted R-squared:  0.3201 
#F-statistic: 4.107 on 5 and 28 DF,  p-value: 0.006367

alpha_div_U2 <- lm(Shannon ~ fire.freq + herbiv, data = shannon_U)
summary(alpha_div_U2)
plot(alpha_div_U2)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.488194   0.207009  26.512   <2e-16 ***
#fire.freqC  -0.513368   0.251689  -2.040   0.0503 .  
#fire.freqT  -0.006779   0.263974  -0.026   0.9797    
#herbivN     -0.447166   0.211461  -2.115   0.0429 *  
#Residual standard error: 0.6165 on 30 degrees of freedom
#Multiple R-squared:  0.2461,	Adjusted R-squared:  0.1707 
#F-statistic: 3.264 on 3 and 30 DF,  p-value: 0.03495

shannon_B <- read.csv("shannon_B.csv")
head(shannon_B)
alpha_div_B <- lm(Shannon ~ fire.freq * herbiv, data = shannon_B)
summary(alpha_div_B)
plot(alpha_div_B)
#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          5.5052     0.1897  29.017   <2e-16 ***
#fire.freqB          -0.1380     0.2683  -0.514    0.612    
#fire.freqC          -0.5206     0.3118  -1.670    0.110    
#herbivN             -0.1557     0.2814  -0.553    0.586    
#fire.freqB:herbivN   0.3768     0.3980   0.947    0.355    
#fire.freqC:herbivN       NA         NA      NA       NA    
#Residual standard error: 0.4647 on 21 degrees of freedom
#Multiple R-squared:  0.2504,	Adjusted R-squared:  0.1077 
#F-statistic: 1.754 on 4 and 21 DF,  p-value: 0.1759

alpha_div_B2 <- lm(Shannon ~ fire.freq + herbiv, data = shannon_B)
summary(alpha_div_B2)
plot(alpha_div_B2)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.41955    0.16638  32.573   <2e-16 ***
#fire.freqB   0.03324    0.19769   0.168   0.8680    
#fire.freqC  -0.62332    0.29156  -2.138   0.0439 *  
#herbivN      0.03268    0.19852   0.165   0.8708    
#Residual standard error: 0.4636 on 22 degrees of freedom
#Multiple R-squared:  0.2185,	Adjusted R-squared:  0.1119 
#F-statistic:  2.05 on 3 and 22 DF,  p-value: 0.1362































##########################################################################
######################## BETA DIVERSITY (between-sample)##################
##########################################################################

#### 4 steps: ;  

# (1) transform; for beta, pre-processing is necessary

# (2) calculate distance: many options of distance measurement = ##bray## (= most common), Unifrac, wunifrac, dpcoa, jsd, manhattan, euclidean, Canberra, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, maximum, binary, minkowski, etc 

# (3) ordinate (many options); ordination methods in plyloseq are:for unconstrained -  "CCA", "DPCoA", "NMDS", "MDS", "PCoA", "DCA", "RDA","CAP", the most commom = ##PCoA and NMDS##.for constrained : "CCA", and specify the model (phyloseq~treatment, "CCA")

# (4) plot (per sample, per OTU, biplot, etc): types of plot : type = c("samples", "sites", "species", "taxa", "biplot", "split", "scree")


# Callahan

log_phyf2 = transform_sample_counts(phyf2, function(x) log(1 + x) )

#katie

total = median(sample_sums(phy))#find median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
M.std = transform_sample_counts(phy, standf)#apply to phyloseq object
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)

#with tree
M.stdtree = transform_sample_counts(phytree, standf)#apply to phyloseq object
M.ftree = filter_taxa(M.std,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)


# two different and widely-used ecological distances: Bray-Curtis and weighted UniFrac, are sensitive to differences in total counts and a deluge of rare taxa. Let's transform to relative abundance, and also apply those filtering criteria that we explored earlier. http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html

mdt = fast_melt(phy)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
# Define taxa to keep.
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
# Define new object with relative abundance
mpra = transform_sample_counts(mp, function(x) x / sum(x))
# Filter this new object
mpraf = prune_taxa(keepTaxa, mpra)



#########################################################
################## UNCONSTRAINED ORDINATION #############



#from Dan Knights video (vegan) (look at the code...)

d.euc <- dist(OTU) #euclidian
pc.euc <- cmdscale(d.euc, k=2) #Principal coordinates

CS <- cca(otu_table(phy)) # chi square. didn't work with cca(OTU)
PC <- CS$CA$u[,1:2]

BC <- vegdist(otu_table(phy)) #bray curtis = default for vegan
PCB <- cmdscale(BC, k=2)
plot(PCB[,1], PCB[,2], col = ..., cex = 3, pch = 16)

  
### PCoA and NMDS - BRAY-CURTIS distances ###

#D <- phyloseq::distance(logt, method = "bray") # distance () didn't work
#O <- ordinate(logt, "PCoA", distance = "D")
#p <- plot_ordination(logt, O, color="...", shape="...")


PCoA_bray_logphyf2 <- ordinate(log_phyf2, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
#evals <- PCoA_bray_logphyf2$values$Eigenvalues

p <- plot_ordination(log_phyf2, PCoA_bray_logphyf2, color="site", shape="treatm", title= "PCoA - Bray curtis distance - fungal populations") + coord_fixed(sqrt(evals[2]/evals[1])) + theme_bw()
p


p <- plot_ordination(log_phyf2, PCoA_bray_logphyf2, color="site", shape="FH", title= "PCoA - Bray curtis distance - fungal populations") + theme_bw()
p + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
p + stat_ellipse(type = "t")


set.seed(2) #For NMDS plots it's important to set a seed since the starting positions of samples in the alogrithm is random.select any start seed, in this case 2 (ensures reproducibility) 
NMDS_bray_logphyf2 <- ordinate(log_phyf2, "NMDS", "bray") #other arguments that might be added: k=2, trymax=100) # stress=0.06
evals <- NMDS_bray_logphyf2$values$Eigenvalues
p <- plot_ordination(log_phyf2, NMDS_bray_logphyf2, color="site", shape="FH", title= "NMDS - Bray curtis distance - fungal populations") + theme_bw()
p



GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100) # stress=0.06

color = c("FH")
shape = c("site")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")

MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"), 
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+ 
  theme_bw() + labs(color=color, shape=shape) + geom_point(size=5)
MDS.1
pdf(paste0(outDir,"/NMDS_tretment_Bray_Curtis.pdf"),8,5)
MDS.1
dev.off()


out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
out.pcoa.M.f <- ordinate(M.f, method = "PCoA", distance = "bray")

evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.M.f, type = "samples", color = c("site"), shape = c("FH"), title = "PCoA of bacterial Communities - BC distance"
) + labs(col = "Treatment (F and H)") + coord_fixed(sqrt(evals[2] / evals[1])) + theme_bw()

# you can also add a facet_wrap(~Phylum, 3) 



### UNIFRAC distance (which takes OTU relatedness into account) ####

GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")
GP.ord.U

color = c("Treatm")
#shape = c("dog)
title=c("NMDS of 16S microbiome, Unifrac distance,k=2")
MDS = plot_ordination(M.f, GP.ord.U, color = color, title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
  theme_bw()+labs(color=color)+geom_point(size=5)
MDS.1

pdf(paste0(outDir,"/NMDS_treatment_Unifrac.pdf"),8,5)
MDS.1
dev.off()


### you can plot all the methods and compare them. http://joey711.github.io/phyloseq/plot_ordination-examples.html

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="site")
}, M.ftree, dist)

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=SampleType, shape=human, fill=SampleType))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
p

# If you want to replot a larger version of an individual plot, you can do by printing from the original plist from which pdataframe was made. Each element of plist is already a ggplot2 graphic. For example, printing the second element of the list.
plist[[2]] 



############ OTU plot and BIPLOT #############################


plot_ordination(M.f, out.pcoa.M.f, type = "taxa", color = c("Phylum"),shape = "site", title = "PCoA of bacterial Communities - BC distance") + labs(col = "Treatment (F and H)") +  theme_bw()

plot_ordination(M.f, out.pcoa.M.f, type = "split", color = c("Phylum"),shape = "site", title = "PCoA of bacterial Communities - BC distance") + labs(col = "Treatment (F and H)") +  theme_bw()

plot_ordination(M.f, out.pcoa.M.f, type = "split", color = c("site"), shape = c("FH"), title = "PCoA of bacterial Communities - BC distance"
) + labs(col = "Treatment (F and H)") + theme_bw()














#########################################################
################## CONSTRAINED ORDINATION #############
##################################################

# to see how environmental variables are associated with these changes in community composition. We constrain the ordination axes to linear combinations of environmental variables. We then plot the environmental scores onto the ordinatio
# most commom = CCA and RDA

####FIRST TEST : 
#remove data points with missing data

M.f.nona <- M.f %>% subset_samples(!is.na(pH_KCl) & ! is.na(S_C_.) & ! is.na(Av_Mg_cmolc_kg) & ! is.na(Av_Na_cmolc_kg)& ! is.na(Av_K_cmolc_kg)& ! is.na(Av_P_cmolc_kg) & ! is.na(Field_Moist_.DW) & ! is.na(Clay_.) & ! is.na(Silt_.) & ! is.na(Sand_.) & ! is.na(Bdfine_g_cm3)) #didn't find how to remove rows with NA's ...


M.f_bray <- phyloseq::distance(physeq = M.f.nona, method = "bray")

# CAP ordinate
cap_ord <- ordinate(
  physeq = M.f.nona, 
  method = "CAP",
  distance = M.f_bray,
  formula = ~ pH_KCl + S_C_. + Av_Mg_cmolc_kg + Av_Na_cmolc_kg + Av_K_cmolc_kg + Av_P_cmolc_kg + Field_Moist_.DW + Clay_. + Silt_. + Sand_. + Bdfine_g_cm3
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = M.f.nona, 
  ordination = cap_ord, 
  color = "site", 
  axes = c(1,2)
) + 
  aes(shape = FH) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) + theme_bw()

#Do a permutational ANOVA on constrained axes used in ordination

anova(cap_ord)





####### TEST 1 CONSTRAINED ORDINATION --- Mix FROM PDF PRETORIA and VEGAN TUTO

# FOR ALL SITES (too big...)

Enviveg <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_tot.csv", header = T)
row.names(Enviveg) <- Enviveg$Sample.ID #define the sample names as row names
str(Enviveg)
head(Enviveg)
Enviveg2 <-  dplyr::select(Enviveg, -Sample.ID, -Sample_name, -site, -site_zone, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -rep.1, -Zone, -site_FH_zone, -site_ttm_zone, -Block, -plot, -L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -Grass_SR, -Forb_SR, -Tree_SR, -Tree_10m_SR, -Total_SR, -Shannon, -Simpson, -spAbundance, -X.bareC, -X.litterC, -X.BasalC, -X.FoliarC, -ndvi, -woody_cover, -V_N_., -V_dN15_14, -V_C_., -V_dC13_12, -V_C.N, -Acid_cmolc_kg, -S_dC13_12)
head(Enviveg2)
str(Enviveg2)
rowstoremove <- c("B37", "B38", "B39", "B40", "U41", "U42", "N91", "N92", "N93", "N94", "N95", "L97", "L98", "L99", "L100", "L101", "L102", "L103", "B12", "L2", "U30", "N54", "U21", "N42") # = samples merged and rows with NA (in the soil mic biomass)
Envi = Enviveg2[!row.names(Enviveg2)%in%rowstoremove,]
rownames(Envi)
str(Envi)

#test correlations:
library(GGally)
ggpairs(EnvU2[,70:80])


# upload ASV table
load("phyf2.RData")
species <- otu_table(phyf2)
species[1:5,1:5]
species <- t(species)
s <- as.data.frame(species)
rowstoremove <- c("B37", "B38", "B39", "B40", "U41", "U42", "N91", "N92", "N93", "N94", "N95", "L97", "L98", "L99", "L100", "L101", "L102", "L103", "B12", "L2", "U30", "N54", "U21", "N42") # = samples merged and rows with NA (in the soil mic biomass)
s = s[!row.names(s)%in%rowstoremove,]

identical(rownames(Envi_norm), rownames(s)) #true
all(rownames(s) %in% rownames(Envi_norm)) #true
rownames(Envi_norm)
rownames(s)
gplots::venn(list(Envi=rownames(Envi_norm), species=rownames(s)))

# Standardise the environmental data : transforms all values for mean of 0 and variance of 1
Envi_norm <- decostand(Envi, "stand")

# normalise the species data
s_norm <- decostand(s, "hellinger")

# Checking for collinearity (excessive correlation) among explanatory variables
species.rda <- rda(s_norm~., data=Envi_norm)
species.rda1 <- rda(s_norm~1, data=Envi_norm)
vif.cca(species.rda)

#S_N_.         S_dN15_14             S_C_.             S_C.N            pH_KCl     Av_Ca_cmol_kg     Av_Mg_cmol_kg       Av_Na_mg_kg 
#1.450382e+03      2.759065e+00      1.404969e+03      1.983029e+01      1.349041e+01      1.011559e+01      1.986155e+01      5.969323e+00 
#Av_K_mg_kg        Av_P_mg_kg   Field_Moist_.DW   Field_Capac_.DW            Clay_.            Silt_.          Vfsand_.           Fsand_. 
#5.038530e+00      9.064386e+00      7.455866e+01      3.975229e+01      4.734129e+01      4.554751e+01      2.749576e+01      4.305480e+01 
#Msand_.           Csand_.          Vcsand_.            Sand_.       Stones_.tot       Bdtot_g_cm3      Bdfine_g_cm3    Bdhybrid_g_cm3 
#2.429286e+01      5.465644e+01                NA                NA      5.613526e+00      3.963972e+01      2.419017e+02      1.246373e+02 
#Residmoist_.airdw          L550_.dw           LOI_.dw     L550.1000_.dw           totMg_.           totAl_.           totSi_.            totP_.
#2.283558e+01      1.108506e+05      1.145156e+05      5.769848e+02      8.576630e+00      4.188762e+01      5.986330e+01      1.273042e+01 
#totS_.           totCl_.            totK_.           totCa_.           totTi_.            totV_.           totCr_.           totMn_. 
#9.529774e+01      4.880015e+00      6.943522e+01      4.874152e+01      2.075469e+01      2.914602e+01      1.870416e+02      1.863820e+01 
#totFe_.           totNi_.           totCu_.           totZn_.           totGa_.           totBr_.           totRb_.           totSr_. 
#2.311541e+02      1.584055e+02      4.561276e+01      4.317303e+01      4.481172e+01      1.226645e+02      5.418415e+01      7.053938e+01 
#totY_.           totZr_.           totNb_.           totBa_.           totHf_.           totTa_.           totPb_.           totTh_. 
#7.820669e+00      2.098627e+01      3.762313e+01      2.310020e+01      8.300308e+00      6.199586e+00      4.581016e+01      1.145096e+01 
#SMBC_ugC_gdw  veg_biomass_g_m2 
#3.559782e+00      6.232264e+00 

# Removal of non-significant explanatory variables... TAKES LOOOONG
step.res <- ordiR2step(species.rda1, scope=(species.rda), perm.max = 999, direction="forward")





####For UKULINGA:

EnvU3 <- dplyr::select(EnvU2, -Sample.ID, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Block, -plot, -Silt_., -V_N_., -V_C.N, -X.bareC, -Acid_cmolc_kg, -Field_Capac_.DW, -Shannon, -Simpson, -Vfsand_., -SMBC_ugC_gdw, -spAbundance, -V_C_., -V_N_., -V_dN15_14, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -S_dN15_14, -V_dC13_12, -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_.,-totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., -totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

#transform the envi table
EnvU_norm <- decostand(EnvU3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvU_norm), species=rownames(OTUU)))

# Checking for collinearity (excessive correlation) among explanatory variables
species.rdaU <- rda(OTUU~., data=EnvU3)
species.rda1U <- rda(OTUU~1, data=EnvU3)
vif.cca(species.rdaU)
# Removal of non-significant explanatory variables
step.resU <- ordistep(species.rda1U, scope=(species.rdaU), direction="forward", pstep=1000)
R2step.resU <- ordiR2step(species.rda1U, scope=(species.rdaU), direction="forward", pstep=1000) # both gave the same model.
#Step: OTUU ~ Av_Mg_cmol_kg + Av_Na_mg_kg + totK_. + totAl_. + Av_P_mg_kg + veg_biomass_g_m2 + Total_SR + X.litterC + X.BasalC + Bdtot_g_cm3 + totBr_. + totCa_. 
#second times (good time) Step: OTUU ~ Av_Mg_cmol_kg + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 +      pH_KCl + X.litterC + Bdtot_g_cm3 + Av_Ca_cmol_kg + X.BasalC +      Total_SR
rda.stepU <- rda(OTUU ~ Av_Mg_cmol_kg + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + pH_KCl + X.litterC + Bdtot_g_cm3 + Av_Ca_cmol_kg + X.BasalC + Total_SR, data = EnvU3)
anova(rda.stepU)
RsquareAdj(rda.stepU)$adj.r.squared	 #0.2798361 #second : [1] 0.2515413
summary(rda.stepU, display = NULL)
vif.cca(rda.stepU) 
screeplot(rda.stepU)
anova.cca(rda.stepU, step = 1000) 
anova.cca(rda.stepU, by = "term", step = 1000)
anova(rda.stepU, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

#rem: the whole procedure with normalized envi table gave the same results
species.rdaUnorm <- rda(OTUU~., data=EnvU_norm)
species.rda1Unorm <- rda(OTUU~1, data=EnvU_norm)
step.resUnorm <- ordistep(species.rda1Unorm, scope=(species.rdaUnorm), direction="forward", pstep=1000)

rda.stepUnorm <- rda(OTUU ~ Av_Mg_cmol_kg + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + pH_KCl + X.litterC + Bdtot_g_cm3 + Av_Ca_cmol_kg + Total_SR + X.BasalC, data = EnvU_norm)

#from package colortools
wheel("#009999")
Square()

#plot in phyloseq

p0 = plot_ordination(U, rda.stepUnorm, type = "biplot", color = "FH")
p0
# Now add the environmental variables as arrows
arrowmat = vegan::scores(rda.stepUnorm, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = RDA1, yend = RDA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * RDA1, y = 1.2 * RDA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.05, "npc"))
p1 = p0 + geom_segment(arrow_map, size = 0.5, data = arrowdf, color = "gray", 
                       arrow = arrowhead) + geom_text(label_map, size = 4, data = arrowdf) + theme_bw()
p1


# plot in Vegan

plot(rda.stepU, type = "n", main = "RDA - fungal population - Ukulinga", xlab = "RDA1 [10.5%]", ylab = "RDA2 [10%]")
points(rda.stepU, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvU2, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvU2, points(rda.stepU, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda.stepU, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
#text(rda.stepU, display = "species", scaling = scl, cex = 0.8, col = "darkcyan") = if I want to add the species text
#legend: each work individually, but not both
#with(EnvU2, legend("topright", legend = levels(FH), bty = "n", col = colvec, pch = 21, pt.bg = colvec))
#legend("topright", legend = "species", cex = 0.8, col = "grey", pch = 3)
legend("topright", legend = c(levels(EnvU2$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))

plot(rda.stepU, type = "n", main = "RDA - fungal population - Ukulinga", xlab = "RDA1 [10.5%]", ylab = "RDA2 [10%]")
points(rda.stepU, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvU2, levels(Treatm)) #"BH" "BN" "CH" "CN" "TH" "TN"
colvec <- c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900")
with(EnvU2, points(rda.stepU, display = "sites", col = colvec[Treatm], pch = 21, bg = colvec[Treatm]))
text(rda.stepU, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("topright", legend = c(levels(EnvU2$Treatm), "species"), col = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"), pch = c(21, 21, 21, 21, 21, 21, 3), pt.bg = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"))


varpartU <- varpart(OTUU, ~ Av_Mg_cmol_kg + Av_Na_mg_kg + Av_P_mg_kg + pH_KCl + Bdtot_g_cm3 + Av_Ca_cmol_kg, ~ veg_biomass_g_m2 + X.litterC + X.BasalC + Total_SR, data=EnvU3)
varpartU


# same but with cca - ukulinga

species.ccaU <- cca(OTUU~., data=EnvU_norm)
species.cca1U <- cca(OTUU~1, data=EnvU_norm)

# Removal of non-significant explanatory variables
step.ccaU <- ordistep(species.cca1U, scope=(species.ccaU), direction="forward", pstep=1000)
rda.ccaU <- cca(OTUU ~ pH_KCl + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + Av_Mg_cmol_kg + X.litterC + Bdtot_g_cm3 + Total_SR, data = EnvU_norm)
#OTUU ~ pH_KCl + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + Av_Mg_cmol_kg + Bdtot_g_cm3 + X.litterC + Total_SR + X.BasalC
#OTUU ~ pH_KCl + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + Av_Mg_cmol_kg + X.litterC + Bdtot_g_cm3 + Total_SR + X.BasalC + V_N_. + Av_K_mg_kg #this one, I realised the list was not perfect: totNb was still there... and I added tot soil N
#OTUU ~ pH_KCl + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + Av_Mg_cmol_kg + X.litterC + Bdtot_g_cm3 + Total_SR
anova(rda.ccaU)
RsquareAdj(rda.ccaU)$adj.r.squared	 #0.2798361 #second : [1] 0.2515413 ... 0.1722828
summary(rda.ccaU, display = NULL)
vif.cca(rda.ccaU) 
screeplot(rda.ccaU)
anova.cca(rda.ccaU, step = 1000) 
anova.cca(rda.ccaU, by = "term", step = 1000)
anova(rda.ccaU, by="axis", perm=500) #Moreover, it is possible to analyse significance of each axis:

# plot in Vegan

plot(rda.ccaU, type = "n", main = "CCA - fungal population - Ukulinga", xlab = "CCA1 [7.9 %]", ylab = "CCA2 [7.3%]")
points(rda.ccaU, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvU2, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvU2, points(rda.ccaU, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda.ccaU, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomleft", legend = c(levels(EnvU2$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))

plot(rda.ccaU, type = "n", main = "CCA - fungal population - Ukulinga", xlab = "CCA1 [7.9 %]", ylab = "CCA2 [7.3%]")
points(rda.ccaU, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvU2, levels(Treatm)) #"BH" "BN" "CH" "CN" "TH" "TN"
colvec <- c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900")
with(EnvU2, points(rda.ccaU, display = "sites", col = colvec[Treatm], pch = 21, bg = colvec[Treatm]))
text(rda.ccaU, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("topright", legend = c(levels(EnvU2$Treatm), "species"), col = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"), pch = c(21, 21, 21, 21, 21, 21, 3), pt.bg = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"))

plot(rda.ccaU, type = "n", main = "CCA - fungal population - Ukulinga", xlab = "CCA1 [7.9 %]", ylab = "CCA2 [7.3%]")
points(rda.ccaU, display = "species", cex = 0.8, scaling = 1, col = "#999999", pch = 3)
with(EnvU2, levels(Treatm)) #"BH" "BN" "CH" "CN" "TH" "TN"
colvec <- c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900")
with(EnvU2, points(rda.ccaU, display = "sites", col = colvec[Treatm], pch = 21, bg = colvec[Treatm]))
text(rda.ccaU, display = "bp", cex = 0.8, col = "black", scaling = 1) #text for the arrows
legend("topright", legend = c(levels(EnvU2$Treatm), "species"), col = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"), pch = c(21, 21, 21, 21, 21, 21, 3), pt.bg = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"))


dbRDA_U <- capscale(OTUU ~ pH_KCl + Av_Na_mg_kg + Av_P_mg_kg + veg_biomass_g_m2 + Av_Mg_cmol_kg + X.litterC + Bdtot_g_cm3 + Total_SR, data = EnvU_norm, distance = "bray")
summary(dbRDA_U, display = NULL)
plot(dbRDA_U, type = "n", main = "db-RDA - fungal population - Ukulinga", xlab = "CAP1 [16.0 %]", ylab = "CAP2 [13.7%]")
points(dbRDA_U, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvU2, levels(Treatm)) #"BH" "BN" "CH" "CN" "TH" "TN"
colvec <- c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900")
with(EnvU2, points(dbRDA_U, display = "sites", col = colvec[Treatm], pch = 21, bg = colvec[Treatm]))
text(dbRDA_U, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("topright", legend = c(levels(EnvU2$Treatm), "species"), col = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"), pch = c(21, 21, 21, 21, 21, 21, 3), pt.bg = c("#990000", "#FF9900", "#4D9900", "#99FF00", "#004D99", "#009999", "#999900", "#CCCCCC"))



varpartU <- varpart(OTUU, ~ Av_Mg_cmol_kg + Av_Na_mg_kg + Av_P_mg_kg + pH_KCl + Bdtot_g_cm3 + Av_Ca_cmol_kg, ~ veg_biomass_g_m2 + X.litterC + X.BasalC + Total_SR, data=EnvU3)
varpartU





#BROTHERTON

EnvB3 <- dplyr::select(EnvB2, -Sample.ID, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Block, -plot, -S_N_., -Silt_., -V_N_., -V_C.N, -Vfsand_., -SMBC_ugC_gdw, -V_C_., V_N_., -V_dN15_14, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -S_dN15_14, -V_dC13_12, -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_.,-totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

EnvB_norm <- decostand(EnvB3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvB3), species=rownames(OTUB)))

# Checking for collinearity (excessive correlation) among explanatory variables
species.rdaB <- rda(OTUB~., data=EnvB_norm)
species.rda1B <- rda(OTUB~1, data=EnvB_norm)

vif.cca(species.rdaB)

# Removal of non-significant explanatory variables
step.resB <- ordistep(species.rda1B, scope=(species.rdaB), direction="forward", pstep=1000)
#First, keeping all significant oneStep: OTUB ~ V_dN15_14 + totAl_. + totV_. + V_N_. + totTi_. 
# second attempt : Step: OTUB ~ Field_Moist_.DW + Av_K_mg_kg + Acid_cmolc_kg + totY_. +      S_C.N 
#last one (good) (OTUB ~ Field_Moist_.DW + Av_K_mg_kg + Acid_cmolc_kg + Av_Na_mg_kg +      S_C.N + Av_Mg_cmol_kg + Stones_.tot)

rda.stepB <- rda(OTUB ~ Field_Moist_.DW + Av_K_mg_kg + Acid_cmolc_kg + Av_Na_mg_kg + S_C.N + Av_Mg_cmol_kg + Stones_.tot, data = EnvB_norm) # = ordination

RsquareAdj(rda.stepB)$adj.r.squared	 #0.2798361 #second : [1] 0.2515413 ] > Second one : 0.1725754... tot Y explained a lot!!!
summary(rda.stepB, display = NULL)
vif.cca(rda.stepB) 
screeplot(rda.stepB)
anova.cca(rda.stepB, step = 1000) 
anova.cca(rda.stepB, by = "term", step = 1000)
anova(rda.stepB, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

# plot in Vegan

plot(rda.stepB, type = "n", main = "RDA - fungal population - Brotherton", xlab = "RDA1 [11.9%]", ylab = "RDA2 [7.0%]")
points(rda.stepB, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvB2, levels(FH)) #"AH" "AN" "BH" "BN" "CN"
colvec <- c("#009999", "#4D0099", "#4D9900")
with(EnvB2, points(rda.stepB, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda.stepB, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomright", legend = c(levels(EnvB2$FH), "species"), col = c("#009999", "#4D0099", "#4D9900", "#999999"), pch = c(21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#4D9900", "#999999"))

plot(rda.stepB, type = "n", main = "RDA - fungal population - Brotherton", xlab = "RDA1 [11.9%]", ylab = "RDA2 [7.0%]")
points(rda.stepB, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvB2, levels(Treatm)) #"AH" "AN" "BH" "BN" "CN"
colvec <- c("#990000", "#FF9900", "#004D99", "#009999", "#4D9900")
with(EnvB2, points(rda.stepB, display = "sites", col = colvec[Treatm], pch = 21, bg = colvec[Treatm]))
text(rda.stepB, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomright", legend = c(levels(EnvB2$Treatm), "species"), col = c("#990000", "#FF9900", "#004D99", "#009999", "#4D9900", "#999999"), pch = c(21, 21, 21, 21, 21, 3), pt.bg = c("#990000", "#FF9900", "#004D99", "#009999", "#4D9900", "#999999"))





#LETABA riparian

EnvLR3 <- dplyr::select(EnvL2R, -Sample.ID, -Zone, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Silt_., -L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -V_N_., -V_C.N, -Vfsand_., -SMBC_ugC_gdw, -V_C_., -V_N_., -V_dN15_14, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -S_dN15_14, -V_dC13_12, -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_., -totMg_., -totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., -totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

colnames(EnvLR3)

EnvLR_norm <- decostand(EnvLR3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvLR_norm), species=rownames(OTULR)))

# Checking for collinearity (excessive correlation) among explanatory variables
rdaLR <- rda(OTULR~., data=EnvLR_norm)
rda1LR <- rda(OTULR~1, data=EnvLR_norm)
step.rdaLR <- ordistep(rda1LR, scope=(rdaLR), direction="forward", pstep=1000)
rda_LR <- rda(OTULR ~ veg_biomass_g_m2 + S_N_. + Stones_.tot + pH_KCl , data = EnvLR_norm) # = ordination

RsquareAdj(rda_LR)$adj.r.squared	 #[1] 0.1513148
summary(rda_LR, display = NULL)
vif.cca(rda_LR) 
screeplot(rda_LR)
anova.cca(rda_LR, step = 1000) 
anova.cca(rda_LR, by = "term", step = 1000)
anova(rda_LR, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

# plot in Vegan

plot(rda_LR, type = "n", main = "RDA - fungal population - Letaba riparian", xlab = "RDA1 [14.5 %]", ylab = "RDA2 [4.9%]")
points(rda_LR, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvL2R, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvL2R, points(rda_LR, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda_LR, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomleft", legend = c(levels(EnvL2R$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))


varpartU <- varpart(OTULR, ~ veg_biomass_g_m2, ~ S_N_. + Stones_.tot + pH_KCl, data=EnvLR_norm)
varpartU

#LETABA CREST

EnvLC3 <- dplyr::select(EnvL2C, -Sample.ID, -Zone, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Silt_., -L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -V_N_., -V_C.N, -Vfsand_., -SMBC_ugC_gdw, -V_C_., -V_N_., -V_dN15_14, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -S_dN15_14, -V_dC13_12, -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_., -totMg_., -totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., -totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

colnames(EnvLC3)

EnvLC_norm <- decostand(EnvLC3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvLC_norm), species=rownames(OTULC)))

# Checking for collinearity (excessive correlation) among explanatory variables
rdaLC <- rda(OTULC~., data=EnvLC_norm)
rda1LC <- rda(OTULC~1, data=EnvLC_norm)
step.rdaLC <- ordistep(rda1LC, scope=(rdaLC), direction="forward", pstep=1000)
rda_LC <- rda(OTULC ~ Field_Moist_.DW + Av_P_mg_kg + S_N_. + veg_biomass_g_m2 + Av_Na_mg_kg , data = EnvLC_norm) # = ordination

RsquareAdj(rda_LC)$adj.r.squared	 #[1] 0.1513148
summary(rda_LC, display = NULL)
vif.cca(rda_LC) 
screeplot(rda_LC)
anova.cca(rda_LC, step = 1000) 
anova.cca(rda_LC, by = "term", step = 1000)
anova(rda_LC, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

# plot in Vegan

plot(rda_LC, type = "n", main = "RDA - fungal population - Letaba Crest", xlab = "RDA1 [5.3 %]", ylab = "RDA2 [4.8%]")
points(rda_LC, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvL2C, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvL2C, points(rda_LC, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda_LC, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomleft", legend = c(levels(EnvL2C$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))


varpartU <- varpart(OTULC, ~ Field_Moist_.DW + Av_P_mg_kg + S_N_. + Av_Na_mg_kg, ~ veg_biomass_g_m2 , data = EnvLC_norm)
varpartU

#NKUHLU riparian

EnvNR3 <- dplyr::select(EnvN2R, -Sample.ID, -Zone, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Silt_., -Vfsand_., -SMBC_ugC_gdw, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_., -totMg_., -totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., -totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

#-L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -V_N_., -V_C.N, -V_C_., -V_N_., -V_dN15_14, -S_dN15_14, -V_dC13_12, 

EnvNR_norm <- decostand(EnvNR3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvNR_norm), species=rownames(OTUNR)))

# Checking for collinearity (excessive correlation) among explanatory variables
rdaNR <- rda(OTUNR~., data=EnvNR_norm)
rda1NR <- rda(OTUNR~1, data=EnvNR_norm)
step.rdaNR <- ordistep(rda1NR, scope=(rdaNR), direction="forward", pstep=1000)
rda_NR <- rda(OTUNR ~ woody_cover + pH_KCl + veg_biomass_g_m2 + Av_K_mg_kg + Av_P_mg_kg, data = EnvNR_norm) # = ordination

RsquareAdj(rda_NR)$adj.r.squared	 #[1]  0.1137933
summary(rda_NR, display = NULL)
vif.cca(rda_NR) 
screeplot(rda_NR)
anova.cca(rda_NR, step = 1000) 
anova.cca(rda_NR, by = "term", step = 1000)
anova(rda_NR, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

# plot in Vegan

plot(rda_NR, type = "n", main = "RDA - fungal population - Nkuhlu riparian", xlab = "RDA1 [10.5 %]", ylab = "RDA2 [5.0 %]")
points(rda_NR, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvN2R, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvN2R, points(rda_NR, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda_NR, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomleft", legend = c(levels(EnvN2R$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))


varpartNR <- varpart(OTUNR, ~ woody_cover + veg_biomass_g_m2, ~ pH_KCl +  Av_K_mg_kg + Av_P_mg_kg, data = EnvNR_norm)
varpartNR

#NKUHLU CREST

EnvNC3 <- dplyr::select(EnvN2C, -Sample.ID, -Zone, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -V_N_., -V_C.N, -V_C_., -V_N_., -V_dN15_14, -S_dN15_14, -V_dC13_12, -Silt_., -Vfsand_., -SMBC_ugC_gdw, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_., -totMg_., -totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., -totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

#-L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -V_N_., -V_C.N, -V_C_., -V_N_., -V_dN15_14, -S_dN15_14, -V_dC13_12, 

EnvNC_norm <- decostand(EnvNC3, "stand")

#the OTU table: take the one from the adonis section: hellinger transformed
gplots::venn(list(Envi=rownames(EnvNC_norm), species=rownames(OTUNC)))

# Checking for collinearity (excessive correlation) among explanatory variables
rdaNC <- rda(OTUNC~., data=EnvNC_norm)
rda1NC <- rda(OTUNC~1, data=EnvNC_norm)
step.rdaNC <- ordistep(rda1NC, scope=(rdaNC), direction="forward", pstep=1000)
rda_NC <- rda(OTUNC ~ Av_Ca_cmol_kg + veg_biomass_g_m2 + pH_KCl + Tree_10m_SR + Av_P_mg_kg + Av_K_mg_kg, data = EnvNC_norm) # = ordination

RsquareAdj(rda_NC)$adj.r.squared	 #[1]  0.09852904
summary(rda_NC, display = NULL)
vif.cca(rda_NC) 
screeplot(rda_NC)
anova.cca(rda_NC, step = 1000) 
anova.cca(rda_NC, by = "term", step = 1000)
anova(rda_NC, by="axis", perm=500) #Moreover, it is possible to analyse signicance of each axis:

# plot in Vegan

plot(rda_NC, type = "n", main = "RDA - fungal population - Nkuhlu Crest", xlab = "RDA1 [6.5 %]", ylab = "RDA2 [4.6 %]")
points(rda_NC, display = "species", cex = 0.8, scaling = 2, col = "#999999", pch = 3)
with(EnvN2C, levels(FH)) #"FH" "FN" "NH" "NN"
colvec <- c("#009999", "#4D0099", "#CC3300", "#4D9900")
with(EnvN2C, points(rda_NC, display = "sites", col = colvec[FH], pch = 21, bg = colvec[FH]))
text(rda_NC, display = "bp", cex = 0.8, col = "black", scaling = 2) #text for the arrows
legend("bottomright", legend = c(levels(EnvN2C$FH), "species"), col = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"), pch = c(21, 21, 21, 21, 3), pt.bg = c("#009999", "#4D0099", "#CC3300", "#4D9900", "#CCCCCC"))


varpartNR <- varpart(OTUNR, ~ woody_cover + veg_biomass_g_m2, ~ pH_KCl +  Av_K_mg_kg + Av_P_mg_kg, data = EnvNR_norm)
varpartNR

####################################################################
####################### STATISTICS #################################
####################################################################

# before doing stats tests: remove rare taxa (present in less than 10% of the samples), because no stat power on them.

####################################################################
################## Differential abundance testing ##################
####################################################################


####################################################################
######### Univariate differential abundance testing##################

#When significant global differences in microbiome composition are detected between groups of samples, a natural question arises: which particular taxa are responsible of that global difference? A common strategy to answer this question is to test every taxa separately for association with the response variable. When the response variable is dichotomous this is known as univariate differential abundance testing.

# First part of the script = from the Dan Knights video. But I should use metagenomeSeq. ttest, to compare alpha diversity between 2 groups, or ANOVA if comparing more than 2 groups. the equivalent nonparametric tests, like the Wilcoxon rank-sum test or the Kruskal-Wallis test, can be applied. However, more powerful parametric approaches are available, such as the Bioconductor packages edgeR [16] and DESeq2 [34], initially proposed for transcriptomics analysis (RNA-Seq data). Both fit a generalized linear model and assume that read counts follow a Negative Binomial distribution. 

#is arthrobacter present ?

taxM <- as.data.frame(tax)
filter(taxM, Genus == "Arthrobacter") # 59 ASV

Arthrobacter <- phy_genus_table[,grep("Arthrobacter", colnames(phy_genus_table))]
hist(Arthrobacter, br=30)
cor.test(Arthrobacter, Env$fire) #doesn't work... must be numeric vector => not with fire
FitArt_fire <- lm(Arthrobacter ~ Env$fire)
FitArt_fire
anovaart_fire <- anova(FitArt_fire)
anovaart_fire
summary(FitArt_fire)
Pval <- anova(FitArt_fire) []

FitArt_FH <- lm(Arthrobacter ~ Env$fire + Env$herbiv)
FitArt_FH
summary(FitArt_FH)
anovaart_FH <- anova(FitArt_FH)
anovaart_FH

Pval <- anova(FitArt_fire) []
Pval
class(Pval)
rownames(Pval)
colnames(Pval)

qqnorm(rstudent(FitArt_fire)); abline(0,1) #very bad !!!

ks.test(rstudent(FitArt_fire), pnorm, mean = mean(rstudent(FitArt_fire)), sd=sd((rstudent(FitArt_fire)))) #very bad!#We can test whether the residuals are normally distributed (then we don't need the data to be normally distributed): Kolmogorov-Smirnov test. If p value >0.05 => we are OK
#One-sample Kolmogorov-Smirnov test
#data:  rstudent(FitArt_fire)
#D = 0.35024, p-value < 2.2e-16

kruskal.test(Arthrobacter ~ Env$fire) #Kruskal-Wallis chi-squared = 0.47474, df = 1, p-value = 0.4908

Fit <- lm(Arthrobacter ~ Env$fire + Env$herbiv + Env$site)
summary(Fit)

#####################################################################
########## make a loop to test if any of the phylums is significantly correlated to Fire and Herbivory

common.ids <- intersect(rownames(Env), rownames(phy_phylum_table)) 
# get just the overlapping samples 
phy_phylum_table2 <- phy_phylum_table[common.ids,] 
Envi <- Env[common.ids,]
dim(phy_phylum_table2)
colnames(phy_phylum_table2)[1:10]

pvals <- numeric(ncol(phy_phylum_table2)) # pvals is a vector initialized with zeroes with enough slots for the different genera
names(pvals) <- colnames(phy_phylum_table2) # "name" the pvalues after the genera
for(i in 1:ncol(phy_phylum_table2)) {
  fit <- lm(phy_phylum_table2[,i] ~ Envi$fire + Envi$herbiv)
  pvals[i] <-  anova(fit)['Envi$fire','Pr(>F)']
}

sort(pvals)[1:10] # print the 10 smallest p-values

#List of all pvalues and see which genus are highly associated with a variable: "apply" with genus, 2 means do something to every column of genus # ("apply" with genus, 1 would mean do something every row) # the last part defines a new function to do to each column, which # will be passed in the temporary variable named "xx" 
pvals <- apply(phy_phylum_table2, 2, function(xx) anova(lm(xx ~ Envi$fire * Envi$herbiv))['Envi$fire','Pr(>F)']) 

# print the 10 smallest p-values: 
sort(pvals)[1:10]

#To correct for multiple test, we can adjust the p-values using the p.adjust function. we are correcting for multiple hypothesis testing use False Discovery Rate (FDR). The adjusted p-values are often called "q-values."

qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10] # print the lowest 10 q-values


########################################################################################
#### test which genus are significantly correlated with Fire, when we correct for herbivory

common.idsG <- intersect(rownames(Env), rownames(phy_genus_table)) 
# get just the overlapping samples 
phy_genus_table2 <- phy_genus_table[common.idsG,] 
Envi <- Env[common.idsG,]
dim(phy_genus_table2)
dim(Env)
colnames(phy_genus_table2)[1:10]

pvals <- numeric(ncol(phy_genus_table2)) # pvals is a vector initialized with zeroes with enough slots for the different genera
names(pvals) <- colnames(phy_genus_table2) # "name" the pvalues after the genera
for(i in 1:ncol(phy_genus_table2)) {
  fit <- lm(phy_genus_table2[,i] ~ Envi$fire + Envi$herbiv)
  pvals[i] <-  anova(fit)[Envi$fire,'Pr(>F)']
}



pvals
qvals <- p.adjust(pvals,'fdr') # correction for multiple hypothesis testing
sort(qvals)[1:10] # print the lowest 10 q-values => nothing is significant now ! 0.2 = a cutoff. all are above

#### test if lm is appropriate (residuals normally distributed)

ks.pvals <- numeric(ncol(phy_genus_table2))
names(ks.pvals) <- colnames(phy_genus_table2)# "name" the pvalues after the genera
options(warn=-1)# turn annoying warnings off 
# Loop through the columns of the genus table, testing each one
for(i in 1:ncol(phy_genus_table2)) {
  fit <- lm(phy_genus_table2[,i] ~ Envi$fire + Envi$herbiv)
  ks.pvals[i] <-  ks.test(rstudent(fit), pnorm, mean=mean(rstudent(fit)), sd=sd(rstudent(fit)),exact=FALSE)$p.value
}
options(warn=0)# turn warnings back on

# Now since we ran 97 tests we should correct for multiple hypothesis testing. What we want: none of these test significant (not significantly not normal). So that we can stick with our linear regression.
ks.qvals <- p.adjust(ks.pvals,'fdr') #here, a lot of them are non normal

# print the lowest 10 q-values
sort(ks.qvals)[1:10]

plot(-log10(seq(0,1,length=98)[-1]), -log10(sort(ks.pvals))); abline(0,1)


######## => we need to use a generalized linear model (negative binomial in edgeR). the negative binomial uses raw counts of sequences (rarefied or not), not the relative abundances. => we need row counts! NO PRETREATMENTS !!!! Good and convenient (from dan knight)


source("C:/Users/VRMMAR024/Dropbox/R code/Wrap.edgeR.R") 

result <- glm.edgeR(x=Envi$fire, Y=phy_genus_table2) # X = independant variable; Y = a matrix, the whole table of genera. and the test is going to be made on each column.

topTags(result) # gives the best hits for genus vs fire. when you look at FRD: even the topone not significant.
#                              logFC   logCPM       LR      PValue       FDR
#Segetibacter             -2.9867956 9.858930 8.212448 0.004160399 0.4882331


result <- glm.edgeR(x=Envi$fire, Y=phy_genus_table2, covariates = Envi[,c("herbiv", "site")])
topTags(result)
#                   logFC    logCPM        LR       PValue         FDR
#Reyranella     -0.4459072 13.088793 18.291769 1.895243e-05 0.009419356
#Segetibacter   -3.0162431  9.858882 15.621211 7.738168e-05 0.019229346
#Devosia         0.6947721 10.861090  9.283722 2.311992e-03 0.383020032
#Luteolibacter  -1.6599354  8.802902  8.281643 4.004787e-03 0.426502438

pvals <- topTags(result,n=Inf)$table[,'PValue']
write.table(topTags(result, n=Inf)$table, file='edgeR_results.txt',sep='\t',quote=FALSE, col.names=NA) # to write the whole table
sum(topTags(result,n=Inf)$table$FDR <= 0.05) # to find how many significantly different: 2 !!!




################# now test if genus are different in different sites

#nkuhlu <- Envi$site == "N"
#resultN <- glm.edgeR(x=nkuhlu, Y=phy_genus_table2, covariates = Envi[,c("herbiv", "site")])
#topTags(result) # highly different !!!!!
#                      logFC    logCPM       LR       PValue          FDR
#Gp2               -4.734747 13.886875 38.97174 4.299857e-10 2.137029e-07
#TT <- topTags(result, n=Inf) # means = i don't care how many genera there are, print all results!
#pvals <- topTags(result,n=Inf)$table[,'PValue']
#plot(-log10(seq(0,1,length=98)[-1]), -log10(sort(pvals))); abline(0,1)



###############################################
##################### non parametric tests

# Mann-Whitney U test, sometimes called the Wilcoxon Signed Rank test or Wilcoxon Rank Sum test.
# =	For testing differences between two categories. Make sure to use the relative abundances, not the absolute abundances.

wilcox.test(prevotella ~ is.USA, exact=FALSE) # get the exact p-value wilcox.test(prevotella ~ is.USA, exact=FALSE)$p.value

# Kruskal-Wallis test. test for differentiation across multiple categories, analogous to ANOVA
kruskal.test(prevotella ~ map$COUNTRY)

# Spearman correlation instead of Pearson correlation.-	Continuous test (For continuous variables)

cor.test(prevotella, map$AGE, method='spearman', exact=FALSE)







####################################################################
######### Multivariate differential abundance testing ##############
####################################################################


#refers to a global test of differences in microbial composition between two or more groups of samples. We can distinguish between distance-based or model-based approaches. PERMANOVA, ANOSIM (analysis of group similarities),  MRPP (Multi-response permutation procedures), and MANTEL (Mantel's test). they perform on distance matrixes




########### ANOSIM ####################


#The anosim function performs a non-parametric test of the significance of the sample-grouping you provide against a permutation-based null distribution, generated by randomly permuting the sample labels many times (999 permutations is the default, used here).

OTU_tab <- otu_table(phyf2)
head(OTU_tab)
OTU_tab[1:10,1:10]
tOTU_tab <- t(OTU_tab) # to have samples as rows
tOTU_tab[1:10,1:10]
dim(tOTU_tab) #192 86998

envi <- sample_data(phyf2)
head(envi)
dim(envi) #192  93
# didn't work : site <- envi[,3] # to extract the parameter to test (dim : #192   1)
site = get_variable(phyf2, "site")
dim(site) #null

#Construct a distance hemi-matrix for biodata

otu_dist<-vegdist(tOTU_tab)
dim(otu_dist)

#Test null hypothesis of no difference in otu_table_transposed composition between the sites

otu_anosim<-anosim(otu_dist, site, permutations = 999, distance = "bray")
otu_anosim
#Call:
#anosim(x = otu_dist, grouping = site, permutations = 999, distance = "bray") 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.8162 
#            The R-statistic in ANOSIM is a ratio between within-group and between-group dissimilarities. 
#            0.75 < R < 1 - highly different
#            0.5 < R < 0.75 - different
#            0.25 < R < 0.5 - different with some overlap
#            0.1 < R < 0.25 - similar with some differences (or high overlap)
#            R < 0.1 - similar

#Significance: 0.001 
#Permutation: free
#Number of permutations: 999


#other way: 

site_anosim = anosim(distance(phyf2, "bray"), site)
site_anosim$signif #0.001
site_anosim$statistic # (R value test for statistical significance between the groups defined) here : 0.816229





#######################################################################
################### PERMANOVA #### ADONIS #############################

#Permutational Multivariate Analysis of Variance Using Distance Matrices, 
#is perhaps the most widely used distance-based method for multivariate community analysis. Implemented in the function "adonis" of the vegan R package,Unfortunately, there are currently no post-hoc tests developed for adonis.


# FIRST: transform the dataset: filter rare taxa (katie's method), because no statistical power on it. Then do a hellinger transformation : square root og the relative abundance. 

load("phyf2.RData")

NITS <- subset_samples(phyf2, site =="N")
NITS
sum(taxa_sums(NITS) ==0) # 38550 => tot = 48448 taxa
#otu_table()   OTU Table:         [ 86998 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 86998 taxa by 7 taxonomic ranks ]

LITS <- subset_samples(phyf2, site =="L")
LITS
sum(taxa_sums(LITS) ==0) #63685 => tot = 23313
#otu_table()   OTU Table:         [ 86998 taxa and 64 samples ]
#sample_data() Sample Data:       [ 64 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 86998 taxa by 7 taxonomic ranks ]

BITS<- subset_samples(phyf2, site =="B")
BITS
sum(taxa_sums(BITS) ==0) #75918 => 11080
#otu_table()   OTU Table:         [ 86998 taxa and 26 samples ]
#sample_data() Sample Data:       [ 26 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 86998 taxa by 7 taxonomic ranks ]

UITS <- subset_samples(phyf2, site =="U")
UITS
sum(taxa_sums(UITS) ==0) #72319 => 14679
#otu_table()   OTU Table:         [ 86998 taxa and 34 samples ]
#sample_data() Sample Data:       [ 34 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 86998 taxa by 7 taxonomic ranks ]





total = median(sample_sums(phyf2))#find median sample read count
phyf_f = filter_taxa(phyf2,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE) #from katie
ntaxa(phyf_f) # we have retained 15195 from the 86998 taxa
save(phyf_f, file = paste0(outDir,"/phyf_f.RData")) 
Helling_phyf_f <- transform_sample_counts(phyf_f, function(x) sqrt(x / sum(x)))


rtreephyf_fb = filter_taxa(rtreephyf2,function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE) #from katie: M.std (here phy.filtered) = standardized to the median
ntaxa(rtreephyf_fb) # (katie's name : M.f)  we have retained 15195 from the 86998 taxa
save(rtreephyf_fb, file = paste0(outDir,"/rtreephyf_fb.RData")) 

phyf2
Helling_phyf_f
#otu_table()   OTU Table:         [ 15195 taxa and 192 samples ]
#sample_data() Sample Data:       [ 192 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]


#extract the individual sites

N <- subset_samples(Helling_phyf_f, site =="N") # select only samples from Nkuhlu
N
#otu_table()   OTU Table:         [ 15195 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]
B <- subset_samples(Helling_phyf_f, site =="B")
B
#otu_table()   OTU Table:         [ 15195 taxa and 26 samples ]
#sample_data() Sample Data:       [ 26 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]
L <- subset_samples(Helling_phyf_f, site =="L")
L
#otu_table()   OTU Table:         [ 15195 taxa and 64 samples ]
#sample_data() Sample Data:       [ 64 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]

U <- subset_samples(Helling_phyf_f, site =="U")
U
#otu_table()   OTU Table:         [ 15195 taxa and 34 samples ]
#sample_data() Sample Data:       [ 34 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]

LC <- subset_samples(L, Zone =="C")
LC
#otu_table()   OTU Table:         [ 15195 taxa and 34 samples ]
#sample_data() Sample Data:       [ 34 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]
LR <- subset_samples(L, Zone =="R")
LR
#otu_table()   OTU Table:         [ 15195 taxa and 30 samples ]
#sample_data() Sample Data:       [ 30 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]
NC <- subset_samples(N, Zone =="C")
NC
#otu_table()   OTU Table:         [ 15195 taxa and 38 samples ]
#sample_data() Sample Data:       [ 38 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]
NR <- subset_samples(N, Zone =="R")
NR
#otu_table()   OTU Table:         [ 15195 taxa and 30 samples ]
#sample_data() Sample Data:       [ 30 samples by 100 sample variables ]
#tax_table()   Taxonomy Table:    [ 15195 taxa by 7 taxonomic ranks ]


# create vegan-compatible OTU
SAMD <-as(sample_data(U), "data.frame")
SAMD <- dplyr::select(SAMD, -site, -site_zone, -Zone, -totMg_.)
SAMDat <- sample_data(SAMD)

OTUU <- otu_table(U)
TAXU <- tax_table(U)
phyloseqU <- phyloseq(OTUU,TAXU, SAMDat)
phyloseqU

OTUtot = otu_table(Helling_phyf_f)
taxa_are_rows(OTUtot) #true
OTUtot = t(OTUtot)

OTUN <- otu_table(N)
taxa_are_rows(OTUN) #true => transpose 
OTUN = t(OTUN)

OTUL <- otu_table(L)
taxa_are_rows(OTUL) #true => transpose 
OTUL = t(OTUL)

OTUB <- otu_table(B)
taxa_are_rows(OTUB) #true => transpose 
OTUB = t(OTUB)

OTUU <- otu_table(U)
taxa_are_rows(OTUU) #true => transpose 
OTUU = t(OTUU)

OTULC <- otu_table(LC)
taxa_are_rows(OTULC) #true => transpose 
OTULC = t(OTULC)
OTULR <- otu_table(LR)
taxa_are_rows(OTULR) #true => transpose 
OTULR = t(OTULR)

OTUNC <- otu_table(NC)
taxa_are_rows(OTUNC) #true => transpose 
OTUNC = t(OTUNC)
OTUNR <- otu_table(NR)
taxa_are_rows(OTUNR) #true => transpose 
OTUNR = t(OTUNR)


#create vegan-compatible envi table. or check that all samples in the newly created OTU table in the envi table

rowstoremove <- c("B37", "B38", "B39", "B40", "U41", "U42", "N91", "N92", "N93", "N94", "N95", "L97", "L98", "L99", "L100", "L101", "L102", "L103") # = samples merged and rows with NA
Env2 = Env[!row.names(Env)%in%rowstoremove,]
rownames(Env2)


rowstoremoveN <- c("N91", "N92", "N93", "N94", "N95") # = samples merged and rows with NA
EnvN2 = EnvN[!row.names(EnvN)%in%rowstoremoveN,]
gplots::venn(list(rownames(EnvN2), rownames(OTUN)))
EnvN2C <- filter(EnvN2, Zone =="C") 
row.names(EnvN2C) <- EnvN2C$Sample.ID#define the sample names as row names
EnvN2R <- filter(EnvN2, Zone =="R")
row.names(EnvN2R) <- EnvN2R$Sample.ID#define the sample names as row names

rowstoremoveL <- c("L97", "L98", "L99", "L100", "L101", "L102", "L103") # = samples merged and rows with NA
EnvL2 = EnvL[!row.names(EnvL)%in%rowstoremoveL,]
rownames(EnvL2)
gplots::venn(list(rownames(EnvL2), rownames(OTUL)))
EnvL2C <- filter(EnvL2, Zone =="C") %>% dplyr::select(-"site_zone")
row.names(EnvL2C) <- EnvL2C$Sample.ID#define the sample names as row names
EnvL2R <- filter(EnvL2, Zone =="R") %>% dplyr::select(-"site_zone")
row.names(EnvL2R) <- EnvL2R$Sample.ID#define the sample names as row names

rowstoremoveU <- c("U41", "U42") # = samples merged and rows with NA
EnvU2 = EnvU[!row.names(EnvU)%in%rowstoremoveU,]
gplots::venn(list(rownames(EnvU2), rownames(OTUU)))

rowstoremoveB <- c("B37", "B38", "B39", "B40") # = samples merged
EnvB2 <-  EnvB[!row.names(EnvB)%in%rowstoremoveB,]
gplots::venn(list(rownames(EnvB2), rownames(OTUB)))




# Calculate bray curtis distance matrix
set.seed(1)

Helling_phyf_f_bray <- phyloseq::distance(Helling_phyf_f, method = "bray")
otu_table(Helling_phyf_f)[1:5, 1:5] # ASV = rows -> transpose to have samples as rows and species as columns
N_Bray <- phyloseq::distance(N, method = "bray")
L_Bray <- phyloseq::distance(L, method = "bray")
U_Bray <- phyloseq::distance(U, method = "bray")
B_Bray <- phyloseq::distance(B, method = "bray")
NC_Bray <- phyloseq::distance(NC, method = "bray")
NR_Bray <- phyloseq::distance(NR, method = "bray")
LR_Bray <- phyloseq::distance(LR, method = "bray")
LC_Bray <- phyloseq::distance(LC, method = "bray")

tot_Bray <- vegdist(OTUtot, method = "bray")
N_Bray <- vegdist(OTUN, method = "bray")
L_Bray <- vegdist(OTUL, method = "bray")
U_Bray <- vegdist(OTUU, method = "bray")
B_Bray <- vegdist(OTUB, method = "bray")
NC_Bray <- vegdist(OTUNC, method = "bray")
NR_Bray <- vegdist(OTUNR, method = "bray")
LR_Bray <- vegdist(OTULR, method = "bray")
LC_Bray <- vegdist(OTULC, method = "bray")

PCoA_bray_U <- ordinate(U, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
Uplot <- plot_ordination(U, PCoA_bray_U, color="Treatm", title= "PCoA - Bray curtis distance - Ukulinga fungal populations") + theme_bw()
Uplot + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
Uplot2 <- plot_ordination(U, PCoA_bray_U, color="FH", title= "PCoA - Bray curtis distance - Ukulinga fungal populations") + theme_bw()
Uplot2 + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
pdf(paste0(outDir,"/PCoA-Bray-U-ITS.pdf"),8,5)
Uplot2
dev.off()

PCoA_bray_B <- ordinate(B, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
Bplot <- plot_ordination(B, PCoA_bray_B, color="Treatm", title= "PCoA - Bray curtis distance - Brotherton fungal populations") + theme_bw()
Bplot + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
Bplot2 <- plot_ordination(B, PCoA_bray_B, color="FH", title= "PCoA - Bray curtis distance - Brotherton fungal populations") + theme_bw()
Bplot2 + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses
pdf(paste0(outDir,"/PCoA-Bray-U-ITS.pdf"),8,5)
Uplot2
dev.off()

PCoA_bray_N <- ordinate(N, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
Bplot <- plot_ordination(N, PCoA_bray_N, color="FH", shape="Zone", title= "PCoA - Bray curtis distance - Nkuhlu fungal populations") + theme_bw()
Bplot + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses

PCoA_bray_L <- ordinate(L, "PCoA", "bray") #other arguments that might be added : k=2, trymax=100) # stress=0.06
Bplot <- plot_ordination(L, PCoA_bray_L, color="FH", shape="Zone", title= "PCoA - Bray curtis distance - Letaba fungal populations") + theme_bw()
Bplot + stat_ellipse(type = "norm", linetype = 2) # generate the confidence ellipses




# run a permanova of each individual column of the envi table

a_tot_fire=adonis(t(otu_table(Helling_phyf_f)) ~  fire, data = data.frame(sample_data(Helling_phyf_f)), method="bray", by="margin")
# exactly the same : adonis(Helling_phyf_f_bray ~ fire, data = data.frame(sample_data(Helling_phyf_f)))
a_tot_fire
a_tot_fire=adonis(OTUtot ~  fire, data = Env2, method="bray", by="margin")
a_tot <- adonis(tot_Bray ~ fire, data = Env2)
a_tot_fire
a_tot

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#fire        1     1.255 1.25505  3.1007 0.01606  0.001 ***
#Residuals 190    76.905 0.40476         0.98394           
#Total     191    78.160                 1.00000  


# if pvalue very small: => we can reject the null hypothesis that our 4 treatments have the same centroid

a_tot_herbiv=adonis(t(otu_table(Helling_phyf_f)) ~  herbiv, data = data.frame(sample_data(Helling_phyf_f)), method="bray", by="margin")
#by = "terms" will assess significance for each term (sequentially from first to last)
#by = "margin" will assess the marginal effects of the terms (each marginal term analysed in a model with all other variables)
#by = NULL will assess the overall significance of all terms together
a_tot_herbiv=adonis(t(otu_table(Helling_phyf_f)) ~  herbiv, data = data.frame(sample_data(Helling_phyf_f)), method="bray", by="terms")
a_tot_herbiv

a_tot_site=adonis(t(otu_table(Helling_phyf_f)) ~  site, data = data.frame(sample_data(Helling_phyf_f)), method="bray", by="margin")
a_tot_site


aS <- adonis(tot_Bray ~ site, Env2)
aS
aFH <- adonis(tot_Bray ~ herbiv * fire, Env2)
aFH
aFH <- adonis(tot_Bray ~ fire * herbiv, Env2)
aFH
aSZFH <- adonis(tot_Bray ~ site_zone * fire * herbiv, Env2)
aSZFH




#make a loop on all columns of the overall dataset

Sample <- data.frame(sample_data(Helling_phyf_f))
dim(Sample) #192 100
S2 <- Env2[,colSums(is.na(Env2)) == 0] # remove all columns with NA's
dim(S2) #192  74
#fungi <- t(otu_table(Helling_phyf_f))
dim(fungi) #192 15195

results <- lapply(colnames(S2), function(x){
  form <- as.formula(paste("Helling_phyf_f_bray", x, sep="~")) 
  z <- adonis(form, data = S2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
results

results <- lapply(colnames(S2), function(x){
  form <- as.formula(paste("tot_Bray", x, sep="~")) 
  z <- adonis(form, data = S2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
results

#make a loop on all columns Letaba

L2 <- EnvL2[,colSums(is.na(EnvL2)) == 0] # remove all columns with NA's
resultsL <- lapply(colnames(L2), function(x){
  form <- as.formula(paste("L_Bray", x, sep="~")) 
  z <- adonis(form, data = L2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsL

aZL <- adonis(L_Bray ~ Zone, EnvL2)
aZL
aFHL <- adonis(L_Bray ~ herbiv * fire, EnvL2)
aFHL
aFHL <- adonis(L_Bray ~ fire * herbiv, EnvL2)
aFHL
aSZFHL <- adonis(L_Bray ~ Zone * fire * herbiv, EnvL2)
aSZFHL


LC2 <- EnvL2C[,colSums(is.na(EnvL2C)) == 0] # remove all columns with NA's
LC2 <- dplyr::select(LC2, -Zone)
resultsLC <- lapply(colnames(LC2), function(x){
  form <- as.formula(paste("LC_Bray", x, sep="~")) 
  z <- adonis(form, data = LC2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsLC

aFHLC <- adonis(LC_Bray ~ herbiv * fire, EnvL2C)
aFHLC
aFHLC <- adonis(LC_Bray ~ fire * herbiv, EnvL2C)
aFHLC


LR2 <- EnvL2R[,colSums(is.na(EnvL2R)) == 0] # remove all columns with NA's
LR2 <- dplyr::select(LR2, -Zone)
resultsLR <- lapply(colnames(LR2), function(x){
  form <- as.formula(paste("LR_Bray", x, sep="~")) 
  z <- adonis(form, data = LR2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsLR

aFHLR <- adonis(LR_Bray ~ herbiv * fire, EnvL2R)
aFHLR
aFHLR <- adonis(LR_Bray ~ fire * herbiv, EnvL2R)
aFHLR




#make a loop on all columns Nkuhlu

N2 <- EnvN2[,colSums(is.na(EnvN2)) == 0] # remove all columns with NA's
resultsN <- lapply(colnames(N2), function(x){
  form <- as.formula(paste("N_Bray", x, sep="~")) 
  z <- adonis(form, data = N2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsN

aZN <- adonis(N_Bray ~ Zone, EnvN2)
aZN
aFHN <- adonis(N_Bray ~ herbiv * fire, EnvN2)
aFHN
aFHN <- adonis(N_Bray ~ fire * herbiv, EnvN2)
aFHN
aSZFHN <- adonis(N_Bray ~ Zone * fire * herbiv, EnvN2)
aSZFHN


NC2 <- EnvN2C[,colSums(is.na(EnvN2C)) == 0] # remove all columns with NA's
NC2 <- dplyr::select(NC2, -Zone)
resultsNC <- lapply(colnames(NC2), function(x){
  form <- as.formula(paste("NC_Bray", x, sep="~")) 
  z <- adonis(form, data = NC2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsNC

aFHNC <- adonis(NC_Bray ~ herbiv * fire, EnvN2C)
aFHNC
aFHNC <- adonis(NC_Bray ~ fire * herbiv, EnvN2C)
aFHNC

NR2 <- EnvN2R[,colSums(is.na(EnvN2R)) == 0] # remove all columns with NA's
NR2 <- dplyr::select(NR2, -Zone)
resultsNR <- lapply(colnames(NR2), function(x){
  form <- as.formula(paste("NR_Bray", x, sep="~")) 
  z <- adonis(form, data = NR2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsNR

aFHNR <- adonis(NR_Bray ~ herbiv * fire, EnvN2R)
aFHNR
aFHNR <- adonis(NR_Bray ~ fire * herbiv, EnvN2R)
aFHNR


#make a loop on all columns Brotherton

B2 <- EnvB2[,colSums(is.na(EnvB2)) == 0] # remove all columns with NA's
resultsB <- lapply(colnames(B2), function(x){
  form <- as.formula(paste("B_Bray", x, sep="~")) 
  z <- adonis(form, data = B2)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame
}
)
resultsB

aFHB <- adonis(B_Bray ~ herbiv * fire, EnvB2)
aFHB
aFHB <- adonis(B_Bray ~ fire * herbiv, EnvB2)
aFHB



#make a loop on all columns Ukulinga

U2 <- EnvU2[,colSums(is.na(EnvU2)) == 0] # remove all columns with NA's
resultsU <- lapply(colnames(U2), function(x){
  form <- as.formula(paste("U_Bray", x, sep="~")) 
  z <- adonis(form, data = U2)
  return(as.data.frame(z$aov.tab))} #convert anova table to a data frame
# return(data.frame(env = rownames(z$aov.tab), Rsq = z$aov.tab$R2,P = z$aov.tab$P))}
)
resultsU

aFHU <- adonis(U_Bray ~ herbiv * fire, EnvU2)
aFHU
aFHU <- adonis(U_Bray ~ fire * herbiv, EnvU2)
aFHU


#oTHER POSSIBILITY :  
#results<-list()
#for (i in colnames(EnvU4)){
#  form <- as.formula(paste("U_Bray", i, sep="~"))
#  results[[i]]<- adonis(form, data=EnvU4)
#} 


#A correlation matrix is a table of correlation coefficients for a set of variables used to determine if a relationship exists between the variables. The coefficient indicates both the strength of the relationship as well as the direction (positive vs. negative correlations).                
U2_num <-  dplyr::select(U2, -Sample.ID, -Sample_name, -fire, -Fire.freq, -herbiv, -Treatm, -FH, -rep, -site_FH_zone, -site_ttm_zone, -Block, -plot)

U2_cor <- cor(U2_num) #This returns a simple correlation matrix showing the correlations between pairs of variables (devices).You can choose the correlation coefficient to be computed using the method parameter. The default method is Pearson, but you can also compute Spearman or Kendall coefficients : mydata.cor = cor(mydata, method = c("spearman"))

#Significance levels (p-values) can also be generated using the rcorr function which is found in the Hmisc package. 
library("Hmisc")
U2_rcorr <- rcorr(as.matrix(U2_num)) #This generates one table of correlation coefficients (the correlation matrix) and another table of the p-values. By default, the correlations and p-values are stored in an object of class type rcorr. To extract the values from this object into a useable data structure, you can use the following syntax:
U2.coeff = U2_rcorr$r
U2.p = U2_rcorr$P

#to visualise: 
library(corrplot)
corrplot(U2_cor)# run the corrplot function providing our original correlation matrix as the data input to the function.
corrplot(U2_cor, method = "circle", tl.col="black", tl.cex=0.5)
corrplot.mixed(U2_cor)
corrplot(U2_cor, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, tl.cex=0.7)
corrplot(U2_cor, order = "hclust", tl.col = "black", tl.cex=0.7)
#palette = colorRampPalette(c("green", "white", "red")) 
#heatmap(x = U2_cor, col = palette, symm = TRUE)
heatmap(x = U2_cor, symm = TRUE)

#other way:
library(GGally)
ggpairs(EnvN2C[,12:20])


#farrar-glaubert test (didn't work...)
library(mctest)
omcdiag(EnvU, SMBC_ugC_gdw)



############################################################
#from Katie : USE PACKAGE VEGAN TO TEST SIGNIFICANCE BETWEEN LOCATIONS BASED ON DISTANCE MATRIX (ADONIS FUNCTION) - 	PERMUTATION MULTIVARIATE ANALYSIS OF VARIANCE USING DISTANCE MATRICES


#1. Using Bray-Curtis distance
diss1 <- distance(M.temp,method = "bray")
a1=adonis(t(otu_table(M.temp)) ~  BV + Nugent.score + Inflammation + location + injectable + BMI + Ethnicity + Age + HPV.risk, data = data.frame(sample_data(M.temp)), method="bray", by="margin")
a1
# a1
# 
# Call:
# adonis(formula = t(otu_table(M.temp)) ~ BV + Nugent.score + Inflammation +      location + injectable + BMI + Ethnicity + Age + HPV.risk,      data = data.frame(sample_data(M.temp)), method = "bray",      by = "margin") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BV             2     8.956  4.4782 19.0905 0.22361  0.001 ***
# Nugent.score   1     0.696  0.6955  2.9650 0.01736  0.016 *  
# Inflammation   1     0.293  0.2930  1.2491 0.00732  0.256    
# location       1     0.531  0.5306  2.2619 0.01325  0.057 .  
# injectable     1     0.218  0.2181  0.9297 0.00544  0.443    
# BMI            1     0.524  0.5239  2.2335 0.01308  0.053 .  
# Ethnicity      8     1.718  0.2148  0.9155 0.04289  0.658    
# Age            1     0.163  0.1633  0.6963 0.00408  0.627    
# HPV.risk       2     1.151  0.5756  2.4536 0.02874  0.014 *  
# Residuals    110    25.804  0.2346         0.64423           
# Total        128    40.054                 1.00000           
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > cat("Synch1485174494534372000\n");
#------

#Notes:
#1. The formula section in the adonis function can be either a data frame with your numeric values where samples = rows and proteins = columns or a dissimilarity matrix.
#2. After the ~ we specify the variables of interest in your case it could be fire, site, pH etc. 
#3. The data term is your dataframe with metadata where samples = rows and variables = columns (your 1st and 3rd steps have to match in terms of samples and you can't have any missing metadata otherwise you'll get an error). If you do have missing metadata you'll have to do something like this:
#remove sample that does not have contraceptives data (otherwise error)
temp.IDs <- rownames(m)[which(!is.na(m$injectable) & !is.na(m$Inflammation) & !is.na(m$location) & !is.na(m$Nugent.score) 
                              & !is.na(m$BMI) & !is.na(m$Ethnicity) & !is.na(m$Age) & !is.na(m$HPV.risk))]
M.temp <- prune_samples(temp.IDs,M.std)

#******
#NB NOTE: Function adonis evaluates terms sequentially. 
#In a model with right-hand-side ~ A + B the effects of A are evaluated first, and the effects of B after removing the effects of A, i.e. order matters!!
#SO ORDER THESE VARS BY MOST TO LEAST LIKELY (FROM WHAT WE ALREADY KNOW) TO INFLUENCE THE MICROBIOME - E.G. FOR MY MICROBIOME DATA IF WE PUT LOCATION BEFORE
#BMI - LOCATION WILL COME UP AS SIGNIFICANT (DUE TO SIGNIFICANT DIFFERENCE IN BMI BY LOCATION), BUT ACTUALLY ITS THE DIFFERENCE IN BMI BETWEEN THE LOCATIONS THAT APPEARS TO BE DRIVING THE DIFFERENCES IN MICROBIOME
#There is another function, adonis2() for which the order does not matter. However you don't get R2 values from adonis2() so run both..
#******
#NB NOTE 2:
#You should check that your data meets the assumption of homogeneous dispersion before you do adonis. For this you can use the betadisper() function:
grouping <- as.factor(unlist(sample_data(M.temp)[,"Inflammation"]))
betadine <- betadisper(diss5, group=grouping)#Where 'diss5' is your distance matrix
anova(betadine)#If NS then adonis meets assumptions of homogeneous dispersion - not sure what group to specify if you're testing multiple effects with adonis though..

#To visualize the contribution of different factors of interest with microbiota composition, you can do a constrained ordination. I used dbRDA on the Bray-Curtis dissimilarity matrix with ggplot, but you already have a constrained ordination plot so maybe don't need more details here.

#########################################################################################################







###############################################
########### MULTITABLE TECHNIQUE###############
#from callahan


Envi_metadata_ITS_num <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_num.csv", row.names=1)
Envnum <-  dplyr::select(Envi_metadata_ITS_num, -L_N_., -L_dN15_14, -L_C_., -L_dC13_12, -L_C.N, -Grass_SR, -Forb_SR, -Tree_SR, -Tree_10m_SR, -Total_SR, -X.BasalC, -X.FoliarC, -ndvi, -woody_cover, -V_N_., -V_dN15_14, -V_C_., -V_dC13_12, -V_C.N, -Acid_cmolc_kg, -S_dC13_12)
head(Envnum)
str(Envnum)
Envnum$Av_Na_mg_kg <- as.numeric(Envnum$Av_Na_mg_kg)
Envnum$Av_P_mg_kg <- as.numeric(Envnum$Av_P_mg_kg)
Envnum$Av_K_mg_kg <- as.numeric(Envnum$Av_K_mg_kg)
rowstoremove <- c("B37", "B38", "B39", "B40", "U41", "U42", "N91", "N92", "N93", "N94", "N95", "L97", "L98", "L99", "L100", "L101", "L102", "L103", "B12", "L2", "U30", "N54", "U21", "N42") # = samples merged and rows with NA
Envnum2 = Envnum[!row.names(Envnum)%in%rowstoremove,]
rownames(Envnum2)

library("genefilter")

Envchar <- read.csv("C:/Users/VRMMAR024/Dropbox/POST-DOC UCT/RESULTS/DNA/ITS/Envi_metadata_ITS_char.csv", header = T)
row.names(Envchar) <- Envchar$Sample.ID#define the sample names as row names
physimple <- phyloseq(OTU,TAX, sample_data(Envchar))
physimple
physimple2 <- filter_taxa(physimple, function (x) {sum(x) > 0}, prune=TRUE)
physimple2
physimple <- prune_taxa(taxa_sums(physimple) > 4, physimple)
physimple <- filter_taxa(physimple, filterfun(kOverA(3, 2)), TRUE) #KOverA(k, A=...) function from genefilter : A filter function for k elements larger than A.
physimple <- subset_samples(physimple, Sample.ID != "B37" & Sample.ID != "B38" & Sample.ID != "B39" & Sample.ID != "B40" & Sample.ID != "U41" & Sample.ID != "U42" & Sample.ID != "N91" & Sample.ID != "N92" & Sample.ID != "N93" & Sample.ID != "N94" & Sample.ID != "N95" & Sample.ID != "L97" & Sample.ID != "L98" & Sample.ID != "L99" & Sample.ID != "L100" & Sample.ID != "L101" & Sample.ID != "L102" & Sample.ID != "L103" & Sample.ID != "B12" & Sample.ID != "L2" & Sample.ID != "U30" & Sample.ID != "N54" & Sample.ID != "U21" & Sample.ID != "N42")

physimple


Envlog <- log(1 + Envnum2, base = 10)
X <- otu_table(physimple)
head(X)
X[X > 50] <- 50
dim(X) #19739   186
rownames(Envlog)
colnames(Envlog)
dim(Envlog) #186  64
 
apply(Envlog, 2, function(x) any(is.na(x))) #OK, aucun NA dans le Dataframe

library(PMA)
cca_res <- CCA(t(X),  Envlog, penaltyx = .15, penaltyz = .15) #sparse CCA. This method compares sets of features across high-dimensional data tables, where there may be more measured features than samples. In the process, it chooses a subset of available features that capture the most covariance - these are the features that reflect signals present across multiple tables. 

head(t(X))















#######################katie ##########################
####################################################################################

### Option 1 : katie. For diferential abundance testing we use raw reads (not standardized) as input because metagenomeSeq has its own internal scaling factor based on total read counts/sample. So lets merge the raw reads at lowest available taxonomy - again we're interested in biologically interpretable units, not strain-level variation

Mraw = tax_glom.kv(phy)
ntaxa(Mraw)#1572
total = median(sample_sums(phy))
Mraw.f = filter_taxa(Mraw,function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(Mraw.f)#1307
ntaxa(phy) #39871

#Let's subset the metagenomeSeq object to eliminate samples from dog K (we would like to compare dogs G & B)
#sub.index <- sample_names(Mraw.f)[sample_data(Mraw.f)[,"Dog"] != "K"]
#phy.temp <- prune_samples(sub.index, Mraw.f)
#nsamples(phy.temp)

#LEt's check if there are ASVs where our subset of samples now have only zero counts (i.e. ASVs that are exclusively found in dog K)
length(which(rowSums(otu_table(Mraw))==0))#0
length(which(rowSums(otu_table(Mraw.cleaned))==1))#0
keep <- which(rowSums(otu_table(Mraw))!=1)
#LEts remove these empty rows
Mraw.cleaned <- prune_taxa(names(keep),Mraw)
ntaxa(Mraw.cleaned)#106

seqps <- data.frame(sample_sums(Mraw.cleaned))

newMraw = subset_samples(Mraw.cleaned, sample_names(Mraw.cleaned) != "L57_S52_L001")

#The function we'll use to do differential abundance testing was built around metagenomeSeq's fitZig() function to build the model and MRfulltable() to extract the results. This customized function also filters the results by fold-change, percent presence across samples and adjusted p-values
#Results are summarized as a heatmap and .csv table
#There are several parameters that we can specify in the super.fitZig.kv() function (including)see microbiom_tutorial_katie for detail):
#Note that we do not HAVE to specify each and every one of these parameters. Some have default values that will be used if we don't specify anything.

b = super.fitZig.kv(physeq = newMraw, factor = "fire", outDir = outDir, FileName =c("1_25FC_0.2_fire"), heatmap.descriptor=c("tax_annot"), main=c("Fire vs no fire, merged taxa"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("site"))

b

############# RANDOM FOREST ################

RF <- RF.k(newMraw, seed = 2,var = "fire", cv.fold=10, outDir = outDir, train.control=FALSE,cv.repeat=3,Nfeatures.validation=NULL, testAndtrain=FALSE,np=50, positive.class = "F", descriptor = c("herbiv"),train.fraction=2/3)

#"THE TOP 50  MOST IMPORTANT FEATURES WERE:"
#predictors MeanDecreaseGini                                     tax
#ASV10         ASV10        1.4121764               Gaiella_occulta(JF423906)
#ASV14         ASV14        1.3620154                             Rhizobiales
#ASV31         ASV31        1.3600860          unidentified_bacterium(X64380)
#ASV70         ASV70        1.3442029          unidentified_bacterium(X64381)
#ASV8           ASV8        1.3121506                            Sphingomonas
#ASV236       ASV236        1.2127377                            Actinomadura
#ASV55         ASV55        1.1682123          unidentified_bacterium(X64382)
#ASV2           ASV2        1.1557335                                Bacillus
#ASV1085     ASV1085        1.0239169                        Chitinophagaceae
#ASV1           ASV1        0.9498887        Bradyrhizobium_retamae(KC247085)

plot(RF)

RF <- RF.k(newMraw, seed = 2,var = "herbiv", cv.fold=10, outDir = outDir, train.control=FALSE,cv.repeat=3,Nfeatures.validation=NULL, testAndtrain=FALSE,np=50, positive.class = "H", descriptor = c("fire"),train.fraction=2/3)



### option 2: differential abundance testing with DESeq2. 

#To test the differences at OTU level between treatments we need to convert the treatment column into factor. Note that we use the data without rarefaction
sample_data(phy)$fire <- as.factor(sample_data(phy)$fire)

#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(phy, ~ fire)
ds = DESeq(ds) #didn't work: Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc: every gene contains at least one zero, cannot compute log geometric means

#Extract the result table from the ds object usind the DESeq2 function results and filter the OTUs using a False Discovery Rate (FDR) cutoff of 0.01. In this example we return the significantly differentially abundant OTU between Fire and no fire:

alpha = 0.01
res = results(ds, contrast=c("Season", "Spring", "Fall"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig

















##################################################################################
############ VEGAN - from pdf pretoria ##################

##Let's check if both the OTU table (the transposed OTU table) and the environmental table are in the same order
rownames(Seqtab)==rownames(test_envitable)

##Alpha diversity (Richness, Shannon, Simpson, inv-Simpson, Pielou evenness)
richness <-rowSums(Seqtab>0) ## richness or number of species per sample

##Mean richness per habitat types
tapply(richness,test_envitable$Herbiv,mean)

##Student T-test: Is there a statistically significant difference in richness between the two habitat types?
t.test(richness ~ test_envitable$Herbiv, paired=T)
H <- diversity(otu_transposed) ## Shannon diversity

##Mean Shannon index per habitat types
tapply(H,env$Type,mean)



#############################drafts

phyf2glom <- tax_glom.kv(phyf2)
filterphyf2glom <- filter_taxa(phyf2glom, function(x) sum(x > 10) > (0.1*length(x)) | sum(x) > 0.001*total, TRUE)

MGS <- make_metagenomeSeq(phyf2glom)
SFK <- super.fitZig.kv(phyf2glom, factor="fire", covariate=NULL, p=0.05, FC = 1.5, perc=0.3, cexCol=1, cexRow=1)

TDF = (data.frame(OTUname = taxa_names(phy), row.names = taxa_names(phy)))
MGS <- make_metagenomeSeq(phy)
if (!taxa_are_rows(phy)) {
  physeq <- t(phy)
}
OTU = as(otu_table(phy), "matrix")
# Convert sample_data to AnnotatedDataFrame
ADF = AnnotatedDataFrame(data.frame(sample_data(phy)))
# define dummy 'feature' data for OTUs, using their name Helps with
# extraction and relating to taxonomy later on.
TDF = AnnotatedDataFrame(data.frame(tax_table(phy)))
# Create the metagenomeSeq object
MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)

# Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
MGS = cumNorm(MGS)
return(MGS)

p = cumNormStatFast(MGS)
#To calculate the scaling factors we simply run cumNorm
lungData = cumNorm(lungData, p = p)

