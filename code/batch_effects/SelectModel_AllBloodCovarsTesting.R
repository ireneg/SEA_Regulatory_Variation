# script written by KSB, 09.04.19

# Load dependencies and set input paths --------------------------

# Load dependencies:
library(edgeR)
library(plyr)
library(NineteenEightyR)
library(RColorBrewer)
library(biomaRt)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(pheatmap)
library(viridis)
library(gplots)
library(circlize)
library(ComplexHeatmap)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Load colour schemes:
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)

# Load log CPM matrix and y object:
# lcpm
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
# y DGE list object
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))

# BEGIN ANALYSIS -------------------------------------------------------------------------------------------------

# First, set up design matrix
# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45

# set up design matrix
design_simple <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$Gran)
design_noCD8T <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD4T + y$samples$NK + y$samples$Bcell + y$samples$Mono + y$samples$Gran)
design_noBcell <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$NK + y$samples$Mono + y$samples$Gran)
design_noNK <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$Bcell + y$samples$Mono + y$samples$Gran)
design_full <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$NK + y$samples$Bcell + y$samples$Mono + y$samples$Gran)

#  we can check what the best model is using limma's selectModel() function. This model uses Akaike's Information Criterion (AIC) to select the best model. i.e., the model with the lowest information criterion score.

designlist=list()
designlist$simple=design_simple
designlist$noCD8T=design_noCD8T
designlist$noBcell=design_noBcell
designlist$noNK=design_noNK
designlist$full=design_ful

out=selectModel(lcpm, designlist)
table(out$pref)

#  simple  noCD8T noBcell    noNK    full 
#   2488    2110    4098    3261    1018 

