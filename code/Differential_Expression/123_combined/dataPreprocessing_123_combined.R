# script created by KSB, 26.07.18
# perform data pre-processing on both batches of indonesian RNA-seq data- third batch

### Last edit: IGR 2019.10.19
### Removed MPI-296 having discovered she is female. 

### 0. Load dependencies and set input paths -------------------------- 
### 1. BEGIN ANALYSIS ----------------------------------------------------------- 
### 2. Covariate setup ------------------------------------------------------ 
### 3. Transformation from the raw scale and filtering of genes --------------------------------------------- 
### 4. Replicate analysis and plotting ---------------------------------------------------------------------- 
### 5. Other basic QC plots --------------------------------------- 
### 6. Normalise gene-expression distribution ----------------------------------------- 
### 7. Testing batch correction with PCA - moved here from the DE script instead



### TO DO:
### Fix some plotting functions
### Fix PCA plotting function for numeric covariates
### Fix malaria information
### Rarefaction curves as png file instead of pdf.

###########################################################################
### 0. Load dependencies and set input paths -------------------------- ###
###########################################################################

# Load dependencies:
library(edgeR) 
library(plyr)
library(NineteenEightyR)
library(openxlsx)
library(RColorBrewer)
library(magrittr)

# Set paths:
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296/" # on server
covariatedir <- "/data/cephfs/punim0586/igallego/indoRNA/"


# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296/"
edaoutput <- paste0(outputdir, "eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)


###################################################################################
# 1. BEGIN ANALYSIS ----------------------------------------------------------- ###
###################################################################################

# load in count data (as DGE list object)
load(paste0(inputdir, "all_samples_unfiltered_DGElistObject.Rda"))

# Create shortened samplenames- insert hyphen and drop suffixes
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# Create the covariate matrix --------------------------------------

covariates <- read.xlsx(paste0(covariatedir, "metadata_RNA_Batch123.xlsx"), sheet=1, detectDates <- T)

# add in blood
blood <- read.table(paste0(covariatedir, "predictedCellCounts_DeconCell.txt"), sep="\t", as.is=T, header=T)
colnames(blood) <- c("Gran","Bcell","CD4T","CD8T","NK","Mono")
blood$ID <- sapply(strsplit(rownames(blood),"[_.]"), `[`, 1)
blood$ID <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", blood$ID)
# add blood data into covariates matrix
covariates <- cbind(covariates, blood[match(covariates$Sample.ID, blood$ID),])

# add in replicate information
covariates$replicate <- duplicated(covariates[,1])

# The date column is being finnicky and importing strangely (i.e., some of the dates are recognised as dates, some are not). To fix this, I have to reformat the dates so that they can be handled correctly for the ANOVA analysis (below).
# first let's check which rows in the sampling date column are not in the correct format
unique(covariates[grep("-", covariates$Sampling.Date),6])
# "2016-03-10"

covariates$Sampling.Date <- gsub("2016-03-10","10/03/2016",covariates$Sampling.Date) %>% as.Date(., tryFormats="%d/%m/%Y")

# Check if any samples in the covariate DF are not in samplenames
covariates[which(!(covariates[,"Sample.ID"] %in% samplenames)),]
# MPI 296 comes up, as expected, so we remove her before we go any further:
covariates <- covariates[covariates$ID != "MPI-296",]
covariates[which(!(covariates[,"Sample.ID"] %in% samplenames)),]
# all good now!

#################################################################################
### 2. Covariate setup ------------------------------------------------------ ###
#################################################################################

# reorder covariates by sequencing batch
covariates <- covariates[order(covariates$Sequencing.Batch),]
# make sure all names in files are within the covariate matrix
all(samplenames == covariates[,"Sample.ID"])
# TRUE

  ###############################################################################
  ### IGR NOTE 2019.04.11 - THIS IS NOT WORKING RIGHT NOW - PATHS TO BE FIXED ###
  ###############################################################################

                                            # add in sick sample information. We can first read in the sick metadata which was generated in the script "indoRNA_DEAnalysis_STAR_HealthyvsSickMappi.R".
                                            # malaria.metadata <- read.table(paste0(inputdir, "DE_Sick/Malaria_summary_table.txt"), sep="\t", header=T)
                                            # # add in malaria information
                                            # covariates$microscopy.pos <- malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"microscopy"]
                                            # covariates$PCR.pos <- malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"PCR"]
                                            # covariates$fract.pfpx.reads <- malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"fract.reads.pfpx"]

# Assign covariates to DGE list
# subtract variables we don't need
subtract <- c("Sample.ID", "Kabisu.Ethnicity","Sequencing.Batch")
# get index of unwanted variables
subtract <- which(colnames(covariates) %in% subtract)
# subtract unwanted variables and add in library size and batch variables
covariate.names  <-  colnames(covariates)[-subtract]

# assign covariates to DGE list. Library size and batch are already assigned so subtract these from the loop.
for (name in covariate.names){
  y$samples[[paste0(name)]]<- covariates[,name]
}

# add in library size and batch infromation to covariate names vector
covariate.names  <-  c(covariate.names,"lib.size", "batch")

# Make sure all covariates in y$samples are the appropriate structure. This is important, as it affects the linear model and ANOVA tests
# Get variables that aren't numeric and transform them into factors. Batch is considered numeric (which we want to change), so we'll also add this into our make.factor vector
which.numeric <- unlist(lapply(y$samples, is.numeric))
make.factor <- c(colnames(y$samples[,!which.numeric]), "batch")
        
# assign variables and make y-samples metadata into appropriate structure
for (name in make.factor){
  if(sum(is.na(y$samples[[name]])) > 0){
    y$samples[[paste0(name)]] <- assign(name, as.factor(addNA(y$samples[[paste0(name)]])))
  }
  else{
    y$samples[[paste0(name)]] <- assign(name, as.factor(y$samples[[paste0(name)]]))
  }
}

# now do this for numeric variables
which.numeric <- unlist(lapply(y$samples, is.numeric))
make.numeric <- colnames(y$samples[,which.numeric])[c(1,3:ncol(y$samples[,which.numeric]))]

# assign numeric information to variables
for (name in make.numeric){
  assign(name, as.numeric(as.character(y$samples[[paste0(name)]])))
}

# Age, RIN, and library size need to be broken up into chunks for easier visualisation of trends (for instance in Age, we want to see high age vs low age rather than the effect of every single age variable)
for (name in c("Age","RIN","lib.size")){
  assign(name, cut(as.numeric(as.character(y$samples[[paste0(name)]])), breaks=5))
}

### IGR NOTE 2019.04.12 - DID IT ONCE NO NEED TO MAKE THEM EVERY TIME. ###

# Perform rarefaction curves for number of expressed genes vs. proportion of pool mRNA (as in Ramskold D, Wang ET, Burge CB, Sandberg R. 2009. An abundance of ubiquitously expressed genes revealed by tissue transcriptome sequence data. PLoS Comput Biol 5:e1000598)
# pdf(paste0(edaoutput, "rarefactionCurves_indoRNA_123Combined.pdf"))
#   for (name in colnames(Filter(is.factor,y$samples))[-c(1,2)]) {
#     plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=name) ## initialize the plot area
#     counter <- 0
#     for (sample in colnames(y)){
#       counter <- counter+1
#       lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=as.numeric(get(name))[counter])
#     }
#     levels <- levels(get(name))
#     levels[which(is.na(levels))]  <-  "NA"
#     legend(x="bottomright", bty="n", col=1:length(levels(get(name))), legend=levels, lty=1, lwd=2)
#   }
#   # add in batch (we do this separately since it has a different colour)
#   plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main="batch") ## initialize the plot area
#   counter <- 0
#   for (sample in colnames(y)){
#     counter <- counter+1
#     lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=batch.col[as.numeric(batch)][counter])
#   }
#   levels <- levels(batch)
#   levels[which(is.na(levels))]  <-  "NA"
#   legend(x = "bottomright", bty="n", col=batch.col[1:length(levels(batch))], legend=levels, lty=1, lwd=2)
# dev.off()


#################################################################################################################
### 3. Transformation from the raw scale and filtering of genes --------------------------------------------- ###
#################################################################################################################

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
indocpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE, prior.count=0.25) 

# Removal of lowly-expressed genes -------------------------------------------------------------------------------
# get histogram of number of genes expressed at log2 cpm > 0.5 and 1 (before filtering)
pdf(paste0(edaoutput, "lcpm_preFiltering_Histogram.pdf_cpm1.pdf"), height=10, width=15)
  par(mfrow=c(1,2))
  hist(rowSums(lcpm>0.5), main= "Genes expressed at log2 cpm over 0.5 \n pre-filtering", xlab="samples", col=4, ylim=c(0,16000))
  hist(rowSums(lcpm>1), main= "Genes expressed at log2 cpm over 1 \n pre-filtering", xlab="samples", col=5, ylim=c(0,16000))
dev.off()

# a gene is only retained if it is expressed at a CPM > 1 in at least half of the group (i.e., in half of MPI, half of SMB, or half of MTW)

# Extract CPM information for each village
cpm_SMB <- indocpm[,grep("SMB",colnames(indocpm))]
cpm_MPI <- indocpm[,grep("SMB",colnames(indocpm))]
cpm_MTW <- indocpm[,grep("MTW",colnames(indocpm))]

# set keep threshold to only retain genes that are expressed at or over a CPPM of one in at least half of the group
keep.SMB <- rowSums(cpm_SMB>1) >= (length(grep("SMB",rownames(y$samples)))*0.5)
keep.MPI <- rowSums(cpm_MPI>1) >= (length(grep("MPI",rownames(y$samples)))*0.5)
keep.MTW <- rowSums(cpm_MTW>1) >= (length(grep("MTW",rownames(y$samples)))*0.5)

table(keep.SMB)
  # keep.SMB
  # FALSE  TRUE 
  # 15076 12337 
table(keep.MPI)
  # keep.MPI
  # FALSE  TRUE 
  # 14447 12966 
table(keep.MTW)
  # keep.MTW
  # FALSE  TRUE 
  # 14997 12416 

# Keep genes that have a CPM >= 1 in at least half of the samples within one of the comparison groups
keep <- keep.SMB | keep.MPI | keep.MTW
yFilt <- y[keep,, keep.lib.sizes=FALSE]
dim(yFilt)
# [1] 12975   123

# eliminate composition bias between libraries by upper quartile normalisation (this will be used in the 'data exploration' step)
yFilt <- calcNormFactors(yFilt, method="TMM")
lcpmFilt <- cpm(yFilt, log=TRUE, prior.count=0.25) 

#################################################################################################################
### 4. Replicate analysis and plotting ---------------------------------------------------------------------- ###
#################################################################################################################

# set replicate names:
allreps <- covariates[,"Sample.ID"][which(covariates$replicate)]
allreps <- unique(allreps)

# Check before filtering:
pdf(paste0(edaoutput, "replicate_comparisons_noFiltering_123Combined.pdf"), height=10, width=10)
  par(mfrow=c(3,3))
  # Since SMB-ANK-027 had two replicates, we need to plot one vs three and two vs three. First, one vs three
  smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[1]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[1]], xlim=c(-5,15), ylim=c(-5,15), main="Technical replicate SMB-ANK-027")
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  # Now, replicates two vs three
  smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[2]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[2]], xlim=c(-5,15), ylim=c(-5,15), main="Technical replicate SMB-ANK-027")
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  for (i in 1:length(allreps)){
    smoothScatter(lcpm[,which(samplenames %in% allreps[i])[1]], lcpm[,which(samplenames %in% allreps[i])[2]], ylab=colnames(lcpm)[which(samplenames %in% allreps[i])[2]], xlab=colnames(lcpm)[which(samplenames %in% allreps[i])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate", allreps[i], sep=" "))
    # best fit/regression line
    abline(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]), col="green")
    # diagonal line
    abline(a=0,b=1,col="red")
  }
dev.off()

# And after:
pdf(paste0(edaoutput, "replicate_comparisons_postFiltering_123Combined.pdf"), height=10, width=10)
  par(mfrow=c(3,3))
  # Since SMB-ANK-027 had two replicates, we need to plot one vs three and two vs three. First, one vs three
  smoothScatter(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[1]], lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpmFilt)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpmFilt)[which(samplenames %in% allreps[1])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[1]]))$r.squared, digits=2)))
  # best fit/regression line
  abline(lm(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[1]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  # Now, replicates two vs three
  smoothScatter(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[2]], lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpmFilt)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpmFilt)[which(samplenames %in% allreps[1])[2]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[2]]))$r.squared, digits=2)))
  # best fit/regression line
  abline(lm(lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpmFilt[,which(samplenames %in% "SMB-ANK-027")[2]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  for (i in 1:length(allreps)){
    smoothScatter(lcpmFilt[,which(samplenames %in% allreps[i])[1]], lcpmFilt[,which(samplenames %in% allreps[i])[2]], ylab=colnames(lcpmFilt)[which(samplenames %in% allreps[i])[2]], xlab=colnames(lcpmFilt)[which(samplenames %in% allreps[i])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate ", allreps[i], "\n", "r2 =",round(summary(lm(lcpmFilt[,which(samplenames %in% allreps[i])[2]]~lcpmFilt[,which(samplenames %in% allreps[i])[1]]))$r.squared, digits=2)))
    # best fit/regression line
    abline(lm(lcpmFilt[,which(samplenames %in% allreps[i])[2]]~lcpmFilt[,which(samplenames %in% allreps[i])[1]]), col="green")
    # diagonal line
    abline(a=0,b=1,col="red")
  }
dev.off()



#######################################################################
### 5. Other basic QC plots --------------------------------------- ###
#######################################################################

# Visualise library size after filtering with barplots
# IGR NOTE 2019.04.12 - this plot looks identical to one that is supposedly made pre-filtering, but looking at column sums there are differences, just very slim. 

pdf(paste0(edaoutput, "librarySize_indoRNA_postFiltering_123Combined.pdf"), height=10, width=15)
  par(oma=c(2,0,0,0))
  barplot(yFilt$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[yFilt$samples$batch], names=samplenames, las=3, ylim=c(0,max(yFilt$samples$lib.size*1e-6)+10), main="Library Size \n Post-filtering")
  # abline(h=10, col="red")
  legend(x="topright", col=batch.col[unique(yFilt$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

pdf(paste0(edaoutput, "nGenes_indoRNA_postFiltering_123Combined.pdf"), height=10, width=15)
  par(oma=c(2,0,0,0))
  barplot(apply(yFilt$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n Post-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[yFilt$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(yFilt$counts, 2, function(c)sum(c!=0)))+3000))
  legend(x="topright", col=batch.col[unique(yFilt$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

pdf(paste0(edaoutput, "libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf"), height=8, width=15)
  plotDensities(lcpm, group=yFilt$samples$batch, main="unfiltered")
  abline(v=0, lty=3)
  plotDensities(lcpmFilt, group=yFilt$samples$batch, main="log2 CPM > 0 in 50% of village")
  abline(v=0, lty=3)
dev.off()

# get histogram of lcpm
pdf(paste0(edaoutput, "lcpm_postFiltering_Histogram.pdf"), height=10, width=15)
  par(mfrow=c(1,2))
  hist(rowSums(lcpmFilt >= 0.5), main= "Genes expressed at log2 cpm over 0.5 \n post-filtering", xlab="samples", col=4)
  hist(rowSums(lcpmFilt >= 1), main= "Genes expressed at log2 cpm over 1 \n post-filtering", xlab="samples", col=5)
dev.off()


###########################################################################################
### 6. Normalise gene-expression distribution ----------------------------------------- ###
###########################################################################################

# First, compare all three methods
pdf(paste0(edaoutput, "NormalisedGeneExpressionDistribution_IndoRNA_all3Methods.pdf"), height=15, width=15)
  par(oma=c(2,0,0,0), mfrow=c(4,1))
  y2 <- yFilt
  y2$samples$norm.factors <- 1
  lcpm.y2 <- cpm(y2, log=TRUE)
  boxplot(lcpm.y2, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="Unnormalised data",ylab="Log-cpm")
  for (testMethods in c("upperquartile", "TMM", "RLE")) {
    y2 <- calcNormFactors(y2, method=testMethods)
    lcpm.y2 <- cpm(y2, log=TRUE)
    boxplot(lcpm.y2, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
    title(main=paste0("Normalised data ", testMethods),ylab="Log-cpm")
  }
dev.off()

# now just plot how well TMM normalisation worked (since wqe like that one the best)
pdf(paste0(edaoutput, "NormalisedGeneExpressionDistribution_IndoRNA_TMM.pdf"), height=15, width=15)
  par(oma=c(2,0,0,0), mfrow=c(2,1))
  y2 <- yFilt
  y2$samples$norm.factors <- 1
  lcpm.y2 <- cpm(y2, log=TRUE)
  boxplot(lcpm.y2, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="A. Unnormalised data",ylab="Log-cpm")
  y2 <- calcNormFactors(y2, method="TMM")
  lcpm.y2 <- cpm(y2, log=TRUE)
  boxplot(lcpm.y2, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="B. Normalised data, TMM",ylab="Log-cpm")
dev.off()

# save data
save(lcpmFilt, file=paste0(outputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
save(yFilt, file=paste0(outputdir, "indoRNA.read_counts.TMM.filtered.Rda"))

# covariate matrix
save(covariates, file=paste0(outputdir, "covariates.Rda"))


####################################################################################
### 7. Testing batch correction with PCA - moved here from the DE script instead ###
####################################################################################

# PCA visualisation after correction and association with covariates ------------------------------------------------------------

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects was given by the lovely John Blischak.
design <- model.matrix(~0 + yFilt$samples$Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)=gsub("yFilt\\$samples\\$", "", colnames(design))
colnames(design)=gsub("West Papua", "Mappi", colnames(design))

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5) to keep this from crashing.
yFilt$samples$Age[which(is.na(yFilt$samples$Age) == T)]=45

batchCorrectedlcpmFilt <- removeBatchEffect(lcpmFilt, batch=yFilt$samples$batch, covariates = cbind(yFilt$samples$Age, yFilt$samples$RIN, yFilt$samples$CD8T, yFilt$samples$CD4T, yFilt$samples$NK, yFilt$samples$Bcell, yFilt$samples$Mono, yFilt$samples$Gran), design=design)

# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        # points(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], col="black", pch=8, cex=2)
        # text(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=3)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(as.numeric(yFilt$samples$batch)), "topright", pch=unique(as.numeric(yFilt$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

# PCA association function
pc.assoc <- function(pca.data){
    all.pcs <- data.frame()
    for (i in 1:ncol(pca.data$x)){
        all.assoc <- vector()
        for (j in 1:ncol(all.covars.df)){
            test.assoc <- anova(lm(pca.data$x[,i] ~ all.covars.df[,j]))[1,5]
            all.assoc <- c(all.assoc, test.assoc)
        }
        single.pc <- c(i, all.assoc)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", colnames(all.covars.df))

    print ("Here are the relationships between PCs and some possible covariates")
    print (all.pcs)
    return (all.pcs)
}

# Prepare covariate matrix
all.covars.df <- yFilt$samples[,covariate.names]

# Plot PCA
for (name in make.factor){
  if (nlevels(get(name)) < 26){
    pdf(paste0(edaoutput,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpmFilt, speciesCol=as.numeric(get(name)),namesPch=as.numeric(yFilt$samples$batch) + 14,sampleNames=get(name))
    pcaresults <- plot.pca(dataToPca=batchCorrectedlcpmFilt, speciesCol=as.numeric(get(name)),namesPch=as.numeric(yFilt$samples$batch) + 14,sampleNames=get(name))

    dev.off()
  } else {
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpmFilt, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    pcaresults <- plot.pca(dataToPca=batchCorrectedlcpmFilt, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# plot batch
pdf(paste0(edaoutput,"pcaresults_batch.pdf"))
  name <- "batch"
  pcaresults <- plot.pca(dataToPca=lcpmFilt, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(yFilt$samples$batch) + 14,sampleNames=batch)
  pcaresults <- plot.pca(dataToPca=batchCorrectedlcpmFilt, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(yFilt$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# # plot numeric variables
# for (name in make.numeric){
#   initial <- .bincode(get(name), breaks=seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
#   bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
#   pdf(paste0(edaoutput,"pcaresults_",name,".pdf"))
#     pcaresults <- plot.pca(dataToPca=lcpmFilt, speciesCol=bloodCol,namesPch=as.numeric(yFilt$samples$batch) + 14, sampleNames=get(name))
#     legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(get(name))], bloodCol[which.min(get(name))]), cex=0.6, title=name, border=F, bty="n")
#     legend(legend=unique(as.numeric(yFilt$samples$batch)), "topright", pch=unique(as.numeric(yFilt$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
#     pcaresults <- plot.pca(dataToPca=batchCorrectedlcpmFilt, speciesCol=bloodCol,namesPch=as.numeric(yFilt$samples$batch) + 14, sampleNames=get(name))
#     legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(get(name))], bloodCol[which.min(get(name))]), cex=0.6, title=name, border=F, bty="n")
#     legend(legend=unique(as.numeric(yFilt$samples$batch)), "topright", pch=unique(as.numeric(yFilt$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
#   dev.off()
# }

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# # plot pca covariates association matrix to illustrate any remaining confounding/batch
# pdf(paste0(outputdir,"significantCovariates_AnovaHeatmap.pdf"))
#   pheatmap(all.pcs[1:5,c(3:20)], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
#   pheatmap(all.pcs[1:5,c(5:ncol(all.pcs)-2)], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
# dev.off()

# Write out the covariates:
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_blood_RNASeqDeconCell.txt"), col.names=T, row.names=F, quote=F, sep="\t")

