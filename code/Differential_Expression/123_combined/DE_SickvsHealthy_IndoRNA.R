# Differential expression analysis for the Indonesian RNA dataset- analysing sick samples vs healthy samplLook at differential expressino between sick and healthy patients
# Code developed by Katalina Bobowik, 11.03.2019

# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: KB 03.04.2019

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
library(openxlsx)
library(magrittr)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/"
covariatedir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/"
housekeepingdir="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick/RandomSample3/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Load colour schemes:
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)

# BEGIN ANALYSIS -------------------------------------------------------------------------------------------------

# load in count data (as DGE list object)
load(paste0(inputdir, "countData/unfiltered_DGElistObject.Rda"))

# Create shortened samplenames- insert hyphen and drop suffixes
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# Create the covariate matrix --------------------------------------

covariates <- read.xlsx(paste0(covariatedir, "metadata_RNA_Batch123.xlsx"), sheet=1, detectDates=T)

# add in blood
blood=read.table(paste0(inputdir, "batchRemoval/predictedCellCounts_DeconCell.txt"), sep="\t", as.is=T, header=T)
colnames(blood)=c("Gran","Bcell","CD4T","CD8T","NK","Mono")
blood$ID=sapply(strsplit(rownames(blood),"[_.]"), `[`, 1)
blood$ID <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", blood$ID)
# add blood data into covariates matrix
covariates=cbind(covariates, blood[match(covariates$Sample.ID, blood$ID),])

# add in replicate information
covariates$replicate=duplicated(covariates[,1])

# The date column is being finnicky and importing strangely (i.e., some of the dates are recognised as dates, some are not). To fix this, I have to reformat the dates so that they can be handled correctly for the ANOVA analysis (below).
# first let's check which rows in the sampling date column are not in the correct format
unique(covariates[grep("-", covariates$Sampling.Date),6])
# "2016-03-10"

covariates$Sampling.Date=gsub("2016-03-10","10/03/2016",covariates$Sampling.Date) %>% as.Date(., tryFormats="%d/%m/%Y")

# Check if any samples in the covariate DF are not in samplenames
covariates[which(!(covariates[,"Sample.ID"] %in% samplenames)),]
# <0 rows> (or 0-length row.names)
# Nothing! So that's good news

# covariate setup ------------------------------------------------------

# reorder covariates by sequencing batch
covariates=covariates[order(covariates$Sequencing.Batch),]
# make sure all names in files are within the covariate matrix
all(samplenames == covariates[,"Sample.ID"])
# TRUE

# add in sick sample information. We can first read in the sick metadata which was generated in the script "indoRNA_DEAnalysis_STAR_HealthyvsSickMappi.R".
malaria.metadata=read.table(paste0(inputdir, "DE_Sick/Malaria_summary_table.txt"), sep="\t", header=T)
# add in malaria information
covariates$microscopy.pos=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"microscopy"]
covariates$PCR.pos=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"PCR"]
covariates$fract.pfpx.reads=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"fract.reads.pfpx"]

# Assign covariates to DGE list
# subtract variables we don't need
subtract=c("Sample.ID", "Kabisu.Ethnicity","Sequencing.Batch")
# get index of unwanted variables
subtract=which(colnames(covariates) %in% subtract)
# subtract unwanted variables and add in library size and batch variables
covariate.names = colnames(covariates)[-subtract]

# assign covariates to DGE list. Library size and batch are already assigned so subtract these from the loop.
for (name in covariate.names){
  y$samples[[paste0(name)]]<- covariates[,name]
}

# add in library size and batch infromation to covariate names vector
covariate.names = c(covariate.names,"lib.size", "batch")

# Make sure all covariates in y$samples are the appropriate structure. This is important, as it affects the linear model and ANOVA tests
# Get variables that aren't numeric and transform them into factors. Batch is considered numeric (which we want to change), so we'll also add this into our make.factor vector
which.numeric=unlist(lapply(y$samples, is.numeric))
make.factor=c(colnames(y$samples[,!which.numeric]), "batch")
        
# assign variables and make y-samples metadata into appropriate structure
for (name in make.factor){
  if(sum(is.na(y$samples[[name]])) > 0){
    y$samples[[paste0(name)]]=assign(name, as.factor(addNA(y$samples[[paste0(name)]])))
  }
  else{
    y$samples[[paste0(name)]]=assign(name, as.factor(y$samples[[paste0(name)]]))
  }
}

# now do this for numeric variables
which.numeric=unlist(lapply(y$samples, is.numeric))
make.numeric=colnames(y$samples[,which.numeric])[c(1,3:ncol(y$samples[,which.numeric]))]

# assign numeric information to variables
for (name in make.numeric){
  assign(name, as.numeric(as.character(y$samples[[paste0(name)]])))
}

# Age, RIN, and library size need to be broken up into chunks for easier visualisation of trends (for instance in Age, we want to see high age vs low age rather than the effect of every single age variable)
for (name in c("Age","RIN","lib.size")){
  assign(name, cut(as.numeric(as.character(y$samples[[paste0(name)]])), breaks=5))
}

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# Removal of lowly-expressed genes -------------------------------------------------------------------------------

# It looks as though MPI-025, MPI-345, SMB-PTB028 all have high reads mapping to Malaria. We'll identify these as sick
# We'll define sick as any sample with over 100000 reads mapping to either of the plasmodium genomes
sick=colnames(y)[which(y$samples$fract.pfpx.reads>0.05)]
y$samples$DiseaseStatus <- colnames(y) %in% sick
y$samples$DiseaseStatus <- as.factor(y$samples$DiseaseStatus %>% gsub("TRUE", "sick", .) %>% gsub("FALSE", "healthy", .))

# randomly sample 3 healthy Mappi samples
sample.subset=rownames(y$samples[which(y$samples$Island == "West Papua" | y$samples$Island == "Sumba"),] %>% .[which(.$DiseaseStatus == "healthy"),] %>% .[which(.$replicate == FALSE),])
sample.subset=sample(sample.subset, 3)

# subtract all other samples from DGE list
keep.samples=c(sample.subset , sick)
keep.samples=which(colnames(y) %in% keep.samples)
y=y[,keep.samples]
# make sure to drop old levels
y$samples=droplevels(y$samples)
dim(y)
# [1] 27413  6

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at a CPM > 1 in at least half of the libraries
keep.expr.sick <- rowSums(cpm[,which(y$samples$DiseaseStatus!="healthy")]>1) >= (length(which(y$samples$DiseaseStatus!="healthy"))*0.5)
keep.expr.healthy <- rowSums(cpm[,which(y$samples$DiseaseStatus=="healthy")]>1) >= (length(which(y$samples$DiseaseStatus=="healthy"))*0.5)
keep=keep.expr.sick|keep.expr.healthy
y <- y[keep,, keep.lib.sizes=FALSE]
dim(y)
# [1] 12979  6

# Compare library size density before and after removing lowly-expressed genes
pdf(paste0(outputdir,"libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf"), height=8, width=15)
nsamples <- ncol(y)
col=brewer.pal(n=6,"Dark2")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.1), las=2, main="", xlab="")
title(main="A. All genes", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=colnames(y), ncol=1, cex=0.8, text.col=col, bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.1), las=2, main="", xlab="")
title(main="B. Filtered genes", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=colnames(y), ncol=1, cex=0.8, text.col=col, bty="n")
dev.off()

# eliminate composition bias between libraries by upper quartile normalisation (this will be used in the 'data exploration' step)
y <- calcNormFactors(y, method="TMM")

# plot how well TMM normalisation worked (since wqe like that one the best)
pdf(paste0(outputdir,"NormalisedGeneExpressionDistribution_IndoRNA_TMM.pdf"), height=15, width=15)
par(oma=c(2,0,0,0), mfrow=c(2,1))
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75, names=colnames(y))
title(main="A. Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2, method="TMM")
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75, names=colnames(y))
title(main="B. Normalised data, TMM",ylab="Log-cpm")
dev.off()

# recalculate cpm
cpm <- cpm(y)
lcpm=cpm(y, log=T)

# Differential expression ----------------------------------------------------------------------

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
# Set up design matrix. Island information isn't necessary here since we're only testing Mappi
# Note: you need a replicate for each covariate incorporated into the model. I had to take out RIN and blood types for this reason, as they all only had one.
# Alternatively, I can put in sample replicates 
design <- model.matrix(~0 + y$samples$DiseaseStatus + y$samples$Island + y$samples$batch)

# rename columns of design matrix
colnames(design)=gsub("DiseaseStatus", "", colnames(design))
colnames(design)=gsub("[\\$]", "", colnames(design))
colnames(design)=gsub("ysamples", "", colnames(design))
colnames(design)=gsub("IslandWest Papua", "Mappi", colnames(design))

# make contrast matrix
contr.matrix <- makeContrasts(HealthyvsSick=healthy - sick, levels=colnames(design))

# Remove heteroscedascity from count data
pdf(paste0(outputdir,"Limma_voom_TMM_cyclicLoess.pdf"))
v <- voom(y, design, plot=TRUE)
dev.off()

# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
dev.off()


dt <- decideTests(efit,p.value=0.01,lfc=1)
summary(dt)

#        HealthyvsSick
# Down               0
# NotSig         12652
# Up                 0

# test different logFC thresholds
logFC.df=matrix(nrow=3,ncol=1)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    logFC.df[counter,]=sum(abs(dt))
}
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file=paste0(outputdir,"logFC_thresholds.txt"))


