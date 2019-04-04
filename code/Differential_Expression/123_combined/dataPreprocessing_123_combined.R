# script created by KSB, 26.07.18
# perform data pre-processing on both batches of indonesian RNA-seq data- third batch

### Last edit: KB 03.04.2019

# Load dependencies and set input paths --------------------------

# Load dependencies:
library(edgeR)
library(plyr)
library(NineteenEightyR)
library(openxlsx)
library(RColorBrewer)

# Set paths:
# inputdir <- "/data/cephfs/punim0586/kbobowik/Sumba/Output/DE_Analysis/123_combined/" # on server
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/" # locally
covariatedir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)

# BEGIN ANALYSIS -----------------------------------------------------------

# load in count data (as DGE list object)
y=readRDS(paste0(inputdir, "countData/unfiltered_DGElistObject.rds"))

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
covariates[,14:19]=blood[match(covariates$Sample.ID, blood$ID),1:6]

# add in replicate information
covariates$replicate=duplicated(covariates[,1])

# The date column is being finnicky and importing strangely (i.e., some of the dates are recognised as dates, some are not). To fix this, I have to reformat the dates in rows 39:70, then change this to a date variable so that it can be handled correctly for the ANOVA analysis (below).
covariates[39:70, 6]="10/03/2016"
covariates[,6]=as.Date(covariates[,6], tryFormats="%d/%m/%Y")

## Check if any samples in the covariate DF are not in samplenames
covariates[which(!(covariates[,1] %in% samplenames)),]
# <0 rows> (or 0-length row.names)
# Nothing! So that's good news

# covariate setup ------------------------------------------------------

# reorder covariates by sequencing batch
covariates=covariates[order(covariates$Sequencing.Batch),]
# make sure all names in files are within the covariate matrix
all(samplenames == covariates[,1])
# TRUE

# add in sick sample information. We can first read in the sick metadata which was generated in the script "indoRNA_DEAnalysis_STAR_HealthyvsSickMappi.R".
malaria.metadata=read.table(paste0(inputdir, "DE_Sick/Malaria_summary_table.txt"), sep="\t", header=T)
# add in malaria information
covariates$microscopy.pos=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"microscopy"]
covariates$PCR.pos=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"PCR"]
covariates$high.pf.reads=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"all.reads.pf"]
covariates$high.vx.reads=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"all.reads.vx"]
covariates$high.pf.reads.unmapped=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"unmapped.reads.pf"]
covariates$high.vx.reads.unmapped=malaria.metadata[match(covariates$Sample.ID, malaria.metadata$sample.ID),"unmapped.reads.vx"]

# Assign covariates to DGE list
# subtract variables we don't need
subtract=c("Sample.ID", "Kabisu.Ethnicity","Sequencing.Batch")
# get index of unwanted variables
subtract=which(colnames(covariates) %in% subtract)
# subtract unwanted variables and add in library size and batch variables
covariate.names = c(colnames(covariates)[-subtract], "lib.size", "batch")
# assign covariates to DGE list. Library size and batch are already assigned so subtract these from the loop.
for (name in head(covariate.names, -2)){
  y$samples[[paste0(name)]]<- covariates[,name]
}

# Make sure all covariates in y$samples are the appropriate structure. This is important, as it affects the linear model and ANOVA tests
# Get variables that aren't numeric and transform them into factors. Batch is considered numeric (which we want to change), so we'll also add this into our make.factor vector
which.numeric=unlist(lapply(y$samples, is.numeric))
make.factor=c(colnames(y$samples[,!which.numeric])[3:ncol(y$samples[,!which.numeric])], "batch")
        
# assign variables and make y-samples metadata into appropriate structure
for (name in make.factor){
  if(sum(is.na(y$samples[[name]])) > 0){
    y$samples[[paste0(name)]]=assign(name, as.factor(addNA(y$samples[[paste0(name)]])))
  }
  else{
    y$samples[[paste0(name)]]=assign(name, as.factor(y$samples[[paste0(name)]]))
  }
}

## Covariate assignment: below, we want to assign all of our covariates to variables.
# Age, RIN, and library size need to be broken up into chunks for easier visualisation of trends (for instance in Age, we want to see high age vs low age rather than the effect of every single age variable)
Age <- cut(as.numeric(as.character(y$samples$Age)), c(14,24,34,44,54,64,74,84), labels=c("15-24","25-34", "35-44", "45-54", "55-64", "65-74", "75-84"))
RIN <- cut(as.numeric(as.character(y$samples$RIN)), c(4.9,5.9,6.9,7.9,8.9), labels=c("5.0-5.9", "6.0-6.9", "7.0-7.9", "8.0-8.9"))
lib.size <- cut(as.numeric(y$samples$lib.size), c(8000000,12000000,16000000,20000000,24000000), labels=c("8e+06-1.2e+07","1.2e+07-1.6e+07", "1.6e+07-2e+07", "2e+07-4.4e+07"))

# assign blood type information to variables
for (name in covariate.names[c(11:16)]){
assign(name, as.numeric(as.character(y$samples[[paste0(name)]])))
}

# assign PF and VX load information
for (name in covariate.names[grep("reads", covariate.names)]){
  assign(name, as.numeric(as.character(y$samples[[paste0(name)]])))
}

# Perform rarefaction curves for number of expressed genes vs. proportion of pool mRNA (as in Ramskold D, Wang ET, Burge CB, Sandberg R. 2009. An abundance of ubiquitously expressed genes revealed by tissue transcriptome sequence data. PLoS Comput Biol 5:e1000598)
pdf(paste0(outputdir, "rarefactionCurves_indoRNA_123Combined.pdf"))
for (name in colnames(Filter(is.factor,y$samples))[-c(1,2)]) {
  plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=name) ## initialize the plot area
  counter=0
  for (sample in colnames(y)){
    counter=counter+1
    lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=as.numeric(get(name))[counter])
  }
  levels=levels(get(name))
  levels[which(is.na(levels))] = "NA"
  legend(x="bottomright", bty="n", col=1:length(levels(get(name))), legend=levels, lty=1, lwd=2)
}
# add in batch (we do this separately since it has a different colour)
plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main="batch") ## initialize the plot area
counter=0
for (sample in colnames(y)){
  counter=counter+1
  lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=batch.col[as.numeric(batch)][counter])
}
levels=levels(batch)
levels[which(is.na(levels))] = "NA"
legend(x="bottomright", bty="n", col=batch.col[1:length(levels(batch))], legend=levels, lty=1, lwd=2)
dev.off()

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# replicate analysis -----------------------------------------------------------------------------------

# set replicate names and rename lcpm columns
allreps=covariates[,1][which(covariates$replicate)]
allreps=unique(allreps)

# look at how well replicate performed
pdf(paste0(outputdir, "replicate_comparisons_noFiltering_123Combined.pdf"), height=10, width=10)
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

# Removal of lowly-expressed genes -------------------------------------------------------------------------------

# get histogram of number of genes expressed at log2 cpm > 0.5 and 1 (before filtering)
pdf(paste0(outputdir, "lcpm_preFiltering_Histogram.pdf_cpm1.pdf"), height=10, width=15)
par(mfrow=c(1,2))
hist(rowSums(lcpm>0.5), main= "Genes expressed at log2 cpm over 0.5 \n pre-filtering", xlab="samples", col=4, ylim=c(0,16000))
hist(rowSums(lcpm>1), main= "Genes expressed at log2 cpm over 1 \n pre-filtering", xlab="samples", col=5, ylim=c(0,16000))
dev.off()

# a gene is only retained if it is expressed at a CPM > 1 in at least half of the group (i.e., in half of MPI, half of SMB, or half of MTW)

# Extract CPM information for each village
cpm_SMB <- cpm[,grep("SMB",colnames(cpm))]
cpm_MPI <- cpm[,grep("SMB",colnames(cpm))]
cpm_MTW <- cpm[,grep("MTW",colnames(cpm))]

# set keep threshold to only retain genes that are expressed at or over a CPPM of one in at least half of the group
keep.SMB <- rowSums(cpm_SMB>1) >= (length(grep("SMB",rownames(y$samples)))*0.5)
keep.MPI <- rowSums(cpm_MPI>1) >= (length(grep("MPI",rownames(y$samples)))*0.5)
keep.MTW <- rowSums(cpm_MTW>1) >= (length(grep("MTW",rownames(y$samples)))*0.5)

# Keep genes that have a CPM >= 1 in at least half of the samples within one of the comparison groups
keep <- keep.SMB | keep.MPI | keep.MTW
y <- y[keep,, keep.lib.sizes=FALSE]
dim(y)
# [1] 12975   123

# Visualise library size after filtering with barplots
pdf(paste0(outputdir, "librarySize_indoRNA_postFiltering_123Combined.pdf"), height=10, width=15)
  par(oma=c(2,0,0,0))
  barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size \n Post-filtering")
  # abline(h=10, col="red")
  legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

pdf(paste0(outputdir, "nGenes_indoRNA_postFiltering_123Combined.pdf"), height=10, width=15)
  par(oma=c(2,0,0,0))
  barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n Post-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
  legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# Compare library size density before and after removing lowly-expressed genes
pdf(paste0(outputdir, "libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf"), height=8, width=15)
  nsamples <- ncol(y)
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=batch.col[y$samples$batch][1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.1), las=2, main="", xlab="")
  title(main="A. All genes", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=batch.col[y$samples$batch][i], lwd=2)
  }
  legend("topright", legend=c("First Batch","Second Batch", "Third Batch"), ncol=1, cex=0.8, text.col=batch.col[unique(y$samples$batch)], bty="n")

  lcpm <- cpm(y, log=TRUE) # why are you defining this again in the middle of the plotting function?
  plot(density(lcpm[,1]), col=batch.col[y$samples$batch][1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.2), las=2, main="", xlab="")
  title(main="B. Filtered genes", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=batch.col[y$samples$batch][i], lwd=2)
  }
  legend("topright", legend=c("First Batch","Second Batch", "Third Batch"), ncol=1, cex=0.8, text.col=batch.col[unique(y$samples$batch)], bty="n")
dev.off()

# get histogram of lcpm
pdf(paste0(outputdir, "lcpm_postFiltering_Histogram.pdf"), height=10, width=15)
  par(mfrow=c(1,2))
  hist(rowSums(lcpm>0.5), main= "Genes expressed at log2 cpm over 0.5 \n post-filtering", xlab="samples", col=4)
  hist(rowSums(lcpm>1), main= "Genes expressed at log2 cpm over 1 \n post-filtering", xlab="samples", col=5)
dev.off()

# replicate analysis after removal of lowly-expressed genes ----------------------------------------------------------------

# get correlation of technical replicates after filtering for lowly expressed genes
pdf(paste0(outputdir, "replicate_comparisons_postFiltering_123Combined.pdf"), height=10, width=10)
  par(mfrow=c(3,3))
  # Since SMB-ANK-027 had two replicates, we need to plot one vs three and two vs three. First, one vs three
  smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[1]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]))$r.squared, digits=2)))
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  # Now, replicates two vs three
  smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[2]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[2]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]))$r.squared, digits=2)))
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")

  for (i in 1:length(allreps)){
    smoothScatter(lcpm[,which(samplenames %in% allreps[i])[1]], lcpm[,which(samplenames %in% allreps[i])[2]], ylab=colnames(lcpm)[which(samplenames %in% allreps[i])[2]], xlab=colnames(lcpm)[which(samplenames %in% allreps[i])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate ", allreps[i], "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]))$r.squared, digits=2)))
    # best fit/regression line
    abline(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]), col="green")
    # diagonal line
    abline(a=0,b=1,col="red")
  }
dev.off()

# normalise gene-expression distribution -------------------------------------------------

# eliminate composition bias between libraries by upper quartile normalisation (this will be used in the 'data exploration' step)
y <- calcNormFactors(y, method="TMM")

# look at performance of normalisation
# First, compare all three methods
pdf(paste0(outputdir, "NormalisedGeneExpressionDistribution_IndoRNA_all3Methods.pdf"), height=15, width=15)
  par(oma=c(2,0,0,0), mfrow=c(4,1))
  y2 <- y
  y2$samples$norm.factors <- 1
  lcpm <- cpm(y2, log=TRUE)
  boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="Unnormalised data",ylab="Log-cpm")
  for (method in c("upperquartile", "TMM", "RLE")) {
    y2 <- calcNormFactors(y2, method="upperquartile")
    lcpm <- cpm(y2, log=TRUE)
    boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
    title(main=paste0("Normalised data ", method),ylab="Log-cpm")
  }
dev.off()

# now just plot how well TMM normalisation worked (since wqe like that one the best)
pdf(paste0(outputdir, "NormalisedGeneExpressionDistribution_IndoRNA_TMM.pdf"), height=15, width=15)
  par(oma=c(2,0,0,0), mfrow=c(2,1))
  y2 <- y
  y2$samples$norm.factors <- 1
  lcpm <- cpm(y2, log=TRUE)
  boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="A. Unnormalised data",ylab="Log-cpm")
  y2 <- calcNormFactors(y2, method="TMM")
  lcpm <- cpm(y2, log=TRUE)
  boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main="B. Normalised data, TMM",ylab="Log-cpm")
dev.off()

# recalculate lcpm
lcpm=cpm(y, log=T)

# save data
save(lcpm, file=paste0(outputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
save(y, file=paste0(outputdir, "indoRNA.read_counts.TMM.filtered.Rda"))
# covariate matrix
save(covariates, file=paste0(outputdir, "covariates.Rda"))


