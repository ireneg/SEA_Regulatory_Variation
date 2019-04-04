# Test out DeconCell for deconvoluting blood using RNASeq data from 123 Indonesian samples
# Code developed by Katalina Bobowik, 26.02.2019
# following the vignette as per: http://htmlpreview.github.io/?https://github.com/molgenis/systemsgenetics/blob/master/Decon2/DeconCell/inst/doc/my-vignette.html


# load packages
library(devtools)
library(DeconCell)
library(edgeR)
library(tidyverse)
library(ghibli)
library(Rcmdr)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
outputdir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/"
deconestimatedir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/"

# load count data
load(paste0(inputdir, "countData/unfiltered_DGElistObject.Rda"))

# assign count table This is y$counts, which has not been filtered for lowly-expressed genes
count.table=y$counts
dim(y$counts)
# [1] 27413   123

# load the model. Models were built using an elastic net and training in 95 healthy dutch volunteers from the 500FG cohort with FACS quantification of 73 circulating cell subpopulations
data("dCell.models")

# Normalize gene expression: approximate the data to a normal-like distribution and account for library sizes by using the dCell.expProcessing function. This function will perform a TMM normalization, a log2(counts+1), and scale (z-transformation) per gene.
dCell.exp <- dCell.expProcessing(count.table, trim = TRUE)
# [INFO]	 Total of 61.13 % genes from dCell are found>

# predict cell counts
prediction <- dCell.predict(dCell.exp, dCell.models, res.type = "median")

# select relevant cell types
predicted.cellcounts <- prediction$dCell.prediction[,c('Granulocytes','B cells (CD19+)','CD4+ T cells','CD8+ T cells','NK cells (CD3- CD56+)','Monocytes (CD14+)')]

# scale to sum to 100
# predicted.cellcounts.scaled <- (predicted.cellcounts/rowSums(predicted.cellcounts))*100

# save table 
write.table(predicted.cellcounts, file=paste0(outputdir,"predictedCellCounts_DeconCell.txt"), sep="\t")

# plot percentages of each cell type
pdf(paste0(outputdir,"DeconCell_RNASeqDeconvolution.pdf"), height=10,width=19)
par(mar=c(14.1,4.1,10.1,2.1),xpd=T)
barplot(t(predicted.cellcounts), col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"), las=3)
legend(123,130,legend=colnames(predicted.cellcounts), col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"), pch=15, cex=0.8)
dev.off()

# Compare old deconvolution data to new deconvolution data ------------------------

# load in deconvoluted matrix, obtained through methylation data by Heini
norm.meth.decon=read.table(paste0(deconestimatedir,"indonesian_cell_counts_rough_estimate_new.txt"), header=T)
# NOTE: this only has 117 rows, as opposed to 120 in the previous metylation deconvolution one by Irene

# make 'not in' operator:
'%!in%' <- function(x,y)!('%in%'(x,y))

blood=read.table(paste0(deconestimatedir,"indonesian_cell_counts_rough_estimate.txt"), sep="\t", as.is=T, header=T)
blood[,9]=gsub("_", "-", blood[,9])
blood[,9]=gsub("-new", "", blood[,9])
blood=blood[,c(2:7,9)]
blood$Sample.ID[which(duplicated(blood$Sample.ID))]
# [1] "MPI-025"     "SMB-ANK-029" "MPI-296"  
# MPI-025, SMB-ANK-029, and MPI-296 are duplicated. I believe these are replicates since they're in two batches in the metylation metadata

# get rid of duplicates
blood=blood[-c(which(duplicated(blood$Sample.ID))),]
# order blood to be same as new methylation data
blood=blood[match(norm.meth.decon$ID,blood$Sample.ID),]

# Now plot the correlation between each blood type from old deconvolution and new deconvolution data
pdf(paste0(outputdir,"unnormalisedVsNormalised_methDeconvolution.pdf"),height=15,width=10)
par(mfrow=c(3,2))
sapply(1:6, function(x) plot(blood[,x],norm.meth.decon[,x], ylab="normalised_meth", xlab="unnormalised_meth", main=paste(colnames(blood)[x],cor(blood[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
dev.off()

# Compare RNA deconvolution data to methylation deconvolution data ------------------------

# define sample names
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

rownames(predicted.cellcounts)=samplenames
# match the row ordering of the RNA deconvolution data to the methylation deconvolution data
RNA.decon=predicted.cellcounts[match(norm.meth.decon$ID,rownames(predicted.cellcounts)),]
# match column ordering
col.order <- c("CD8+ T cells","CD4+ T cells","NK cells (CD3- CD56+)","B cells (CD19+)","Monocytes (CD14+)","Granulocytes")
RNA.decon=RNA.decon[,col.order]
# multiply methylation deconvolution data to get numbers on same scale
norm.meth.decon[,c(1:6)]=norm.meth.decon[,c(1:6)]*100

# now test correlation
pdf(paste0(outputdir,"correlation_deconvolution_RNAvsMeth.pdf"))
# sapply(1:6, function(x) plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
par(mfrow=c(2,3))
for (x in 1:6){
	plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],"\n r = ",round(cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),digits=2),sep=""))
	# best fit regression line
	abline(lm(norm.meth.decon[,x]~RNA.decon[,x]), col="red")
}
dev.off()


