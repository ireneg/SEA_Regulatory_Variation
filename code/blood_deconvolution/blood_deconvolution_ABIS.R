# Test out ABIS for deconvoluting blood using RNASeq data from 123 Indonesian samples
# Code developed by Katalina Bobowik, 20.03.2019
# using the shiny app from the paper from Monaco et al, 2019: https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30059-2


# load packages
library(shiny)
library(dplyr)
library(tidyverse)
library(scater)

# load count data
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")

# setwd
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval")

# ABIS needs the following data and format: Gene expression values in TPM values, gene names must be gene symbols, files must be in tab-delimited format, and files cannot be larger than 5MB
# First, convert expression data to counts per million
cpm <- cpm(y)
# Rownames must be HGNC symbols so we'll use Biomart to convert Ensembl IDs to HGNC symbols
all.gene.IDs=getBM(attributes = c("ensembl_gene_id",'hgnc_symbol'), mart = ensembl.mart.90,values=rownames(cpm), filters="ensembl_gene_id")
new.rownames=all.gene.IDs$hgnc_symbol[match(rownames(cpm),all.gene.IDs$ensembl_gene_id)]
rownames(cpm)=new.rownames
# get rid of emty and NA values (note: filtering doesn't actually do anything to the results of the samples, but it's nice to haev an organised list of actual gene names)
cpm=cpm[-which(rownames(cpm)==""),] %>% .[-which(is.na(rownames(.))),]
dim(cpm)
# 20934   123

# save expression data. Unfortunately, saving the whole df with all decimal points is too large so we have to save it in chunks (with reduced decimal points)
write.table(round(as.data.frame(cpm[,1:50]),2), file="cpm_Samples1to50.txt", sep="\t")

# Now repeat the same thing for the rest of the samples
write.table(round(as.data.frame(cpm[,51:100]),2), file="cpm_Samples51to100.txt", sep="\t")
write.table(round(as.data.frame(cpm[,101:123]),2), file="cpm_Samples101to123.txt", sep="\t")

# Now read in data to ABIS shiny app
install.packages(c("shiny", "MASS"), dependencies = TRUE)
runGitHub("ABIS", user="giannimonaco")

# read in cpm_Samples1to50.txt from launched shiny app and download
# note: got weird error message: Warning in rlm.default(as.matrix(sigRNAseq[genes, ]), x) : 'rlm' failed to converge in 20 steps
ABIS.decon1.50=read.table("cpm_Samples1to50.txt_deconvolution.txt", header=T, row.names=1)
ABIS.decon51.100=read.table("cpm_Samples51to100.txt_deconvolution.txt", header=T, row.names=1)
ABIS.decon101.123=read.table("cpm_Samples101to123.txt_deconvolution.txt", header=T, row.names=1)
# merge all columns together
ABIS.decon=cbind(ABIS.decon1.50,ABIS.decon51.100,ABIS.decon101.123)
ABIS.decon=t(ABIS.decon)

# select relevant cell types - might also need to add  "Monocytes NC+I" in there as well
predicted.cellcounts <- as.data.frame(ABIS.decon[,c("Neutrophils LD","Basophils LD","B Memory","T CD4 Memory","T CD8 Memory","NK","Monocytes C")])
# take out spaces form column names
colnames(predicted.cellcounts)=gsub(" ", ".",colnames(predicted.cellcounts))
# merge Neutrophils LD and Basophils LD and make into Granulocytes
predicted.cellcounts=cbind(Granulocytes=predicted.cellcounts$Neutrophils.LD + predicted.cellcounts$Basophils.LD, predicted.cellcounts)
predicted.cellcounts=predicted.cellcounts[,-c(2,3)]
# scale to sum to 100 (rowSums will add to 100)
predicted.cellcounts.scaled <- (predicted.cellcounts/rowSums(predicted.cellcounts))*100

# barplot
pdf("ABIS_RNASeqDeconvolution.pdf", height=10,width=19)
par(mar=c(14.1,4.1,10.1,2.1),xpd=T)
barplot(t(predicted.cellcounts.scaled), las=3, col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"))
legend(123,130,pch=15, cex=0.8,legend=rownames(t(predicted.cellcounts.scaled)), col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"))
dev.off()

# Compare ABIS deconvolution method to Deconcell ------------------------

# read in deconCell data
decon=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/predictedCellCounts_DeconCell.txt")
# rename decon column names
colnames(decon)=c("Granulocytes", "B.cells", "CD4.T.cells", "CD8.T.cells","NK.cells", "Monocytes")
# rename predicted cell counts (ABIS) df to ABIS (to avoid confusion)
ABIS=as.matrix(predicted.cellcounts.scaled)
rownames(ABIS)=gsub("\\.","-",rownames(ABIS))
identical(rownames(ABIS), rownames(decon))
# TRUE

# now test correlation
pdf("correlation_deconvolution_DeconCellvsABIS.pdf")
# sapply(1:6, function(x) plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
par(mfrow=c(2,3))
for (x in 1:6){
	plot(ABIS[,x],decon[,x], main=paste(colnames(decon)[x],"\n r = ",round(cor(ABIS[,x],decon[,x],method="pearson"),digits=2),sep=""))
	# best fit regression line
	abline(lm(decon[,x]~ABIS[,x]), col="red")
}
dev.off()

# Compare ABIS deconvolution method to Methylation ---------------------------------------

# load in deconvoluted matrix, obtained through methylation data by Heini
norm.meth.decon=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/indonesian_cell_counts_rough_estimate_new.txt", header=T)

# rename rownames in ABIS matrix to match methylation matrix
rownames(ABIS)=samplenames
# match the row ordering of ABIS data to the methylation deconvolution data
ABIS.decon=ABIS[match(norm.meth.decon$ID,rownames(ABIS)),]
# match column ordering
col.order <- c("T.CD8.Memory","T.CD4.Memory","NK","B.Memory","Monocytes.C","Granulocytes")
ABIS.decon=ABIS.decon[,col.order]
# multiply methylation deconvolution data to get numbers on same scale
norm.meth.decon[,c(1:6)]=norm.meth.decon[,c(1:6)]*100

# now test correlation
pdf("correlation_deconvolution_ABISvsMeth_noLogTransform.pdf")
# sapply(1:6, function(x) plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
par(mfrow=c(2,3))
for (x in 1:6){
	plot(ABIS.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],"\n r = ",round(cor(ABIS.decon[,x],norm.meth.decon[,x],method="pearson"),digits=2),sep=""))
	# best fit regression line
	abline(lm(norm.meth.decon[,x]~ABIS.decon[,x]), col="red")
}
dev.off()


# now try using tpm ----------------------------------------

# load in gene lengths df
gene.lengths=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/length_gcContent_allSamples.txt")
y$lengths=gene.lengths[,"Length"][match(rownames(y$counts), rownames(gene.lengths))]
tpm=calculateTPM(y$counts, effective_length=y$lengths)

# Rownames must be HGNC symbols so we'll use Biomart to convert Ensembl IDs to HGNC symbols
new.rownames=all.gene.IDs$hgnc_symbol[match(rownames(tpm),all.gene.IDs$ensembl_gene_id)]
rownames(tpm)=new.rownames

# save expression data. Unfortunately, saving the whole df with all decimal points is too large so we have to save it in chunks (with reduced decimal points)
write.table(round(as.data.frame(tpm[,1:39]),2), file="tpm_Samples1to39.txt", sep="\t")
# Now repeat the same thing for the rest of the samples
write.table(round(as.data.frame(tpm[,40:79]),2), file="tpm_Samples40to79.txt", sep="\t")
write.table(round(as.data.frame(tpm[,80:119]),2), file="tpm_Samples80to119.txt", sep="\t")
write.table(round(as.data.frame(tpm[,120:123]),2), file="tpm_Samples120to123.txt", sep="\t")


# Now read in data to ABIS shiny app
runGitHub("ABIS", user="giannimonaco")

# read in cpm_Samples1to50.txt from launched shiny app and download
# note: got weird error message: Warning in rlm.default(as.matrix(sigRNAseq[genes, ]), x) : 'rlm' failed to converge in 20 steps
ABIS.decon1.39=read.table("tpm_Samples1to39.txt_deconvolution.txt", header=T, row.names=1)
ABIS.decon40.79=read.table("tpm_Samples40to79.txt_deconvolution.txt", header=T, row.names=1)
ABIS.decon80.119=read.table("tpm_Samples80to119.txt_deconvolution.txt", header=T, row.names=1)
ABIS.decon120.123=read.table("tpm_Samples120to123.txt_deconvolution.txt", header=T, row.names=1)

# merge all columns together
ABIS.decon=cbind(ABIS.decon1.39,ABIS.decon40.79,ABIS.decon80.119,ABIS.decon120.123)
ABIS.decon=t(ABIS.decon)

# select relevant cell types - might also need to add  "Monocytes NC+I" in there as well
predicted.cellcounts <- as.data.frame(ABIS.decon[,c("Neutrophils LD","Basophils LD","B Memory","T CD4 Memory","T CD8 Memory","NK","Monocytes C")])
# take out spaces form column names
colnames(predicted.cellcounts)=gsub(" ", ".",colnames(predicted.cellcounts))
# merge Neutrophils LD and Basophils LD and make into Granulocytes
predicted.cellcounts=cbind(Granulocytes=predicted.cellcounts$Neutrophils.LD + predicted.cellcounts$Basophils.LD, predicted.cellcounts)
predicted.cellcounts=predicted.cellcounts[,-c(2,3)]
# scale to sum to 100 (rowSums will add to 100)
predicted.cellcounts.scaled <- (predicted.cellcounts/rowSums(predicted.cellcounts))*100

# barplot
pdf("ABIS_RNASeqDeconvolution_TPM.pdf", height=10,width=19)
par(mar=c(14.1,4.1,10.1,2.1),xpd=T)
barplot(t(predicted.cellcounts.scaled), las=3, col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"))
legend(123,130,pch=15, cex=0.8,legend=rownames(t(predicted.cellcounts.scaled)), col=c("#ffd92f","#e78ac3","#fc8d62","#66c2a5","#8da0cb","#a6d854"))
dev.off()

# Compare ABIS TPM deconvolution method to Deconcell ------------------------

# read in deconCell data
decon=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/predictedCellCounts_DeconCell.txt")
# rename decon column names
colnames(decon)=c("Granulocytes", "B.cells", "CD4.T.cells", "CD8.T.cells","NK.cells", "Monocytes")
# rename predicted cell counts (ABIS) df to ABIS (to avoid confusion)
ABIS=as.matrix(predicted.cellcounts.scaled)
rownames(ABIS)=gsub("\\.","-",rownames(ABIS))
identical(rownames(ABIS), rownames(decon))
# TRUE

# now test correlation
pdf("correlation_deconvolution_DeconCellvsABIS_TPM.pdf")
# sapply(1:6, function(x) plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
par(mfrow=c(2,3))
for (x in 1:6){
	plot(ABIS[,x],decon[,x], main=paste(colnames(decon)[x],"\n r = ",round(cor(ABIS[,x],decon[,x],method="pearson"),digits=2),sep=""))
	# best fit regression line
	abline(lm(decon[,x]~ABIS[,x]), col="red")
}
dev.off()

# Compare ABIS TPM deconvolution method to Methylation ------------------------

# load in deconvoluted matrix, obtained through methylation data by Heini
norm.meth.decon=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/indonesian_cell_counts_rough_estimate_new.txt", header=T)

# rename rownames in ABIS matrix to match methylation matrix
rownames(ABIS)=samplenames
# match the row ordering of ABIS data to the methylation deconvolution data
ABIS.decon=ABIS[match(norm.meth.decon$ID,rownames(ABIS)),]
# match column ordering
col.order <- c("T.CD8.Memory","T.CD4.Memory","NK","B.Memory","Monocytes.C","Granulocytes")
ABIS.decon=ABIS.decon[,col.order]
# multiply methylation deconvolution data to get numbers on same scale
norm.meth.decon[,c(1:6)]=norm.meth.decon[,c(1:6)]*100

# now test correlation
pdf("correlation_deconvolution_ABISvsMeth_TPM.pdf")
# sapply(1:6, function(x) plot(RNA.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],cor(RNA.decon[,x],norm.meth.decon[,x],method="pearson"),sep="\n")))
par(mfrow=c(2,3))
for (x in 1:6){
	plot(ABIS.decon[,x],norm.meth.decon[,x], main=paste(colnames(norm.meth.decon)[x],"\n r = ",round(cor(ABIS.decon[,x],norm.meth.decon[,x],method="pearson"),digits=2),sep=""))
	# best fit regression line
	abline(lm(norm.meth.decon[,x]~ABIS.decon[,x]), col="red")
}
dev.off()



