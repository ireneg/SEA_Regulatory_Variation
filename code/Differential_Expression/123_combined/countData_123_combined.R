# Load in FeatureCounts data, set up covariate matrices, and get initial statistics from first and second batch of RNASeq Indonesian DE analysis
# Code developed by Katalina Bobowik, 06.06.2018


# load packages
library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(Homo.sapiens)
library(limma)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(readr)
library(openxlsx)
library(pheatmap)
library(devtools)
library(ggbiplot)
library(biomaRt)
library(biomartr)
library(gplots)
library(sva)
library(magrittr)
library(dendextend)
library(qvalue)
library(rowr)
library(reshape2)
library(RUVSeq)
library(doParallel)
library(car)
library(ggpubr)
library(GO.db)
library(goseq)
library(ggplot2)
library(ggsignif)
library(wesanderson)
library(treemap)
library(NineteenEightyR)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(vioplot)
library(ReactomePA)


# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
dev.off()
# set up colour palette for batch
batch.col=electronic_night(n=3)

# Set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/countData")

# read in count files from featureCounts. Here, I'm loading in files for all three batches.
files=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba", pattern="Filter", full.names=T)
files.secondbatch=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/second_batch", pattern="Filter", full.names=T)
files.thirdbatch=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/third_batch", pattern="Filter", full.names=T)

# set up DGE matrix combining first, second, and third batch
y <- readDGE(c(files, files.secondbatch, files.thirdbatch), columns=c(1,3)) 
# Organise gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(y)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes

# Trim file names into shorter sample names and apply to column names
colnames(y)[grep("second_batch",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("second_batch",colnames(y))],"[_.]"), `[`, 11), "secondBatch", sep="_")
colnames(y)[grep("third_batch",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("third_batch",colnames(y))],"[_.]"), `[`, 9), "thirdBatch", sep="_")
# make sure to run in this order, as grepping by 'trimmed output' only works after applying the first two commands
colnames(y)[grep("trimmedOutput",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("trimmedOutput",colnames(y))],"[_.]"), `[`, 10), "firstBatch", sep="_")

# sample MPI-336 from batch three is actually sample MPI-381. Let's make sure to change this
colnames(y)[103] <- "MPI-381_thirdBatch"

# create samplenames vector without batch information
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)
# now add a hyphen inbetween the tribe abbreviation and number to keep consistency
samplenames[104:123]=sapply(samplenames[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))

# assign batch
y$samples$batch <- c(rep(1, length(grep("firstBatch",colnames(y)))), rep(2, length(grep("secondBatch",colnames(y)))), rep(3, length(grep("thirdBatch",colnames(y)))))

# Create covariate matrix
covariates = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/metadata_RNA_Batch123.xlsx",sheet=1, detectDates=T)

# add in blood
blood=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/predictedCellCounts_DeconCell.txt", sep="\t", as.is=T, header=T)
colnames(blood)=c("Gran","Bcell","CD4T","CD8T","NK","Mono")
blood$ID=sapply(strsplit(rownames(blood),"[_.]"), `[`, 1)
blood$ID[104:123]=sapply(blood$ID[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))

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

## Get initial statistics before pre-processing

# Visualise library size
pdf("librarysizeIndoRNA_preFiltering.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size")
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf("nGenesIndoRNA_preFilterin.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n pre-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# save unfiltered counts file
saveRDS(y$counts, file = "unfiltered_counts.rds")
