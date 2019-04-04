# Load in FeatureCounts data, set up covariate matrices, and get initial statistics from first and second batch of RNASeq Indonesian DE analysis
# Code developed by Katalina Bobowik, 06.06.2018

### Last edit: 02.04.2019
### IGR edit read ins, clean package loading list and get rid of hardcoding of paths and indexes 

# load packages
# library(Rsubread) #not needed
library(RColorBrewer)
library(edgeR)
library(plyr)
# library(Homo.sapiens) # not needed
# library(limma) (loaded through edgeR)
# library(ensembldb) # not needed
# library(EnsDb.Hsapiens.v75) # not needed
# library(AnnotationDbi) # not needed
library(readr)
library(openxlsx)
# library(pheatmap)
# library(devtools)
library(ggbiplot)
library(biomaRt)
# library(biomartr)
library(gplots)
# library(sva) # not needed
library(magrittr)
library(dendextend)
library(qvalue)
library(rowr)
library(reshape2)
library(RUVSeq)
# library(doParallel)
library(car)
library(ggpubr)
# library(GO.db)
library(goseq)
library(ggplot2)
library(ggsignif)
# library(wesanderson) # not needed
library(treemap)
library(NineteenEightyR)
# library(ComplexHeatmap)
library(circlize)
library(viridis)
library(vioplot)
library(ReactomePA)


# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)

# Set paths:
inputdir <- "/data/cephfs/punim0586/kbobowik/Sumba/"
covariatedir <- "~/"
blooddir <- "/data/cephfs/punim0586/kbobowik/Sumba/Output/DE_Analysis/123_combined/batchRemoval/"

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA_testing/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# BEGIN ANALYSIS

# read in count files from featureCounts. Here, I'm loading in files for all three batches.
files=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA/sample_counts"), pattern="Filter", full.names=T)
files.secondbatch=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA_second_batch/sample_counts"), pattern="Filter", full.names=T)
files.thirdbatch=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA_third_batch/sample_counts"), pattern="Filter", full.names=T)

featureCountsOut <- lapply(c(files, files.secondbatch, files.thirdbatch), read.delim)

indoReads <- data.frame(t(ldply(featureCountsOut, "[",,3)))

# set col and row names:
names1st <- paste(sapply(strsplit(files, "[_.]"), `[`, 10), "firstBatch", sep="_")
names2nd <- paste(sapply(strsplit(files.secondbatch, "[_.]"), `[`, 12), "secondBatch", sep="_")
names3rd <- paste(sapply(strsplit(files.thirdbatch, "[_.]"), `[`, 10), "thirdBatch", sep="_")

row.names(indoReads) <- featureCountsOut[[1]]$Geneid
names(indoReads) <- c(names1st, names2nd, names3rd) #need to manually confirm ordering. 

colnames(indoReads) <- gsub("MPI-336_thirdBatch", "MPI-381_thirdBatch", colnames(indoReads)) #not hard coded. 

# To avoid hardcoding while standardising, fix samplenames object first, then drop the suffixes:
samplenames <- as.character(colnames(indoReads))
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenamesSlim <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# Into DGEList:
y <- DGEList(indoReads, genes=rownames(indoReads), samples=samplenames)

rm(featureCountsOut) # clean up, big object

# assign batch
y$samples$batch <- c(rep(1, length(grep("firstBatch",colnames(y)))), rep(2, length(grep("secondBatch",colnames(y)))), rep(3, length(grep("thirdBatch",colnames(y)))))

# Create covariate matrix
covariates <- read.xlsx(paste0(covariatedir, "metadata_RNA_Batch123.xlsx"), sheet=1, detectDates=T)

# add in blood
blood=read.table(paste0(blooddir, "predictedCellCounts_DeconCell.txt"), sep="\t", as.is=T, header=T)
colnames(blood)=c("Gran","Bcell","CD4T","CD8T","NK","Mono")
blood$ID=sapply(strsplit(rownames(blood),"[_.]"), `[`, 1)
blood$ID <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", blood$ID)

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
pdf(paste0(outputdir, "librarysizeIndoRNA_preFiltering.pdf"), height=10, width=15)
par(oma=c(2,0,0,0))
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size")
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf(paste0(outputdir, "nGenesIndoRNA_preFilterin.pdf"), height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n pre-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# save unfiltered counts file
saveRDS(y$counts, file = paste0(outputdir, "unfiltered_counts.rds"))
