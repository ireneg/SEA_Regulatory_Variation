# Load in FeatureCounts data, set up covariate matrices, and get initial statistics from first and second batch of RNASeq Indonesian DE analysis
# Code developed by Katalina Bobowik, 06.06.2018

### Last edit: 02.04.2019
### IGR edit read ins, clean package loading list and get rid of hardcoding of paths and indexes 

# load packages
library(RColorBrewer)
library(edgeR)
library(plyr)
library(NineteenEightyR)
# library(Homo.sapiens)

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)

# Set paths:
inputdir <- "/data/cephfs/punim0586/kbobowik/Sumba/" # on server 

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/"
edaoutput <- paste0(outputdir, "/eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

# BEGIN ANALYSIS -----------------------------------------------------------

# read in count files from featureCounts. Here, I'm loading in files for all three batches.
files=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA/sample_counts"), pattern="Filter", full.names=T)
files.secondbatch=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA_second_batch/sample_counts"), pattern="Filter", full.names=T)
files.thirdbatch=list.files(path=paste0(inputdir, "FeatureCounts/indoRNA_third_batch/sample_counts"), pattern="Filter", full.names=T)

# read files into a list
featureCountsOut <- lapply(c(files, files.secondbatch, files.thirdbatch), read.delim)
# transform list into a dataframe
indoReads <- data.frame(t(ldply(featureCountsOut, "[",,3)))

# set col and row names:
names1st <- paste(sapply(strsplit(basename(files), "[_.]"), `[`, 9), "firstBatch", sep="_")
names2nd <- paste(sapply(strsplit(basename(files.secondbatch), "[_.]"), `[`, 9), "secondBatch", sep="_")
names3rd <- paste(sapply(strsplit(basename(files.thirdbatch), "[_.]"), `[`, 7), "thirdBatch", sep="_")

row.names(indoReads) <- featureCountsOut[[1]]$Geneid
names(indoReads) <- c(names1st, names2nd, names3rd) #need to manually confirm ordering --> KB: what do you mean by this?

# rename MPI-336 to MPI-381 (mixup in naming when perfomring sequencing)
colnames(indoReads) <- gsub("MPI-336_thirdBatch", "MPI-381_thirdBatch", colnames(indoReads)) 

# To avoid hardcoding while standardising, fix samplenames object first, then drop the suffixes:
samplenames <- as.character(colnames(indoReads))
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenamesSlim <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# Into DGEList:
y <- DGEList(indoReads, genes=rownames(indoReads), samples=samplenames)

rm(featureCountsOut) # clean up, big object


# assign batch
y$samples$batch <- c(rep(1, length(grep("firstBatch",colnames(y)))), rep(2, length(grep("secondBatch",colnames(y)))), rep(3, length(grep("thirdBatch",colnames(y)))))

# Get initial statistics before pre-processing -----------------------

# Visualise library size
pdf(paste0(edaoutput, "librarysizeIndoRNA_preFiltering.pdf"), height=10, width=15)
    par(oma=c(2,0,0,0))
    barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size")
    legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf(paste0(edaoutput, "nGenesIndoRNA_preFilterin.pdf"), height=10, width=15)
    par(oma=c(2,0,0,0))
    barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n pre-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
    legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# save unfiltered counts file
saveRDS(y$counts, file = paste0(outputdir, "all_samples_unfiltered_counts.rds"))
# save whole DGE list as output for downstream analysis
save(y, file = paste0(outputdir, "all_samples_unfiltered_DGElistObject.Rda"))
