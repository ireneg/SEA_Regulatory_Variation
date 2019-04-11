# Differential expression analysis for the Indonesian RNA dataset- analysing sick samples vs healthy samplLook at differential expressino between sick and healthy patients
# Code developed by Katalina Bobowik, 11.03.2019

# Load dependencies and set input paths --------------------------

# Load dependencies:
library(edgeR)
library(plyr)
library(openxlsx)
library(RColorBrewer)
library(magrittr)

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
dev.off()

# Set paths:
# inputdir <- "/data/cephfs/punim0586/kbobowik/Sumba/Output/DE_Analysis/123_combined/" # on server
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/" # locally
FeatureCountsDir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/FeatureCounts/"
refdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Read in human count data --------------------------------------------------------------------

# load in count data (as DGE list object)
load(paste0(inputdir, "countData/unfiltered_DGElistObject.Rda"))

# Create shortened samplenames- insert hyphen and drop suffixes
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# load covariate matrix
load(paste0(inputdir, "dataPreprocessing/covariates.Rda"))

# Read in plasmodium count data and identify sick samples for IndoRNA -----------------------------------------

# read in count files from featureCounts for plasmodium falciparum and vivax and assign to DGElist object
files.pfpx=list.files(path=paste0(FeatureCountsDir,"PFPX_Combined_IndoSamples/"), full.names=T)

# read files into a list
featureCountsOut.plasmo <- lapply(files.pfpx, read.delim)
# transform list into a dataframe
plasmoReads.indo <- data.frame(t(ldply(featureCountsOut.plasmo, "[",,3)))

# set col and row names:
names1st <- paste(sapply(strsplit(basename(files.pfpx[grep("first",basename(files.pfpx))]), "[_.]"), `[`, 8), "firstBatch", sep="_")
names2nd <- paste(sapply(strsplit(basename(files.pfpx[grep("second",basename(files.pfpx))]), "[_.]"), `[`, 8), "secondBatch", sep="_")
names3rd <- paste(sapply(strsplit(basename(files.pfpx[grep("third",basename(files.pfpx))]), "[_.]"), `[`, 8), "thirdBatch", sep="_")

row.names(plasmoReads.indo) <- featureCountsOut.plasmo[[1]]$Geneid
names(plasmoReads.indo) <- c(names1st, names2nd, names3rd) #need to manually confirm ordering --> KB: what do you mean by this?

# rename MPI-336 to MPI-381 (mixup in naming when perfomring sequencing)
colnames(plasmoReads.indo) <- gsub("MPI-336_thirdBatch", "MPI-381_thirdBatch", colnames(plasmoReads.indo)) 

# Into DGEList:
pfpx <- DGEList(plasmoReads.indo, genes=rownames(plasmoReads.indo), samples=colnames(plasmoReads.indo))
dim(pfpx)
# [1] 32232   123

# Create shortened samplenames- insert hyphen and drop suffixes
samplenames.plas <- as.character(pfpx$samples$samples)
samplenames.plas <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames.plas)
samplenames.plas <- sapply(strsplit(samplenames.plas, "[_.]"), `[`, 1)

rm(featureCountsOut.plasmo) # clean up, big object

# plot and see which samples are infected
pdf(paste0(outputdir,"Plasmodium_LibrarySize.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(pfpx$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames.plas, col=as.numeric(as.factor(covariates$Island)), las=3, ylim=c(0,max(pfpx$samples$lib.size*1e-3)+200), main="Number of Raw Reads \nCombined Falc & Vivax")
legend(x="topright", col=unique(as.numeric(as.factor(covariates$Island))), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf(paste0(outputdir,"nGenes_plasmodium.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(apply(pfpx$counts, 2, function(c)sum(c!=0)),main="Number of Genes \nCombined Falc & Vivax", ylab="n Genes", cex.names=0.75, col=as.numeric(as.factor(covariates$Island)), names=samplenames.plas, las=3, ylim=c(0,max(apply(pfpx$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(as.numeric(as.factor(covariates$Island))), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# plot normalised count data
# read in unmapped reads file
number.unmapped.reads.df=read.table(paste0(refdir,"unmappedReads_Counts.txt"))
# first 123 rows are duplicated. get rid of these
number.unmapped.reads.df=number.unmapped.reads.df[124:246,]
# set column names
colnames(number.unmapped.reads.df)=c("Sample.ID","Unmapped.Reads")
number.unmapped.reads.df$Sample.ID=gsub("first_batch", "firstBatch", number.unmapped.reads.df$Sample.ID) %>% gsub("second_batch", "secondBatch", .) %>% gsub("third_batch", "thirdBatch", .)
# rename MPI-336 to MPI-381 (mixup in naming when perfomring sequencing)
number.unmapped.reads.df$Sample.ID=gsub("MPI-336_thirdBatch", "MPI-381_thirdBatch", number.unmapped.reads.df$Sample.ID)
# See if the total unmapped reads and PFPX file are the same
identical(number.unmapped.reads.df$Sample.ID, colnames(pfpx))
# [1] TRUE

pdf(paste0(outputdir,"FractionReadsMappingToPlasmodium_NormalisedByUnmappedReads.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(pfpx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, ylab="Fraction of Total Unmapped Reads", cex.names=0.75,names=samplenames.plas, col=as.numeric(as.factor(covariates$Island)), las=3, main="Fraction of Combined Falc & Vivax Reads\n Normalised By Total Reads")
dev.off()

# make a table of all samples that have malaria, confirmed by PCR and number of reads
samplesheet = read.xlsx(paste0(refdir,"/indoRNA_SequencingFiles/SampleList_120417_v2d.xlsx"),sheet=1, detectDates=T)
samplesheet$Sample.ID=gsub(" ", "-",samplesheet$Sample.ID)

# assign variables and make DF
sample.ID=samplenames
microscopy=samplesheet$Malaria.microscopy[match(samplenames,samplesheet$Sample.ID)]
PCR=samplesheet$Malaria.PCR[match(samplenames,samplesheet$Sample.ID)]
fract.reads.pfpx=pfpx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads

malaria.summary=data.frame(sample.ID,microscopy,PCR,fract.reads.pfpx)
write.table(malaria.summary, file=paste0(outputdir,"Malaria_summary_table.txt"), quote=F, row.names=F, col.names=T, sep="\t")


# Identify sick Yamagishi samples ------------------------------

# read in count files from featureCounts for plasmodium falciparum and vivax and assign to DGElist object
pfpx.yam.sick=list.files(path=paste0(FeatureCountsDir,"PFPX_Combined_Yamagishi/Sick/"), full.names=T)
pfpx.yam.controls=list.files(path=paste0(FeatureCountsDir,"PFPX_Combined_Yamagishi/Controls/"), full.names=T)

# read files into a list
featureCountsOut.yam.all <- lapply(c(pfpx.yam.sick,pfpx.yam.controls), read.delim)

# transform list into a dataframe
plasmoReads.yam <- data.frame(t(ldply(featureCountsOut.yam.all, "[",,3)))

# set col and row names:
namesSick <- paste(sapply(strsplit(basename(pfpx.yam.sick), "[_.]"), `[`, 10), "Sick", sep="_") %>% gsub("Aligned", "", .)
namesControls <- paste(sapply(strsplit(basename(pfpx.yam.controls), "[_.]"), `[`, 10), "Controls", sep="_") %>% gsub("Aligned", "", .)

row.names(plasmoReads.yam) <- featureCountsOut.yam.all[[1]]$Geneid
names(plasmoReads.yam) <- c(namesSick, namesControls) #need to manually confirm ordering --> KB: what do you mean by this?

# Into DGEList:
pfpx.yam <- DGEList(plasmoReads.yam, genes=rownames(plasmoReads.yam), samples=colnames(plasmoReads.yam))
dim(pfpx.yam)
# [1] 32232   147

rm(featureCountsOut.yam.all) # clean up, big object

# set colours for sick and control samples
cols=c(rep(7,length(grep("Sick", colnames(pfpx.yam)))), rep(8,length(grep("Controls", colnames(pfpx.yam)))))

# plot and see which samples are infected
pdf(paste0(outputdir,"Plasmodium_LibrarySize.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(pfpx.yam$samples$lib.size, ylab="Library size (millions)", cex.names=0.75,names=colnames(pfpx.yam), col=cols, las=3, main="Number of Raw Reads \nCombined Falc & Vivax")
legend(x="topright", col=unique(cols), legend=c("Sick","Control"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf(paste0(outputdir,"nGenes_plasmodium.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(apply(pfpx$counts, 2, function(c)sum(c!=0)),main="Number of Genes\nCombined Falc & Vivax", ylab="n Genes", cex.names=0.75, col=cols, names=samplenames.plas, las=3, ylim=c(0,max(apply(pfpx$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(as.numeric(as.factor(covariates$Island))), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# plot normalised count data
# read in unmapped reads file
number.unmapped.reads.df.yam=read.table(paste0(refdir,"unmappedReads_Counts_Yamagishi.txt"))
# set column names
colnames(number.unmapped.reads.df.yam)=c("Sample.ID","Unmapped.Reads")
# Reorder the unmapped reads df
number.unmapped.reads.df.yam=number.unmapped.reads.df.yam[match(colnames(pfpx.yam), number.unmapped.reads.df.yam$Sample.ID),]
identical(colnames(pfpx.yam), as.character(number.unmapped.reads.df.yam$Sample.ID))

pdf(paste0(outputdir,"FractionReadsMappingToPlasmodium_NormalisedByUnmappedReads.pdf"), height=10, width=15)
par(oma=c(5,0,0,0))
barplot(pfpx.yam$samples$lib.size/number.unmapped.reads.df.yam$Unmapped.Reads, ylab="Fraction of Total Unmapped Reads", cex.names=0.75,names=colnames(pfpx.yam), col=cols, las=3, main="Fraction of Combined Falc & Vivax Reads\n Normalised By Total Reads")
dev.off()

fract.reads.pfpx.yam=pfpx.yam$samples$lib.size/number.unmapped.reads.df.yam$Unmapped.Reads
malaria.summary.yam=data.frame(colnames(pfpx.yam),fract.reads.pfpx.yam)
write.table(malaria.summary.yam, file=paste0(outputdir,"Malaria_summary_table_Yamagishi.txt"), quote=F, row.names=F, col.names=T, sep="\t")


