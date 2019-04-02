# Differential expression analysis for the Indonesian RNA dataset- analysing sick samples vs healthy samplLook at differential expressino between sick and healthy patients
# Code developed by Katalina Bobowik, 11.03.2019

# load packages and set up colour palettes -----------------

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
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick")

# Read in human count data --------------------------------------------------------------------

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
dim(y$genes)
# [1] 27413  3

# Trim file names into shorter sample names and apply to column names
colnames(y)[grep("second_batch",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("second_batch",colnames(y))],"[_.]"), `[`, 11), "secondBatch", sep="_")
colnames(y)[grep("third_batch",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("third_batch",colnames(y))],"[_.]"), `[`, 9), "thirdBatch", sep="_")
colnames(y)[grep("trimmedOutput",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("trimmedOutput",colnames(y))],"[_.]"), `[`, 10), "firstBatch", sep="_")

# sample MPI-336 from batch three is actually sample MPI-381. Let's make sure to change this
colnames(y)[103] <- "MPI-381_thirdBatch"

# create samplenames vector without batch information
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)
samplenames[104:123]=sapply(samplenames[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))

# assign batch
y$samples$batch <- c(rep(1, length(grep("firstBatch",colnames(y)))), rep(2, length(grep("secondBatch",colnames(y)))), rep(3, length(grep("thirdBatch",colnames(y)))))

# Create covariate matrix
covariates = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/metadata_RNA_Batch123.xlsx",sheet=1, detectDates=T)

# add in blood
blood=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/indonesian_cell_counts_rough_estimate_new.txt", sep="\t", as.is=T, header=T)
# blood[,9]=gsub("_", "-", blood[,9])
# blood[,9]=gsub("-new", "", blood[,9])
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

# Read in plasmodium count data and identify sick samples -----------------------------------------

# read in count files from featureCounts for plasmodium falciparum and vivax and assign to DGElist object
files.pf=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/pFalciparum", full.names=T)
files.vx=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/pVivax", full.names=T)

pf <- readDGE(files.pf, columns=c(1,3))
dim(pf)
# [1] 14495   123
vx <- readDGE(files.vx, columns=c(1,3))
dim(vx)
# [1] 17737   123

# rename samplenames
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)

# assign samplenames
samplenames.plas <- sapply(strsplit(colnames(pf),"[_.]"), `[`, 10)
# samplenames and samplenames.plas should be the same so match and make this the new samplename.plas name
samplenames.plas=colnames(y)[match(samplenames.plas, samplenames)]
# sample MPI-336 from batch three is actually sample MPI-381 and SMB-ANK-027_FirstBatch shoudl be changed to second batch
samplenames.plas[c(12,84)] <- c("MPI-381_thirdBatch","SMB-ANK-027_secondBatch")

# assign samplenames as column names
colnames(pf) <- samplenames.plas
colnames(vx) <- samplenames.plas

# keep ordering the same as human data
pf=pf[,order(match(colnames(pf), colnames(y)))]
dim(pf)
#[1] 14495   123
vx=vx[,order(match(colnames(vx), colnames(y)))]
dim(vx)
# [1] 17737   123

# rename samplenames.plas after reordering
samplenames.plas=colnames(pf)

# plot and see which samples are infected
pdf("Plasmodium_LibrarySize.pdf", height=10, width=15)
par(oma=c(5,0,0,0))
barplot(pf$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames.plas, las=3, ylim=c(0,max(pf$samples$lib.size*1e-3)+100), col=Island, main="Number of Raw Reads \n P. Falciparum")
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(vx$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames.plas, las=3, ylim=c(0,max(vx$samples$lib.size*1e-3)+100), col=Island, main="Number of Raw Reads \n P. Vivax")
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf("nGenes_plasmodium.pdf", height=10, width=15)
par(oma=c(5,0,0,0))
barplot(apply(pf$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Falciparum", ylab="n Genes", cex.names=0.75, col=Island, names=samplenames.plas, las=3, ylim=c(0,max(apply(pf$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(apply(vx$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Vivax", ylab="n Genes", cex.names=0.75, col=Island, names=samplenames.plas, las=3, ylim=c(0,max(apply(vx$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# plot normalised count data
# read in unmapped reads file
number.unmapped.reads.df=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/unmappedReads_Counts.txt")
# first 123 rows are duplicated. get rid of these
number.unmapped.reads.df=number.unmapped.reads.df[124:246,]
# set column names
colnames(number.unmapped.reads.df)=c("Sample.ID","Unmapped.Reads")
number.unmapped.reads.df$Sample.ID=gsub("first_batch", "firstBatch", number.unmapped.reads.df$Sample.ID) %>% gsub("second_batch", "secondBatch", .) %>% gsub("third_batch", "thirdBatch", .)
number.unmapped.reads.df$Sample.ID=gsub("MPI-336_thirdBatch", "MPI-381_thirdBatch", number.unmapped.reads.df$Sample.ID)
# match the ordering of plasmodium and unmapped reads df
number.unmapped.reads.df=number.unmapped.reads.df[order(match(number.unmapped.reads.df$Sample.ID, colnames(pf))),]

# First, divide by total number of reads in human data
infected=c("MPI-025_firstBatch", "MPI-334_secondBatch", "MPI-345_secondBatch", "SMB-PTB028_thirdBatch")
pf.val=sort(pf$samples$lib.size/y$samples$lib.size, decreasing=T)
pf.names=samplenames.plas[sort(pf$samples$lib.size/y$samples$lib.size, index.return=T, decreasing=T)$ix]
vx.val=sort(vx$samples$lib.size/y$samples$lib.size, decreasing=T)
vx.names=samplenames.plas[sort(vx$samples$lib.size/y$samples$lib.size, index.return=T, decreasing=T)$ix]

pdf("FractionReadsMappingToPlasmodium_AllHumanReads.pdf", height=10, width=15)
par(oma=c(5,0,0,0))
b=barplot(pf.val,ylab="PF reads/all reads", cex.names=0.75,names=pf.names, las=3, col=Island[sort(pf$samples$lib.size/y$samples$lib.size, index.return=T, decreasing=T)$ix], main="Fraction of Reads Mapping to Plasmodium Falciparum \nUnmapped Reads", ylim=c(0,0.03))
text(b,pf.val,ifelse(pf.names %in% infected,"*",""),pos=3,cex=1.5,xpd=NA)
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
b=barplot(vx.val,ylab="PF reads/all reads", cex.names=0.75,names=vx.names, las=3, col=Island[sort(vx$samples$lib.size/y$samples$lib.size, index.return=T, decreasing=T)$ix], main="Fraction of Reads Mapping to Plasmodium Falciparum \nUnmapped Reads", ylim=c(0,0.1))
text(b,vx.val,ifelse(vx.names %in% infected,"*",""),pos=3,cex=1.5,xpd=NA)
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# Now divide by total number of unmapped reads
pf.val.unmapped=sort(pf$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, decreasing=T)
pf.names.unmapped=samplenames.plas[sort(pf$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, index.return=T, decreasing=T)$ix]
vx.val.unmapped=sort(vx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, decreasing=T)
vx.names.unmapped=samplenames.plas[sort(vx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, index.return=T, decreasing=T)$ix]

pdf("FractionReadsMappingToPlasmodium_UnmappedReads.pdf", height=10, width=15)
par(oma=c(5,0,0,0))
b=barplot(pf.val.unmapped,ylab="PF reads/unmapped reads", cex.names=0.75,names=pf.names.unmapped, las=3, col=Island[sort(pf$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, index.return=T, decreasing=T)$ix], main="Fraction of Reads Mapping to Plasmodium Falciparum \nUnmapped Reads", ylim=c(0,0.4))
text(b,pf.val.unmapped,ifelse(pf.names.unmapped %in% infected,"*",""),pos=3,cex=1.5,xpd=NA)
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
b=barplot(vx.val.unmapped,ylab="VX reads/unmapped reads", cex.names=0.75,names=vx.names.unmapped, las=3, col=Island[sort(vx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads, index.return=T, decreasing=T)$ix], main="Fraction of Reads Mapping to Plasmodium Vivax \nUnmapped Reads", ylim=c(0,0.4))
text(b,vx.val.unmapped,ifelse(vx.names.unmapped %in% infected,"*",""),pos=3,cex=1.5,xpd=NA)
legend(x="topright", col=unique(Island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# make a table of all samples that have malaria, confirmed by PCR and number of reads
samplesheet = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/SampleList_120417_v2d.xlsx",sheet=1, detectDates=T)
samplesheet$Sample.ID=gsub(" ", "-",samplesheet$Sample.ID)

# insert hyphen into samplenames that are merged together
samplenames[104:123]=sapply(samplenames[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))

# assign variables and make DF
sample.ID=samplenames
microscopy=samplesheet$Malaria.microscopy[match(samplenames,samplesheet$Sample.ID)]
PCR=samplesheet$Malaria.PCR[match(samplenames,samplesheet$Sample.ID)]
all.reads.pf=pf$samples$lib.size/y$samples$lib
all.reads.vx=vx$samples$lib.size/y$samples$lib
unmapped.reads.pf=pf$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads
unmapped.reads.vx=vx$samples$lib.size/number.unmapped.reads.df$Unmapped.Reads

malaria.summary=data.frame(sample.ID,microscopy,PCR,all.reads.pf,all.reads.vx,unmapped.reads.pf,unmapped.reads.vx)
write.table(malaria.summary, file="Malaria_summary_table.txt", quote=F, row.names=F, col.names=T, sep="\t")

# covariate setup ------------------------------------------------------

# reorder covariates by sequencing batch
covariates=covariates[order(covariates$Sequencing.Batch),]
# make sure all names in files are within the covariate matrix
all(samplenames == covariates[,1])
# TRUE

# add in sick sample information. We can first read in the sick metadata which was generated in the script "indoRNA_DEAnalysis_STAR_HealthyvsSickMappi.R".
malaria.metadata=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick/Malaria_summary_table.txt", sep="\t", header=T)
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

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)


# Removal of lowly-expressed genes -------------------------------------------------------------------------------

# It looks as though MPI-025, MPI-345, SMB-PTB028 all have high reads mapping to Malaria. We'll identify these as sick
# We'll define sick as any sample with over 100000 reads mapping to either of the plasmodium genomes
sick=colnames(y)[which(y$samples$high.pf.reads>0.01)]
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

# Compare library size density before and after removing lowly-expressed genes
pdf("libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf", height=8, width=15)
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
pdf("NormalisedGeneExpressionDistribution_IndoRNA_TMM.pdf", height=15, width=15)
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
design <- model.matrix(~0 + y$samples$DiseaseStatus + y$samples$Age + y$samples$Island + y$samples$batch)

# rename columns of design matrix
colnames(design)=gsub("DiseaseStatus", "", colnames(design))
colnames(design)=gsub("[\\$]", "", colnames(design))
colnames(design)=gsub("y2samples", "", colnames(design))
colnames(design)[4]="Mappi"

# make contrast matrix
contr.matrix <- makeContrasts(HealthyvsSick=healthy - sick, levels=colnames(design))

# Remove heteroscedascity from count data
pdf("DE/Limma_Voom/Limma_voom_TMM_cyclicLoess.pdf")
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

