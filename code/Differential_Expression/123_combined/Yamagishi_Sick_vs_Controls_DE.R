# by KSB, 09.22.17

# create an RNASeq pipeline for the Yamagishi et al 2014 reads from the paper: Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum. 


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

# set colour palette
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
dev.off()

# Set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/Yamagishi")

# read in count files for sick samples
files.sick=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi", pattern="Filter", full.names=T)
files.controls=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi/HealthyControls", pattern="Filter", full.names=T)

# set up DGE matrix combining health and sick samples
y <- readDGE(c(files.sick, files.controls), columns=c(1,3)) 
# Organising gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(y)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes

# assign healthy and control samples
y$samples$diseaseStatus[grep("Controls", colnames(y))]="control"
y$samples$diseaseStatus[grep("Controls", colnames(y), invert = T)]="malaria"
# make disease status into factor
y$samples$diseaseStatus=as.factor(y$samples$diseaseStatus)

# Trim file names into shorter sample names and apply to column names
samplenames <- sapply(1:length(colnames(y)), function(x) tail(strsplit(colnames(y)[x],"_")[[1]],1))
colnames(y) <- samplenames

# Get initial statistics before pre-processing
# Visualise library size
pdf("librarysizeYamagishi_preFiltering.pdf", height=10, width=15)
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
abline(h=6, col="red")
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf("totalGenesYamagishi_preFiltering.pdf", height=10, width=15)
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+2000))
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# Load in Yamagishi et al 2014 supplementary table 11 (with patient information)
sup11.sick=read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/Supplemental_Table_11.xlsx", sheet=1)
sup11.controls=read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/Supplemental_Table_11.xlsx", sheet=3)
sup11.all=rbind(sup11.sick, sup11.controls)

# load in SRA run table info for sick samples and controls
sra.sick=read.delim("/Users/katalinabobowik/Documents/Singapore_StemCells/Projects/Sumba/Papers/SupplementaryMaterials/SraRunTable.txt", as.is=T, strip.white=T)
sra.controls=read.delim("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/SraRunTable_Controls.txt", as.is=T, strip.white=T)
sra.sick=sra.sick[,c("Library_Name_s","Run_s")]
sra.controls=sra.controls[,c("Library_Name","Run")]
colnames(sra.sick)=colnames(sra.controls)
sra.all=rbind(sra.sick,sra.controls)

# Create empty matrix and fill in with patient names, Age, Sex, and SRA Run ID
sample_summary=matrix("NA",nrow=nrow(sup11.all),ncol=7)
sample_summary[,1:5]=c(sup11.all$Patient_Name,sup11.all$Age, sup11.all$Sex, sup11.all[,"%Pf_tags"], sup11.all$From)
colnames(sample_summary)=c("Patient_Name", "Age", "Sex", "PF_Tags","From", "sample_ID", "SRA_ID")
# convert to data frame
sample_summary=as.data.frame(sample_summary)

# There is a discrepancy between the SRA table and Supplementary 11 table- one says "malaria7#09" and one says "malaria7#009". We need to change "malaria7#09" to "malaria7#009"
sample_summary$Patient_Name=gsub("malaria7#09", "malaria7#009", sample_summary$Patient_Name)

# Match patient names with the patient names in the SRA table, then grab the corresponding SRA run IDs
sra.all[match(sample_summary$Patient_Name, sra.all$Library_Name),"Run"]
# we seem to have NA values in the dataframe. Which ones are these?
sample_summary$Patient_Name[which(is.na(match(sample_summary$Patient_Name, sra.all$Library_Name)))]
# [1] "malaria11#1" "malaria11#2" "malaria11#3" "malaria11#4" "malaria11#5"
# [6] "malaria11#6" "malaria11#7" "malaria11#8" "malaria11#9"

# when checking the sra run table sheet, these samples have a 0 in front of them. Let's take this out to keep sample names the same.
sra.all$Library_Name=gsub("malaria11#0","malaria11#", sra.all$Library_Name)
sample_summary[,"SRA_ID"]=sra.all[match(sample_summary$Patient_Name, sra.all$Library_Name),"Run"]

# Assign gender, age, location, and PF load
y$samples$gender <- as.factor(sample_summary[match(colnames(y), sample_summary$SRA_ID),"Sex"])
y$samples$age <- as.numeric(sample_summary[match(colnames(y), sample_summary$SRA_ID),"Age"])
y$samples$location <- as.factor(sample_summary[match(colnames(y), sample_summary$SRA_ID),"From"])
y$samples$PFload <- as.numeric(as.character(sample_summary[match(colnames(y), sample_summary$SRA_ID),"PF_Tags"]))

# let's also add our own PF load made by getting the total number of unmapped reads that mapped to the combined PFPX genome
alignedPFPXreads=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick/Malaria_summary_table_Yamagishi.txt", header=T)
colnames(alignedPFPXreads)[1]="SampleID"
alignedPFPXreads$SampleID=gsub("_Controls","",alignedPFPXreads$SampleID) %>% gsub("_Sick","",.)
identical(alignedPFPXreads$SampleID, colnames(y))
# [1] TRUE
y$samples$alignedPFPX=alignedPFPXreads$fract.reads.pfpx.yam

#assign covariates to variables
gender <- addNA(y$samples$gender)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(0,1000000,10000000,20000000,40000000, 60000000), labels=c("0-1","1-10", "10-20", "20-40", "40-60"))
PFload=addNA(cut(y$samples$PFload, c(-1,20,40,60,80), labels=c("0-20","20-40", "40-60", "60-80")))
diseaseStatus=y$samples$diseaseStatus
location=addNA(y$samples$location)
alignedPFPX=y$samples$alignedPFPX

# assign covariate names
covariate.names=c("gender", "age", "lib.size", "PFload","diseaseStatus","location")

pdf("rarefactionCurves.pdf")
for (name in covariate.names) {
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
dev.off

# Data pre-processing ------------------------------------------------------------------------

# Filter out samples with library size <5 million and samples with no gender information and assign covariate names
y=y[,which(y$samples$lib.size >= 9000000)]
dim(y)
# [1] 27413  124

# Visualise library size after filtering
pdf("librarysizeYamagishi_postFiltering.pdf", height=10, width=15)
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# n genes
pdf("totalGenesYamagishi_postFiltering.pdf", height=10, width=15)
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+2000))
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# get histogram of number of genes expressed at log2 cpm > 0.5 and 1 (before filtering)
hist(rowSums(cpm>0.5), main= "n Genes expressed at cpm > 0.5 \n pre-filtering", xlab="samples", col=4)
hist(rowSums(cpm>1), main= "n Genes expressed at cpm > 1 \n pre-filtering", xlab="samples", col=4)

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs <- rowSums(cpm>1) >= (nrow(y$samples)*0.5)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)
# [1] 11096   124

# Compare library sizes before and after removing lowly-expressed genes
nsamples <- ncol(y)
col <- as.numeric(factor(y$samples$diseaseStatus))
pdf("libraryDensity_afterFiltering.pdf", height=10, width=15)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", levels(gender), ncol=1, cex=0.6, text.col=unique(col), bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")
dev.off()

# Normalise gene expression distributions (i.e., no bias introduced during sample preparation/sequencing)
y <- calcNormFactors(y, method = "TMM")

# Duplicate data, set normalisation back to 1, and plot difference between normalised and non-normalised data
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
pdf("NormalisedGeneExpressionDistribution_YamagishiHealthyvsControls.pdf", height=15, width=15)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

# get density plot after normalisation
pdf("libraryDensity_afterFilteringandNormalisation.pdf", height=10, width=15)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")
dev.off()

# Data exploration --------------------------------------------------------------------

# Reassign covariates
gender <- addNA(y$samples$gender)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(0,1000000,10000000,20000000,40000000, 60000000), labels=c("0-1","1-10", "10-20", "20-40", "40-60"))
PFload=addNA(cut(y$samples$PFload, c(-1,20,40,60,80), labels=c("0-20","20-40", "40-60", "60-80")))
diseaseStatus=y$samples$diseaseStatus
location=addNA(y$samples$location)
alignedPFPX=y$samples$alignedPFPX


# plot MDS
for (name in covariate.names) {
    plotMDS(lcpm, labels=get(name), col=as.numeric(get(name)))
    title(main=name)
}


# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
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
all.covars.df <- y$samples[,c(covariate.names, "alignedPFPX")] 

# Plot PCA
for (name in covariate.names){
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=15,sampleNames=get(name))
    dev.off()
}

# plot numeric variables
name="alignedPFPX"
initial = .bincode(alignedPFPX, breaks=seq(min(alignedPFPX, na.rm=T), max(alignedPFPX, na.rm=T), len = 80),include.lowest = TRUE)
bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
pdf("pcaresults_alignedPlasmodiumLoad.pdf")
pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=bloodCol,namesPch=15,sampleNames=alignedPFPX)
legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(alignedPFPX)], bloodCol[which.min(alignedPFPX)]), cex=0.6, title=name, border=F, bty="n")
dev.off()

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf("significantCovariates_AnovaHeatmap.pdf")
pheatmap(all.pcs[1:5,c(covariate.names,"alignedPFPX")], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file="pca_covariates_significanceLevels.txt", col.names=T, row.names=F, quote=F, sep="\t")

# DE analysis ------------

# Set up design matrix
design <- model.matrix(~0 + y$samples$diseaseStatus + y$samples$alignedPFPX)
colnames(design)=gsub("diseaseStatus", "", colnames(design))
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))

# set up contrast matrix
contr.matrix <- makeContrasts(HealthyvsSick=control - malaria, levels=colnames(design))

# Remove heteroscedascity from count data
pdf("Limma_voom_TMM_cyclicLoess.pdf")
v <- voom(y, design, plot=TRUE)
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
dev.off()


dt <- decideTests(efit,p.value=0.01,lfc=1)
summary(dt)

#        HealthyvsSick
#Down             1407
#NotSig          9613
#Up               974

# look at different DE thresholds
# test different logFC thresholds
logFC.df=matrix(nrow=3,ncol=1)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    logFC.df[counter,]=sum(abs(dt))
}
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file="logFC_thresholds.txt")

# get top genes
toptable=topTable(efit, coef=1, p.value=0.01, n=Inf, sort.by="p")
write.table(toptable,file="topTable_healthyvsSick.txt")


