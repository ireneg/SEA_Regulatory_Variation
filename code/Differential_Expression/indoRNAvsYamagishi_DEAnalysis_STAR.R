# Differential expression analysis for the Indonesian RNA dataset
# Code developed by Katalina Bobowik, 9.22.2017


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
library(RSkittleBrewer)
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

# set up  color palette
require(graphics)
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "darkgreen", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
dev.off()

#######################################################
# 1. Read in count data and create covariate matrices #
#######################################################

# Set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/indoRNA/indoRNAvsYamagishi")

# read in count files from featureCounts
files=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba", pattern="Filter", full.names=T)
files.secondbatch=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/second_batch", pattern="Filter", full.names=T)
files.yamagishi=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi", pattern="Filter", full.names=T)

# set up DGE matrix, combining both first and second batch
y <- readDGE(c(files, files.secondbatch, files.yamagishi), columns=c(1,3)) 
# Organise gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(y)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes

# Organise sample names by accession number
colnames(y)[grep("uniquelyMapped",colnames(y),invert=T)] <- sapply(strsplit(colnames(y)[grep("uniquelyMapped",colnames(y),invert=T)],"[_.]"), `[`, 7)
colnames(y)[grep("second_batch",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("second_batch",colnames(y))],"[_.]"), `[`, 11), "secondBatch", sep="_")
colnames(y)[grep("uniquelyMapped",colnames(y))] <- paste(sapply(strsplit(colnames(y)[grep("uniquelyMapped",colnames(y))],"[_.]"), `[`, 10), "firstBatch", sep="_")

# Organise sample names by accession number
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)

# assign batch
y$samples$batch <- c(rep(1, length(grep("firstBatch",colnames(y)))), rep(2, length(grep("secondBatch",colnames(y)))),rep(3, length(grep("DRR",colnames(y)))))

# Create covariate matrix for Indonesian samples
sample_list = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/SampleList_120417_v2d.xlsx",sheet=1)
sequencing_pool = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/Sequencing_multiplexes_HiSeq2500.xlsx", sheet=1)
sequencing_pool_newSamples = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/Sequencing_pools_extra_23.xlsx", sheet=1)
rin = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/FullQCResults_RIN.xlsx", sheet=1)
rin_newSamples = read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/FullQCResults_RIN_23NewSamples.xlsx", sheet=1, colNames=FALSE)
# Reformat data frames so they all match by changing column names and omitting spaces in samplenames
sample_list[,7]=convertToDate(sample_list[,7])
# get rid of unnecessary columns in seq pool df for merge
sample_list=subset(sample_list, select = -c(1,5,9:15,17:23))
sequencing_pool=subset(sequencing_pool, select = -c(2))
colnames(rin)[1]= "Sample.ID"
colnames(sequencing_pool)[1]="Sample.ID"
colnames(sequencing_pool_newSamples)[2]="Sample.ID"
sample_list[,1] = gsub(" ", "-", sample_list[,1])
sequencing_pool_newSamples[,2] = gsub(" ", "-", sequencing_pool_newSamples[,2])
# merge all dataframes together
first_merge = merge(sample_list, sequencing_pool, by="Sample.ID", all=T)
covariates = merge(first_merge, rin, by="Sample.ID", all=T)
# add in RIN and sequencing pool information for second batch of samples
covariates[match(rin_newSamples[,1], covariates[,1]),9] = rin_newSamples[,2]
covariates[match(sequencing_pool_newSamples[,2], covariates[,1]),8] = sequencing_pool_newSamples[,3]
# Check if any samples in the covariate DF are not in samplenames
covariates[which(!(covariates[,1] %in% samplenames)),]

# 69 SMB-ANK-002
# 77 SMB-ANK-012

# remove SMB-ANK-002 and SMB-ANK-012 from covariates dataframe and rename sequencing pool column
covariates = covariates[-c(69,77), ]
colnames(covariates)[8]="Sequencing.pool"

# remove duplicated covariate row
covariates = covariates[-c(which(duplicated(covariates[,1]))),]

# Add in information on healthy and sick samples
sick=c("MPI-025","MPI-345","MPI-334","MPI-376","MPI-061","MPI-296")
covariates$DiseaseStatus <- covariates[,1] %in% sick
covariates$DiseaseStatus <- covariates$DiseaseStatus %>% gsub("TRUE", "sick", .) %>% gsub("FALSE", "healthy", .)

# We know that "MPI-025","MPI-345","MPI-334", and "MPI-376" are all infected with P. Falciparum and "MPI-061" and "MPI-296" with P. Vivax so let's add this info to make it more granular.
# Make rownames the same as samplenames
rownames(covariates)=covariates[,1]
covariates[c("MPI-025","MPI-345","MPI-334","MPI-376"),10]="p.falciparum"
covariates[c("MPI-061","MPI-296"),10]="p.vivax"

## Load in Yamagishi covariate matrices

# Load in Yamagishi et al 2014 supplementary table 11 (with patient information) and SRA run table
sup11=read.delim("/Users/katalinabobowik/Documents/Singapore_StemCells/Projects/Sumba/Papers/SupplementaryMaterials/Yamagishi2014_Supplemental_Table_11.txt", as.is=T, strip.white=T)
sra=read.delim("/Users/katalinabobowik/Documents/Singapore_StemCells/Projects/Sumba/Papers/SupplementaryMaterials/SraRunTable.txt", as.is=T, strip.white=T)

# Create empty matrix and fill in with patient names, Age, Sex, and SRA Run ID
sample_summary=matrix("NA",nrow=nrow(sup11),ncol=5)
sample_summary[,1:3]=c(sup11$Patient_Name,sup11$Age, sup11$Sex)
colnames(sample_summary)=c("Patient_Name", "Age", "Sex", "SRA_ID", "sample_ID")

# Discrepancy in SRA table and Supplementary 11 table (i.e., one says "malaria7#09" and one says "malaria7#009"
sample_summary[103,1]="malaria7#009"

# Match patient names with the patient names in the SRA table, then grab the corresponding SRA run IDs
sample_summary[,4]=sra[match(sample_summary[,1], sra[,4]),7]
sample_summary[,5]=paste("aligned_human_",sample_summary[,4],".sam", sep="")

##########################
# 2. Data Pre-processing #
##########################

# Assign covariates to DGE list
covariate.names = c(colnames(covariates)[c(2,3,6:10)], "lib.size", "batch", "sex")
for (name in covariate.names[1:7]){
  y$samples[[paste0(name)]]<- c(covariates[match(samplenames[1:102], covariates[,1]),name],rep("NA", length(grep("DRR", colnames(y)))))
}

# add in Yamagishi age info
y$samples$Age[grep("DRR", colnames(y))] <- as.factor(sample_summary[match(colnames(y)[grep("DRR", colnames(y))], sample_summary[,4]),2])
y$samples$DiseaseStatus[grep("DRR", colnames(y))] <- "p.falciparum"
y$samples$Island[grep("DRR", colnames(y))] <- "Sulawesi"
y$samples$Sampling.Site[grep("DRR", colnames(y))] <- "Manado"

# Assign gender
y$samples$sex <- sample_summary[match(colnames(y), sample_summary[,4]),3]
y$samples$sex[grep("Batch", colnames(y))] <- "M"

# make variables factors and numeric
y$samples[,c(5:7,9:10,12:13)]=lapply(y$samples[,c(5:7,9:10,12:13)], as.factor)
y$samples[,c(8,11)]=lapply(y$samples[,c(8,11)], as.numeric)

# include NA's into data
y$samples[,c(8:11,13)]=lapply(y$samples[,c(8:11,13)], addNA)

# Filter out samples with library size <8 million
y=y[,which(y$samples$lib.size >= 8000000)]

# reassign samplenames after filtering
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)

# Transform from the raw scale to CPM and log-CPM values (prior count for logCPM = 0.25)
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
# keep.exprs <- rowSums(lcpm>1) >= (nrow(y$samples)*0.5)
keep.expr.sick <- rowSums(lcpm[,which(y$samples$DiseaseStatus!="healthy")]>1) >= (length(which(y$samples$DiseaseStatus!="healthy"))*0.5)
y <- y[keep.expr.sick,, keep.lib.sizes=FALSE]
lcpm <- cpm(y, log=TRUE)
keep.expr.healthy <- rowSums(lcpm[,which(y$samples$DiseaseStatus=="healthy")]>1) >= (length(which(y$samples$DiseaseStatus=="healthy"))*0.5)
y <- y[keep.expr.healthy,, keep.lib.sizes=FALSE]

# After removing lowly expressed genes, assign covariates

# cut Age, RIN, and library size and assign to variables
Age <- addNA(cut(as.numeric(as.character(y$samples$Age)), c(0,14,24,34,44,54,64,74,84), labels=c("0-14","15-24","25-34", "35-44", "45-54", "55-64", "65-74", "75-84")))
RIN <- addNA(cut(as.numeric(as.character(y$samples$RIN)), c(4.9,5.9,6.9,7.9,8.9), labels=c("5.0-5.9", "6.0-6.9", "7.0-7.9", "8.0-8.9")))
lib.size <- cut(as.numeric(y$samples$lib.size), c(50000,8000000,12000000,16000000,20000000,24000000,50000000), labels=c("5e+05-8e+06","8e+06-1.2e+07","1.2e+07-1.6e+07", "1.6e+07-2e+07", "2e+07-2.4e+07","2.4e+07-5e+07"))
        
# assign values to other covariate names
for (name in covariate.names[c(1,2,4,5,7,9,10)]){
  if(sum(is.na(y$samples[[name]])) > 0){
    assign(name, addNA(y$samples[[paste0(name)]]))
  }
  else{
    assign(name, as.factor(y$samples[[paste0(name)]]))
  }
}

# Visualise library size after filtering with barplots
pdf("librarySize_indoRNA_postFiltering.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=y$samples$batch, names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size \n Post-filtering")
abline(h=10, col="red")
legend(x="topright", col=unique(y$samples$batch), legend=c("first batch", "second batch", "Yamagishi"), pch=15, cex=0.8)
dev.off()
pdf("nGenes_indoRNA_postFiltering.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n Post-Filtering", ylab="n Genes", cex.names=0.75, col=y$samples$batch, names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
legend(x="topright", col=unique(y$samples$batch), legend=c("first batch", "second batch", "Yamagishi"), pch=15, cex=0.8)
dev.off()

# Compare library size density before and after removing lowly-expressed genes
pdf("libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf")
nsamples <- ncol(y)
plot(density(lcpm[,1]), col=as.numeric(DiseaseStatus)[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.1), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=as.numeric(DiseaseStatus)[i], lwd=2)
}
legend("topright", legend=c("sick","healthy"), ncol=1, cex=0.8, text.col=unique(as.numeric(DiseaseStatus)), bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=as.numeric(DiseaseStatus)[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=as.numeric(DiseaseStatus)[i], lwd=2)
}
legend("topright", legend=c("sick","healthy"), ncol=1, cex=0.8, text.col=unique(as.numeric(DiseaseStatus)), bty="n")
dev.off()

# Normalise gene expression distributions (i.e., correctino for batch effects)
y <- calcNormFactors(y, method = "TMM")

# Duplicate data, set normalisation back to 1, and plot difference between normalised and non-normalised data
pdf("NormalisedGeneExpressionDistribution_IndoRNA.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=y$samples$batch, main="", cex.axis=0.75, names=samplenames)
title(main="A. Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=y$samples$batch, main="", cex.axis=0.75, names=samplenames)
title(main="B. Normalised data, TMM",ylab="Log-cpm")
dev.off()

# get density plot after normalisation
pdf("densityPlot_NormalisedGeneExpressionDistribution_IndoRNA.pdf")
plot(density(lcpm[,1]), col=y$samples$batch, lwd=2, ylim=c(0,max(density(lcpm)$y)+.2), las=2, main="", xlab="")
title(main="B. Filtered data \n Post-Normalisation", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=y$samples$batch[i], lwd=2)
}
legend("topright", legend=c("First Batch","Second Batch", "Yamagishi"), ncol=1, cex=0.8, text.col=as.numeric(unique(y$samples$batch)), bty="n")
dev.off()

# recalculate lcpm after normalisation
lcpm <- cpm(y, log=TRUE)

######################################
# 4. Covariate and sample analysis #
######################################

# MDS
pdf("indoRNA_MDS.pdf", height=10, width=10)
for (name in covariate.names) {
    plotMDS(lcpm, labels=get(name), col=as.numeric(get(name)))
    title(main=name)
}
dev.off()

# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        points(pca$x[,pca_axis1]["SMB-ANK-027_firstBatch"], pca$x[,pca_axis2]["SMB-ANK-027_firstBatch"], col="black", pch=16, cex=2)
        points(pca$x[,pca_axis1]["SMB-ANK-027_secondBatch"], pca$x[,pca_axis2]["SMB-ANK-027_secondBatch"], col="black", pch=17, cex=2)
        #legend(legend=unique(sampleNames), col=unique(speciesCol), pch=unique(namesPch), x="bottomright", cex=0.6)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(y$samples$batch), "topright", pch=unique(as.numeric(y$samples$batch)) + 15, title="Batch", cex=0.6, border=F, bty="n")
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
all.covars.df <- y$samples[,c(3,5:13)] 

# Plot PCA
for (name in covariate.names){
  if (nlevels(get(name)) < 26){
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=as.numeric(y$samples$batch) + 15,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}


# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf("significantCovariates_AnovaHeatmap.pdf")
pheatmap(all.pcs[1:5,c(2:10)], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file="pca_covariates_blood.txt", col.names=T, row.names=F, quote=F, sep="\t")

# Get relationship of all covariates
new_cov=apply(covariates[2:10], 2, FUN=function(x){x=as.numeric(as.factor(x))})
pdf("covariateHeatmap.pdf", height=10, width=15)
pheatmap(t(new_cov),cluster_col=FALSE,cluster_rows=FALSE, labels_col=covariates[,1], colorRampPalette(rev(brewer.pal(n=10,name="Spectral")))(100),cex=1.1, main="Covariates")
dev.off()

# set colnames to samplenames to conserve space
colnames(lcpm)=make.unique(samplenames)

# Dissimilarity matrix with euclidean distances
pdf("SampleDistances.pdf", height=10, width=15)
par(mar=c(6.1,4.1,4.1,2.1))
eucl.distance <- dist(t(lcpm), method = "euclidean")
eucl.cluster <- hclust(eucl.distance, method = "complete")
dend.eucl=as.dendrogram(eucl.cluster)
labels_colors(dend.eucl)=as.numeric(Island)[order.dendrogram(dend.eucl)]
plot(dend.eucl, main="log2-CPM \n Euclidean Distances")

# Manhattan distance
par(mar=c(6.1,4.1,4.1,2.1))
manht.distance <- dist(t(lcpm), method = "manhattan")
manht.cluster <- hclust(manht.distance, method = "complete" )
dend.manht=as.dendrogram(manht.cluster)
labels_colors(dend.manht)=as.numeric(Island)[order.dendrogram(dend.manht)]
plot(dend.manht, main="log2-CPM \n Manhattan Distances")
dev.off()

# lcpm distances
pdf("lcpmCorrelationHeatmaps.pdf", height=10, width=15)
heatmap.2(cor(lcpm,method="spearman"), labCol=samplenames, margins = c(7, 7), main="Spearman Correlation \n log2-CPM")
dev.off()


######################################################################################
# 6. Differential expression analysis setup and batch effect correction confirmation #
######################################################################################

### set up design

# design matrix
# design <- model.matrix(~0 + batch + Sampling.Site + Island + as.numeric(y$samples$Age) + as.numeric(y$samples$RIN) + DiseaseStatus)

#rename columns to exclude spaces and unrecognised characters
# colnames(design)[11:13]=c("Papua","Age", "RIN") 
# contr.matrix <- makeContrasts(MurrayVsYamagishi=(batch1+batch2)/2 - batch3, levels=colnames(design))


