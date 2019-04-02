# by KSB, 09.22.17

# create an RNASeq pipeline for the Yamagishi et al 2014 reads from the paper: Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum. 


library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(Homo.sapiens)
library(limma)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(readr)
library(xlsx)
library(pheatmap)
library(devtools)
library(ggbiplot)
library(biomaRt)
library(biomartr)
library(gplots)
library(sva)
library(topGO)
library(magrittr)
library(dendextend)
library(qvalue)

# read in count files from featureCounts (see readMe "featureCounts")
files=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi", pattern="Filter", full.names=T)
y <- readDGE(files, columns=c(1,3)) #DGEList(counts=y$counts, genes=y$annotation[,c("GeneID","Length")]
# Organising gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(y)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes
# Organise sample names by accession number
samplenames <- sapply(strsplit(colnames(y),"[_.]"), `[`, 7)
colnames(y) <- samplenames

# Visualise library size
cols <- brewer.pal(12,"Set3")
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=cols, names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
abline(h=10, col="red")

# Total number of genes
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=cols, names=colnames(y), las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+2000))

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

# Assign gender and age
y$samples$group <- as.factor(sample_summary[match(colnames(y), sample_summary[,4]),3])
y$samples$age <- sample_summary[match(colnames(y), sample_summary[,4]),2]

#assign covariates to variables
group <- addNA(y$samples$group)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(1000000,10000000,20000000,40000000, 60000000, 80000000), labels=c("1-10", "10-20", "20-40", "40-60", "60-80"))

# assign colors
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group),"Set1")
col.age <- age
levels(col.age) <- brewer.pal(nlevels(col.age),"Set1")
col.lib.size=lib.size
levels(col.lib.size) <- brewer.pal(nlevels(col.lib.size),"Set1")

for (covariate in c("group", "age", "lib.size")) {
    plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=covariate) ## initialize the plot area
    counter=0
    for (sample in colnames(y)){
        counter=counter+1
        lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=get(paste("col",covariate, sep="."))[counter])
    }
    levels=levels(get(covariate))
    levels[which(is.na(levels))] = "NA"
    legend(x="bottomright", bty="n", col=1:length(levels(get(covariate))), legend=levels, lty=1, lwd=2)

 # Filter out samples with library size <10 million and samples with no gender information and assign covariate names
y=y[,which(y$samples$lib.size >= 10000000)]
y=y[,which(y$samples$group != "NA")]

# Visualise library size after filtering
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=cols, names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
# n genes
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=cols, names=colnames(y), las=3)

# Transform from the raw scale to CPM and log-CPM values
# Prior count for logCPM = 0.25
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# get histogram of number of genes expressed at log2 cpm > 0.5 and 1 (before filtering)
hist(rowSums(lcpm>0.5), main= "n Genes expressed at log2 cpm over 0.5 \n pre-filtering", xlab="samples", col=5)
hist(rowSums(lcpm>1), main= "n Genes expressed at log2 cpm over 1 \n pre-filtering", xlab="samples", col=5)

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs <- rowSums(lcpm>1) >= (nrow(y$samples)*0.5)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

# Compare library sizes before and after removing lowly-expressed genes
nsamples <- ncol(y)
col <- colorRampPalette(brewer.pal(nsamples, "Paired")) (nsamples)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")

# get histogram of lcpm after filtering
hist(rowSums(lcpm>0.5), main= "n Genes expressed at log2 cpm over 1 \n post-filtering", xlab="samples", col=5)
hist(rowSums(lcpm>1), main= "n Genes expressed at log2 cpm over 1 \n post-filtering", xlab="samples", col=5)

# Normalise gene expression distributions (i.e., no bias introduced during sample preparation/sequencing)
y <- calcNormFactors(y, method = "TMM")

# Duplicate data, set normalisation back to 1, and plot difference between normalised and non-normalised data
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="B. Example: Normalised data",ylab="Log-cpm")

# get density plot after normalisation
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")

# Reassign covariates
group <- addNA(y$samples$group)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(1000000,10000000,20000000,40000000, 60000000, 80000000), labels=c("1-10", "10-20", "20-40", "40-60", "60-80"))

# Assign colors
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group),"Set1")
col.group <- as.character(col.group)
col.age <- age
levels(col.age) <- brewer.pal(nlevels(col.age),"Set1")
col.age <- as.character(col.age)
col.lib.size=lib.size
levels(col.lib.size) <- brewer.pal(nlevels(col.lib.size),"Set1")
col.lib.size=as.character(col.lib.size)

# MDS
lcpm <- cpm(y, log=TRUE)
plotMDS(lcpm, labels=y$samples$group, col=c("#c51b7d", "#4d9221"))
title(main="Sample Gender")
plotMDS(lcpm, labels=lib.size, col=col.lib.size)
title(main="Sample Library Size")
plotMDS(lcpm, labels=age, col=col.age)
title(main="Sample Age")

# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    plot(pca$x[,1], pca$x[,2], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC1 (", round(pca.var[1]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), main="PCA")
    plot(pca$x[,2], pca$x[,3], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""))
    plot(pca$x[,3], pca$x[,4], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC4 (", round(pca.var[4]*100, digits=2), "% of variance)", sep=""))

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
all.covars.df <- y$samples[,c(2,3,5)] 
all.covars.df$group <- factor(all.covars.df$group)
all.covars.df$age <- factor(all.covars.df$age)

# Plot gender
pcaresults <- plot.pca(lcpm, col.group, 20, colnames(y))

    
# Get PCA associations
# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    plot(pca$x[,1], pca$x[,2], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC1 (", round(pca.var[1]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), main="PCA")
    plot(pca$x[,2], pca$x[,3], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""))
    plot(pca$x[,3], pca$x[,4], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC4 (", round(pca.var[4]*100, digits=2), "% of variance)", sep=""))

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
all.covars.df <- y$samples[,c(2,3,5)] 
all.covars.df$group <- factor(all.covars.df$group)
all.covars.df$age <- factor(all.covars.df$age)

# Plot gender
pcaresults <- plot.pca(lcpm, col.group, 20, colnames(y))

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs
