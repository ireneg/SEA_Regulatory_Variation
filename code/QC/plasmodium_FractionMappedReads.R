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
library(topGO)
library(magrittr)
library(dendextend)
library(qvalue)
library(rowr)
library(reshape2)
library(RSkittleBrewer)
library(RUVSeq)
library(variancePartition)
library(doParallel)

# set up  color palette
require(graphics)
tropical=RSkittleBrewer('tropical')
palette(c(tropical, brewer.pal(8,"Dark2")))
dev.off()

#########################################
# 1. Read in p. Falciparum count data #
#########################################

# Set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/MappedPlasmodiumReads")

# read in count files from featureCounts for plasmodium falciparum and vivax
files.pf=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/pFalciparum", pattern="Filter", full.names=T)
files.vx=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Sumba/pVivax", pattern="Filter", full.names=T)

pf <- readDGE(files.pf, columns=c(1,3))
vx <- falciparum <- readDGE(files.vx, columns=c(1,3))

# assign samplenames (both PF and VX DGE list objects are the same)
samplenames.pf <- sapply(strsplit(colnames(pf),"[_.]"), `[`, 10)

# Organise gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
# geneid <- rownames(y)
# geneid <- gsub("exon_", "", geneid)
# geneid <- sapply(strsplit(geneid,"[-.]"), `[`, 1)
# genes <- select(org.Pf.plasmo.db, keys=geneid, columns=c("ORF", "GENENAME"), keytype="SYMBOL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
# genes <- genes[!duplicated(genes$SYMBOL),]
# y$genes <- genes

# assign colours to island and visualise library size
island=c(rep(1,length(grep("MPI", samplenames))), rep(2,length(grep("MTW", samplenames))), rep(3,length(grep("SMB", samplenames))))
pdf("librarysize_plasmodium.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(pf$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames, las=3, ylim=c(0,max(pf$samples$lib.size*1e-3)+100), col=island, main="P. Falciparum")
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(vx$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames, las=3, ylim=c(0,max(vx$samples$lib.size*1e-3)+100), col=island, main="P. Vivax")
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf("nGenes_plasmodium.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(pf$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Falciparum", ylab="n Genes", cex.names=0.75, col=island, names=samplenames, las=3, ylim=c(0,max(apply(pf$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(apply(vx$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Vivax", ylab="n Genes", cex.names=0.75, col=island, names=samplenames, las=3, ylim=c(0,max(apply(vx$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# filter out lowly expressed genes

falc <- rowSums(pf$counts)
isexpr.pf <- falc > 50
pf <- pf[isexpr, , keep.lib.size = FALSE]
viv <- rowSums(vx$counts)
isexpr.vx <- viv > 50
vx <- vx[isexpr, , keep.lib.size = FALSE]

# get library size after filterig
pdf("librarysize_plasmodium_postFiltering.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(pf$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames, las=3, ylim=c(0,max(pf$samples$lib.size*1e-3)+100), col=island, main="P. Falciparum")
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(vx$samples$lib.size*1e-3, ylab="Library size (thousands)", cex.names=0.75,names=samplenames, las=3, ylim=c(0,max(vx$samples$lib.size*1e-3)+100), col=island, main="P. Vivax")
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()

# Total number of genes after filtering
pdf("nGenes_plasmodium_postFiltering.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(pf$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Falciparum", ylab="n Genes", cex.names=0.75, col=island, names=samplenames, las=3, ylim=c(0,max(apply(pf$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
barplot(apply(vx$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n P. Vivax", ylab="n Genes", cex.names=0.75, col=island, names=samplenames, las=3, ylim=c(0,max(apply(vx$counts, 2, function(c)sum(c!=0)))+1000))
legend(x="topright", col=unique(island), legend=c("Mappi", "Mentawai", "Sumba"), pch=15, cex=0.8)
dev.off()






