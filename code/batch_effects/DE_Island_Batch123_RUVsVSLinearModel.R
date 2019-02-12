# script created by KSB, 06.06.18
# Perform DE analysing using RUVs vs a traditional linear modelling approach

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVvsLinearModel")
library(NineteenEightyR)
require(VennDiagram)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(eulerr)

## setup ---------------------------------------

colors <- electronic_night(n=5)
# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# First, construct a matrix specifying the replicates. 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}
genes <- rownames(y)

## RUVs ----------------------------------------

# First, set up pheno data with all known factors of unwanted variation
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))
# normalise with upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")
set1 <- RUVs(set, genes, k=5, replicates)

# create design matrix
design.RUV <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design.RUV)=gsub("Island", "", colnames(design.RUV))
colnames(design.RUV)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")

# set up gene names
geneid <- rownames(z)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
z$genes <- genes

# make contrast matrix
contr.matrix.RUVs <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design.RUV))
v.RUV <- voom(z, design.RUV, plot=FALSE)
# fit linear models
vfit.RUV <- lmFit(v.RUV, design.RUV)
vfit.RUV <- contrasts.fit(vfit.RUV, contrasts=contr.matrix.RUVs)
efit.RUV <- eBayes(vfit.RUV)
# for now, we'll continue sticking with a k of 5. Let's see how many DE genes we get setting a lfc of 1 and 0.01 pvalue
dt <- decideTests(efit.RUV, p.value=0.01, lfc=1)
summary(dt)

#       SMBvsMTW SMBvsMPI MTWvsMPI
# Down          6       33       43
# NotSig    11450    11349    11332
# Up           14       88       95

# which genes are these?
mycol <- colorpanel(1000,"blue","white","red")
colors=c("steelblue","goldenrod1","red2")
cc=colors[Island]
pdf("topGenes_Heatmap_RUVs.pdf", height=15, width=15)
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    if (nrow(topTable) > 0){
        index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:100])
        heatmap.2(v$E[index,], scale="row",labRow=v$genes$SYMBOL[index], labCol=Island, col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column", ColSideColors=cc, keysize=1, cexRow=1.5, main=colnames(efit)[i])
    }
    # write out topTable genes for the whole gene list
    write.table(topTable, file=paste0("topTable_",colnames(efit)[i],".txt"))
}
dev.off()


# now let's do the same on a regular limma LM pipeline ----------------------------------------------------------------------------------------------------

# set up design
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design.lm <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design.lm)=gsub("Island", "", colnames(design.lm))
#rename columns to exclude spaces and unrecognised characters
colnames(design.lm)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
# set up contrast matrix using nested design
contr.matrix.limma <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design.lm))
v.lm <- voom(y, design.lm, plot=FALSE)
# fit linear models
vfit.lm <- lmFit(v.lm, design.lm)
vfit.lm <- contrasts.fit(vfit.lm, contrasts=contr.matrix.limma)
efit.lm <- eBayes(vfit.lm)
dt.lm <- decideTests(efit.lm, p.value=0.01, lfc=1)
summary(dt.lm)

#       SMBvsMTW SMBvsMPI MTWvsMPI
# Down          1       60       63
# NotSig    11455    11273    11309
# Up           14      137       98


# LM vs RUV similiarities --------------------------------

# Get which genes overlap between RUVs and linear model
pdf("VennDiagram_LMvsRUVsComparison.pdf", height=4, width=12)
# sapply(1:3, function(x) vennDiagram(cbind(dt[,x],dt.lm[,x]), circle.col=c("red","blue"), names=c("RUVs", "LinearModel"), main=colnames(dt.lm)[x]))
par(mfrow=  c(1,3),mar = c(0,0,0,0)) 
for (x in c(1:3)){
    vennDiagram(cbind(dt[,x],dt.lm[,x]), circle.col=c("red","blue"), names=c("RUVs", "LinearModel"))
    title(colnames(dt.lm)[x], line=-6)
}
dev.off()

for (i in 1:3){
	commonGenes=dt[,i] & dt.lm[,i]
	commonGenes=which(commonGenes == TRUE)
	commonGenes.results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=names(commonGenes), filters="ensembl_gene_id")
	write.table(commonGenes.results,file=paste0("allCommonGenes_",colnames(dt)[i],".txt"))
}

# test top genes -----------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 

# analysis of top gene pathways through Camera
pdf("Top10CameraGeneSets_VennDiagram.pdf", height=12, width=14)
for (i in 1:3){
    camera.lm.matrix=camera(v.lm,idx,design.lm,contrast=contr.matrix.limma[,i])
    camera.RUV.matrix=camera(v.RUV,idx,design.RUV,contrast=contr.matrix.RUVs[,i])
    LM=rownames(camera.lm.matrix[1:100,])
    RUV=rownames(camera.RUV.matrix[1:100,])
    intersect=length(intersect(LM, RUV))
    fit2 <- euler(c(LM=length(LM)-intersect,RUV=length(RUV)-intersect,"LM&RUV"=intersect))
    assign(colnames(efit.lm)[i],plot(fit2, fills = c("dodgerblue4","darkgoldenrod1"),edges = FALSE,fontsize = 8,quantities = list(fontsize = 10, col="white"), alpha=0.8, main=colnames(efit.lm)[i], cex=1, counts=T))
    write.table(camera.lm.matrix, file=paste0("CAMERA_LM_",colnames(efit)[i],".txt"))
    write.table(camera.lm.matrix, file=paste0("CAMERA_RUV_",colnames(efit)[i],".txt"))
}
grid.arrange(SMBvsMTW, SMBvsMPI, MTWvsMPI, ncol=3, widths=c(1.3,1.1,1))
# you can also plot this with cowplot: plot_grid(SMBvsMTW, SMBvsMPI, MTWvsMPI,rel_heights = c(1/2, 1/4, 1/4), align="h", ncol=3). Source: https://stackoverflow.com/questions/36198451/specify-widths-and-heights-of-plots-with-grid-arrange
dev.off()

# look at the top genes using top table
for (i in 1:3){
    topTable.lm <- topTable(efit.RUV, coef=i, p.value=0.01, n=Inf, sort.by="P")
    topTable.RUV <- topTable(efit.lm, coef=i, p.value=0.01, n=Inf, sort.by="P")
    LM=topTable.lm$SYMBOL[1:100]
    RUV=topTable.RUV$SYMBOL[1:100]
    intersect=length(intersect(LM, RUV))
    fit2 <- euler(c(LM=length(LM)-intersect,RUV=length(RUV)-intersect,"LM&RUV"=intersect))
    assign(colnames(efit.lm)[i],plot(fit2, fills = c("dodgerblue4","darkgoldenrod1"),edges = FALSE,fontsize = 8,quantities = list(fontsize = 10, col="white"), alpha=0.8, main=colnames(efit.lm)[i], cex=1, counts=T))
}
grid.arrange(SMBvsMTW, SMBvsMPI, MTWvsMPI, ncol=3, widths=c(1.3,1.1,1))

# plot pca and associations of both -------------------------------------------

# rename column names of lcpm to sample names (so that they are shorter and easier to read)
colnames(lcpm)=samplenames

# first set up batch-corrected lm data
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[(3)]=c("Mappi")
# here, we're using log2-cpm values of the TMM-normalised data. This seems to be suffiecient (i.e., voom-corrected output is not necessary). See link here: https://support.bioconductor.org/p/76837/.
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$sample$CD8T, y$sample$CD4T, y$sample$NK, y$sample$Bcell, y$sample$Monoy$sample$Gran),design=design)

# This is our PCA plotting function. We might need it later to explore differnet dimensions of the PCA
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        # points(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], col="black", pch=8, cex=2)
        # text(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=3)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
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
all.covars.df <- y$samples[,c(3,5:22)] 

# Let's just look at the forst pca for now and make a function to do this (without any labels to make the visualisation of batch clear)
one.dimension=function(data) {
    pca <- prcomp(t(data), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    plot(pca$x[,1], pca$x[,2], col=batch.col[as.numeric(batch)], pch=as.numeric(y$samples$batch) + 14, cex=2, xlab=paste0("PC", 1, " (", round(pca.var[1]*100, digits=2), "% of variance)"), ylab=paste0("PC", 2, " (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""))
    legend(legend=unique(batch), pch=unique(as.numeric(y$samples$batch)) + 14, x="bottomright", col=unique(batch.col[as.numeric(batch)]), cex=0.6, title="batch", border=F, bty="n")
}

# Compare the first PCA of the LM-corrected method to the RUVs corrected method
pdf("PCA_RUVvsLM_FirstDimension.pdf", height=8, width=15)
par(mfrow=c(1,2))
one.dimension(batch.corrected.lcpm)
title("Limma-LM")
# here, we're looking at the RUVs-corrected, log2-transformed output
one.dimension(cpm(normCounts(set1), log=T))
title("RUVs")
dev.off()

# We also want to look at a heatmap of the significant covariates compared to each dimension of the PCA

# first for limma-corrected data
pcaresults.limma <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
all.pcs.limma <- pc.assoc(pcaresults.limma)
all.pcs.limma$Variance <- pcaresults.limma$sdev^2/sum(pcaresults.limma$sdev^2)
# we can save the table to look at the association in all of the dimensions
write.table(all.pcs.limma, file="PC_associations_limmaCorrectedData.txt")

# now for the RUVs-corrected data
pcaresults.RUVs <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
all.pcs.RUVs <- pc.assoc(pcaresults.RUVs)
all.pcs.RUVs$Variance <- pcaresults.RUVs$sdev^2/sum(pcaresults.RUVs$sdev^2)
# save the associations for the RUVs data
write.table(all.pcs.RUVs, file="PC_associations_RUVsCorrectedData.txt")

# set up heatmap information
# first limma
limma.pc=all.pcs.limma[1:5,c(3:20)]
limma.hm=melt(limma.pc)
limma.hm$Dimension=rep(1:5,ncol(limma.pc))
colnames(limma.hm)[c(1,2)]=c("Covariate","Association")

# RUVs
RUV.pc=all.pcs.RUVs[1:5,c(3:20)]
RUV.hm=melt(RUV.pc)
RUV.hm$Dimension=rep(1:5,ncol(RUV.pc))
colnames(RUV.hm)[c(1,2)]=c("Covariate","Association")

# set up colour palette
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
# Save output of both plots
limma=ggplot(data = limma.hm, aes(x = Covariate, y = Dimension)) + geom_tile(aes(fill = Association))  + scale_fill_gradientn(colours = hm.palette(100)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Limma") 
RUV=ggplot(data = RUV.hm, aes(x = Covariate, y = Dimension)) + geom_tile(aes(fill = Association))  + scale_fill_gradientn(colours = hm.palette(100)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("RUVs") 
# use ggarange to plot both heatmpas
pdf("Heatmap_SigCovar_RUVvsLM.pdf", height=7, width=15)
ggarrange(limma,RUV)
dev.off()

# Compare p-value distributions of both -----------------------------------------------

pdf("VolcanoPlot_RUVvsLM.pdf", height=10, width=15)
layout(matrix(c(1,2,3,4,5,6), nrow=3,byrow = F))
for (i in 1:ncol(efit.lm)){
    plot(efit.lm$coef[,i], -log10(as.matrix(efit.lm$p.value)[,i]), pch=20, main=paste("LM",colnames(efit)[i],sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit.lm$coef[,i][which(names(efit.lm$coef[,i]) %in% hkControls)], -log10(as.matrix(efit.lm$p.value)[,i][which(names(efit.lm$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
for (i in 1:ncol(efit.RUV)){
    plot(efit.RUV$coef[,i], -log10(as.matrix(efit.RUV$p.value)[,i]), pch=20, main=paste("RUV",colnames(efit)[i],sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit.RUV$coef[,i][which(names(efit.RUV$coef[,i]) %in% hkControls)], -log10(as.matrix(efit.RUV$p.value)[,i][which(names(efit.RUV$coef[,i]) %in% hkControls)]) , pch=20, col=2, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=2)
    abline(v=c(-1,1))
}
dev.off()

# See which genes are not in common and why ----------------------------------------------------------------------

# reset wd to keep things organised
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LMvsRUVs_allCommonGenes")

# set up not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# see which genes aren't common and why this might be
pdf("genesNotInCommon_RUVsVsLM.pdf", height=10, width=15)
for (i in 1:ncol(efit)){
    RUVs=dt[,i][which(abs(dt[,i])==T)]
    lm=dt.lm[,i][which(abs(dt.lm[,i])==T)]
    notInLM=names(RUVs[which(names(RUVs) %!in% names(lm))])
    notInRUV=names(lm[which(names(lm) %!in% names(RUVs))])
    topTable.lm <- topTable(efit.lm, coef=i, n=Inf)
    topTable.RUV <- topTable(efit, coef=i, n=Inf)

    plot(topTable.lm[notInLM,]$adj.P.Val, abs(topTable.lm[notInLM,]$logFC), pch=16, main=paste("Top DE Genes Not in Common", colnames(efit)[i], sep="\n"), xlab="adjusted p.value", ylab="logFC")
    points(topTable.RUV[notInRUV,]$adj.P.Val, abs(topTable.RUV[notInRUV,]$logFC), col="red", pch=16)
    text(topTable.lm[notInLM,]$adj.P.Val, abs(topTable.lm[notInLM,]$logFC), labels=topTable.lm[notInLM,]$SYMBOL, pos=1)
    text(topTable.RUV[notInRUV,]$adj.P.Val, abs(topTable.RUV[notInRUV,]$logFC), labels=topTable.RUV[notInRUV,]$SYMBOL, col="red", pos=2)
    legend("topright", c("genes not in lm", "genes not in RUVs"), pch=16, col=c("black", "red"))
    abline(v=0.01, lty=2)
    abline(h=1, lty=2)
}
dev.off()









