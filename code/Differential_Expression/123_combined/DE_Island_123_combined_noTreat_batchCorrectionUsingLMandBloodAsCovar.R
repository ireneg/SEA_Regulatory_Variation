# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: IGR 03.04.2019
### Clean up paths, remove hardcoding, test for stand-alone

# first load dependencies and set input paths(i.e., count data)
### IGR: should ditch the source command and instead load saved output from object and re-reading covariate matrices etc. Right now a lot of variables from the first script get invisibly recreated in this one if they are run separately - explicitly restating paths here: 

# Set paths:
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA_testing/dataPreprocessing" #output from previous script. 
# covariatedir <- "~/"
# blooddir <- "/data/cephfs/punim0586/kbobowik/Sumba/Output/DE_Analysis/123_combined/batchRemoval/"
#repodir <- "~/repos/SEA_Regulatory_Variation/"
#source(paste0(repodir, "code/Differential_Expression/123_combined/countData_123_combined.R"))

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA_testing/DE_testing"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Load colour schemes:
KAT TO INSERT THEM HERE # This line purposefully throws an error.


# Load all packages here:
library(RColorBrewer)
library(edgeR)
library(plyr)
# library(readr)
# library(openxlsx)
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


### BEGIN ANALYSIS

# Load log CPM matrix and y object:
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))


# Removing heteroscedascity with voom and fitting linear models -----------------------------------------------------------------------

# First, set up design matrix
# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45

# We also need to create a new variable with individual IDs:
y$samples$ind <- sapply(strsplit(as.character(y$samples$samples), "[_.]"), `[`, 1)

# design <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$NK + y$samples$Bcell + y$samples$Mono + y$samples$Gran)
# colnames(design)=gsub("Island", "", colnames(design))
# #rename columns to exclude spaces and unrecognised characters
# # colnames(design)[c(1:13)]=c("Mentawai", "Sumba", "Mappi", "Age", "batch2", "batch3", "RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# FOR TESTING ONLY:
design <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$Gran)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[c(1:8)]=c("Mentawai", "Sumba", "Mappi", "Age", "batch2", "batch3", "RIN", "Gran")

# set up contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))

# plot mean-variance trend for different normalisation methods. This is helpfule for understandimg what should be done with estimating the dispersion: https://support.bioconductor.org/p/77664/
# pdf(paste0(outputdir, "EstimatingDispersion_allNormalisationMethods.pdf"), height=15, width=15)
# par(mfrow=c(4,3))

# Fails below for statistical reasons, but runs to here as a stand-alone using the load commands!

# first perform voom on unnormalised data
# y$samples$norm.factors <- 1
# estimateDisp <- estimateDisp(y, design, robust=TRUE)
# plotBCV(estimateDisp)
# title(paste0("None","\nDispersion Range = ",round(min(estimateDisp$prior.df), 2), "-", round(max(estimateDisp$prior.df), 2)))
# v <- voom(y, design, plot=TRUE)
# title(main="None", line=0.5)
# # fit linear models
# vfit <- lmFit(v, design)
# vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
# efit <- eBayes(vfit, robust=T)
# plotSA(efit, main="Mean-variance trend elimination \n None")

# # now perform normalisation on all the others
# for (method in c("upperquartile", "TMM", "RLE")) {
#     y <- calcNormFactors(y, method=method)
#     estimateDisp <- estimateDisp(y, design, robust=TRUE)
#     plotBCV(estimateDisp)
#     title(paste0(method,"\nDispersion Range = ",round(min(estimateDisp$prior.df), 2), "-", round(max(estimateDisp$prior.df), 2)))
#     v <- voom(y, design, plot=TRUE)
#     title(main=method, line=0.5)
#     # fit linear models
#     vfit <- lmFit(v, design)
#     vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#     efit <- eBayes(vfit, robust=T)
#     plotSA(efit, main=paste("Mean-variance trend elimination",method,sep="\n"))
# }
# dev.off()

# we'll go ahead an stick with TMM normalisation
y <- calcNormFactors(y, method="TMM")

# now go ahead with voom normalisation
pdf("Limma_voom_TMMNormalisation.pdf", height=8, width=12)
par(mfrow=c(ncol=1,nrow=2))
v <- voom(y, design, plot=TRUE)
dupcor <- duplicateCorrelation(v, design, block=y$samples$ind) # 46 non-convergences
dupcor$consensus # sanity check
# [1] 0.7236156
median(v$weights) # another sanity check:
# [1] 19.03065

vDup <- voom(y, design, plot=TRUE, block=y$samples$ind, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(vDup, design, block=y$samples$ind) # 28 non-convergences
dupcor$consensus # sanity check pt 2
# [1] 0.7239365
median(vDup$weights) # another sanity check, pt 2 - small change, so it didn't matter too much, but it is good to have done it. 
# [1] 18.71378

# fit linear models
# With duplicate correction and blocking:
    voomDupVfit <- lmFit(vDup, design, block=y$samples$ind, correlation=dupcor$consensus)
    voomDupVfit <- contrasts.fit(voomDupVfit, contrasts=contr.matrix)
    voomDupEfit <- eBayes(voomDupVfit, robust=T)

# And without:
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit, robust=T)

plotSA(voomDupEfit, main="Mean-variance trend elimination with duplicate correction")
plotSA(efit, main="Mean-variance trend elimination without duplicate correction")
dev.off()

### IGR: Quick DE qc after fitting lms:
# With dup correction:
    voomDupTopTableSMB.MTW <- topTable(voomDupEfit, coef=1, p.value=0.01, n=Inf, sort.by="p")
    voomDupTopTableSMB.MPI <- topTable(voomDupEfit, coef=2, p.value=0.01, n=Inf, sort.by="p")
    voomDupTopTableMTW.MPI <- topTable(voomDupEfit, coef=3, p.value=0.01, n=Inf, sort.by="p")

    summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down       1393     2158     1590
    # NotSig    10575     8760     9867
    # Up         1007     2057     1518

    summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down        152      574      432
    # NotSig    12631    11695    12041
    # Up          192      706      502

    summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down         14       81       96
    # NotSig    12942    12706    12744
    # Up           19      188      135

# And without:
    topTableSMB.MTW <- topTable(efit, coef=1, p.value=0.01, n=Inf, sort.by="p")
    topTableSMB.MPI <- topTable(efit, coef=2, p.value=0.01, n=Inf, sort.by="p")
    topTableMTW.MPI <- topTable(efit, coef=3, p.value=0.01, n=Inf, sort.by="p")

    summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down       1447     2320     1728
    # NotSig    10439     8483     9641
    # Up         1089     2172     1606

    summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down        128      603      441
    # NotSig    12610    11624    12026
    # Up          237      748      508

    summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
    #        SMBvsMTW SMBvsMPI MTWvsMPI
    # Down         11       82       97
    # NotSig    12935    12680    12741
    # Up           29      213      137


# And finally, correlation between those two measurements - sort by gene first, then cor test
    MTW.MPI <- join(voomDupTopTableMTW.MPI, topTableMTW.MPI, by="genes")
    cor(MTW.MPI[,6], MTW.MPI[,12], method="spearman", use="complete")
    # [1] 0.969342

    SMB.MPI <- join(voomDupTopTableSMB.MPI, topTableSMB.MPI, by="genes")
    cor(SMB.MPI[,6], SMB.MPI[,12], method="spearman", use="complete")
    # [1] 0.929687

    SMB.MTW <- join(voomDupTopTableSMB.MTW, topTableSMB.MTW, by="genes")
    cor(SMB.MTW[,6], SMB.MTW[,12], method="spearman", use="complete")
    # [1] 0.8754186

# Didn't break anything! But comparisons Sumba look notably different, which is interesting.

# Back to Kat's regularly scheduled code:

# We can view how UQ normalisation performed using MD plots. This visualizes the library size-adjusted log-fold change between
# two libraries (the difference) against the average log-expression across those libraries (themean). MD plots are generated by comparing sample 1 against an artificial
# library constructed from the average of all other samples. Ideally, the bulk of genes should be centred at a log-fold change of zero.  This indicates
# that any composition bias between libraries has been successfully removed

pdf("MDPlots_TMM_Normalisation_OutliersCheck.pdf", height=15, width=15)
par(mfrow=c(4,4))
for (i in 1:ncol(y)){
  plotMD(cpm(y, log=TRUE), column=i)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

# QC after fitting linear models --------------------------------------------------------------------------------------

# check to see p-value distribution is normal
pdf("PvalueDist_NotAdjusted.pdf", height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    hist(efit$p.value[,i], main=colnames(efit)[i], ylim=c(0,max(table(round(efit$p.value[,i], 1)))+1000), xlab="p-value")
}
dev.off()

# check p-value distribution for adjusted p-values
pdf("PvalueDist_Adjusted.pdf", height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, n=Inf)
    histData <- hist(topTable$adj.P.Val, main=colnames(efit)[i], xlab="p-value")
    hist(topTable$adj.P.Val, main=colnames(efit)[i], ylim=c(0,max(histData$counts)+1000), xlab="p-value")
}
dev.off()

# Verify that control housekeeping genes are not significantly DE. Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# Volcano plot with points of housekeeping genes
pdf("VolcanoPlots.pdf", height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    plot(efit$coef[,i], -log10(as.matrix(efit$p.value)[,i]), pch=20, main=colnames(efit)[i], xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit$coef[,i][which(names(efit$coef[,i]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
dev.off()

# PCA visualisation after correction and association with covariates ------------------------------------------------------------

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects was given by the lovely John Blischak.
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[(3)]=c("Mappi")
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$sample$CD8T, y$sample$CD4T, y$sample$NK, y$sample$Bcell, y$sample$Monoy$sample$Gran),design=design)

# rename column names of lcpm to sample names (so that they are shorter and easier to read)
colnames(lcpm)=samplenames

# get which samples are replicates
allreplicated=as.factor(samplenames %in% allreps)

# PCA plotting function
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

# Plot PCA
for (name in covariate.names[c(1:10,17:18)]){
  if (nlevels(get(name)) < 26){
    pdf(paste0("batchCorrected_pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# plot batch (this has to be done separately since it ha sa different colour scheme)
pdf(paste0("batchCorrected_pcaresults_batch.pdf"))
name="batch"
pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# plot blood - this one also has a differnet gradient colour scheme
for (name in covariate.names[c(11:16)]){
    initial <- cut(get(name), breaks = seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0("batchCorrected_pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
}

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any remaining confounding/batch
pdf("significantCovariates_AnovaHeatmap.pdf")
pheatmap(all.pcs[1:5,c(3:20)], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file="pca_covariates_blood.txt", col.names=T, row.names=F, quote=F, sep="\t")

# Summary and visualisation of gene trends ---------------------------------------------------------------------------

# first see which logFC threshold is best
logFC.df=matrix(nrow=3,ncol=3)
colnames(logFC.df)=colnames(efit)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    values=c(sum(abs(dt[,1])), sum(abs(dt[,2])), sum(abs(dt[,3])))
    logFC.df[counter,]=values
}
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file="logFC_thresholds.txt")

dt <- decideTests(efit, p.value=0.01, lfc=1)
# get summary of decide tests statistics
write.table(summary(dt), file="numberSigDEgenes_voom_efit.txt")
# write out top table results 
write.fit(efit, dt, file="toptable_SigGenes_voom_efit.txt")

# graphical representation of DE results through MD plot
pdf("MD_Plots_pval01_lfc1.pdf", height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    o <- which(names(efit$Amean) %in% names(which(abs(dt[,i])==1)))
    x <- efit$Amean
    z <- efit$coefficients[,i]
    t=which(names(efit$coefficients[,i]) %in% names(which(abs(dt[,i])==1)))
    G <- efit$genes[names(which(abs(dt[,i])==1)),]$SYMBOL
    plotMD(efit, column=i, status=dt[,i], main=colnames(efit)[i], hl.col=c("blue","red"), values=c(-1,1))
    abline(h=c(1,-1), lty=2)
    legend(legend=paste(names(summary(dt)[,i]), summary(dt)[,i], sep="="), x="bottomright", border=F, bty="n")
    text(x[o], z[t], labels=G)
}
dev.off()

# plot log2 fold change between islands
pdf("log2FC_IslandComparisons_pval01.pdf")
# note 'p.value' is the cutoff value for adjusted p-values
topTable <- topTable(efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=9, xlim=c(-2,2), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,-0.5,0.5,1), lty=3)
counter=0
for (i in 2:ncol(efit)){
	counter=counter+1
    topTable <- topTable(efit, coef=i, n=Inf, p.value=0.01)
    lines(density(topTable$logFC), col=9+counter, xlim=c(-2,2))
}
legend(x="topright", bty="n", col=9:11, legend=colnames(efit), lty=1, lwd=2)
dev.off()

# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
# first, make a heatmap of all top genes in one pdf

# reset ensemble row names to gene symbols
rownames(v$E)=v$genes$SYMBOL

# set up annotation
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

df1=data.frame(island = as.character(Island))
df2=data.frame(batch = as.numeric(batch))
ha1 = HeatmapAnnotation(df = df1, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3)))

pdf("HeatmapAllPops.pdf", height=15, width=15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
# set up layout row position
layout.row=c(1,1,2)
# set up column position
layout.col=c(1,2,1)

for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf, sort.by="p")
    index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
    pushViewport(viewport(layout.pos.row = layout.row[i], layout.pos.col = layout.col[i]))
    draw(Heatmap(t(scale(t(v$E[index,]))), col=col_fun, column_title = colnames(efit)[i], top_annotation = ha1, show_row_names = T, show_heatmap_legend = F, show_column_names = F, name = "Z-Score"),show_annotation_legend = FALSE,newpage=F)
    upViewport()

}

lgd = Legend(at = c(-4,0,4), title = "Row Z-Score", col_fun = col_fun, grid_height = unit(1, "cm"), grid_width = unit(10, "mm"))
lgd2 = Legend(at = c(1,2,3), legend_gp = gpar(fill = 1:3), labels=c("Mentawai", "Sumba","West Papua"),title = "Island", grid_height = unit(1, "cm"), grid_width = unit(10, "mm"))
pd = packLegend(lgd, lgd2, direction = "horizontal")

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(pd)
upViewport()
dev.off()

# We can also make individual pdfs of the top genes
island1=c("Sumba","Mentawai","West Papua")
island2=c("Sumba","Mentawai","West Papua")

counter=0

for (i1 in island1){
    island2=island2[-1]
    for (i2 in island2){
        counter=counter+1
        topTable <- topTable(efit, coef=counter, p.value=0.01, lfc=1, n=Inf, sort.by="p")
        index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
        df=data.frame(island = as.character(Island[grep(paste(i1,i2,sep="|"), Island)]))
        ha =  HeatmapAnnotation(df = df, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3)))
        pdf(paste("HeatmapTopeGenes",i1,i2,".pdf",sep="_"), height=10, width=15)
        draw(Heatmap(t(scale(t(v$E[index,])))[,grep(paste(i1,i2,sep="|"), Island)], col=col_fun, column_title = colnames(efit)[counter], top_annotation = ha, show_row_names = T, show_heatmap_legend = T, show_column_names = F, name = "Z-Score"),show_annotation_legend = TRUE,newpage=F)
        dev.off()
    }
}

# Get top DE genes through topTable (FDR 0.01) with log fold change of one and save gene information to file
for (i in 1:ncol(efit)){
    # note 'p.value' is the cutoff value for adjusted p-values
    topTable <- topTable(efit, coef=i, p.value=0.01, n=Inf, lfc=1, sort.by="p")
    topGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006', "interpro","interpro_description"), mart = ensembl.mart.90,values=topTable$ENSEMBL, filters="ensembl_gene_id")
    write.table(topGenes, file=paste0("top_genes_",colnames(efit)[i],".txt"), quote=F, row.names=F)
    write.table(topTable, file=paste0("top_Table_",colnames(efit)[i],".txt"), quote=F, row.names=F)
}

# show the number of DE genes between all islands
pdf("vennDiagram_allSigDEGenes_pval01_FDR1.pdf", height=15, width=15)
vennDiagram(dt[,1:3], circle.col=c(9,10,11))
dev.off()

# get DE genes in common with populations compared to Mappi, i.e., SMBvsMPI and MTWvsMPI (since we think this is an interesting island comparison)
de.common.MPI = which(dt[,2]!=0 & dt[,3]!=0)
# get what these genes are doing and save them to a file
commonGenes.MPI <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common.MPI), filters="ensembl_gene_id")
write.table(de.common.MPI, file="allCommonGenes_MPI.txt")
# save the common gene names 
de.common.MPI=efit$genes[names(de.common.MPI),]

# now plot the common genes to see if they're being regulated in the same direction
tt.SMBvsMPI=topTable(efit, coef=2, p.value=0.01, n=Inf, lfc=1, sort.by="p")
tt.MTWvsMPI=topTable(efit, coef=3, p.value=0.01, n=Inf, lfc=1, sort.by="p")
pdf("logFC_commonMPIgenes.pdf")
plot(tt.SMBvsMPI[rownames(de.common.MPI),"logFC"], tt.MTWvsMPI[rownames(de.common.MPI),"logFC"], xlab="logFC SMBvsMPI", ylab="logFC MTWvsMPI", pch=20, main="Common DE Genes", xlim=c(-5,5), ylim=c(-6,6))
text(tt.SMBvsMPI[rownames(de.common.MPI),"logFC"], tt.MTWvsMPI[rownames(de.common.MPI),"logFC"], labels=tt.SMBvsMPI[rownames(de.common.MPI),"SYMBOL"], pos=3)
abline(h=0,v=0, lty=2)
dev.off()

# top ranked genes -----------------------------------------------------------------------------------------

# Let's see how the expression levels of all of the significantly DE genes in population comparisons with Mappi are distributed within each island. First, assign our top genes and ensembl IDs to variables
topGenes=de.common.MPI[,2]
topEnsembl=de.common.MPI[,1]

# To visualise distributions, we'll be making violin plots using ggpubr which needs p-value labels. Let's go ahead and make a matrix to input this into ggpubr
# first set up matrix
topGenes.pvalue=matrix(nrow=length(topEnsembl), ncol=ncol(efit))
rownames(topGenes.pvalue)=topEnsembl
colnames(topGenes.pvalue)=colnames(efit)
for (i in 1:ncol(efit)){
    # get significant genes over a logFC of 1 for all Island comparisons
    topTable <- topTable(efit, coef=i, n=Inf)
    for(j in topEnsembl){
        # input the adjusted p.value for each gene
        topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]
    }
}

# make pvalues into scientific notation with max 3 digits
topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# convert e notation to base 10 notation
topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

# We can make the violin plots using ggpubr
pdf("TopGenes_ggboxplot_Island.pdf", height=8, width=10)
counter=0
for(ensembl in topEnsembl){
    counter=counter+1
    gene.df <- data.frame(v$E[which(v$genes$ENSEMBL==ensembl),],Island)
    colnames(gene.df)=c("CPM", "Island")
    annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
    print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=topGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
}
dev.off()

# after analysing the distributions and reading up on some of the genes, my three favourite genes are Siglec6, Siglec7, and MARCO. Lets plot out the distribution solely for these three genes
favGenes=c("SIGLEC6","SIGLEC7","MARCO","RSAD2","AIM2","TNFSF4")
favEnsembl=de.common.MPI[,1][sapply(1:length(favGenes), function(x)grep(favGenes[x],de.common.MPI[,2]))]

# set up pvalue matrix
topGenes.pvalue=matrix(nrow=length(favEnsembl), ncol=ncol(efit))
rownames(topGenes.pvalue)=favEnsembl
colnames(topGenes.pvalue)=colnames(efit)
for (i in 1:ncol(efit)){
    # get significant genes over a logFC of 1 for all Island comparisons
    topTable <- topTable(efit, coef=i, n=Inf)
    for(j in favEnsembl){
        # input the adjusted p.value for each gene
        topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]
    }
}

# make pvalues into scientific notation with max 3 digits
topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# convert e notation to base 10 notation
topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

# We can make the violin plots using ggpubr
counter=0
for(ensembl in favEnsembl){
    counter=counter+1
    # pdf(paste0("FavouriteGenes_ggboxplot_",favGenes[counter],".pdf"), height=8, width=10)
    gene.df <- data.frame(v$E[which(v$genes$ENSEMBL==ensembl),],Island)
    colnames(gene.df)=c("CPM", "Island")
    annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
    # print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=favGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
    assign(favGenes[counter], ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=favGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 3, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
}

pdf("favouriteTopGenes_distribution_Island.pdf", height=12, width=15)
ggarrange(SIGLEC6,SIGLEC7,MARCO,AIM2,TNFSF4,RSAD2)
dev.off()


