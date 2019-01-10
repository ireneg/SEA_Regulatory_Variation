# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood")

# Removing heteroscedascity with voom and fitting linear models -----------------------------------------------------------------------

# First, set up design matrix
# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# set up contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))

# plot mean-variance trend for different normalisation methods
pdf("Limma_voom_allNormalisationMethods.pdf", height=8, width=8)
par(mfrow=c(2,2))
y2$samples$norm.factors <- 1
v <- voom(y, design, plot=TRUE)
title(main="None", line=0.5)
for (method in c("upperquartile", "TMM", "RLE")) {
    y <- calcNormFactors(y, method=method)
    v <- voom(y, design, plot=TRUE)
    title(main=method, line=0.5)
}
dev.off()

# we'll go ahead an stick with TMM normalisation
y <- calcNormFactors(y, method="TMM")

# now go ahead with voom normalisation
pdf("Limma_voom_TMMNormalisation.pdf", height=8, width=12)
par(mfrow=c(ncol=1,nrow=2))
v <- voom(y, design, plot=TRUE)
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Mean-variance trend elimination")
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

# verify that my control genes are not significantly DE

# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# volcano plot with points of housekeeping genes and male-specific, y chromsome genes
pdf("VolcanoPlots.pdf", height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    plot(efit$coef[,i], -log10(as.matrix(efit$p.value)[,i]), pch=20, main=colnames(efit)[i], xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit$coef[,i][which(names(efit$coef[,i]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    # text(efit$coef[,i][which(names(efit$coefficients[,i]) %in% names(which(abs(dt[,i])==1)))], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coefficients[,i]) %in% names(which(abs(dt[,i])==1)))]), labels=names(which(abs(dt[,i])==1)))
    # points(efit$coef[,i][which(names(efit$coef[,i]) %in% Ycontrols)], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coef[,i]) %in% Ycontrols)]) , pch=20, col=2, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
dev.off()

# PCA visualisation after correction and association with covariates ------------------------------------------------------------

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects given by John Blischak.
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[(3)]=c("Mappi")
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$sample$CD8T, y$sample$CD4T, y$sample$NK, y$sample$Bcell, y$sample$Monoy$sample$Gran),design=design)

# rename column names of lcpm
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

# plot batch
pdf(paste0("batchCorrected_pcaresults_batch.pdf"))
name="batch"
pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# plot blood
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

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf("significantCovariates_AnovaHeatmap.pdf")
pheatmap(all.pcs[1:5,c(3:20)], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file="pca_covariates_blood.txt", col.names=T, row.names=F, quote=F, sep="\t")

# Summary and visualisation of gene trends ---------------------------------------------------------------------------

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
pdf("log2FC_IslandComparisons_pval01.pdf", height=15, width=10)
# note 'p.value' is the cutoff value for adjusted p-values
topTable <- topTable(efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=1, xlim=c(-2,2), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,-0.5,0.5,1), lty=3)
for (i in 2:ncol(efit)){
    topTable <- topTable(efit, coef=i, n=Inf, p.value=0.01)
    lines(density(topTable$logFC), col=i, xlim=c(-2,2))
}
legend(x="topright", bty="n", col=1:ncol(efit), legend=colnames(efit), lty=1, lwd=2)
dev.off()

# Get top ten DE genes through topTable (FDR 0.01) with log fold change of one and save gene information to file
for (i in 1:ncol(efit)){
    # note 'p.value' is the cutoff value for adjusted p-values
    topTable <- topTable(efit, coef=i, p.value=0.01, n=10, lfc=1, sort.by="p")
    topGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006'), mart = ensembl.mart.90,values=topTable$ENSEMBL, filters="ensembl_gene_id")
    write.table(topGenes, file=paste0("top10_genes_",colnames(efit)[i],".txt"), quote=F, row.names=F)
}

# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
mycol <- colorpanel(1000,"blue","white","red")
colors=c("steelblue","goldenrod1","red2")
cc=colors[Island]
for (i in 1:ncol(efit)){
    pdf(paste0("topGenes_Heatmap_",colnames(efit)[i],".pdf"), height=15, width=15)
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    if (nrow(topTable) > 0){
        index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
        heatmap.2(v$E[index,], scale="row",labRow=v$genes$SYMBOL[index], labCol="", col=mycol, trace="none", density.info="none", margin=c(10,10), lhei=c(2,10), dendrogram="column", ColSideColors=cc, keysize=1, cexRow=1.5, main=colnames(efit)[i])
    }
    # write out topTable genes for the whole gene list
    write.table(topTable, file=paste0("topTable_",colnames(efit)[i],".txt"))
    dev.off()
}

# show the number of DE genes between all islands
pdf("vennDiagram_allSigDEGenes_pval01_FDR1.pdf", height=15, width=15)
vennDiagram(dt[,1:3], circle.col=c("red", "blue", "green"))
dev.off()

# get DE genes in common with all islands
de.common <- which(dt[,1]!=0 & dt[,2]!=0 & dt[,3]!=0)
de.common.MPI = which(dt[,2]!=0 & dt[,3]!=0)
commonGenes.MPI <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common.MPI), filters="ensembl_gene_id")
write.table(de.common.MPI, file="allCommonGenes_MPI.txt")
common.genes=efit$genes[names(de.common.MPI),]

# top ranked genes -----------------------------------------------------------------------------------------

# since populations compared to Mappi are the only populations with DE genes over a logFC of 1, we'll apply topTreat to get our top-ranked genes
SMBvsMTW <- topTable(efit, coef=1, p.value=0.01, lfc=1, n=Inf)
SMBvsMPI <- topTable(efit, coef=2, p.value=0.01, lfc=1, n=Inf)
MTWvsMPI <- topTable(efit, coef=3, p.value=0.01, lfc=1, n=Inf)

# Let's find out what these genes are doing
SMBvsMTW.topGenes = getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006', "interpro","interpro_description"), mart = ensembl.mart.90,values=SMBvsMTW$ENSEMBL, filters="ensembl_gene_id")
write.table(SMBvsMTW.topGenes, file="SMBvsMPI_topGenes.txt")
SMBvsMPI.topGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006', "interpro","interpro_description"), mart = ensembl.mart.90,values=SMBvsMPI$ENSEMBL, filters="ensembl_gene_id")
write.table(SMBvsMPI.topGenes, file="SMBvsMPI_topGenes.txt")
MTWvsMPI.topGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006', "interpro","interpro_description"), mart = ensembl.mart.90,values=MTWvsMPI$ENSEMBL, filters="ensembl_gene_id")
write.table(MTWvsMPI.topGenes, file="MTWvsMPI_topGenes.txt")

# Let's see how the expression levels of some of our favourite/most significant genes are distributed within each island. First, assign our top genes and ensembl IDs to variables
topGenes=common.genes[,2]
topEnsembl=common.genes[,1]

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

# finally, let's make the violin plots using fancy ggpubr
pdf("TopGenes_ggboxplot_Island.pdf", height=8, width=10)
counter=0
for(ensembl in topEnsembl){
    counter=counter+1
    gene.df <- data.frame(cpm(y$counts[which(y$genes$ENSEMBL==ensembl),],log=T),Island)
    colnames(gene.df)=c("CPM", "Island")
    annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
    print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=topGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
}
dev.off()


# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# reset design matrix
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 
for (i in 1:ncol(contr.matrix)){
    camera.matrix=camera(v,idx,design,contrast=contr.matrix[,i])
    write.table(camera.matrix, file=paste0("cameraMatrix_",colnames(contr.matrix)[i],".txt"))
}

# gene set testing with goSeq
for(pop in 1:ncol(efit)){
    for(pval in c(0.05, 0.01)){
        topTable <- topTable(efit, coef=pop, n=Inf, p.value=pval, lfc=1)
        gene.vector=as.integer(rownames(y) %in% rownames(topTable))
        names(gene.vector)=rownames(y)

        # set the probability weighting fcuntion, i.e., implement a weight for each gene dependent on its length
        pwf=nullp(gene.vector,"hg19","ensGene")
        # use  the  default  method  to  calculate  the  over  and  under  expressed  GO categories among DE genes
        GO.wall=goseq(pwf,"hg19","ensGene",use_genes_without_cat=TRUE)

        # now let's interpret the results. First we need to apply a multiple hypothesis testing correction set at 5% (BH method)
        enriched.GO=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO, file=paste("enrichedGOterms",colnames(efit)[pop],pval,".txt", sep="_"), quote=F, row.names=F)
        if(nrow(enriched.GO)>0){
            zz=file(paste("topTen_enrichedGOterms",pop,pval,".txt", sep="_"), open="wt")
            sink(zz)
            sink(zz, type = "message")
            # get GO terms for top ten enriched GO terms - write output with the sink() function
            for(go in 1:length(enriched.GO$category)){
                print(GOTERM[[enriched.GO$category[go]]])
            }
        }
        sink(type = "message")
        sink()

        # KEGG pathway analysis
        en2eg=as.list(org.Hs.egENSEMBL2EG)
        # Get the mapping from Entrez 2 KEGG
        eg2kegg=as.list(org.Hs.egPATH)
        # Define a function which gets all unique KEGG IDs
        # associated with a set of Entrez IDs
        grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
        # Apply this function to every entry in the mapping from
        # ENSEMBL 2 Entrez to combine the two maps
        kegg=lapply(en2eg,grepKEGG,eg2kegg)

        # produce PWF as before
        pwf=nullp(gene.vector,"hg19","ensGene")
        KEGG=goseq(pwf,gene2cat=kegg, use_genes_without_cat=TRUE)
        enriched.GO.kegg=KEGG[p.adjust(KEGG$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO.kegg, file=paste("enrichedGOkegg",pop,pval,".txt", sep="_"))
    }
}

# Enrichment Visualisation ----------------------------------------------------------------------------------------------------


# Having played around with many visualisation options, it seems as though using Revigo is the best and most user-friendly. Using the output from the GoSeq enriched Go results (FDR <0.05), we can plug our GO IDs into Revigo (http://revigo.irb.hr/) and then output this as an R script.
# First save  output (remember that SMBvsMTW has no significantly enriched GO pathways so we'll only do this for SMBvsMPI and MTWvsMPI)  

# Revigo for SMB vs MPI (pval of 0.05 FDR)
# source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/REVIGO/123_combined/REVIGO-SMBvsMPI.r")

# Revigo for SMB vs MPI
# source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/REVIGO/123_combined/REVIGO-MTWvsMPI.r")

# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

