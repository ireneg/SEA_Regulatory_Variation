# script created by KSB, 06.06.18
# data exploration (MDS, PCA, and sample outlier analysis) of Indonesia RNA-seq data (both batches)

# load dependencies: human count data and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataExploration")

# MDS
pdf("indoRNA_MDS_123Combined_allCovariates.pdf", height=30, width=25)
par(mfrow=c(5,4))
for (name in covariate.names[c(1:10,17:18)]) {
    plotMDS(lcpm, labels=get(name), col=as.numeric(get(name)))
    title(main=name)
}
# plot blood separately
plotMDS(lcpm, labels=batch, col=batch.col[as.numeric(batch)])
title(main="batch")
for (name in covariate.names[c(11:16)]) {
    # assign gradient colours for blood cell types
    initial <- cut(get(name), breaks = seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    plotMDS(lcpm, labels=get(name), col=bloodCol)
    title(main=name)
}
dev.off()

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
        points(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], col="black", pch=8, cex=2)
        text(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=3)
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
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# plot batch
pdf(paste0("pcaresults_batch.pdf"))
name="batch"
pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# plot blood
for (name in covariate.names[c(11:16)]){
    initial <- cut(get(name), breaks = seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0("pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
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

# Getting the loadings and top genes:
geneLoadings <- as.data.frame(pcaresults$rotation)
# get loadings for PCAs 1:10
for (i in 1:10){
    absLoadings=geneLoadings[order(-abs(geneLoadings[i])),][1:100,]
    assign(paste0("topGenes_",i), rownames(absLoadings))
    write.table(absLoadings, file=paste0("Top_100_absolute_loadings_PCA",i,".txt"), quote=F, row.names=T, col.names=T)
}

# Get relationship of all covariates
new_cov=apply(covariates[2:11], 2, FUN=function(x){x=as.numeric(as.factor(x))})
pdf("covariateHeatmap.pdf", height=10, width=15)
pheatmap(t(new_cov),cluster_col=FALSE,cluster_rows=FALSE, labels_col=covariates[,1], colorRampPalette(rev(brewer.pal(n=10,name="Spectral")))(100),cex=1.1, main="Covariates")
dev.off()

# set colnames to samplenames to conserve space
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)
samplenames[104:123]=sapply(samplenames[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))
colnames(lcpm)=make.unique(samplenames)

# Dissimilarity matrix with euclidean distances
pdf("SampleDistances.pdf", height=10, width=15)
par(mar=c(6.1,4.1,4.1,2.1), mfrow=c(2,1))
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

## Heatmap annotation

# old annotation
# df=data.frame(island = as.character(Island), batch = as.numeric(batch))
# ha = HeatmapAnnotation(df = df, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3), batch = c("1" = 5, "2" = 6, "3" = 7)))
# Heatmap(cor(lcpm,method="spearman"), col=viridis(100), column_title = "Spearman Correlation \n log2-CPM", name="Corr Coefficient", top_annotation = ha)

# set up annotation
df1=data.frame(island = as.character(Island))
df2=data.frame(batch = as.numeric(batch))
ha1 = HeatmapAnnotation(df = df1, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3)))
ha2 = rowAnnotation(df = df2, col= list(batch=c("1" = batch.col[1], "2" = batch.col[2], "3" = batch.col[3])))

# lcpm distances
pdf("lcpmCorrelationHeatmaps.pdf", height=10, width=15)
Heatmap(cor(lcpm,method="pearson"), col=magma(100), column_title = "Pearson Correlation \n log2-CPM", name="Corr Coeff", top_annotation = ha1, show_row_names = FALSE, column_names_gp=gpar(fontsize = 8)) + ha2
dev.off()

# analyse sample variation
island=c("MTW", "SMB","MPI")
village=c("MPI", "MDB", "TLL", "ANK", "WNG")

#for (method in c("spearman", "pearson")){
#    correlation=cor(lcpm,method=method)
#    # samples from same village
#    withinVillage=list()
#    for (sample in village){
#        corLcpm=correlation[grepl(sample,rownames(correlation)),grepl(sample,colnames(correlation))]
#        diag(corLcpm)=NA
#        allCPM=melt(corLcpm)
#        allCPM$value[duplicated(allCPM$value)]=NA
#        withinVillage[[sample]]=allCPM[,3]
#    }
#    pdf(paste0("IslandVariation_",method,".pdf"),height=15, width=15)
#    par(mai=rep(0.5, 4))
#    layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE))
#    withinIsland=list()
#    variation=list()
#    for (sample in island){
#        corLcpm=correlation[grepl(sample,rownames(correlation)),grepl(sample,colnames(correlation))]
#        diag(corLcpm)=NA
#        allCPM=melt(corLcpm)
#        allCPM$value[duplicated(allCPM$value)]=NA
#        withinIsland[[sample]]=allCPM[,3]
#        sampVar=allCPM[which(abs(scale(allCPM[,3])) > 3),]
#        if (nrow(sampVar) > 0) {
#            variation[[sample]]=transform(sampVar, names=paste(sapply(strsplit(as.character(sampVar[,1]),"[-.]"), `[`, 3), sapply(strsplit(as.character(sampVar[,2]),"[-.]"), `[`, 3), sep="vs"))
#        }
#    }
#    boxplot(withinIsland, main=paste("Within-Island Variation",method,sep="\n"))
#    stripchart(withinIsland, vertical=T, method = "jitter", add = TRUE, pch = 20, cex=2, col=c(1,2,3))
#    sapply(1:2, function(x) text(x=x, y=variation[[x]][,3], labels=as.character(variation[[x]][1:nrow(variation[[x]]),4]), cex=0.8, pos=1))
#        
#    # Inter island
#    interIsland <- list()
#    variation=list()
#    island2 <- island
#    for (i1 in island) {
#        island2 <- island2[-1]
#        for (i2 in island2) {
#            corLcpm=correlation[grepl(i1,rownames(correlation)),grepl(i2,colnames(correlation))]
#            diag(corLcpm)=NA
#            allCPM=melt(corLcpm)
#            allCPM$value[duplicated(allCPM$value)]=NA
#            interIsland[[paste(i1,"vs",i2,sep="_")]]=allCPM[,3]
#            sampVar=allCPM[which(abs(scale(allCPM[,3])) > 3),]
#            if (nrow(sampVar) > 0) {
#                if (i2=="MPI"){
#                    variation[[paste(i1,"vs",i2,sep="_")]]=transform(sampVar, names=paste(sapply(strsplit(as.character(sampVar[,1]),"[-.]"), `[`, 3), sapply(strsplit(as.character(sampVar[,2]),"[-.]"), `[`, 2), sep="vs"))
#                } else {
#                    variation[[paste(i1,"vs",i2,sep="_")]]=transform(sampVar, names=paste(sapply(strsplit(as.character(sampVar[,1]),"[-.]"), `[`, 3), sapply(strsplit(as.character(sampVar[,2]),"[-.]"), `[`, 3), sep="vs"))
#        }
#        }
#        }
#    }
#    boxplot(interIsland, main=paste("Inter-Island Variation",method,sep="\n"))
#    stripchart(interIsland, vertical=T, method = "jitter", add = TRUE, pch = 20, cex=2, col=c(4,5,7))
#    sapply(1:3, function(x) text(x=x, y=variation[[x]][,3], labels=as.character(variation[[x]][1:nrow(variation[[x]]),4]), cex=0.8, pos=1))
##
#   # melt all three dataframes
#    withinIsland.melt=melt(withinIsland)[,1]
#    interIsland.melt=melt(interIsland)[,1]

#    meta=list(withinVillage.melt, withinIsland.melt, interIsland.melt)
#    names(meta)=c("withinVillage", "withinIsland", "interIsland")
#    boxplot(meta, main=paste("Sample Correlation",method,sep="\n"))
#    dev.off()
#}

# look for sample outliers from PCA
pca.outliers.final=matrix(nrow=0, ncol=3)

for (i in 1:ncol(pcaresults$x)){
    pca.dim=c()
    outlier.sample=c()
    outlier.zscore=c()
    zscore=scale(pcaresults$x[,i])
    outliers=which(abs(zscore) >= 3)
    if (length(outliers) > 0){
        pca.dim=c(pca.dim, i)
        outlier.sample=c(outlier.sample, names(pcaresults$x[,i][outliers]))
        outlier.zscore=c(outlier.zscore, zscore[outliers])
        pca.outliers=matrix(c(rep(pca.dim,length(outlier.sample)), outlier.sample, outlier.zscore), nrow=length(outlier.sample), ncol=3)
        pca.outliers.final=rbind(pca.outliers.final, pca.outliers)
    }
}
colnames(pca.outliers.final)=c("Pca.dim", "Samples", "Z.score")
write.table(pca.outliers.final, file="sample_outliersInPCA.txt", quote=F, row.names=F)

# Analyse what might be driving variation
pdf("CovariateOutliers_SamplingSite.pdf", height=10, width=15)
for (covariate in colnames(y$samples)[c(3,10,13,16:21)]){
    Boxplot(get(covariate)~Island,data=y$samples, main=covariate, col=1:5)
}
dev.off()

# get total amount of variation for each gene and compare islands
island=c("MPI","MTW","SMB")
rv=list()
i=0
pdf("geneVariance.pdf")
plot(0,0,type="n",xlim=c(0.5,3.5), ylim=c(0,20),  xaxt = 'n', xlab ="", ylab = "log2-CPM",  main ="Gene Varaince by Island")
for (sample in island){
    i=i+1
    rv[[sample]]=rowVars(lcpm[,grepl(sample, colnames(lcpm))])
    names(rv[[sample]])=rownames(lcpm)
    vioplot(rv[[sample]], at=i, add=T, col = c(1:3)[i])
}
dev.off()

# let's look at the correlation of gene variance by all population pairs
pdf("geneVariance_correlationbyIsland_Pearson.pdf")
upper.panel<-function(x, y){
  points(x,y, pch=19)
  r <- round(cor(x, y, method="pearson"), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
}

pairs(rv, lower.panel = NULL, upper.panel = upper.panel, main="Gene Variance")
dev.off()








