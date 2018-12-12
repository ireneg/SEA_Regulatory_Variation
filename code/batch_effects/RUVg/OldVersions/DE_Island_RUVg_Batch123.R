# script created by KSB, 06.06.18
# Perform DE analysing relationship between islands using RUVSeq

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVg")
library(NineteenEightyR)

# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))

# set up colors
colors <- electronic_night(n=5)

# exploratory analyis pf results before normalisation
pdf("no_Normalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F)
legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 14, title="Batch", cex=0.6, border=F, bty="n", col=unique(colors[batch]))
dev.off()

# Normalisation can be performed using median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
set <- betweenLaneNormalization(set, which="upper")

# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# exploratory analysis after normalisation
pdf("upperQuartileNOrmalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F)
legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 14, title="Batch", cex=0.6, border=F, bty="n", col=unique(colors[batch]))
legend(legend=unique(allreplicated), pch=16, x="bottomright", col=unique(colors[allreplicated]), cex=0.6, title="replicate", border=F, bty="n")
dev.off()

# First, construct a matrix specifying the replicates. 
# create matrix filled with -1s 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}
genes <- rownames(y)

# set housekeeping genes and perform RUVg. 
# first, set up control genes using housekeeping genes from the dataset taken from Eisenberg and Levanon, 2003
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

pdf("BatchCorrected_PCA_469HousekeepingGenes.pdf")
for (i in c(1,5,10,50)){
	set1 <- RUVg(set, hkControls, k=i)
	plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[batch], main=i)
	plotPCA(set1, col=colors[allreplicated], cex=1.2, main=i)
}
dev.off()

# now do empirically-derived genes

# set up duplicated dataframe and look at how well replicates sit next to each other
# let's see how the PCAs look going through all the covariates
for (name in covariate.names[c(1:10,17:19)]){
	pdf(paste0("RUVsNormalisation_k5_",name,".pdf"))
	plotPCA(set1, col=as.numeric(get(name)), cex=1.2, pch=as.numeric(batch) + 14, labels=F, main=name)
	legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 14, title="Batch", cex=0.6, border=F, bty="n")
	legend(legend=unique(get(name)), pch=16, x="bottomright", col=unique(as.numeric(get(name))), cex=0.6, title="replicate", border=F, bty="n")
	dev.off()
}


# EdgeR pipeline --------------------------------------------------------------------------------------------------

# now do DE by using the empiricall-derived control set of genes
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
# conversely, you can set this up as: design <- model.matrix(~Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
write.table(topTags(lrt), file="topTags_edgeR.txt")

# Limma pipeline --------------------------------------------------------------------------------------------------

contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

# NOw that we have our DE genes, let's see what the best k is with a Limma pipeline
deGenes=vector()
k=vector()
pdf("k_deGenes_RUVs.pdf")
for (j in 1:ncol(efit)){
	counter=0
	for (i in c(1,2,3,4,5,10,20,50,100)){
		counter=counter+1
		set1 <- RUVs(set, genes, k=i, replicates)
		design <- cbind(model.matrix(~0 + Island), pData(set1)[2:(i+1)])
		colnames(design)=gsub("Island", "", colnames(design))
		colnames(design)[3]="Mappi"
		contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
		v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
		vfit <- lmFit(v, design)
		vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
		efit <- eBayes(vfit)
		topTable <- topTable(efit, coef=j, p.value=0.01, lfc=1, n=Inf)
		deGenes[counter]=nrow(topTable)
		k[counter]=i
	}
	plot(k, deGenes, main=colnames(efit)[j])
}
dev.off()

# for now, we'll continue sticking with a k of 5. Let's see how many DE genes we get setting a lfc of 1 and 0.01 pvalue
dt <- decideTests(efit, p.value=0.01, lfc=1)
summary(dt)

#        SMBvsMTW SMBvsMPI MTWvsMPI
# Down          5       35       52
# NotSig    11356    11274    11255
# Up           17       69       71

# Finally, let's visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects given by John Blischak.
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[(3)]=c("Mappi")
batch.corrected.lcpm <- removeBatchEffect(lcpm, covariates=cbind(pData(set1)$W_1 + pData(set1)$W_2 + pData(set1)$W_3 + pData(set1)$W_4 + pData(set1)$W_5),design=design)

# rename column names of lcpm
colnames(lcpm)=samplenames

# plot and see how it looks. Technical replicates should ideally group together
# define pca plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        points(pca$x[,pca_axis1][which(colnames(lcpm) %in% allreps)], pca$x[,pca_axis2][which(colnames(lcpm) %in% allreps)], col="black", pch=8, cex=2)
        text(pca$x[,pca_axis1][which(colnames(lcpm) %in% allreps)], pca$x[,pca_axis2][which(colnames(lcpm) %in% allreps)], labels=samplenames[which(colnames(lcpm) %in% allreps)], pos=3)
        #legend(legend=unique(sampleNames), col=unique(speciesCol), pch=unique(namesPch), x="bottomright", cex=0.6)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 15, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

for (name in covariate.names[c(1:10,17,18,19)]){
  if (nlevels(get(name)) < 26){
    pdf(paste0("batchCorrectedPCA_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=y$samples$batch + 15,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("batchCorrectedPCA_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# this result is a bit strange. The replicates are definitely not as close as in the RUV EDASeq plots