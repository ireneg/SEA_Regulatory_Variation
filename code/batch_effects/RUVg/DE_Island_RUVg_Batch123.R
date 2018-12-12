# script created by KSB, 06.06.18
# Perform DE analysing relationship between islands using RUVg

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVg")
library(NineteenEightyR)

# set up colors
colors <- electronic_night(n=5)
lightcolors=c("thistle", "darkblue")


# housekeeping genes --------------------------------------------------------------------

# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))
# Normalisation can be performed using median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
set <- betweenLaneNormalization(set, which="upper")

# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/batchEffects/HouseKeepingGenes_RevisitedPaper.txt",header=T)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])

# get housekeeping genes(file sent from Ramyar)
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# initially, we should see clear grouping my batch (if the housekeeping genes are neutral)
pdf("hKcontrols_PCA_uncorrected.pdf")
plotPCA(set[hkControls,],col=as.numeric(batch))
dev.off()

# let's decide on which k we need to use by checking the variabilty in each RLE plot between varying amounts of k, as well as how closely our replicates sit to one another in PCA plots
pdf("hKcontrols_choosingK.pdf")
for (i in c(1,2,3,4,5,10)){
	set1 <- RUVg(set, hkControls, k=i)
	plotRLE(set1, main=paste("Housekeeping Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
	# plot PCA to see how closely replicates sit
	plotPCA(set1, col=lightcolors[allreplicated], cex=1.2, main=paste(i, "replicates", sep="\n"))
	# plot pca to see if batch effect is eliminated/minimized
	plotPCA(set1, col=as.numeric(batch), cex=1.2, main="Batch")
	# lastly, plot housekeeping genes to see if batch effect is minimized
	plotPCA(set1[hkControls,], col=as.numeric(batch), cex=1.2, main=paste(i, "housekeeping Genes by batch", sep="\n"))
}
dev.off()

# we can also look at the p-value distributions 
deGenes=vector()
k=vector()
for (j in 1:3){
	counter=0
	pdf(paste0("pvalueDist_choosingK_RUVg_housekeeping_",colnames(efit)[j],".pdf"))
	for (i in c(1,2,3,4,5,10,20)){
		counter=counter+1
		set1 <- RUVg(set, hkControls, k=i)
		design <- cbind(model.matrix(~0 + Island), pData(set1)[2:(i+1)])
		colnames(design)=gsub("Island", "", colnames(design))
		colnames(design)[3]="Mappi"
		z <- DGEList(counts=counts(set1), group=Island)
		z <- calcNormFactors(z, method="upperquartile")
		contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
		v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
		vfit <- lmFit(v, design)
		vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
		efit <- eBayes(vfit)
		topTable <- topTable(efit, coef=j, p.value=0.01, lfc=1, n=Inf)
		deGenes[counter]=nrow(topTable)
		k[counter]=i
		# get p-value distribution from raw pvalues
		hist(efit$p.value[,j], main=paste(colnames(efit)[j],i,sep="\n"), ylim=c(0,max(table(round(efit$p.value[,j], 1)))+1000), xlab="p-value")
		# get volcano plots
		plot(efit$coef[,j], -log10(as.matrix(efit$p.value)[,j]), pch=20, main=paste(colnames(efit)[j],i,sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    	points(efit$coef[,j][which(names(efit$coef[,j]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,j][which(names(efit$coef[,j]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    	legend("topleft", c("genes", "hk genes"),fill=c("black",4))
    	abline(v=c(-1,1), lty=2)
	}
	dev.off()
	pdf(paste("numberDeGenes_choosingK_RUVg_housekeeping_",colnames(efit)[j],".pdf"))
	plot(k, deGenes, main=colnames(efit)[j])
	dev.off()
}

# we'll set k to 5 for now
set1 <- RUVg(set, hkControls, k=5)

# let's do DE with our housekeeping genes and see what we get
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit.hk <- eBayes(vfit)
# get number of DE genes that pass our theshold
dt.hk <- decideTests(efit.hk, p.value=0.01, lfc=1)
# get summary of decide tests statistics
write.table(summary(dt.hk), file="numberSigDEgenes_RUVg_housekeepingGenes_voom.txt")
# for all contrasts, get names of significant genes
for (i in 1:ncol(efit.hk)){
	topTable <- topTable(efit.hk, coef=i, p.value=0.01, lfc=1, n=Inf)
	sig.hkGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
	write.table(sig.hkGenes, file=paste0("significantGenes_",colnames(efit.hk)[i],".txt"))
	# make mdplot
	o <- which(names(efit.hk$Amean) %in% names(which(abs(dt.hk[,i])==1)))
    x <- efit.hk$Amean
    z <- efit.hk$coefficients[,i]
    t=which(names(efit.hk$coefficients[,i]) %in% names(which(abs(dt.hk[,i])==1)))
    G <- y$genes[names(which(dt.hk[,i]==1)),]$SYMBOL
    pdf(paste0("MDPlot_housekeeping_TopGenes_FDRpval01_LFC01_",colnames(efit.hk)[i],".pdf"), height=15,width=13)
    plotMD(efit.hk, column=i, status=dt.hk[,i], main=colnames(efit.hk)[i], hl.col=c("blue","red"), values=c(-1,1))
    abline(h=c(1,-1), lty=2)
    legend(legend=paste(names(summary(dt.hk)[,i]), summary(dt.hk)[,i], sep="="), x="bottomright", border=F, bty="n")
    text(x[o], z[t], labels=G)
    dev.off()
}

# We can also visualise how RUVg batch correction looks like with a PCA. First, turn the seqExpressionSet object into a summarized experiment

se.hk=makeSummarizedExperimentFromExpressionSet(set1)

# set up PCa [plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(assay(dataToPca)), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        points(pca$x[,pca_axis1][which(colnames(lcpm) %in% allreps)], pca$x[,pca_axis2][which(colnames(lcpm) %in% allreps)], col="black", pch=8, cex=2)
        text(pca$x[,pca_axis1][which(samplenames %in% allreps)], pca$x[,pca_axis2][which(samplenames %in% allreps)], labels=samplenames[which(samplenames %in% allreps)], pos=3)
        #legend(legend=unique(sampleNames), col=unique(speciesCol), pch=unique(namesPch), x="bottomright", cex=0.6)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 15, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

for (name in c(covariate.names[c(1:10,18)], "allreplicated")){
  if (nlevels(get(name)) < 26){
    pdf(paste0("batchCorrectedPCA_RUVg_housekeeping_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=se.hk, speciesCol=as.numeric(get(name)),namesPch=y$samples$batch + 15,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("batchCorrectedPCA_RUVg_housekeeping_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=se.hk, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# empirically-derived control genes---------------------------------------------------------

# Let's see where our replicates sit when we use empirically-derived control genes. We can use empirically-derived control genes by using all but the top 5000 genes. First, we need to set up our design matrix and perform DE in order to find the least significant DE genes
# reset design back to the normalised seq expression set
design <- model.matrix(~Island, data=pData(set))
z <- DGEList(counts=counts(set), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
# set up top tags for all 3 coefficients
lrt1 <- glmLRT(fit, coef=1)
lrt2 <- glmLRT(fit, coef=2)
lrt3 <- glmLRT(fit, coef=3)
top1 <- topTags(lrt1 , n=nrow(set))$table
top2 <- topTags(lrt2 , n=nrow(set))$table
top3 <- topTags(lrt3 , n=nrow(set))$table
a=rownames(top1)[1:7000]
b=rownames(top2)[1:7000]
c=rownames(top2)[1:7000]
mergedTop <- intersect(intersect(a,b),c)
empirical <- rownames(set)[which(!(rownames(set) %in% mergedTop))]


# define new empirical gene set after getting top DE genes. First, let's see what k should be. 
pdf("empiricalControls_choosingK.pdf")
for (i in c(1,2,3,4,5,10)){
	set2 <- RUVg(set, empirical, k=i)
	plotRLE(set2, main=paste("Empirical Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
	plotPCA(set1, col=lightcolors[allreplicated], cex=1.2, main=i)
}
dev.off()

# we can also look at the p-value distributions 
deGenes=vector()
k=vector()
for (j in 1:3){
	counter=0
	pdf(paste0("pvalueDist_choosingK_RUVg_empirical_",colnames(efit)[j],".pdf"))
	for (i in c(1,2,3,4,5,10,20)){
		counter=counter+1
		set2 <- RUVg(set, empirical, k=i)
		design <- cbind(model.matrix(~0 + Island), pData(set2)[2:(i+1)])
		colnames(design)=gsub("Island", "", colnames(design))
		colnames(design)[3]="Mappi"
		z <- DGEList(counts=counts(set2), group=Island)
		z <- calcNormFactors(z, method="upperquartile")
		contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
		v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
		vfit <- lmFit(v, design)
		vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
		efit <- eBayes(vfit)
		topTable <- topTable(efit, coef=j, p.value=0.01, lfc=1, n=Inf)
		deGenes[counter]=nrow(topTable)
		k[counter]=i
		# get p-value distribution from raw pvalues
		hist(efit$p.value[,j], main=paste(colnames(efit)[j],i,sep="\n"), ylim=c(0,max(table(round(efit$p.value[,j], 1)))+1000), xlab="p-value")
		# get volcano plots
		plot(efit$coef[,j], -log10(as.matrix(efit$p.value)[,j]), pch=20, main=paste(colnames(efit)[j],i,sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    	points(efit$coef[,j][which(names(efit$coef[,j]) %in% empirical)], -log10(as.matrix(efit$p.value)[,j][which(names(efit$coef[,j]) %in% empirical)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    	legend("topleft", c("genes", "empirical genes"),fill=c("black",4))
    	abline(v=c(-1,1), lty=2)
	}
	dev.off()
	pdf(paste("numberDeGenes_choosingK_RUVg_empirical_",colnames(efit)[j],".pdf"))
	plot(k, deGenes, main=colnames(efit)[j])
	dev.off()
}

# a k of 5 looks to be the best
set2 <- RUVg(set, empirical, k=5)

# let's do DE with our housekeeping genes and see what we get
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set2))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set2), group=Island)
z <- calcNormFactors(z, method="upperquartile")
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit.emp <- eBayes(vfit)
# get number of DE genes that pass our theshold
dt.emp <- decideTests(efit.emp, p.value=0.01, lfc=1)
# get summary of decide tests statistics
write.table(summary(dt.emp), file="numberSigDEgenes_RUVg_empiricalGenes_voom.txt")
# for all contrasts, get names of significant genes
for (i in 1:ncol(efit.emp)){
	topTable <- topTable(efit.emp, coef=i, p.value=0.01, lfc=1, n=Inf)
	sig.empGenes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
	write.table(sig.empGenes, file=paste0("significantGenes_",colnames(efit.emp)[i],".txt"))
	# make mdplot
	o <- which(names(efit.emp$Amean) %in% names(which(abs(dt.emp[,i])==1)))
    x <- efit.emp$Amean
    z <- efit.emp$coefficients[,i]
    t=which(names(efit.emp$coefficients[,i]) %in% names(which(abs(dt.emp[,i])==1)))
    G <- y$genes[names(which(dt.emp[,i]==1)),]$SYMBOL
    pdf(paste0("MDPlot_housekeeping_TopGenes_FDRpval01_LFC01_",colnames(efit.emp)[i],".pdf"), height=15,width=13)
    plotMD(efit.emp, column=i, status=dt.emp[,i], main=colnames(efit.emp)[i], hl.col=c("blue","red"), values=c(-1,1))
    abline(h=c(1,-1), lty=2)
    legend(legend=paste(names(summary(dt.emp)[,i]), summary(dt.emp)[,i], sep="="), x="bottomright", border=F, bty="n")
    text(x[o], z[t], labels=G)
    dev.off()
}

# we can also visualise empirical control genes with a PCA
se.emp=makeSummarizedExperimentFromExpressionSet(set2)

for (name in c(covariate.names[c(1:10,18)], "allreplicated")){
  if (nlevels(get(name)) < 26){
    pdf(paste0("batchCorrectedPCA_RUVg_empirical",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=se.emp, speciesCol=as.numeric(get(name)),namesPch=y$samples$batch + 15,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("batchCorrectedPCA_RUVg_empirical",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=se.emp, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# Get which genes overlap between RUVs and linear model
pdf("VennDiagram_hkVSempComparison.pdf")
sapply(1:3, function(x) vennDiagram(cbind(dt.hk[,x],dt.emp[,x]), circle.col=c("red","blue"), names=c("hk", "empirical"), main=colnames(dt.hk)[x]))
dev.off()

decideTestsDGE

# get all of the genes that inersect
for (i in 1:3){
	commonGenes=dt.hk[,i] & dt.emp[,i]
	commonGenes=which(commonGenes == TRUE)
	commonGenes.results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=names(commonGenes), filters="ensembl_gene_id")
	write.table(commonGenes.results,file=paste0("allCommonGenes_empVsHk",colnames(dt.emp)[i],".txt"))
}

