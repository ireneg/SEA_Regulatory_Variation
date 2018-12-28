# script created by KSB, 08.08.18
# Analyse batch correction by RUVs and set up data for downstream DE analysis

# load dependencies: libraries, human count data and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/RUVs_Setup")

# set up colors
colors <- electronic_night(n=5)
lightcolors=c("thistle", "darkblue")

# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))

# exploratory analyis of results before normalisation
pdf("no_Normalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F)
legend(legend=unique(as.numeric(batch)), "topright", pch=unique(as.numeric(batch)) + 14, title="Batch", cex=0.6, border=F, bty="n", col=unique(colors[batch]))
dev.off()

# Normalisation can be performed using "median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
set <- betweenLaneNormalization(set, which="upper")

# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# exploratory analysis after normalisation
pdf("upperQuartileNOrmalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F, main="Batch")
# plot PCA to see how closely replicates sit
plotPCA(set, col=lightcolors[allreplicated], cex=1.2, main="replicates")
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

# set up genes used for RUVs
genes <- rownames(y)

# set RUVs. We'll first set this to 5 and see how it looks
pdf("RUVsNormalisation_choosingK.pdf")
for (i in c(1,2,3,4,5,10)){
    set1 <- RUVs(set, genes, k=i, replicates)
    plotRLE(set1, main=paste("Housekeeping Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
    # plot PCA to see how closely replicates sit
    plotPCA(set1, col=lightcolors[allreplicated], cex=1.2, main=paste(i, "replicates", sep="\n"))
    # plot pca to see if batch effect is eliminated/minimized
    plotPCA(set1, col=as.numeric(batch), cex=1.2, main="Batch")
}
dev.off()


# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003) for volcano plots
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# look at the p-value distributions and choose the best k
deGenes=vector()
k=vector()
for (j in 1:3){
    counter=0
    pdf(paste0("pvalueDist_choosingK_RUVs_",j,".pdf"))
    for (i in c(1:6)){
        counter=counter+1
        set1 <- RUVs(set, genes, k=i, replicates)
        design <- cbind(model.matrix(~0 + Island), pData(set1)[2:(i+1)])
        colnames(design)=gsub("Island", "", colnames(design))
        colnames(design)[3]="Mappi"
        z <- DGEList(counts=counts(set1), group=Island)
        # here, we use uppserquarile normalisation since "upper" is used for between lane normalisation of our 'set' object
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
    pdf(paste("numberDeGenes_choosingK_RUVs_",colnames(efit)[j],".pdf"))
    plot(k, deGenes, main=colnames(efit)[j])
    dev.off()
}

# we'll set k to 5 for now
set1 <- RUVs(set, genes, k=5, replicates)

# make a matrix for labels of replicates
xcord=c(-0.14,-0.03, 0.02, 0.02, 0.05)
ycord=c(0.1,-0.09,0.0,0.087,0.03)
textcords=cbind(xcord,ycord)

# now plot all of the covariates and see how the batch correction worked
for (name in c(covariate.names[c(1:10,18)], "allreplicated")){
    pdf(paste0("batchCorrectedPCA_RUVs_",name,".pdf"))
    plotPCA(set1, labels=F, pch=as.numeric(batch) + 14, col=as.numeric(get(name)), main=name)
    # text(textcords, c("SMB-WNG-021", "MPI-381", "SMB-ANK-016", "SMB-ANK-027", "MTW-013"), cex=0.8)
    legend(legend=unique(get(name)), pch=16, x="bottomright", col=unique(as.numeric(get(name))), cex=0.6, title=name, border=F, bty="n")
    legend(legend=unique(y$samples$batch), "topright", pch=unique(as.numeric(batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    dev.off()
}