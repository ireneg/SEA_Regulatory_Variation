# script created by KSB, 06.06.18
# Perform batch correction testing whether RUVSeq or linear model method works better

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/lmRUVcomp")

# Batch Removal ----------------------------------------------------------------------------------------------------

# first, set up regular linear modelling pipeline
# set up design
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7)]=c("Mappi","Age","RIN")
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
# It looks like cyclicloess normalisation outputs data with the least amount of varaince so we'll go ahead and use this as our normalisation method
v.lm <- voom(y, design, plot=FALSE, normalize.method = "cyclicloess")
# fit linear models
vfit.lm <- lmFit(v.lm, design)
vfit.lm <- contrasts.fit(vfit.lm, contrasts=contr.matrix)
efit.lm <- eBayes(vfit.lm)
dt.linmod <- decideTests(efit.lm,adjust.method="BH",p.value=0.05)

## RUVSeq

# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])

# first, set up control genes using housekeeping genes from the dataset taken from Eisenberg and Levanon, 2003
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]
# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))
# Normalisation can be performed using median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
set <- betweenLaneNormalization(set, which="upper")

# we'll set k to 5 for now
housekeepingGenes <- RUVg(set, hkControls, k=5)

# set up ruvseq DE function:
ruvseq.de=function(set){
    design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set))
    colnames(design)=gsub("Island", "", colnames(design))
    colnames(design)[3]="Mappi"
    z <- DGEList(counts=counts(set), group=Island)
    z <- calcNormFactors(z, method="upperquartile")
    z <- estimateGLMCommonDisp(z, design)
    z <- estimateGLMTagwiseDisp(z, design)
    contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
    v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
    # fit linear models
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    dt <- decideTests(efit, adjust.method="BH", p.value=0.05)
    return(dt)
}

# find top DE genes for hk genes
dt.ruv.hk <- ruvseq.de(housekeepingGenes)

# Now let's get empirical control genes- first we'll test it out with RUVs way of getting empirical genes, which is performing DE and getting all biut the 5000 top DE genes
design <- model.matrix(~Island, data=pData(set))
z <- DGEList(counts=counts(set), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
top <- topTags(lrt, n=nrow(set))$table
empirical.leastDE <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]


# Choosing a k of 5 results in clear grouping by batch, whereas a k of 50 removes this grouping and technical replicates sit much closer together. We'll choose a k of 10 for now.
leastDEgenes <- RUVg(set, empirical.leastDE, k=5)
# now do DE by using the empiricall-derived control set of genes
dt.ruv.leastDE <- ruvseq.de(leastDEgenes)

# now let's do DE with empirical controls that are obtained using the least variable genes
design <- model.matrix(~Island, data=pData(set))
z <- DGEList(counts=counts(set), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
lcpm=cpm(z, log=TRUE)

# get 6000 least variable genes
batch1.var=sort(apply(lcpm[,grep("firstBatch", colnames(lcpm))], 1, var))[1:6000]
batch2.var=sort(apply(lcpm[,grep("secondBatch", colnames(lcpm))], 1, var))[1:6000]
batch3.var=sort(apply(lcpm[,grep("thirdBatch", colnames(lcpm))], 1, var))[1:6000]
batch.variability=list(batch1=names(batch1.var), batch2=names(batch2.var), batch3=names(batch3.var))
pdf("venn_all3batches.pdf")
venn(batch.variability)
dev.off()

# set empirical control genes derived from least variable genes from all three batches
empirical.var=Reduce(intersect, batch.variability)

# apply RUVSeq
empirical <- RUVg(set, empirical.var, k=5)

# now do DE 
dt.ruv.var <- ruvseq.de(empirical)

# now compare all three methods with vennDiagram
pdf("VennDiagram_allDE_lmRUVhk.pdf")
sapply(1:3, function(x) vennDiagram(cbind(dt.linmod[,x],dt.ruv.hk[,x], dt.ruv.leastDE[,x], dt.ruv.var[,x]), circle.col=c("red", "green", "blue", "yellow"), names=c("linearModel", "housekeeping", "leastDE.emp", "leastVar.emp"), main=colnames(dt.linmod)[x]))
dev.off()

# find out how the logFold change looks like after subtracting common genes 
ruvseq.logFC.plot=function(set, main=main) {
    design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set))
    colnames(design)=gsub("Island", "", colnames(design))
    colnames(design)[3]="Mappi"
    z <- DGEList(counts=counts(set), group=Island)
    z <- calcNormFactors(z, method="upperquartile")
    z <- estimateGLMCommonDisp(z, design)
    z <- estimateGLMTagwiseDisp(z, design)
    contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
    v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
    # fit linear models
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    topTable <- topTable(efit, coef=1, n=Inf, p.value=0.05)
    a=dt.linmod[,1] & dt.ruv.var[,1] & dt.ruv.leastDE[,1] & dt.ruv.hk[,1]
    a=which(a == TRUE)
    topTable=topTable[-c(which(rownames(topTable) %in% names(a))),]
    plot(density(topTable$logFC), col=1, main=main, xlab="LogFC", ylab="Density")
    abline(v=c(-1,-0.5,0.5,1), lty=3)
    for (i in 2:ncol(fit)) {
        topTable <- topTable(efit, coef=i, n=Inf, p.value=0.05)
        a=dt.linmod[,i] & dt.ruv.var[,i] & dt.ruv.leastDE[,i] & dt.ruv.hk[,i]
        a=which(a == TRUE)
        topTable=topTable[-c(which(rownames(topTable) %in% names(a))),]
        lines(density(topTable$logFC), col=i, xlim=c(-2,2))
    }   
    legend(x="topright", bty="n", col=1:ncol(efit), legend=colnames(efit), lty=1, lwd=2)
}

for (method in c("housekeepingGenes", "leastDEgenes", "empirical")){
    pdf(paste0("logFC_minusCoreGenes_",method,".pdf"))
    ruvseq.logFC.plot(get(method), main=method)
    dev.off()
}

# also plot for linear model method
topTable <- topTable(efit.lm, coef=1, n=Inf, p.value=0.05)
a=dt.linmod[,1] & dt.ruv.var[,1] & dt.ruv.leastDE[,1] & dt.ruv.hk[,1]
a=which(a == TRUE)
topTable=topTable[-c(which(rownames(topTable) %in% names(a))),]
pdf("logFC_minusCoreGenes_linearModel.pdf")
plot(density(topTable$logFC), col=1, main="linear model", xlab="LogFC", ylab="Density")
abline(v=c(-1,-0.5,0.5,1), lty=3)
for (i in 2:ncol(efit.lm)) {
    topTable <- topTable(efit.lm, coef=i, n=Inf, p.value=0.05)
    a=dt.linmod[,i] & dt.ruv.var[,i] & dt.ruv.leastDE[,i] & dt.ruv.hk[,i]
    a=which(a == TRUE)
    topTable=topTable[-c(which(rownames(topTable) %in% names(a))),]
    lines(density(topTable$logFC), col=i, xlim=c(-2,2))
}   
legend(x="topright", bty="n", col=1:ncol(efit.lm), legend=colnames(efit.lm), lty=1, lwd=2)
dev.off()


