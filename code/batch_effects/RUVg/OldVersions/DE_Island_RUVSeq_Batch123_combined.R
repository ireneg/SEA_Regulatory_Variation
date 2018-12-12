# script created by KSB, 06.06.18
# Perform DE analysing relationship between islands using RUVSeq

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval")

# Batch Removal ----------------------------------------------------------------------------------------------------

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
# work RUV's magic! First, let's decide on which k we need to use
pdf("hKcontrols_choosingK.pdf")
for (i in c(1,2,5,10,50)){
	set1 <- RUVg(set, hkControls, k=i)
	plotRLE(set1, main=paste("Housekeeping Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
}
dev.off()

# set up new colour scheme with batch and replicate information
batch.replicate=as.numeric(batch)
replicate.index=grep("MPI-381_firstBatch|MPI-381_thirdBatch|SMB-ANK-027_firstBatch|SMB-ANK027_thirdBatch", colnames(y))
batch.replicate[replicate.index]="black"

# let's set the k to 1,5, 10, and and see how it looks
pdf("BatchCorrected_PCA_469HousekeepingGenes.pdf")
for (i in c(1,5,10,50)){
	set1 <- RUVg(set, hkControls, k=i)
	plotPCA(set1, col=batch.replicate, main=i)
}
dev.off()

# we'll set k to 5 for now
set1 <- RUVg(set, hkControls, k=5)

# let's dp DE with our housekeeping genes and see what we get
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)

# first, try method with edgeR
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
top.housekeeping=topTags(lrt)

# adjusting k certainly makes a difference. We can see that the technical replicates sit closer together and the amount of variation explained by the first dimeniosn of the PCA changes dependent on k. Is this potentially masking real biology or getting rid of unwanted variation?

# Let's see where our replicates sit when we use empirically-derived control genes. We can use empirically-derived control genes by using all but the top 5000 genes. First, we need to set up our design matrix and perform DE in order to find the least significant DE genes
design <- model.matrix(~Island, data=pData(set))
z <- DGEList(counts=counts(set), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
# TODO: question: which coefficient should I be using??
lrt <- glmLRT(fit, coef=3)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
# define new empirical gene set after getting top DE genes. First, let's see what k should be. 
pdf("empiricalControls_choosingK.pdf")
for (i in c(1,2,5,10,50)){
	set2 <- RUVg(set, empirical, k=i)
	plotRLE(set2, main=paste("Empirical Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
}
dev.off()

# Let's see what a k of 50 and 5 look like in the PCA
pdf("BatchCorrected_PCA_empiricallyDerivedControlGenes.pdf")
for (i in c(1,5,10,50)){
	set2 <- RUVg(set, empirical, k=i)
	plotPCA(set2, col=batch.replicate, cex=1.2, main=i)
}
dev.off()

# Choosing a k of 5 results in clear grouping by batch, whereas a k of 50 removes this grouping and technical replicates sit much closer together. We'll choose a k of 10 for now.
set2 <- RUVg(set, empirical, k=5)

# now do DE by using the empiricall-derived control set of genes
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set2))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"

z <- DGEList(counts=counts(set2), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)

# first, try method with edgeR
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
top.empirical=topTags(lrt)

# let's see how if the top genes in our housekeeping control set vs the emoirically-derived control genes come up with the same top genes:
rownames(top.housekeeping) %in% rownames(top.empirical)

# [1]  TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE

# Looks as if some are and some aren't (60%). So in summary, it looks as though adjusting your control genes makes a difference, as well as adjusting k. 

# let's also find out what these genes are doing
topGenes.edgeR <- getBM(attributes = c('external_gene_name', 'description'), mart = ensembl.mart.90,values=rownames(top.empirical), filters="ensembl_gene_id")
print(topGenes.function[,1])

#  [1] "MARCO"      "MYO1B"      "ADAMTS1"    "PCDH1"      "CAPN12"    
#  [6] "TCL6"       "ZNF677"     "AC104389.1" "AC092802.3" "AL691447.2"

# Limma pipeline --------------------------------------------------------------------------------------------------

# Let's try the method now with Limma
# note that batch effect removel using RUVSeq with limma has not been tested, as per this Bioconductor thread: https://support.bioconductor.org/p/86461/

# set up contrast matrix first with empirical genes
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set2))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=TRUE, normalize.method = "quantile")
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
dt <- decideTests(efit)
topTable <- topTable(efit, coef=3, p.value=0.05, n=10)

# do we see the same results in Limma's top genes as in EdgeR's top empirical genes?
print(rownames(topTable) %in% rownames(top.empirical))

# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

# let's get the top DE gene names:
topGenes.voom.empirical <- getBM(attributes = c('external_gene_name', 'description'), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
print(topGenes.voom.empirical[,1])

#  [1] "STAT4"     "NEDD1"     "UBASH3B"   "SAMD3"     "SRRM2"     "RSL1D1"   
#  [7] "TNFRSF10D" "NSUN3"     "MTF1"      "ARHGAP19" 

# now let's test housekeeping genes
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=TRUE, normalize.method = "quantile")
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
dt <- decideTests(efit)
topTable <- topTable(efit, coef=3, p.value=0.05, n=10)

# do we see the same results in Limma's top genes as in EdgeR's top hk genes?
print(rownames(topTable) %in% rownames(top.housekeeping))

# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

# let's get the top DE gene names:
topGenes.voom.hk <- getBM(attributes = c('external_gene_name', 'description'), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
print(topGenes.voom.hk[,1])

#  [1] "MATK"     "PLEK"     "SYTL2"    "STAT4"    "GRAP"     "AUTS2"   
#  [7] "SAMD3"    "SIGLEC7"  "ARHGAP19" "APOBEC3G"

# Doesn't seem like a limma-voom pipeline yields the same results as EdgeR. What if we tried a more simple Limma-trend pipeline (wthout the weights)

# limma trend method
logCPM <- cpm(z, log=TRUE, prior.count=0.25)
fit <- lmFit(logCPM, design)
fit2 <- contrasts.fit(fit, contrasts=contr.matrix)
trend <- eBayes(fit2, trend=TRUE)
topTable <- topTable(efit, coef=3, p.value=0.05, n=10)

# do we see any similarities between trend and EdgeR
rownames(topTable) %in% rownames(top.housekeeping)

# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

# get top trend DE gene names
topGenes.trend <- getBM(attributes = c('external_gene_name', 'description'), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
print(topGenes.trend[,1])

# [1] "MATK"     "PLEK"     "SYTL2"    "STAT4"    "GRAP"     "AUTS2"   
# [7] "SAMD3"    "SIGLEC7"  "ARHGAP19" "APOBEC3G"


