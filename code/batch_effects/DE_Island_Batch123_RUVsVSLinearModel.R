# script created by KSB, 06.06.18
# Perform DE analysing using RUVs vs a traditional linear modelling approach

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVs/RUVvsLinearModel")
library(NineteenEightyR)

## setup ----------------

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

## Analysis -----------------

# First, set up pheno data with all known factors of unwanted variation
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))
# normalise with upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")
set1 <- RUVs(set, genes, k=5, replicates)

# create design matrix
design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")

# set up gene names
geneid <- rownames(z)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
z$genes <- genes

# make contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v <- voom(z, design, plot=FALSE, normalize.method = "cyclicloess")
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
# for now, we'll continue sticking with a k of 5. Let's see how many DE genes we get setting a lfc of 1 and 0.01 pvalue
dt <- decideTests(efit, p.value=0.01, lfc=1)
summary(dt)

#        SMBvsMTW SMBvsMPI MTWvsMPI
# Down          5       35       52
# NotSig    11356    11274    11255
# Up           17       69       71

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
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
# set up contrast matrix using nested design
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v.lm <- voom(y, design, plot=FALSE, normalize.method = "quantile")
# fit linear models
vfit.lm <- lmFit(v.lm, design)
vfit.lm <- contrasts.fit(vfit.lm, contrasts=contr.matrix)
efit.lm <- eBayes(vfit.lm)
dt.lm <- decideTests(efit.lm, p.value=0.01, lfc=1)
summary(dt.lm)

#        SMBvsMTW SMBvsMPI MTWvsMPI
# Down          1       62       70
# NotSig    11364    11203    11237
# Up           13      113       71

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


# EdgeR pipeline --------------------------------------------------------------------------------------------------
# now let's see what EdgeR is doing


design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))

conversely, you can set this up as: design <- model.matrix(~Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
z <- calcNormFactors(z, method="upperquartile")
z <- estimateGLMCommonDisp(z, design)
z <- estimateGLMTagwiseDisp(z, design)
fit <- glmFit(z, design)
lrt <- glmLRT(fit, coef=3)
write.table(topTags(lrt), file="topTags_edgeR.txt")

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









