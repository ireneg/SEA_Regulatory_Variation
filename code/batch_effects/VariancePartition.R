# script created by KSB, 08.08.18
# Analyse the amount of variance each covariate contributes to the design matrix

### Last edit: KB 05.04.2019

# Load dependencies and set input paths --------------------------

# load libraries-
library(variancePartition)
library(pvca)
library(doParallel)
library(edgeR)
library(magrittr)
library(RUVSeq)
library(Homo.sapiens)
library(ggubr)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVvsLinearModel/"
# Load the DGE list object. y:
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))

# set up parallel processing
cl <- makeCluster(4)
registerDoParallel(cl)

# Variance partition per variable ---------------------------------
# First, using Limma implementing known covariates 

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
# set up design
design <- model.matrix(~0 + y$samples$Island + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$NK + y$samples$Bcell + y$samples$Mono + y$samples$Gran)
colnames(design)=gsub("Island", "", colnames(design)) %>% gsub("[\\y$]", "", .) %>% gsub("samples", "", .) %>% gsub("West Papua", "Mappi", .)
v <- voom(y, design, plot=F)

# Specify variables to consider
# Age is continuous so model it as a fixed effect. Individual and Tissue are both categorical, so model them as random effects.
# Note the syntax used to specify random effects
form1 <- ~ Age + RIN + CD8T +  CD4T + NK + Bcell + Mono + Gran + (1|Island) + (1|batch)

# Fit model and extract results
# 1) fit linear mixed model on gene expression. If categorical variables are specified,
# a linear mixed model is used. If all variables are modeled as fixed effects,
# a linear model is used. Each entry results in a regression model fit on a single gene.
# 2) extract variance fractions from each model fit for each gene, the fraction of variation attributable
# to each variable is returned
# Interpretation: the variance explained by each variables after correcting for all other variables
varPart1 <- fitExtractVarPartModel(v, form1, y$samples)

# sort variables (i.e. columns) by median fraction of variance explained
vp1 <- sortCols(varPart1)

# violin plot of contribution of each variable to total variance
fig=plotVarPart(vp1, main="Limma")
ggsave(file=paste0(outputdir,"totalVarianceContribution_lmModel_allVariables.pdf"), fig)

write.table(summary(vp1), file=paste0(outputdir,"summary_lmModel_VariancePartition.txt"))

# Now look at RUVs ----------------------------------------------------------------------------------------------------

# identify which samples are replicated
load(paste0(inputdir, "covariates.Rda"))
allreps=covariates[,1][which(covariates$replicate)]
allreps=unique(allreps)
allreplicated=as.factor(samplenames %in% allreps)

allreplicated=as.factor(samplenames %in% allreps)

# define sample names
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# construct a matrix specifying the replicates. 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}
genes <- rownames(y)

# set up pheno data with all known factors of unwanted variation
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(y$samples$Island, row.names=colnames(y)))
# normalise with upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")
set1 <- RUVs(set, genes, k=5, replicates)

# create design matrix
design <- model.matrix(~0 + y.samples.Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("y.samples.Island", "", colnames(design))
colnames(design)=gsub("West Papua", "Mappi", colnames(design))

z <- DGEList(counts=counts(set1), group=y$samples$Island)
z <- calcNormFactors(z, method="upperquartile")

# set up gene names
geneid <- rownames(z)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
z$genes <- genes

# make contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
v.ruv <- voom(z, design, plot=FALSE, normalize.method = "cyclicloess")

# Define formula
form2 <- ~ W_1 + W_2 + W_3 + W_4 + W_5 + (1|y.samples.Island)
# fit modelf
varPart2 <- fitExtractVarPartModel(v.ruv, form2, pData(set1))

# sort variables (i.e. columns) by median fraction of variance explained
vp2 <- sortCols(varPart2)

# violin plot of contribution of each variable to total variance
fig2=plotVarPart(vp2, main="RUVs")
ggsave(file="totalVarianceContribution_RUVs_allVariables.pdf", fig2)
write.table(summary(vp2), file="summary_RUVs_VariancePartition.txt")

# Merge both and plot
pdf(file="VarianceExplained_LMvsRUVs.pdf", height=10, width=15)
ggarrange(fig, fig2, labels=c("A","B"))
dev.off()

# Now estimate the total variablity due to batch effects using PVCA -----------------------------------------------------

#data.lm=AnnotatedDataFrame(y$samples[,c(5,6,10,13,16:21)])
#eset.X.lm=new("ExpressionSet", exprs=lcpm, phenoData=data.lm)
#pvcaObj.lm <- pvcaBatchAssess(eset.X.lm, colnames(data.lm), 0.1)

#data.RUV=AnnotatedDataFrame(pData(set1))
#eset.X.RUV=new("ExpressionSet", exprs=normCounts(set1), phenoData=data.RUV)
#pvcaObj.RUVs <- pvcaBatchAssess(eset.X.RUV, c("Island", "W_1"), 0.6)


#bp <- barplot(pvcaObj.before$dat,  xlab = "Effects",
#         ylab = "Weighted average proportion variance",
#         ylim= c(0,1.1),col = c("blue"), las=2,
#         main="LM Variance Estimation")
##axis(1, at = bp, labels = pvcaObj.before$label, xlab = "Effects", cex.axis = 0.5, las=2)
#values = pvcaObj$dat
#new_values = round(values , 3)
#text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)

