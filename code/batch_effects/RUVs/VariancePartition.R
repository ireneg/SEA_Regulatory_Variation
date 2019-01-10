source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# load libraries-
library(variancePartition)
library(pvca)

# set up parallel processing
cl <- makeCluster(4)
registerDoParallel(cl)

setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/lmRUVcomparison")

# First, we need to test the regular linear model method implementing known covariates -----------------------------------

# Perform TMM normalization
y <- calcNormFactors(y, method="TMM")

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
# set up design
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
# Estimate precision weights for each gene and sample
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
fig=plotVarPart(vp1)
ggsave(file="totalVarianceContribution_lmModel_allVariables.pdf", fig)

write.table(summary(vp1), file="summary_lmModel_VariancePartition.txt")

# Now look at RUVs ----------------------------------------------------------------------------------------------------

# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# construct a matrix specifying the replicates. 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}
genes <- rownames(y)

# set up pheno data with all known factors of unwanted variation
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
v.ruv <- voom(z, design, plot=FALSE, normalize.method = "cyclicloess")

# Define formula
form2 <- ~ W_1 + W_2 + W_3 + W_4 + W_5 + (1|Island)
# fit modelf
varPart2 <- fitExtractVarPartModel(v.ruv, form2, pData(set1))

# sort variables (i.e. columns) by median fraction of variance explained
vp2 <- sortCols(varPart2)

# violin plot of contribution of each variable to total variance
fig2=plotVarPart(vp2)
ggsave(file="totalVarianceContribution_RUVs_allVariables.pdf", fig2)

write.table(summary(vp2), file="summary_RUVs_VariancePartition.txt")

# Now estimate the total variablity due to batch effects using PVCA -----------------------------------------------------
# batch.dataframe=y$samples[,c(5,6,10,13,16:21)]
# batch.dataframe[,3:10] <- lapply(batch.dataframe[,3:10], factor)

# a=AnnotatedDataFrame(batch.dataframe)
# b=new("ExpressionSet", exprs=lcpm, phenoData=a)
# batch.factors <- colnames(a)
# pct_threshold <- 0.6
# pvcaObj <- pvcaBatchAssess(b, batch.factors, pct_threshold) 

# note numeric doesn't work, must do as.factor

# bp <- barplot(pvcaObj$dat,  xlab = "Effects",
#         ylab = "Weighted average proportion variance",
#         ylim= c(0,1.1),col = c("blue"), las=2,
#         main="LM Variance Estimation")
#axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.5, las=2)
#values = pvcaObj$dat
#new_values = round(values , 3)
#text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)


# now with -------------------------------
