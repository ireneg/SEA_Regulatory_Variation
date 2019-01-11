# script created by KSB, 08.08.18
# Look at effectiveness of batch correction by regressing out covariates

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/linearModelling")

# DE analysis ----------------------------------------------------------------------------------------------------

# set up design
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# As per JOhn Blischak's recommendations, see if my variable of interest (Island) is confounded with any of technical variables:
pdf("Correlation_AgeIsland.pdf")
boxplot(y$samples$Age ~ y$samples$Island, main=c("Age", paste0("r2 = ", round(summary(lm(y$samples$Age ~ y$samples$Island))$adj.r.squared, digits=2))))
dev.off()

# get r2 between variables
summary(lm(y$samples$Age ~ y$samples$Island))$adj.r.squared

# 0.2153212 - looks like the correlation is pretty low

pdf("Correlation_RinIsland.pdf")
boxplot(y$samples$RIN ~ y$samples$Island, main=c("RIN", paste0("r2 = ", round(summary(lm(y$samples$RIN ~ y$samples$Island))$adj.r.squared, digits=2))))
dev.off()

summary(lm(y$samples$RIN ~ y$samples$Island))$adj.r.squared

# 0.1862963 - once again, low correlation

table(batch, y$samples$Island)

# batch 1 and 2 look fine, however batch 3 is confounded with Sumba

fisher.test(batch, y$samples$Island)

# p-value = 3.628e-06
# alternative hypothesis: two.sided

# set up contrast matrix using nested design
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))

# It looks like cyclicloess normalisation outputs data with the least amount of varaince so we'll go ahead and use this as our normalisation method
v <- voom(y, design, plot=FALSE, normalize.method = "cyclicloess")
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

# if p-value distributions look realistic (which they seem to), continue on with DE analysis
dt <- decideTests(efit)

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects given by John Blischak.
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[3]=c("Mappi")
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=batch, covariates = cbind(y$samples$Age, y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran),design=design)

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

for (name in covariate.names[c(1:10,17:19)]){
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

# run removeBatchEffect without protecting variable of interest:
x2 <- removeBatchEffect(batch.corrected.lcpm, batch=batch,  covariates = cbind(y$samples$Age, y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran))

for (name in covariate.names[c(1:10,17:19)]){
  if (nlevels(get(name)) < 26){
    pdf(paste0("batchCorrectedPCA_",name,"_unprotectedVariable.pdf"))
    pcaresults <- plot.pca(dataToPca=x2, speciesCol=as.numeric(get(name)),namesPch=y$samples$batch + 15,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0("batchCorrectedPCA_",name,"_unprotectedVariable.pdf"))
    pcaresults <- plot.pca(dataToPca=x2, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}








