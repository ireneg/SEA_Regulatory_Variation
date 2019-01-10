# script created by KSB, 26.07.18
# perform data pre-processing on both batches of indonesian RNA-seq data- third batch

# first laod dependencies (i.e., count data)
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")

# Set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing")


# covariate setup ------------------------------------------------------

# reorder covariates by sequencing batch
covariates=covariates[order(covariates$Sequencing.Batch),]
all(samplenames == covariates[,1])
# TRUE

# Assign covariates to DGE list
covariate.names = c(colnames(covariates)[c(2,3,5:12,14:20)], "lib.size", "batch")
for (name in covariate.names[1:17]){
  y$samples[[paste0(name)]]<- covariates[,name]
}
# make sequencing pool into a factor instead of numeric
y$samples$Sequencing.pool=as.factor(y$samples$Sequencing.pool)

## Covariate assignment: below, we want to assign all of our covariates to variables.
# Age, RIN, and library size need to be broken up into chunks for easier visualisation of trends (for instance in Age, we want to see high age vs low age rather than the effect of every single age variable)
Age <- cut(as.numeric(as.character(y$samples$Age)), c(14,24,34,44,54,64,74,84), labels=c("15-24","25-34", "35-44", "45-54", "55-64", "65-74", "75-84"))
RIN <- cut(as.numeric(as.character(y$samples$RIN)), c(4.9,5.9,6.9,7.9,8.9), labels=c("5.0-5.9", "6.0-6.9", "7.0-7.9", "8.0-8.9"))
lib.size <- cut(as.numeric(y$samples$lib.size), c(8000000,12000000,16000000,20000000,24000000), labels=c("8e+06-1.2e+07","1.2e+07-1.6e+07", "1.6e+07-2e+07", "2e+07-4.4e+07"))
# assign blood type information to variables
for (name in covariate.names[c(11:16)]){
assign(name, as.numeric(as.character(y$samples[[paste0(name)]])))
}

# Now do the rest- get which samples are not numeric in y$samples and assign values to variables
nums=unlist(lapply(y$samples, is.numeric))
# this looks confusing, but here we're getting variables that aren't numeric and subtracting the first two variables, which are files and group (which we don't need)
covar.variables=c(colnames(y$samples[,!nums])[3:ncol(y$samples[,!nums])], "batch")
        
# assign values to other covariate names
for (name in covar.variables){
  if(sum(is.na(y$samples[[name]])) > 0){
    y$samples[[paste0(name)]]=assign(name, as.factor(addNA(y$samples[[paste0(name)]])))
  }
  else{
    y$samples[[paste0(name)]]=assign(name, as.factor(y$samples[[paste0(name)]]))
  }
}

# Perform rarefaction curves for number of expressed genes vs. proportion of pool mRNA (as in Ramskold D, Wang ET, Burge CB, Sandberg R. 2009. An abundance of ubiquitously expressed genes revealed by tissue transcriptome sequence data. PLoS Comput Biol 5:e1000598)
pdf("rarefactionCurves_indoRNA_123Combined.pdf")
for (name in covariate.names[c(1:10,17,18)]) {
  plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=name) ## initialize the plot area
  counter=0
  for (sample in colnames(y)){
    counter=counter+1
    lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=as.numeric(get(name))[counter])
  }
  levels=levels(get(name))
  levels[which(is.na(levels))] = "NA"
  legend(x="bottomright", bty="n", col=1:length(levels(get(name))), legend=levels, lty=1, lwd=2)
}
# add in batch
plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main="batch") ## initialize the plot area
counter=0
for (sample in colnames(y)){
  counter=counter+1
  lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=batch.col[as.numeric(batch)][counter])
}
levels=levels(batch)
levels[which(is.na(levels))] = "NA"
legend(x="bottomright", bty="n", col=batch.col[1:length(levels(batch))], legend=levels, lty=1, lwd=2)
dev.off()

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

# replicate analysis -----------------------------------------------------------------------------------

# reset samplenames 
samplenames=sapply(strsplit(colnames(y),"[_.]"), `[`, 1)
samplenames[104:123]=sapply(samplenames[104:123],function(x)sub("([[:digit:]]{3,3})$", "-\\1", x))

# set replicate names and rename lcpm columns
allreps=covariates[,1][which(covariates$replicate)]
allreps=unique(allreps)

# look at how well replicate performed
pdf("replicate_comparisons_noFiltering_123Combined.pdf", height=10, width=10)
par(mfrow=c(3,3))
# Since SMB-ANK-027 had two replicates, we need to plot one vs three and two vs three. First, one vs three
smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[1]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[1]], xlim=c(-5,15), ylim=c(-5,15), main="Technical replicate SMB-ANK-027")
# best fit/regression line
abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]), col="green")
# diagonal line
abline(a=0,b=1,col="red")

# Now, replicates two vs three
smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[2]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[2]], xlim=c(-5,15), ylim=c(-5,15), main="Technical replicate SMB-ANK-027")
# best fit/regression line
abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]), col="green")
# diagonal line
abline(a=0,b=1,col="red")

for (i in 1:length(allreps)){
  smoothScatter(lcpm[,which(samplenames %in% allreps[i])[1]], lcpm[,which(samplenames %in% allreps[i])[2]], ylab=colnames(lcpm)[which(samplenames %in% allreps[i])[2]], xlab=colnames(lcpm)[which(samplenames %in% allreps[i])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate", allreps[i], sep=" "))
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")
}
dev.off()

# Removal of lowly-expressed genes -------------------------------------------------------------------------------

# a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs <- rowSums(lcpm>1) >= (nrow(y$samples)*0.5)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

# Visualise library size after filtering with barplots
pdf("librarySize_indoRNA_postFiltering_123Combined.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10), main="Library Size \n Post-filtering")
# abline(h=10, col="red")
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()
pdf("nGenes_indoRNA_postFiltering_123Combined.pdf", height=10, width=15)
par(oma=c(2,0,0,0))
barplot(apply(y$counts, 2, function(c)sum(c!=0)),main="Number of Genes \n Post-Filtering", ylab="n Genes", cex.names=0.75, col=batch.col[y$samples$batch], names=samplenames, las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+3000))
legend(x="topright", col=batch.col[unique(y$samples$batch)], legend=c("first batch", "second batch", "third batch"), pch=15, cex=0.8)
dev.off()

# Compare library size density before and after removing lowly-expressed genes
pdf("libraryDensity_afterFiltering_afterNormalization_indoRNA.pdf", height=8, width=15)
nsamples <- ncol(y)
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=batch.col[y$samples$batch][1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.1), las=2, main="", xlab="")
title(main="A. All genes", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=batch.col[y$samples$batch][i], lwd=2)
}
legend("topright", legend=c("First Batch","Second Batch", "Third Batch"), ncol=1, cex=0.8, text.col=batch.col[unique(y$samples$batch)], bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=batch.col[y$samples$batch][1], lwd=2, ylim=c(0,max(density(lcpm)$y)+.2), las=2, main="", xlab="")
title(main="B. Filtered genes", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=batch.col[y$samples$batch][i], lwd=2)
}
legend("topright", legend=c("First Batch","Second Batch", "Third Batch"), ncol=1, cex=0.8, text.col=batch.col[unique(y$samples$batch)], bty="n")
dev.off()

# replicate analysis after removal of lowl-expressed genes ----------------------------------------------------------------

# get correlation of technical replicates after filtering for lowly expressed genes
pdf("replicate_comparisons_postFiltering_123Combined.pdf", height=10, width=10)
par(mfrow=c(3,3))
# Since SMB-ANK-027 had two replicates, we need to plot one vs three and two vs three. First, one vs three
smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[1]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]))$r.squared, digits=2)))
# best fit/regression line
abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[1]]), col="green")
# diagonal line
abline(a=0,b=1,col="red")

# Now, replicates two vs three
smoothScatter(lcpm[,which(samplenames %in% "SMB-ANK-027")[2]], lcpm[,which(samplenames %in% "SMB-ANK-027")[3]], ylab=colnames(lcpm)[which(samplenames %in% allreps[1])[3]], xlab=colnames(lcpm)[which(samplenames %in% allreps[1])[2]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate SMB-ANK-027", "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]))$r.squared, digits=2)))
# best fit/regression line
abline(lm(lcpm[,which(samplenames %in% "SMB-ANK-027")[3]]~lcpm[,which(samplenames %in% "SMB-ANK-027")[2]]), col="green")
# diagonal line
abline(a=0,b=1,col="red")

for (i in 1:length(allreps)){
  smoothScatter(lcpm[,which(samplenames %in% allreps[i])[1]], lcpm[,which(samplenames %in% allreps[i])[2]], ylab=colnames(lcpm)[which(samplenames %in% allreps[i])[2]], xlab=colnames(lcpm)[which(samplenames %in% allreps[i])[1]], xlim=c(-5,15), ylim=c(-5,15), main=paste("Technical replicate ", allreps[i], "\n", "r2 =",round(summary(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]))$r.squared, digits=2)))
  # best fit/regression line
  abline(lm(lcpm[,which(samplenames %in% allreps[i])[2]]~lcpm[,which(samplenames %in% allreps[i])[1]]), col="green")
  # diagonal line
  abline(a=0,b=1,col="red")
}
dev.off()

# normalise gene-expression distribution -------------------------------------------------

# eliminate composition bias between libraries by upper quartile normalisation (this will be used in the 'data exploration' step)
y <- calcNormFactors(y, method="TMM")

# look at performance of normalisation
# First, compare all three methods
pdf("NormalisedGeneExpressionDistribution_IndoRNA_all3Methods.pdf", height=15, width=15)
par(oma=c(2,0,0,0), mfrow=c(4,1))
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
title(main="Unnormalised data",ylab="Log-cpm")
for (method in c("upperquartile", "TMM", "RLE")) {
  y2 <- calcNormFactors(y2, method="upperquartile")
  lcpm <- cpm(y2, log=TRUE)
  boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
  title(main=paste0("Normalised data ", method),ylab="Log-cpm")
}
dev.off()

# now just plot how well TMM normalisation worked (since wqe like that one the best)
pdf("NormalisedGeneExpressionDistribution_IndoRNA_TMM.pdf", height=15, width=15)
par(oma=c(2,0,0,0), mfrow=c(2,1))
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
title(main="A. Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2, method="TMM")
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=batch.col[as.numeric(batch)], main="", cex.axis=0.75, names=samplenames)
title(main="B. Normalised data, TMM",ylab="Log-cpm")
dev.off()

# recalculate lcpm
lcpm=cpm(y, log=T)



