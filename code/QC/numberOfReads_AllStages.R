# Script for creating barplot of number of reads at all different processing stages
# Code developed by Katalina Bobowik, 10.12.2018

# load library and colour palette
library(viridis)
library(Rcmdr)

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Paper_Figures")

# read in summary file and tidy up 
a=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Thesis_Writing/Tables/nReads_AllSamples.txt", header=T, sep="\t")
a=a[1:123,1:5]
rownames(a)=make.unique(as.character(a[,1]))
a$Sample = NULL

# plot
pdf("TotalandPercentageOfReads_allFilteringStages.pdf", height=20,width=17)
par(mar=c(6.1,4.1,10.1,2.1), xpd=T, mfrow=c(2,1))
# first plot total number of reads
barplot(as.matrix(t(a)*1e-6), col=viridis(4), cex.names=0.75, las=3, ylim=c(0,200), main="Number of Reads \nAll filtering stages", ylab="Number of Reads (millions)")
legend(135,210, col=viridis(4), legend=colnames(a), pch=15, cex=0.9)
# now plot reads as a percentage
#barplot(colPercents(t(a))[1:4,]/100, col=viridis(4), cex.names=0.75, las=3, main="Percentage of Reads \nAll filtering stages", ylab="Percentage of Reads")
barplot(t(apply(t(a), 1, "/", a[,1])*100), las=3, ylim=c(0,400),col=viridis(4), cex.names=0.75,main="Percentage of Reads \nAll filtering stages", ylab="Percentage of Reads")
legend(135,410, col=viridis(4), legend=colnames(a), pch=15, cex=0.9)
dev.off()

# pdf("TotalReads_allFilteringStages.pdf", height=10,width=15)
#par(mar=c(6.1,4.1,10.1,2.1), xpd=T)
#barplot(as.matrix(t(a)*1e-6), col=viridis(4), cex.names=0.75, las=3, ylim=c(0,200), main="Number of Reads \n All filtering stages")
#legend(135,210, col=viridis(4), legend=colnames(a), pch=15, cex=0.9)
#dev.off()
