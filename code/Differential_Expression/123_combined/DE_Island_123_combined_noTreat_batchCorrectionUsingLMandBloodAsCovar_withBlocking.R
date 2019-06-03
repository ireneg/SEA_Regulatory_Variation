# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: IGR 2019.06.03 
### Changed references to Mappi into Korowai, kept West_Papua in there for downstream uses and better nomenclature.

### 0. Load dependencies and functions and set input paths -------------------------- ###
### 1. Begin analyses and initial QC ---------------------------------------------------------------------------------- ###
### 2. DE testing with duplicate correlation and blocking ----------------------------------------------------- ###
### 3. DE testing without duplicate correlation ------------------------------------------------- ###
### 4. Visual QC of duplicate correlation voom output after fitting linear models ---------------- ###
### 5. Summary and visualisation of gene trends --------------------------------------------- ###
### 6. Looking at the top ranked genes ------------------------------------------------- ###

### TO DO:
### Fix everything that's commented out (just figures)


#########################################################################################
### 0. Load dependencies and functions and set input paths -------------------------- ###
#########################################################################################

# Load dependencies:
library(edgeR)
library(plyr)
library(NineteenEightyR)
library(RColorBrewer)
# library(biomaRt)
# library(ggpubr)
library(ggplot2)
library(ggsignif)
# library(pheatmap)
library(viridis)
# library(gplots)
library(circlize)
library(ComplexHeatmap)
library(VennDiagram)
library(wesanderson)

# Set paths:
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/" # on server
covariatedir <- "/data/cephfs/punim0586/igallego/indoRNA/"

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/"
edaoutput <- paste0(outputdir, "/eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

# Load colour schemes:
#Colour schemes:
korowai <- wes_palette("Zissou1", 20, type = "continuous")[20]
mentawai <- wes_palette("Zissou1", 20, type = "continuous")[1]
sumba <- wes_palette("Zissou1", 20, type = "continuous")[11]

smb_mtw <- wes_palette("Darjeeling1", 9, type = "continuous")[3]
smb_kor <- wes_palette("Darjeeling1", 9, type = "continuous")[7]
mtw_kor <- "darkorchid4"
# set up colour palette for batch
batch.col=electronic_night(n=3)

# Venn diagram plotting figure:
make.venn.triple <- function(geneset1, geneset2, geneset3, prefix, geneset1.label, geneset2.label, geneset3.label, universe){
    universe$g1 <- universe$genes %in% geneset1
    universe$g2 <- universe$genes %in% geneset2
    universe$g3 <- universe$genes %in% geneset3
    pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
    venn.placeholder <- draw.triple.venn(length(geneset1), length(geneset2), length(geneset3), dim(universe[universe$g1 == T & universe$g2 == T,])[1], dim(universe[universe$g2 == T & universe$g3 == T,])[1], dim(universe[universe$g1 == T & universe$g3 == T,])[1], dim(universe[universe$g1 == T & universe$g2 == T & universe$g3 == T,])[1], c(geneset1.label, geneset2.label, geneset3.label), fill=c("goldenrod", "plum4", "steelblue3"), alpha=c(0.5, 0.5, 0.5),col=NA, euler.d=T)
    complement.size <- dim(universe[universe$g1 == F & universe$g2 == F & universe$g3 == F,][1])
    grid.text(paste(complement.size, " not in\nany", sep=""), x=0.1, y=0.1)
    dev.off()
    print(paste("Genes in a: ", length(geneset1), sep=""))
    print(paste("Genes in b: ", length(geneset2), sep=""))
    print(paste("Genes in c: ", length(geneset3), sep=""))
    print(paste("Common genes: ", dim(universe[universe$g1 == T & universe$g2 == T & universe$g3 == T,])[1], sep=""))
}

# Load log CPM matrix and y object:
# lcpm
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
# y DGE list object
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))


###########################################################################################################################
### 1. Begin analyses and initial QC ---------------------------------------------------------------------------------- ###
###########################################################################################################################

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
yFilt$samples$Age[which(is.na(yFilt$samples$Age) == T)]=45

# Set up design matrix
design <- model.matrix(~0 + yFilt$samples$Island + yFilt$samples$Age + yFilt$samples$batch + yFilt$samples$RIN + yFilt$samples$CD8T + yFilt$samples$CD4T + yFilt$samples$NK + yFilt$samples$Bcell + yFilt$samples$Mono + yFilt$samples$Gran)
colnames(design)=gsub("Island", "", colnames(design))

# rename columns to exclude spaces and unrecognised characters
colnames(design)=gsub("yFilt\\$samples\\$", "", colnames(design))
colnames(design)=gsub("West Papua", "West_Papua", colnames(design))

# Rename Mappi to Korowai for downstream processing:
yFilt$samples$Sampling.Site <- gsub("Mappi", "Korowai", yFilt$samples$Sampling.Site)

# set up contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsKOR=Sumba - West_Papua, MTWvsKOR=Mentawai - West_Papua, levels=colnames(design))

yFilt <- calcNormFactors(yFilt, method="TMM")

# We can view how TMM normalisation performed using MD plots. This visualizes the library size-adjusted log-fold change between
# two libraries (the difference) against the average log-expression across those libraries (themean). MD plots are generated by comparing sample 1 against an artificial
# library constructed from the average of all other samples. Ideally, the bulk of genes should be centred at a log-fold change of zero.  This indicates
# that any composition bias between libraries has been successfully removed
pdf(paste0(outputdir,"MDPlots_TMM_Normalisation_OutliersCheck.pdf"), height=15, width=15)
par(mfrow=c(4,4))
for (i in 1:ncol(yFilt)){
  plotMD(cpm(yFilt, log=TRUE), column=i)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

###################################################################################################################
### 2. DE testing with duplicate correlation and blocking ----------------------------------------------------- ###
###################################################################################################################

# create a new variable for blocking using sample IDs
yFilt$samples$ind <- sapply(strsplit(as.character(yFilt$samples$samples), "[_.]"), `[`, 1)

# First, we need to perform voom normalisation

# No normalisation between samples beyond tmm and voom:
    voomNoNorm <- voom(yFilt, design, normalize.method="none", plot=F) 
    dupcorNone <- duplicateCorrelation(voomNoNorm, design, block=yFilt$samples$ind) # 23 or more non-convergences
    # The value dupcor$consensus estimates the average correlation within the blocks and should be positive
    dupcorNone$consensus # sanity check
    # [1] 0.6796663
    median(voomNoNorm$weights) # another sanity check:
    # [1] 22.8338
    save(voomNoNorm, file=paste0(outputdir, "voomNoNorm.tmm.filtered.indoRNA.Rda"))

    # Second round:
    voomNoNormDup <- voom(yFilt, design, plot=TRUE, block=yFilt$samples$ind, correlation=dupcorNone$consensus)
    dupcorNoneDup <- duplicateCorrelation(voomNoNormDup, design, block=yFilt$samples$ind) # 24 non convergences
    dupcorNoneDup$consensus # sanity check pt 2
    # [1] 0.6796721
    median(voomNoNormDup$weights) # another sanity check, pt 2 
    # [1] 22.41583

    pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.densities.pdf"))
        plotDensities(voomNoNormDup, group=yFilt$samples$batch)
        plotDensities(voomNoNormDup, group=yFilt$samples$Island)
    dev.off()
    save(voomNoNormDup, file=paste0(outputdir, "voomNoNorm.tmm.filtered.duplicate_corrected.indoRNA.Rda"))

    # DE testing:
    # the inter-subject correlation is input into the linear model fit
    voomNoNormDupVfit <- lmFit(voomNoNormDup, design, block=yFilt$samples$ind, correlation=dupcorNoneDup$consensus)
    voomNoNormDupVfit <- contrasts.fit(voomNoNormDupVfit, contrasts=contr.matrix)
    voomNoNormDupEfit <- eBayes(voomNoNormDupVfit, robust=T)

    pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.mean-variance-trend.pdf"))
        plotSA(voomNoNormDupEfit, main="Mean-variance trend elimination with duplicate correction")
    dev.off()

    # get top genes using toptable
    voomNoNormDupTopTableSMB.MTW <- topTable(voomNoNormDupEfit, coef=1, n=Inf, sort.by="p")
    voomNoNormDupTopTableSMB.KOR <- topTable(voomNoNormDupEfit, coef=2, n=Inf, sort.by="p")
    voomNoNormDupTopTableMTW.KOR <- topTable(voomNoNormDupEfit, coef=3, n=Inf, sort.by="p")

# Quantile normalisation between samples and tmm and voom:
    voomQuant <- voom(yFilt, design, normalize.method="quantile", plot=F) 
    dupcorQuant <- duplicateCorrelation(voomQuant, design, block=yFilt$samples$ind) # 46 warnings
    # The value dupcor$consensus estimates the average correlation within the blocks and should be positive
    dupcorQuant$consensus # sanity check
    # [1] 0.6888087
    median(voomQuant$weights) # another sanity check:
    # [1] 22.94226
    save(voomQuant, file=paste0(outputdir, "voomQuant.tmm.filtered.indoRNA.Rda"))

    # Second round:
    voomQuantDup <- voom(yFilt, design, plot=TRUE, block=yFilt$samples$ind, correlation=dupcorQuant$consensus)
    dupcorQuantDup <- duplicateCorrelation(voomQuantDup, design, block=yFilt$samples$ind) # 48 warnings
    dupcorQuantDup$consensus # sanity check pt 2
    # [1] 0.6796766
    median(voomQuantDup$weights) # another sanity check, pt 2 
    # [1] 22.37909

    pdf(file=paste0(edaoutput, "voomQuant.tmm.filtered.indoRNA.densities.pdf"))
        plotDensities(voomQuantDup, group=yFilt$samples$batch)
        plotDensities(voomQuantDup, group=yFilt$samples$Island)
    dev.off()
    save(voomQuantDup, file=paste0(outputdir, "voomQuant.tmm.filtered.duplicate_corrected.indoRNA.Rda"))

    # DE testing:
    # the inter-subject correlation is input into the linear model fit
    voomQuantDupVfit <- lmFit(voomQuantDup, design, block=yFilt$samples$ind, correlation=dupcorQuantDup$consensus)
    voomQuantDupVfit <- contrasts.fit(voomQuantDupVfit, contrasts=contr.matrix)
    voomQuantDupEfit <- eBayes(voomQuantDupVfit, robust=T)

    pdf(file=paste0(edaoutput, "voomQuant.tmm.filtered.indoRNA.mean-variance-trend.pdf"))
        plotSA(voomQuantDupEfit, main="Mean-variance trend elimination with duplicate correction")
    dev.off()

    # get top genes using toptable
    voomQuantDupTopTableSMB.MTW <- topTable(voomQuantDupEfit, coef=1, n=Inf, sort.by="p")
    voomQuantDupTopTableSMB.KOR <- topTable(voomQuantDupEfit, coef=2, n=Inf, sort.by="p")
    voomQuantDupTopTableMTW.KOR <- topTable(voomQuantDupEfit, coef=3, n=Inf, sort.by="p")


# Loess normalisation between samples and tmm and voom:
    voomLoess <- voom(yFilt, design, normalize.method="cyclicloess", plot=F) 
    dupcorLoess <- duplicateCorrelation(voomLoess, design, block=yFilt$samples$ind) # 46 warnings
    # The value dupcor$consensus estimates the average correlation within the blocks and should be positive
    dupcorLoess$consensus # sanity check
    # [1] 0.6825221
    median(voomLoess$weights) # another sanity check:
    # [1] 23.18414
    save(voomLoess, file=paste0(outputdir, "voomLoess.tmm.filtered.indoRNA.Rda"))

    # Second round:
    voomLoessDup <- voom(yFilt, design, plot=TRUE, block=yFilt$samples$ind, correlation=dupcorLoess$consensus)
    dupcorLoessDup <- duplicateCorrelation(voomLoessDup, design, block=yFilt$samples$ind) # 48 warnings
    dupcorLoessDup$consensus # sanity check pt 2
    # [1] 0.6796738
    median(voomLoessDup$weights) # another sanity check, pt 2 
    # [1] 22.40257

    pdf(file=paste0(edaoutput, "voomLoess.tmm.filtered.indoRNA.densities.pdf"))
        plotDensities(voomLoessDup, group=yFilt$samples$batch)
        plotDensities(voomLoessDup, group=yFilt$samples$Island)
    dev.off()
    save(voomLoessDup, file=paste0(outputdir, "voomLoess.tmm.filtered.duplicate_corrected.indoRNA.Rda"))

    # DE testing:
    # the inter-subject correlation is input into the linear model fit
    voomLoessDupVfit <- lmFit(voomLoessDup, design, block=yFilt$samples$ind, correlation=dupcorLoessDup$consensus)
    voomLoessDupVfit <- contrasts.fit(voomLoessDupVfit, contrasts=contr.matrix)
    voomLoessDupEfit <- eBayes(voomLoessDupVfit, robust=T)

    pdf(file=paste0(edaoutput, "voomLoess.tmm.filtered.indoRNA.mean-variance-trend.pdf"))
        plotSA(voomLoessDupEfit, main="Mean-variance trend elimination with duplicate correction")
    dev.off()

    # get top genes using toptable
    voomLoessDupTopTableSMB.MTW <- topTable(voomLoessDupEfit, coef=1, n=Inf, sort.by="p")
    voomLoessDupTopTableSMB.KOR <- topTable(voomLoessDupEfit, coef=2, n=Inf, sort.by="p")
    voomLoessDupTopTableMTW.KOR <- topTable(voomLoessDupEfit, coef=3, n=Inf, sort.by="p")

# LOL omg what was the point? On the basis of this, going with no norm so there's no need to think about it any deeper.
summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down        898     2325     2102
# NotSig    11430     8479     8887
# Up          647     2171     1986
summary(decideTests(voomQuantDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down        896     2324     2102
# NotSig    11432     8480     8887
# Up          647     2171     1986
summary(decideTests(voomLoessDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down        898     2325     2102
# NotSig    11430     8479     8887
# Up          647     2171     1986

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down         96      606      536
# NotSig    12661    11577    11958
# Up          218      792      481
summary(decideTests(voomQuantDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down         96      606      536
# NotSig    12661    11577    11958
# Up          218      792      481
summary(decideTests(voomLoessDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down         96      606      536
# NotSig    12661    11577    11958
# Up          218      792      481

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down          6       87       96
# NotSig    12940    12662    12748
# Up           29      226      131

summary(decideTests(voomQuantDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down          6       87       96
# NotSig    12940    12662    12748
# Up           29      226      131

summary(decideTests(voomLoessDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#        SMBvsMTW SMBvsKOR MTWvsKOR
# Down          6       87       96
# NotSig    12940    12662    12748
# Up           29      226      131

write.table(voomNoNormDupTopTableSMB.MTW, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.SMB-MTW.txt"))
write.table(voomNoNormDupTopTableSMB.KOR, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.SMB-KOR.txt"))
write.table(voomNoNormDupTopTableMTW.KOR, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.MTW-KOR.txt"))

# For easily combining with the villages:
save(voomNoNormDupEfit, file=paste0(outputdir, "voomNoNorm.tmm.filtered.dup_corrected.Efit_object.Rda"))

#####################################################################################################
### 3. DE testing without duplicate correlation ------------------------------------------------- ###
#####################################################################################################

### Only doing the nonorm one, because why bother with the others?

    voomNoNormVfit <- lmFit(voomNoNorm, design)
    voomNoNormVfit <- contrasts.fit(voomNoNormVfit, contrasts=contr.matrix)
    voomNoNormEfit <- eBayes(voomNoNormVfit, robust=T)

    # voomQuantVfit <- lmFit(voomQuant, design)
    # voomQuantVfit <- contrasts.fit(voomQuantVfit, contrasts=contr.matrix)
    # voomQuantEfit <- eBayes(voomQuantVfit, robust=T)

    # voomLoessVfit <- lmFit(voomLoess, design)
    # voomLoessVfit <- contrasts.fit(voomLoessVfit, contrasts=contr.matrix)
    # voomLoessEfit <- eBayes(voomLoessVfit, robust=T)

pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.no_dup_correction.mean-variance-trend.pdf"))
    plotSA(voomNoNormEfit, main="Mean-variance trend elimination without duplicate correction")
dev.off()

# get top genes using toptable
topTableSMB.MTW <- topTable(voomNoNormEfit, coef=1, n=Inf, sort.by="p")
topTableSMB.KOR <- topTable(voomNoNormEfit, coef=2, n=Inf, sort.by="p")
topTableMTW.KOR <- topTable(voomNoNormEfit, coef=3, n=Inf, sort.by="p")

# no LFC threshold
summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#       SMBvsMTW SMBvsKOR MTWvsKOR
#Down       1032     2569     2228
#NotSig    11183     8106     8669
#Up          760     2300     2078

# LFC of 0.05
summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#       SMBvsMTW SMBvsKOR MTWvsKOR
#Down        112      654      558
#NotSig    12603    11492    11920
#Up          260      829      497

# LFC of 1
summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#       SMBvsMTW SMBvsKOR MTWvsKOR
#Down          6       93      103
#NotSig    12923    12643    12737
#Up           46      239      135

# Let's check the correlation between those two approaches - sort by gene first, then cor test on adjusted p-value
MTW.KOR <- join(voomNoNormDupTopTableMTW.KOR, topTableMTW.KOR, by="genes")
cor(MTW.KOR[,6], MTW.KOR[,12], method="spearman"    )
# [1] 0.9884076

SMB.KOR <- join(voomNoNormDupTopTableSMB.KOR, topTableSMB.KOR, by="genes")
cor(SMB.KOR[,6], SMB.KOR[,12], method="spearman", use="complete")
# [1] 0.9771854

SMB.MTW <- join(voomNoNormDupTopTableSMB.MTW, topTableSMB.MTW, by="genes")
cor(SMB.MTW[,6], SMB.MTW[,12], method="spearman", use="complete")
# [1] 0.955927

write.table(topTableSMB.MTW, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.not_dup_corrected.SMB-MTW.txt"))
write.table(topTableSMB.KOR, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.not_dup_corrected.SMB-KOR.txt"))
write.table(topTableMTW.KOR, file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.not_dup_corrected.MTW-KOR.txt"))


######################################################################################################
### 4. Visual QC of duplicate correlation voom output after fitting linear models ---------------- ###
######################################################################################################

# check to see p-value distribution is normal
pdf(paste0(edaoutput,"PvalueDist_NotAdjusted_dupCor.pdf"), height=15, width=10)
    par(mfrow=c(3,1))
    for (i in 1:ncol(voomNoNormDupEfit)){
        hist(voomNoNormDupEfit$p.value[,i], main=colnames(voomNoNormDupEfit)[i], ylim=c(0,max(table(round(voomNoNormDupEfit$p.value[,i], 1)))+1000), xlab="p-value")
    }
dev.off()

# check p-value distribution for adjusted p-values
pdf(paste0(edaoutput,"PvalueDist_Adjusted_dupCor.pdf"), height=15, width=10)
    par(mfrow=c(3,1))
    for (i in 1:ncol(voomNoNormDupEfit)){
        topTable <- topTable(voomNoNormDupEfit, coef=i, n=Inf)
        histData <- hist(topTable$adj.P.Val, main=colnames(voomNoNormDupEfit)[i], xlab="p-value")
        hist(topTable$adj.P.Val, main=colnames(voomNoNormDupEfit)[i], ylim=c(0,max(histData$counts)+1000), xlab="p-value")
    }
dev.off()

##########################################################################################################
### IGR NOTE 2019.04.12 - NEEDS FIXING WITH NEW GENE LIST, CAN THIS EVEN RUN ON SPARTAN COMPUTE NODES? ###
##########################################################################################################

                                            # # Verify that control housekeeping genes are not significantly DE. Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
                                            # housekeeping <- read.table(paste0(housekeepingdir,"Housekeeping_ControlGenes.txt"), as.is=T, header=F) 
                                            # # if this is broken, use host = "uswest.ensembl.org"
                                            # ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
                                            # biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
                                            # hkGenes <- as.vector(biomart.results.table[,1])
                                            # hkControls <- hkGenes[which(hkGenes %in% rownames(yFilt$counts))]

                                            # # Volcano plot with points of housekeeping genes
                                            # pdf(paste0(outputdir,"VolcanoPlots_dupCorEfit.pdf"), height=15, width=10)
                                            #     par(mfrow=c(3,1))
                                            #     for (i in 1:ncol(voomDupEfit)){
                                            #         plot(voomDupEfit$coef[,i], -log10(as.matrix(voomDupEfit$p.value)[,i]), pch=20, main=colnames(voomDupEfit)[i], xlab="log2FoldChange", ylab="-log10(pvalue)")
                                            #         points(voomDupEfit$coef[,i][which(names(voomDupEfit$coef[,i]) %in% hkControls)], -log10(as.matrix(voomDupEfit$p.value)[,i][which(names(voomDupEfit$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
                                            #         legend("topleft", "genes", "hk genes",fill=4)
                                            #         abline(v=c(-1,1))
                                            #     }
                                            # dev.off()

#################################################################################################
### 5. Summary and visualisation of gene trends --------------------------------------------- ###
#################################################################################################

# plot log2 fold change between islands
pdf(paste0(edaoutput,"log2FC_IslandComparisons_pval01_dupCor.pdf"))
# note 'p.value' is the cutoff value for adjusted p-values
    topTable <- topTable(voomNoNormDupEfit, coef=1, n=Inf, p.value=0.01)
    plot(density(topTable$logFC), col=9, xlim=c(-2,2), main="LogFC Density", xlab="LogFC", ylab="Density", lwd=3)
    abline(v=c(-1,-0.5,0.5,1), lty=3)
    counter=0
    for (i in 2:ncol(voomNoNormDupEfit)){
        counter=counter+1
        topTable <- topTable(voomNoNormDupEfit, coef=i, n=Inf, p.value=0.01)
        lines(density(topTable$logFC), col=9+counter, xlim=c(-2,2), lwd=3)
    }
    legend(x="topright", bty="n", col=9:11, legend=colnames(voomNoNormDupEfit), lty=1, lwd=2)
dev.off()

# graphical representation of DE results through MD plot
pdf(paste0(edaoutput,"MD_Plots_pval01_lfc1_dupCor.pdf"))
    plotMD(voomNoNormDupEfit, column = 1, array = NULL, xlab = "Average log-expression", ylab = "Expression log-ratio",
       main = colnames(voomNoNormDupEfit)[1], status=voomNoNormDupEfit$genes$Status, zero.weights = FALSE)
    abline(h=c(-1,-0.5,0.5,1), lty=3)
    plotMD(voomNoNormDupEfit, column = 2, array = NULL, xlab = "Average log-expression", ylab = "Expression log-ratio",
       main = colnames(voomNoNormDupEfit)[2], status=voomNoNormDupEfit$genes$Status, zero.weights = FALSE)
    abline(h=c(-1,-0.5,0.5,1), lty=3)
    plotMD(voomNoNormDupEfit, column = 3, array = NULL, xlab = "Average log-expression", ylab = "Expression log-ratio",
       main = colnames(voomNoNormDupEfit)[3], status=voomNoNormDupEfit$genes$Status, zero.weights = FALSE)
    abline(h=c(-1,-0.5,0.5,1), lty=3)
dev.off()

# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering

#######################################################################################
### IGR NOTE 2019.04.12 - THIS TO BE REPLACED WITH A CLEANER CALL TO COMPLEXHEATMAP ###
#######################################################################################


                    # first, make a heatmap of all top genes in one pdf

                    # reset ensemble row names to gene symbols
                    # rownames(vDup$E)=vDup$genes$SYMBOL

                    # # set up annotation
                    # col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

                    # df1=data.frame(island = as.character(Island))
                    # df2=data.frame(batch = as.numeric(batch))
                    # ha1 = HeatmapAnnotation(df = df1, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3)))

                    # pdf(paste0(outputdir,"HeatmapAllPops_dupCor.pdf"), height=15, width=15)
                    # grid.newpage()
                    # pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
                    # # set up layout row position
                    # layout.row=c(1,1,2)
                    # # set up column position
                    # layout.col=c(1,2,1)

                    # for (i in 1:ncol(voomNoNormDupEfit)){
                    #     topTable <- topTable(voomNoNormDupEfit, coef=i, p.value=0.01, lfc=1, n=Inf, sort.by="p")
                    #     index <- which(voomNoNormDupEfit$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
                    #     pushViewport(viewport(layout.pos.row = layout.row[i], layout.pos.col = layout.col[i]))
                    #     draw(Heatmap(t(scale(t(vDup$E[index,]))), col=col_fun, column_title = colnames(voomNoNormDupEfit)[i], top_annotation = ha1, show_row_names = T, show_heatmap_legend = F, show_column_names = F, name = "Z-Score"),show_annotation_legend = FALSE,newpage=F)
                    #     upViewport()

                    # }

                    # lgd = Legend(at = c(-4,0,4), title = "Row Z-Score", col_fun = col_fun, grid_height = unit(1, "cm"), grid_width = unit(10, "mm"))
                    # lgd2 = Legend(at = c(1,2,3), legend_gp = gpar(fill = 1:3), labels=c("Mentawai", "Sumba","West Papua"),title = "Island", grid_height = unit(1, "cm"), grid_width = unit(10, "mm"))
                    # pd = packLegend(lgd, lgd2, direction = "horizontal")

                    # pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
                    # grid.draw(pd)
                    # upViewport()
                    # dev.off()

# We can also make individual pdfs of the top genes
island1 <- c("Sumba","Mentawai","West Papua")
island2 <- c("Sumba","Mentawai","West Papua")

counter <- 0

# for (i1 in island1){
#     island2=island2[-1]
#     for (i2 in island2){
#         counter=counter+1
#         topTable <- topTable(voomNoNormDupEfit, coef=counter, p.value=0.01, lfc=1, n=Inf, sort.by="p")
#         index <- which(vDup$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
#         df=data.frame(island = as.character(Island[grep(paste(i1,i2,sep="|"), Island)]))
#         ha =  HeatmapAnnotation(df = df, col = list(island = c("Mentawai" =  1, "Sumba" = 2, "West Papua" = 3)))
#         pdf(paste0(outputdir,"HeatmapTopeGenes_",i1,"_vs_",i2,"_dupCor.pdf"), height=10, width=15)
#         draw(Heatmap(t(scale(t(vDup$E[index,])))[,grep(paste(i1,i2,sep="|"), Island)], col=col_fun, column_title = colnames(voomNoNormDupEfit)[counter], top_annotation = ha, show_row_names = T, show_heatmap_legend = T, show_column_names = F, name = "Z-Score"),show_annotation_legend = TRUE,newpage=F)
#         dev.off()
#     }
# }



# Some Venn diagrams
make.venn.triple(voomNoNormDupTopTableSMB.MTW[voomNoNormDupTopTableSMB.MTW$adj.P.Val <= 0.01,]$genes, voomNoNormDupTopTableSMB.KOR[voomNoNormDupTopTableSMB.KOR$adj.P.Val <= 0.01,]$genes, voomNoNormDupTopTableMTW.KOR[voomNoNormDupTopTableMTW.KOR$adj.P.Val <= 0.01,]$genes, paste0(edaoutput, "all_islands.fdr_0.01"), "Sumba vs\nMentawai", "Sumba vs\nKorowai", "Mentawai\nvs Korowai", voomNoNormDupTopTableSMB.MTW)

make.venn.triple(voomNoNormDupTopTableSMB.MTW[voomNoNormDupTopTableSMB.MTW$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableSMB.MTW$logFC)>= 0.5,]$genes, voomNoNormDupTopTableSMB.KOR[voomNoNormDupTopTableSMB.KOR$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableSMB.KOR$logFC)>= 0.5,]$genes, voomNoNormDupTopTableMTW.KOR[voomNoNormDupTopTableMTW.KOR$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableMTW.KOR$logFC)>= 0.5,]$genes, paste0(edaoutput, "all_islands.fdr_0.01.logfc_0.5"), "Sumba vs\nMentawai", "Sumba vs\nKorowai", "Mentawai\nvs Korowai", voomNoNormDupTopTableSMB.MTW)

make.venn.triple(voomNoNormDupTopTableSMB.MTW[voomNoNormDupTopTableSMB.MTW$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableSMB.MTW$logFC)>= 1,]$genes, voomNoNormDupTopTableSMB.KOR[voomNoNormDupTopTableSMB.KOR$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableSMB.KOR$logFC)>= 1,]$genes, voomNoNormDupTopTableMTW.KOR[voomNoNormDupTopTableMTW.KOR$adj.P.Val <= 0.01 & abs(voomNoNormDupTopTableMTW.KOR$logFC)>= 1,]$genes, paste0(edaoutput, "all_islands.fdr_0.01.logfc_1"), "Sumba vs\nMentawai", "Sumba vs\nKorowai", "Mentawai\nvs Korowai", voomNoNormDupTopTableSMB.MTW)


    ###################################################################################################
    ### IGR NOTE 2019.04.12 - I BELIEVE HEINI IS NOW MAKING THIS FIGURE, WILL REVISIT IT OTHERWISE. ###
    ###################################################################################################

# get DE genes in common with populations compared to Korowai, i.e., SMBvsKOR and MTWvsKOR (since we think this is an interesting island comparison)
allGenes <- merge(voomNoNormDupTopTableSMB.MTW, voomNoNormDupTopTableSMB.KOR, by.x="genes", by.y="genes", suffixes=c(".SMB.MTW", ".SMB.KOR"))
allGenes <- merge(allGenes, voomNoNormDupTopTableMTW.KOR, by.x="genes", by.y="genes")
names(allGenes)[13:19] <- paste0(names(allGenes)[13:19], ".MTW.KOR")

deSummaryAll <- decideTests(voomNoNormDupEfit, p.value=0.01)
deSummary05 <- decideTests(voomNoNormDupEfit, p.value=0.01, lfc=0.5)
deSummary1 <- decideTests(voomNoNormDupEfit, p.value=0.01, lfc=1)

deCommonKOR = which(deSummaryAll[,2]!=0 & deSummaryAll[,3]!=0)
deCommonKOR05 = which(deSummary05[,2]!=0 & deSummary05[,3]!=0)
deCommonKOR1 = which(deSummary1[,2]!=0 & deSummary1[,3]!=0)

# get what these genes are doing and save them to a file
# commonGenes.KOR <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common.KOR), filters="ensembl_gene_id")
# write.table(de.common.KOR, file=paste0(outputdir,"allCommonGenes_KOR_dupcor.txt"))
# # save the common gene names 
# de.common.KOR=voomNoNormDupEfit$genes[names(de.common.KOR),]

# now plot the common genes to see if they're being regulated in the same direction
pdf(paste0(edaoutput,"logFC_commonKORgenes_dupCor.pdf"))
    plot(voomNoNormDupTopTableSMB.KOR[rownames(deCommonKOR), "logFC"], voomNoNormDupTopTableMTW.KOR[rownames(deCommonKOR), "logFC"], xlab="logFC SMBvsKOR", ylab="logFC MTWvsKOR", pch=20, main="Common DE Genes", xlim=c(-5,5), ylim=c(-6,6))
    # text(tt.SMBvsKOR[rownames(de.common.KOR),"logFC"], tt.MTWvsKOR[rownames(de.common.KOR),"logFC"], labels=tt.SMBvsKOR[rownames(de.common.KOR),"SYMBOL"], pos=3)
    abline(h=0,v=0, lty=2)
dev.off()


############################################################################################
### 6. Looking at the top ranked genes ------------------------------------------------- ###
############################################################################################

# # Let's see how the expression levels of all of the significantly DE genes in population comparisons with Korowai are distributed within each island. First, assign our top genes and ensembl IDs to variables
# topGenes=de.common.KOR[,2]
# topEnsembl=de.common.KOR[,1]

# # To visualise distributions, we'll be making violin plots using ggpubr which needs p-value labels. Let's go ahead and make a matrix to input this into ggpubr
# # first set up matrix
# topGenes.pvalue=matrix(nrow=length(topEnsembl), ncol=ncol(voomNoNormDupEfit))
# rownames(topGenes.pvalue)=topEnsembl
# colnames(topGenes.pvalue)=colnames(voomNoNormDupEfit)
# for (i in 1:ncol(voomNoNormDupEfit)){
#     # get significant genes over a logFC of 1 for all Island comparisons
#     topTable <- topTable(voomNoNormDupEfit, coef=i, n=Inf)
#     for(j in topEnsembl){
#         # input the adjusted p.value for each gene
#         topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]
#     }
# }

# # make pvalues into scientific notation with max 3 digits
# topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# # convert e notation to base 10 notation
# topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

# # We can make the violin plots using ggpubr
# pdf(paste0(outputdir,"TopGenes_ggboxplot_Island.pdf"), height=8, width=10)
# counter=0
# for(ensembl in topEnsembl){
#     counter=counter+1
#     gene.df <- data.frame(vDup$E[which(vDup$genes$ENSEMBL==ensembl),],Island)
#     colnames(gene.df)=c("CPM", "Island")
#     annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
#     print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=topGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
# }
# dev.off()

# # after analysing the distributions and reading up on some of the genes, my three favourite genes are Siglec6, Siglec7, and MARCO. Lets plot out the distribution solely for these three genes
# favGenes=c("SIGLEC6","SIGLEC7","MARCO")
# #"RSAD2","AIM2","TNFSF4")
# favEnsembl=de.common.KOR[,1][sapply(1:length(favGenes), function(x) grep(favGenes[x],de.common.KOR[,2]))]

# # set up pvalue matrix
# topGenes.pvalue=matrix(nrow=length(favEnsembl), ncol=ncol(voomNoNormDupEfit))
# rownames(topGenes.pvalue)=favEnsembl
# colnames(topGenes.pvalue)=colnames(voomNoNormDupEfit)
# for (i in 1:ncol(voomNoNormDupEfit)){
#     # get significant genes over a logFC of 1 for all Island comparisons
#     topTable <- topTable(voomNoNormDupEfit, coef=i, n=Inf)
#     for(j in favEnsembl){
#         # input the adjusted p.value for each gene
#         topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]
#     }
# }

# # make pvalues into scientific notation with max 3 digits
# topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# # convert e notation to base 10 notation
# topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

# # We can make the violin plots using ggpubr
# counter=0
# for(ensembl in favEnsembl){
#     counter=counter+1
#     # pdf(paste0("FavouriteGenes_ggboxplot_",favGenes[counter],".pdf"), height=8, width=10)
#     gene.df <- data.frame(vDup$E[which(vDup$genes$ENSEMBL==ensembl),],Island)
#     colnames(gene.df)=c("CPM", "Island")
#     annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
#     # print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=favGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
#     assign(favGenes[counter], ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=favGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 3, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
# }

# pdf(paste0(outputdir,"favouriteTopGenes_distribution_Island.pdf"), height=12, width=15)
# ggarrange(SIGLEC6,SIGLEC7,MARCO)
# #AIM2,TNFSF4,RSAD2)
# dev.off()

# # finally, get logFC thresholds
# logFC.df=matrix(nrow=3,ncol=3)
# colnames(logFC.df)=colnames(voomNoNormDupEfit)
# counter=0
# for (number in c(0,0.5,1)){
#     counter=counter+1
#     dt <- decideTests(voomNoNormDupEfit, p.value=0.01, lfc=number)
#     values=c(sum(abs(dt[,1])), sum(abs(dt[,2])), sum(abs(dt[,3])))
#     logFC.df[counter,]=values
# }
# logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
# write.table(logFC.df, file=paste0(outputdir,"logFC_thresholds.txt"))



