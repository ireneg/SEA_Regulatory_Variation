# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: IGR 2019.10.19
### Changed paths to deal with removal of MPI-296

### 0. Load dependencies and functions and set input paths ###
### 1. Begin analyses and initial QC ###
### 2. DE testing with duplicate correlation and blocking ###
### 3. DE testing without duplicate correlation ###
### 4. Visual QC of duplicate correlation voom output after fitting linear models ###
### 5. Summary and visualisation of gene trends  ###
### 6. UpsetR plots and merging island and village. ###
### 7. Looking at the top ranked genes ###
### 8. Heatmaps of top DE genes between some villages but not others. ###
### 9. Quick check of variance by village, to see what drives the weird distribution of DE results. ###
### 10. Good old plot of pairwise correlations within each village and level etc etc... ###

### TO DO:
### Fix everything that's commented out (just figures)

##############################################################
### 0. Load dependencies and functions and set input paths ###
##############################################################

# Load dependencies:
library(edgeR)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(ggsignif)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(VennDiagram)
library(UpSetR)
library(matrixStats)
library(reshape)
library(wesanderson)
library(dplyr)
library(ggpubr)
#library(GenomicRanges)


# Set paths
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296/" # on server
covariatedir <- "/data/cephfs/punim0586/igallego/indoRNA/"
# inputdir <- "~/Desktop/indoRNA_temp/de_testing/"
# covariatedir <- "~/Desktop/indoRNA_temp/"

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296/"
edaoutput <- paste0(outputdir, "eda/")
# outputdir <- "~/Desktop/indoRNA_temp/de_testing/"
# edaoutput <- paste0(outputdir, "local_eda/")

if (file.exists(edaoutput) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

#Colour schemes:
korowai <- wes_palette("Zissou1", 20, type = "continuous")[20]
mentawai <- wes_palette("Zissou1", 20, type = "continuous")[1]
sumba <- wes_palette("Zissou1", 20, type = "continuous")[11]

smb_mtw <- wes_palette("Darjeeling1", 9, type = "continuous")[3]
smb_kor <- wes_palette("Darjeeling1", 9, type = "continuous")[7]
mtw_kor <- "darkorchid4"

# Load log CPM matrix and y object:
# lcpm
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
# y DGE list object
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))


########################################
### 1. Begin analyses and initial QC ###
########################################

# Rename Mappi to Korowai for downstream processing:
yFilt$samples$Sampling.Site <- gsub("Mappi", "Korowai", yFilt$samples$Sampling.Site)

# First, remove samples that have less than ten individuals per village
table(yFilt$samples$Sampling.Site)
# Anakalung    Bilarenge    Hupu Mada      Korowai      Madobag  Padira Tana Patiala Bawa        Rindi     Taileleu        Wunga   Wura Homba 
#       20            1            5           21           17            3            1            5           32           17            1 

# remove Bilarenge, Hupu Mada, Padira Tana, Patiala Bawa, Rindi, and Wura Homba
yVillage <- yFilt[,-grep("Bilarenge|Hupu Mada|Padira Tana|Patiala Bawa|Rindi|Wura Homba", yFilt$samples$Sampling.Site)]
# drop unused levels
yVillage$samples <- droplevels(yVillage$samples)

# Set up design matrix
design <- model.matrix(~0 + yVillage$samples$Sampling.Site + yVillage$samples$Age + yVillage$samples$batch + yVillage$samples$RIN + yVillage$samples$CD8T + yVillage$samples$CD4T + yVillage$samples$NK + yVillage$samples$Bcell + yVillage$samples$Mono + yVillage$samples$Gran)
# rename columns to exclude spaces and unrecognised characters
colnames(design)=gsub("yVillage\\$samples\\$", "", colnames(design))
colnames(design)=gsub("Sampling.Site", "", colnames(design))

# set up contrast matrix
contr.matrix <- makeContrasts(ANKvsMDB=Anakalung-Madobag, ANKvsKOR=Anakalung-Korowai, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, MDBvsKOR=Madobag-Korowai, MDBvsTLL=Madobag-Taileleu, WNGvsMDB=Wunga-Madobag, TLLvsKOR=Taileleu-Korowai, WNGvsKOR=Wunga-Korowai, WNGvsTLL=Wunga-Taileleu, levels=colnames(design)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

yVillage <- calcNormFactors(yVillage, method="TMM")

#############################################################
### 2. DE testing with duplicate correlation and blocking ###
#############################################################

# create a new variable for blocking using sample IDs
yVillage$samples$ind <- sapply(strsplit(as.character(yVillage$samples$samples), "[_.]"), `[`, 1)

# First, we need to perform voom normalisation

# No normalisation between samples beyond tmm and voom:
    voomNoNorm <- voom(yVillage, design, normalize.method="none", plot=F) 
    dupcorNone <- duplicateCorrelation(voomNoNorm, design, block=yVillage$samples$ind) # 24 non-convergences
    # The value dupcor$consensus estimates the average correlation within the blocks and should be positive
    dupcorNone$consensus # sanity check
    # [1] 0.6835068
    median(voomNoNorm$weights) # another sanity check:
    # [1] 23.90951
    save(voomNoNorm, file=paste0(outputdir, "voomNoNorm.tmm.filtered.indoRNA.village.Rda"))

    # Second round:
    voomNoNormDup <- voom(yVillage, design, plot=TRUE, block=yVillage$samples$ind, correlation=dupcorNone$consensus)
    dupcorNoneDup <- duplicateCorrelation(voomNoNormDup, design, block=yVillage$samples$ind) # 25 non convergences
    dupcorNoneDup$consensus # sanity check pt 2
    # [1] 0.6836133
    median(voomNoNormDup$weights) # another sanity check, pt 2 
    # [1] 23.35687

    pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.densities.villages.pdf"))
        plotDensities(voomNoNormDup, group=yVillage$samples$batch)
        plotDensities(voomNoNormDup, group=yVillage$samples$Island)
    dev.off()
    save(voomNoNormDup, file=paste0(outputdir, "voomNoNorm.tmm.filtered.duplicate_corrected.indoRNA.village.Rda"))

    # DE testing:
    # the inter-subject correlation is input into the linear model fit
    voomNoNormDupVfit <- lmFit(voomNoNormDup, design, block=yVillage$samples$ind, correlation=dupcorNoneDup$consensus)
    voomNoNormDupVfit <- contrasts.fit(voomNoNormDupVfit, contrasts=contr.matrix)
    voomNoNormDupEfit <- eBayes(voomNoNormDupVfit, robust=T)

    pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.mean-variance-trend.village.pdf"))
        plotSA(voomNoNormDupEfit, main="Mean-variance trend elimination with duplicate correction")
    dev.off()

    # get top genes using toptable
    allDEresults <- list()

    for(i in 1:10){
        allDEresults[[i]] <- topTable(voomNoNormDupEfit, coef=i, n=Inf, sort.by="p")
    }

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down         45     1078      407        1      405      137      482     2098     2039      647
# NotSig    12852    10847    12106    12973    12026    12483    12034     8985     8907    11774
# Up           78     1050      462        1      544      355      459     1892     2029      554

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
# Down         15      303       78        0      163        8      138      711      697      136
# NotSig    12906    12112    12685    12974    12489    12924    12659    11818    11358    12606
# Up           54      560      212        1      323       43      178      446      920      233

 summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down          4       47        9        0       52        2       17      121      125        8
# NotSig    12947    12717    12911    12974    12800    12969    12906    12723    12546    12913
# Up           24      211       55        1      123        4       52      131      304       54

for (i in 1:10){
    write.table(allDEresults[[i]], file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.village.", colnames(contr.matrix)[i], ".txt"))
}


###################################################
### 3. DE testing without duplicate correlation ###
###################################################

### Only doing the nonorm one, because why bother with the others?
    voomNoNormVfit <- lmFit(voomNoNorm, design)
    voomNoNormVfit <- contrasts.fit(voomNoNormVfit, contrasts=contr.matrix)
    voomNoNormEfit <- eBayes(voomNoNormVfit, robust=T)

pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.no_dup_correction.mean-variance-trend.village.pdf"))
    plotSA(voomNoNormEfit, main="Mean-variance trend elimination without duplicate correction")
dev.off()

    # get top genes using toptable
    allDEresultsNoDup <- list()

    for(i in 1:10){
        allDEresultsNoDup[[i]] <- topTable(voomNoNormEfit, coef=i, n=Inf, sort.by="p")
    }

summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down         67     1441      483        1      517      133      564     2161     2340      792
# NotSig    12801    10212    11959    12973    11849    12511    11844     8905     8413    11435
# Up          107     1322      533        1      609      331      567     1909     2222      748
 
summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down         19      408       93        1      191       10      159      719      784      166
# NotSig    12886    11911    12637    12973    12432    12924    12597    11789    11205    12494
# Up           70      656      245        1      352       41      219      467      986      315

summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down          3       58        9        1       61        2       21      119      142       13
# NotSig    12946    12674    12902    12973    12782    12969    12881    12727    12501    12891
# Up           26      243       64        1      132        4       73      129      332       71

for (i in 1:10){
    write.table(allDEresultsNoDup[[i]], file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.not_dup_corrected.village.", colnames(contr.matrix)[i], ".txt"))
}

for (i in 1:10){
    withWithout <- join(allDEresults[[i]], allDEresultsNoDup[[i]], by="genes")
    print(paste0("Spearman correlation between methods in ", colnames(contr.matrix)[i], ":"))
    print(cor(withWithout[,6], withWithout[,12], method="spearman"))
}

# [1] "Spearman correlation between methods in ANKvsMDB:"
# [1] 0.9723409
# [1] "Spearman correlation between methods in ANKvsKOR:"
# [1] 0.9721176
# [1] "Spearman correlation between methods in ANKvsTLL:"
# [1] 0.9588628
# [1] "Spearman correlation between methods in ANKvsWNG:"
# [1] 0.9775763
# [1] "Spearman correlation between methods in MDBvsKOR:"
# [1] 0.9868406
# [1] "Spearman correlation between methods in MDBvsTLL:"
# [1] 0.9800908
# [1] "Spearman correlation between methods in WNGvsMDB:"
# [1] 0.9866404
# [1] "Spearman correlation between methods in TLLvsKOR:"
# [1] 0.984165
# [1] "Spearman correlation between methods in WNGvsKOR:"
# [1] 0.9826641
# [1] "Spearman correlation between methods in WNGvsTLL:"
# [1] 0.9747141



#####################################################################################
### 4. Visual QC of duplicate correlation voom output after fitting linear models ###
#####################################################################################

# check to see p-value distribution is normal
pdf(paste0(edaoutput,"PvalueDist_NotAdjusted_dupCor.village.pdf"), height=15, width=10)
    par(mfrow=c(3,1))
    for (i in 1:ncol(voomNoNormDupEfit)){
        hist(voomNoNormDupEfit$p.value[,i], main=colnames(voomNoNormDupEfit)[i], ylim=c(0,max(table(round(voomNoNormDupEfit$p.value[,i], 1)))+1000), xlab="p-value")
    }
dev.off()

# check p-value distribution for adjusted p-values
pdf(paste0(edaoutput,"PvalueDist_Adjusted_dupCor.village.pdf"), height=15, width=10)
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
                                            # hkControls <- hkGenes[which(hkGenes %in% rownames(yVillage$counts))]

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

####################################################
### 5. Summary and visualisation of gene trends. ###
####################################################

# plot log2 fold change between islands
pdf(paste0(edaoutput,"log2FC_VillageComparisons_pval01_dupCor.village.pdf"))
# note 'p.value' is the cutoff value for adjusted p-values
    topTable <- topTable(voomNoNormDupEfit, coef=1, n=Inf, p.value=0.01)
    plot(density(topTable$logFC), col=9, xlim=c(-2,2), main="LogFC Density", xlab="LogFC", ylab="Density", lwd=3, ylim=c(0,1))
    abline(v=c(-1,-0.5,0.5,1), lty=3)
    counter=0
    for (i in 2:ncol(voomNoNormDupEfit)){
        counter=counter+1
        topTable <- topTable(voomNoNormDupEfit, coef=i, n=Inf, p.value=0.01)
        lines(density(topTable$logFC), col=9+counter, xlim=c(-2,2), lwd=3, ylim=c(0,1))
    }
    legend(x="topright", bty="n", col=9:18, legend=colnames(voomNoNormDupEfit), lty=1, lwd=2)
dev.off()

# graphical representation of DE results through MD plot
pdf(paste0(edaoutput,"MD_Plots_pval01_lfc1_dupCor.village.pdf"))
    for(i in 1:ncol(voomNoNormDupEfit)){
        plotMD(voomNoNormDupEfit, column = i, array = NULL, xlab = "Average log-expression", ylab = "Expression log-ratio",
           main = colnames(voomNoNormDupEfit)[i], status=voomNoNormDupEfit$genes$Status, zero.weights = FALSE)
        abline(h=c(-1,-0.5,0.5,1), lty=3)
    }
dev.off()


#######################################################
### 6. UpsetR plots and merging island and village. ###
#######################################################

# Bring in the by Island tests:
# First some temporary renaming so things don't become messed up
villageVoomNoNormDupEfit <- voomNoNormDupEfit
rm(voomNoNormDupEfit)
load(paste0(outputdir, "voomNoNorm.tmm.filtered.dup_corrected.Efit_object.Rda"))
islandVoomNoNormDupEfit <- voomNoNormDupEfit
rm(voomNoNormDupEfit)

byVillages <- decideTests(villageVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01)
byVillages05 <- decideTests(villageVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5)
byVillages1 <- decideTests(villageVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1)

byIslands <- decideTests(islandVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01)
byIslands05 <- decideTests(islandVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5)
byIslands1 <- decideTests(islandVoomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1)

#Now order them both by the same so upsetR can judge them both:
table(rownames(byIslands) == rownames(byVillages))
table(rownames(byIslands05) == rownames(byVillages05))
table(rownames(byIslands1) == rownames(byVillages1))

allTogether <- data.frame(byVillages, byIslands)
allTogether05 <- data.frame(byVillages05, byIslands05)
allTogether1 <- data.frame(byVillages1, byIslands1)

# Look at which genes are in common using UpsetR
# IGR 2019.06.14: Cleaned most of these plots out, they were just cluttering the script. Check repo for past versions.
### Consider using the group.by approach if figures become too messy, but this is neat.

# Making the final figure:
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_fc05_testers.pdf"), width=12)
    # Sort by number of intersects without the island level comparisons, out of curiosity, include all pairwise comparisons
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 1.5, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T) 
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=30,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 1.5, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)
dev.off()

# And now let's give up and make a four way Venn diagram instead:
### UPDATE 2019.06.13: We hate Venn Diagrams. Dropped this section. See below for more exploration.
# Venn diagrams only for the supplementary figure, with three comparisons they convey information better than the equivalent upsetR plot:

# Triple:
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

allTogether05Venn <- allTogether05
allTogether05Venn$genes <- rownames(allTogether05)

    make.venn.triple(allTogether05Venn[allTogether05Venn$ANKvsKOR != 0,]$genes, allTogether05Venn[allTogether05Venn$WNGvsKOR != 0,]$genes, allTogether05Venn[allTogether05Venn$SMBvsKOR != 0,]$genes, paste0(edaoutput, "venn_SMB_vs_KOR_all_levels_fc05"), "ANKvsKOR", "WNGvsKOR", "SMBvsKOR", allTogether05Venn)

    make.venn.triple(allTogether05Venn[allTogether05Venn$MDBvsKOR != 0,]$genes, allTogether05Venn[allTogether05Venn$TLLvsKOR != 0,]$genes, allTogether05Venn[allTogether05Venn$MTWvsKOR != 0,]$genes, paste0(edaoutput, "venn_MTW_vs_KOR_all_levels_fc05"), "MDBvsKOR", "TLLvsKOR", "MTWvsKOR", allTogether05Venn)


##########################################
### 7. Looking at the top ranked genes ###
##########################################

# Let's look at signal across some of the genes that are DE between WNG and KOR, and between SMB and KOR and ANK and KOR, to check what's going there with the villages

# Define plotting function:
singleVillageGenes <- function(singleVillageDF, nGenes, comp1, comp2){
    for (i in 1:nGenes){
        cpmSingle <- voomNoNormDup$E[singleVillageDF$genes[i],]
        singleGene <- data.frame(cpmSingle, yVillage$samples$Sampling.Site)
        names(singleGene)[2] <- "Sampling.Site"

        singleGenePlot <- ggplot(singleGene, aes(x=Sampling.Site, y=cpmSingle, fill=Sampling.Site)) +
            geom_violin(trim=T) + 
            geom_boxplot(width=0.1, fill="white") + 
            geom_jitter(colour = "black", width = 0.2) +
            scale_fill_manual(values=c(sumba, mentawai, korowai, mentawai, sumba)) + 
            scale_x_discrete(labels=c("Anakalung", "Madobag", "Korowai", "Taileleu", "Wunga")) +
            theme_bw() + 
            labs(title=paste0(singleVillageDF$genes[i], ": p = ", signif(singleVillageDF[i,6], 3), " ", comp1, ",\np = ", signif(singleVillageDF[i,12], 3), " ", comp2), y="log CPM", x="") + 
            theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F)

        print(singleGenePlot)
    }
}

# Sumba vs Korowai:
    ankKor <- allDEresults[[2]]
    wngKor <- allDEresults[[9]]

    smbVillageKor <- merge(ankKor, wngKor, by.x="genes", by.y="genes", suffixes=c(".ank", ".wng"))
    smbVillageKor <- smbVillageKor[order(smbVillageKor$adj.P.Val.wng),]

    wngOnly <- smbVillageKor[smbVillageKor$adj.P.Val.wng <= 0.01 & abs(smbVillageKor$logFC.wng) >= 0.5 & (smbVillageKor$adj.P.Val.ank > 0.01 | abs(smbVillageKor$logFC.ank) < 0.5) ,]
    # This one actually recapitulates the results of the venn diagram perfectly, which is what I want - either it is not DE, or the log FC is too low
    wngOnly <- wngOnly[order(wngOnly$adj.P.Val.wng),] # Too lazy to order inside function...

    ankOnly <- smbVillageKor[smbVillageKor$adj.P.Val.ank <= 0.01 & abs(smbVillageKor$logFC.ank) >= 0.5 & (smbVillageKor$adj.P.Val.wng > 0.01 | abs(smbVillageKor$logFC.wng) < 0.5),]
    ankOnly <- ankOnly[order(ankOnly$adj.P.Val.ank),] # Too lazy to order inside function...

    pdf(file=paste0(edaoutput, "wng_kor_only_top_genes.pdf"))
        singleVillageGenes(wngOnly, 30, "ANKvsKOR", "WNGvsKOR")
    dev.off()

    pdf(file=paste0(edaoutput, "ank_kor_only_top_genes.pdf"))
        singleVillageGenes(ankOnly, 30, "ANKvsKOR", "WNGvsKOR")
    dev.off()

# Mentawai Korowai:
    tllKor <- allDEresults[[8]]
    mdbKor <- allDEresults[[5]]

    mtwVillageKor <- merge(tllKor, mdbKor, by.x="genes", by.y="genes", suffixes=c(".tll", ".mdb"))
    mtwVillageKor <- mtwVillageKor[order(mtwVillageKor$adj.P.Val.mdb),]

    mdbOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.mdb <= 0.01 & abs(mtwVillageKor$logFC.mdb) >= 0.5 & (mtwVillageKor$adj.P.Val.tll > 0.01 | abs(mtwVillageKor$logFC.tll) < 0.5) ,]
    # This one actually recapitulates the results of the venn diagram perfectly, which is what I want - either it is not DE, or the log FC is too low
    mdbOnly <- mdbOnly[order(mdbOnly$adj.P.Val.mdb),] # Too lazy to order inside function...

    tllOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.tll <= 0.01 & abs(mtwVillageKor$logFC.tll) >= 0.5 & (mtwVillageKor$adj.P.Val.mdb > 0.01 | abs(mtwVillageKor$logFC.mdb) < 0.5),]
    tllOnly <- tllOnly[order(tllOnly$adj.P.Val.tll),] # Too lazy to order inside function...

    pdf(file=paste0(edaoutput, "mdb_kor_only_top_genes.pdf"))
        singleVillageGenes(mdbOnly, 30, "TLLvsKOR", "MDBvsKOR")
    dev.off()

    pdf(file=paste0(edaoutput, "tll_kor_only_top_genes.pdf"))
        singleVillageGenes(tllOnly, 30, "TLLvsKOR", "MDBvsKOR")
    dev.off()

# What's the rank correlation across the two villages in each island? does it get worse as you go down quintiles/deciles?
    smbKorIslandDE <- topTable(islandVoomNoNormDupEfit, coef=2, n=Inf, sort.by="p")
    mtwKorIslandDE <- topTable(islandVoomNoNormDupEfit, coef=3, n=Inf, sort.by="p")

    # Merge island info with village info, pt 2:
    smbAllKor <- merge(smbVillageKor, smbKorIslandDE, by.x="genes", by.y="genes")
    smbAllKor <- smbAllKor[order(smbAllKor$adj.P.Val),]

    mtwAllKor <- merge(mtwVillageKor, mtwKorIslandDE, by.x="genes", by.y="genes")
    mtwAllKor <- mtwAllKor[order(mtwAllKor$adj.P.Val),]

    cor(smbAllKor[,c(6,12,18)], method="spearman")
    #               adj.P.Val.ank adj.P.Val.wng adj.P.Val
    # adj.P.Val.ank     1.0000000     0.7113743 0.8436221
    # adj.P.Val.wng     0.7113743     1.0000000 0.9171904
    # adj.P.Val         0.8436221     0.9171904 1.0000000

    cor(mtwAllKor[,c(6,12,18)], method="spearman")
    #               adj.P.Val.tll adj.P.Val.mdb adj.P.Val
    # adj.P.Val.tll     1.0000000     0.4735451 0.9053022
    # adj.P.Val.mdb     0.4735451     1.0000000 0.6829288
    # adj.P.Val         0.9053022     0.6829288 1.0000000

    # And now, by quintiles/deciles, determined on the island-wide p-value? mean expression?:

    mtwAllKor$quintile <- ntile(mtwAllKor$adj.P.Val, 5)
    smbAllKor$quintile <- ntile(smbAllKor$adj.P.Val, 5)

    by(mtwAllKor, mtwAllKor$quintile, function(x) cor(x[,c(6,12,18)], method="spearman")) # These are DISMAL
    by(smbAllKor, smbAllKor$quintile, function(x) cor(x[,c(6,12,18)], method="spearman")) # These are not much better

    # What about by logFC, is that better?
    # Condition on DE (look at top quintile etc), consider log FC direction in general, 

    cor(smbAllKor[,c(2,8,14)], method="spearman")
    #           logFC.ank logFC.wng     logFC
    # logFC.ank 1.0000000 0.8805480 0.9392975
    # logFC.wng 0.8805480 1.0000000 0.9701557
    # logFC     0.9392975 0.9701557 1.0000000
    cor(mtwAllKor[,c(2,8,14)], method="spearman")
    #           logFC.tll logFC.mdb     logFC
    # logFC.tll 1.0000000 0.7821874 0.9671235
    # logFC.mdb 0.7821874 1.0000000 0.8878783
    # logFC     0.9671235 0.8878783 1.0000000

    by(mtwAllKor, mtwAllKor$quintile, function(x) cor(x[,c(2,8,14)], method="spearman")) # These are much improved
    by(smbAllKor, smbAllKor$quintile, function(x) cor(x[,c(2,8,14)], method="spearman")) # Yes indeed

# Labels for the plot
quintileLabels <- c("0-20% inter-island\np-value", "20-40%", "40-60%", "60-80%", "80-100% inter-island\np-value")
names(quintileLabels) <- c("1", "2", "3", "4", "5")

    pdf(file=paste0(edaoutput, "village_correlations_by_quintile.pdf"), width=18, height=4)
        # define plotting function:
        ggplot(mtwAllKor, aes(y=-log10(adj.P.Val.tll), x=-log10(adj.P.Val.mdb), group=quintile)) +
            geom_hline(yintercept=-log10(0.1), linetype=2, colour="grey60") +
            geom_vline(xintercept=-log10(0.1), linetype=2, colour="grey60") +
            geom_point(size= 0.7, alpha=0.5) + 
            theme_bw() + 
            labs(title="", y="-log10 pval TLLvsKOR", x="-log10 pval MDBvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ quintile, labeller = labeller(quintile = quintileLabels))

        ggplot(smbAllKor, aes(y=-log10(adj.P.Val.wng), x=-log10(adj.P.Val.ank), group=quintile)) +
            geom_hline(yintercept=-log10(0.1), linetype=2, colour="grey60") +
            geom_vline(xintercept=-log10(0.1), linetype=2, colour="grey60") +
            geom_point(size= 0.7, alpha=0.5) + 
            theme_bw() + 
            labs(title="", y="-log10 pval WNGvsKOR", x="-log10 pval ANKvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ quintile, labeller = labeller(quintile = quintileLabels))

        ggplot(mtwAllKor, aes(y=logFC.tll, x=logFC.mdb, group=quintile)) +
            geom_hline(yintercept=0, linetype=3, colour="grey30") +
            geom_vline(xintercept=0, linetype=3, colour="grey30") +
            geom_hline(yintercept=-0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=-0.5, linetype=2, colour="grey60") +
            geom_hline(yintercept=0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=0.5, linetype=2, colour="grey60") +
            geom_point(size= 0.7, alpha=0.5) + 
            theme_bw() + 
            labs(title="", y="log FC TLLvsKOR", x="log FC MDBvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ quintile, labeller = labeller(quintile = quintileLabels))

        ggplot(smbAllKor, aes(y=logFC.wng, x=logFC.ank, group=quintile)) +
            geom_hline(yintercept=0, linetype=3, colour="grey30") +
            geom_vline(xintercept=0, linetype=3, colour="grey30") +
            geom_hline(yintercept=-0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=-0.5, linetype=2, colour="grey60") +
            geom_hline(yintercept=0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=0.5, linetype=2, colour="grey60") +
            geom_point(size= 0.7, alpha=0.5) + 
            theme_bw() + 
            labs(title="", y="log FC WNGvsKOR", x="log FC ANKvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ quintile, labeller = labeller(quintile = quintileLabels))

    dev.off()    




#########################################################################
### 8. Heatmaps of top DE genes between some villages but not others. ###
#########################################################################

### Decided on phone call that we want to pull maybe top 100 genes or so showing DE between some villages but not others, so testing this out here now:
# pull top 100 genes... but how do we define them? Can do it at the island level, or at the focal village level

allCPM <- voomNoNormDup$E

# Sumba vs Korowai
    # Remove Mentawai
    wngKorCPM <- allCPM[wngOnly$genes[1:100],]
    wngKorCPM <- wngKorCPM[,grepl("MPI|SMB", colnames(wngKorCPM))]

    ankKorCPM <- allCPM[ankOnly$genes[1:86],] # only 86 genes
    ankKorCPM <- ankKorCPM[,grepl("MPI|SMB", colnames(ankKorCPM))]

    smbKorCPM <- allCPM[smbKorIslandDE$genes[1:100],] # This one isn't filtered, but it's ok because the figure is useless.
    smbKorCPM <- smbKorCPM[,grepl("MPI|SMB", colnames(smbKorCPM))]

    # Column annotation - same for all plots
    colMetadata <- yVillage$samples[,c(4,6,7)]
    islandCols <- c("West Papua" = korowai, "Sumba" = sumba, "Mentawai" = mentawai)
    villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod")

    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|SMB", colMetadata$samples),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Rows: 
    # Wunga-centric:
    wngOnly <- merge(wngOnly, smbKorIslandDE, by.x="genes", by.y="genes", all=F, sort=F)
    rowMetadataWng <- data.frame(rownames(wngKorCPM), -log10(wngOnly$adj.P.Val.wng[1:100]), -log10(wngOnly$adj.P.Val.ank[1:100]), -log10(wngOnly$adj.P.Val[1:100]))
    names(rowMetadataWng) <- c("genes", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsWng <- HeatmapAnnotation(rowMetadataWng[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black"))), which="row")

    # Anakalung
    ankOnly <- merge(ankOnly, smbKorIslandDE, by.x="genes", by.y="genes", all=F, sort=F)
    rowMetadataAnk <- data.frame(rownames(ankKorCPM), -log10(ankOnly$adj.P.Val.wng[1:86]), -log10(ankOnly$adj.P.Val.ank[1:86]), -log10(ankOnly$adj.P.Val[1:86]))
    names(rowMetadataAnk) <- c("genes", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsAnk <- HeatmapAnnotation(rowMetadataAnk[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black"))), which="row")

    # For the island-wide one
    rowMetadataSmb <- data.frame(rownames(smbKorCPM), -log10(smbAllKor$adj.P.Val.wng[1:100]), -log10(smbAllKor$adj.P.Val.ank[1:100]), -log10(smbAllKor$adj.P.Val[1:100]))
    names(rowMetadataSmb) <- c("genes", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsSmb <- HeatmapAnnotation(rowMetadataSmb[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black"))), which="row")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "smb_kor_DE_heatmaps.pdf"), height=8, width=6)
        wngMap <- Heatmap(t(scale(t(wngKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(wngKorCPM)), cluster_columns=F)
        draw(rowColsWng + wngMap, row_dend_side = "left")

        ankMap <- Heatmap(t(scale(t(ankKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(ankKorCPM)), cluster_columns=F)
        draw(rowColsAnk + ankMap, row_dend_side = "left")

        smbMap <- Heatmap(t(scale(t(smbKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(smbKorCPM)), cluster_columns=F)
        draw(rowColsSmb + smbMap, row_dend_side = "left")
    dev.off()

### Now Mentawai vs Korowai:
    # Remove Sumba
    tllKorCPM <- allCPM[tllOnly$genes[1:100],]
    tllKorCPM <- tllKorCPM[,grepl("MPI|MTW", colnames(tllKorCPM))]

    mdbKorCPM <- allCPM[mdbOnly$genes[1:100],]
    mdbKorCPM <- mdbKorCPM[,grepl("MPI|MTW", colnames(mdbKorCPM))]

    mtwKorCPM <- allCPM[mtwKorIslandDE$genes[1:100],]
    mtwKorCPM <- mtwKorCPM[,grepl("MPI|MTW", colnames(mtwKorCPM))]

    # Column annotation - same for all plots
    colMetadata <- yVillage$samples[,c(4,6,7)]
    islandCols <- c("West Papua" = korowai, "Sumba" = sumba, "Mentawai" = mentawai)
    villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod")

    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|MTW", colMetadata$samples),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Rows: 
    # Taileleu-centric:
    tllOnly <- merge(tllOnly, mtwKorIslandDE, by.x="genes", by.y="genes", all=F, sort=F)
    rowMetadataTll <- data.frame(rownames(tllKorCPM), -log10(tllOnly$adj.P.Val.tll[1:100]), -log10(tllOnly$adj.P.Val.mdb[1:100]), -log10(tllOnly$adj.P.Val[1:100]))
    names(rowMetadataTll) <- c("genes", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsTll <- HeatmapAnnotation(rowMetadataTll[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black"))), which="row")

    # Madobag
    mdbOnly <- merge(mdbOnly, mtwKorIslandDE, by.x="genes", by.y="genes", all=F, sort=F)
    rowMetadataMdb <- data.frame(rownames(mdbKorCPM), -log10(mdbOnly$adj.P.Val.tll[1:100]), -log10(mdbOnly$adj.P.Val.mdb[1:100]), -log10(mdbOnly$adj.P.Val[1:100]))
    names(rowMetadataMdb) <- c("genes", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsMdb <- HeatmapAnnotation(rowMetadataMdb[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black"))), which="row")

    # Island-wide
    rowMetadataMtw <- data.frame(rownames(mtwKorCPM), -log10(mtwAllKor$adj.P.Val.tll[1:100]), -log10(mtwAllKor$adj.P.Val.mdb[1:100]), -log10(mtwAllKor$adj.P.Val[1:100]))
    names(rowMetadataMtw) <- c("genes", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsMtw <- HeatmapAnnotation(rowMetadataMtw[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black"))), which="row")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "mtw_kor_DE_heatmaps.pdf"), height=8, width=6)
        tllMap <- Heatmap(t(scale(t(tllKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(tllKorCPM)), cluster_columns=F)
        draw(rowColsTll + tllMap, row_dend_side = "left")

        mdbMap <- Heatmap(t(scale(t(mdbKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mdbKorCPM)), cluster_columns=F)
        draw(rowColsMdb + mdbMap, row_dend_side = "left")

        mtwMap <- Heatmap(t(scale(t(mtwKorCPM))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mtwKorCPM)), cluster_columns=F)
        draw(rowColsMtw + mtwMap, row_dend_side = "left")
    dev.off()


#######################################################################################################
### 9. Quick check of variance by village, to see what drives the weird distribution of DE results. ###
#######################################################################################################

# load(paste0(outputdir, "voomNoNorm.tmm.filtered.duplicate_corrected.indoRNA.Rda"))

# Define CoV function:
calcCoV <- function(x){
     (sd(x)/mean(x) )
 }

# Can't use CoV with the log transformation, need to undo it: (wikipedia said so, and yes, the negative numbers were probably messing things up)
normLCPM <- data.frame(t(voomNoNormDup$E))
normCPM <- normLCPM^2

perVillageCoV <- t(ddply(normCPM, .(yVillage$samples$Sampling.Site), function(x) apply(as.matrix(x), 2, function(x) calcCoV(x)))) # That's a lot of transposing, but I checked it manually.

# Clean up all that messy output...
perVillageCoVDF <- data.frame(perVillageCoV[-1,], stringsAsFactors=F)
names(perVillageCoVDF) <- perVillageCoV[1,]
perVillageCoVDF <- as.data.frame(lapply(perVillageCoVDF, as.numeric)) # Worked fine when checking lines manually.
summary(perVillageCoVDF)

dataForPlotting <- melt(perVillageCoVDF)

pdf(paste0(edaoutput, "cov_by_village_densities.pdf"))
    ggplot(dataForPlotting, aes(x=value, fill=variable)) +
        geom_density(alpha=0.3) + 
        # geom_boxplot(width=0.05, fill="white") + 
        scale_fill_manual(values=c(sumba, mentawai, korowai, mentawai, sumba)) + 
        # scale_x_discrete(labels=c("Anakalung", "Madobag", "Korowai", "Taileleu", "Wunga")) +
        theme_bw() + 
        labs(title="", x="CoV all genes", y="Density") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(fill=F)
dev.off()

# Well that doesn't look too different... so why that DE trend?

# Lots of t.tests:

tTestOut <- matrix(0, nrow=5, ncol=5)
colnames(tTestOut) <- names(perVillageCoVDF)
rownames(tTestOut) <- names(perVillageCoVDF)

for (i in 1:5){
    for (j in 1:5){
        if(j > i){
            tTestOut[i,j] <- t.test(perVillageCoVDF[,i], perVillageCoVDF[,j])$p.value
        }
    }
}

# They are different, as they were of course going to be, but not THAT different, surely, to explain the difference in power? Also from the plot the effect size does not go in the direction you would expect. The things that are different are not different 
signif(tTestOut, digits=3)
#           Anakalung  Madobag    Korowai Taileleu    Wunga
# Anakalung         0 2.73e-08 1.26e-02 2.66e-01 4.36e-04
# Madobag           0 0.00e+00 2.29e-16 9.70e-06 4.30e-02
# Korowai             0 0.00e+00 0.00e+00 2.89e-04 1.07e-09
# Taileleu          0 0.00e+00 0.00e+00 0.00e+00 1.66e-02
# Wunga             0 0.00e+00 0.00e+00 0.00e+00 0.00e+00

# Similar plots of pairwise correlations within each village, to see if anything is as noisy as Korowai. But then how do you reconcile the CoV observations?


# And then, for a bit of overkill, plots of CoVs across all villages subset by genes that are only DE in a single village.
# Add row names so you can subset by genes:
rownames(perVillageCoVDF) <- rownames(perVillageCoV)[-1]

villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod")

# Define plotting function
plotCoV <- function(inputDF, comparison){
    dataForPlotting <- melt(inputDF)

    covOverkill <- ggplot(dataForPlotting, aes(x=value, fill=variable)) +
            geom_density(alpha=0.5) + 
            scale_fill_manual(name = "",values = villageCols) +
            theme_bw() + 
            labs(title=paste0("DE ", comparison, " (", nrow(inputDF), " genes)"), x="CoV by gene", x="") + 
            theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
    return(covOverkill)

    # print(covOverkill)
}

WngNoAnk <- plotCoV(perVillageCoVDF[wngOnly$genes,c("Wunga", "Anakalung")], "WNGvsKOR not ANKvsKOR")
AnkNoWng <- plotCoV(perVillageCoVDF[ankOnly$genes,c("Wunga", "Anakalung")], "ANKvsKOR not WNGvsKOR")
MdbNoTll <- plotCoV(perVillageCoVDF[mdbOnly$genes,c("Madobag", "Taileleu")], "MDBvsKOR not TLLvsKOR")
TllNoMdb <- plotCoV(perVillageCoVDF[tllOnly$genes,c("Madobag", "Taileleu")], "TLLvsKOR not MDBvsKOR")

pdf(paste0(edaoutput, "cov_by_single_village_DE_densities.pdf"), height=8, width=8)
    # plotCoV(perVillageCoVDF[wngOnly$genes,c("Wunga", "Anakalung")], "WNGvsKOR not ANKvsKOR")
    # plotCoV(perVillageCoVDF[ankOnly$genes,c("Wunga", "Anakalung")], "ANKvsKOR not WNGvsKOR")
    # plotCoV(perVillageCoVDF[mdbOnly$genes,c("Madobag", "Taileleu")], "MDBvsKOR not TLLvsKOR")
    # plotCoV(perVillageCoVDF[tllOnly$genes,c("Madobag", "Taileleu")], "TLLvsKOR not MDBvsKOR")
    ggarrange(WngNoAnk, AnkNoWng, MdbNoTll, TllNoMdb, ncol=2, nrow=2, common.legend=T)
dev.off()


# These need the t-test

###########################################################################################
### 10. Good old plot of pairwise correlations within each village and level etc etc... ###
#############################################################################################################

# Define hideous function
plot.reproducibility <- function(data.to.test, metadata, method){
    corMat <- cor(data.to.test, method=method, use="pairwise.complete.obs")

    indRep <- vector()
    villageBatchRep <- vector()
    islandBatchRep <- vector()
    batchRep <- vector()
    villageNoBatchRep <- vector()
    islandNoBatchRep <- vector()
    betweenBatch <- vector()

    for (i in 1:ncol(data.to.test)){
        for (j in 1:ncol(data.to.test)){
            if (j > i){
                if (metadata$ID[i] == metadata$ID[j]) {
                    indRep <- c(indRep, corMat[i,j])
                } else if (metadata$batch[i] == metadata$batch[j]){
                    if (metadata$Sampling.Site[i] == metadata$Sampling.Site[j]){
                        villageBatchRep <- c(villageBatchRep, corMat[i,j])
                    } else if (metadata$Island[i] == metadata$Island[j]){
                        islandBatchRep <- c(islandBatchRep, corMat[i,j])
                    } else
                        batchRep <- c(batchRep, corMat[i,j])
                } else if (metadata$batch[i] != metadata$batch[j]){
                    if (metadata$Sampling.Site[i] == metadata$Sampling.Site[j]){
                        villageNoBatchRep <- c(villageNoBatchRep, corMat[i,j])
                    } else if (metadata$Island[i] == metadata$Island[j]){
                        islandNoBatchRep <- c(islandNoBatchRep, corMat[i,j])
                    } else {betweenBatch <- c(betweenBatch, corMat[i,j])}
                }
            }
        }
    }

    forPlot <- melt(list(indRep, villageBatchRep, villageNoBatchRep, islandBatchRep, islandNoBatchRep, batchRep, betweenBatch))
    forPlot$L1 <- as.factor(forPlot$L1)

    ggplot(forPlot, aes(x=L1, y=value, fill=L1)) +
        geom_violin(trim=T) + 
        geom_boxplot(width=0.05, fill="white") + 
        scale_x_discrete(labels=c("Individual reps", "within village\nand batch", "within island\nand batch", "diff island\nwithin batch", "within village\ndiff batch", "within island\ndiff batch", "diff island\ndiff batch")) +
        theme_bw() + 
        labs(title="", y=paste0(method, " pairwise correlation"), x="") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(fill=F)

}

# And now, by village and by island plots (separed by batch)
plotWithinSite <- function(data.to.test, metadata, method){
    corMat <- cor(data.to.test, method=method, use="pairwise.complete.obs")

    indRep <- vector()
    smbBatchRep <- vector()
    wngBatchRep <- vector()
    ankBatchRep <- vector()
    mtwBatchRep <- vector()
    mdbBatchRep <- vector()
    tllBatchRep <- vector()
    korBatchRep <- vector()
    smbNoBatchRep <- vector()
    wngNoBatchRep <- vector()
    ankNoBatchRep <- vector()
    mtwNoBatchRep <- vector()
    mdbNoBatchRep <- vector()
    tllNoBatchRep <- vector()
    korNoBatchRep <- vector()

    for (i in 1:ncol(data.to.test)){
        for (j in 1:ncol(data.to.test)){
            if (j > i){
                if (metadata$ID[i] == metadata$ID[j]) {
                    indRep <- c(indRep, corMat[i,j])
                } else if (metadata$batch[i] == metadata$batch[j]){
                    if (metadata$Sampling.Site[i] == metadata$Sampling.Site[j]){
                        if(metadata$Sampling.Site[i] == "Wunga"){
                            wngBatchRep <- c(wngBatchRep, corMat[i,j])
                            smbBatchRep <- c(smbBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Anakalung"){
                            ankBatchRep <- c(ankBatchRep, corMat[i,j])
                            smbBatchRep <- c(smbBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Madobag"){
                            mdbBatchRep <- c(mdbBatchRep, corMat[i,j])
                            mtwBatchRep <- c(mtwBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Taileleu"){
                            tllBatchRep <- c(tllBatchRep, corMat[i,j])
                            mtwBatchRep <- c(mtwBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Korowai"){
                            korBatchRep <- c(korBatchRep, corMat[i,j])
                        }
                    }
                } else if (metadata$batch[i] != metadata$batch[j]){
                    if (metadata$Sampling.Site[i] == metadata$Sampling.Site[j]){
                        if(metadata$Sampling.Site[i] == "Wunga"){
                            wngNoBatchRep <- c(wngNoBatchRep, corMat[i,j])
                            smbNoBatchRep <- c(smbNoBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Anakalung"){
                            ankNoBatchRep <- c(ankNoBatchRep, corMat[i,j])
                            smbNoBatchRep <- c(smbNoBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Madobag"){
                            mdbNoBatchRep <- c(mdbNoBatchRep, corMat[i,j])
                            mtwNoBatchRep <- c(mtwNoBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Taileleu"){
                            tllNoBatchRep <- c(tllNoBatchRep, corMat[i,j])
                            mtwNoBatchRep <- c(mtwNoBatchRep, corMat[i,j])
                        } else if(metadata$Sampling.Site[i] == "Korowai"){
                            korNoBatchRep <- c(korNoBatchRep, corMat[i,j])
                        }
                    }
                }
            }
        }
    }

    forPlot <- melt(list(ankBatchRep, ankNoBatchRep, wngBatchRep, wngNoBatchRep, mdbBatchRep, mdbNoBatchRep, tllBatchRep, tllNoBatchRep, korBatchRep, korNoBatchRep))
    forPlot$L1 <- as.factor(forPlot$L1)

    forPlotIsland <- melt(list(smbBatchRep, mtwBatchRep, korBatchRep, smbNoBatchRep, mtwNoBatchRep, korNoBatchRep))
    forPlotIsland$L1 <- as.factor(forPlotIsland$L1)

    byVillagePlot <- ggplot(forPlot, aes(x=L1, y=value, fill=L1)) +
        geom_violin(trim=T) + 
        geom_boxplot(width=0.05, fill="white") + 
        scale_x_discrete(labels=c("ANK within batch", "ANK bw batch", "WNG within batch", "WNG bw batch", "MDB within batch", "MDB bw batch", "TLL within batch", "TLL bw batch", "KOR within batch", "KOR bw batch")) +
        theme_bw() + 
        labs(title="", y=paste0(method, "pairwise correlation"), x="") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(fill=F)

    byIslandPlot <- ggplot(forPlotIsland, aes(x=L1, y=value, fill=L1)) +
        geom_violin(trim=T) + 
        geom_boxplot(width=0.05, fill="white") + 
        scale_x_discrete(labels=c("SMB within batch", "MTW within batch", "KOR within batch", "SMB bw batch", "MTW bw batch", "KOR bw batch")) +
        theme_bw() + 
        labs(title="", y=paste0(method, " pairwise correlation"), x="") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(fill=F)

    print(byVillagePlot)
    print(byIslandPlot)

}

pdf(file=paste0(edaoutput, "reproducibility_by_levels.pdf"))
    plot.reproducibility(voomNoNormDup$E, yVillage$samples, "spearman")
    plot.reproducibility(voomNoNormDup$E, yVillage$samples, "pearson")
dev.off()

pdf(file=paste0(edaoutput, "correlation_within_sites.pdf"))
    plotWithinSite(voomNoNormDup$E, yVillage$samples, "spearman")
    plotWithinSite(voomNoNormDup$E, yVillage$samples, "pearson")
dev.off()

### Final thought: should we sample everyone down to 15 inds, without reps? 
