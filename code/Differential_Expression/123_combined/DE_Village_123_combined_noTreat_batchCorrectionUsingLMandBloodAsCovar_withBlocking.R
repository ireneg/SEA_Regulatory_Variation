# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: IGR 2019.06.03 
### Changed MPI to KOR, Mappi to Korowai throughout - see line 81. Reran all figures; file names will differ between old and new versions but the contents should be identical. 

### 0. Load dependencies and functions and set input paths -------------------------- ###
### 1. Begin analyses and initial QC ---------------------------------------------------------------------------------- ###
### 2. DE testing with duplicate correlation and blocking ----------------------------------------------------- ###
### 3. DE testing without duplicate correlation ------------------------------------------------- ###
### 4. Visual QC of duplicate correlation voom output after fitting linear models ---------------- ###
### 5. Summary and visualisation of gene trends --------------------------------------------- ###
### 6. Looking at the top ranked genes ------------------------------------------------- ###
### 7. Quick check of variance by village, to see what drives the weird distribution of DE results. ------ ###
### 8. Good old plot of pairwise correlations within each village and level etc etc... ------------------ ###


### TO DO:
### Fix everything that's commented out (just figures)
### Triple check all numbers.

#########################################################################################
### 0. Load dependencies and functions and set input paths -------------------------- ###
#########################################################################################

# Load dependencies:
library(edgeR)
library(plyr)
# library(NineteenEightyR)
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
library(UpSetR)
library(matrixStats)
library(reshape)
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


###########################################################################################################################
### 1. Begin analyses and initial QC ---------------------------------------------------------------------------------- ###
###########################################################################################################################

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

###################################################################################################################
### 2. DE testing with duplicate correlation and blocking ----------------------------------------------------- ###
###################################################################################################################

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


#####################################################################################################
### 3. DE testing without duplicate correlation ------------------------------------------------- ###
#####################################################################################################

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



######################################################################################################
### 4. Visual QC of duplicate correlation voom output after fitting linear models ---------------- ###
######################################################################################################

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

#################################################################################################
### 5. Summary and visualisation of gene trends --------------------------------------------- ###
#################################################################################################

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

# # We can also make individual pdfs of the top genes
# island1 <- c("Sumba","Mentawai","West Papua")
# island2 <- c("Sumba","Mentawai","West Papua")

# counter <- 0

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

# Bring in the by Island tests:
# First some temporary renaming so things don't become messed up:

#    load(voomNoNormDupEfit, file=paste0(outputdir, "voomNoNorm.tmm.filtered.duplicate_corrected.indoRNA.Rda"))


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
### Consider using the group.by approach if figures become too messy, but this is neat.

pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_by_village_allfcs.pdf"), width=12)
    upset(as.data.frame(abs(byVillages)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(byVillages05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(byVillages1)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
dev.off()

pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_by_island_allfcs.pdf"), width=12)
    upset(as.data.frame(abs(allTogether)), sets = c("SMBvsMTW", "SMBvsKOR", "MTWvsKOR"), sets.bar.color = c(smb_mtw, smb_kor, mtw_kor), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsMTW", "SMBvsKOR", "MTWvsKOR"), sets.bar.color = c(smb_mtw, smb_kor, mtw_kor), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether1)), sets = c("SMBvsMTW", "SMBvsKOR", "MTWvsKOR"), sets.bar.color = c(smb_mtw, smb_kor, mtw_kor), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
dev.off()

pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_allfcs.pdf"), width=12)
    upset(as.data.frame(abs(allTogether)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), rep(smb_kor, 3), rep(mtw_kor, 3), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), rep(smb_kor, 3), rep(mtw_kor, 3), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), rep(smb_kor, 3), rep(mtw_kor, 3), sumba, mentawai), nintersects=40,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), rep(smb_kor, 3), rep(mtw_kor, 3), sumba, mentawai), nintersects=30,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether1)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), rep(smb_kor, 3), rep(mtw_kor, 3), sumba, mentawai), nintersects=30,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
dev.off()

### And now, focusing only on each inter-island comparison:
### SMB-MTW
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_allfcs_SMB_MTW.pdf"), width=12)
    upset(as.data.frame(abs(allTogether)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether1)), sets = c("SMBvsMTW", "ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,5), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
dev.off()

### SMB-KOR
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_allfcs_SMB_KOR.pdf"), width=12)
    upset(as.data.frame(abs(allTogether)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "ANKvsWNG"), sets.bar.color = c(rep(smb_kor, 3), sumba), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "ANKvsWNG"), sets.bar.color = c(rep(smb_kor, 3), sumba), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
    upset(as.data.frame(abs(allTogether1)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "ANKvsWNG"), sets.bar.color = c(rep(smb_kor, 3), sumba), nintersects=50,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2)
dev.off()

###MTW-KOR
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_allfcs_MTW_KOR.pdf"), width=12)
    upset(as.data.frame(abs(allTogether)), sets = c("MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "MDBvsTLL"), sets.bar.color = c(rep(mtw_kor, 3), mentawai), nintersects=100,  order.by = "freq", keep.order=T)
    upset(as.data.frame(abs(allTogether05)), sets = c("MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "MDBvsTLL"), sets.bar.color = c(rep(mtw_kor, 3), mentawai), nintersects=100,  order.by = "freq", keep.order=T)
    upset(as.data.frame(abs(allTogether1)), sets = c("MTWvsKOR", "MDBvsKOR", "TLLvsKOR", "MDBvsTLL"), sets.bar.color = c(rep(mtw_kor, 3), mentawai), nintersects=100,  order.by = "freq", keep.order=T)
dev.off()


# Making the final figure:
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_fc05_testers.pdf"), width=12)
    # Sort by Frequency of set, only
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=10,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3))
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=20,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3))
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=30,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3))

    # Sort by number of intersects
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=10,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=20,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=30,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)

    # Sort by number of intersects, decreasing = F - here the set number has to be higher, because 10 is useless
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=10,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=F)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=20,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=F)
    upset(as.data.frame(abs(allTogether05)), sets = c("SMBvsKOR", "ANKvsKOR", "WNGvsKOR", "MTWvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 3), rep(mtw_kor, 3)), nintersects=30,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=F)


    # Sort by number of intersects without the island level comparisons, out of curiosity (16 total comparisons, so show only that one).
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)

    # Sort by both things without the island level comparisons, out of curiosity (16 total comparisons, so lose the 30)
    # Dropped because for some reason it replicates some of the intersects
    # upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)

    # Sort by number of intersects without the island level comparisons, out of curiosity, decreasing = F
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = "degree", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=F)

    # Group by sets, no island-level:
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  group.by = "sets", cutoff=5, keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3))

    # Sort by number of intersects and height of bars (decreasing doesn't work here, which makes it look terrible.)
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = c("degree", "freq"), keep.order=T, number.angles = 30, point.size = 3.5, line.size = 2, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3))
dev.off()



    ###################################################################################################
    ### IGR NOTE 2019.04.12 - I BELIEVE HEINI IS NOW MAKING THIS FIGURE, WILL REVISIT IT OTHERWISE. ###
    ###################################################################################################

# # get DE genes in common with populations compared to Korowai, i.e., SMBvsKOR and MTWvsKOR (since we think this is an interesting island comparison)
# allGenes <- merge(voomNoNormDupTopTableSMB.MTW, voomNoNormDupTopTableSMB.KOR, by.x="genes", by.y="genes", suffixes=c(".SMB.MTW", ".SMB.KOR"))
# allGenes <- merge(allGenes, voomNoNormDupTopTableMTW.KOR, by.x="genes", by.y="genes")
# names(allGenes)[13:19] <- paste0(names(allGenes)[13:19], ".MTW.KOR")

# deSummaryAll <- decideTests(voomNoNormDupEfit, p.value=0.01)
# deSummary05 <- decideTests(voomNoNormDupEfit, p.value=0.01, lfc=0.5)
# deSummary1 <- decideTests(voomNoNormDupEfit, p.value=0.01, lfc=1)

# deCommonKOR = which(deSummaryAll[,2]!=0 & deSummaryAll[,3]!=0)
# deCommonKOR05 = which(deSummary05[,2]!=0 & deSummary05[,3]!=0)
# deCommonKOR1 = which(deSummary1[,2]!=0 & deSummary1[,3]!=0)

# # get what these genes are doing and save them to a file
# # commonGenes.KOR <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common.KOR), filters="ensembl_gene_id")
# # write.table(de.common.KOR, file=paste0(outputdir,"allCommonGenes_KOR_dupcor.txt"))
# # # save the common gene names 
# # de.common.KOR=voomNoNormDupEfit$genes[names(de.common.KOR),]

# # now plot the common genes to see if they're being regulated in the same direction
# pdf(paste0(edaoutput,"logFC_commonKORgenes_dupCor.pdf"))
#     plot(voomNoNormDupTopTableSMB.KOR[rownames(deCommonKOR), "logFC"], voomNoNormDupTopTableMTW.KOR[rownames(deCommonKOR), "logFC"], xlab="logFC SMBvsKOR", ylab="logFC MTWvsKOR", pch=20, main="Common DE Genes", xlim=c(-5,5), ylim=c(-6,6))
#     # text(tt.SMBvsKOR[rownames(de.common.KOR),"logFC"], tt.MTWvsKOR[rownames(de.common.KOR),"logFC"], labels=tt.SMBvsKOR[rownames(de.common.KOR),"SYMBOL"], pos=3)
#     abline(h=0,v=0, lty=2)
# dev.off()


############################################################################################
### 6. Looking at the top ranked genes ------------------------------------------------- ###
############################################################################################

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
    AnkKor <- allDEresults[[2]]
    WngKor <- allDEresults[[9]]

    smbVillageKor <- merge(AnkKor, WngKor, by.x="genes", by.y="genes", suffixes=c(".ank", ".wng"))
    smbVillageKor <- smbVillageKor[order(smbVillageKor$adj.P.Val.wng),]

    # cor(smbVillageKor[,6], smbVillageKor[,12], method="pearson")
    # # [1] 0.5513258
    # cor(smbVillageKor[,6], smbVillageKor[,12], method="spearman")
    # # [1] 0.7113743

    wngOnly <- smbVillageKor[smbVillageKor$adj.P.Val.wng <= 0.01 & smbVillageKor$adj.P.Val.ank > 0.01,]
    wngOnly <- wngOnly[order(wngOnly$adj.P.Val.wng),] # Too lazy to order inside function...

    ankOnly <- smbVillageKor[smbVillageKor$adj.P.Val.ank <= 0.01 & smbVillageKor$adj.P.Val.wng > 0.01,]
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

    # cor(mtwVillageKor[,6], mtwVillageKor[,12], method="pearson")
    # # [1] 0.5513258
    # cor(mtwVillageKor[,6], mtwVillageKor[,12], method="spearman")
    # # [1] 0.7113743

    mdbOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.mdb <= 0.01 & mtwVillageKor$adj.P.Val.tll > 0.01,]
    mdbOnly <- mdbOnly[order(mdbOnly$adj.P.Val.mdb),] # Too lazy to order inside function...

    tllOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.tll <= 0.01 & mtwVillageKor$adj.P.Val.mdb > 0.01,]
    tllOnly <- tllOnly[order(tllOnly$adj.P.Val.tll),] # Too lazy to order inside function...

    pdf(file=paste0(edaoutput, "mdb_kor_only_top_genes.pdf"))
        singleVillageGenes(mdbOnly, 30, "TLLvsKOR", "MDBvsKOR")
    dev.off()

    pdf(file=paste0(edaoutput, "tll_kor_only_top_genes.pdf"))
        singleVillageGenes(tllOnly, 30, "TLLvsKOR", "MDBvsKOR")
    dev.off()


##############################################################################################################
### 7. Quick check of variance by village, to see what drives the weird distribution of DE results. ------ ###
##############################################################################################################

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

pdf(paste0(edaoutput, "cov_by_village.pdf"))
    ggplot(dataForPlotting, aes(x=variable, y=value, fill=variable)) +
        geom_violin(trim=T) + 
        geom_boxplot(width=0.05, fill="white") + 
        scale_fill_manual(values=c(sumba, mentawai, korowai, mentawai, sumba)) + 
        scale_x_discrete(labels=c("Anakalung", "Madobag", "Korowai", "Taileleu", "Wunga")) +
        theme_bw() + 
        labs(title="", y="CoV CPM", x="") + 
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

# Define plotting function
plotCoV <- function(inputDF, comparison){
    dataForPlotting <- melt(inputDF)

    covOverkill <- ggplot(dataForPlotting, aes(x=variable, y=value, fill=variable)) +
            geom_violin(trim=T) + 
            geom_boxplot(width=0.05, fill="white") + 
            scale_fill_manual(values=c(sumba, mentawai, korowai, mentawai, sumba)) + 
            scale_x_discrete(labels=c("Anakalung", "Madobag", "Korowai", "Taileleu", "Wunga")) +
            theme_bw() + 
            labs(title=paste0("DE ", comparison, " (", nrow(inputDF), " genes)"), y="CoV CPM", x="") + 
            theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F)
    print(covOverkill)
}

pdf(paste0(edaoutput, "cov_by_single_village_DE.pdf"))
    plotCoV(perVillageCoVDF[wngOnly$genes,], "WNGvsKOR not ANKvsKOR")
    plotCoV(perVillageCoVDF[ankOnly$genes,], "ANKvsKOR not WNGvsKOR")
    plotCoV(perVillageCoVDF[mdbOnly$genes,], "MDBvsKOR not TLLvsKOR")
    plotCoV(perVillageCoVDF[tllOnly$genes,], "TLLvsKOR not MDBvsKOR")
dev.off()


#############################################################################################################
### 8. Good old plot of pairwise correlations within each village and level etc etc... ------------------ ###
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
