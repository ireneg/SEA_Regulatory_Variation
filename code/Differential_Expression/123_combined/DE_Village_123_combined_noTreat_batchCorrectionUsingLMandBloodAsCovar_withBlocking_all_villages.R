# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

### Last edit: IGR 2019.10.19
### Changed paths to deal with removal of MPI-296

### 0. Load dependencies and functions and set input paths -------------------------- 
### 1. Begin analyses and initial QC ---------------------------------------------------------------------------------- 
### 2. DE testing with duplicate correlation and blocking ----------------------------------------------------- 
### 3. And now with random subsetting for power reasons... -----------------------------

### TO DO:
### Fix everything that's commented out (just figures)
### Triple check all numbers.

### Repeat this 1000 times.

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


# Set paths:
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296" # on server
covariatedir <- "/data/cephfs/punim0586/igallego/indoRNA/"

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/de_testing/no_mpi296"
edaoutput <- paste0(outputdir, "/eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

# Load colour schemes:
mappi <- wes_palette("Zissou1", 20, type = "continuous")[20]
mentawai <- wes_palette("Zissou1", 20, type = "continuous")[1]
sumba <- wes_palette("Zissou1", 20, type = "continuous")[11]

smb_mtw <- wes_palette("Darjeeling1", 9, type = "continuous")[3]
smb_mpi <- wes_palette("Darjeeling1", 9, type = "continuous")[7]
mtw_mpi <- "darkorchid4"


# Load log CPM matrix and y object:
# lcpm
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
# y DGE list object
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))


########################################
### 1. Begin analyses and initial QC ###
########################################

# First, remove samples that have less than ten individuals per village
table(yFilt$samples$Sampling.Site)
    #    Anakalung    Bilarenge    Hupu Mada      Madobag        Mappi  Padira Tana 
    #           20            1            5           17           21            3 
    # Patiala Bawa        Rindi     Taileleu        Wunga   Wura Homba 
    #            1            5           32           17            1 

# remove Bilarenge, Hupu Mada, Padira Tana, Patiala Bawa, Rindi, and Wura Homba
yVillage <- yFilt[,-grep("Bilarenge|Patiala Bawa|Wura Homba", yFilt$samples$Sampling.Site)]
# drop unused levels
yVillage$samples <- droplevels(yVillage$samples)

# Set up design matrix
design <- model.matrix(~0 + yVillage$samples$Sampling.Site + yVillage$samples$Age + yVillage$samples$batch + yVillage$samples$RIN + yVillage$samples$CD8T + yVillage$samples$CD4T + yVillage$samples$NK + yVillage$samples$Bcell + yVillage$samples$Mono + yVillage$samples$Gran)
# rename columns to exclude spaces and unrecognised characters
colnames(design)=gsub("yVillage\\$samples\\$", "", colnames(design))
colnames(design)=gsub("Sampling.Site", "", colnames(design))
colnames(design)=gsub(" ", "_", colnames(design))

# set up contrast matrix
contr.matrix <- makeContrasts(  ANKvsMDB=Anakalung-Madobag, ANKvsMPI=Anakalung-Mappi, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, ANKvsRIN = Anakalung-Rindi, ANKvsHPM = Anakalung-Hupu_Mada, ANKvsPDT = Anakalung-Padira_Tana,
    WNGvsMDB=Wunga-Madobag, WNGvsMPI=Wunga-Mappi, WNGvsTLL=Wunga-Taileleu, WNGvsRIN = Wunga-Rindi, WNGvsHPM = Wunga-Hupu_Mada, WNGvsPDT = Wunga-Padira_Tana,
    RINvsMDB= Rindi-Madobag, RINvsTLL= Rindi-Taileleu, RINvsMPI= Rindi-Mappi, RINvsHPM = Rindi-Hupu_Mada, RINvsPDT=Rindi-Padira_Tana,
    HPMvsMDB= Hupu_Mada-Madobag, HPMvsTLL= Hupu_Mada-Taileleu, HPMvsMPI= Hupu_Mada-Mappi, HPMvsPDT=Hupu_Mada-Padira_Tana,
    PDTvsMDB = Padira_Tana-Madobag, PDTvsTLL = Padira_Tana-Taileleu, PDTvsMPI = Padira_Tana-Mappi,
    MDBvsMPI=Madobag-Mappi, MDBvsTLL=Madobag-Taileleu,   
    TLLvsMPI=Taileleu-Mappi, 
    levels=colnames(design)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

yVillage <- calcNormFactors(yVillage, method="TMM")

#############################################################
### 2. DE testing with duplicate correlation and blocking ###
#############################################################

# create a new variable for blocking using sample IDs
yVillage$samples$ind <- sapply(strsplit(as.character(yVillage$samples$samples), "[_.]"), `[`, 1)

# First, we need to perform voom normalisation

# No normalisation between samples beyond tmm and voom:
    voomNoNorm <- voom(yVillage, design, normalize.method="none", plot=F) 
    dupcorNone <- duplicateCorrelation(voomNoNorm, design, block=yVillage$samples$ind) # 19 non-convergences
    # The value dupcor$consensus estimates the average correlation within the blocks and should be positive
    dupcorNone$consensus # sanity check
    # [1] 0.6866513
    median(voomNoNorm$weights) # another sanity check:
    # [1] 23.6336
    save(voomNoNorm, file=paste0(outputdir, "voomNoNorm.tmm.filtered.indoRNA.village.all_villages.Rda"))

    # Second round:
    voomNoNormDup <- voom(yVillage, design, plot=TRUE, block=yVillage$samples$ind, correlation=dupcorNone$consensus)
    dupcorNoneDup <- duplicateCorrelation(voomNoNormDup, design, block=yVillage$samples$ind) # 19 non convergences
    dupcorNoneDup$consensus # sanity check pt 2
    # [1] 0.6870139
    median(voomNoNormDup$weights) # another sanity check, pt 2 
    # [1] 23.13211

    pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.densities.village.all_villages.pdf"))
        plotDensities(voomNoNormDup, group=yVillage$samples$batch)
        plotDensities(voomNoNormDup, group=yVillage$samples$Island)
    dev.off()
    save(voomNoNormDup, file=paste0(outputdir, "voomNoNorm.tmm.filtered.duplicate_corrected.indoRNA.village.all_villages.Rda"))

    # DE testing:
    # the inter-subject correlation is input into the linear model fit
    voomNoNormDupVfit <- lmFit(voomNoNormDup, design, block=yVillage$samples$ind, correlation=dupcorNoneDup$consensus)
    voomNoNormDupVfit <- contrasts.fit(voomNoNormDupVfit, contrasts=contr.matrix)
    voomNoNormDupEfit <- eBayes(voomNoNormDupVfit, robust=T)

    # pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.mean-variance-trend.village.all_villages.pdf"))
    #     plotSA(voomNoNormDupEfit, main="Mean-variance trend elimination with duplicate correction")
    # dev.off()

    # get top genes using toptable
    allDEresults <- list()

    for(i in 1:28){
        allDEresults[[i]] <- topTable(voomNoNormDupEfit, coef=i, n=Inf, sort.by="p")
    }

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))

summary(decideTests(voomNoNormDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))

# for (i in 1:28){
#     write.table(allDEresults[[i]], file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.village.all_villages.", colnames(contr.matrix)[i], ".txt"))
# }


##############################################################
### 3. And now with random subsetting for power reasons... ###
##############################################################

# Let's drop PadiraTana because even I agree that three samples is ridiculous. 
yVillage5 <- yVillage[,-grep("Padira Tana", yVillage$samples$Sampling.Site)]
yVillage5$samples <- droplevels(yVillage5$samples) # drop unused levels

# Since individual seeds were so dubious, I'm running this using a counter - 1000 random draws without replicates in there. 
# God bless the internet:
# https://stackoverflow.com/questions/20507247/r-repeat-function-until-condition-met


deTestingSubsets <- function(testingData) {
    # draw the sample
    toKeep <- by(testingData$samples, testingData$samples$Sampling.Site, function(x) sample(x$samples, 5, replace=F))
    yVillageSub <- testingData[,grepl(paste(unlist(toKeep), collapse= '|'), testingData$samples$samples)]
    cleanSubset <- length(unique(yVillageSub$samples$ID)) == 35
    
    while (!cleanSubset) {
        toKeep <- by(testingData$samples, testingData$samples$Sampling.Site, function(x) sample(x$samples, 5, replace=F))
        yVillageSub <- testingData[,grepl(paste(unlist(toKeep), collapse= '|'), testingData$samples$samples)]
        cleanSubset <- length(unique(yVillageSub$samples$ID)) == 35
    }

    # # A sanity check for debugging - output the list of samples used:
    # print("doing DE testing with the following samples:")
    # print(yVillageSub$samples$ID)

    # Once we are happy with our subset, set up design and contrast matrices:
    design5 <- model.matrix(~0 + yVillageSub$samples$Sampling.Site + yVillageSub$samples$Age + yVillageSub$samples$batch + yVillageSub$samples$RIN + yVillageSub$samples$CD8T + yVillageSub$samples$CD4T + yVillageSub$samples$NK + yVillageSub$samples$Bcell + yVillageSub$samples$Mono + yVillageSub$samples$Gran)
    # rename columns to exclude spaces and unrecognised characters
    colnames(design5)=gsub("yVillageSub\\$samples\\$", "", colnames(design5))
    colnames(design5)=gsub("Sampling.Site", "", colnames(design5))
    colnames(design5)=gsub(" ", "_", colnames(design5))

    contr.matrix5 <- makeContrasts(  ANKvsMDB=Anakalung-Madobag, ANKvsMPI=Anakalung-Mappi, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, ANKvsRIN = Anakalung-Rindi, ANKvsHPM = Anakalung-Hupu_Mada,
        WNGvsMDB=Wunga-Madobag, WNGvsMPI=Wunga-Mappi, WNGvsTLL=Wunga-Taileleu, WNGvsRIN = Wunga-Rindi, WNGvsHPM = Wunga-Hupu_Mada,
        RINvsMDB= Rindi-Madobag, RINvsTLL= Rindi-Taileleu, RINvsMPI= Rindi-Mappi, RINvsHPM = Rindi-Hupu_Mada,
        HPMvsMDB= Hupu_Mada-Madobag, HPMvsTLL= Hupu_Mada-Taileleu, HPMvsMPI= Hupu_Mada-Mappi,
        MDBvsMPI=Madobag-Mappi, MDBvsTLL=Madobag-Taileleu,   
        TLLvsMPI=Taileleu-Mappi, 
        levels=colnames(design5)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

    yVillageSub <- calcNormFactors(yVillageSub, method="TMM")

    # There are not enough reps in here to use blocking, so we'll do it without the blocking:
        voomNoNorm5 <- voom(yVillageSub, design5, normalize.method="none", plot=F) 
        voomNoNormVfit5 <- lmFit(voomNoNorm5, design5)
        voomNoNormVfit5 <- contrasts.fit(voomNoNormVfit5, contrasts=contr.matrix5)
        voomNoNormEfit5 <- eBayes(voomNoNormVfit5, robust=T)

    # get top genes using toptable
    allDEresultsNoDup5 <- list()

    for(i in 1:21){
        allDEresultsNoDup5[[i]] <- topTable(voomNoNormEfit5, coef=i, n=Inf, sort.by="p")
    }

    deTesting <- summary(decideTests(voomNoNormEfit5, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
    return(deTesting)
}

set.seed(110584) # Let's go with this one. 
iterations <- 1000
allDESummaries <- list()
allDESummaries <- replicate(iterations, deTestingSubsets(yVillage5), simplify=FALSE)

# Now some analyses... 
# Every entry is three rows, so first is to sum DE genes for every test:
    
deGenesTables <- ldply(allDESummaries, function(x) colSums(x[c(1,3),]))
summary(deGenesTables) # Yeah this is hard to interpret, so probably worth putting aside for now. What about n = 15 instead, and repeating with those? There's no variability in Madobag, though, so maybe 10 is a good compromise

deGenesPlotting <- melt(deGenesTables)
deGenesPlotting[deGenesPlotting == 0] <- NA

pdf(paste0(edaoutput, "subsampling_densities.pdf"), width=18)
    ggplot(deGenesPlotting, aes(x=value, colour=variable, fill=variable, group=variable)) +
        geom_density() +
        xlim(0,200) + 
        theme_bw() + 
        # facet_wrap(. ~ variable, nrow = 5, ncol = 5) + 
        labs(title="", y="density", x="DE genes") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(colour=F, fill=F)
dev.off()

# OK and now that we know what village is weird, let's just plot the RINs out by village, because this is soooo bizarre.

pdf(paste0(edaoutput, "RIN_by_village.pdf"))
    ggplot(yFilt$samples, aes(x=Sampling.Site, y=RIN, fill=Sampling.Site)) +
        geom_violin(trim=T) +
        geom_boxplot(width=0.05, fill="white") + 
        theme_bw() + 
        labs(title="", y="RIN", x="") + 
        theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(fill=F)
dev.off()

##########################################
### 4. Random subsetting to 10 inds... ###
##########################################

# Let's drop small villages: 
yVillage10 <- yVillage[,-grep("Padira Tana|Rindi|Hupu Mada", yVillage$samples$Sampling.Site)]
yVillage10$samples <- droplevels(yVillage10$samples) # drop unused levels
table(yVillage10$samples$Sampling.Site)
# Anakalung   Madobag     Mappi  Taileleu     Wunga 
#        20        17        21        32        17 

deTestingSubsets10 <- function(testingData) {
    # draw the sample
    toKeep <- by(testingData$samples, testingData$samples$Sampling.Site, function(x) sample(x$samples, 10, replace=F))
    yVillageSub <- testingData[,grepl(paste(unlist(toKeep), collapse= '|'), testingData$samples$samples)]
    cleanSubset <- length(unique(yVillageSub$samples$ID)) == 50
    
    while (!cleanSubset) {
        toKeep <- by(testingData$samples, testingData$samples$Sampling.Site, function(x) sample(x$samples, 10, replace=F))
        yVillageSub <- testingData[,grepl(paste(unlist(toKeep), collapse= '|'), testingData$samples$samples)]
        cleanSubset <- length(unique(yVillageSub$samples$ID)) == 50
    }

    # # A sanity check for debugging - output the list of samples used:
    print("doing DE testing with the following samples:")
    print(yVillageSub$samples$ID)

    # Once we are happy with our subset, set up design and contrast matrices:
    design <- model.matrix(~0 + yVillageSub$samples$Sampling.Site + yVillageSub$samples$Age + yVillageSub$samples$batch + yVillageSub$samples$RIN + yVillageSub$samples$CD8T + yVillageSub$samples$CD4T + yVillageSub$samples$NK + yVillageSub$samples$Bcell + yVillageSub$samples$Mono + yVillageSub$samples$Gran)
    # rename columns to exclude spaces and unrecognised characters
    colnames(design)=gsub("yVillageSub\\$samples\\$", "", colnames(design))
    colnames(design)=gsub("Sampling.Site", "", colnames(design))
    colnames(design)=gsub(" ", "_", colnames(design))

    contr.matrix <- makeContrasts(  ANKvsMDB=Anakalung-Madobag, ANKvsMPI=Anakalung-Mappi, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga,
        WNGvsMDB=Wunga-Madobag, WNGvsMPI=Wunga-Mappi, WNGvsTLL=Wunga-Taileleu,
        MDBvsMPI=Madobag-Mappi, MDBvsTLL=Madobag-Taileleu,   
        TLLvsMPI=Taileleu-Mappi, 
        levels=colnames(design)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

    yVillageSub <- calcNormFactors(yVillageSub, method="TMM")

    # There are not enough reps in here to use blocking, so we'll do it without the blocking:
        voomNoNorm <- voom(yVillageSub, design, normalize.method="none", plot=F) 
        voomNoNormVfit <- lmFit(voomNoNorm, design)
        voomNoNormVfit <- contrasts.fit(voomNoNormVfit, contrasts=contr.matrix)
        voomNoNormEfit <- eBayes(voomNoNormVfit, robust=T)

    # get top genes using toptable
    allDEresultsNoDup <- list()

    for(i in 1:10){
        allDEresultsNoDup[[i]] <- topTable(voomNoNormEfit, coef=i, n=Inf, sort.by="p")
    }

    deTesting <- summary(decideTests(voomNoNormEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
    return(deTesting)
}

set.seed(110584) # Let's go with this one. 
iterations <- 1000
allDESummaries10 <- list()
allDESummaries10 <- replicate(iterations, deTestingSubsets10(yVillage10), simplify=FALSE)

# Now some analyses... 
# Every entry is three rows, so first is to sum DE genes for every test:
    
deGenesTables10 <- ldply(allDESummaries10, function(x) colSums(x[c(1,3),]))
summary(deGenesTables10) # Yeah this is hard to interpret, so probably worth putting aside for now. What about n = 15 instead, and repeating with those? There's no variability in Madobag, though, so maybe 10 is a good compromise

write.table(deGenesTables10, file=paste0(outputdir, "DE_subsampling_10_inds.txt"), quote=F, row.names=F, col.names=T, sep="\t", eol="\n")

#And now... some t.tests:
t.test(deGenesTables10$ANKvsMPI, deGenesTables10$WNGvsMPI)$p.value
# [1] 1.290068e-191

t.test(deGenesTables10$MDBvsMPI, deGenesTables10$TLLvsMPI)$p.value
# [1] 1.727248e-212