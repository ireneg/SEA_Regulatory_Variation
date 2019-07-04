### Generate figures and analyses for methylation that are akin to those from figure 4 in the RNA data and main text.
### IGR 2019.06.24

### Last edit: IGR 2019.06.25
### Fully compliant, replicates H. Natri's analyses and adds my own.

### 0. Load dependencies and functions and set input paths ###
### 1. Wrangling it using her code on the repo ###
### 2. Model matrices and limma testing. ### 
### 3. UpsetR plots and merging island and village. ###
### 4. Looking at correlations across sites ###
### 5. Heatmaps of top DM genes between some villages but not others. ###
### 6. Heatmaps, but with beta values instead of m-values. ###
### 7. Plotting siglec7 because I am curious. ###
### 8. Good old plot of pairwise correlations within each village and level etc etc... ###


##############################################################
### 0. Load dependencies and functions and set input paths ###
##############################################################

# Load dependencies:
library(edgeR)
library(openxlsx)
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
library(hexbin)
library(data.table)

# Set paths:
inputdir <- "/data/cephfs/punim0586/igallego/indoRNA/dm_testing/" # on server
covariatedir <- "/data/cephfs/punim0586/igallego/indoRNA/"

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/dm_testing/"
edaoutput <- paste0(outputdir, "eda/")

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

# Load cleaned mvalues from Heini:
# y methylation list object
mval <- data.frame(fread(paste0(inputdir, "m_vals.tsv"))) # data.table makes strange life choices...
rownames(mval) <- mval$V1
mval <- mval[,-1] 
head(mval)

##################################################
### 1. Wrangling it using her code on the repo ###
##################################################

# Read covariate matrix
covariates <- read.xlsx(paste0(covariatedir, "metadata_RNA_Batch123.xlsx"), sheet=1, detectDates <- T)

# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
covariates$Age[which(is.na(covariates$Age) == T)]=45

# Sorting to match the order of samples in the methylation dataframe
samplenames <- colnames(mval)
samplenames <- gsub("\\.methyl1", "", samplenames)
samplenames <- gsub("\\.methyl2", "", samplenames)
samplenames <- gsub("\\.", "-", samplenames)
# samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

covariates <- covariates[match(samplenames, covariates$Sample.ID),]
dim(covariates)

covariates$Island <- gsub("West Papua", "Korowai", covariates$Island)
covariates$Sampling.Site <- gsub("Mappi", "Korowai", covariates$Sampling.Site)
covariates$Island <- as.factor(covariates$Island)
covariates$Sequencing.Batch <- as.factor(covariates$Sequencing.Batch)

# Adding methylation batch
covariates$methyl_id <- colnames(mval)
covariates$methyl_exp_id <- with(covariates, paste0(methyl_id, "-rna", Sequencing.Batch))
covariates$methyl_batch <- as.factor(as.numeric(grepl('methyl1', covariates$methyl_id, ignore.case=T))) # Fixed this, hope it replicates. 

# Adding blood deconvolution estimates
blood <- read.table(paste0(covariatedir, "predictedCellCounts_DeconCell.txt"), sep="\t", as.is=T, header=T)
colnames(blood) <- c("Gran","Bcell","CD4T","CD8T","NK","Mono")
blood$ID <- sapply(strsplit(rownames(blood),"[_.]"), `[`, 1)
blood$ID <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", blood$ID)
head(blood)

covariates[,17:23] <- blood[match(covariates$Sample.ID, blood$ID),1:7]
head(covariates)
covariates$replicate <- duplicated(covariates[,1])
unique(covariates[grep("-", covariates$Sampling.Date),6])
covariates$Sampling.Date <- gsub("2016-03-10","10/03/2016",covariates$Sampling.Date) %>% as.Date(., tryFormats="%d/%m/%Y")
covariates[which(!(covariates[,"Sample.ID"] %in% samplenames)),]

# Updating:
colnames(mval) <- covariates$Sample.ID

############################################
### 2. Model matrices and limma testing. ### 
############################################

# Creating a design matrix for islands
designIsland <- model.matrix(~0 + covariates$Island + covariates$Age + covariates$methyl_batch + covariates$CD8T + covariates$CD4T + covariates$NK + covariates$Bcell + covariates$Mono + covariates$Gran)
colnames(designIsland) <- c("Korowai", "Mentawai", "Sumba", "Age", "Methyl_batch", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# Defining pairwise comparisons
contrastsIsland <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsKOR=Sumba - Korowai, MTWvsKOR=Mentawai - Korowai, levels=colnames(designIsland))
all_comparisonsIsland <-colnames(data.frame(contrasts))

# Fitting the linear model
fitIsland <- lmFit(mval, designIsland)
vfitIsland <- contrasts.fit(fitIsland, contrasts=contrastsIsland)
efitIsland <- eBayes(vfitIsland)

# save(efitIsland, file=paste0(outputdir, "efitIsland.Rda"))

allIslandresults <- list()

for(i in 1:3){
    allIslandresults[[i]] <- topTable(efitIsland, coef=i, n=Inf, sort.by="p")
}

summary(decideTests(efitIsland, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        SMBvsMTW SMBvsMPI MTWvsMPI
# Down       2596    18099    10362
# NotSig   855457   837215   845236
# Up         1351     4090     3806

colSums(summary(decideTests(efitIsland, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))[c(1,3),])
# SMBvsMTW SMBvsMPI MTWvsMPI 
#     3947    22189    14168  # matches Heini's numbers perfectly.

# for (i in 1:3){
#     write.table(allIslandresults[[i]], file=paste0(outputdir,"topTable.filtered.village.", colnames(contrastsIsland)[i], ".txt"))
# }

# Now for villages:

# First, remove samples that have less than ten individuals per village, to keep it consistent with the gene expression testing.
table(covariates$Sampling.Site)
#    Anakalung    Bilarenge    Hupu Mada      Madobag        Korowai  Padira Tana 
#           18            1            5           17           22            3 
# Patiala Bawa        Rindi     Taileleu        Wunga   Wura Homba 
#            1            5           31           16            1 

covarVillage <- covariates[-grep("Bilarenge|Patiala Bawa|Wura Homba|Hupu Mada|Padira Tana|Rindi", covariates$Sampling.Site),]
table(covarVillage$Sampling.Site)
# Anakalung   Korowai   Madobag  Taileleu     Wunga 
#        18        22        17        31        16 

mvalVillage <- mval[,-grep("Bilarenge|Patiala Bawa|Wura Homba|Hupu Mada|Padira Tana|Rindi", covariates$Sampling.Site)]

# Set up design matrix
designVillage <- model.matrix(~0 + covarVillage$Sampling.Site + covarVillage$Age + covarVillage$methyl_batch + covarVillage$CD8T + covarVillage$CD4T + covarVillage$NK + covarVillage$Bcell + covarVillage$Mono + covarVillage$Gran)
# rename columns to exclude spaces and unrecognised characters
colnames(designVillage)=gsub("covarVillage\\$", "", colnames(designVillage))
colnames(designVillage)=gsub("Sampling.Site", "", colnames(designVillage))
colnames(designVillage)=gsub(" ", "_", colnames(designVillage))

# Defining pairwise comparisons
# set up contrast matrix
contrastsVillage <- makeContrasts(ANKvsMDB=Anakalung-Madobag, ANKvsKOR=Anakalung-Korowai, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, MDBvsKOR=Madobag-Korowai, MDBvsTLL=Madobag-Taileleu, WNGvsMDB=Wunga-Madobag, TLLvsKOR=Taileleu-Korowai, WNGvsKOR=Wunga-Korowai, WNGvsTLL=Wunga-Taileleu, levels=colnames(designVillage)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

all_comparisonsVillage <- colnames(data.frame(contrastsVillage))

# Fitting the linear model
fitVillage <- lmFit(mvalVillage, designVillage)
vfitVillage <- contrasts.fit(fitVillage, contrasts=contrastsVillage)
efitVillage <- eBayes(vfitVillage)

# save(efitVillage, file=paste0(outputdir, "efitVillage.Rda"))
# load (paste0(outputdir, "efitVillage.Rda"))

# get top genes using toptable
allVillageResults <- list()

for(i in 1:10){
    allVillageResults[[i]] <- topTable(efitVillage, coef=i, n=Inf, sort.by="p")
}

summary(decideTests(efitVillage, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#        ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL
# Down        623    15331     1190       11     4887       38     2515     6959    18884     4962
# NotSig   858090   840741   856949   859386   852122   859348   854799   849773   834847   850776
# Up          691     3332     1265        7     2395       18     2090     2672     5673     3666

colSums(summary(decideTests(efitVillage, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))[c(1,3),])
# ANKvsMDB ANKvsKOR ANKvsTLL ANKvsWNG MDBvsKOR MDBvsTLL WNGvsMDB TLLvsKOR WNGvsKOR WNGvsTLL 
#     1314    18663     2455       18     7282       56     4605     9631    24557     8628 

# for (i in 1:10){
#     write.table(allVillageResults[[i]], file=paste0(outputdir,"topTable.filtered.village.", colnames(contrastsVillage)[i], ".txt"))
# }

#######################################################
### 3. UpsetR plots and merging island and village. ###
#######################################################

byVillages05 <- decideTests(efitVillage, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5)
byIslands05 <- decideTests(efitIsland, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5)

#Now order them both by the same so upsetR can judge them both:
table(rownames(byIslands05) == rownames(byVillages05))

allTogether05 <- data.frame(byVillages05, byIslands05)

# Look at which genes are in common using UpsetR

# Making the final figure:
pdf(paste0(edaoutput, "UpsetR_SamplingSiteComparison_all_levels_fc05_testers.pdf"), width=12)
    # Sort by number of intersects without the island level comparisons, out of curiosity, include all pairwise comparisons
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=20,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 1.5, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T) 
    upset(as.data.frame(abs(allTogether05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsKOR", "WNGvsKOR", "MDBvsKOR", "TLLvsKOR"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_kor, 2), rep(mtw_kor, 2)), nintersects=30,  order.by = "freq", keep.order=T, number.angles = 30, point.size = 3.5, line.size = 1.5, sets.x.label="DEG", mainbar.y.label="DEG in common", text.scale=c(1.6, 1.6, 1.6, 1.6, 1.6, 1.3), decreasing=T)
dev.off()

# And now let's make a three way Venn diagram instead:

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


###############################################
### 4. Looking at correlations across sites ###
###############################################

# Let's look at signal across some of the genes that are DE between WNG and KOR, and between SMB and KOR and ANK and KOR, to check what's going there with the villages

# Sumba vs Korowai:
    ankKor <- allVillageResults[[2]]
    ankKor$probe <- rownames(ankKor)
    wngKor <- allVillageResults[[9]]
    wngKor$probe <- rownames(wngKor)

    smbVillageKor <- merge(ankKor, wngKor, by.x="probe", by.y="probe", suffixes=c(".ank", ".wng"))
    smbVillageKor <- smbVillageKor[order(smbVillageKor$adj.P.Val.wng),]

    wngOnly <- smbVillageKor[smbVillageKor$adj.P.Val.wng <= 0.01 & smbVillageKor$adj.P.Val.ank > 0.01 & abs(smbVillageKor$logFC.wng) >= 0.5,] # was forgetting to filter by logFC here!
    wngOnly <- wngOnly[order(wngOnly$adj.P.Val.wng),] # Too lazy to order inside function...

    ankOnly <- smbVillageKor[smbVillageKor$adj.P.Val.ank <= 0.01 & smbVillageKor$adj.P.Val.wng > 0.01 & abs(smbVillageKor$logFC.ank) >= 0.5,]
    ankOnly <- ankOnly[order(ankOnly$adj.P.Val.ank),]

# Mentawai Korowai:
    tllKor <- allVillageResults[[8]]
    tllKor$probe <- rownames(tllKor)
    mdbKor <- allVillageResults[[5]]
    mdbKor$probe <- rownames(mdbKor)

    mtwVillageKor <- merge(tllKor, mdbKor, by.x="probe", by.y="probe", suffixes=c(".tll", ".mdb"))
    mtwVillageKor <- mtwVillageKor[order(mtwVillageKor$adj.P.Val.mdb),]

    mdbOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.mdb <= 0.01 & mtwVillageKor$adj.P.Val.tll > 0.01 & abs(mtwVillageKor$logFC.mdb) >= 0.5,] # was forgetting to filter by logFC here!
    mdbOnly <- mdbOnly[order(mdbOnly$adj.P.Val.mdb),] # Too lazy to order inside function...

    tllOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.tll <= 0.01 & mtwVillageKor$adj.P.Val.mdb > 0.01 & abs(mtwVillageKor$logFC.tll) >= 0.5,]
    tllOnly <- tllOnly[order(tllOnly$adj.P.Val.tll),]

# What's the rank correlation across the two villages in each island? does it get worse as you go down quintiles/deciles?
smbKorIsland <- allIslandresults[[2]]
smbKorIsland$probe <- rownames(smbKorIsland)
mtwKorIsland <- allIslandresults[[3]]
mtwKorIsland$probe <- rownames(mtwKorIsland)

# Merge island info with village info, pt 2:
smbAllKor <- merge(smbVillageKor, smbKorIsland, by.x="probe", by.y="probe")
smbAllKor <- smbAllKor[order(smbAllKor$adj.P.Val),]

mtwAllKor <- merge(mtwVillageKor, mtwKorIsland, by.x="probe", by.y="probe")
mtwAllKor <- mtwAllKor[order(mtwAllKor$adj.P.Val),]

# And now, by quintiles/deciles, determined on the island-wide p-value? mean expression?:
mtwAllKor$decile <- ntile(mtwAllKor$adj.P.Val, 10)
smbAllKor$decile <- ntile(smbAllKor$adj.P.Val, 10)

# By logFC, is that better?
# Condition on DE (look at top quintile etc), consider log FC direction in general, 
# But there's going to be sooooooooo many ties at the bottom... so by decile or quintile makes a lot more sense here.

cor(smbAllKor[,c(2,8,14)], method="spearman")
#           logFC.ank logFC.wng     logFC
# logFC.ank 1.0000000 0.7717367 0.8930302
# logFC.wng 0.7717367 1.0000000 0.9125966
# logFC     0.8930302 0.9125966 1.0000000
cor(mtwAllKor[,c(2,8,14)], method="spearman")
#           logFC.tll logFC.mdb     logFC
# logFC.tll 1.0000000 0.7717367 0.6525039
# logFC.mdb 0.7717367 1.0000000 0.6112676
# logFC     0.6525039 0.6112676 1.0000000

by(mtwAllKor, mtwAllKor$decile, function(x) cor(x[,c(2,8,14)], method="spearman")) # These are excellent in the top deciles
by(smbAllKor, smbAllKor$decile, function(x) cor(x[,c(2,8,14)], method="spearman")) # Yes indeed - tho this is better than MTW

# Labels for the plot
decileLabels <- c("0-10% inter-island\np-value", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100% inter-island\np-value")
names(decileLabels) <- seq(1,10,1)

    pdf(file=paste0(edaoutput, "village_correlations_by_decile.pdf"), width=18, height=4)
        ggplot(mtwAllKor, aes(y=logFC.tll, x=logFC.mdb, group=decile)) +
            geom_hline(yintercept=0, linetype=3, colour="grey30") +
            geom_vline(xintercept=0, linetype=3, colour="grey30") +
            geom_hline(yintercept=-0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=-0.5, linetype=2, colour="grey60") +
            geom_hline(yintercept=0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=0.5, linetype=2, colour="grey60") +
            geom_hex(bins=50) + 
            theme_bw() + 
            labs(title="", y="log FC TLLvsKOR", x="log FC MDBvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ decile, labeller = labeller(decile = decileLabels))

        ggplot(smbAllKor, aes(y=logFC.wng, x=logFC.ank, group=decile)) +
            geom_hline(yintercept=0, linetype=3, colour="grey30") +
            geom_vline(xintercept=0, linetype=3, colour="grey30") +
            geom_hline(yintercept=-0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=-0.5, linetype=2, colour="grey60") +
            geom_hline(yintercept=0.5, linetype=2, colour="grey60") +
            geom_vline(xintercept=0.5, linetype=2, colour="grey60") +
            geom_hex(bins=50) + 
            theme_bw() + 
            labs(title="", y="log FC WNGvsKOR", x="log FC ANKvsKOR") + 
            theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            guides(fill=F) +
            coord_equal(ratio=1) +
            facet_grid(. ~ decile, labeller = labeller(decile = decileLabels))
    dev.off()    


#########################################################################
### 5. Heatmaps of top DE genes between some villages but not others. ###
#########################################################################

# Instead of top 100 genes let's do top 1000 genes.

# Sumba vs Korowai
    # Remove Mentawai
    wngKorMval <- mvalVillage[wngOnly$probe[1:1000],]
    wngKorMval <- wngKorMval[,grepl("MPI|SMB", colnames(wngKorMval))]

    ankKorMval <- mvalVillage[ankOnly$probe[1:1000],]
    ankKorMval <- ankKorMval[,grepl("MPI|SMB", colnames(ankKorMval))]

    smbKorMval <- mvalVillage[smbKorIsland$probe[1:1000],]
    smbKorMval <- smbKorMval[,grepl("MPI|SMB", colnames(smbKorMval))]

    # Column annotation - same for all plots
    colMetadata <- covarVillage[,c(1,2,3)]
    islandCols <- c("Korowai" = korowai, "Sumba" = sumba, "Mentawai" = mentawai)
    villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod")

    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|SMB", colMetadata$Sample),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Rows: 
    # Wunga-centric:
    wngOnly <- merge(wngOnly, smbKorIsland, by.x="probe", by.y="probe", all=F, sort=F)
    rowMetadataWng <- data.frame(rownames(wngKorMval), -log10(wngOnly$adj.P.Val.wng[1:1000]), -log10(wngOnly$adj.P.Val.ank[1:1000]), -log10(wngOnly$adj.P.Val[1:1000]))
    names(rowMetadataWng) <- c("probe", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsWng <- HeatmapAnnotation(rowMetadataWng[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataWng[,2:4]), max(rowMetadataWng[,2:4])), c("white", "black"))), which="row")

    # Anakalung
    ankOnly <- merge(ankOnly, smbKorIsland, by.x="probe", by.y="probe", all=F, sort=F)
    rowMetadataAnk <- data.frame(rownames(ankKorMval), -log10(ankOnly$adj.P.Val.wng[1:1000]), -log10(ankOnly$adj.P.Val.ank[1:1000]), -log10(ankOnly$adj.P.Val[1:1000]))
    names(rowMetadataAnk) <- c("probe", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsAnk <- HeatmapAnnotation(rowMetadataAnk[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataAnk[,2:4]), max(rowMetadataAnk[,2:4])), c("white", "black"))), which="row")

    # For the island-wide one
    rowMetadataSmb <- data.frame(rownames(smbKorMval), -log10(smbAllKor$adj.P.Val.wng[1:1000]), -log10(smbAllKor$adj.P.Val.ank[1:1000]), -log10(smbAllKor$adj.P.Val[1:1000]))
    names(rowMetadataSmb) <- c("probe", "wngKorpval", "ankKorpval", "smbKorpval")
    rowColsSmb <- HeatmapAnnotation(rowMetadataSmb[,2:4], col = 
        list(wngKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black")), 
             ankKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black")), 
             smbKorpval=colorRamp2(c(min(rowMetadataSmb[,2:4]), max(rowMetadataSmb[,2:4])), c("white", "black"))), which="row")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "smb_kor_DE_heatmaps.pdf"), height=9, width=6)
        wngMap <- Heatmap(t(scale(t(wngKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(wngKorMval)), cluster_columns=F)
        draw(rowColsWng + wngMap, row_dend_side = "left")

        ankMap <- Heatmap(t(scale(t(ankKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(ankKorMval)), cluster_columns=F)
        draw(rowColsAnk + ankMap, row_dend_side = "left")

        smbMap <- Heatmap(t(scale(t(smbKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(smbKorMval)), cluster_columns=F)
        draw(rowColsSmb + smbMap, row_dend_side = "left")
    dev.off()

### Now Mentawai vs Korowai:
    # Remove Sumba
    tllKorMval <- mvalVillage[tllOnly$probe[1:1000],]
    tllKorMval <- tllKorMval[,grepl("MPI|MTW", colnames(tllKorMval))]

    mdbKorMval <- mvalVillage[mdbOnly$probe[1:1000],]
    mdbKorMval <- mdbKorMval[,grepl("MPI|MTW", colnames(mdbKorMval))]

    mtwKorMval <- mvalVillage[mtwKorIsland$probe[1:1000],]
    mtwKorMval <- mtwKorMval[,grepl("MPI|MTW", colnames(mtwKorMval))]

    # Column annotation - same for all plots
    colMetadata <- covarVillage[,c(1,2,3)]
    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|MTW", colMetadata$Sample),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Rows: 
    # Taileleu-centric:
    tllOnly <- merge(tllOnly, mtwKorIsland, by.x="probe", by.y="probe", all=F, sort=F)
    rowMetadataTll <- data.frame(rownames(tllKorMval), -log10(tllOnly$adj.P.Val.tll[1:1000]), -log10(tllOnly$adj.P.Val.mdb[1:1000]), -log10(tllOnly$adj.P.Val[1:1000]))
    names(rowMetadataTll) <- c("probe", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsTll <- HeatmapAnnotation(rowMetadataTll[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataTll[,2:4]), max(rowMetadataTll[,2:4])), c("white", "black"))), which="row")

    # Anakalung
    mdbOnly <- merge(mdbOnly, mtwKorIsland, by.x="probe", by.y="probe", all=F, sort=F)
    rowMetadataMdb <- data.frame(rownames(mdbKorMval), -log10(mdbOnly$adj.P.Val.tll[1:1000]), -log10(mdbOnly$adj.P.Val.mdb[1:1000]), -log10(mdbOnly$adj.P.Val[1:1000]))
    names(rowMetadataMdb) <- c("probe", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsMdb <- HeatmapAnnotation(rowMetadataMdb[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataMdb[,2:4]), max(rowMetadataMdb[,2:4])), c("white", "black"))), which="row")

    # For the island-wide one
    rowMetadataMtw <- data.frame(rownames(mtwKorMval), -log10(mtwAllKor$adj.P.Val.tll[1:1000]), -log10(mtwAllKor$adj.P.Val.mdb[1:1000]), -log10(mtwAllKor$adj.P.Val[1:1000]))
    names(rowMetadataMtw) <- c("probe", "tllKorpval", "mdbKorpval", "mtwKorpval")
    rowColsMtw <- HeatmapAnnotation(rowMetadataMtw[,2:4], col = 
        list(tllKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black")), 
             mdbKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black")), 
             mtwKorpval=colorRamp2(c(min(rowMetadataMtw[,2:4]), max(rowMetadataMtw[,2:4])), c("white", "black"))), which="row")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "mtw_kor_DE_heatmaps.pdf"), height=8, width=6)
        tllMap <- Heatmap(t(scale(t(tllKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(tllKorMval)), cluster_columns=F)
        draw(rowColsTll + tllMap, row_dend_side = "left")

        mdbMap <- Heatmap(t(scale(t(mdbKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mdbKorMval)), cluster_columns=F)
        draw(rowColsMdb + mdbMap, row_dend_side = "left")

        mtwMap <- Heatmap(t(scale(t(mtwKorMval))), col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="expression Z-score", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mtwKorMval)), cluster_columns=F)
        draw(rowColsMtw + mtwMap, row_dend_side = "left")
    dev.off()


##############################################################
### 6. Heatmaps, but with beta values instead of m-values. ###
##############################################################

rm(mvalVillage) # clean up a bit

beta <- data.frame(fread(paste0(inputdir, "beta_vals.tsv")))
rownames(beta) <- beta$V1
beta <- beta[,-1] 
head(beta)

colnames(beta) <- covariates$Sample.ID
betaVillage <- beta[,-grep("Bilarenge|Patiala Bawa|Wura Homba|Hupu Mada|Padira Tana|Rindi", covariates$Sampling.Site)]

# Now the heatmaps - just need to pull the new beta values and replot, the row and colour info shouldn't change. 

# Sumba vs Korowai
    # Remove Mentawai
    wngKorbeta <- betaVillage[wngOnly$probe[1:100],]
    wngKorbeta <- wngKorbeta[,grepl("MPI|SMB", colnames(wngKorbeta))]

    ankKorbeta <- betaVillage[ankOnly$probe[1:100],]
    ankKorbeta <- ankKorbeta[,grepl("MPI|SMB", colnames(ankKorbeta))]

    smbKorbeta <- betaVillage[smbKorIsland$probe[1:100],]
    smbKorbeta <- smbKorbeta[,grepl("MPI|SMB", colnames(smbKorbeta))]

    # The metadata is the same as before, so no need to redefine it beyond this...
    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|SMB", colMetadata$Sample),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "smb_kor_DE_heatmaps_betas.pdf"), height=9, width=6)
        wngMap <- Heatmap(wngKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(wngKorbeta)), cluster_columns=F)
        draw(rowColsWng + wngMap, row_dend_side = "left")

        ankMap <- Heatmap(ankKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(ankKorbeta)), cluster_columns=F)
        draw(rowColsAnk + ankMap, row_dend_side = "left")

        smbMap <- Heatmap(smbKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(smbKorbeta)), cluster_columns=F)
        draw(rowColsSmb + smbMap, row_dend_side = "left")
    dev.off()

### Now Mentawai vs Korowai:
    # Remove Sumba
    tllKorbeta <- betaVillage[tllOnly$probe[1:100],]
    tllKorbeta <- tllKorbeta[,grepl("MPI|MTW", colnames(tllKorbeta))]

    mdbKorbeta <- betaVillage[mdbOnly$probe[1:100],]
    mdbKorbeta <- mdbKorbeta[,grepl("MPI|MTW", colnames(mdbKorbeta))]

    mtwKorbeta <- betaVillage[mtwKorIsland$probe[1:100],]
    mtwKorbeta <- mtwKorbeta[,grepl("MPI|MTW", colnames(mtwKorbeta))]

    # Same here...
    colCols <- HeatmapAnnotation(df = colMetadata[grepl("MPI|MTW", colMetadata$Sample),2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

    # Heatmaps:
    pdf(file=paste0(edaoutput, "mtw_kor_DE_heatmaps_betas.pdf"), height=8, width=6)
        tllMap <- Heatmap(tllKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value", show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(tllKorbeta)), cluster_columns=F)
        draw(rowColsTll + tllMap, row_dend_side = "left")

        mdbMap <- Heatmap(mdbKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value",  show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mdbKorbeta)), cluster_columns=F)
        draw(rowColsMdb + mdbMap, row_dend_side = "left")

        mtwMap <- Heatmap(mtwKorbeta, col=colorRampPalette(brewer.pal(9, "PuOr"))(100), name="beta value",  show_column_names=FALSE, show_row_names=FALSE, top_annotation = colCols, column_order = sort(colnames(mtwKorbeta)), cluster_columns=F)
        draw(rowColsMtw + mtwMap, row_dend_side = "left")
    dev.off()

rm(betaVillage)

#################################################
### 7. Plotting siglec7 because I am curious. ###
#################################################

# read in the manifest:
manifest <- fread(paste0(inputdir, "MethylationEPIC_v-1-0_B4.csv"), sep=",", skip=7)
siglec7 <- manifest[grepl("SIGLEC7", manifest$GencodeBasicV12_NAME),]
siglec7Alt <- manifest[grepl("SIGLEC7", manifest$UCSC_RefGene_Name),]

siglecBetas <- beta[siglec7$Name,]
siglecAltBetas <- beta[siglec7Alt$Name,]

#get missing pvals
smbMtwIsland <- allIslandresults[[1]]
smbMtwIsland$probe <- rownames(smbMtwIsland)

siglecPvals <- merge(smbMtwIsland[siglec7$Name,], smbKorIsland[siglec7$Name,], by.x="probe", by.y="probe", suffixes=c(".smbmtw", ".smbkor"))
siglecPvals <- merge(siglecPvals, mtwKorIsland[siglec7$Name,], by.x="probe", by.y="probe")

siglecAltPvals <- merge(smbMtwIsland[siglec7Alt$Name,], smbKorIsland[siglec7Alt$Name,], by.x="probe", by.y="probe", suffixes=c(".smbmtw", ".smbkor"))
siglecAltPvals <- merge(siglecAltPvals, mtwKorIsland[siglec7Alt$Name,], by.x="probe", by.y="probe")


# Column annotation
colMetadata <- covariates[,c(1,2,3)]
islandCols <- c("Korowai" = korowai, "Sumba" = sumba, "Mentawai" = mentawai)
villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod", "Bilarenge" = "gold4", "Hupu Mada" = "gold4", "Padira Tana" = "gold4", "Patiala Bawa" = "gold4", "Wura Homba" = "gold4", "Rindi" = "gold4")

colCols <- HeatmapAnnotation(df = colMetadata[,2:3], col = list(Island = islandCols, Sampling.Site = villageCols), which="col")

# Read in CPM data to make a nice extra info column
# coming soon

rowMetadataSiglec <- data.frame(rownames(siglecPvals), -log10(siglecPvals[,6]), -log10(siglecPvals[,12]), -log10(siglecPvals[,18]))
names(rowMetadataSiglec) <- c("probe", "SMB_MTW_pval", "SMB_KOR_pval", "MTW_KOR_pval")
rowColsSiglec <- HeatmapAnnotation(rowMetadataSiglec[,2:4], col = 
    list(SMB_MTW_pval=colorRamp2(c(min(rowMetadataSiglec[,2:4]), max(rowMetadataSiglec[,2:4])), c("white", "black")), 
         SMB_KOR_pval=colorRamp2(c(min(rowMetadataSiglec[,2:4]), max(rowMetadataSiglec[,2:4])), c("white", "black")), 
         MTW_KOR_pval=colorRamp2(c(min(rowMetadataSiglec[,2:4]), max(rowMetadataSiglec[,2:4])), c("white", "black"))), which="row")

rowMetadatasiglecAlt <- data.frame(rownames(siglecAltPvals), -log10(siglecAltPvals[,6]), -log10(siglecAltPvals[,12]), -log10(siglecAltPvals[,18]))
names(rowMetadatasiglecAlt) <- c("probe", "SMB_MTW_pval", "SMB_KOR_pval", "MTW_KOR_pval")
rowColssiglecAlt <- HeatmapAnnotation(rowMetadatasiglecAlt[,2:4], col = 
    list(SMB_MTW_pval=colorRamp2(c(min(rowMetadatasiglecAlt[,2:4]), max(rowMetadatasiglecAlt[,2:4])), c("white", "black")), 
         SMB_KOR_pval=colorRamp2(c(min(rowMetadatasiglecAlt[,2:4]), max(rowMetadatasiglecAlt[,2:4])), c("white", "black")), 
         MTW_KOR_pval=colorRamp2(c(min(rowMetadatasiglecAlt[,2:4]), max(rowMetadatasiglecAlt[,2:4])), c("white", "black"))), which="row")

# Sort the data... 
siglecBetas <- siglecBetas[sort(rownames(siglecBetas)),]
siglecAltBetas <- siglecAltBetas[sort(rownames(siglecAltBetas)),]

pdf(file=paste0(edaoutput, "siglec7_heatmaps_betas.pdf"))
    siglecMap <- Heatmap(siglecBetas, col=colorRampPalette(brewer.pal(9, "PuOr"))(130), name="beta value", show_column_names=FALSE, show_row_names=TRUE, top_annotation = colCols, cluster_columns=T)
    draw(rowColsSiglec + siglecMap, row_dend_side = "left")
    siglecAltMap <- Heatmap(siglecAltBetas, col=colorRampPalette(brewer.pal(9, "PuOr"))(130), name="beta value", show_column_names=FALSE, show_row_names=TRUE, top_annotation = colCols, cluster_columns=T)
    draw(rowColssiglecAlt + siglecAltMap, row_dend_side = "left")
dev.off()

# OK that's done and it is ugly and underwhelming


##########################################################################################
### 8. Good old plot of pairwise correlations within each village and level etc etc... ###
##########################################################################################

# This is so slow that it absolutely goes last, and we clean up a bit first
rm(beta)
rm(manifest)

# Define hideous functions
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
                } else if (metadata$methyl_batch[i] == metadata$methyl_batch[j]){
                    if (metadata$Sampling.Site[i] == metadata$Sampling.Site[j]){
                        villageBatchRep <- c(villageBatchRep, corMat[i,j])
                    } else if (metadata$Island[i] == metadata$Island[j]){
                        islandBatchRep <- c(islandBatchRep, corMat[i,j])
                    } else
                        batchRep <- c(batchRep, corMat[i,j])
                } else if (metadata$methyl_batch[i] != metadata$methyl_batch[j]){
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
                } else if (metadata$methyl_batch[i] == metadata$methyl_batch[j]){
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
                } else if (metadata$methyl_batch[i] != metadata$methyl_batch[j]){
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
    plot.reproducibility(mval, covariates, "spearman")
    plot.reproducibility(mval, covariates, "pearson")
dev.off()

pdf(file=paste0(edaoutput, "correlation_within_sites.pdf"))
    plotWithinSite(mval, covariates, "spearman")
    plotWithinSite(mval, covariates, "pearson")
dev.off()


