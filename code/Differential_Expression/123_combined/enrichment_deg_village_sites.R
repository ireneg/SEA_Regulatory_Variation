### Script from Pai Kusuma, modified a little bit by IGR to detect enrichment in genes DE between some comparisons but not others

### Last update: 2019.07.05 IGR
### Fixed a bug in KEGG testing where the background was not getting set correctly, cleaned up.

### 0. Load the relevant packages etc. ###
### 1. Define the GO/KEGG enrichment testing function. ###
### 2. Testing the island and village-level genes. ###
### 3. Testing the single-village DE genes. ### 
### 4. Annotating for Denisovan introgression and t test of frequencies. ###


##########################################
### 0. Load the relevant packages etc. ###
##########################################

library(tidyverse)
library(stringr)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(tidyr)

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

files <- list.files(path=inputdir, 
                    pattern="topTable.*\\.dup_corrected.*txt", full.names=TRUE, recursive=FALSE)

mtwKor <- read.table(files[1], header=T, stringsAsFactors=F)
smbKor <- read.table(files[3], header=T, stringsAsFactors=F)
smbMtw <- read.table(files[5], header=T, stringsAsFactors=F)

ankKor <- read.table(files[6], header=T, stringsAsFactors=F)
wngKor <- read.table(files[16], header=T, stringsAsFactors=F)
smbVillageKor <- merge(ankKor, wngKor, by.x="genes", by.y="genes", suffixes=c(".ank", ".wng"))
wngOnly <- smbVillageKor[smbVillageKor$adj.P.Val.wng <= 0.01 & smbVillageKor$adj.P.Val.ank > 0.01,]
ankOnly <- smbVillageKor[smbVillageKor$adj.P.Val.ank <= 0.01 & smbVillageKor$adj.P.Val.wng > 0.01,]

tllKor <- read.table(files[14], header=T, stringsAsFactors=F)
mdbKor <- read.table(files[11], header=T, stringsAsFactors=F)
mtwVillageKor <- merge(tllKor, mdbKor, by.x="genes", by.y="genes", suffixes=c(".tll", ".mdb"))
mdbOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.mdb <= 0.01 & mtwVillageKor$adj.P.Val.tll > 0.01,]
tllOnly <- mtwVillageKor[mtwVillageKor$adj.P.Val.tll <= 0.01 & mtwVillageKor$adj.P.Val.mdb > 0.01,]

methylMix <- read.table(paste0(inputdir, "methylmix_genename_ensemblid.csv"), sep=",", header=T)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://grch37.ensembl.org/", verbose=T) # super slow if we use the default mirrors...


##########################################################
### 1. Define the GO/KEGG enrichment testing function. ###
##########################################################

degEnrichment <- function(DEG_pop, DEG_comparison) {
    # outfile_annot <- paste0(edaoutput, DEG_comparison, "annot.txt")
    outfile_go <- paste0(edaoutput, DEG_comparison, "GO.txt")
    outfile_kegg <- paste0(edaoutput, DEG_comparison, "KEGG.txt")

    message("")
    message("####")
    message(paste("Annotating ", DEG_comparison, " ...", sep = ""))

    #listFilters(ensembl) %>% filter(str_detect(name, "ensembl"))
    Type_ensembl <- "ensembl_gene_id"
    Values <- DEG_pop$genes
    attribute_Names = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'gene_biotype')

    annot_DEG_pop <- getBM(attributes = attribute_Names, filters = Type_ensembl, values = Values, mart = ensembl)
    tab_DEG_pop <- left_join(DEG_pop, annot_DEG_pop, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)

    # For debugging
    print(head(tab_DEG_pop))

    ## enrichment
    message(paste("Performing GO and KEGG enrichment for ", DEG_comparison, "...", sep = ""))

    ego_or_ALL_pop <- enrichGO(gene = tab_DEG_pop$genes, universe = background$genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff  = 0.05, readable = TRUE, pool=TRUE)

    write_tsv(filter(ego_or_ALL_pop@result, p.adjust <= 0.05), outfile_go)

    # For KEGG we need to annotate the background too, because otherwise there's no way to get results, since the background cannot be mapped properly!
    # Because this is slow we do it outside the function, once per background set, instead of every time.
    kegg_or_ALL_pop <- enrichKEGG(gene = tab_DEG_pop$entrezgene_id, universe = background$entrezgene_id, organism = "hsa", pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff = 0.05, use_internal_data=F) 
    kegg_readable <- setReadable(kegg_or_ALL_pop, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

    write_tsv(filter(kegg_or_ALL_pop@result, p.adjust <= 0.05), outfile_kegg)

    # I learned to make figures!
    message(paste("Plotting GO and KEGG enrichment for ", DEG_comparison, "...", sep = ""))

    pdf(file=paste0(edaoutput, DEG_comparison, "plots.pdf"), height=9, width=9)
        try(print(cnetplot(ego_or_ALL_pop)))
        try(print(heatplot(ego_or_ALL_pop)))
        try(print(emapplot(ego_or_ALL_pop)))
        try(print(cnetplot(kegg_readable)))
        try(print(heatplot(kegg_readable)))
        try(print(emapplot(kegg_readable)))
    dev.off()

    # For debugging:
    return(kegg_readable)

    message(paste("Done"))
    message("####")
}


######################################################
### 2. Testing the island and village-level genes. ###
######################################################

# Define the background for all of these tests:
background <- ankKor # doesn't matter, just need all tested genes. But because annotation is slow, we do it outside the function:
    annot_background <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand', 'kegg_enzyme'), filters = 'ensembl_gene_id', values = background$genes, mart = ensembl)
    background <- left_join(background, annot_background, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)  
    background$entrezgene_id <- as.character(background$entrezgene_id)

# We'll need these later, so save the output, but first filter to genes with FDR < 0.01
# mtwKor <- mtwKor[mtwKor$adj.P.Val <= 0.01,]
# smbKor <- smbKor[smbKor$adj.P.Val <= 0.01,]
# smbKor <- smbKor[smbKor$adj.P.Val <= 0.01,]

mtwKorKEGG <- degEnrichment(mtwKor[mtwKor$adj.P.Val <= 0.01,], "MTW_vs_KOR.") # 0 / 2 
smbKorKEGG <- degEnrichment(smbKor[smbKor$adj.P.Val <= 0.01,], "SMB_vs_KOR.") # 0 / 2 
smbMtwKEGG <- degEnrichment(smbMtw[smbMtw$adj.P.Val <= 0.01,], "SMB_vs_MTW.") # 0 / 2 

# Village level:
# Don't forget the p-value cutoff if you're going to run these!
# degEnrichment(ankKor, "ANK_vs_KOR")  
# degEnrichment(wngKor, "WNG_vs_KOR")  

# degEnrichment(tllKor, "TLL_vs_KOR")  
# degEnrichment(mdbKor, "MDB_vs_KOR") 

# Now the island-unique ones. No need to redefine background for this first set:
degEnrichment(ankOnly, "ankOnly.") #   0 /  1
degEnrichment(wngOnly, "wngOnly.") #  30 /  8
degEnrichment(tllOnly, "tllOnly.") # 217 / 30
degEnrichment(mdbOnly, "mdbOnly.") #   4 /  0


### Finally, enrichment of methylMix genes vs all genes - would be nice to test genes against particular comparisons, but it doesn't seem to have that resolution. 
names(methylMix)[2] <- "genes"
methylMix <- methylMix[!duplicated(methylMix$genes),]
dim(methylMix)
# [1] 1261    3 # A lot more than if we do it by HGNC names, which is worrisome.

methylMixKEGG <- degEnrichment(methylMix, "methylMix.all.") # 0 / 2 

### Now we gotta integrate DE with methylMix output


###############################################
### 3. Testing the single-village DE genes. ### 
###############################################

# Now we do redefine the background although it's unclear what the merits are:
# NB this is my ideal test case for the bug, since the background would be bigger than the list of tested genes for the GO stuff unless I have fixed it properly.
# Edit: I was correct - the old version tested against all known KEGG annotations, the background annotation code below tests against only the background set.  

# This version considers all genes that are DE in either village against Korowai. This is the best compromise, because the island level testing is missing some genes that are DE at a single village level. Could also go crazy and incorporate those, but don't see the need.
background <- smbVillageKor[smbVillageKor$adj.P.Val.wng <= 0.01 | smbVillageKor$adj.P.Val.ank <= 0.01,]
    annot_background <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand', 'kegg_enzyme'), filters = 'ensembl_gene_id', values = background$genes, mart = ensembl)
    background <- left_join(background, annot_background, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)  
    background$entrezgene_id <- as.character(background$entrezgene_id)

degEnrichment(ankOnly, "ankOnly.islandBG.") # 0 / 3
degEnrichment(wngOnly, "wngOnly.islandBG.") # 0 / 8

background <- mtwVillageKor[mtwVillageKor$adj.P.Val.mdb <= 0.01 | mtwVillageKor$adj.P.Val.tll <= 0.01,]
    annot_background <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand', 'kegg_enzyme'), filters = 'ensembl_gene_id', values = background$genes, mart = ensembl)
    background <- left_join(background, annot_background, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)  
    background$entrezgene_id <- as.character(background$entrezgene_id)

degEnrichment(tllOnly, "tllOnly.islandBG.") # 98 / 30
degEnrichment(mdbOnly, "mdbOnly.islandBG.") #  0 /  0


############################################################################
### 4. Annotating for Denisovan introgression and t test of frequencies. ###
############################################################################

## Denisovan introgression and any excess thereof:
denisovan <- read.table(paste0(covariatedir, "denisovan_genes_jacobs_et_al.txt"), sep = "\t", header = T, quote = NULL)
# But of course that file is a mess because every row can have one gene or many, or none, so we use tidyr. 
denisovan2 <- separate_rows(denisovan, GENES, sep="\\,") #OMG 

# Executive decision: if a gene is seen multiple times keep the highest frequency
# In which case we can sort by frequency, then filter duplicate gene names:
denisovan2 <- denisovan2[order(-denisovan2$PROPORTION_NG),]
denisovanFilt <- denisovan2[!duplicated(denisovan2$GENES),]
dim(denisovanFilt)
# [1] 3224    6

denisovanAnnot <- function(DEG_pop, DEG_name) {
    outfile_annot <- paste0(edaoutput, DEG_name, "annot.txt")
    message("")
    message("####")
    message(paste("Annotating ", DEG_name, " ...", sep = ""))

    Type_ensembl <- "ensembl_gene_id"
    Values <- DEG_pop$genes
    attribute_Names = c('ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')

    annot_DEG_pop <- getBM(attributes = attribute_Names, filters = Type_ensembl, values = Values, mart = ensembl)
    tab_DEG_pop <- left_join(DEG_pop, annot_DEG_pop, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)

    tab_DEG_pop_deni <- denisovanFilt %>% inner_join(tab_DEG_pop, by = c("GENES"="external_gene_name")) 
    ## save annotation
    write_tsv(tab_DEG_pop_deni, outfile_annot)

    return(tab_DEG_pop_deni)
    message(paste("Done"))
    message("####")
}

# This is probably a bit fraught because different haplotypes are different lengths and we need to control for that and gene lengths, but I am being sloppy. 

smbKorAnnot <- denisovanAnnot(smbKor, "Smb-Kor.")
mtwKorAnnot <- denisovanAnnot(mtwKor, "Mtw-Kor.")
smbMtwAnnot <- denisovanAnnot(smbMtw, "Smb-Mtw.")

# Two tailed or one-tailed?
t.test(smbKorAnnot[smbKorAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG, mtwKorAnnot[mtwKorAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG) # two tailed, we think they're equal #NB that it wouldn't be significant even if it was one tailed. 
# data:  smbKorAnnot[smbKorAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG and mtwKorAnnot[mtwKorAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG
# t = -0.06737, df = 1432.7, p-value = 0.9463
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.00971487  0.00906973
# sample estimates:
#  mean of x  mean of y 
# 0.07547917 0.07580174 

t.test(smbKorAnnot[smbKorAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG, smbMtwAnnot[smbMtwAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG, alternative="g") # greater, we think kor should contribute
# data:  smbKorAnnot[smbKorAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG and smbMtwAnnot[smbMtwAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG
# t = 1.8188, df = 618.62, p-value = 0.03471
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#  0.0009350708          Inf
# sample estimates:
#  mean of x  mean of y 
# 0.07547917 0.06556098 

t.test(mtwKorAnnot[mtwKorAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG, smbMtwAnnot[smbMtwAnnot$adj.P.Val <= 0.01,]$PROPORTION_NG, alternative="g") # and here too. 
# data:  mtwKorAnnot[mtwKorAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG and smbMtwAnnot[smbMtwAnnot$adj.P.Val <= 0.01, ]$PROPORTION_NG
# t = 1.8244, df = 659.83, p-value = 0.03427
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#  0.0009948329          Inf
# sample estimates:
#  mean of x  mean of y 
# 0.07580174 0.06556098 

