### Script from Pai Kusuma, modified a little bit by IGR to detect enrichment in genes DE between some comparisons but not others
### 

library(tidyverse)
library(stringr)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

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

background <- ankKor # doesn't matter, just need all tested genes.

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://grch37.ensembl.org/", verbose=T) # super slow if we use the default mirrors...

degEnrichment <- function(DEG_pop, DEG_comparison) {
    # outfile_annot <- paste0(edaoutput, DEG_comparison, "annot.txt")
    outfile_go <- paste0(edaoutput, DEG_comparison, "GO.txt")
    outfile_kegg <- paste0(edaoutput, DEG_comparison, "KEGG.txt")

    message("")
    message("####")
    message(paste("Annotating ", DEG_comparison, " ...", sep = ""))

    # DEG_pop <- read.table(infile, header = TRUE)

    #listFilters(ensembl) %>% filter(str_detect(name, "ensembl"))
    Type_ensembl <- "ensembl_gene_id"
    Values <- DEG_pop$genes
    attribute_Names = c('ensembl_gene_id', 'entrezgene', 'external_gene_name', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')

    annot_DEG_pop <- getBM(attributes = attribute_Names, filters = Type_ensembl, values = Values, mart = ensembl)
    tab_DEG_pop <- left_join(DEG_pop, annot_DEG_pop, by = c("genes"="ensembl_gene_id")) %>% distinct(genes, .keep_all = TRUE)

    # ## deni
    # dat <- read.table("~/Dropbox/Work/projects/rna/deni_introg/QT.txt", sep = "\t", header = T, quote = NULL)
    # QT <- dat %>% fill(QUANTILE) %>% distinct(GENES, .keep_all = TRUE)

    # tab_DEG_pop_deni <- tab_DEG_pop %>% left_join(QT, by = c("external_gene_name"="GENES")) %>% dplyr::rename("DENI_intro"="QUANTILE")

    # ## save annotation
    # write_tsv(tab_DEG_pop_deni, outfile_annot)

    ## enrichment
    message(paste("Performing GO and KEGG enrichment for ", DEG_comparison, "...", sep = ""))

    ego_or_ALL_pop <- enrichGO(gene = tab_DEG_pop$genes, universe = background$genes, OrgDb = org.Hs.eg.db, 
                               keyType = "ENSEMBL", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.01, 
                               qvalueCutoff  = 0.05, readable = TRUE, pool=TRUE)

    write_tsv(filter(ego_or_ALL_pop@result, p.adjust <= 0.1), outfile_go)

    kegg_or_ALL_pop <- enrichKEGG(gene = as.character(na.omit(tab_DEG_pop$entrezgene)), 
                                  universe = as.character(na.omit(background$entrezgene)), organism = "hsa", 
                                  pAdjustMethod = "fdr", pvalueCutoff = 0.01, qvalueCutoff = 0.05) %>%
      setReadable(., OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

    write_tsv(filter(kegg_or_ALL_pop@result, p.adjust <= 0.1), outfile_kegg)

    message(paste("Done"))
    message("####")
}

degEnrichment(ankOnly, "ankOnly")
degEnrichment(wngOnly, "wngOnly")
degEnrichment(tllOnly, "tllOnly")
degEnrichment(mdbOnly, "mdbOnly")