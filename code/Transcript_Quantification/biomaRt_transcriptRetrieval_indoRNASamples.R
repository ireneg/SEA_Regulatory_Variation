# created on 09.04.2017 by KSB
# create a filtered GTF file for GRCh38.p10 transcript annotation via biomaRt R package

library(biomaRt)

# select for transcripts from protein-coding genes                                                
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='aug2017.archive.ensembl.org')
ensembl.biotype <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'gene_biotype', 'external_gene_name', 'transcript_tsl', 'transcript_gencode_basic'), mart = ensembl.mart.90)

# filter for GENCODE basic annotation and trancsript support levels 1-3
ensembl90transcripts = ensembl.biotype[(ensembl.biotype$transcript_tsl %in% c('tsl1','tsl2','tsl3') & ensembl.biotype$transcript_gencode_basic == 'GENCODE basic'),]

# save whole table
write.table(ensembl90transcripts, file="~/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Ensembl90_Transcriptome/biomaRt_GTF_GencodeBasic_TSL123.txt", sep="\t", col.names=T)

# save transcript names
write.table(unique(ensembl90transcripts$ensembl_transcript_id), file="~/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Ensembl90_Transcriptome/biomaRt_transcriptNames_GTF_GencodeBasic_TSL123.txt", col.names=F, row.names=F, quote=F)

# save gene names
write.table(unique(ensembl90transcripts$ensembl_gene_id), file="~/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Ensembl90_Transcriptome/biomaRt_geneNames_GTF_GencodeBasic_TSL123.txt", col.names=F, row.names=F, quote=F)
