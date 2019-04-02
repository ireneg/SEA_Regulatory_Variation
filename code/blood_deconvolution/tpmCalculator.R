# Compute gene length (and GC content) from GTF file
# script taken from: https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
# this was found following the thread: https://bioinformatics.stackexchange.com/questions/2567/how-can-i-calculate-gene-length-for-rpkm-calculation-from-counts-data

library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

# Get gene length and GC content on cluster --------------------------

# computing gene length is performed by getting the "union gene model". In this method, the non-duplicated exons for each gene are simply summed up ("non-duplicated" in that no genomic base is double counted)

# load in GTF and fasta files
GTFfile = "/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf"
FASTAfile = "/data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa"

# Load the annotation and reduce it
GTF <- import.gff(GTFfile)
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
reducedGTF$gene_id <- rep(names(grl), elementNROWS(grl))

# Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

# Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

# Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")


