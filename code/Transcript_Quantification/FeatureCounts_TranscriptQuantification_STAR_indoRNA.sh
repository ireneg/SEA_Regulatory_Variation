# created by KSB, 09.13.2017

# quantify the abundance of transcripts from the STAR alignment of Indonesian reads using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

#download GTF file from Ensembl and unzip
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
gunzip /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf.gz
 
#intersect filtered transcriptome names, gathered from biomaRt R package
# fgrep -f /vlsci/SG0008/kbobowik/Sumba/Ensembl90_Transcriptome/biomaRt_transcriptNames_GTF.txt /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf > Homo_sapiens.GRCh38.90.filteredTranscripts.gtf
fgrep -f /vlsci/SG0008/kbobowik/Sumba/Ensembl90_Transcriptome/biomaRt_transcriptNames_GTF_GencodeBasic_TSL123.txt /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf > /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf

#check number of transcripts
cut -d ';' -f3 /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf | sort | uniq | wc -l

## 54209 (which means we're missing 4181)

# save transcripts that are in GTF
cut -d ';' -f3 /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf | sort | uniq | cut -d ' ' -f 3 > /vlsci/SG0008/kbobowik/Sumba/Ensembl90_Transcriptome/output.transcriptIDs_GencodeBasic_TSL123.txt

### FeatureCounts pipeline -------------------------

# execute FeatureCounts 
# here are the following flags I used:
# -T: Number of the threads.  The value should be between 1 and 32.  1 by default.
# -s: Indicate if strand-specific read counting should be performed. Acceptable  values:  0  (unstranded),  1  (stranded)  and  2  (re-versely stranded).  0 by default.  For paired-end reads, strand of the first read is taken as the strand of the whole fragment. FLAG field is used to tell if a read is first or second read in a pair.
# -p: If specified, fragments (or templates) will be counted instead of reads.  This option is only applicable for paired-end reads.
# -a: Give the name of an annotation file
# - t: Specify the feature type.  Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting.  ‘exon’ by default.
# -o: Give the name of the output file
# -g: Specify the attribute type used to group features (eg.  exons) into meta-features (eg.  genes)

for file in /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass/sample_output/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  featureCounts -T 6 -s 2 -p -a /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf -t exon -g gene_id -o ~/Sumba/FeatureCounts/indoRNA/sample_counts/${sample}.txt $file &> ~/Sumba/FeatureCounts/indoRNA/sample_counts/${ID}_OutputSummary.txt
  cat ~/Sumba/FeatureCounts/indoRNA/sample_counts/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ~/Sumba/FeatureCounts/indoRNA/sample_counts/Filter_GeneLengthCount_${sample}.txt
done

# make summary file
for file in /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass/sample_output/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  tail -n +2 ~/Sumba/FeatureCounts/indoRNA/sample_counts/${sample}.txt.summary | cut -f2 | tr -s '\n' '\t' |  awk -v prefix="$ID\t" '{print prefix $0}' >> ~/Sumba/FeatureCounts/indoRNA/sample_counts/FeatureCounts_summaryFile.txt
  done
# add header to file
a=`tail -n +2 ~/Sumba/FeatureCounts/indoRNA/sample_counts/${sample}.txt.summary | cut -f1 | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" ~/Sumba/FeatureCounts/indoRNA/sample_counts/FeatureCounts_summaryFile.txt
