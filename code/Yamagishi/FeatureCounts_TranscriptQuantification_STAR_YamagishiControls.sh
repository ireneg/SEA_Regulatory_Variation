# created by KSB on 09.13.2017
# quantify the abundance of transcripts from the STAR alignment of Yamagishi reads to the human genome (hg38) using featureCounts

input_dir="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls"
GTF="/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf"
output_dir="/data/cephfs/punim0586/kbobowik/Sumba/FeatureCounts/Yamagishi/Controls"
module load Subread
module load SAMtools

# execute FeatureCounts 
# here are the following flags I used:
# -T: Number of the threads.  The value should be between 1 and 32.  1 by default.
# -s: Indicate if strand-specific read counting should be performed. Acceptable  values:  0  (unstranded),  1  (stranded)  and  2  (re-versely stranded).  0 by default.  For paired-end reads, strand of the first read is taken as the strand of the whole fragment. FLAG field is used to tell if a read is first or second read in a pair.
# -p: If specified, fragments (or templates) will be counted instead of reads.  This option is only applicable for paired-end reads.
# -a: Give the name of an annotation file
# - t: Specify the feature type.  Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting.  ‘exon’ by default.
# -o: Give the name of the output file
# -g: Specify the attribute type used to group features (eg.  exons) into meta-features (eg.  genes)

# Execute featureCounts with an array script
for file in ${input_dir}/STAR*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  # -T = threads; -a = GTF file; -t = feature type; -g = attribute type used to group features; -o = output file
  echo featureCounts -T 12 -a $GTF -t exon -g gene_id -o ${output_dir}/${sample}.txt $file
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/FeatureCounts_ArrayTable_YamagishiControls.txt

for file in ${input_dir}/STAR*.bam; do
	sample=`basename $file Aligned.sortedByCoord.out.bam`
	ID=${sample##*_}
	cat ${output_dir}/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${output_dir}/Filter_GeneLengthCount_${sample}.txt
done

