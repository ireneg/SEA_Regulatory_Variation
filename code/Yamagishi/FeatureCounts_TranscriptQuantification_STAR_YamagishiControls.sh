# created by KSB on 09.13.2017
# quantify the abundance of transcripts from the STAR alignment of Yamagishi reads to the human genome (hg38) using featureCounts

input_dir="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls"
GTF="/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.filteredTranscripts.GencodeBasic_TSL123.gtf"
output_dir="/data/cephfs/punim0586/kbobowik/Sumba/FeatureCounts/Yamagishi/Controls"
module load Subread

# Execute featureCounts with an array script
for file in ${input_dir}/uniquely*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  # -T = threads; -a = GTF file; -t = feature type; -g = attribute type used to group features; -o = output file
  echo featureCounts -T 12 -a $GTF -t exon -g gene_id -o ${output_dir}/${sample}.txt $file
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/FeatureCounts_ArrayTable_YamagishiControls.txt

for file in ${input_dir}/uniquely*.bam; do
	sample=`basename $file Aligned.sortedByCoord.out.bam`
	ID=${sample##*_}
	cat ${output_dir}/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${output_dir}/Filter_GeneLengthCount_${sample}.txt
done
