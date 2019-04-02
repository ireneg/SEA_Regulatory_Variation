# created by KSB, 29.03.2019
# use the genome aligner STAR to map all Yamagishi healthy samples and controls

# load modules
module load STAR

genomeDir=/data/cephfs/punim0586/kbobowik/genomes
mainDir=/data/cephfs/punim0586/kbobowik
star=/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass

# combine PFalciparum and PVivax genomes
cat ${genomeDir}/PlasmoDB_36_Pfalciparum3D7_Genome.fasta ${genomeDir}/PlasmoDB_36_PvivaxP01_Genome.fasta > ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta

# Generate genome indices
# Note: jobs will get killed of you don't ask for enough memory. I found asking for at least 12 CPUs and 200,000 MB was sufficient

# the '--sjdbGTFtagExonParentTranscript Parent' flag makes STAR work with gff files
STAR --runMode genomeGenerate --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX --genomeFastaFiles ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta --sjdbGTFfile ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff --sjdbOverhang 35 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# First Pass - yamagishi
for file in ${star}/Yamagishi/{Controls,Sick}/*.fastq; do
  sample=`basename $file _trimmed.fastq.gz`
  healthStatus=`echo $file | cut -d/ -f 9`
  if [[ $healthStatus == "Controls" ]]; then
  	echo STAR --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX --readFilesIn $file --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix ${mainDir}/STAR/Yamagishi_CombinedPFPX/Controls/STAR_PXPFCombined_${healthStatus}_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
  else
  	echo STAR --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX --readFilesIn $file --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix ${mainDir}/STAR/Yamagishi_CombinedPFPX/Sick/STAR_PXPFCombined_${healthStatus}_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
  fi
done > ${mainDir}/Sumba/scripts/Array_Scripts/STARArrayTable_YamagishiSamples_PXPFCombined.txt
# run as sbatch script 

----
# combine all of the SJ.out.tab files from from all files created in first pass into one file 
cat ${mainDir}/STAR/Yamagishi_CombinedPFPX/{Controls,Sick}/*.out.tab > ${mainDir}/STAR/Yamagishi_CombinedPFPX/Sick/STAR_PFPX_allFilesSJ.out.tab

# generate genomice indices for the second time
STAR --runMode genomeGenerate --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX_SecondPass --genomeFastaFiles ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta --sjdbFileChrStartEnd ${mainDir}/STAR/Yamagishi_CombinedPFPX/Sick/STAR_PFPX_allFilesSJ.out.tab --sjdbGTFfile ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff --sjdbOverhang 35 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# Map the reads for the second pass
for file in ${star}/Yamagishi/{Controls,Sick}/*.fastq; do
  sample=`basename $file _trimmed.fastq.gz`
  healthStatus=`echo $file | cut -d/ -f 9`
  if [[ $healthStatus == "Controls" ]]; then
  	echo STAR --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX_SecondPass --readFilesIn $file --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix ${mainDir}/STAR/Yamagishi_CombinedPFPX_SecondPass/Controls/STAR_PXPFCombined_${healthStatus}_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
  else
  	echo STAR --genomeDir ${mainDir}/STAR/Yamagishi_CombinedPFPX_SecondPass --readFilesIn $file --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix ${mainDir}/STAR/Yamagishi_CombinedPFPX_SecondPass/Sick/STAR_PXPFCombined_${healthStatus}_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
  fi
done > ${mainDir}/Sumba/scripts/Array_Scripts/STARArrayTable_YamagishiSamples_PXPFCombined_SecondPass.txt

---#---
# filter the bam file so that only uniquely mapped reads remain
for pathname in /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls/STAR_*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # "-q 255" = unique reads
    echo samtools view -q 255 $pathname '>' ${directory}/uniquelyMapped_${filename} 
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/uniquelyMappedReads_ArrayTable_YamagishiControls.txt

