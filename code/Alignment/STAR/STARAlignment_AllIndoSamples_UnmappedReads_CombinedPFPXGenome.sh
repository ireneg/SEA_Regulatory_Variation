# created by KSB, 29.03.2019
# use the genome aligner STAR to map all Indonesian samples and controls

# load modules
module load STAR

genomeDir=/data/cephfs/punim0586/kbobowik/genomes
mainDir=/data/cephfs/punim0586/kbobowik/
star="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass"

# combine PFalciparum and PVivax genomes
cat ${genomeDir}/PlasmoDB_36_Pfalciparum3D7_Genome.fasta ${genomeDir}/PlasmoDB_36_PvivaxP01_Genome.fasta > ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta

# Generate genome indices
STAR --runMode genomeGenerate --genomeDir ${mainDir}/STAR/IndoRNA_CombinedPFPX --genomeFastaFiles ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta --sjdbGTFfile ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff --sjdbOverhang 100 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# see what maximum intron value is for mapping
cut -f 4,5 ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff | awk '{a=$2-$1;print a;}' | awk 'BEGIN{a=0}{if ($1>0+a) a=$1} END{print a}'
# 34385

# First Pass- align all human reads from the trimmomatic output to the Plasmodium falciparum genome
for file in ${star}/{sample_output,second_batch/sample_output,third_batch/sample_output}/R1_*.fastq; do
	readDir=`dirname $file`
	batch=`echo $file | cut -d/ -f 8`
	# the first batch of sequences isn't labelled as 'first batch' (in the 7th field separator) so we need to add this information in. 
	if [[ $batch == "sample_output" ]]; then
		batch="first_batch"
	fi
	f1=`basename $file`
	f2=R2${f1#"R1"}
	shortenedFile=`basename $f1 Aligned.sortedByCoord.out.fastq`
	sampleID=${shortenedFile##*_}
	# --alignMatesGapMax = maximum genomic distance between mates. Note that for PF, the maximum annotated intron in the GTF is ~2.5kb, so it makes sense to limit the max gap to 10-20kb.
	# --genomeSAindexNbases = For small genomes, the parameter--genomeSAindexNbases needs to be scaled down, with a typicalvalue ofmin(14, log2(GenomeLength)/2 - 1).  For example, for 1 megaBase genome, this is equalto 9, for 100 kiloBase genome, this is equal to 7.
	# The PX genome is 22MB and the PF genome is 22.9MB. Let's just say the max MB is 23--> log2(23) = 24.5/2 = 12.25 - 1 = 11.25. We'll use 11 as the --genomeSAindexNbases
	echo STAR --genomeDir ${mainDir}/STAR/IndoRNA_CombinedPFPX --readFilesIn ${readDir}/${f1} ${readDir}/${f2} --runThreadN 12 --sjdbOverhang 100 --outFileNamePrefix ${mainDir}/STAR/IndoRNA_CombinedPFPX/IndoOutput/STAR_PXPFCombined_${batch}_${sampleID} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
done > ${mainDir}/Sumba/scripts/Array_Scripts/STARArrayTable_IndoRNA_PXPFCombined.txt
# saved as shell file: STAR_FirstPass_PFPXCombined_indoRNAallBatches.sh

# combine all of the SJ.out.tab files from from all files created in first pass into one file 
cat ${mainDir}/STAR/IndoRNA_CombinedPFPX/IndoOutput/*.out.tab > ${mainDir}/STAR/IndoRNA_CombinedPFPX/IndoOutput/STAR_PFPX_allFilesSJ.out.tab

# generate genomice indices for the second time
STAR --runMode genomeGenerate --genomeDir ${mainDir}/STAR/IndoRNA_CombinedPFPX_SecondPass --genomeFastaFiles ${genomeDir}/combined_PFalc3D7_PvivaxP01_Genome.fasta --sjdbFileChrStartEnd ${mainDir}/STAR/IndoRNA_CombinedPFPX/IndoOutput/STAR_PFPX_allFilesSJ.out.tab --sjdbGTFfile ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff --sjdbOverhang 100 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# Second pass-
for file in ${star}/{sample_output,second_batch/sample_output,third_batch/sample_output}/R1_*.fastq; do
	readDir=`dirname $file`
	batch=`echo $file | cut -d/ -f 8`
	# the first batch of sequences isn't labelled as 'first batch' (in the 7th field separator) so we need to add this information in. 
	if [[ $batch == "sample_output" ]]; then
		batch="first_batch"
	fi
	f1=`basename $file`
	f2=R2${f1#"R1"}
	shortenedFile=`basename $f1 Aligned.sortedByCoord.out.fastq`
	sampleID=${shortenedFile##*_}
	# --alignMatesGapMax = maximum genomic distance between mates. Note that for PF, the maximum annotated intron in the GTF is ~2.5kb, so it makes sense to limit the max gap to 10-20kb.
	# --genomeSAindexNbases = For small genomes, the parameter--genomeSAindexNbases needs to be scaled down, with a typicalvalue ofmin(14, log2(GenomeLength)/2 - 1).  For example, for 1 megaBase genome, this is equalto 9, for 100 kiloBase genome, this is equal to 7.
	# The PX genome is 22MB and the PF genome is 22.9MB. Let's just say the max MB is 23--> log2(23) = 24.5/2 = 12.25 - 1 = 11.25. We'll use 11 as the --genomeSAindexNbases
	echo STAR --genomeDir ${mainDir}/STAR/IndoRNA_CombinedPFPX_SecondPass --readFilesIn ${readDir}/${f1} ${readDir}/${f2} --runThreadN 12 --sjdbOverhang 100 --outFileNamePrefix ${mainDir}/STAR/IndoRNA_CombinedPFPX_SecondPass/IndoOutput/STAR_PXPFCombined_${batch}_${sampleID} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 35000 --alignMatesGapMax 35000
done > ${mainDir}/Sumba/scripts/Array_Scripts/STARArrayTable_SecondPass_IndoRNA_PXPFCombined.txt
# saved as shell file: STAR_SecondPass_PFPXCombined_indoRNAallBatches.sh

# filter the bam file so that only uniquely mapped reads remain
for pathname in /data/cephfs/punim0586/kbobowik/STAR/IndoRNA_CombinedPFPX_SecondPass/IndoOutput/STAR_*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # "-q 255" = unique reads
    samtools view -q 255 $pathname > ${directory}/uniquelyMapped_${filename} 
done
