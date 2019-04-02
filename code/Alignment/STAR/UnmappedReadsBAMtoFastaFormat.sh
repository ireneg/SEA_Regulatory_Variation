# created by KSB, 14.03.2019
# Get total number of unmapped read pairs from Indonesian samples after STAR alignment

# Load modules -------------------

# load samtools 1.9
module load SAMtools

# specify directory
rootdir=/data/cephfs/punim0586/kbobowik

# Filter reads for only unmapped reads ---------------

# the samtools -f flag outputs alignmemnts only matching this flag. The '4' flag (or 00100) corresponds to all unmapped reads (as per the Samtools manual: http://www.htslib.org/doc/samtools-1.2.html and the Picard 'explain samtools flags' website: http://broadinstitute.github.io/picard/explain-flags.html).
# This is the command for one file
samtools view -f 4 alignment.bam > unmapped.bam

# We need to execute this for all of our files. It takes a while to exceute this for one file (around 3 minutes) so we can make an array script for each file to speed up the process. 
for file in ${rootdir}/STAR/Hg38_SecondPass/sample_output/STAR*.bam; do
	filename=`basename $file`
	path=`dirname $file`
	echo samtools view -f 4 $file '>' ${path}/unmapped_${filename}
done > unmapped_array.txt

# Convert files from BAM to Fasta format ---------------

# set dir name
star="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass"

# conversion of BAM files to FASTA format was followd using the following pipeline: http://www.metagenomics.wiki/tools/samtools/converting-bam-to-fastq
for file in ${star}/{sample_output,second_batch/sample_output,third_batch/sample_output}/unmapped_*.bam; do
	noExtension=${file%.bam}
	samtools bam2fq $file > ${noExtension}.fastq
done
                                                                 
# Split reads into two separate files- extract reads ending with '/1' or '/2'
for file in ${star}/{sample_output,second_batch/sample_output,third_batch/sample_output}/unmapped_*.fastq; do
	filename=`basename $file`
	path=`dirname $file`
	cat $file | grep '^@.*/1$' -A 3 --no-group-separator > ${path}/R1_${filename}
	cat $file | grep '^@.*/2$' -A 3 --no-group-separator > ${path}/R2_${filename}
done



