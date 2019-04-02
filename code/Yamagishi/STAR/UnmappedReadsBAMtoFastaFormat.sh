# created by KSB, 14.03.2019
# Get total number of unmapped read pairs of Yamagishi samples (indlucing controls) after STAR alignemnt

# Load modules -------------------

# load samtools 1.9
module load SAMtools

# specify directory
rootdir=/data/cephfs/punim0586/kbobowik
yamagishi=/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Sick
scriptdir=/data/cephfs/punim0586/kbobowik/Sumba/scripts

# Filter reads for only unmapped reads ---------------

# the samtools -f flag outputs alignmemnts only matching this flag. The '4' flag (or 00100) corresponds to all unmapped reads (as per the Samtools manual: http://www.htslib.org/doc/samtools-1.2.html and the Picard 'explain samtools flags' website: http://broadinstitute.github.io/picard/explain-flags.html).
# note: this happens automatically when using the --outSAMunmapped Within flag in STAR

# We need to execute this for all of our files. It takes a while to exceute this for one file (around 3 minutes) so we can make an array script for each file to speed up the process. 
for file in ${rootdir}/STAR/Hg38_SecondPass/Yamagishi/{Controls,Sick}/STAR*.bam; do
	filename=`basename $file`
	path=`dirname $file`
	echo samtools view -f 4 $file '>' ${path}/unmapped_${filename}
done > ${scriptdir}/Array_Scripts/unmapped_YamagshiSamples_array.txt

# Convert files from BAM to Fasta format ---------------

# set dir name
star="/data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass"

# conversion of BAM files to FASTA format was followd using the following pipeline: http://www.metagenomics.wiki/tools/samtools/converting-bam-to-fastq
for file in ${star}/Yamagishi/{Controls,Sick}/unmapped_*.bam; do
	noExtension=${file%.bam}
	samtools bam2fq $file > ${noExtension}.fastq
done


