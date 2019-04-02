# Created 14.02.19
# Find out what unmapped sequences from STAR bam files are mapping to
 
# first, laod the SRA-Toolkit module
module load SRA-Toolkit
module load web_proxy

# set outdir
outdir=/data/cephfs/punim0586/kbobowik/Yamagishi/data

# First, download accession list (SRR_Acc_List.txt under 'Download') from ncbi: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRP001953

# download each sample by reading line by line from accession list file
for i in `cat SRR_Acc_List.txt`; do
	echo fastq-dump --gzip --outdir $outdir $i
done > Yamagishi_Controls_Array.txt
