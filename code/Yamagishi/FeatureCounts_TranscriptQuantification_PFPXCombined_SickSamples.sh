# created by KSB, 13.03.2018

# quantify the abundance of transcripts from the STAR alignment of plasmodium using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

# load modules
module load Subread

# set directory name
dir=/data/cephfs/punim0586/kbobowik
genomeDir=/data/cephfs/punim0586/kbobowik/genomes
stardir=/data/cephfs/punim0586/kbobowik/STAR/Yamagishi_CombinedPFPX_SecondPass/Sick

# execute FeatureCounts 
# here are the following flags I used:
# -T: Number of the threads.  The value should be between 1 and 32.  1 by default.
# -s: Indicate if strand-specific read counting should be performed. Acceptable  values:  0  (unstranded),  1  (stranded)  and  2  (re-versely stranded).  0 by default.  For paired-end reads, strand of the first read is taken as the strand of the whole fragment. FLAG field is used to tell if a read is first or second read in a pair.
# -p: If specified, fragments (or templates) will be counted instead of reads.  This option is only applicable for paired-end reads.
# -a: Give the name of an annotation file
# - t: Specify the feature type.  Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting.  ‘exon’ by default.
# -o: Give the name of the output file
# -g: Specify the attribute type used to group features (eg.  exons) into meta-features (eg. genes)

# execute FeatureCounts
for file in ${stardir}/STAR*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  # note!!!! Subread version 1.6.1 seemed to break on me so I'm using the full path below (version 1.5.3) which seems to work
  /data/cephfs/punim0586/kbobowik/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 12 -a ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff -t exon -g gene_id -o ${dir}/Sumba/FeatureCounts/PFPX_Combined_Yamagishi/Sick/${sample}.txt $file --verbose -B -C
  cat ${dir}/Sumba/FeatureCounts/PFPX_Combined_Yamagishi/Sick/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${dir}/Sumba/FeatureCounts/PFPX_Combined_Yamagishi/Sick/Filter_GeneLengthCount_${sample}.txt
done
