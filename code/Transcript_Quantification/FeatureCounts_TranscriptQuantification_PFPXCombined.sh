# created by KSB, 13.03.2018

# quantify the abundance of transcripts from the STAR alignment of plasmodium using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

# set directory name
dir=/data/cephfs/punim0586/kbobowik
genomeDir=/data/cephfs/punim0586/kbobowik/genomes
stardir=/data/cephfs/punim0586/kbobowik/STAR/IndoRNA_CombinedPFPX_SecondPass/IndoOutput

# execute FeatureCounts
for file in ${stardir}/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  # note!!!! Subread version 1.6.1 seemed to break on me so I'm using the full path below (version 1.5.3) which seems to work
  /data/cephfs/punim0586/kbobowik/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 6 -s 2 -p -a ${genomeDir}/combined_PFalc3D7_PvivaxP01_GFF.gff -t exon -g gene_id -o ${dir}/Sumba/FeatureCounts/PFPX_Combined_IndoSamples/${sample}.txt $file --verbose -B -C
  cat ${dir}/Sumba/FeatureCounts/PFPX_Combined_IndoSamples/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${dir}/Sumba/FeatureCounts/PFPX_Combined_IndoSamples/Filter_GeneLengthCount_${sample}.txt
done
