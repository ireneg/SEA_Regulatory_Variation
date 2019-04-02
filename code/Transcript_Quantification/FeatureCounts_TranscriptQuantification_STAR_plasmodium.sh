# created by KSB, 13.03.2018

# quantify the abundance of transcripts from the STAR alignment of plasmodium using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

# set directory name
dir="/data/cephfs/punim0586/kbobowik"

# pFalciparum ------------

# fix last column in gff file by switching 'ID' with 'gene_id'
cp PlasmoDB-36_Pfalciparum3D7.gff PlasmoDB-36_Pfalciparum3D7_GeneID.gff
sed -i 's/ID=/gene_id /g' PlasmoDB-36_Pfalciparum3D7_GeneID.gff

# execute FeatureCounts
for file in ${dir}/STAR/pFalciparum_SecondPass/sample_output_allBatches/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  # note!!!! Subread version 1.6.1 seemed to break on me so I'm using the full path below (version 1.5.3) which seems to work
  /data/cephfs/punim0586/kbobowik/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 6 -s 2 -p -a ${dir}/genomes/PlasmoDB-36_Pfalciparum3D7_GeneID.gff -t exon -g gene_id -o ${dir}/Sumba/FeatureCounts/pFalciparum/sample_counts_allSamples/${sample}.txt $file --verbose -B -C
  cat ${dir}/Sumba/FeatureCounts/pFalciparum/sample_counts_allSamples/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${dir}/Sumba/FeatureCounts/pFalciparum/sample_counts_allSamples/Filter_GeneLengthCount_${sample}.txt

done

# pVivax ----------------------------------------------------

# fix last column in gff file by switching 'ID' with 'gene_id'
cp PlasmoDB-36_PvivaxP01.gff PlasmoDB-36_PvivaxP01_GeneID.gff
sed -i 's/ID=/gene_id /g' PlasmoDB-36_PvivaxP01_GeneID.gff

# execute FeatureCounts
for file in ${dir}/STAR/pVivax_SecondPass/sample_output_allBatches/uniquelyMapped*.bam; do
	sample=`basename $file Aligned.sortedByCoord.out.bam`
	ID=${sample##*_}
	# note!!!! Subread version 1.6.1 seemed to break on me so I'm using the full path below (version 1.5.3) which seems to work
	/data/cephfs/punim0586/kbobowik/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts -T 6 -s 2 -p -a ${dir}/genomes/PlasmoDB-36_PvivaxP01_GeneID.gff -t exon -g gene_id -o ${dir}/Sumba/FeatureCounts/pVivax/sample_counts_allSamples/${sample}.txt $file --verbose -B -C
    cat ${dir}/Sumba/FeatureCounts/pVivax/sample_counts_allSamples/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${dir}/Sumba/FeatureCounts/pVivax/sample_counts_allSamples/Filter_GeneLengthCount_${sample}.txt
done
