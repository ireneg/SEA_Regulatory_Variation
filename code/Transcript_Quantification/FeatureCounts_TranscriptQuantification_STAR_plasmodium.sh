# created by KSB, 13.03.2018

# quantify the abundance of transcripts from the STAR alignment of plasmodium using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

################
# pFalciparum #
################

# fix last column in gff file by switching 'ID' with 'gene_id'
cp PlasmoDB-36_Pfalciparum3D7.gff PlasmoDB-36_Pfalciparum3D7_GeneID.gff
sed -i 's/ID=/gene_id /g' PlasmoDB-36_Pfalciparum3D7_GeneID.gff

# execute FeatureCounts
for file in /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum_SecondPass/sample_output/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  featureCounts -T 6 -s 2 -p -a /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_Pfalciparum3D7_GeneID.gff -t exon -g gene_id -o ~/Sumba/FeatureCounts/pFalciparum/sample_counts/${sample}.txt $file --verbose -B -C
  cat ~/Sumba/FeatureCounts/pFalciparum/sample_counts/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ~/Sumba/FeatureCounts/pFalciparum/sample_counts/Filter_GeneLengthCount_${sample}.txt

done

################
# pVivax #
################

# fix last column in gff file by switching 'ID' with 'gene_id'
cp PlasmoDB-36_PvivaxP01.gff PlasmoDB-36_PvivaxP01_GeneID.gff
sed -i 's/ID=/gene_id /g' PlasmoDB-36_PvivaxP01_GeneID.gff

# execute FeatureCounts
for file in /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax_SecondPass/sample_output/uniquelyMapped*.bam; do
  sample=`basename $file Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  featureCounts -T 6 -s 2 -p -a /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_PvivaxP01_GeneID.gff -t exon -g gene_id -o ~/Sumba/FeatureCounts/pVivax/sample_counts/${sample}.txt $file --verbose -B -C
    cat ~/Sumba/FeatureCounts/pVivax/sample_counts/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ~/Sumba/FeatureCounts/pVivax/sample_counts/Filter_GeneLengthCount_${sample}.txt

done
