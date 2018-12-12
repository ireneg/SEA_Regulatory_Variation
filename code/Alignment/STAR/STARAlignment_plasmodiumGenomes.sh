# created by KSB on 09.11.2017
# use the genome aligner STAR to align indoRNA reads to plasmodium falciparum genome


#########################
# P. Falciparum Genome #
########################

# genome downloaded and GFF files from Plasmos DB, http://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-36_Pfalciparum3D7_Genome.fasta

# clean up file by removing everything after space
awk '{print $1}' Homo_sapiens.GRCh38.dna.primary_assembly.fa > Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa

# Generate genome indices
/vlsci/SG0008/shared/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum --genomeFastaFiles /vlsci/SG0008/kbobowik/genomes/PlasmoDB_36_Pfalciparum3D7_Genome.fasta --sjdbGTFfile /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_Pfalciparum3D7.gff --sjdbOverhang 100 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# First Pass- align all trimmed reads from the trimmomatic output to the human genome fasta file
for file in /vlsci/SG0008/shared/all_studyFiles/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  echo STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum --readFilesIn /vlsci/SG0008/shared/all_studyFiles/${f1} /vlsci/SG0008/shared/all_studyFiles/${f2} --readFilesCommand zcat --runThreadN 26 --sjdbOverhang 100 --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum/sample_output/STAR_pFalciparum_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 10000 --alignMatesGapMax 10000
done > STARArrayTable_FirstPass_pFalciparum.txt

# combine all of the SJ.out.tab files from from all files created in first pass into one file
cat /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum/sample_output/*.out.tab > /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum/sample_output/STAR_pFalciparum_allFilesSJ.out.tab

# Second pass- generate genomice indices
/vlsci/SG0008/shared/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum_SecondPass --genomeFastaFiles /vlsci/SG0008/kbobowik/genomes/PlasmoDB_36_Pfalciparum3D7_Genome.fasta --sjdbFileChrStartEnd /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum/sample_output/STAR_pFalciparum_allFilesSJ.out.tab --sjdbOverhang 100 --runThreadN 12 --sjdbGTFfile /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_Pfalciparum3D7.gff --sjdbGTFtagExonParentTranscript Parent --limitSjdbInsertNsj 2000000

# Map the reads
for file in /vlsci/SG0008/shared/all_studyFiles/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  echo STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum_SecondPass --readFilesIn /vlsci/SG0008/shared/all_studyFiles/${f1} /vlsci/SG0008/shared/all_studyFiles/${f2} --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum_SecondPass/sample_output/STAR_pFalciparum_second_${sample} --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 10000 --alignMatesGapMax 10000
done > STARArrayTable_SecondPass_pFalciparum.txt

# filter the bam file so that only uniquely mapped reads remain
for pathname in /vlsci/SG0008/kbobowik/Sumba/STAR/pFalciparum_SecondPass/sample_output/*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # this filters out singletons as well as the 255 flag
    samtools view -q 255 -f 3 $pathname > ${directory}/uniquelyMapped_${filename}
done


###################
# P. Vivax Genome #
###################

# genome downloaded and GFF files from Plasmos DB, http://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-36_PvivaxP01_Genome.fasta
# clean up file by removing everything after space
awk '{print $1}' PlasmoDB-36_PvivaxP01_Genome.fasta > PlasmoDB_36_PvivaxP01_Genome.fasta

# Generate genome indices
/vlsci/SG0008/shared/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax --genomeFastaFiles /vlsci/SG0008/kbobowik/genomes/PlasmoDB_36_PvivaxP01_Genome.fasta --sjdbGTFfile /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_PvivaxP01.gff --sjdbOverhang 100 --runThreadN 12 --sjdbGTFtagExonParentTranscript Parent

# First Pass- align all trimmed reads from the trimmomatic output to the human genome fasta file
for file in /vlsci/SG0008/shared/all_studyFiles/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  echo STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax --readFilesIn /vlsci/SG0008/shared/all_studyFiles/${f1} /vlsci/SG0008/shared/all_studyFiles/${f2} --readFilesCommand zcat --runThreadN 6 --sjdbOverhang 100 --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax/sample_output/STAR_pVivax_${sample} --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 10000 --alignMatesGapMax 10000
done > STARArrayTable_FirstPass_pVivax.txt

# combine all of the SJ.out.tab files from from all files created in first pass into one file
cat /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax/sample_output/*.out.tab > /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax/sample_output/STAR_pVivax_allFilesSJ.out.tab

# Second pass- generate genomice indices
/vlsci/SG0008/shared/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax_SecondPass --genomeFastaFiles /vlsci/SG0008/kbobowik/genomes/PlasmoDB_36_PvivaxP01_Genome.fasta --sjdbFileChrStartEnd /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax/sample_output/STAR_pVivax_allFilesSJ.out.tab --sjdbOverhang 100 --runThreadN 12 --sjdbGTFfile /vlsci/SG0008/kbobowik/genomes/PlasmoDB-36_PvivaxP01.gff --sjdbGTFtagExonParentTranscript Parent --limitSjdbInsertNsj 2000000

# Map the reads
for file in /vlsci/SG0008/shared/all_studyFiles/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  echo STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax_SecondPass --readFilesIn /vlsci/SG0008/shared/all_studyFiles/${f1} /vlsci/SG0008/shared/all_studyFiles/${f2} --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax_SecondPass/sample_output/STAR_pVivax_second_${sample} --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 11 --alignIntronMax 10000 --alignMatesGapMax 10000
done > STARArrayTable_SecondPass_pVivax.txt

# filter the bam file so that only uniquely mapped reads remain
for pathname in /vlsci/SG0008/kbobowik/Sumba/STAR/pVivax_SecondPass/sample_output/*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # this filters out singletons as well as the 255 flag
    samtools view -q 255 -f 3 $pathname > ${directory}/uniquelyMapped_${filename}
done
