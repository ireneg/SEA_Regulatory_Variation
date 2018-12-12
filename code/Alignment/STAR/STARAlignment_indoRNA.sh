# created by KSB on 09.11.2017
# use the genome aligner STAR to align the Indonesian reads


# Get genome from Ensembl (release 90, GenBank Assembly ID GCA_000001405.25)
wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip, clean up file by removing everything after space. Rename to p10 and ensembl v90 to indicate build and source
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
awk '{print $1}' Homo_sapiens.GRCh38.dna.primary_assembly.fa > Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa

# Download and extract GTF file
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
gunzip /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf

# Generate genome indices
/vlsci/SG0008/shared/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38 --genomeFastaFiles /vlsci/SG0008/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbGTFfile /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 100 --runThreadN 12

# First Pass- align all trimmed reads from the trimmomatic output to the human genome fasta file
for file in /vlsci/SG0008/shared/indoRNA/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38 --readFilesIn /vlsci/SG0008/shared/indoRNA/${f1} /vlsci/SG0008/shared/indoRNA/${f2} --readFilesCommand zcat --runThreadN 6 --sjdbOverhang 100 --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38/sample_output/STAR_Hg38_${sample} --outSAMtype BAM SortedByCoordinate
done

# combine all of the SJ.out.tab files from from all files created in first pass into one file
cat /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38/sample_output/*.out.tab > /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38/sample_output/STAR_Hg38_allFilesSJ.out.tab

# Second pass- generate genomice indices
STAR --runMode genomeGenerate --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass --genomeFastaFiles /vlsci/SG0008/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbFileChrStartEnd /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38/sample_output/STAR_Hg38_allFilesSJ.out.tab --sjdbOverhang 100 --runThreadN 6 --sjdbGTFfile /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --limitSjdbInsertNsj 2000000

# Map the reads
for file in /vlsci/SG0008/shared/indoRNA/paired_trimmedOutput_*_1.fastq.gz; do
  f1=`basename $file`
  f2=`basename $file _1.fastq.gz`"_2.fastq.gz"
  sample=${f1%_*}
  STAR --genomeDir /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass --readFilesIn /vlsci/SG0008/shared/indoRNA/${f1} /vlsci/SG0008/shared/indoRNA/${f2} --outFileNamePrefix /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass/sample_output/STAR_Hg38_second_${sample} --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within
done

# filter the bam file so that only uniquely mapped reads remain
for pathname in /vlsci/SG0008/kbobowik/Sumba/STAR/Hg38_SecondPass/sample_output/*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # this filters out singletons as well as the 255 flag
    samtools view -b -q 255 -f 3 $pathname '>' ${directory}/uniquelyMapped_${filename}
done
