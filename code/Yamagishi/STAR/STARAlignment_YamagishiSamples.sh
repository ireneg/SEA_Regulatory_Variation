# created by KSB, 09.11.2017
# use the genome aligner STAR to align the 122 Yamagishi samples

# download and tidy up genome ----------------------

# Get genome from Ensembl (release 90, GenBank Assembly ID GCA_000001405.25)
wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename to p10 and ensembl v90 to indicate build and source
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa.gz

# Download GTF file
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz

# Extract genome and GTF files
gunzip /data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa
gunzip /data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf


# Run STAR -------------------------------------------

# set dir names
rootdir=/data/cephfs/punim0586/kbobowik
stardir=${rootdir}/STAR/Hg38/Yamagishi
stardirSecond=${rootdir}/STAR/Hg38_SecondPass/Yamagishi
genomedir=/data/cephfs/punim0586/shared/genomes/hg38
trimmedDataDir=/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/Sick
scriptdir=/data/cephfs/punim0586/kbobowik/Sumba/scripts

# Generate genome indices
# Note: jobs will get killed of you don't ask for enough memory. I found asking for at least 12 CPUs and 200,000 MB was sufficient
${rootdir}/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir ${rootdir}/STAR_Yamagishi/Hg38/Genome --genomeFastaFiles ${genomedir}/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbGTFfile ${genomedir}/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 35 --runThreadN 12

# First Pass
for file in ${trimmedDataDir}/*trimmed.fastq.gz; do
  sample=`basename $file _trimmed.fastq.gz`
  echo STAR --genomeDir ${rootdir}/STAR_Yamagishi/Hg38/Genome --sjdbGTFfile ${genomedir}/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --readFilesIn $file --readFilesCommand zcat --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix ${stardir}/Sick/STAR_Hg38_${sample} --outSAMtype BAM SortedByCoordinate
done > ${scriptdir}/Array_Scripts/STARArrayTable_FirstPass_YamagishiSick.txt
# run as sbatch script STAR_FirstPass_Yamagishi.sh

# combine all of the SJ.out.tab files from from all files created in first pass into one file 
cat ${stardir}/Sick/*.out.tab > ${stardir}/Sick/STAR_PFPX_allFilesSJ.out.tab

# generate genomice indices for the second time
${rootdir}/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir ${rootdir}/STAR_Yamagishi/Hg38_SecondPass/Sick/Genome --genomeFastaFiles ${genomedir}/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbFileChrStartEnd ${stardir}/Sick/STAR_PFPX_allFilesSJ.out.tab --sjdbGTFfile ${genomedir}/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 35 --runThreadN 12

# Map the reads for the second pass
for file in ${trimmedDataDir}/*trimmed.fastq.gz; do
  sample=`basename $file _trimmed.fastq.gz`
  echo STAR --genomeDir ${rootdir}/STAR_Yamagishi/Hg38_SecondPass/Sick/Genome --readFilesIn $file --outFileNamePrefix ${stardirSecond}/Sick/STAR_Hg38_second_${sample} --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within
done > ${scriptdir}/Array_Scripts/STARArrayTable_SecondPass_YamagishiSick.txt
# run as sbatch script STAR_FirstPass_Yamagishi.sh

# filter the bam file so that only uniquely mapped reads remain
for pathname in ${stardirSecond}/Sick/STAR_*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # "-q 255" = unique reads
    echo samtools view -q 255 $pathname '>' ${directory}/uniquelyMapped_${filename} 
done > ${scriptdir}/Array_Scripts/uniquelyMappedReads_ArrayTable_Yamagishi.txt

