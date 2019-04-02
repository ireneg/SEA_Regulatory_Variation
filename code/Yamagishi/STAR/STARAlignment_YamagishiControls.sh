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

# Generate genome indices
# Note: jobs will get killed of you don't ask for enough memory. I found asking for at least 12 CPUs and 200,000 MB was sufficient
/data/cephfs/punim0586/kbobowik/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Controls/Hg38/Genome --genomeFastaFiles /data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbGTFfile /data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 35 --runThreadN 12

# First Pass
for file in /data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/*trimmed.fastq.gz; do
  sample=`basename $file _trimmed.fastq.gz`
  echo STAR --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Controls/Hg38/Genome --sjdbGTFfile /data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --readFilesIn $file --readFilesCommand zcat --runThreadN 12 --sjdbOverhang 35 --outFileNamePrefix /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/STAR_Hg38_${sample} --outSAMtype BAM SortedByCoordinate
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/STARArrayTable_FirstPass_YamagishiControls.txt
# run as sbatch script 

# combine all of the SJ.out.tab files from from all files created in first pass into one file 
cat /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/*.out.tab > /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/STAR_PFPX_allFilesSJ.out.tab

# generate genomice indices for the second time
/data/cephfs/punim0586/kbobowik/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Yamagishi/Hg38_SecondPass/Controls/Genome --genomeFastaFiles /data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbFileChrStartEnd /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/STAR_PFPX_allFilesSJ.out.tab --sjdbGTFfile /data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 35 --runThreadN 12

# Map the reads for the second pass
# note unmapped reads get output automatically as a bam file when using the --outSAMunmapped Within flag
for file in /data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/*trimmed.fastq.gz; do
  sample=`basename $file _trimmed.fastq.gz`
  echo STAR --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Yamagishi/Hg38_SecondPass/Controls/Genome --readFilesIn $file --outFileNamePrefix /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls/STAR_Hg38_second_${sample} --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/STARArrayTable_SecondPass_YamagishiControls.txt

# filter the bam file so that only uniquely mapped reads remain
for pathname in /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls/STAR_*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # "-q 255" = unique reads
    echo samtools view -q 255 $pathname '>' ${directory}/uniquelyMapped_${filename} 
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/uniquelyMappedReads_ArrayTable_YamagishiControls.txt
