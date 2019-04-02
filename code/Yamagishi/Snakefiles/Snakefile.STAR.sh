
"""
A simple Snakemake's workflow to process RNA-Seq using STAR
taken from https://www.biostars.org/p/219748/
and: https://github.com/saketkc/rna-seq-snakemake/blob/master/Snakefile_sra_pe
"""
import os

configfile: "cluster-configs/config.yaml"

samples = ['DRR017562', 'DRR017563', 'DRR017564', 'DRR017565','DRR017566','DRR017567','DRR017568','DRR017569','DRR017570','DRR017571','DRR017572','DRR017573','DRR017574','DRR017575','DRR017576','DRR017577','DRR017578','DRR017579','DRR017580','DRR017581','DRR017582','DRR017583', 'DRR017584', 'DRR017585', 'DRR017586', 'DRR017587', 'DRR017588', 'DRR017589', 'DRR017590','DRR017591','DRR017592','DRR017593','DRR017594','DRR017595','DRR017596','DRR017597','DRR017598','DRR017599','DRR017600','DRR017601','DRR017602','DRR017603','DRR017604','DRR017605','DRR017606','DRR017607','DRR017608','DRR017609','DRR017610','DRR017611','DRR017612','DRR017613','DRR017614','DRR017615','DRR017616','DRR017617','DRR017618','DRR017619','DRR017620','DRR017621','DRR017622']

rule all:
    input: '/scratch/punim0586/kat/STAR_Controls/Hg38'

rule index_star:
    input: fa = "/data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa",
           gtf = "/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf"
    output: '/scratch/punim0586/kat/STAR_Controls/Hg38'
    threads: 12
    shell:
        r'''mkdir -p {output} && STAR --runMode genomeGenerate        \
                   --readFilesCommand zcat         \
                   --runThreadN {threads}          \
                   --genomeFastaFiles {input.fa}   \
                   --sjdbGTFfile {input.gtf}       \
                   --sjdbOverhang 35               \
                   --genomeDir {output}'''

# run using snakemake --resources mem_mb=100000
------
# example 1
rule align_star:
    """
    Align sequencing reads using STAR.
    """
    input: fq1 = '{prefix}_1.fq.gz',
           fq2 = '{prefix}_2.fq.gz',
           idx = rules.index_star.output.Genome
    output: bam = '{prefix}.star.bam',
            sj = '{prefix}.star.sj'
    threads: 8
    run:
        res_dir = os.path.dirname(output.bam)
        idx_dir = os.path.dirname(input.idx)
        shell("""
              STAR --runThreadN {threads}                 \
                   --runMode alignReads                   \
                   --readFilesCommand pigz -d -c          \
                   --outSAMtype BAM Unsorted              \
                   --genomeDir {idx_dir}/                 \
                   --outFileNamePrefix {res_dir}/         \
                   --readFilesIn {input.fq1} {input.fq2}
              mv {res_dir}/Aligned.out.bam {output.bam}
              mv {res_dir}/SJ.out.tab {output.sj}   
              """)

# example #2
rule map_star:
    input:
        R1='preprocessed/{sample}/{sample}_1_val_1.fq',
        R2='preprocessed/{sample}/{sample}_2_val_2.fq',
        index=STAR_INDEX
    output: 'mapped/bams/star/{sample}.bam'
    params:
        prefix = 'mapped/bams/star/{sample}',
        unmapped = 'unmapped/fastq/star/{sample}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STAR --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --outFilterMatchNmin 50\
             --outFilterMismatchNmax 100\
             --outReadsUnmapped {params.unmapped} && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
'''

-----------

rule sort_bam:
    """
    Sort bam by coordinates using sambamba.
    """
    input: bam = '{prefix}.bam'
    output: bam = '{prefix}.sorted.bam'
    params: mem = '35G'
    threads: 8
    run:
        tmp_dir = os.path.dirname(output.bam)
        shell('sambamba sort --tmpdir {tmp_dir} -t {threads} -m {params.mem} -o {output.bam} {input.bam}')

--------------

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
cat /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/*.out.tab > /data/cephfs/punim0586/kbobowik/STAR/Hg38/Yamagishi/Controls/STAR_Hg38_allFilesSJ.out.tab

# generate genomice indices
/data/cephfs/punim0586/kbobowik/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --readFilesCommand zcat --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Controls/Hg38_SecondPass/Genome --genomeFastaFiles /data/cephfs/punim0586/shared/genomes/hg38/Homo_sapiens.GRCh38.p10.ensemblv90.dna.primary_assembly.fa --sjdbGTFfile /data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 35 --runThreadN 12

# Map the reads for the second pass
for file in /data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/*trimmed.fastq.gz; do
  sample=`basename $file _trimmed.fastq.gz`
  echo STAR --genomeDir /data/cephfs/punim0586/kbobowik/STAR_Controls/Hg38_SecondPass/Genome --readFilesIn $file --outFileNamePrefix /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls/STAR_Hg38_second_${sample} --runThreadN 12 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/STARArrayTable_SecondPass_YamagishiControls.txt

# filter the bam file so that only uniquely mapped reads remain
for pathname in /data/cephfs/punim0586/kbobowik/STAR/Hg38_SecondPass/Yamagishi/Controls/STAR_*.bam; do
    filename=`basename $pathname`
    directory=`dirname $pathname`
    # "-q 255" = unique reads; -b = bam file output
    echo samtools view -b -q 255 $pathname '>' ${directory}/uniquelyMapped_${filename} 
done > /data/cephfs/punim0586/kbobowik/Sumba/scripts/Array_Scripts/uniquelyMappedReads_ArrayTable_YamagishiControls.txt

