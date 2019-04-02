from itertools import chain
from os.path import join
import os
import glob

SAMPLES = []
for path in glob.glob('{}/*.fastq.gz'.format(/data/cephfs/punim0586/kbobowik/Yamagishi/data/Sick)):
    dir, filename = os.path.split(path)
SAMPLES.append(filename.replace('.fastq.gz',''))

rule all:
    input:
        expand("/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.html", sample=SAMPLES),
        expand("/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/Sick/{sample}_trimmed.fastq.gz", sample=SAMPLES),
        expand("/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/Sick/SecondRun/{sample}_fastqc.html", sample=SAMPLES)

rule fastqc:
    input:
         "/data/cephfs/punim0586/kbobowik/Yamagishi/data/Sick/{sample}.fastq.gz"
    output:
        zip  = "/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.zip",
        html = "/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.html",
    shell:
        "fastqc {input} --outdir=/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/"

rule trimmomatic_se:
	input:
		"/data/cephfs/punim0586/kbobowik/Yamagishi/data/Sick/{sample}.fastq.gz"
	output:
		"/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/Sick/{sample}_trimmed.fastq.gz"
	log:
        "/data/cephfs/punim0586/kbobowik/Sumba/QC/Trimmomatic/Yamagishi/logs/Sick/{sample}.log"
    shell:
		"java -jar /data/cephfs/punim0586/kbobowik/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -trimlog {log} {input} {output} SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:32"

rule fastqc_trimmed:
    input:
         "/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/Sick/{sample}_trimmed.fastq.gz"
    output:
        zip  = "/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/Sick/SecondRun/{sample}_fastqc.zip",
        html = "/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/Sick/SecondRun/{sample}_fastqc.html",
    shell:
        "fastqc {input} --outdir=/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/Sick/SecondRun/"
