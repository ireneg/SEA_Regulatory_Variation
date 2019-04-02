SAMPLES_FASTQC = ['DRR017562', 'DRR017563', 'DRR017564', 'DRR017565','DRR017566','DRR017567','DRR017568','DRR017569','DRR017570','DRR017571','DRR017572','DRR017573','DRR017574','DRR017575','DRR017576','DRR017577','DRR017578','DRR017579','DRR017580','DRR017581','DRR017582','DRR017583', 'DRR017584', 'DRR017585', 'DRR017586', 'DRR017587', 'DRR017588', 'DRR017589', 'DRR017590','DRR017591','DRR017592','DRR017593','DRR017594','DRR017595','DRR017596','DRR017597','DRR017598','DRR017599','DRR017600','DRR017601','DRR017602','DRR017603','DRR017604','DRR017605','DRR017606','DRR017607','DRR017608','DRR017609','DRR017610','DRR017611','DRR017612','DRR017613','DRR017614','DRR017615','DRR017616','DRR017617','DRR017618','DRR017619','DRR017620','DRR017621','DRR017622']

rule all:
    input:
        expand("/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.html", sample=SAMPLES_FASTQC),
        expand("/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/{sample}_trimmed.fastq.gz", sample=SAMPLES_FASTQC),
        expand("/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/controls/SecondRun/{sample}_fastqc.html", sample=SAMPLES_FASTQC)

rule fastqc:
    input:
         "/data/cephfs/punim0586/kbobowik/Yamagishi/data/{sample}.fastq.gz"
    output:
        zip  = "/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.zip",
	 	html = "/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/{sample}_fastqc.html",
    shell:
        "fastqc {input} --outdir=/data/cephfs/punim0586/kbobowik/Yamagishi/FastQC/reports/"

rule trimmomatic_se:
	input:
		"/data/cephfs/punim0586/kbobowik/Yamagishi/data/{sample}.fastq.gz"
	output:
		"/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/{sample}_trimmed.fastq.gz"
	log:
        "/data/cephfs/punim0586/kbobowik/Sumba/QC/Trimmomatic/Yamagishi/logs/controls/{sample}.log"
    shell:
		"java -jar /data/cephfs/punim0586/kbobowik/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -trimlog {log} {input} {output} SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:32"

rule fastqc_trimmed:
    input:
         "/data/cephfs/punim0586/kbobowik/processed_data/Yamagishi/controls/{sample}_trimmed.fastq.gz"
    output:
        zip  = "/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/controls/SecondRun/{sample}_fastqc.zip",
        html = "/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/controls/SecondRun/{sample}_fastqc.html",
    shell:
        "fastqc {input} --outdir=/data/cephfs/punim0586/kbobowik/Sumba/QC/FastQC/Yamagishi/controls/SecondRun/"
