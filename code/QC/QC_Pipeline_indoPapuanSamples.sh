# created 19.09.2017 by KSB
# perform QC on indonesian and Papuan reads

# Globin-depleted (primary) samples
for path in /vlsci/SG0008/shared/raw_data/indoRNA/*.fastq.gz; do
  file=`basename $path`
  sample=`basename $path .fastq.gz`
  fastqc $path --extract --outdir=~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/FirstRun
  #transpose rows and add in column for sample number
  cut -f1 ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/FirstRun/${sample}_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="$sample\t" '{print prefix $0}' >> ~/QC/FastQC/IndoPapua/Globin_Depleted/FirstRun/allRuns_summary_file.txt
done
# add header to file
a=`cut -f2 ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/FirstRun/${sample}_fastqc/summary.txt | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/FirstRun/allRuns_summary_file.txt

# Globin Samples
for path in /vlsci/SG0008/shared/raw_data/indoRNA_no_globin/*.fastq.gz; do
  file=`basename $path`
  sample=`basename $path .fastq.gz`
  fastqc $path --extract --outdir=/vlsci/SG0008/kbobowik/Sumba/QC/FastQC/IndoPapua/Globin
  #transpose rows and add in column for sample number
  cut -f1 ~/Sumba/QC/FastQC/IndoPapua/Globin/${sample}_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="$sample\t" '{print prefix $0}' >> ~/Sumba/QC/FastQC/IndoPapua/Globin/allRuns_summary_file.txt
done
# add header to file
a=`cut -f2 ~/Sumba/QC/FastQC/IndoPapua/Globin/${sample}_fastqc/summary.txt | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" ~/Sumba/QC/FastQC/IndoPapua/Globin/allRuns_summary_file.txt

# Leading and trailing bases are removed if below a quality of 20
for file in /vlsci/SG0008/shared/raw_data/indoRNA/*_1.fastq.gz; do
  path=`dirname $file`
  f1=`basename $file`
  f2=${f1%%_*}"_2.fastq.gz"
  java -jar /vlsci/SG0008/shared/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 -trimlog ~/Sumba/QC/Trimmomatic/IndoPapua/${f1%%_*}.fastq.trim.log $path/${f1} $path/${f2} /vlsci/SG0008/shared/indoRNA/paired_trimmedOutput_${f1} /vlsci/SG0008/shared/indoRNA/unpaired_trimmedOutput_${f1} /vlsci/SG0008/shared/indoRNA/paired_trimmedOutput_${f2} /vlsci/SG0008/shared/indoRNA/unpaired_trimmedOutput_${f2} LEADING:20 TRAILING:20 MINLEN:90 &> ~/Sumba/QC/Trimmomatic/IndoPapua/summary/summary_${f1%%_*}.txt
done

# Run FastQC a second time on globin-depleted (primary) samples to check the quality of the trimming
for path in /vlsci/SG0008/shared/indoRNA/paired_trimmedOutput_*.fastq.gz; do
  file=`basename $path`
  sample=`basename $path .fastq.gz`
  fastqc $path --extract --outdir=/vlsci/SG0008/kbobowik/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/SecondRun
  #transpose rows and add in column for sample number
  cut -f1 ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/SecondRun/${sample}_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="$sample\t" '{print prefix $0}' >> ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/SecondRun/allRuns_summary_file.txt
  done
# add header to file
a=`cut -f2 ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/SecondRun/${sample}_fastqc/summary.txt | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" ~/Sumba/QC/FastQC/IndoPapua/Globin_Depleted/SecondRun/allRuns_summary_file.txt




