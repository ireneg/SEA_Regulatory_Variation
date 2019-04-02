# created 18.03.2019
# run FastQc and Trimmomatic on Yamagishi control samples

# initial FastQC
for i in `cat /vlsci/SG0008/kbobowik/SRA/SRR_Acc_List.txt`; do
  fastqc /scratch/SG0008/kbobowik/SRA_Files/${i}.fastq.gz --extract --outdir=/vlsci/SG0008/kbobowik/FastQC
  #transpose rows and add in column for sample number
cut -f1 ${i}_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="$i\t" '{print prefix $0}' >> ~/FastQC/allRuns_summary_file.txt
done
# add header to file
a=`cut -f2 ${i}_fastqc/summary.txt | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" allRuns_summary_file.txt

# trimmomatic- set minimum length to 32bp with a sliding window of 4bp averaging to a minimum phred score of 15. Leading and trailing bases are removed if below a quality of 3
for i in `cat /vlsci/SG0008/kbobowik/SRA/SRR_Acc_List.txt`; do  
  java -jar /vlsci/SG0008/shared/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -trimlog /vlsci/SG0008/kbobowik/QC/Trimmomatic/${i}.fastq.trim.log /scratch/SG0008/kbobowik/SRA_Files/${i}.fastq.gz /scratch/SG0008/kbobowik/SRA_Files/${i}_trimmed.fastq.gz SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:32 &> /vlsci/SG0008/kbobowik/QC/Trimmomatic/summary/summary_${i}.txt
done > /vlsci/SG0008/kbobowik/QC/Trimmomatic/Trimmomatic_ArrayTable.txt

# run FastQC a second time to check quality of the trimming
for i in `cat /vlsci/SG0008/kbobowik/SRA/SRR_Acc_List.txt`; do
  fastqc /scratch/SG0008/kbobowik/SRA_Files/${i}_trimmed.fastq.gz --extract --outdir=/vlsci/SG0008/kbobowik/QC/FastQC/SecondRun/
  #transpose rows and add in column for sample number
  cut -f1 /vlsci/SG0008/kbobowik/QC/FastQC/SecondRun/${i}_trimmed_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="$i\t" '{print prefix $0}' >> ~/QC/FastQC/allRuns_secondRun_summary_file.txt
done
# add header to file
a=`cut -f2 ${i}_fastqc/summary.txt | tr -s '\n' '\t'`
sed -i "1s/^/$a\n/" allRuns_summary_file.txt

# After viewing the files, we can see that two files, DRR006486 and DRR006487 have lost 75% and 82% of their reads, respectively. # remove DRR006486 and DRR006487, as they lost 75% and 82% of their reads, respectively. Three other reads, DRR006483-DRR006485, have poor quality scores past 35bp so will be cropped
for i in "DRR006483" "DRR006484" "DRR006485"; do
  echo $i
	java -jar /vlsci/SG0008/shared/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -trimlog /vlsci/SG0008/kbobowik/QC/Trimmomatic/${i}.fastq.trim.log /scratch/SG0008/kbobowik/SRA_Files/${i}.fastq.gz /scratch/SG0008/kbobowik/SRA_Files/${i}_trimmed.fastq.gz SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:32 CROP:35 &> /vlsci/SG0008/kbobowik/QC/Trimmomatic/summary/summary_${i}.txt
done

# Make sample summary for all trimmed reads
for i in `cat /vlsci/SG0008/kbobowik/SRA/SRR_Acc_List.txt`; do
    echo "$i `grep "Input" summary_${i}.txt`"
done > /vlsci/SG0008/kbobowik/QC/Trimmomatic/summary/TrimmomaticSampleSummary_AllSamples_minlength32_cropped.txt

# rerun FastQC on the three samples to see if quality has improved
for i in "DRR006483" "DRR006484" "DRR006485"; do
  fastqc /scratch/SG0008/kbobowik/SRA_Files/${i}_trimmed.fastq.gz --extract --outdir=/vlsci/SG0008/kbobowik/QC/FastQC/SecondRun/
  #transpose rows and add in column for sample number
  cut -f1 /vlsci/SG0008/kbobowik/QC/FastQC/SecondRun/${i}_trimmed_fastqc/summary.txt | tr -s '\n' '\t' | awk -v prefix="${i}_SecondTrim \t" '{print prefix $0}' >> ~/QC/FastQC/allRuns_secondRun_summary_file.txt
done


