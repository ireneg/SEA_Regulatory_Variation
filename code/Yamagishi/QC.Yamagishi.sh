# run trimmomatic on Yamagishi samples

dir=/data/cephfs/punim0586/kbobowik
for sample in `cat SRR_Acc_List_MalariaSamples.txt`; do
    echo java -jar ${dir}/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -trimlog ${dir}/Sumba/QC/Trimmomatic/Yamagishi/logs/Sick/${sample}.log ${dir}/Yamagishi/data/Sick/${sample}.fastq.gz ${dir}/processed_data/Yamagishi/Sick/${sample}_trimmed.fastq.gz SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:32
done > ${dir}/Sumba/scripts/Array_Scripts/Trimmomatic_ArrayTable_Yamagishi.txt
