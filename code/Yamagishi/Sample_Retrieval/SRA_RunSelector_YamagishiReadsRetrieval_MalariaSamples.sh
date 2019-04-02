# CReated 28.03.18 by KSB - all previoulsy-downlaoded Yamagishi files were on scratch, which is now gone. Fantastic.
# retrieve samples without any extra flags using the SRA toolkit
# this is taken from SRA study DRP000987. The accession list file was downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRP000987&go=go

module load SRA-Toolkit
module load web_proxy

ref_dir=/data/cephfs/punim0586/kbobowik/Yamagishi/data
out_dir=/data/cephfs/punim0586/kbobowik/Yamagishi/data/Sick
main_dir=/data/cephfs/punim0586/kbobowik

for i in `cat ${ref_dir}/SRR_Acc_List_MalariaSamples.txt`; do
	echo fastq-dump --gzip --outdir $out_dir $i
done > ${main_dir}/Sumba/scripts/Array_Scripts/Yamagishi_MalariaSamples_Array.txt

# execute sbatch script 

# see if all files completed successfully
vdb-dump --info DRR006383

