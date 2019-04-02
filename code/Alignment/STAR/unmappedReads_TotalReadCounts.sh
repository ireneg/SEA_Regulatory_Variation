# created by KSB, 14.03.2019

# Get total number of unmapped read pairs after STAR alignemnt

# set directory name
dir="/data/cephfs/punim0586/kbobowik"

# load modules
module load SAMtools

# Quantify number of unmapped reads by counting total number of reads/2 (so we get read pairs rather than reads)
for file in ${dir}/STAR/Hg38_SecondPass/{sample_output,second_batch/sample_output,third_batch/sample_output}/unmapped*.bam; do
	# get batch information ftom the 8th field seperator in the sequence file
	batch=`echo $file | cut -d/ -f 8` sample=`basename $file Aligned.sortedByCoord.out.bam`
	# the first batch of sequences isn't labelled as 'first batch' (there is onlu 'sample_output' in teh 8th field separator) so we need to add this information in. 
	if [[ $batch == "sample_output" ]]; then
		batch="first_batch"
	fi
  	shortenedFile=`basename $file Aligned.sortedByCoord.out.bam`
	# ${variable##pattern} is like $variable, minus the longest matching pattern from front-end
	sampleID=${shortenedFile##*_}
	batchPlusID="${sampleID}_${batch}"
	samtools view $file | echo $batchPlusID `expr $(wc -l) / 2` >> ${dir}/STAR/unmappedReads_Hg38/unmappedReads_Counts.txt
	echo ${sampleID} done
done
