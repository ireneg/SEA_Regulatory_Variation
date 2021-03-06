# created by KSB, 14.03.2019

# Get total number of unmapped read pairs after STAR alignemnt

# set directory name
dir="/data/cephfs/punim0586/kbobowik"

# load modules
module load SAMtools

# Quantify number of unmapped reads by counting total number of reads/2 (so we get read pairs rather than reads)
for file in ${dir}/STAR/Hg38_SecondPass/Yamagishi/{Controls,Sick}/unmapped*.bam; do
	# get batch information ftom the 8th field seperator in the sequence file
	healthStatus=`echo $file | cut -d/ -f 9` 
	# the first batch of sequences isn't labelled as 'first batch' (there is onlu 'sample_output' in teh 8th field separator) so we need to add this information in. 
	if [[ $batch == "sample_output" ]]; then
		batch="first_batch"
	fi
  	shortenedFile=`basename $file Aligned.sortedByCoord.out.bam`
	# ${variable##pattern} is like $variable, minus the longest matching pattern from front-end
	sampleID=${shortenedFile##*_}
	healthStatusPlusID="${sampleID}_${healthStatus}"
	samtools view $file | echo $healthStatusPlusID `expr $(wc -l) / 2` >> ${dir}/STAR/unmappedReads_Hg38/unmappedReads_Counts_Yamagishi.txt
	echo ${sampleID} done
done
