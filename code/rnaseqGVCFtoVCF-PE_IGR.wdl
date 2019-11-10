## rnaseq GATK best practices pipeline
## Bobbie Shaban
## 09/04/2018
## Melbourne Integrative genomics
##
## This GATK pipeline will take a number of bam files
## and convert them to GVCF files and the output
## it will employ scatter gather and hpc
## instructions to make are here: https://gatkforums.broadinstitute.org/wdl/discussion/6716/scatter-gather-parallelism
## and here: https://gatkforums.broadinstitute.org/wdl/discussion/7614/4-howto-use-scatter-gather-to-joint-call-genotypes

##Import tasks
#possible to reuse but will do this later

#import "./tasks/starAlignment.wdl" as starAlignmentTask
import "./tasks/combineGVCFs.wdl" as combineGVCFsTask
#import "./tasks/createBamIndex.wdl" as createBamIndexTask
#import "./tasks/createRefIndex.wdl" as createRefIndexTask
import "./tasks/GenotypeGVCFs.wdl" as GenotypeGVCFsTask      
# import "./tasks/indexCalibratedBam.wdl" as indexCalibratedBamTask   
#import "./tasks/picard.wdl" as picardTask      
# import "./tasks/sortBam.wdl" as sortBamTask          
# import "./tasks/baseRecalibrator1.wdl" as baseRecalibrator1Task
# import "./tasks/baseRecalibrator2.wdl" as baseRecalibrator2Task
#import "./tasks/convertSamToBam.wdl" as convertSamToBamTask
#import "./tasks/createPicardBamIndex.wdl" as createPicardBamIndexTask
# import "./tasks/generatePlots.wdl" as generatePlotsTask
import "./tasks/HaplotypeCallerERC.wdl" as HaplotypeCallerERCTask
#import "./tasks/picardMarkDuplicates.wdl" as picardMarkDuplicatesTask
# import "./tasks/printReads.wdl" as printReadsTask
#import "./tasks/splitNCigarReads.wdl" as splitNCigarReadsTask
import "./tasks/VariantFiltration.wdl" as VariantFiltrationTask
import "./tasks/copyOutput.wdl" as copyOutputTask

workflow rnaseqGVCFtoGVCF {

File inputSamplesFile
Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
File refFasta
File refDict
File refIndex
String picardLocation
String gatkLocation
String genomeDir
String outputDir
File dbSnpVcf
File dbSnpVcfIndex
File knownVcfs
File knownVcfsIndices

  # scatter (sample in inputSamples) {

    # call starAlignmentTask.starAlignment_task {
      # Int starAlignmentRunThreads
      # Int starAlignmentRunMinutes
      # Int starAlignmentRunMem
      # input:
          # genomeDir=genomeDir,
          # inputFastqRead1=sample[1],
          # inputFastqRead2=sample[2],
          # sampleName=sample[0]
    # }

    # call convertSamToBamTask.convertSamToBam_task {
      # Int convertSamToBamRunThreads
      # Int convertSamToBamRunMinutes
      # Int convertSamToBamRunMem
      # input:
        # alignmentSam=starAlignment_task.outputSam,
        # sampleName=sample[0]
    # }

    # call picardTask.picard_task {
      # String typeARRG
      # String sortOrder
      # String readGroupID
      # String readGroupLibrary
      # String readGroupPlatform
      # String readGroupPlatformBarcode
      # Int picardRunMinutes
      # Int picardThreads
    # input:
        # picardLocation=picardLocation,
        # picardInputBam=convertSamToBam_task.outputBam,
        # sampleName=sample[0]
    # }

    # call picardMarkDuplicatesTask.picardMarkDuplicates_task {
      # String typeMD
      # String createIndex
      # String validationStringency
      # String outputMetrics
      # String picardMarkDuplicatesRunMinutes
      # String picardMarkDuplicatesThreads

      # input:
        # picardLocation=picardLocation,
        # picardMDInputBam=picard_task.picardOutputBam,
        # sampleName=sample[0]
    # }

    # call createPicardBamIndexTask.createPicardBamIndex_task {
      # Int createPicardBamIndexRunMinutes
      # Int createPicardBamIndexThreads
      # Int createPicardBamIndexMem

      # input:
        # picardMDBamToBeIndexed=picardMarkDuplicates_task.picardDeduppedBam
    # }

    # call createRefIndexTask.createRefIndex_task{
      # Int createRefIndexRunMinutes
      # Int createRefIndexRunThreads
      # Int createRefIndexMem
      # input:
        # refFasta=refFasta
    # }

    # Split intron-spanning reads with N cigar strings
    # call splitNCigarReadsTask.splitNCigarReads_task{
      # String splitCigars
      # Int RF
      # Int RMQF
      # String RMQT
      # String U
      # Int splitNCigarReadsRunMinutes  
      # Int splitNCigarReadsThreads
      # Int splitNCigarReadsMem

      # input:
        # refFasta=refFasta,
        # refFastaIndex=createRefIndex_task.refFastaIndex,
        # refDictionary=refDict,
        # sampleName=sample[0],
        # gatkLocation=gatkLocation,
        # splitCigarsInputBam=picardMarkDuplicates_task.picardDeduppedBam,
        # splitCigarsInputBamIndex=createPicardBamIndex_task.mDBamIndex
    # }

    # Sorts and indexes bam file.
    # call sortBamTask.sortBam_task {
    #   input:
    #     bam2sort=sample[1],
    #     sampleName=sample[0]
    # }

    # # Index splitNcigar bamfile, pass this as output to next stage. Is there a reason this should be sorted? 
    # call createBamIndexTask.createBamIndex_task{
    #   Int createBamIndexRunMinutes
    #   Int createBamIndexThreads
    #   Int createBamIndexMem

    #   input:
    #     # bamToBeIndexed=splitNCigarReads_task.splitCigarsBamOutput
    #     # bamToBeIndexed=sortBam_task.sortedB
    # }

    # call baseRecalibrator1Task.baseRecalibrator1_task {
    #   input:
    #     gatkLocation=gatkLocation,
    #     sortedBam=sortBam_task.sortedBam,
    #     dbsnp=dbSnpVcf,
    #     goldStandard=knownVcfs,
    #     sampleName=sample[0],
    #     ref_fasta=refFasta,
    #     ref_fasta_index=refIndex,
    #     ref_dict=refDict,
    #     outputSortedBamIndex=sortBam_task.outputSortedBam
    # }

    # call baseRecalibrator2Task.baseRecalibrator2_task {
    #   input:
    #     gatkLocation=gatkLocation,
    #     sortedBam=sortBam_task.sortedBam,
    #     dbsnp=dbSnpVcf,
    #     goldStandard=knownVcfs,
    #     sampleName=sample[0],
    #     ref_fasta=refFasta,
    #     ref_fasta_index=refIndex,
    #     ref_dict=refDict,
    #     outputSortedBamIndex=sortBam_task.outputSortedBam,
    #     calibratedGrp=baseRecalibrator1_task.calibratedFile1
    # }

    # # QC plots
    # call generatePlotsTask.generatePlots_task {
    #   input:
    #     gatkLocation=gatkLocation,
    #     ref_fasta=refFasta,
    #     ref_fasta_index=refIndex,
    #     sampleName=sample[0],
    #     ref_dict=refDict,
    #     calibratedFile1=baseRecalibrator1_task.calibratedFile1,
    #     calibratedFile2=baseRecalibrator2_task.calibratedFile2
    # }

    # # Update quality scores and print new bam files.
    # call printReadsTask.printReads_task {
    #   input:
    #     gatkLocation=gatkLocation,
    #     ref_fasta=refFasta,
    #     ref_fasta_index=refIndex,
    #     ref_dict=refDict,
    #     sortedBam=sortBam_task.sortedBam,
    #     calibratedFile1=baseRecalibrator1_task.calibratedFile1,
    #     sortedBamIndex=sortBam_task.outputSortedBam,
    #     sampleName=sample[0]
    # }

  #   call indexCalibratedBamTask.indexCalibratedBam_task {
  #     input:
  #       sampleName=sample[0],
  #       refFasta=refFasta,
  #       calBam=sample[1]
  #   }

  #   call HaplotypeCallerERCTask.HaplotypeCallerERC_task {
  #     Int haplotypeCallerRunMinutes
  #     Int haplotypeCallerThreads
  #     Int haplotypeCallerMem

  #     input: GATK=gatkLocation, 
  #       RefFasta=refFasta, 
  #       RefIndex=refIndex,
  #       RefDict=refDict, 
  #       sampleName=sample[0],
  #       bamFile=indexCalibratedBam_task.calibratedBam, 
  #       bamIndex=indexCalibratedBam_task.calibratedBamIndex
  #   }
  # }

  # End of scattered jobs, end of first pipeline.

  call combineGVCFsTask.combineGVCFs_task {
    Int combineRunMinutes
    Int combineRunThreads
    Int combineRunMem

    input: 
      GATK=gatkLocation,
      RefFasta=refFasta,
      RefIndex=refIndex,
      RefDict=refDict,
      GVCFs=sample[1],
      sampleName="indoRNA_test6_combinedGVCFs"  
  }

  call GenotypeGVCFsTask.GenotypeGVCFs_task {
    Int genotypeRunMinutes
    Int genotypeThreads
    Int genotypeMem

    input: 
      GATK=gatkLocation, 
      RefFasta=refFasta, 
      RefIndex=refIndex, 
      RefDict=refDict, 
      sampleName="indoRNA_test6", 
      combinedVCF=combineGVCFs_task.combinedOutput
  }

  call VariantFiltrationTask.VariantFiltration_task {
    Int variantFilterRunMinutes
    Int variantFilterThreads
    Int variantFilterMem

    input:
      input_vcf = GenotypeGVCFs_task.rawVCF,
      input_vcf_index = GenotypeGVCFs_task.rawVCFidx,
      base_name = "indoRNA_test6",
      ref_fasta = refFasta,
      ref_fasta_index = refIndex,
      ref_dict = refDict,
      gatk_path=gatkLocation
  }

  call copyOutputTask.copyOutput_task {
    Int copyOutputRunThreads
    Int copyOutputRunMinutes
    Int copyOutputRunMem

    input:
      outputDir=outputDir,
      variantFiles=VariantFiltration_task.output_vcf,
      variantFilesIndex=VariantFiltration_task.output_vcf_index,
      variantFilesArray=GenotypeGVCFs_task.variantFiles,
      variantFilesIndexArray=GenotypeGVCFs_task.variantFilesIndex,
      haplotypeFiles=HaplotypeCallerERC_task.GVCF,
      generatedPlots=generatePlots_task.calibratedPlots
  }

}
#end workflow calls
