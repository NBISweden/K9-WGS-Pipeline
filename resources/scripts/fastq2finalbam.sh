#!/bin/bash
# To run
#		./finalbam.sh finalbam.list			output in finalbam
#
# Developer: Jessika Nordin							Date: 2016-08-29
#
#	Change "Variables to change" if necessary 
#	and change to your email in the #SBATCH field below.
#

while read line
do

######## Variables to change #########
name=${line:16}	#Usually 16
# Ex: line=Sample_AID_H_JM_SpA-CA01_tag4=>name=SpA-CA01_tag4
folder=$line
# Which folder is used?
location=/proj/snic2017-1-403/Raw_data
# Path to the sample in the finalbam.list	
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa
# The reference you want to use (can be re-indexed with modifications of the script for bwa version wanted: /proj/snic2017-1-403/Scripts/Index_ref.sh)
known=/proj/snic2017-1-403/CanineGenomeResources/sorted_canine_snplist2012_canfam3.bed
# # Known sites for the recalibration with GATK (row 111)
# bait=/proj/b2013060/ReferenceTools/bait_hg19.bed
# # The bait region for HS metrics (row 139)
# target=/proj/b2013060/ReferenceTools/target_hg19.bed	
# # The target region for HS metrics (row 139)

echo "#!/bin/bash -l" > Finalbam.sc;
echo "#SBATCH -A b2013060" >> Finalbam.sc;
echo "#SBATCH -p core -n 3" >> Finalbam.sc;
echo "#SBATCH -J $name-160729" >> Finalbam.sc;
echo "#SBATCH -t 96:00:00" >> Finalbam.sc;
echo "#SBATCH -o /proj/b2013060/nobackup/output/$name.161005.output" >> Finalbam.sc;
echo "#SBATCH -e /proj/b2013060/nobackup/output/$name.161005.error" >> Finalbam.sc;
echo "#SBATCH --mail-user jessika.nordin@imbim.uu.se" >> Finalbam.sc;
echo "#SBATCH --mail-type=ALL" >> Finalbam.sc;

# load some modules
echo "module load bioinfo-tools" >> Finalbam.sc;
echo "module load samtools" >> Finalbam.sc;
echo "module load bwa/0.7.12" >> Finalbam.sc;
echo "module load tabix" >> Finalbam.sc;
#echo "module load R/3.1.0" >> Finalbam.sc;

# Variables, same as above #
echo "ref=$ref" >> Finalbam.sc;
echo "location=$location" >> Finalbam.sc;
echo "line=$line" >> Finalbam.sc;
echo "folder=$folder" >> Finalbam.sc;
echo "name=$name" >> Finalbam.sc;
echo "known=$known" >> Finalbam.sc;
echo "bait=$bait" >> Finalbam.sc;
echo "target=$target" >> Finalbam.sc;

# Where should the program create all files?
echo "cd /proj/b2013060/nobackup/finalbam/running" >> Finalbam.sc;

# Go from fastQ to BAM #
echo "bwa mem -t 6 -M -R '@RG\tID:$name\tSM:bar' $ref $location/$line/*_R1_* $location/$line/*_R2_* > \$SNIC_TMP/$name.sam" >> Finalbam.sc;
 
#Use PICARD, start with the SAM file and ask it to simultaneously add read groups,sort the file, and spit it out as BAM. 
echo "java -Xmx15g -jar /sw/apps/bioinfo/picard/1.92/milou/AddOrReplaceReadGroups.jar \
 INPUT=\$SNIC_TMP/$name.sam OUTPUT=$name.bam SORT_ORDER=coordinate RGID=$name-id RGLB=$name-lib \
 RGPL=ILLUMINA RGPU=$name-01 RGSM=$name VALIDATION_STRINGENCY=LENIENT" >> Finalbam.sc;

# Index the Picard derived BAM
echo "java -Xmx23g -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar INPUT=$name.bam \
  VALIDATION_STRINGENCY=LENIENT" >> Finalbam.sc;

echo "rm $name.sam" >> Finalbam.sc;
# 
# #Process reads with GATK
#Step 1 - realign locally around potential indels. First, identify possible sites to realign.
echo "java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -I $name.bam -R $ref \
 -T RealignerTargetCreator -o $name.intervals" >> Finalbam.sc;

#Feed the intervals file back into GATK with a different argument to actually do the realignments
echo "java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -I $name.bam -R $ref \
 -T IndelRealigner -o $name.realigned.bam -targetIntervals $name.intervals" >> Finalbam.sc;

echo "rm   $name.ba* $name.intervals" >> Finalbam.sc;
 
#Use Picard to mark duplicate reads
echo "java -Xmx23g -jar /sw/apps/bioinfo/picard/1.92/milou/MarkDuplicates.jar INPUT=$name.realigned.bam \
 OUTPUT=$name.realigned.marked.bam METRICS_FILE=$name.realigned.marked.metrics VALIDATION_STRINGENCY=LENIENT" >> Finalbam.sc;
 
echo "rm $name.realigned.ba*" >> Finalbam.sc;

#Reindexing the picard file
echo "java -Xmx23g -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar INPUT=$name.realigned.marked.bam \
 VALIDATION_STRINGENCY=LENIENT" >> Finalbam.sc;

#Perform quality recalibration with GATK. #
#Do this last, as we want all the data to be as clean as possible. Two Step Process. 

#Step 1 - We need a list of known sites otherwise, GATK will think all the real SNPs in our data are errors.
#Failure to remove real SNPs from the recalibration will result in globally lower quality scores
echo "java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T BaseRecalibrator -I $name.realigned.marked.bam \
 -R $ref -knownSites $known -o $name.recal_data.table -rf BadCigar" >> Finalbam.sc; 
#If needed use "-rf BadCigar"
#Step2 - 		#Changes have been made
echo "java -Xmx6g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T BaseRecalibrator -I $name.realigned.marked.bam \
 -R $ref -knownSites $known -BQSR $name.recal_data.table -o $name.post_recal_data.table " >> Finalbam.sc;

#Not gotten step 3 to work. Creates a before and after plot when recalibrating. 
#Step 3 - Generate before and after QS plots
# echo "java -Xmx6g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T AnalyzeCovariates \
#     -R $ref -before $name.recal_data.table -after $name.post_recal_data.table \
#     -plots ../quality-test/$name.recalibration_plots.pdf" >> Finalbam.sc;
    
#Step4 - 		#Changes have been made
echo "java -Xmx6g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T PrintReads -I $name.realigned.marked.bam \
 -R $ref -BQSR $name.recal_data.table -o $name.realigned.marked.recalibrated.bam" >> Finalbam.sc;


echo "rm $name.realigned.marked.ba* $name.realigned.marked.metrics" >> Finalbam.sc;

#Quality test #
#Doing flagstat
echo "samtools flagstat $name.realigned.marked.recalibrated.bam > ../quality-test/$name.flagstat" >> Finalbam.sc;

#  #HSmetrics
# echo "java -Xmx6g -jar /sw/apps/bioinfo/picard/1.92/milou/CalculateHsMetrics.jar R=$ref BAIT_INTERVALS=$bait \
#  TARGET_INTERVALS=$target INPUT=$name.realigned.marked.recalibrated.bam OUTPUT=../quality-test/$name.hybrid_selection_metrics \
#  VALIDATION_STRINGENCY=LENIENT" >> Finalbam.sc;

# Delete files #
# Delete files to free up space
echo "rm core* $name.recal_data.table $name.post_recal_data.table" >> Finalbam.sc;

#Move BAM to backupspace #
echo "cd /proj/b2013060/nobackup/finalbam/" >> Finalbam.sc;
echo "mv /proj/b2013060/nobackup/finalbam/running/$name.realigned.marked.recalibrated.bam $name.final.bam" >> Finalbam.sc;
echo "mv /proj/b2013060/nobackup/finalbam/running/$name.realigned.marked.recalibrated.bai $name.final.bai" >> Finalbam.sc;

echo "cd /proj/b2013060/nobackup/haplotypeCaller/" >> Finalbam.sc;

### HaplotypeCaller for all samples ###
# Males should also have their chromosome X genotyped with haplotypecaller_ploidy1.sh
echo "java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref \
  -I /proj/b2013060/nobackup/finalbam/$name.final.bam --emitRefConfidence GVCF --variant_index_type LINEAR \
  --variant_index_parameter 128000 -o $name.g.vcf -rf BadCigar" >> Finalbam.sc;
 
echo "bgzip $name.g.vcf" >> Finalbam.sc;
echo "tabix $name.g.vcf.gz" >> Finalbam.sc;

#let me know it is done
echo "echo "Job is done"" >> Finalbam.sc;

echo "date" >> Finalbam.sc;

sbatch Finalbam.sc;

done <finalbam.list
