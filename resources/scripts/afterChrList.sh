#! /bin/bash -l
#
# This script is used after using a the "chromosome.list". 
# Takes all chromosomes.vcf and make them to one 
# To run
#		sbatch --dependency=afterok:8497735:8497735 afterChrList.sh			output in haplotypeCaller
#	

#SBATCH -A snic2017-1-403
#SBATCH -p core -n 3
#SBATCH -J afterChrList
#SBATCH -t 196:00:00
#SBATCH -o /proj/b2013060/nobackup/output/afterChrList.output
#SBATCH -e /proj/b2013060/nobackup/output/afterChrList.error
#SBATCH --mail-user jessika.nordin@imbim.uu.se
#SBATCH --mail-type=ALL

# load some modules
module load bioinfo-tools
module load vcftools
module load tabix 

## Variables ##
# CombineGvcfs groups 
name=group1
# Genotyping 
# name=all_samples_genotyping
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa

#Which directory are you working in?
cd /proj/snic2017-1-403/nobackup/


# Combine the vcfs
java -Xmx23g -cp /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
  -R $ref  -V chr1-$name.vcf.gz -V chr2-$name.vcf.gz -V chr3-$name.vcf.gz -V chr4-$name.vcf.gz -V chr5-$name.vcf.gz \
  -V chr6-$name.vcf.gz -V chr7-$name.vcf.gz -V chr8-$name.vcf.gz -V chr9-$name.vcf.gz -V chr10-$name.vcf.gz \
  -V chr11-$name.vcf.gz -V chr12-$name.vcf.gz -V chr13-$name.vcf.gz -V chr14-$name.vcf.gz -V chr15-$name.vcf.gz \
  -V chr16-$name.vcf.gz -V chr17-$name.vcf.gz -V chr18-$name.vcf.gz -V chr19-$name.vcf.gz -V chr20-$name.vcf.gz \
  -V chr21-$name.vcf.gz -V chr22-$name.vcf.gz -V chrX-$name.vcf.gz -V chrY-$name.vcf.gz -V chrM-$name.vcf.gz \
  -out $name.vcf.gz -assumeSorted


#let me know it is done
echo "Job is done"

date