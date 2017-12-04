#! /bin/bash -l
#SBATCH -A b2013060
#SBATCH -p core -n 3
#SBATCH -J super_filter
#SBATCH -t 96:00:00
#SBATCH -o /proj/b2013060/nobackup/output/s_filter.output
#SBATCH -e /proj/b2013060/nobackup/output/s_filter.error
#SBATCH --mail-user jessika.nordin@imbim.uu.se
#SBATCH --mail-type=ALL

# load some modules
module load bioinfo-tools
module load vcftools
module load tabix
module load samtools
module load R/3.0.1
module load snpEff/4.0
module load plinkseq

## Variables ##
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa
 
# #Do we have a canine "good" set to use in VQSR???
# hapmap=/proj/b2013060/ReferenceTools/hapmap.vcf 
# omni=/proj/b2013060/ReferenceTools/omni.vcf
# snphc=/proj/b2013060/ReferenceTools/1000Gsnp.vcf  
# dbsnp=/proj/b2013060/ReferenceTools/dbSNP_147.vcf
# mills=/proj/b2013060/ReferenceTools/mills.vcf 


cd /proj/b2013060/nobackup/haplotypeCaller

### VQSR ###

# SNP #
java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T VariantRecalibrator \
 -R $ref --input ploidy_one.vcf.gz \
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
 -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $snphc \
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
 -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 \
 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal \
 -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R --maxGaussians 4 

#Use recalibration model on snp call set
java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref \
 --input ploidy_one.vcf.gz -mode SNP --ts_filter_level 99.0 \
 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o SNP_raw_indels.vcf 
 
bgzip SNP_raw_indels.vcf 

# Get only positions that PASSED VQSR and are polymorphic
# Final_SpA-set = should contain the names on samples taken forward (e.g. not the family)	--keep Final_SpA-set
vcftools --gzvcf SNP_raw_indels.vcf.gz --keep-filtered PASS --out VQSR_SNPs --remove-filtered-geno-all  \
 --remove-filtered-all --recode --recode-INFO-all

 