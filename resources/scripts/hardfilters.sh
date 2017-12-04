#! /bin/bash -l
#SBATCH -A b2013060
#SBATCH -p core -n 1
#SBATCH -J Filter_chrX_ploidy1
#SBATCH -t 100:00:00
#SBATCH -o /proj/b2013060/nobackup/output/filter_ploidy1.output
#SBATCH -e /proj/b2013060/nobackup/output/filter_ploidy1.error
#SBATCH --mail-user jessika.nordin@imbim.uu.se
#SBATCH --mail-type=ALL

# load some modules
module load bioinfo-tools
module load vcftools


# Choose reference according to if it is chrX or all others that are going to be filtered.
 ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa 
# ref=/proj/snic2017-1-403/CanineGenomeResources/Autosome_fas/chrX.fa

cd /proj/b2013060/nobackup/finalbam/gvcf/

java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants -R ref \
	-V minDP_GQ_1.recode.vcf -selectType SNP -o raw_SNP_1.vcf

java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_SNP_1.vcf \
	--filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
    --filterExpression "FS > 60.0" --filterName "FS_greater_than_60_filter" \
    --filterExpression "SOR > 3.0" --filterName "SOR_greater_than_3_filter" \
    --filterExpression "MQ < 40.0" --filterName "MQ_less_than_40_filter" \
    --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum_less_than_-12.5_filter" \
    --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
    -o filtered_snps_1.vcf \

vcftools --vcf filtered_snps_1.vcf --keep-filtered PASS --out pass_SNP_1 --remove-filtered-geno-all  \
 --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992
    
java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T SelectVariants -R ref \
	-V minDP_GQ_1.recode.vcf -selectType INDEL -o raw_INDEL_1.vcf

java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_INDEL_1.vcf \
    --filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
	--filterExpression "FS > 200.0" --filterName "FS_greater_than_60_filter" \
    --filterExpression "SOR > 10.0" --filterName "SOR_greater_than_3_filter" \
    --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
    -o filtered_indels_1.vcf 

vcftools --vcf filtered_indels_1.vcf  --keep-filtered PASS --out pass_INDEL_1 --remove-filtered-geno-all  \
 --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992

#let me know it is done
echo "Job is done"

date

