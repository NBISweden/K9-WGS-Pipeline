#! /bin/bash -l
#SBATCH -A b2014341
#SBATCH -p core -n 1
#SBATCH -J Genotyping_chrX_ploidy1
#SBATCH -t 200:00:00
#SBATCH -o /proj/b2013060/nobackup/output/genotyping_ploidy1.output
#SBATCH -e /proj/b2013060/nobackup/output/genotyping_ploidy1.error
#SBATCH --mail-user jessika.nordin@imbim.uu.se
#SBATCH --mail-type=ALL

# load some modules
module load bioinfo-tools

cd /proj/b2013060/nobackup/finalbam/gvcf/

java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.5.0/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /proj/snic2017-1-403/CanineGenomeResources/Autosome_fas/chrX.fa --sample_ploidy 1 -V ploidy1_samples.list \
-o ploidy_1.vcf.gz

# Possible options
# -L /proj/b2013060/nobackup/finalbam/chrX_pos.bed
# --includeNonVariantSites

#
# Example of ploidy1_samples.list:
# SpA-CA01_tag1-ploidy1.g.vcf.gz
# SpA-CA01_tag2-ploidy1.g.vcf.gz
# SpA-CA01_tag3-ploidy1.g.vcf.gz
# 


#let me know it is done
echo "Job is done"

date