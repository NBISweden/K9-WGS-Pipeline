#! /bin/bash -l
#SBATCH -A b2014341
#SBATCH -p core -n 1
#SBATCH -J Genotyping
#SBATCH -t 200:00:00
#SBATCH -o /proj/b2013060/nobackup/output/genotyping.output
#SBATCH -e /proj/b2013060/nobackup/output/genotyping.error
#SBATCH --mail-user jessika.nordin@imbim.uu.se
#SBATCH --mail-type=ALL

# load some modules
module load bioinfo-tools
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa

cd /proj/b2013060/nobackup/finalbam/gvcf/

java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R $ref -V samples.list -o all_samples_genotyping.vcf.gz

#
# Example of samples.list:
# SpA-CA01_tag1-ploidy1.g.vcf.gz
# SpA-CA01_tag2-ploidy1.g.vcf.gz
# SpA-CA01_tag3-ploidy1.g.vcf.gz
# 


#let me know it is done
echo "Job is done"

date