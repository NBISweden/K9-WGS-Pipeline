#!/bin/bash

while read line
do

######## Variables to change #########
name=$line #Where does the name start in the list? Usually like 29 or 33
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa
# The reference you want to use (indexed with script: /proj/b2013060/Script/Index_ref.sh)

echo "#!/bin/bash -l" > combineGVCF.sc;
echo "#SBATCH -A snic2017-1-403" >> combineGVCF.sc;
echo "#SBATCH -p core -n 1" >> combineGVCF.sc;
echo "#SBATCH -J $name-combineGVCFs1-160812" >> combineGVCF.sc;
echo "#SBATCH -t 96:00:00" >> combineGVCF.sc; #Changed the time
echo "#SBATCH -o /proj/b2013060/nobackup/output/$name-combineGVCFs1-160812.output" >> combineGVCF.sc;
echo "#SBATCH -e /proj/b2013060/nobackup/output/$name-combineGVCFs1-160812.error" >> combineGVCF.sc;
echo "#SBATCH --mail-user jessika.nordin@imbim.uu.se" >> combineGVCF.sc;
echo "#SBATCH --mail-type=ALL" >> combineGVCF.sc;

# load some modules
echo "module load bioinfo-tools" >> combineGVCF.sc;
echo "module load tabix" >> combineGVCF.sc;

# Variables, same as above #
echo "ref=$ref" >> combineGVCF.sc;
echo "name=$name" >> combineGVCF.sc;
echo "line=$line" >> combineGVCF.sc;


# Where should the program create all files?
echo "cd /proj/snic2017-1-403/nobackup/" >> combineGVCF.sc;

echo "java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T CombineGVCFs -R $ref \
-V ploidy1_samples.list -o $name-group1.vcf -L $name" >> combineGVCF.sc;

#echo "java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -R $ref \
# -T ValidateVariants --validationTypeToExclude ALL --variant $name-group1.vcf" >> combineGVCF.sc;
 
echo "bgzip $name-group1.vcf" >> combineGVCF.sc;

#let me know it is done
echo "echo "Job is done"" >> combineGVCF.sc;

echo "date" >> combineGVCF.sc;

sbatch combineGVCF.sc;

done <chromosome.list

