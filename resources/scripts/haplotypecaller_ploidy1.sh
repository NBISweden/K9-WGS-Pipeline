#!/bin/bash

while read line
do

######## Variables to change #########
name=$line	#Usually 16
# Ex: line=Sample_AID_H_JM_SpA-CA01_tag4=>name=SpA-CA01_tag4	
ref=/proj/snic2017-1-403/CanineGenomeResources/canFam3.fa 
# The reference you want to use (indexed with script: /proj/b2013060/Script/Index_ref.sh)

echo "#!/bin/bash -l" > Finalbam.sc;
echo "#SBATCH -A b2013060" >> Finalbam.sc;
echo "#SBATCH -p core -n 2" >> Finalbam.sc;
echo "#SBATCH -J $name-sex-170814" >> Finalbam.sc;
echo "#SBATCH -t 24:00:00" >> Finalbam.sc;
#echo "#SBATCH -o /proj/b2013060/nobackup/finalbam/$name-sex.170814.output" >> Finalbam.sc;
#echo "#SBATCH -e /proj/b2013060/nobackup/finalbam/$name-sex.170814.error" >> Finalbam.sc;
echo "#SBATCH --mail-user jessika.nordin@imbim.uu.se" >> Finalbam.sc;
#echo "#SBATCH --mail-type=ALL" >> Finalbam.sc;

# load some modules
echo "module load bioinfo-tools" >> Finalbam.sc;
echo "module load tabix" >> Finalbam.sc;

# Variables, same as above #
echo "ref=$ref" >> Finalbam.sc;
echo "name=$name" >> Finalbam.sc;

# Where should the program create all files?
echo "cd /proj/snic2017-1-403/CanineGenomeResources/" >> Finalbam.sc;

echo "java -Xmx15g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref --sample_ploidy 1 \
	-I $name.final.bam --emitRefConfidence GVCF --variant_index_type LINEAR \
	--variant_index_parameter 128000 -o gvcf/$name-ploidy1.g.vcf -L chrX -rf BadCigar" >> Finalbam.sc;  


echo "bgzip $name.g.vcf" >> Finalbam.sc;
echo "bgzip $name-ploidy1.g.vcf" >> Finalbam.sc;
echo "tabix $name-ploidy1.g.vcf.gz" >> Finalbam.sc;


#let me know it is done
echo "echo "Job is done"" >> Finalbam.sc;

echo "date" >> Finalbam.sc;

sbatch Finalbam.sc;

done <finalbam.list
