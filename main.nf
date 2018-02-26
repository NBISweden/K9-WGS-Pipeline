
if (params.help) {
    exit 0, usageMessage()
}

checkInputParams()

reference = file(params.reference)
refdir = file(reference.getParent())
referenceBwaIndex = file("${reference}.{amb,ann,bwt,pac,sa}")
referenceFaIndex = file("${reference}.fai")
referenceDict = file("${refdir}/${reference.getBaseName()}.dict")
known     = file(params.known)
outdir    = params.out

if ( refdir.name.contains('test-data-tiny') ) {
    chromosomes = ['chr38'] 
    }
else if ( refdir.name.contains('test-data-small') ) {
    chromosomes = (36..38).collect {"chr${it}_1000000_1030000"} + ['chrX_1000000_1030000']
    }
else {
    chromosomes = (1..38).collect {"chr${it}"} + ['chrX', 'chrY', 'chrM']
    }

fastqFiles = Channel.fromFilePairs(params.fastqDir + '/*R{1,2}.fq.gz')
fastqFiles.into { fastq_qc; fastq_bwa }


infoMessage()


process fastqc {
    input:
        set val(key), file(fastqs) from fastq_qc

    output:
        file "*_fastqc.{zip,html}" into fastQCreport

    publishDir "out"


    script:
    """
    fastqc -q $fastqs
    """
}


process bwa {
    input:
        set val(key), file(fastqs) from fastq_bwa
        set file(reference), file(bwaindex) from Channel.value([reference, referenceBwaIndex])
    output:
        set val(key), file("${key}.bam"), file("${key}.bam.bai") into mapped_reads

    publishDir "out"

    tag "$key"


    script:
    readGroup = "@RG\\tID:Sample_79162\\tSM:bar\\tPL:ILLUMINA"
    """
    bwa mem -t $task.cpus -M -R \'$readGroup\' $reference $fastqs |
        samtools sort --threads $task.cpus -m 4G > ${key}.bam
    samtools index ${key}.bam
    """
}


process gatk_realign {
    input:
        set val(key), file(bamfile), file(bamix) from mapped_reads
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        file('name.intervals')
        set val(key), file("${key}.realigned.bam"), file("${key}.realigned.bai") into realigned_reads

    publishDir 'out'


    script:
    """
    ls
    java -jar /usr/GenomeAnalysisTK.jar \
        -I $bamfile \
        -R $reference \
        -T RealignerTargetCreator \
        -o name.intervals

    java -jar /usr/GenomeAnalysisTK.jar \
        -I $bamfile \
        -R $reference \
        -T IndelRealigner \
        -targetIntervals name.intervals \
        -o ${key}.realigned.bam

    """
}


process mark_duplicates {
    input:
        set val(key), file(bamfile), file(bamix) from realigned_reads
    output:
        set val(key), file("${key}.realigned.marked.bam"), file("${key}.realigned.marked.bai") into marked_reads

    publishDir 'out'


    script:
    """
    picard MarkDuplicates.jar \
        INPUT=$bamfile \
        OUTPUT=${key}.realigned.marked.bam \
        METRICS_FILE=${key}.realigned.marked.metrics \
        VALIDATION_STRINGENCY=LENIENT

    picard BuildBamIndex.jar \
        INPUT=${key}.realigned.marked.bam \
        VALIDATION_STRINGENCY=LENIENT
    """
}


process quality_recalibration {
    input:
        set val(key), file(bamfile), file(bamix) from marked_reads
        file known
    output:
        file("${key}.recal_data.table")
        file("${key}.post_recal_data.table")
        set val(key), file("${key}.recalibrated.bam"), file("${key}.recalibrated.bai") into recalibrated_bam
        file("${key}.recalibration_plots.pdf")

    publishDir 'out'


    script:
    """
    # Step 1 - We need a list of known sites otherwise, GATK will think all the
    #          real SNPs in our data are errors. Failure to remove real SNPs
    #          from the recalibration will result in globally lower quality
    #          scores
    java -jar /usr/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -I $bamfile \
        -R $reference \
        -knownSites $known \
        -o ${key}.recal_data.table \
        -rf BadCigar

    # Step 2
    java -jar /usr/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -I $bamfile \
        -R $reference \
        -knownSites $known \
        -BQSR ${key}.recal_data.table \
        -o ${key}.post_recal_data.table

    # Step 3 - Create before and after plots
    java -jar /usr/GenomeAnalysisTK.jar \
        -T AnalyzeCovariates \
        -R $reference \
        -before ${key}.recal_data.table \
        -after ${key}.post_recal_data.table \
        -plots ${key}.recalibration_plots.pdf

    # Step 4 - Base calibration
    java -jar /usr/GenomeAnalysisTK.jar \
        -T PrintReads \
        -I $bamfile \
        -R $reference \
        -BQSR ${key}.recal_data.table \
        -o ${key}.recalibrated.bam
    """
}


recalibrated_bam.tap { recalibrated_bam_flagstats; recalibrated_bam_hsmetrics; recalibrated_bam_haplotype }


process flagstats {
    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_flagstats
    output:
        file("${key}.flagstat")

    publishDir 'out'


    script:
    """
    samtools flagstat $bamfile > ${key}.flagstat
    """
}


process haplotypeCaller {
    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_haplotype
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        set val(key), file("${key}.g.vcf") into haplotype_caller


    script:
    """
    java -jar /usr/GenomeAnalysisTK.jar \
        -T HaplotypeCaller\
        -R $reference \
        -I $bamfile \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -o ${key}.g.vcf \
        -rf BadCigar
    """
}


process haplotypeCallerCompress {
    input:
        set val(key), file(vcffile) from haplotype_caller
    output:
        set val(key), file("${key}.g.vcf.gz"), file("${key}.g.vcf.gz.tbi") into compress_haplocalled

    publishDir 'out'


    script:
    """
    bgzip $vcffile
    tabix ${vcffile}.gz
    """
}


process hsmetrics {
    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_hsmetrics
        file reference

    output:
        file("${key}.hybridd_selection_metrics")

    publishDir 'out'

    when: false

    script:
    """
    picard CalculateHsMetrics.jar \
        R=$reference \
        BAIT_INTERVALS=$bait \
        TARGET_INTERVALS=$target \
        INPUT=$bamfile \
        OUTPUT=${key}.hybrid_selection_metrics \
        VALIDATION_STRINGENCY=LENIENT
    """
}

/* COMBINE GVCF
This needs to collect 150/200 samples from samples.list
Unsure how to get those
Script also runs one per chromosome it seems, could use each
*/

compress_haplocalled
  .reduce([keys:[], vcfs:[], ixvcfs:[]]) { a, b -> a.keys.add(b[0]); a.vcfs.add(b[1]); a.ixvcfs.add(b[2]); return a }
  .map { it -> [it.keys, it.vcfs, it.ixvcfs] }
  .set { collect_haplovcfs }


process gVCFCombine {

    input:
    set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    set val(keys), file(vcfs), file(ix_vcfs) from collect_haplovcfs
    each chrom from chromosomes
     
    output:
    file "${chrom}.vcf" into combined

    script:
    """
    #count=1; for k in ${keys.join(' ')}; do ln -s `pwd`/vcf\$count \$k.vcf.gz; ln -s `pwd`/ix_vcf\$count \$k.vcf.gz.tbi; ((count++)); done
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -V ${vcfs.join(' -V ')} \
        -R $reference \
        -o ${chrom}.vcf -L $chrom
    """
}


process bgZipCombinedGVCF {

    input:
    file combined_gvcf from combined

    output:
    file "${combined_gvcf}.gz" into compressed_comb_gvcf

    """
    bgzip $combined_gvcf
    """
}

compressed_comb_gvcf
  .collect()
  .set { combgvcfs }


process afterChrList {
    input:
    file reference
    file chromosomes from combgvcfs

    output:
    file 'chromosomes.vcf.gz' into combined_chromosomes

    script:
    """
    java -Xmx23g -cp /usr/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
        -R $reference \
        -V ${chromosomes.join(' -V ')} \
        -out chromosomes.vcf.gz -assumeSorted
    """
}


/*
process genotype {
    input:
        set file('file.vcf.gz'), file('file.vcf.gz.tbi') from input
        file reference
    output:
        file 'genotyped.vcf.gz'


    script:
    """
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R $reference \
        -V samples.list \
        -o all_samples_genotyping.vcf.gz
    """
}
*/


////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////


def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false
    if ( ! params.fastqDir ) {
        log.warn("You need to provide a fastqDir (--fastqDir)")
        fatal_error = true
    }
    if ( ! params.reference ) {
        log.warn("You need to provide a genome reference (--reference)")
        fatal_error = true
    }
    if ( ! params.known ) {
        log.warn("You need to provide a file of known sites (--known)")
        fatal_error = true
    }

    if (fatal_error) {
        log.warn "See --help for more information"
        exit 1
    }
}

def usageMessage() {
    log.info """\
    Usage:
        nextflow run main.nf --fastqDir <directory>
    Options:
        --fastqDir <Dir>
           Directory containing fastq samples (.fq.gz)
        --reference <file>
           Genome reference file (has to be indexed)
        --known <file>
           File with known sites for quality calibration.
        --out <dir>
           Directory for output files
    """
}

def infoMessage() {
    log.info """\
*** K9 WGS Pipeline ***
Configuration environemnt:
    Out directory: $params.out
    """
}
