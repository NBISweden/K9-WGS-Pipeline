
if (params.help) {
    exit 0, usageMessage()
}

checkInputParams()

reference         = file(params.reference)
refdir            = file(reference.getParent())
referenceBwaIndex = file("${reference}.{amb,ann,bwt,pac,sa}")
referenceFaIndex  = file("${reference}.fai")
referenceDict     = file("${refdir}/${reference.getBaseName()}.dict")
known             = file(params.known)
outdir            = params.out

chromosomes = (1..38).collect {"chr${it}"} + ['chrX', 'chrY', 'chrM']

fastqFiles = Channel.fromFilePairs(params.fastqDir + '/*R{1,2}.fq.gz')
fastqFiles.into { fastq_qc; fastq_bwa }

bamFiles = Channel.fromPath(params.bamDir + '*.bam')
bamFiles
  .map { [it.baseName, it, infer_bam_index_from_bam(it)] }
  .set { readyBamFiles }


infoMessage()


process fastqc {
    tag "$key"

    input:
        set val(key), file(fastqs) from fastq_qc
    output:
        file "*fastqc.{html,zip}"

    publishDir "${params.out}/report", mode: 'copy'


    when: params.fastqDir

    script:
    """
    fastqc -q $fastqs
    """
}


process bwa {
    tag "$key"

    input:
        set val(key), file(fastqs) from fastq_bwa
        set file(reference), file(bwaindex) from Channel.value([reference, referenceBwaIndex])
    output:
        set val(key), file("${key}.bam"), file("${key}.bam.bai") into mapped_reads


    when: params.fastqDir && ! params.bamDir

    script:
    readGroup = "@RG\\tID:$key\\tSM:$key\\tPL:ILLUMINA"
    """
    bwa mem -t $task.cpus -M -R \'$readGroup\' $reference $fastqs |
        samtools sort --threads $task.cpus -m 4G > ${key}.bam
    samtools index ${key}.bam
    """
}


process gatk_realign {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from mapped_reads
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        set val(key), file("${key}.realigned.bam"), file("${key}.realigned.bai") into realigned_reads


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
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from realigned_reads
    output:
        set val(key), file("${key}.marked.bam"), file("${key}.marked.bai") into marked_reads
        file("*.metrics")

    publishDir "${params.out}/report", mode: 'copy', saveAs: { it =~ /metrics$/ ? it : null }


    script:
    """
    picard MarkDuplicates.jar \
        INPUT=$bamfile \
        OUTPUT=${key}.marked.bam \
        METRICS_FILE=${key}.marked.metrics \
        VALIDATION_STRINGENCY=LENIENT

    picard BuildBamIndex.jar \
        INPUT=${key}.marked.bam \
        VALIDATION_STRINGENCY=LENIENT
    """
}


process quality_recalibration {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from marked_reads
        file known
    output:
        file("${key}.*recal_data.table")
        set val(key), file("${key}.recalibrated.bam"), file("${key}.recalibrated.bai") into recalibrated_bam

    publishDir "${params.out}", mode: 'copy', saveAs: {
            type = it =~ /(ba(m|i))$/;
            type ? "bam/${key}.${type[0][1]}" : "report/$it"
        }


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

    # Step 3 - Base calibration
    java -jar /usr/GenomeAnalysisTK.jar \
        -T PrintReads \
        -I $bamfile \
        -R $reference \
        -BQSR ${key}.recal_data.table \
        -o ${key}.recalibrated.bam
    """
}


recalibrated_bam.concat( readyBamFiles ).tap { recalibrated_bam_flagstats; recalibrated_bam_hsmetrics; recalibrated_bam_haplotype }


process flagstats {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_flagstats
    output:
        file("${key}.*stat*")

    publishDir "${params.out}/report", mode: 'copy'


    script:
    """
    samtools flagstat $bamfile > ${key}.flagstat
    samtools stats    $bamfile > ${key}.stats
    """
}


process haplotypeCaller {
    tag "$key"

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
    tag "$key"

    input:
        set val(key), file(vcffile) from haplotype_caller
    output:
        set val(key), file("${key}.g.vcf.gz"), file("${key}.g.vcf.gz.tbi") into compress_haplocalled

    publishDir "${params.out}/haplotypeCaller", mode: 'copy'


    script:
    """
    bgzip $vcffile
    tabix ${vcffile}.gz
    """
}


process hsmetrics {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_hsmetrics
        file reference
    output:
        file("${key}.hybridd_selection_metrics")

    publishDir "${params.out}/report"


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

compress_haplocalled.toList().transpose().toList().set { collect_haplovcfs }

gVCFCombine_ch = Channel.create()
genotyping = Channel.create()
if ( params.combineByChromosome ) {
    collect_haplovcfs.set { gVCFCombine_ch }
    genotyping.close()
}
else {
    collect_haplovcfs.set { genotyping }
    gVCFCombine_ch.close()
}


process grab_chromosome_names_from_reference {
    input:
        file reference from Channel.value(reference)
    output:
        file 'chromosomes.names' into chromosome_names

    executor 'local'

    script:
    """
    grep '^>' $reference | sed 's/^>//' > chromosomes.names
    """
}

process filter_chromosomes {
    tag "$chrom"

    input:
        file chrom_names from chromosome_names
        each chrom from chromosomes
    output:
        val chrom into chromosomes_existing

    errorStrategy 'ignore'
    executor 'local'

    script:
    """
    grep -q '^$chrom\$' $chrom_names
    """
}


process gVCFCombine {
    tag "$chrom"

    input:
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
        set val(key), file(vcfs), file(ix_vcfs) from gVCFCombine_ch
        each chrom from chromosomes_existing
    output:
        set val(chrom), file("${chrom}.vcf") into combined


    when params.combineByChromosome

    script:
    """
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -V ${vcfs.join(' -V ')} \
        -R $reference \
        -o ${chrom}.vcf -L $chrom
    """
}


process bgZipCombinedGVCF {
    tag "$chrom"

    input:
        set val(chrom), file(combined_gvcf) from combined
    output:
        set val(chrom), file("${combined_gvcf}.gz"), file("*.gz.tbi") into compressed_comb_gvcf


    """
    bgzip $combined_gvcf
    tabix ${combined_gvcf}.gz
    """
}


genotyping.mix( compressed_comb_gvcf )
          .map { it[0] instanceof List ? ['all', it[1], it[2]] : it }
          .set { genotyping }


process genotype {
    tag "$key"

    input:
        set val(key), file(vcfs), file(ix_vcfs) from genotyping
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        set val(key), file("${key}_genotyping.vcf.gz"), file("${key}_genotyping.vcf.gz.tbi") into genotyped

    publishDir "${params.out}/genotype", mode: 'copy'


    script:
    """
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R $reference \
        -V ${vcfs.join(' -V ')} \
        -o ${key}_genotyping.vcf.gz
    """
}


hardfilters.into { hardfilters_snp; hardfilters_indel }


process hardfilters_snp {
    input:
        set file(vcf), file(index) from hardfilters_snp
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        set file('*SNP*vcf'), file('filtered_snps*vcf')

    publishDir "${params.out}/genotype", mode: 'copy'


    script:
    """
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T SelectVariants -R $reference \
        -V $vcf -selectType SNP -o raw_SNP_1.vcf

    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V raw_SNP_1.vcf \
        --filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
        --filterExpression "FS > 60.0" --filterName "FS_greater_than_60_filter" \
        --filterExpression "SOR > 3.0" --filterName "SOR_greater_than_3_filter" \
        --filterExpression "MQ < 40.0" --filterName "MQ_less_than_40_filter" \
        --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum_less_than_-12.5_filter" \
        --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
        -o filtered_snps_1.vcf \

    vcftools --vcf filtered_snps_1.vcf --keep-filtered PASS --out pass_SNP_1 --remove-filtered-geno-all  \
     --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992
     """
}


process hardfilters_indel {
    input:
        set file(vcf), file(index) from hardfilters_indel
        set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
    output:
        set file('*INDEL*vcf'), file('filtered_indels*vcf')

    publishDir "${params.out}/genotype", mode: 'copy'


    script:
    """
    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T SelectVariants -R $reference \
        -V $vcf -selectType INDEL -o raw_INDEL_1.vcf

    java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V raw_INDEL_1.vcf \
        --filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
        --filterExpression "FS > 200.0" --filterName "FS_greater_than_60_filter" \
        --filterExpression "SOR > 10.0" --filterName "SOR_greater_than_3_filter" \
        --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
        -o filtered_indels_1.vcf

    vcftools --vcf filtered_indels_1.vcf  --keep-filtered PASS --out pass_INDEL_1 --remove-filtered-geno-all  \
     --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992
    """
}


////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////


def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false
    if ( ! params.fastqDir && ! params.bamDir ) {
        log.warn("You need to provide a fastqDir (--fastqDir) or a bamDir (--bamDir)")
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
        --bamDir <Dir>
           Instead of --fastqDir, directory containing bam and indexed bam sample files (.bam, .bam.bai)
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

def infer_bam_index_from_bam(f) {
    // If the ".bam.bai" file does not exist, try ".bai" without ".bam"
    return infer_filepath(f, /$/, '.bai')
        ?: infer_filepath(f, /.bam$/, '.bai')
        ?: filepath_from(f, /$/, '.bai') // Default filename if none exist
}

def filepath_from(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
    return path
}

def infer_filepath(from, match, replace) {
    path = file( from.toString().replaceAll(match, replace) )
    if (path.exists()) {
        return path
    }
    return false
}
