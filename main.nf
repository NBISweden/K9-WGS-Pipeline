
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

chromosomes = params.chromosomes.split(',').collect { it.trim() }

fastqFiles = Channel.fromFilePairs(params.fastqDir + '/*R{1,2}*.f*q.gz')
fastqFiles.into { fastq_qc; fastq_bwa }

bdir = params.bamDir + '/*.bam'
bamFiles = Channel.fromPath(bdir)
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

    publishDir "${params.outdir}/reports", mode: 'copy'


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


    // errorStrategy 'ignore'

    publishDir "${params.outdir}/bam-temp"

    when: params.fastqDir && ! params.bamDir

    script:
    readGroup = "@RG\\tID:$key\\tSM:$key\\tPL:ILLUMINA"
    """
    bwa mem -t $task.cpus -M -R \'$readGroup\' $reference $fastqs > temp.sam

    samtools sort --threads ${task.cpus - 2} -m ${params.singleCPUMem.toMega() - 1000}M < temp.sam > ${key}.bam

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
    java -Xmx${params.singleCPUMem.toMega() - 500}m \
        -jar /usr/GenomeAnalysisTK.jar \
        -I $bamfile \
        -R $reference \
        -T RealignerTargetCreator \
        -o name.intervals

    java -Xmx${params.singleCPUMem.toMega() - 500}m \
        -jar /usr/GenomeAnalysisTK.jar \
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

    publishDir "${params.outdir}/reports", mode: 'copy', saveAs: { it =~ /metrics$/ ? it : null }


    script:
    """
    picard -Xmx${params.singleCPUMem.toMega() - 500}m \
        MarkDuplicates \
        INPUT=${bamfile} \
        OUTPUT=${key}.marked.bam \
        METRICS_FILE=${key}.marked.metrics \
        TMP_DIR=. \
        CREATE_INDEX=TRUE \
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

    publishDir "${params.outdir}", mode: 'copy', saveAs: {
            type = it =~ /(ba(m|i))$/;
            type ? "bam/${key}.${type[0][1]}" : "reports/$it"
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


recalibrated_bam.concat( readyBamFiles ).tap { recalibrated_bam_flagstats; recalibrated_bam_wgsmetrics; recalibrated_bam_haplotype }


process flagstats {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_flagstats
    output:
        file("${key}.*stat*")

    publishDir "${params.outdir}/reports", mode: 'copy'


    script:
    """
    samtools flagstat $bamfile > ${key}.flagstat
    samtools stats    $bamfile > ${key}.stats
    """
}


process wgsmetrics {
    tag "$key"

    input:
        set val(key), file(bamfile), file(bamix) from recalibrated_bam_wgsmetrics
        file reference
    output:
        file("${key}.wgs_metrics")

    publishDir "${params.outdir}/reports", mode: 'copy'


    script:
    """
    picard CollectWgsMetrics \
        R=$reference \
        INPUT=$bamfile \
        OUTPUT=${key}.wgs_metrics
    """
}


if ( !params.onlyMap ) {
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

        publishDir "${params.outdir}/haplotypeCaller", mode: 'copy'


        script:
        """
        bgzip $vcffile
        tabix ${vcffile}.gz
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


    process gVCFCombine {
        tag "$chrom"

        input:
            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
            set val(key), file(vcfs), file(ix_vcfs) from gVCFCombine_ch
            each chrom from chromosomes
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


        script:
        """
        java -Xmx7g -jar /usr/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R $reference \
            -V ${vcfs.join(' -V ')} \
            -o ${key}_genotyping.vcf.gz
        """
    }


    genotyped.toList().transpose().toList().set { comb_input }


    process combineChrVCFs {
        input:
            set val(keys), file(vcf), file(idx) from comb_input
        output:
            set val('all'), file('all.vcf.gz'), file('all.vcf.gz.tbi') into hardfilters

        publishDir "${params.outdir}/genotype", mode: 'copy'


        script:
        """
        java -Xmx23g -cp /usr/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
            -R $reference \
            -V ${vcf.join(' -V ')} \
            -out all.vcf.gz -assumeSorted
        """
    }



    hardfilters.into { hardfilters_snp; hardfilters_indel }


    process hardfilters_snp {
        input:
            set val(key), file(vcf), file(index) from hardfilters_snp
            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
        output:
            set file('*SNP*vcf'), file('*filtered_snps*vcf')

        publishDir "${params.outdir}/genotype", mode: 'copy'


        script:
        """
        java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T SelectVariants -R $reference \
            -V $vcf -selectType SNP -o ${key}_raw_SNP_1.vcf

        java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V ${key}_raw_SNP_1.vcf \
            --filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
            --filterExpression "FS > 60.0" --filterName "FS_greater_than_60_filter" \
            --filterExpression "SOR > 3.0" --filterName "SOR_greater_than_3_filter" \
            --filterExpression "MQ < 40.0" --filterName "MQ_less_than_40_filter" \
            --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum_less_than_-12.5_filter" \
            --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
            -o ${key}_filtered_snps_1.vcf \

        vcftools --vcf ${key}_filtered_snps_1.vcf --keep-filtered PASS --out ${key}_pass_SNP_1 --remove-filtered-geno-all  \
            --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992
         """
    }


    process hardfilters_indel {
        input:
            set val(key), file(vcf), file(index) from hardfilters_indel
            set file(reference), file(refindex), file(refdict) from Channel.value([reference, referenceFaIndex, referenceDict])
        output:
            set file('*INDEL*vcf'), file('*filtered_indels*vcf')

        publishDir "${params.outdir}/genotype", mode: 'copy'


        script:
        """
        java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T SelectVariants -R $reference \
            -V $vcf -selectType INDEL -o ${key}_raw_INDEL_1.vcf

        java -Xmx7g -jar /usr/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V ${key}_raw_INDEL_1.vcf \
            --filterExpression "QD < 2.0" --filterName "QD_less_than_2_filter" \
            --filterExpression "FS > 200.0" --filterName "FS_greater_than_60_filter" \
            --filterExpression "SOR > 10.0" --filterName "SOR_greater_than_3_filter" \
            --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum_less_than_-8_filter" \
            -o ${key}_filtered_indels_1.vcf

        vcftools --vcf ${key}_filtered_indels_1.vcf  --keep-filtered PASS --out ${key}_pass_INDEL_1 --remove-filtered-geno-all  \
            --remove-filtered-all --recode --recode-INFO-all --maf 0.00001 --max-maf 0.99992
        """
    }
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
    if ( params.onlyMap && !params.fastqDir ) {
        log.warn("You need to specify --fastqDir when doing the mapping (--onlyMap)")
        fatal_error = true
    }
    if ( params.onlyMap && params.bamDir ) {
        log.warn("You need to specify --fastqDir not --bamDir when doing only doing the mapping (--onlyMap)")
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
        --help
           Print this help message
        --fastqDir <Dir>
           Directory containing fastq samples (.fq.gz)
        --bamDir <Dir>
           Instead of --fastqDir, directory containing bam and indexed bam sample files (.bam, .bam.bai)
        --reference <file>
           Genome reference file (has to be indexed)
        --known <file>
           File with known sites for quality calibration.
        --outdir <dir>
           Directory for output files
        --onlyMap
           Only run the mapping steps
        --project
           Slurm project to run with
    """
}

def infoMessage() {
    log.info """\
*** K9 WGS Pipeline ***
Configuration environemnt:
    Out directory:   $params.outdir
    Fastq directory: $params.fastqDir
    Bam directory:   $params.bamDir
    Reference:       $params.reference
    OnlyMap?:        $params.onlyMap
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
