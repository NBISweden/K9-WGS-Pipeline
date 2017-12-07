fastqFiles = Channel.fromFilePairs('test-data/test-data-tiny/*R{1,2}.fq.gz')

reference = file(params.reference)
known     = file(params.known)


fastqFiles.into { fastq_qc; fastq_bwa }

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
    output:
        set file("file.bam"), file("file.bam.bai") into mapped_reads


    publishDir "out"


    script:
    readGroup = "@RG\\tID:Sample_79162\\tSM:bar\\tPL:ILLUMINA"
    """
    bwa mem -t $task.cpus -M -R \'$readGroup\' $reference $fastqs |
        samtools sort --threads $task.cpus -m 4G > file.bam
    samtools index file.bam
    """
}


process gatk_realign {
    input:
        set file("file.bam"), file("file.bam.bai") from mapped_reads

    output:
        file('name.intervals')
        set file('file.realigned.bam'), file('file.realigned.bai') into realigned_reads


    publishDir 'out'


    script:
    """
    java -jar /usr/GenomeAnalysisTK.jar \
        -I file.bam \
        -R $reference \
        -T RealignerTargetCreator \
        -o name.intervals

    java -jar /usr/GenomeAnalysisTK.jar \
        -I file.bam \
        -R $reference \
        -T IndelRealigner \
        -targetIntervals name.intervals \
        -o file.realigned.bam
    """
}


process mark_duplicates {
    input:
        set file('file.bam'), file('file.bam.bai') from realigned_reads

    output:
        set file('file.realigned.marked.bam'), file('file.realigned.marked.bai') into marked_reads


    publishDir 'out'


    script:
    """
    picard MarkDuplicates.jar \
        INPUT=file.bam \
        OUTPUT=file.realigned.marked.bam \
        METRICS_FILE=file.realigned.marked.metrics \
        VALIDATION_STRINGENCY=LENIENT

    picard BuildBamIndex.jar \
        INPUT=file.realigned.marked.bam \
        VALIDATION_STRINGENCY=LENIENT
    """
}


process quality_recalibration {
    input:
        set file('file.bam'), file('file.bam.bai') from marked_reads
    output:
        file('file.recal_data.table')
        file('file.post_recal_data.table')
        file('file.recalibrated.bam')
        file('file.recalibration_plots.pdf')


    publishDir 'out'


    script:
    """
    # Step 1 - We need a list of known sites otherwise, GATK will think all the
    #          real SNPs in our data are errors. Failure to remove real SNPs
    #          from the recalibration will result in globally lower quality
    #          scores
    java -jar /usr/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -I file.bam \
        -R $reference \
        -knownSites $known \
        -o file.recal_data.table \
        -rf BadCigar

    # Step 2
    java -jar /usr/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -I file.bam \
        -R $reference \
        -knownSites $known \
        -BQSR file.recal_data.table \
        -o file.post_recal_data.table

    # Step 3 - Create before and after plots
    java -jar /usr/GenomeAnalysisTK.jar \
        -T AnalyzeCovariates \
        -R $reference \
        -before file.recal_data.table \
        -after file.post_recal_data.table \
        -plots file.recalibration_plots.pdf

    # Step 4 - Base calibration
    java -jar /usr/GenomeAnalysisTK.jar \
        -T PrintReads \
        -I file.bam \
        -R $reference \
        -BQSR file.recal_data.table \
        -o file.recalibrated.bam
    """
}
