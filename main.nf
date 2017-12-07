fastqFiles = Channel.fromFilePairs('test-data/test-data-tiny/*R{1,2}.fq.gz')

reference = file(params.reference)

process MapReads {
    input:
        set val(key), file(fastqs) from fastqFiles
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
