fastqFiles = Channel.fromFilePairs('test-data/test-data-tiny/*R{1,2}.fq.gz')

reference = file(params.reference)

process MapReads {
    input:
        set val(key), file(fastqs) from fastqFiles
    output:
        file("file.bam")

    publishDir "out"

    script:
    readGroup = "@RG\\tID:Sample_79162\\tSM:bar"
    """
    bwa mem -t $task.cpus -M -R \'$readGroup\' $reference $fastqs |
        samtools sort --threads $task.cpus -m 4G > file.bam
    """
}
