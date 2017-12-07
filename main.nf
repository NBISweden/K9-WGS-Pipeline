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
