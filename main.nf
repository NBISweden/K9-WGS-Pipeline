fastqFiles = Channel.fromFilePairs('test-data/test-data-tiny/R{1,2}.fq')
params.genomeFile = 'test-data/test-data-tiny/reference_chr38-1000000-1010000.fa'

process MapReads {
    input:
        set val(key), file(fastqs) from fastqFiles
        file(genomeFile) from Channel.fromPath(params.genomeFile)
    output:
        file("file.bam")

    publishDir "out"

    script:
    readGroup = "@RG\\tID:Sample_79162\\tSM:bar"
    """
    bwa index $genomeFile
    bwa mem -t $task.cpus -M -R \'$readGroup\' $genomeFile $fastqs |
        samtools sort --threads $task.cpus -m 4G > file.bam
    """
}
