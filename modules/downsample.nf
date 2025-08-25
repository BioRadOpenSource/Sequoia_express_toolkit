/*
Downsample FASTQ files to a specific number of reads
*/

process downsample {
    tag "downsample on $sample_id"
    container "bioraddbg/omnition-core-dev:v1.0.0-alpha.98"
    publishDir "${params.outDir}/downsampled_fastq", pattern: '*.fastq.gz', mode: 'copy', overwrite: true
    memory '15 GB'
    cpus '8'

    input:
    tuple val(sampleId), path(reads, stageAs: 'raw/*')

    output:
    tuple val(sample_id), path("*.fastq.gz") emit: fastq

    script:
    """
    # Run downsample script
    downsample.py --I ${reads[0]} ${reads[1]} --O ./ --size 2M
    """
}