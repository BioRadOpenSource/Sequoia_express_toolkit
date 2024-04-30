process fastQc {
    tag "FASTQC on $sample_id"
    label 'micro_cpu' 
    publishDir "${params.outDir}/Sample_Files/$sample_id/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: report_fastqc

    script:
    """
    fastqc $reads
    """
}