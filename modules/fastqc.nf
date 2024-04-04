process fastQc {
    tag "FASTQC on $sample_id"
    label 'micro_cpu' 
    publishDir "${params.outDir}/Sample_Files/$sample_id/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set sample_id, file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results, report_fastqc

    script:
    """
    fastqc $reads
    """
}