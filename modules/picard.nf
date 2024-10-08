process picardAlignSummary {
    label 'low_memory'
    label 'micro_cpu'
    tag "picardAlignSummary on $name"
    publishDir "${params.outDir}/Sample_Files/$name/picardAlignSummary", mode: 'copy'

    input:
    tuple val(name), path(bams)
    path(refFlatFile)
    path(ribosomalIntervalFile)

    output:
    path('rna_metrics.txt.*'), emit: report_picard

    script:
    (bam, bai) = bams
    strand = params.reverseStrand ? "SECOND_READ_TRANSCRIPTION_STRAND" : "FIRST_READ_TRANSCRIPTION_STRAND"
    //check format of memory to avoid breaking picard tools command line vs config
    if(task.memory.toGiga() ==0){
	    picard_mem = task.memory.toBytes()
    }
    else{
	    picard_mem = task.memory.toGiga() 
    }
    """
    picard CollectRnaSeqMetrics -I $bam \
    -O rna_metrics.txt \
    -REF_FLAT $refFlatFile \
    -STRAND $strand \
    -RIBOSOMAL_INTERVALS $ribosomalIntervalFile \
    -Xmx${picard_mem}g
    cp rna_metrics.txt rna_metrics.txt.$name
    """
}