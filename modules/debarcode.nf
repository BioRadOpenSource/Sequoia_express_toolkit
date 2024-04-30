process debarcode{
label 'micro_cpu'
label 'low_memory'
tag "debarcode DEAD on $sample_id"
publishDir "${params.outDir}/Sample_Files/$sample_id/debarcode", mode: 'copy'
input:
tuple val(sample_id), path(reads)
output:
tuple val(sample_id), path("*.fastq.gz"), emit: debarcoded_ch
path("*stats.tsv"), emit: report_debarcode  
script:
"""
export RUST_LOG=info	
dead -c /opt/biorad/src/2dcomplete.json -a DefaultParser -o ./ -i $sample_id $reads
"""
}