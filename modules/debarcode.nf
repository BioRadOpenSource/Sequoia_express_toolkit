process debarcode{
label 'micro_cpu'
label 'low_memory'
tag "debarcode DEAD on $sample_id"
publishDir "${params.outDir}/Sample_Files/$sample_id/debarcode", mode: 'copy'
input:
set sample_id, file(reads) from raw_reads
output:
set val(sample_id), file("*.fastq.gz") into debarcoded_ch
file("*stats.tsv") into report_debarcode  
script:
"""
export RUST_LOG=info	
dead -c /opt/biorad/src/2dcomplete.json -a DefaultParser -o ./ -i $sample_id $reads
"""
}