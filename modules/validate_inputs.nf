process validateInputs {
    tag "Validation on $sample_id"
    publishDir "${params.outDir}/Sample_Files/$sample_id/validation", mode: 'copy'
    input:
    set sample_id, file(reads) from raw_reads_validation
    script:
    """
    python3 /opt/biorad/src/validate.py $reads 2>&1 > validation.log
    """
}
