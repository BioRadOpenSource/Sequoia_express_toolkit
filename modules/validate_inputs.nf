process validateInputs {
    tag "Validation on $sample_id"
    publishDir "${params.outDir}/Sample_Files/$sample_id/validation", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    script:
    """
    python3 /opt/biorad/src/validate.py $reads 2>&1 > validation.log
    """
}
