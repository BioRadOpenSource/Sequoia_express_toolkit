process count_rna {
    label 'low_cpu'
    label 'low_memory'
    tag "countLongRNA on $name"
    publishDir "${params.outDir}/Sample_Files/$name/RNACounts", mode: 'copy'

    input:
    tuple val(name), path(bam)
    path(longRNAgtfFile)
    val(geneId)
    
    output:
    tuple val(name),path ("gene_counts_longRNA.$name"), emit: counts_ch
    path("*.$name"), emit: report_longRNACounts
    path("*.featureCounts.bam")

    script:
    strand = params.reverseStrand ? "-s 2" : "-s 1"
    just_bam = bam[0]
    paired = "-p" 
    if (params.seqType == "SE") {
	paried = ""
	}
    """
    featureCounts $paired -T $task.cpus --primary -M -t exon -g $geneId $strand -Q $params.minMapqToCount \
    -a $longRNAgtfFile \
    -o ./gene_counts_longRNA \
    -R BAM $just_bam
    mv gene_counts_longRNA.summary gene_counts_longRNA.summary.$name
    mv gene_counts_longRNA gene_counts_longRNA.$name
    """
}

process calcRPKMTPM {
    label 'low_memory'
    label 'low_cpu'
    tag "calcRPKMTPM on $name"
    publishDir "${params.outDir}/Sample_Files/$name/calcRPKMTPM", mode: 'copy'
    input:
    tuple val(name), path(counts)

    output:
    tuple val(name),path('gene_counts_rpkmtpm.txt'), emit: rpkm_tpm_ch

    script:
    """
    python3 /opt/biorad/src/calc_rpkm_tpm.py $counts ./gene_counts_rpkmtpm.txt
    """
}

process thresholdResults{
        label 'low_memory'
        label 'micro_cpu'
        tag 'thresholdGenes'
        publishDir "${params.outDir}/Sample_Files/$name/RNACounts", mode:'copy'
        // take in user specified cutoff and type and generate appropriate report
        // should also include biotype 
        input:
        tuple val(name), path(rpkm) ,path(counts) 
        path(annoDirPath)
        output:
        path("Full_count_table.csv")
        path("Filter_count_table.csv")
        path("Filter_count_table.csv.$name"), emit: threshold_ch
        script:
        """
        mkdir -p ./out/
        mv $rpkm ./out/
        mv $counts ./out/gene_counts_longRNA
        mkdir -p ./tmp
        cp /opt/biorad/src/threshold_results.R ./tmp/threshold_results.R
        Rscript ./tmp/threshold_results.R "${params.minGeneType}" "${params.minGeneCutoff}" \$(readlink -f ./out) \$(readlink -f ./tmp)  \$(readlink -f $annoDirPath)
        cp Filter_count_table.csv Filter_count_table.csv.$name                
        """
}

