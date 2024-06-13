process starAlign {
    label 'high_memory'
    label 'mid_cpu'
    tag "starAlign on $name"
    publishDir "${params.outDir}/Sample_Files/$name/star", mode: 'copy'


    input:
    tuple val(name), path(input) 
    path(genomeDirPath)
    path(sjdbGTFFile)

    output:
    tuple val(name), path("Aligned.sortedByCoord.out.bam*"), emit: sorted_bam
    path ("Unmapped.out.mate*")
    path( "Log.final.out.*"), emit: report_star

    script:
           """
    STAR \
        --readFilesIn $input \
        --readFilesCommand zcat \
        --genomeDir $genomeDirPath \
        --runThreadN $task.cpus \
        --sjdbGTFfile $sjdbGTFFile \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMatchNmin 15 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 \
        --outMultimapperOrder Random \
        --runRNGseed 1234 \
        --outFileNamePrefix ./ > star_log.txt 2>&1
    rm -rf _STARgenome
    sambamba index -t $task.cpus Aligned.sortedByCoord.out.bam
    cp Log.final.out Log.final.out.$name
    """ 
}