process starAlign {
    label 'high_memory'
    label 'mid_cpu'
    tag "starAlign on $name"
    publishDir "${params.outDir}/Sample_Files/$name/star", mode: 'copy'


    input:
    set val(name), file(input) from trimmed_ch
    file genomeDirPath
    file sjdbGTFFile

    output:
    set val(name), file("Aligned.sortedByCoord.out.bam*") into umiTagging_ch, picardBam_ch
    file ("Unmapped.out.mate*")
    val name into meta_names_star
    file "Log.final.out.*" into meta_star, report_star

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