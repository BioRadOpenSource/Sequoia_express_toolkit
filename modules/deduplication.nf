process umiTagging {
    label 'low_cpu'
    label 'mid_memory'
    tag "umiTagging on $name"
    //publishDir "${params.outDir}/$name/umiTagging", mode: 'copy'
    input:
    tuple val(name) , path(bams)
    output:
    tuple val(name), path("Aligned.sortedByCoord.tagged.bam*"), emit: dedup_in_ch
    
    script:
    (bam, bai) = bams
    """
    samtools idxstats $bam | cut -f 1 | grep -E 'chr|ERCC-*' > ./Aligned.sortedByCoord.idxstats.txt
    #samtools idxstats $bam | cut -f 1 | uniq > ./Aligned.sortedByCoord.idxstats.txt
    mkdir -p ./tmp/
    python3 /opt/biorad/src/tagBamFile.py $bam ./Aligned.sortedByCoord.idxstats.txt ./tmp/ $task.cpus
    sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep -E 'chr|ERCC-*')
    #sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep .bam)
    sambamba index -t $task.cpus ./Aligned.sortedByCoord.tagged.bam
    rm -r ./tmp/
    """
}

process deduplication {
label 'high_memory'
label 'mid_cpu'
tag "rumi dedup on $name"
publishDir "${params.outDir}/Sample_Files/$name/dedup", mode: 'copy'
	
input:
tuple val(name), path(bams)
output:	
tuple val(name), path("rumi_dedup.sort.bam*"), emit:  bamLong_ch
path('dedup.log.*'), emit: report_dedup
script:
(bam, bai) = bams
"""
export RUST_LOG=info 
rumi --is_paired $bam --output rumi_dedup.bam --umi_tag XU #>dedup.log
sambamba sort -t $task.cpus ./rumi_dedup.bam -o rumi_dedup.sort.bam 
sambamba index -t $task.cpus ./rumi_dedup.sort.bam
touch dedup.log
printf "Reads In: " >>./dedup.log; samtools view $bam |cut -f1| sort -u | wc -l >> ./dedup.log
printf "Reads Out: " >>./dedup.log; samtools view ./rumi_dedup.sort.bam |cut -f1| sort -u | wc -l >> ./dedup.log  
#printf "unique_input_umi: " >> ./dedup.log; samtools view $bam | cut -f16 | sort -u | wc -l >> ./dedup.log
printf "No Mate: " >> ./dedup.log;samtools view -F 2 $bam | cut -f1,16 | sort -u| wc -l >> ./dedup.log
printf "unique_input_reads: " >> ./dedup.log; samtools view -f 2 $bam | cut -f1,16 | sort -u| wc -l >> ./dedup.log
printf "unique_umi: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f16 | sort -u | wc -l >> ./dedup.log
printf "unique_output_reads: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f1,16 | sort -u | wc -l >> ./dedup.log
cp dedup.log dedup.log.$name
"""
}