process cutAdapt {
    label 'low_cpu'
    label 'mid_memory'
    tag "cutAdapt on $name"
    publishDir "${params.outDir}/Sample_Files/$name/cutAdapt", mode: 'copy'

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path( "trimmed_*.fastq.gz"), emit: trimmed_ch
    path("trimlog.log.*"), emit: report_trim

    script:
	cutter = "-u 1"
	if(params.noTrim){
		cutter = ""
	}

	if (params.seqType == "SE") {
	read1 = reads[0]
    	"""
   	 cutadapt $cutter -m ${params.minBp} -j $task.cpus \
             -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
             -o trimmed_R1.fastq.gz $read1 1> trimlog.log
    	mv trimlog.log trimlog.log.$name
    	"""
	}
	else{
	//paired end 
	read1 = reads[0]
	read2 = reads[1] 
        cutter = "-u 1"
	if(params.skipUmi){
		cutter = cutter+" -U 8"
	}
	if(params.noTrim){
		cutter = ""
	}
	"""
    	cutadapt $cutter -m ${params.minBp} -j $task.cpus \
        -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
        -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz $read1 $read2 1> trimlog.log
    	mv trimlog.log trimlog.log.$name
    	"""

	}
}