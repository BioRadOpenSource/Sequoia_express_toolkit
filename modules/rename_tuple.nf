process rename_tuple{
	input:
	set sid, file(fastq) from read_files
	output:
	set stdout, file(reads) into raw_reads_fastqc, raw_reads, raw_reads_validation
	script:
	reads = fastq
	"""
	echo ${sid} | sed "s/\\(_\\)R1.*//g" | xargs printf
	"""	
}