process rename_tuple{
	input:
	tuple val(sid), path(fastq)
	output:
	tuple stdout, path(reads), emit: reads_tuple
	script:
	reads = fastq
	"""
	echo ${sid} | sed "s/\\(_\\)R1.*//g" | xargs printf
	"""	
}