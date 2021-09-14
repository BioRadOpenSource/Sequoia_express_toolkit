
def paramsWithUsage = readParamsFromJsonSettings()

// Constants
acceptableGenomes = ["rnor6", "hg38", "mm10"]
allowedSpikes     = ["ercc"]

// Show help emssage
if (params.help){
    helpMessage(paramsWithUsage)
    exit 0
}

// Validate that genome is in correct set
if ( !acceptableGenomes.contains(params.genome) ) {
    log.error "$params.genome not in acceptable genomes $acceptableGenomes"
    exit 1
}

if (params.spikeType != "NONE" && !allowedSpikes.contains(params.spikeType)) {
    log.error "$params.spikeType is not in $allowedSpikes"
    exit 1
}

// Make all the genome related things file resources
genomeDirPath        = file(params.genomes[params.genome][params.spikeType].genomeDir)
annoDirPath          = file(params.genomes[params.genome][params.spikeType].annoDir)
geneId               = params.genomes[params.genome][params.spikeType].geneId
sjdbGTFFile          = file(params.genomes[params.genome][params.spikeType].sjdbGTFFile)
refFlatFile          = file(params.genomes[params.genome][params.spikeType].refFlatFile)
ribosomalIntervalFile = file(params.genomes[params.genome][params.spikeType].ribosomalIntervalFile)
longRNAgtfFile       = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)
sizesFile            = file(params.genomes[params.genome][params.spikeType].sizesFile)

// Check that inputs are viable / R1 only vs both
if (params.reads == "NOINPUT") {
    log.error "No input reads were supplied"
    exit 1
}


// Create Summary
def summary = [:]
summary['Run Name'] = workflow.runName
summary['Reads'] = params.reads
summary['Genome'] = params.genome
summary['Spike Type'] = params.spikeType
summary['Reference Dir'] = params.genomes_base
summary['Annotations Dir'] = annoDirPath
summary['Skip UMI?'] = params.skipUmi
summary['Min MAPQ To Count'] = params.minMapqToCount
summary['Output Dir'] = params.outDir
summary['Trace Dir'] = params.tracedir
/*summary['Max Cores'] = task.cpus*/
summary['geneId'] = geneId
summary['sjdb GTF File'] = sjdbGTFFile
summary['ref Flat File'] = refFlatFile
summary['Ribosomal Intervals'] = ribosomalIntervalFile
summary['Long RNA GTF File'] = longRNAgtfFile
summary['Sizes File'] = sizesFile
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
log.info bioradHeader()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"



reads = params.reads+"*{R1,R2}*"

Channel
    .fromFilePairs( reads, size:-1)
    .ifEmpty { exit 1, "Cannot find any reads in illumina format in dir: $reads\nIf not using R2 please use --seqType SE" }
    .set { read_files}

read_files.into{ raw_reads_fastqc; raw_reads; raw_reads_validation; raw_reads_dead}
// Begin Processing

if (params.validateInputs) {
    process validateInputs {
        tag "Validation on $sample_id"
        publishDir "${params.outDir}/$sample_id/validation", mode: 'copy'

        input:
        set sample_id, file(reads) from raw_reads_validation

        script:
        """
        python3 /opt/biorad/src/validate.py $reads 2>&1 > validation.log
        """
    }
}
process fastQc {
    tag "FASTQC on $sample_id"
    publishDir "${params.outDir}/$sample_id/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set sample_id, file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results, report_fastqc

    script:
    """
    fastqc $reads
    """
}
//TODO: this will need an update based on where they put the UMI
// Only extract barcodes if umiAware
if (!params.skipUmi) {
    
    process debarcode{
	label 'mid_cpu'
	label 'mid_mem'
	tag "debarcode DEAD on $sample_id"
	publishDir "${params.outDir}/$sample_id/debarcode", mode: 'copy'

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
} else {
    debarcoded_ch = raw_reads
    report_debarcode = Channel.empty()
}

process cutAdapt {
    label 'mid_cpu'
    tag "cutAdapt on $name"
    publishDir "${params.outDir}/$name/cutAdapt", mode: 'copy'

    input:
    set val(name), file(reads) from debarcoded_ch

    output:
    set val(name), file( "trimmed_*.fastq.gz") into trimmed_ch
    file "trimlog.log.*" into report_trim

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
        cutter = "-u 1 -U 8"
	"""
    	cutadapt $cutter -m ${params.minBp} -j $task.cpus \
        -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
        -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz $read1 $read2 1> trimlog.log
    	mv trimlog.log trimlog.log.$name
    	"""

	}
}

process starAlign {
    label 'high_memory'
    label 'mid_cpu'
    tag "starAlign on $name"
    publishDir "${params.outDir}/$name/star", mode: 'copy'


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

process picardAlignSummary {
    label 'low_memory'
    tag "picardAlignSummary on $name"
    publishDir "${params.outDir}/$name/picardAlignSummary", mode: 'copy'

    input:
    set val(name), file(bams) from picardBam_ch
    file refFlatFile
    file ribosomalIntervalFile

    output:
    val name into meta_names_picard
    file 'rna_metrics.txt.*' into meta_picard, report_picard

    script:
    (bam, bai) = bams
    strand = params.reverseStrand ? "SECOND_READ_TRANSCRIPTION_STRAND" : "FIRST_READ_TRANSCRIPTION_STRAND"
    """
    java -jar /opt/picard/picard.jar CollectRnaSeqMetrics I=$bam \
    O=rna_metrics.txt \
    REF_FLAT=$refFlatFile \
    STRAND=$strand \
    RIBOSOMAL_INTERVALS=$ribosomalIntervalFile
    cp rna_metrics.txt rna_metrics.txt.$name
    """
}

if (!params.skipUmi) {
    process umiTagging {
        label 'mid_cpu'
        label 'mid_memory'
        tag "umiTagging on $name"
        publishDir "${params.outDir}/$name/umiTagging", mode: 'copy'
        input:
        set val(name) , file(bams) from umiTagging_ch

        output:
        set val(name), file("Aligned.sortedByCoord.tagged.bam*") into dedup_in_ch
	set val(cutname), file("cutAligned.sortedByCoord.tagged.bam*") into cutdedup_in_ch
        
        script:
        (bam, bai) = bams
	cutname = "cut$name"
        """
        samtools idxstats $bam | cut -f 1 | grep -E 'chr|ERCC-*' > ./Aligned.sortedByCoord.idxstats.txt
        #samtools idxstats $bam | cut -f 1 | uniq > ./Aligned.sortedByCoord.idxstats.txt
        mkdir -p ./tmp/
        python3 /opt/biorad/src/tagBamFile.py $bam ./Aligned.sortedByCoord.idxstats.txt ./tmp/ $task.cpus
        sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep -E 'chr|ERCC-*')
        #sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep .bam)
        sambamba index -t $task.cpus ./Aligned.sortedByCoord.tagged.bam
        rm -r ./tmp/
	
	cat <(samtools view -H $bam) <(samtools view $bam | cut -c2-) | samtools view -Sb > cut.bam
	sambamba index -t $task.cpus cut.bam
	samtools idxstats cut.bam | cut -f 1 | grep -E 'chr|ERCC-*' > ./Aligned.sortedByCoord.idxstats.txt
	mkdir -p ./tmp/
	python3 /opt/biorad/src/tagBamFile.py cut.bam ./Aligned.sortedByCoord.idxstats.txt ./tmp/ $task.cpus
	sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep -E 'chr|ERCC-*')
	sambamba index -t $task.cpus ./cutAligned.sortedByCoord.tagged.bam
        rm -r ./tmp/
        """

    }
    process deduplication {
	label 'high_memory'
	label 'mid_cpu'
	tag "rumi dedup on $name"
	publishDir "${params.outDir}/$name/dedup", mode: 'copy'
	
	input:
	set val(name), file(bams) from dedup_in_ch
	output:	
	set val(name), file("rumi_dedup.sort.bam*") into BamLong_ch
	file 'dedup.log.*' into report_dedup, meta_dedup

	script:
	(bam, bai) = bams
	"""
	export RUST_LOG=info 
	rumi --is_paired $bam --output rumi_dedup.bam --umi_tag XU >dedup.log
	sambamba sort -t $task.cpus ./rumi_dedup.bam -o rumi_dedup.sort.bam 
	sambamba index -t $task.cpus ./rumi_dedup.sort.bam
	#printf "unique_input_umi: " >> ./dedup.log; samtools view $bam | cut -f1 | sort -u | wc -l >> ./dedup.log
	printf "unique_input_reads: " >> ./dedup.log; samtools view $bam | cut -f1,10 | wc -l >> ./dedup.log
	printf "unique_umi: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f1 | sort -u | wc -l >> ./dedup.log
	printf "unique_output_reads: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f1,10 | sort -u | wc -l >> ./dedup.log
	cp dedup.log dedup.log.$name
	"""
	}
	
	process cutdeduplication {
	// dead end process - does not pipe into subsequent processes
	label 'high_memory'
	label 'mid_cpu'
	tag "rumi dedup on $name"
	publishDir "${params.outDir}/$name/dedup", mode: 'copy'
	
	input:
	set val(name), file(bams) from cutdedup_in_ch
	
	output:	
	set val(name), file("rumi_dedup.sort.bam*")
	file 'dedup.log.*'

	script:
	(bam, bai) = bams
	"""
	export RUST_LOG=info 
	rumi --is_paired $bam --output rumi_dedup.bam --umi_tag XU >dedup.log
	sambamba sort -t $task.cpus ./rumi_dedup.bam -o rumi_dedup.sort.bam 
	sambamba index -t $task.cpus ./rumi_dedup.sort.bam
	#printf "unique_input_umi: " >> ./dedup.log; samtools view $bam | cut -f1 | sort -u | wc -l >> ./dedup.log
	printf "unique_input_reads: " >> ./dedup.log; samtools view $bam | cut -f1,10 | wc -l >> ./dedup.log
	printf "unique_umi: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f1 | sort -u | wc -l >> ./dedup.log
	printf "unique_output_reads: " >> ./dedup.log; samtools view ./rumi_dedup.sort.bam | cut -f1,10 | sort -u | wc -l >> ./dedup.log
	cp dedup.log dedup.log.$name
	"""
	}
	
} else {
    umiTagging_ch.into { BamLong_ch} 
    report_dedup = Channel.empty()
    meta_dedup = Channel.empty()
}


process count_rna {
    label 'mid_cpu'
    label 'low_memory'
    tag "countLongRNA on $name"
    publishDir "${params.outDir}/$name/RNACounts", mode: 'copy'

    input:
    set val(name), file(bam) from BamLong_ch
    file longRNAgtfFile
    
    output:
    val(name) into (counts_name, xls_name, threshold_name)
    file "gene_counts_longRNA*" into (counts_ch, counts_xls, count_threshold_ch)
    file "*.$name" into report_longRNACounts

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
    cp gene_counts_longRNA.summary gene_counts_longRNA.summary.$name
    cp gene_counts_longRNA gene_counts_longRNA.$name
    """
}

process calcRPMKTPM {
    label 'low_memory'
    label 'low_cpu'
    tag "calcRPMKTPM on $name"
    publishDir "${params.outDir}/$name/calcRPMKTPM", mode: 'copy'
    input:
    val name from counts_name
    file(counts) from counts_ch

    output:
    file 'gene_counts_rpkmtpm.txt' into rpkm_tpm_ch, normalize_xls, rpkm_threshold_ch
    val name into thresh_ch

    script:
    """
    python3 /opt/biorad/src/calc_rpkm_tpm.py gene_counts_longRNA ./gene_counts_rpkmtpm.txt
    """

}
if(params.minGeneType != "none"){
        process thresholdResults{
                label 'low_memory'
                tag 'thresholdGenes'
                publishDir "${params.outDir}/$name/RNACounts", mode:'copy'

                // take in user specified cutoff and type and generate appropriate report
                // should also include biotype 
                input:
		val name from thresh_ch
		file ("./out/") from rpkm_threshold_ch 
		file ("./out/") from count_threshold_ch	

                output:
                file "Full_count_table.csv"
                file "Filter_count_table.csv"
                file "Filter_count_table.csv.$name" into threshold_ch

                script:
                """
		mkdir -p ./tmp
                cp /opt/biorad/src/threshold_report.R ./tmp/threshold_report.R
                Rscript ./tmp/threshold_report.R "${params.minGeneType}" "${params.minGeneCutoff}" \$(readlink -f ./out) \$(readlink -f ./tmp)  \$(readlink -f $annoDirPath)
		cp Filter_count_table.csv Filter_count_table.csv.$name                
                """
        }
}
else{
        threshold_ch = Channel.empty()
}


process assembleReport {
    label 'low_memory'
    tag "assembleReport"
    publishDir "${params.outDir}/reports", mode: 'copy' // TODO: Filter down the outputs since so much stuff will be in this dir

    input:
    file annoDirPath
    file(fastqc: "out/fastqc/*") from report_fastqc.collect()
    file("out/debarcode/*") from report_debarcode.collect().ifEmpty([]) // optional
    file("out/cutAdapt/*") from report_trim.collect()
    file("out/star/*") from report_star.collect() 
    file("out/star/*") from report_picard.collect() // Goes into star for reasons
    file("out/umitools/*") from report_dedup.collect().ifEmpty([]) // optional
    file("out/counts/*") from report_longRNACounts.collect()
	
    output:
    file '*_htmlReport.html'
    file '*_pdfReport.pdf'
    file '*_csvReport.csv'
    
    script:
    """
    mkdir -p ./tmp
    cp /opt/biorad/src/htmlReport.R ./tmp/htmlReport.R
    cp /opt/biorad/src/pdfReport.R ./tmp/pdfReport.R
    cp /opt/biorad/src/csvReport.R ./tmp/csvReport.R
    Rscript /opt/biorad/src/generateRmdReport.R \$(readlink -f ./out) \$(readlink -f ./tmp)  \$(readlink -f $annoDirPath)
    cp ./tmp/*_htmlReport.html ./
    cp ./tmp/*_pdfReport.pdf ./
    cp ./tmp/*_csvReport.csv ./
    """
}
process combinedXLS{
	label 'low_memory'
	tag "countsAsXls"
	publishDir "${params.outDir}/$name/calcRPMKTPM", mode:'copy'

	input:
	file rpkm from normalize_xls
	val(name) from xls_name 
        file(count_file) from counts_xls

	output:
	file "readcount_report.xlsx"
	
	script:
	"""
	python3 /opt/biorad/src/converge_xls.py gene_counts_longRNA $rpkm 
	"""

}
process metaReport{
	tag "Overall Batch Summary"
	publishDir "${params.outDir}/", mode:'copy'
	// generate a high level summary of batch run
	
	input:
	file("out/star/") from meta_star.collect()
	file("out/picard/") from meta_picard.collect()
	file("out/dedup/") from meta_dedup.collect().ifEmpty([])

	output:
	file 'batch_summary.csv'
	file 'batch_summary.html'
	file 'batch_summary.pdf'
	val "done" into clean_workspace
	script:	
	"""
	mkdir -p tmp/
	cp /opt/biorad/src/batch_html.R ./tmp/batch_html.R
	cp /opt/biorad/src/batch_pdf.R ./tmp/batch_pdf.R
	Rscript /opt/biorad/src/meta_report.R \$(readlink -f ./out) \$(readlink -f ./tmp)
	cp ./tmp/batch_html.html ./batch_summary.html
	cp ./tmp/batch_pdf.pdf ./batch_summary.pdf
	"""
}
if(params.tidy == true){
	process tidy_up{
		input: val x from clean_workspace

		
		"""
		 rm -r ${workflow.workDir}
		""" 

	}
}

/* Helper Functions */
def readParamsFromJsonSettings() {
    List paramsWithUsage
    try {
        paramsWithUsage = tryReadParamsFromJsonSettings()
    } catch (Exception e) {
        println "Could not read parameters settings from json. $e"
        pramsWithUsage = Collections.emptyMap()
    }
    return paramsWithUsage
}

def tryReadParamsFromJsonSettings() throws Exception {
    def paramsContent = new File(config.params_description.path).text
    def paramsWithUsage = new groovy.json.JsonSlurper().parseText(paramsContent)
    return paramsWithUsage.get('parameters')
}

String prettyFormatParamsWithPaddingAndIndent(List paramGroup, Integer padding=2, Integer indent=4) {
    def fields = ["name", "usage", "type", "default_value", "pattern", "choices"]
    def maxFields = fields.collectEntries { String field ->
        [(field): paramGroup.collect {
            def val = it.get(field)
            val ? val.toString().size() : 1
        }.max()]
    }
    def formatter = {param -> 
        sprintf("%${indent}s%-${maxFields.name}s (%-${maxFields.type}s) %-${maxFields.default_value}s %-${maxFields.usage}s %-${maxFields.choices}s %-${maxFields.pattern}s\n", "",
                                param.name ?: "", param.type ?: "",  param.default_value ?: "", param.usage ?: "", param.choices ?: "", param.pattern ?: "")
    }
    def requiredParamsFormattedList = paramGroup.sort { it.name }.findAll { it.required }.collect { Map param -> formatter(param) }
    def optionalParamsFormattedList = paramGroup.sort { it.name }.findAll { !it.required }.collect { Map param -> formatter(param) }
    return String.format("REQUIRED:\n%s\nOPTIONAL:\n%s\n", requiredParamsFormattedList.join(), optionalParamsFormattedList.join())
}

def helpMessage(paramsWithUsage) {
    def helpMessage = String.format(
    """\
    %s
    
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run . --reads './tests/*_R{1,2}.fastq.gz' --genome hg38 --outDir /data/out --skipUmi --genomes_base /mnt/genome-annotations

    Args:

    %s
    """.stripIndent(), bioradHeader(), prettyFormatParamsWithPaddingAndIndent(paramsWithUsage, 8, 4))
    log.info helpMessage
}

def bioradHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    return """
    ${c_reset}
    ${c_green}
    ${c_green}/-----------------------------------------------------------\\ 
    ${c_green}| __________.__                __________             .___  |
    ${c_green}|  \\_____   \\__|____           \\______   \\____      __| _/  |
    ${c_green}|   |  |  _/  |/  _ \\   ______   |     _/\\__  \\   / __ |    |
    ${c_green}|   |  |   \\  (  <_> ) /_____/   |  |   \\ / __ \\_/ /_/ |    |
    ${c_green}|   |____  /__|\\____/            |__|_  /(____  /\\____ |    |
    ${c_green}|        \\/                           \\/      \\/      \\/    |
    ${c_green}\\___________________________________________________________/
    ${c_reset}
    """.stripIndent()
}
