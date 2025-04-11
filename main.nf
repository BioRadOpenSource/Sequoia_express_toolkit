nextflow.enable.dsl = 2
def paramsWithUsage = readParamsFromJsonSettings()

// Constants
acceptableGenomes = ["rnor6", "hg38", "mm10","tair10","sacCer3","dm6","danRer11","ce11"]
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
genomeDirPath         = file(params.genomes[params.genome][params.spikeType].genomeDir)
annoDirPath           = file(params.genomes[params.genome][params.spikeType].annoDir)
geneId                = params.genomes[params.genome][params.spikeType].geneId
sjdbGTFFile           = file(params.genomes[params.genome][params.spikeType].sjdbGTFFile)
refFlatFile           = file(params.genomes[params.genome][params.spikeType].refFlatFile)
ribosomalIntervalFile = file(params.genomes[params.genome][params.spikeType].ribosomalIntervalFile)
longRNAgtfFile        = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)
sizesFile             = file(params.genomes[params.genome][params.spikeType].sizesFile)

//import modules for processes
include {debarcode}                                 from './modules/debarcode'
include {cutAdapt}                                  from './modules/cutadapt'
include {umiTagging; deduplication}                 from './modules/deduplication'
include {fastQc}                                    from './modules/fastqc'
include {count_rna; calcRPKMTPM; thresholdResults} from './modules/genecount'
include {picardAlignSummary}                        from './modules/picard'
include {assembleReport; combinedXLS; metaReport}   from './modules/report'
include {starAlign}                                 from './modules/star'
include {rename_tuple}                              from './modules/rename_tuple'
include {validateInputs}                            from './modules/validate_inputs'

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
summary['Seq Type'] = params.seqType
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
log.info summary.collect { k,v -> "${k.padRight(19)}: $v" }.join("\n")
log.info "----------------------------------------------------"

workflow{

    if(params.seqType =="SE"){
    	reads = params.reads+"*_R1*"
    	params.skipUmi =true
    }
    else{
    	reads = params.reads+"*{R1,R2}*"
    }

    read_files = Channel
        .fromFilePairs( reads, size:-1)
        .ifEmpty { exit 1, "Cannot find any reads in illumina format in dir: $reads\nIf not using R2 please use --seqType SE" }


    if(params.seqType =="SE"){
        rename_tuple(read_files)
        raw_reads = rename_tuple.out.reads_tuple
    }
    else{
    	raw_reads = read_files
    }
    // Begin Processing

    if (params.validateInputs) {
        validateInputs(raw_reads)
    }

    fastQc(raw_reads)
    // Only extract barcodes if umiAware
    if (!params.skipUmi && params.seqType=="PE") {
        debarcode(raw_reads)
        debarcoded_ch = debarcode.out.debarcoded_ch
        report_debarcode = debarcode.out.report_debarcode

    } else {
        debarcoded_ch = raw_reads
        report_debarcode = Channel.empty().ifEmpty([])
    }

    //trimming
    cutAdapt(debarcoded_ch)
    //alignment
    starAlign(cutAdapt.out.trimmed_ch, genomeDirPath, sjdbGTFFile)
    //alignment stats
    picardAlignSummary(starAlign.out.sorted_bam, refFlatFile, ribosomalIntervalFile)

    //process umis and deduplicate reads
    if (!params.skipUmi && params.seqType=="PE") {
        umiTagging(starAlign.out.sorted_bam)
        deduplication(umiTagging.out.dedup_in_ch)
        bamLong_ch = deduplication.out.bamLong_ch
        deduplication_report = deduplication.out.report_dedup

    } else {
        bamLong_ch = starAlign.out.sorted_bam
        deduplication_report = Channel.empty()
    }

    //feature counting 
    count_rna(bamLong_ch, longRNAgtfFile, geneId)

    //normalized read counts 
    calcRPKMTPM(count_rna.out.counts_ch)

    if(params.minGeneType != "none"){
            thresholdResults(
                calcRPKMTPM.out.rpkm_tpm_ch.combine(
                count_rna.out.counts_ch, by: 0),
                annoDirPath)
    }
    else{
            threshold_ch = Channel.empty().ifEmpty([])
    }

    //individual sample reports 
    assembleReport(
        annoDirPath,
        fastQc.out.report_fastqc.collect(),
        report_debarcode.collect().ifEmpty([]),
        cutAdapt.out.report_trim.collect(),
        starAlign.out.report_star.collect(),
        picardAlignSummary.out.report_picard.collect(),
        deduplication_report.collect().ifEmpty([]),
        count_rna.out.report_longRNACounts.collect()
        )

    //gene counts and output as spreadsheet 
    combinedXLS(count_rna.out.counts_ch.join(calcRPKMTPM.out.rpkm_tpm_ch))

    //generate batch report
    metaReport(
        starAlign.out.report_star.collect(), 
        picardAlignSummary.out.report_picard.collect(), 
        deduplication_report.collect().ifEmpty([])
        )


    //end workflow 
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
    def paramsContent = new File(params.settings).text
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
    nextflow run Sequoia_express_toolkit/main.nf  --outDir ./output/ --reads '~/read/express/' --genome hg38 --genomes_base ./genomes/

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
