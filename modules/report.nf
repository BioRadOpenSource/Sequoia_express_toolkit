process assembleReport {
    label 'low_memory'
    label 'micro_cpu'
    tag "assembleReport"
    publishDir "${params.outDir}/report", mode: 'copy' // TODO: Filter down the outputs since so much stuff will be in this dir

    input:
    path(annoDirPath)
    path(fastqc_collect), 		stageAs: "out/fastqc/*"
    path(debarcode_collect), 	stageAs: "out/debarcode/*"  // optional
    path(cutadapt_collect), 	stageAs: "out/cutAdapt/*" 
    path(star_collect), 		stageAs: "out/star/*"
    path(picard_collect), 		stageAs: "out/star/*"		//Goes into star for reasons
    path(dedup_collect), 		stageAs: "out/umitools/*"  	// optional
    path(count_collect), 		stageAs: "out/counts/*"
	
    output:
    path('*_htmlReport.html')
    path('*_pdfReport.pdf')
    path('*_csvReport.csv')
    
    script:
    """
    mkdir -p ./tmp
    cp /opt/biorad/src/htmlReport.R ./tmp/htmlReport.R
    cp /opt/biorad/src/pdfReport.R ./tmp/pdfReport.R
    cp /opt/biorad/src/csvReport.R ./tmp/csvReport.R
    Rscript /opt/biorad/src/generateRmdReport.R \$(readlink -f ./out) \$(readlink -f ./tmp)  \$(readlink -f $annoDirPath) ${params.genome} ${params.spikeType}
    cp ./tmp/*_htmlReport.html ./
    cp ./tmp/*_pdfReport.pdf ./
    cp ./tmp/*_csvReport.csv ./
    """
}
process combinedXLS{
	label 'low_memory'
	label 'micro_cpu'
	tag "countsAsXls on $name"
	publishDir "${params.outDir}/Sample_Files/$name/calcRPKMTPM", mode:'copy'

	input:
	tuple val(name), path(count_file), path(rpkm)

	output:
	path("readcount_report.xlsx")
	
	script:
	"""
	python3 /opt/biorad/src/converge_xls.py $count_file $rpkm 
	"""

}
process metaReport{
	tag "Overall Batch Summary"
	label 'micro_cpu'
	publishDir "${params.outDir}/report", mode:'copy'
	// generate a high level summary of batch run
	
	input:
	path(star_collect), stageAs: "out/star/*"
	path(picard_collect), stageAs: "out/picard/*"
	path(dedup_collect), stageAs: "out/dedup/*"

	output:
	path('batch_summary.csv')
	path('batch_summary.html')
	path('batch_summary.pdf')
	
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
