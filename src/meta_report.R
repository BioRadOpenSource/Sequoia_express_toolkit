library(knitr)
library(rmarkdown)
#get dir from command line arg

args = commandArgs(trailingOnly=T)

outPath = args[1]
#go recursive to get all star outfiles 

alignmentDir <- unique(dirname(list.files(outPath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="rna_metrics.txt.*")))
alignmentDirExists <- length(alignmentDir) >= 1
write(paste("alignmentDirExists: ", alignmentDirExists), stderr())

starDir <- unique(dirname(list.files(outPath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="Log.final.out.*")))
starDirExists <- length(starDir) >= 1
write(paste("logDirExists: ", starDirExists), stderr())

files = list.files(alignmentDir)
metrics = lapply(files, function(x) read.table(paste0(alignmentDir,"/",x), header=T, nrows=1,fill=T))
for(l in 1:length(metrics)){
	m = metrics[[l]]
	m[grep("PCT", names(m))]<- paste0((m[grep("PCT", names(m))]*100),"%")
	metrics[[l]] = m
}


metrics = lapply(metrics, function(x) data.frame(names(x), as.character(x[1,]), stringsAsFactors=F))
sample_names = gsub("rna_metrics.txt.","",files)
print(metrics)
Metric = metrics[[1]][,1]
#rename cols of each dataframe
align_frame = data.frame(Metric, stringsAsFactors=F)
for(i in 1:length(metrics)){
	align_frame = cbind(align_frame, metrics[[i]][,2])
}
colnames(align_frame) = c("Metric", sample_names)


files = list.files(starDir)
star_names = gsub("Log.final.out.","",files)
reports = lapply(files, function(x) read.table(paste0(starDir,"/",x), header=T, sep="|",fill=T))
labels = reports[[1]][,1]
reports = lapply(reports, function(x) x[,2] <-gsub("\t","", x[,2]))
report_frame = data.frame(labels)
for(i in 1:length(reports)){
	report_frame = cbind(report_frame, reports[[i]])
}
colnames(report_frame) = c("Metric", star_names)

#deduplication
dedupDir <- unique(dirname(list.files(outPath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="dedup.log.*")))
dedupDirExists <- length(dedupDir) == 1
write(paste("dedupDirExists: ", dedupDirExists), stderr())
if(dedupDirExists){
files = list.files(dedupDir)
dedup_frame = data.frame(Metrics=c("Total input alignments",
				   "Total output alignments",
				   "Unique UMIs observed",
				   "Average UMIs per position",
				   "Maximum UMIs per position",
				   "Unique Input Reads",
				   "Unique Output Reads",
				   "% PCR Duplicates"
				))
for( f in files){
	write(paste0("working on reading file ",f), stderr())
	file_loc = paste0(dedupDir,"/",f)
	umisObserved <- as.numeric(system(paste('grep -F "#umis"', file_loc, "| cut -d' ' -f5"), intern=T))
	inputAlignments <- as.numeric(system(paste('grep "Input Reads"', file_loc, "| cut -d' ' -f7|sed 's/,//g'"), intern=T))
	outputAlignments <- as.numeric(system(paste('grep "reads out"', file_loc, "| cut -d: -f4"), intern=T))
	meanUmiPerPos <- as.numeric(system(paste('grep "Mean number of unique UMIs per position"', file_loc, "| cut -d: -f4"), intern=T))
	maxUmiPerPos <- as.numeric(system(paste('grep "Max. number of unique UMIs per position"', file_loc, "| cut -d: -f4"), intern=T))
	uniqInputReads <- as.numeric(system(paste('grep "unique_input_reads"', file_loc, "| cut -d ' ' -f2"), intern=T))
	uniqOutputReads <- as.numeric(system(paste('grep "unique_output_reads"', file_loc, "| cut -d ' ' -f2"), intern=T))
	pcr_dup = round(((1 - (uniqOutputReads / uniqInputReads)) * 100),digits=2)
	df = data.frame(Results=c(
			inputAlignments,
			outputAlignments,
			umisObserved,
			meanUmiPerPos,
			maxUmiPerPos,
			uniqInputReads,
			uniqOutputReads,
			pcr_dup), check.names=F)
	colnames(df) = unlist(strsplit(f,"\\."))[3]

	dedup_frame = cbind(dedup_frame,df)
}
}

write("Dataframes complete generating CSV",stderr())
sink("batch_summary.csv")
cat("\n")
cat("Alignment Stats\n")
write.csv(align_frame)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Alignment Report\n")
write.csv(report_frame)
cat("___________________________________")
cat("\n")
cat("\n")
if(dedupDirExists){
cat("Dedupilcation Report\n")
write.csv(dedup_frame)
cat("___________________________________")
cat("\n")
cat("\n")
}
sink()

temp_dir = "./tmp"
write("CSV complete generating HTML", stderr())
#generate HTML version to rmd
html_batch = paste(temp_dir, "batch_html.R", sep="/")
htmlRmd =paste(temp_dir, "batch_html.Rmd", sep="/")
spin(html_batch, knit=F)
rmarkdown::render(htmlRmd)
write("HTML complete generating pdf")
#gnerate pdf version 
pdf_batch = paste(temp_dir, "batch_pdf.R", sep="/")
pdfRmd =paste(temp_dir, "batch_pdf.Rmd", sep="/")
spin(pdf_batch, knit=F)
rmarkdown::render(pdfRmd)
