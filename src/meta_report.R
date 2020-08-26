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
metrics = lapply(metrics, function(x) data.frame(names(x), as.character(x[1,])))
sample_names = gsub("rna_metrics.txt.","",files)
print(metrics)
Metric = metrics[[1]][,1]
#rename cols of each dataframe
align_frame = data.frame(Metric)
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


meta_report = rbind(align_frame, report_frame)

write.table(meta_report, "batch_summary.csv", sep=",", row.names=F)

temp_dir = "./tmp"

#generate HTML version to rmd
html_batch = paste(temp_dir, "batch_html.R", sep="/")
htmlRmd =paste(temp_dir, "batch_html.Rmd", sep="/")
spin(html_batch, knit=F)
rmarkdown::render(htmlRmd)

#gnerate pdf version 
pdf_batch = paste(temp_dir, "batch_pdf.R", sep="/")
pdfRmd =paste(temp_dir, "batch_pdf.Rmd", sep="/")
spin(pdf_batch, knit=F)
rmarkdown::render(pdfRmd)
