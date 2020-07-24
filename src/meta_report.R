#get dir from command line arg

args = commandArgs(trailingOnly=T)

outPath = args[1]

#go recursive to get all star outfiles 

alignmentDir <- unique(dirname(list.files(outPath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="rna_metrics.txt")))
alignmentDirExists <- length(alignmentDir) >= 1
write(paste("alignmentDirExists: ", alignmentDirExists), stderr())

starDir <- unique(dirname(list.files(outPath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="Log.final.out")))
starDirExists <- length(starDir) >= 1
write(paste("logDirExists: ", starDirExists), stderr())

metrics = lapply(alignmentDir, function(x) read.table(paste0(x,"/rna_metrics.txt"), header=T, nrows=1,fill=T))
metrics = lapply(metrics, function(x) data.frame(names(x), as.character(x[1,])))
sample_names = unlist(lapply(alignmentDir, function(x) strsplit(x, "/")[[1]][3]))
Metric = metrics[[1]][,1]
#rename cols of each dataframe
align_frame = data.frame(Metric)
for(i in 1:length(metrics)){
	align_frame = cbind(align_frame, metrics[[i]][,2])
}
colnames(align_frame) = c("Metric", sample_names)

reports = lapply(starDir, function(x) read.table(paste0(x,"/Log.final.out"), header=T, sep="|",fill=T))
labels = reports[[1]][,1]
reports = lapply(reports, function(x) x[,2] <-gsub("\t","", x[,2]))
report_frame = data.frame(labels)
for(i in 1:length(reports)){
	report_frame = cbind(report_frame, reports[[i]])
}
colnames(report_frame) = c("Metric", unlist(lapply(starDir, function(x) strsplit(x, "/")[[1]][3])))


meta_report = rbind(align_frame, report_frame)

write.table(meta_report, "batch_summary.csv", sep=",", row.names=F)
