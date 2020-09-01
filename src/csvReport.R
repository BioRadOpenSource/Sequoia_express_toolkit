library(fastqcr)


qcFiles <- list.files(fastqcDir, full.names=TRUE,recursive=T)

r1 <- qcFiles[grepl(paste0(n,"_R1_.*zip"),qcFiles)]
r2 <- qcFiles[grepl(paste0(n,"_R2_.*zip"),qcFiles)]

r1Qc <- sapply(r1, parseQcFile, simplify=F)
r1df = as.data.frame(r1Qc)
colnames(r1df) = paste0("R1_",gsub(".*zip.","",colnames(r1df)))
if(length(r2) >0){
	r2Qc <- sapply(r2, parseQcFile, simplify=F)
	r2df = as.data.frame(r2Qc)
	colnames(r2df) = paste0("R2_",gsub(".*zip.","",colnames(r2df)))
	if(dim(r1df)[1] >= dim(r2df)[1]){
		diff = dim(r1df)[1] - dim(r2df)[1]
		empty = data.frame(matrix(NA,nrow=diff,ncol= 7))
		colnames(empty) = colnames(r2df) 
		r2df = rbind(r2df, empty)
		r1df = cbind(r1df,r2df)
	}
}

####r1df is not either one of both the fastq files 


####trimming 
trim <- data.frame(" "="Reads",
	"Reads Input" = rt$Value[grep("Total read.* processed",rt$Metric)],
	"Reads Too Short" = rt$Value[grep(".* that were too short", rt$Metric)],
	"Reads Written" = rt$Value[grep(".* written \\(passing filters\\)",rt$Metric)], check.names=FALSE)

###alignment
#Show threshold genes (will be gene count folder) 

#show threshold with read count filter
colnames(env) <- c("Parameter","Value")
#combine each to one csv
sink("csvReport.csv")
#load.image("/opt/biorad/src/vendor-logo.png")
cat("\n")
cat("Fastqc Report\n")
write.csv(r1df)
cat("___________________________________")
cat("\n")
cat("\n")
if(debarcodeDirExists){
cat("Debarcode Report\n")
write.csv(deb_df)
cat("___________________________________")
cat("\n")
cat("\n")
}
cat("Trimming Report\n")
write.csv(rt)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Alignment Stats\n")
write.csv(align_df)
cat("___________________________________")
cat("\n")
cat("\n")
if(dedupDirExists){
cat("Dedupilcation Report\n")
write.csv(dedup_df)
cat("___________________________________")
cat("\n")
cat("\n")
}
cat("Gene Count Summary\n")
write.csv(longRNAcounts)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Count table with biotypes\n")
write.csv(countByBiotype)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Enviorment Metadata\n")
write.csv(env,row.names=F)
cat("___________________________________")
cat("\n")
cat("\n")
sink()
