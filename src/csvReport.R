library(fastqcr)

# get each report 
comArgs <- commandArgs(TRUE)

base_dir <- comArgs[1]
temp_dir <- comArgs[2]
anno_dir <- comArgs[3]

#fastqc
fastqcDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="_fastqc.html")))
fastqcDirExists <- length(fastqcDir) == 1
write(paste("fastqcDirExists: ", fastqcDirExists), stderr())

#debarcoding
debarcodeDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="debarcode_stats.txt")))
debarcodeDirExists <- length(debarcodeDir) == 1
write(paste("debarcodeDirExists: ", debarcodeDirExists), stderr())

#trimming
trimDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="trimlog.log")))
trimDirExists <- length(trimDir) == 1
write(paste("trimDirExists: ", trimDirExists), stderr())

#alignments
alignmentDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="rna_metrics.txt")))
alignmentDirExists <- length(alignmentDir) == 1
write(paste("alignmentDirExists: ", alignmentDirExists), stderr())

starDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="Log.final.out")))
starDirExists <- length(starDir) == 1
write(paste("alignmentDirExists: ", starDirExists), stderr())

#deduplication
dedupDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="dedup.log")))
dedupDirExists <- length(dedupDir) == 1
write(paste("dedupDirExists: ", dedupDirExists), stderr())

#count summary
countsDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="longRNA.summary")))
countsDirExists <- length(countsDir) == 1
write(paste("countsDirExists: ", countsDirExists), stderr())
#actual gene counts
gene_counts = read.table(paste0(countsDir,"/gene_counts_longRNA"), sep='\t', comment='#', header=T)
colnames(gene_counts)[7] ="Counts"
gene_counts = gene_counts[,c("Geneid","Counts")]

#create combined table 
#### fastqc files #####
parseQcFile <- function(fileName)
{
	  qc <- qc_read(fileName, modules="Per base sequence quality")
  return(qc$per_base_sequence_quality)
}

qcFiles <- list.files(fastqcDir, full.names=TRUE,recursive=T)

r1 <- qcFiles[grepl("_R1_.*zip",qcFiles)]
r2 <- qcFiles[grepl("_R2_.*zip",qcFiles)]

r1Qc <- sapply(r1, parseQcFile, simplify=F)
r1df = as.data.frame(r1Qc)
colnames(r1df) = paste0("R1_",gsub(".*zip.","",colnames(r1df)))
if(length(r2) >0){
	r2Qc <- sapply(r2, parseQcFile, simplify=F)
	r2df = as.data.frame(r2Qc)
	colnames(r2df) = paste0("R2_",gsub(".*zip.","",colnames(r2df)))
	if(dim(r1df)[1] > dim(r2df)[1]){
		diff = dim(r1df)[1] - dim(r2df)[1]
		empty = data.frame(matrix(NA,nrow=diff,ncol= 7))
		colnames(empty) = colnames(r2df) 
		r2df = rbind(r2df, empty)
		r1df = cbind(r1df,r2df)
	}
}

####r1df is not either one of both the fastq files 

###Debarcode 
if(debarcodeDirExists){
deb <- read.table(paste(debarcodeDir,"debarcode_stats.txt", sep="/"), fill=T)
inputReads <- as.numeric(as.character(deb$V3[1]))
validBcReads <- as.numeric(as.character(deb$V3[2]))
invalidBcReads <- inputReads-validBcReads

dBarcode <- data.frame(
	Metric = c("Input Reads", "Reads with Valid UMI", "% Reads with Valid UMI"),
	Value = prettyNum(c(inputReads, validBcReads, signif(validBcReads/inputReads, 3) * 100), big.mark = ",", scientific = F),
	stringsAsFactors = FALSE)
}
####trimming 
trim <- read.table(paste(trimDir, "trimlog.log", sep="/"), skip=7, nrows=7, fill=T, sep=":") 
names(trim) <- c("Metric","Value")
trim$Value <- as.numeric(gsub(",","",unlist(lapply(strsplit(as.character(trim$Value), split="\\s+"), `[[`, 2)))) #this is gross, i'm sorry for nesting 6 functions

###alignment

#Using a subset of PICARD outputs
picard <- read.table(paste(alignmentDir, "rna_metrics.txt", sep="/"), nrows=1, header=T, fill=T)
picard_df = data.frame(names(picard),as.character(picard[1,]))
names(picard_df) = c("Metric", "Value")		       
#alignement report
starReport <- read.table(paste(starDir, "Log.final.out", sep="/"), sep="|", fill=T)
starReport$V2 <- gsub("\t","",starReport$V2)
names(starReport) = c("Metric","Value")

###gene counts
longRNAcounts <- read.table(paste(countsDir, "gene_counts_longRNA.summary", sep="/"), skip=1)
colnames(longRNAcounts) <- c("Result", "Count")
longRNAcounts$Result <- gsub("_", " ", longRNAcounts$Result)
longRNAcounts <- rbind(data.frame(Result="Total Alignments", Count=sum(longRNAcounts$Count)), longRNAcounts)
countLong = countLong <- read.table(paste(countsDir, "gene_counts_longRNA", sep="/"), sep="\t", header=T, col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
longRNAcounts <- rbind(longRNAcounts,data.frame(Result="Genes with >0 reads", Count=dim(countLong[countLong$Count >1,])[1]))

### biotypes w/ threshold 
biotypes <- read.table(paste(anno_dir,"gene_biotypes.tsv", sep="/"), sep="\t", header=T)
countLong = countLong <- read.table(paste(countsDir, "gene_counts_longRNA", sep="/"), sep="\t", header=T, col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
#countLong <- left_join(countLong, biotypes, by = c("Gene" = "gene_id"))
#counts with biotypes as one table 
colnames(biotypes)[1] = "Gene"
countsWbiotype = merge(countLong, biotypes[biotypes[,"Gene"] %in% countLong[,"Gene"],], by = c("Gene"), all=T) 


#enviormental metadata from other reports
env <- Sys.getenv(c("FASTQC_VERSION","STAR_VERSION","BEDTOOLS_VERSION","PICARD_VERSION","UMI_TOOLS_VERSION","SUBREAD_VERSION","SAMBAMBA_VERSION"))
env <- as.data.frame(env, stringsAsFactors=FALSE)
umi_tools_version <- system("umi_tools --version", intern=T)
umi_tools_version <- strsplit(umi_tools_version, ":")[[1]][2]
env[which(env$rowname=="UMI_TOOLS_VERSION"), 2] = gsub(" ", "", umi_tools_version)
containerInfo <- read.table("/opt/biorad/imageInfo.txt", stringsAsFactors=FALSE)
containerInfo <- data.frame("rowname" = paste(containerInfo[,1], containerInfo[,2]), env = containerInfo[,3], stringsAsFactors=FALSE)
containerInfo[3,2] <- substr(containerInfo[3,2], 1,7)
anno_path <- unlist(strsplit(anno_dir,"/"))
referenceGenome <- anno_path[grepl("hg38|mm10|rnor6", anno_path)]
if(length(referenceGenome) == 0)
{
	  referenceGenome = "NA"
}
isErcc <- any(grepl("ercc", anno_path))
anno_version <- read.table(paste(anno_dir,"annotation_version.txt", sep="/"), comment.char="", fill=T, sep=",")
anno_source <- gsub("#!annotation-source ", "", anno_version$V1[grep("annotation-source", anno_version$V1)])
localVars <- data.frame(rowname = c("Reference Genome", "Annotation Source", "UMI Aware", "ERCC"), env = c(referenceGenome, anno_source, dedupDirExists, isErcc), stringsAsFactors=FALSE)
env <- rbind(containerInfo, localVars, env)
env[nrow(env) + 1,] = list("Report Generated", paste(as.character(Sys.time()), "UTC"))
colnames(env) <- NULL

#Show threshold genes (will be gene count folder) 


#show threshold with read count filter


#combine each to one csv
sink("Sequoia_express_report.csv")
load.image("/opt/biorad/src/vendor-logo.png")
cat("\n")
cat("Fastqc Report\n")
write.csv(r1df)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Debarcode Report\n")
write.csv(debarcode)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Trimming Report\n")
write.csv(trim)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Alignment Stats\n")
write.csv(picard)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Alignement Report\n")
write.csv(starReport)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Gene Count Summary\n")
write.csv(longRNAcounts)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Count table with biotypes\n")
write.csv(countsWbiotype)
cat("___________________________________")
cat("\n")
cat("\n")
cat("Enviorment Metadata\n")
write.csv(env)
cat("___________________________________")
cat("\n")
cat("\n")
sink()
