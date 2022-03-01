library(rmarkdown)
library(knitr)
library(kableExtra)
library(dplyr)
library(data.table)
library(ggplot2)
library(tibble)
library(plotly)
library(fastqcr)
library(rlist)
options(warn = -1)
options(tinytex.verbose = TRUE)
comArgs <- commandArgs(TRUE)

base_dir <- comArgs[1]
temp_dir <- comArgs[2]
anno_dir <- comArgs[3]

names = system("ls ./out/counts/",intern=T ) 
names = unlist(lapply(names, function(x) unlist(strsplit(x, "\\."))[3]))
names = names[!is.na(names)]

write(paste("anno_dir: ", anno_dir), stderr())

#fastqc
fastqcDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="_fastqc.html")))
fastqcDirExists <- length(fastqcDir) == 1
write(paste("fastqcDirExists: ", fastqcDirExists), stderr())

#debarcoding
debarcodeDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="_barcode_stats.tsv")))
debarcodeDirExists <- length(debarcodeDir) == 1
write(paste("debarcodeDirExists: ", debarcodeDirExists), stderr())

#trimming
trimDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="trimlog.log.*")))
trimDirExists <- length(trimDir) == 1
write(paste("trimDirExists: ", trimDirExists), stderr())

#alignments
alignmentDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="rna_metrics.txt.*")))
alignmentDirExists <- length(alignmentDir) == 1
write(paste("alignmentDirExists: ", alignmentDirExists), stderr())

#deduplication
dedupDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="dedup.log.*")))
dedupDirExists <- length(dedupDir) == 1
write(paste("dedupDirExists: ", dedupDirExists), stderr())

#counts
countsDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="gene_counts_longRNA.summary.*")))
countsDirExists <- length(countsDir) == 1
write(paste("countsDirExists: ", countsDirExists), stderr())

parseQcFile <- function(fileName)
{
	  qc <- qc_read(fileName, modules="Per base sequence quality")
  return(qc$per_base_sequence_quality)
}


setwd(temp_dir)
for(n in names){
	#set vraibles to be used by all to lower processing time 

	if(debarcodeDirExists){
		deb <- read.table(paste0(debarcodeDir,"/",n,"_R1_barcode_stats.tsv"), fill=T,sep="\t",header=T)
		inputReads <- as.numeric(as.character(deb$count[1]))
		validBcReads <- as.numeric(as.character(deb$count[2]))
		invalidBcReads <- inputReads-validBcReads

		#create data frame
		deb_df <- data.frame(
		   Metric = c("Input Reads", "Reads with Valid UMI", "% Reads with Valid UMI"),
		   Value = c(inputReads, validBcReads, signif(validBcReads/inputReads, 4) * 100),    
		   stringsAsFactors = FALSE
		   )

	}

	if(trimDirExists){
		rt <- read.table(paste0(trimDir, "/trimlog.log.",n), skip=7, fill=T, sep=":", stringsAsFactors=F)
		names(rt) <- c("Metric","Value")
		rt$Value <- as.numeric(gsub(",","",unlist(lapply(strsplit(as.character(rt$Value), split="\\s+"), `[[`, 2)))) #this is gross, i'm sorry for nesting 6 functions
		if(length(grep("Read 2", rt$Metric))>0){
			index = grep("Read [1|2]$", rt$Metric)
			rt$Metric[index] = paste0(rt$Metric[index]," basepairs")
		}	
	}

	if(alignmentDirExists){
		#Using a subset of PICARD outputs
		rna <- read.table(paste0(alignmentDir, "/rna_metrics.txt.",n), nrows=1, header=T, fill=T)
		starReport <- read.table(paste0(alignmentDir, "/Log.final.out.",n), sep="|", fill=T)
		starReport$V2 <- gsub("\t","",starReport$V2)

		#parse out relevant stuff from STAR report
		inputReads <- as.numeric(starReport$V2[grep("Number of input reads",starReport$V1)])
		uniqueMapped <- as.numeric(starReport$V2[grep("Uniquely mapped reads number",starReport$V1)])
		multiMapped <- as.numeric(starReport$V2[grep("Number of reads mapped to multiple loci",starReport$V1)])
		tooManyMapped <- as.numeric(starReport$V2[grep("Number of reads mapped to too many loci",starReport$V1)])
		unmapped <- inputReads-uniqueMapped-multiMapped

		align_df <- data.frame("Reads Input" = inputReads,
				 "Uniquely Mapped Reads" = uniqueMapped,
				 "Multi-mapped Reads" = multiMapped,
      				 "Reads mapped to too many loci" = tooManyMapped,
     				 "Unmapped Reads" = unmapped,
				 "PF Bases" = rna$PF_BASES,
				 "PF Aligned Bases" = rna$PF_ALIGNED_BASES,
				 "Coding Bases" = rna$CODING_BASES,
				 "UTR Bases" = rna$UTR_BASES,
				 "Intronic Bases" = rna$INTRONIC_BASES,
				 "Intergenic Bases" = rna$INTERGENIC_BASES,
				 "Ribosomal Bases" = rna$RIBOSOMAL_BASES,
				 "Median CV Coverage" = rna$MEDIAN_CV_COVERAGE,
				 "Median 5' Bias" = rna$MEDIAN_5PRIME_BIAS,
				 "Median 3' Bias" = rna$MEDIAN_3PRIME_BIAS,
				 "Median 5' to 3' Bias" = rna$MEDIAN_5PRIME_TO_3PRIME_BIAS,
				 "% Stranded" = rna$PCT_CORRECT_STRAND_READS * 100,
				 "% rRNA bases" = rna$PCT_RIBOSOMAL_BASES *100,
				 check.names= F)
		align_df <- as.data.frame(t(align_df)) %>% rownames_to_column()
		names(align_df) <- c("Metric", "Value")
	}
	
	if(dedupDirExists){
		file_loc = paste0(dedupDir, "/dedup.log.",n)
		umisObserved <- as.numeric(system(paste('grep -F "unique_umi"', file_loc, "| cut -d' ' -f2"), intern=T))
		inputAlignments <- as.numeric(system(paste('grep "Reads In"', file_loc, "| cut -d' ' -f3"), intern=T))
		outputAlignments <- as.numeric(system(paste('grep "Reads Out"', file_loc, "| cut -d' ' -f3"), intern=T))
		uMate <- as.numeric(system(paste('grep "No Mate"', file_loc, "| cut -d' ' -f3"), intern=T))
		#unpaired <- as.numeric(system(paste('grep "Reads Unpaired"', file_loc, "| cut -d' ' -f3"), intern=T))
		#meanUmiPerPos <- as.numeric(system(paste('grep "Mean number of unique UMIs per position"', file_loc, "| cut -d: -f4"), intern=T))
		#maxUmiPerPos <- as.numeric(system(paste('grep "Max. number of unique UMIs per position"', file_loc, "| cut -d: -f4"), intern=T))
		#chimera <- as.numeric(system(paste('grep "Reads Chimeric"', file_loc, "| cut -d ' ' -f3"), intern=T))
		uniqInputReads <- as.numeric(system(paste('grep "unique_input_reads"', file_loc, "| cut -d ' ' -f2"), intern=T))
		uniqOutputReads <- as.numeric(system(paste('grep "unique_output_reads"', file_loc, "| cut -d ' ' -f2"), intern=T))

		dedup_df <- data.frame("Total input alignments" = inputAlignments,
				 "Total output alignments" = outputAlignments,
		                 "Unique UMIs observed" = umisObserved,
				 #"Average UMIs per position" = meanUmiPerPos,
				 #"Maximum UMIs per position" = maxUmiPerPos,
				 #"Chimeric Reads" =chimera,
				 "Reads with unpaired mate" = uMate,
				 "Unique Input Reads" = uniqInputReads,
				 "Unique Output Reads" = uniqOutputReads,
				 "% PCR Duplicates" = (1 - (uniqOutputReads / uniqInputReads)) * 100,
				 check.names= F)
		dedup_df <- as.data.frame(t(dedup_df)) %>% rownames_to_column()
		names(dedup_df) <- c("Metric", "Value")
	}

	if(countsDirExists){
		write("Processing gene_counts_longRNA.summary", stderr())
		longRNAcounts <- read.table(paste0(countsDir, "/gene_counts_longRNA.summary.",n), skip=1)
		colnames(longRNAcounts) <- c("Result", "Count")
		longRNAcounts$Result <- gsub("_", " ", longRNAcounts$Result)
		longRNAcounts <- rbind(data.frame(Result="Total Alignments", Count=sum(longRNAcounts$Count)), longRNAcounts)
		
		#handle biotypes
		write("Processing gene_counts_longRNA", stderr())
		countLong <- read.table(paste0(countsDir, "/gene_counts_longRNA.",n), header=T, sep="\t", col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
		write("Processing gene_biotypes.tsv", stderr())
		biotypes <- read.table(paste(anno_dir,"gene_biotypes.tsv", sep="/"), sep="\t", header=T)
		#if("rRNA" %in% biotypes[,"gene_biotype"]){
		#biotypes[biotypes["gene_biotype"] == "rRNA",]["gene_biotype"] <- "mitochondrial_rRNA"
		#}
		countLong <- left_join(countLong, biotypes, by = c("Gene" = "gene_id"))
		countAll <- countLong

		countLong = countLong[,c("Gene","Length","Count","gene_biotype")]
		colnames(countLong) = c("Gene","Length","Count","Biotype")

		longRNAcounts <- rbind(longRNAcounts,data.frame(Result="Genes with >0 Counts", Count=dim(countLong[countLong$Count >0,])[1]))
		countByBiotype <- countAll %>% filter(!is.na(gene_biotype)) %>% group_by(gene_biotype) %>% summarise(count = sum(Count)) %>% arrange(-count)
		countByBiotype$gene_biotype <- factor(countByBiotype$gene_biotype, levels = unique(countByBiotype$gene_biotype)[order(countByBiotype$count, decreasing = TRUE)])
		countByBiotype <- countByBiotype %>% filter(count > 0) #filter biotypes with no counts

	}

	env <- Sys.getenv(c("Software Name","FASTQC_VERSION","STAR_VERSION","PICARD_VERSION","RUMI_VERSION","SUBREAD_VERSION"))
	env <- as.data.frame(env, stringsAsFactors=FALSE) %>% tibble::rownames_to_column()
	umi_tools_version <- system("rumi -V", intern=T)
	umi_tools_version <- strsplit(umi_tools_version, " ")[[1]][2]
	env[which(env$rowname=="RUMI_VERSION"), 2] = gsub(" ", "", umi_tools_version)
	write("Preparing to read imageInfo.txt", stderr())
	containerInfo <- read.table("/opt/biorad/imageInfo.txt", stringsAsFactors=FALSE,sep=":")
	write("Read imageInfo.txt", stderr())
	containerInfo <- data.frame("rowname" = containerInfo[,1], env = containerInfo[,2], stringsAsFactors=FALSE)
	containerInfo[3,2] <- substr(containerInfo[3,2], 1,8)
	anno_path <- unlist(strsplit(anno_dir,"/"))
	referenceGenome <- anno_path[grepl("hg38|mm10|rnor6|tair10|sacCer3|dm6|danRer11|ce11", anno_path)]
	if(length(referenceGenome) == 0)
	{
		        referenceGenome = "NA"
	}
	isErcc <- any(grepl("ercc", anno_path))
	write(paste("Preparing to read: ", paste(anno_dir,"annotation_version.txt", sep="/")), stderr())
	anno_version <- read.table(paste(anno_dir,"annotation_version.txt", sep="/"), comment.char="", fill=T, sep=",")
	write("Read annotation_version", stderr())
	anno_source <- gsub("#!annotation-source ", "", anno_version$V1[grep("annotation-source", anno_version$V1)])
	localVars <- data.frame(rowname = c("Sample Name","Reference Genome", "Annotation Source", "UMI Aware", "ERCC"), env = c(n,referenceGenome, anno_source, dedupDirExists, isErcc), stringsAsFactors=FALSE)
	env <- rbind(containerInfo, localVars, env)
	env[nrow(env) + 1,] = list("Report Generated", paste(as.character(Sys.time()), "UTC"))
	colnames(env) <- NULL

	#placeholder for biotype plot
	biotype_pl = NULL
	
	debar <- data.frame(
		"Input Reads"	= 		"The number of reads in the input FASTQ files", 
		"Reads with Valid UMI" = 	"The number of R1 reads with a valid UMI" ,
		"% Reads with Valid UMI" =	"The percentage of input reads with a valid UMI",
		check.names=F)
		glossary <- data.frame(
		"Reads Input" = 		"The number of reads input to alignment.",
	      	"Uniquely Mapped Reads" = 	"The number of reads mapped to a unique genetic locus.",
	      	"Multi-mapped Reads" = 		"The number of reads mapping to multiple genetic loci.",
	      	"Reads mapped to too many loci"="The number of reads mapping to too many genetic loci for further processing (a subset of the value listed above).",
	      	"Unmapped Reads" = 		"The number of reads that did not map to the reference genome.",
	      	"PF Bases" = 			"The total number of PF bases including non-aligned reads.",
	      	"PF Aligned Bases" = 		"The total number of aligned PF bases. Non-primary alignments are not counted. Bases in aligned reads that do not correspond to reference (e.g. soft clips, insertions) are not counted.",
	     	"Coding Bases" = 		"Number of bases in primary alignments that align to a non-UTR coding base for some gene, and not ribosomal sequence.",
	      	"UTR Bases" = 			"Number of bases in primary alignments that align to a UTR base for some gene, and not a coding base.",
	      	"Intronic Bases" = 		"Number of bases in primary alignments that align to an intronic base for some gene, and not a coding or UTR base.",
	      	"Intergenic Bases" =		"Number of bases in primary alignments that do not align to any gene." ,
	      	"Ribosomal Bases" = 		"Number of bases in primary alignments that align to ribosomal sequence.",
	      	"Median CV Coverage" = 		"The median coefficient of variation (CV) or stdev/mean for coverage values of the 1000 most highly expressed transcripts. Ideal value = 0.",
	      	"Median 5' Bias" = 		"The median 5 prime bias of the 1000 most highly expressed transcripts. The 5 prime bias is calculated per transcript as: mean coverage of the 5 prime-most 100 bases divided by the mean coverage of the whole transcript.",
	      	"Median 3' Bias" = 		"The median 3 prime bias of the 1000 most highly expressed transcripts, where 3 prime bias is calculated per transcript as: mean coverage of the 3 prime-most 100 bases divided by the mean coverage of the whole transcript.",
	      	"Median 5' to 3' Bias" = 	"The ratio of coverage at the 5 prime end to the 3 prime end based on the 1000 most highly expressed transcripts.",
	      	"% Stranded" = 			"The percentage of reads corresponding to transcripts which map to the correct strand of a reference genome ",
	      	"% rRNA bases" = 		"Percent of aligned bases that mapped to regions encoding ribosomal RNA",
		check.names=F)
	dedup_terms = data.frame(	
		
		"Total input alignments" = "The total number of alignments passed into deduplication",
		"Total output alignments" = "The total number of alignments output after deduplication",
		"Unique UMIs observed" = 	"The total number of unique UMIs observed",
		"Reads with unpaired mate" = 	"Number of reads from input that have no paired mate post alignment",
		"Unique Input Reads" = 		"The number of unique input reads passed into deduplication",
		"Unique Output Reads" = 	"The number of unique output reads after deduplication",
		"% PCR Duplicates" = 		"The percentage of reads that are PCR duplicates ((1 - (Unique Output Reads / Unique Input Reads)) * 100)",
		check.names= F
	)
	glossary = as.data.frame(t(glossary))
	colnames(glossary) <-NULL
	dedup_terms = as.data.frame(t(dedup_terms))
	colnames(dedup_terms) <-NULL
	debar = as.data.frame(t(debar))
	colnames(debar) <-NULL

	#Breaker ------
	htmlReport <- paste(temp_dir, "htmlReport.R", sep="/")
	pdfReport <- paste(temp_dir, "pdfReport.R", sep="/")
	htmlReportRmd <- paste(temp_dir, "htmlReport.Rmd", sep="/")
	pdfReportRmd <- paste(temp_dir, "pdfReport.Rmd", sep="/")
	write("Spinning up the html report", stderr())
	spin(htmlReport, knit=FALSE)
	write("html report created", stderr())
	write("Spinning up the pdfreport", stderr())
	spin(pdfReport, knit=FALSE)
	write("pdf report created", stderr())
	write("Rendering html report", stderr())
	rmarkdown::render(htmlReportRmd)
	write("html report rendered")
	write("Rendering pdf report", stderr())
	rmarkdown::render(pdfReportRmd)
	write("pdf report rendered")
	write("generating CSV report")
	source("./csvReport.R")
	write("CSV report complete")
	#change name of pdf and html reports so nothing is over written
	system(paste0("mv pdfReport.pdf ",n,"_pdfReport.pdf"))
	system(paste0("mv htmlReport.html ",n,"_htmlReport.html"))
	system(paste0("mv csvReport.csv ",n,"_csvReport.csv"))

}
unlink(htmlReportRmd)
unlink(htmlReport)
unlink(pdfReportRmd)
unlink(pdfReport)

