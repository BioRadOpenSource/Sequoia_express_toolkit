library(knitr)
library(rmarkdown)
options(warn = -1)
options(tinytex.verbose = TRUE)
comArgs <- commandArgs(TRUE)

base_dir <- comArgs[1]
temp_dir <- comArgs[2]
anno_dir <- comArgs[3]

system("mv ./out ./cache")
system("mkdir ./out")
system("mkdir -p ./out/{star,fastqc,counts,cutAdapt}")
names = system("ls ./cache/counts/",intern=T ) 
names = unlist(lapply(names, function(x) unlist(strsplit(x, "\\."))[3]))

write(paste("anno_dir: ", anno_dir), stderr())

setwd(temp_dir)

for(n in names){
	print(n)
	if(dir.exists("./cache/debarcode")){
		system("mkdir ../out/debarcode")
		system(paste0("cp ./cache/debarcode/debarcode_sats.txt.",n," ./out/star/debarcode_stats.txt"))
	}
	if(dir.exists("./cahce/umiTools")){
		system("mkdir ../out/umiTools")
		system(paste0("cp ./cache/umiTools/dedup.log.",n," ./out/umiTools/dedup.log"))
	}
	#move files to proper place in out folder
	system(paste0("cp ../cache/fastqc/",n,"* ../out/fastqc/"))
	system(paste0("cp ../cache/counts/*.",n," ../out/counts/gene_counts_longRNA.summary"))
	system(paste0("cp ../cache/star/Log.final.out.",n," ../out/star/Log.final.out"))
	system(paste0("cp ../cache/star/rna_metrics.txt.",n," ../out/star/rna_metrics.txt"))
	system(paste0("cp ../cache/cutAdapt/trimlog.log.",n," ../out/cutAdapt/trimlog.log"))

	htmlReport <- "htmlReport.R"
	pdfReport <- "pdfReport.R"
	htmlReportRmd <- "htmlReport.Rmd"
	pdfReportRmd <- "pdfReport.Rmd"
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
	unlink(htmlReportRmd)
	unlink(htmlReport)
	unlink(pdfReportRmd)
	unlink(pdfReport)

	#change name of pdf and html reports so nothing is over written
	system(paste0("mv pdfReport.pdf ",n,"_pdfReport.pdf"))
	system(paste0("mv htmlReport.html ",n,"_htmlReport.html"))
	#clean fastqc folder
	system("rm ./out/fastqc/*")
}
unlink(htmlReportRmd)
unlink(htmlReport)
unlink(pdfReportRmd)
unlink(pdfReport)

