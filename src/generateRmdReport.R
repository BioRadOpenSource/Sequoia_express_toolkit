library(knitr)
library(rmarkdown)
options(warn = -1)
options(tinytex.verbose = TRUE)
comArgs <- commandArgs(TRUE)

base_dir <- comArgs[1]
temp_dir <- comArgs[2]
anno_dir <- comArgs[3]

names = system("ls ./out/counts/",intern=T ) 
names = unlist(lapply(names, function(x) unlist(strsplit(x, "\\."))[3]))

write(paste("anno_dir: ", anno_dir), stderr())

setwd(temp_dir)

for(n in names){
	print(n)
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
	unlink(htmlReportRmd)
	unlink(htmlReport)
	unlink(pdfReportRmd)
	unlink(pdfReport)

	#change name of pdf and html reports so nothing is over written
	system(paste0("mv pdfReport.pdf ",n,"_pdfReport.pdf"))
	system(paste0("mv htmlReport.html ",n,"_htmlReport.html"))
}
unlink(htmlReportRmd)
unlink(htmlReport)
unlink(pdfReportRmd)
unlink(pdfReport)

