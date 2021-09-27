#' ---
#' title: |
#'   | \vspace{2cm} \LARGE{SEQuoia Express Analysis Report}
#' header-includes:
#' - \usepackage{titling}
#' - \usepackage{fontspec}
#' - \setmainfont{FreeSans}
#' - \usepackage{booktabs}
#' - \usepackage[table]{xcolor}
#' - \pretitle{\vspace{5cm}\begin{center}\LARGE\includegraphics[width=12cm]{/opt/biorad/src/vendor-logo.png}\\[\bigskipamount]}
#' - \posttitle{\end{center}\newpage}
#' output: 
#'   pdf_document:
#'     latex_engine: xelatex
#'     toc: true
#' always_allow_html: yes
#' ---


#' \centering
#' \raggedright
#' \newpage

#+ setup, include=FALSE

#basics
library(knitr)
library(kableExtra)
library(dplyr) 
library(data.table)
library(ggplot2)
library(tibble)
library(plotly)
library(fastqcr)
library(rlist)

#muting warnings
options(warn=-1)

#setwd("/")
#dictate kable styling options
kableStyle <- c("striped", "condensed", "hover", "responsive")


#' `r if(fastqcDirExists) { "# Read QC" }`
#+ eval=fastqcDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F

qcFiles <- list.files(fastqcDir, full.names=TRUE)[grepl("zip",list.files(fastqcDir))]

r1 <- qcFiles[grepl(paste0(n,"_R1_"),qcFiles)]
r2 <- qcFiles[grepl(paste0(n,"_R2_"),qcFiles)]

r1Qc <- sapply(r1, parseQcFile, simplify=F)

r1QcDf <- list.cbind(r1Qc)

p1 <- plot_ly(width = 700) %>% 
  layout(
    yaxis = list(title = "Q Score", range = c(0, 40)),
    xaxis = list(title = "Position in read (bp)",
                 tickvals = seq(1, dim(r1QcDf)[1], by=2),
                 tickmode = "array",
                 ticktext = r1Qc[[1]]$Base[seq(1, dim(r1QcDf)[1], by=2)],
                 tickangle = 90,
                 ticks = "outside"))

for(i in grep("Mean",names(r1QcDf)))
{
  p1 <- add_trace(p1, x = ~c(1:(length(r1QcDf[,1]))), y = r1QcDf[,i], type = 'scatter', mode ='lines', name = gsub("_fastqc.zip.Mean", " ", basename(names(r1QcDf[i]))))
}

p1 <- p1 %>% layout(
  shapes = list(
    list(type = "rect",
         fillcolor = "red",
         line = list(color="red"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=0, y1=20, yref="y",
         layer = "below"),
    list(type = "rect",
         fillcolor = "yellow",
         line = list(color="yellow"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=20, y1=28, yref="y",
         layer = "below"),
    list(type = "rect",
         fillcolor = "green",
         line = list(color="green"),
         opacity = 0.3,
         x0=0, x1=length(r1QcDf[,1]), xref="x",
         y0=28, y1=40, yref="y",
         layer = "below")
    ),
    legend = list(
      orientation = 'h',
      x = 10,
      y = -0.15)
)

if (length(r2) > 0) {

  r2Qc <- sapply(r2, parseQcFile, simplify=F)
  r2QcDf <- list.cbind(r2Qc)
  p2 <- plot_ly(width = 700) %>% 
    layout(
      yaxis = list(title = "Q Score", range = c(0, 40)),
      xaxis = list(title = "Position in read (bp)",
       tickvals = seq(1, dim(r2QcDf)[1], by=2),
       tickmode = "array",
       ticktext = r1Qc[[1]]$Base[seq(1, dim(r2QcDf)[1], by=2)],
       tickangle = 90,
       ticks = "outside"))

  for(i in grep("Mean", names(r2QcDf)))
  {
    p2 <- add_trace(p2, x = ~c(1:(length(r2QcDf[,1]))), y = r2QcDf[,i], type = 'scatter', mode ='lines', name = gsub("_fastqc.zip.Mean", " ", basename(names(r2QcDf[i]))))
  }

  p2 <- p2 %>% layout(
    shapes = list(
      list(type = "rect",
     fillcolor = "red",
     line = list(color="red"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=0, y1=20, yref="y",
     layer = "below"),
      list(type = "rect",
     fillcolor = "yellow",
     line = list(color="yellow"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=20, y1=28, yref="y",
     layer = "below"),
      list(type = "rect",
     fillcolor = "green",
     line = list(color="green"),
     opacity = 0.3,
     x0=0, x1=length(r2QcDf[,1]), xref="x",
     y0=28, y1=40, yref="y",
     layer = "below")
    ),
    legend = list(
      orientation = 'h',
      x = 10,
      y = -0.15)
  )
}

#' `r if(exists("p1")) { "## Read 1" }`
#+ eval=exists("p1"), echo=FALSE, fig.asp=1, fig.align="center", message=F
p1

#' `r if(exists("p2")) { "## Read 2" }`
#+ eval=exists("p2"), echo=FALSE, fig.asp=1, fig.align="center", message=F
p2

#' \newpage

#' `r if(debarcodeDirExists) { "# UMI Parsing" }`
#+ eval=debarcodeDirExists, echo=FALSE, fig.width=8, fig.height=5, fig.align="left"
if(debarcodeDirExists){

#create kable output
kable(deb_df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#plotly
dfx <- data.frame(" "="barcode", "valid"=validBcReads, "invalid"=inputReads-validBcReads)
pl <- plot_ly(dfx, x = ~valid, y = ~" ", type = "bar", name="Valid", orientation = "h", hoverinfo = 'text', text = ~paste("Valid Reads: ", valid)) %>% 
  add_trace(x = ~invalid, name = "Invalid", hoverinfo = 'text', text = ~paste("Invalid Reads: ", invalid)) %>% 
  layout(barmode = 'stack', xaxis = list(title = "Reads"), yaxis = list(title = ""))

#create plot
pl
}
#' \newpage

#' `r if(trimDirExists) { "# Read Trimming" }`
#+ eval=trimDirExists, echo=FALSE, fig.asp=0.75, fig.align="center"

#Import metadata from cutadapt and format it

#generate table for plotting
dfx <- data.frame(" "="Reads",
	"Reads Input" = rt$Value[grep("Total read.* processed",rt$Metric)],
	"Reads Too Short" = rt$Value[grep(".* that were too short", rt$Metric)],
	"Reads Written" = rt$Value[grep(".* written \\(passing filters\\)",rt$Metric)], check.names=FALSE)

pl <- plot_ly(dfx, x = ~`Reads Too Short`, y = ~" ", type = "bar", name = "Reads Too Short", orientation = "h", hoverinfo = "text", text = ~paste("Reads Too Short: ", dfx$`Reads Too Short`)) %>%
  add_trace(x = ~`Reads Written`, name = "Reads Written", hoverinfo = "text", text = ~paste("Reads Written: ", dfx$`Reads Written`)) %>%
  layout(barmode = "stack", xaxis = list(title = "Reads"), yaxis = list(title = ""), legend = list(orientation = 'h', x=0.3, y=-0.4))

#write table
rt$Value <- prettyNum(rt$Value, big.mark = ",", scientific = F) #this is gross, i'm sorry for nesting 6 functions
kable(rt, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


#' \newpage

#' `r if(alignmentDirExists) { "# Alignment" }`
#+ eval=alignmentDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F

#Using a subset of PICARD outputs
dfx <- data.frame(" "="Bases",
                  "Coding" = align_df$Value[which(align_df$Metric=="Coding Bases")],
                  "UTR" = align_df$Value[which(align_df$Metric=="UTR Bases")],
                  "Intronic" = align_df$Value[which(align_df$Metric=="Intronic Bases")],
                  "Intergenic" = align_df$Value[which(align_df$Metric=="Intergenic Bases")],
                  "Ribosomal" = align_df$Value[which(align_df$Metric=="Ribosomal Bases")], check.names=F)


pl <- plot_ly(dfx, x = ~`Coding`, y= ~" ", type = "bar", name="Coding", orientation = "h", hoverinfo = 'text', text = ~paste("Coding: ", dfx$`Coding`)) %>% 
  add_trace(x = ~`UTR`, name = "UTR", hoverinfo = 'text', text = ~paste("UTR: ", round(as.numeric(dfx$`UTR`)))) %>%
  add_trace(x = ~`Intronic`, name = "Intronic", hoverinfo = 'text', text = ~paste("Intronic: ", round(as.numeric(dfx$`Intronic`)))) %>%
  add_trace(x = ~`Intergenic`, name = "Intergenic", hoverinfo = 'text', text = ~paste("Intergenic: ", round(as.numeric(dfx$`Intergenic`)))) %>%
  add_trace(x = ~`Ribosomal`, name = "Ribosomal", hoverinfo = 'text', text = ~paste("Ribosomal: ", round(as.numeric(dfx$`Ribosomal`)))) %>%
  layout(barmode = 'stack', xaxis = list(title = "Aligned Bases"), yaxis = list(title = ""), legend = list(orientation = 'h', x=0.3, y=-0.4))

#create table
align_df$Value <- prettyNum(align_df$Value, big.mark = ",", scientific = F)
kable(align_df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

rt_cov <- read.table(paste0(alignmentDir, "/rna_metrics.txt.",n), skip = 10, header=T, fill=T)

cov <- plot_ly(width = 700) %>% 
  layout(
    yaxis = list(title = "Normalized Coverage", range = c(0, max(rt_cov$All_Reads.normalized_coverage+0.1*rt_cov$All_Reads.normalized_coverage))),
    xaxis = list(title = "Normalized Position",
                 tickvals = seq(0, 100, by=2),
                 tickmode = "array",
                 ticktext = seq(0, 100, by=2),
                 tickangle = 90,
                 ticks = "outside"))

cov <- add_trace(cov, x = ~c(1:101), y = rt_cov$All_Reads.normalized_coverage, type = 'scatter', mode ='lines')

#' \newpage

#' `r if(exists("cov")) { "## Transcript Coverage" }`
#+ eval=exists("cov"), echo=FALSE, fig.asp=1, fig.align="center", message=F
cov

#' \newpage

#' `r if(exists("pl")) { "## Base Distribution" }`
#+ eval=exists("pl"), echo=FALSE, fig.asp=1, fig.align="center", message=F
pl

#' \newpage

#' `r if(dedupDirExists) { "# Deduplication" }`
#+ eval=dedupDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
dedup_df$Value <- prettyNum(dedup_df$Value, big.mark = ",", scientific=FALSE)
kable(dedup_df, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage

#' `r if(countsDirExists) { "# Transcriptome" }`
#+ eval=countsDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
#parse genecounts summary for longRNA
#create plot with labels above bars; plotly handles autoscaling

write("Plotting counts information", stderr())
#countByBiotype$Biotype <- factor(countByBiotype$Biotype, levels = unique(countByBiotype$Biotype)[order(countByBiotype$Count, decreasing = TRUE)]) 

pl <- plot_ly(countByBiotype,
	      x=~Biotype,
	      y=~as.numeric(gsub(",","",Count)),
	      text=~prettyNum(Count, big.mark = ",", scientific=FALSE),
	      textposition='outside',
	      type='bar')
plot_by_biotype = countByBiotype
plot_by_biotype$Count = as.numeric(gsub(",","",plot_by_biotype$Count))
pl <- plot_ly(plot_by_biotype, x=~Biotype, y=~Count, type='bar')

#countByBiotype$count <- prettyNum(countByBiotype$count, big.mark = ",", scientific=FALSE)

#render tables and plots on separate pages

#' `r if(exists("longRNAcounts")) { "## longRNA Counts" }`
#+ eval=exists("longRNAcounts"), echo=FALSE, fig.asp=1, fig.align="center", message=F
kable(longRNAcounts, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage

#' `r if(exists("countByBiotype")) { "## Gene Biotypes" }`
#+ eval=exists("countByBiotype"), echo=FALSE, fig.asp=1, fig.align="center", message=F
kable(countByBiotype, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#render plot
pl

#' \newpage

#' `r if(TRUE) { "# Pipeline Metadata" }`
#+ eval=TRUE, echo=FALSE, fig.asp=1, fig.align="center", message=F, results="asis", warn=F
kable(env, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


#' \newpage
#' `r if(TRUE) { "# Glossary of Terms" }`
#+ eval=TRUE, echo=FALSE, fig.asp=1, fig.align="center", message=F, results="asis", warn=F
kable(glossary, "latex", booktabs = T) %>%
	  kable_styling(latex_options = c("striped", "hold_position"))
