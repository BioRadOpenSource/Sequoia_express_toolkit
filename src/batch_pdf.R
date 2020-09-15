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
#' `r {"# Summary Alignment Statistics"}`
#+ echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
kable(align_frame,  "latex", booktabs = T) %>% 
	kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage
#' `r {"# Summary Alignment Report"}`
#+ echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
kable(report_frame,  "latex", booktabs = T) %>%
	kable_styling(latex_options = c("striped", "hold_position"))

#' \newpage
#' `r if(dedupDirExists) {"# Summary Alignment Report"}`
#+ eval=dedupDirExists, echo=FALSE, fig.asp=1, fig.align="center", message=F, warn=F
kable(dedup_frame,  "latex", booktabs = T) %>% 
	kable_styling(latex_options = c("striped", "hold_position"))

