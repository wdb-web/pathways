% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEA.R
\name{GSEA}
\alias{GSEA}
\alias{GESA}
\title{GSE}
\usage{
GESA(
  data,
  org = "10090",
  exponent = 1,
  minGSSize = 2,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)
}
\arguments{
\item{data}{clusterprofiler GSEA data}

\item{org}{tai_id 10090 is mouse}
}
\value{
GSEA
}
\description{
GSE
}
\examples{

\dontrun{ 
  system.file("data", "Tes.RData", package = "pathways")
n@data$a[1,]\%>\%c\%>\%unlist\%>\%abs->data2
gsea<-GESA(data=data2,org="10090")
# or
}
}
