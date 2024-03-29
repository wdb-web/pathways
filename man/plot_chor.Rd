% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/char_plot.R
\name{plot_chor}
\alias{plot_chor}
\title{plot_chor is plot metaProfiler data}
\usage{
plot_chor(
  f = f,
  metaProfiler,
  scale_height = 0.6,
  begin = 0,
  pathways_scale_color = structure(c(-log(0.05), -log(0.01)), names = c("red", "blue")),
  directional = 0,
  annotationTrack = NULL,
  line_scale_col = NULL,
  ...
)
}
\arguments{
\item{f}{using dplyr to filter data}

\item{metaProfiler}{metaProfiler data}

\item{pathways_scale_color}{pathways bar color}

\item{directional}{https://jokergoo.github.io/circlize_book/book/}

\item{annotationTrack}{https://jokergoo.github.io/circlize_book/book/}

\item{line_scale_col}{line color}

\item{...}{https://jokergoo.github.io/circlize_book/book/}
}
\value{
a list of length ,you can using packLegend and draw to add length
}
\description{
plot_chor is plot metaProfiler data
}
\examples{
\dontrun{
load(system.file("data", "Tes.Rdata",package = "pathways"))
 n@kegg_analyst$compareClusterResult@result\%>\%dplyr::filter(
Description\%in\%c((n@kegg_analyst$compareClusterResult\%>\%group_by(Description)\%>\%
                    summarise(n=n()>1))\%>\%dplyr::filter(n==T)\%>\%.$Description)|qvalue     <0.05
)\%>\%as.data.frame()->f
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,start.degree=90)
plot_chor(f,metaProfiler = n)->v
#looking for ComplexHeatmap[https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html]
pd = packLegend(list =  v)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
}
}
