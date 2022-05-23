## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",out.height = "50%",
  message = F,warning = F,tidy = T
)

## ----cars---------------------------------------------------------------------
library(pathways)
library(dplyr)
library(clusterProfiler)
library(pathways)
library(tidyverse)
library(igraph)
library(ggraph)
load(system.file("data", "Tes.Rdata",package = "pathways"))
 n@kegg_analyst$compareClusterResult%>%clusterProfiler::filter(org=="KEGG")%>%
 clusterProfiler::filter(
(Description%in%c((n@kegg_analyst$compareClusterResult%>%group_by(Description)%>%
                    summarise(n=n()>1))%>%
                    dplyr::filter(n==T)%>%.$Description)|qvalue     <0.05)
)%>%as.data.frame()->f
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,start.degree=90)
plot_chor(f,metaProfiler = n)->v
pd = packLegend(list =  v)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))


## ----ds-----------------------------------------------------------------------
comp_dotplot(n@kegg_analyst$enrichKEGG$a)
n@kegg_analyst$enrichKEGG$a->lllll
lllll->da
da@result%>% dplyr::filter(org=="KEGG") ->y
#plot_funmap(y)

