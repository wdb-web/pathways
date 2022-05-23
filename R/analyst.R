#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @importMethodsFrom clusterProfiler enrichResult compareClusterResult
#' @import purrr
#' @import dplyr
#' @importFrom stringr str_split str_locate_all str_sub str_remove str_match
#' @importFrom methods setClass news 
setClass("metaProfiler",
         slots = c(
           kegg_analyst = "list",
           data_mixOmics_analyst = "list",
           data_mean="list",
           org="data.frame",
           org_organism="character",
           data="list",
           group="data.frame"
         ),
         prototype = list(
           kegg_analyst = list(),
           data_mixOmics_analyst = list(),
           data_mean=list(),
           org=data.frame(),
           org_organism=character(),
           data=list(),
           group=data.frame()
         )
)
#data%>%purrr::map(colnames)%>%.[[2]]->data
#' pathways_analy to pathways and analy
#'
#' @param data a list of data
#' @param group data group
#' @param org looking for https://www.genome.jp/kegg/catalog/org_list.html,using fread to keggAPI
#' @param fisher.alternative alternative  of fisher.test
#' @param scale scale of data ,It well to .get_data_analyst_plot
#' @param p.adjust.methods  looking for  p.adjust
#' @param import_model import_model=c("betweenness","degree")
#' @return metaProfiler data
#' @export
#'
#' @examples
#'
##' \dontrun{
##' 
##' 
##' x<-readxl::read_xlsx(system.file("data", "Col-top50.xlsx", package = "pathways"))
##' y<-readxl::read_xlsx(system.file("data", "metadata.xlsx", package = "pathways"))%>% as.data.frame()%>% na.omit

##' rownames(y)<-y$Sample%>%str_extract_all("[0-9]+[0-9]")%>% unlist %>% paste("cid:",., sep = "")
##' y[,-1]->y
##' library(biomaRt)
##' mart <- useMart("ensembl","mmusculus_gene_ensembl")
##' gene_name<-getBM(attributes=c("ensembl_gene_id","entrezgene_id"),
##'                  filters = "ensembl_gene_id",values = x$Gene, mart = mart)
##' inner_join(gene_name,x,by=c("ensembl_gene_id"="Gene"))%>% na.omit->getkegg
##' getkegg$ensembl_gene_id[!is.na(getkegg$entrezgene_id)]<-getkegg$entrezgene_id[!is.na(getkegg$entrezgene_id)]
##' getkegg[,-2]->genes
##' rownames(genes)<-genes$ensembl_gene_id
##' genes[,-1]->genes
##' group<-c(rep("A",4),rep("B",5))
##' list(a=y%>%t%>%as.data.frame(),b=genes%>%t%>%as.data.frame())->data
##' data%>% setNames(c("a","b"))->data
##' names<-match(data[[1]]%>% rownames,data[[2]]%>% rownames)
##' data[[2]]<-data[[2]][names,]
##' pathways_analy(data = data,group = group,org = "10090",scale = T)->n

##' 
##' # or
##' load(system.file("data", "Tes.Rdata",package = "pathways"))
##' }
##' 
pathways_analy <- function(data,group,org="hsa",p_model=c("phyper","fisher"),
                           scale=T,
                 p.adjust.methods="holm",import_model=c("betweenness","degree")) {
  if(is.list(data)){
    data%>%purrr::map(colnames)%>%
    kegg_pathway1(data = . ,org = org,p.adjust.methods=p.adjust.methods,p_model = p_model,
                  import_model=import_model)->y}else{
    kegg_pathway1(data =colnames(data) ,org = org,p.adjust.methods=p.adjust.methods,p_model = p_model,
                       import_model=import_model)->y}
  .get_data_analyst_plot( data=data,group = group,scale = scale)->h
  slot(y, "data_mixOmics_analyst") <- list(h$analyst_data)
  slot(y, "data_mean") <-h$men
  slot(y, "data") <-data
  slot(y, "group") <-group%>%as.data.frame()
  
  return(y)
}
