library(clusterProfiler)
library(pathways)
library(tidyverse)
library(igraph)
library(ggraph)
load("D:\\wdb\\packages\\kinase_metabolism-master\\pathways\\pathways.RData")
comp_dotplot(da@kegg_analyst$enrichKEGG$a)
da->lllll
lllll->da
da@kegg_analyst$enrichKEGG$a@result%>% dplyr::filter(org=="KEGG") ->y
da@kegg_analyst$enrichKEGG$d->d23
plot_funmap(y)



load(system.file("data", "Tes.Rdata",package = "pathways"))
 n@kegg_analyst$compareClusterResult%>%clusterProfiler::filter(org=="KEGG")%>%
 clusterProfiler::filter(
(Description%in%c((n@kegg_analyst$compareClusterResult%>%group_by(Description)%>%
                    summarise(n=n()>1))%>%dplyr::filter(n==T)%>%.$Description)|qvalue     <0.05)
)%>%as.data.frame()->f
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,start.degree=90)
plot_chor(f,metaProfiler = n)->v
#looking for ComplexHeatmap[https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html]
pd = packLegend(list =  v)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))











x<-readxl::read_xlsx(system.file("data", "Col-top50.xlsx", package = "pathways"))
y<-readxl::read_xlsx(system.file("data", "metadata.xlsx", package = "pathways"))%>% as.data.frame()%>% na.omit
rownames(y)<-y$Sample%>%str_extract_all("[0-9]+[0-9]")%>% unlist %>% paste("cid:",., sep = "")
y[,-1]->y
library(biomaRt)
mart <- useMart("ensembl","mmusculus_gene_ensembl")
gene_name<-getBM(attributes=c("ensembl_gene_id","entrezgene_id"),
                 filters = "ensembl_gene_id",values = x$Gene, mart = mart)
inner_join(gene_name,x,by=c("ensembl_gene_id"="Gene"))%>% na.omit->getkegg
getkegg$ensembl_gene_id[!is.na(getkegg$entrezgene_id)]<-getkegg$entrezgene_id[!is.na(getkegg$entrezgene_id)]
getkegg[,-2]->genes
rownames(genes)<-genes$ensembl_gene_id
genes[,-1]->genes
group<-c(rep("A",4),rep("B",5))
list(a=y%>%t%>%as.data.frame(),b=genes%>%t%>%as.data.frame())->data
data%>% setNames(c("a","b"))->data
names<-match(data[[1]]%>% rownames,data[[2]]%>% rownames)
data[[2]]<-data[[2]][names,]
pathways_analy(data = data,group = group,org = "10090",scale = T,import_model=c("degree"))->n
#save(n,file="Tes.Rdata")
# or



