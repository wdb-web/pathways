# pathways
easy to analy pathways
You'll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

``` r
library(clusterProfiler)
library(pathways)
library(tidyverse)
library(igraph)
library(ggraph)

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

```


``` r
# kegg_pathway1(data=c("cid:5997","cid:65094","cid:5280335"))->da
# kegg_pathway1(data=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1"))->da
# kegg_pathway1(data=c("3.4.22.15","3.4.22.69","3.6.4.12","3.6.4.13"))->da
# kegg_pathway1(c("6921","6923","8453","8883","9039","9978","79699"))->da
# kegg_pathway1(a=c("6921","6923","8453","8883","9039","9978","79699"),b=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1"))->da

```


``` r
# load(system.file('data', 'Tes.Rdata',package = 'pathways'))
n@kegg_analyst$compareClusterResult %>%
    clusterProfiler::filter(Description %in% c((n@kegg_analyst$compareClusterResult %>%
        group_by(Description) %>%
        summarise(n = n() > 1)) %>%
        dplyr::filter(n == T) %>%
        .$Description) | qvalue < 0.05) %>%
    as.data.frame() -> f
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE, start.degree = 90)
plot_chor(f, metaProfiler = n) -> v
# looking for
# ComplexHeatmap[https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html]
pd = packLegend(list = v)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
```

<img src="vignettes/man/figures/README-cars-1.png" width="100%" height="50%"/>

In that case, don't forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.


``` r
comp_dotplot(n@kegg_analyst$enrichKEGG$a)
```

<img src="vignettes/man/figures/README-ds-1.png" width="100%" height="50%"/>


``` r
n@kegg_analyst$enrichKEGG$a@result%>% dplyr::filter(org=="KEGG") ->y
n@kegg_analyst$enrichKEGG$d@result->d23
plot_funmap(d23)
```

<img src="vignettes/man/figures/funmap.png" width="100%" height="50%"/>


``` r
easy.clusterProfiler(n@kegg_analyst$enrichKEGG$a@result)
```

<img src="docs/reference/easy.clusterProfiler-1.png" width="100%" height="50%"/>

``` r
 clusterProfiler::dotplot(n@kegg_analyst$enrichKEGG$a)
```


<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" height="50%"/>

``` r
enrichplot::emapplot(enrichplot::pairwise_termsim(n@kegg_analyst$compareClusterResult))
```

<img src="/man/figures/README-unnamed-chunk-3-3.png" width="100%" height="50%"/>
