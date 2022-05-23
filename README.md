# pathways
easy to analy pathways
You'll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

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

<img src="help/figures/README-pressure-1.png" width="100%" height="50%"/>

In that case, don't forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

``` r
easy.clusterProfiler(n@kegg_analyst$enrichKEGG$a@result)
```

<img src="help/figures/README-unnamed-chunk-3-1.png" width="100%" height="50%"/>


    clusterProfiler::dotplot(n@kegg_analyst$enrichKEGG$a)

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" height="50%"/>

``` r
enrichplot::emapplot(enrichplot::pairwise_termsim(n@kegg_analyst$compareClusterResult))
```

<img src="/man/figures/README-unnamed-chunk-3-3.png" width="100%" height="50%"/>
