
## instead of `cellName`, users can use other features (e.g. `cancerType`)
get_GSEA<-function(data2,n,exponent = 1, minGSSize = 2, maxGSSize = 500,
                                      eps  = 1e-10,
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      verbose = TRUE,
                                      seed = FALSE,
                                      by = 'fgsea') {
  data <- n  %>%
  dplyr::group_by(ID, Compound_ID) %>%
  dplyr::summarise(Compound_ID = str_split(Compound_ID, '\\|')%>%unlist) %>% dplyr::ungroup()
GSEA_internal <- DOSE:::GSEA_internal(data2%>%.[order(- . )],  USER_DATA     = DOSE:::build_Anno(data, NA),exponent = 1,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps  = eps,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      verbose = verbose,
                                      seed = seed,
                                      by = by)

  GSEA_internal@result%>% dplyr::left_join(n[,1:4],by="ID")%>% dplyr::mutate(Description=name)%>% dplyr::select(!name)->GSEA_internal@result

return(GSEA_internal)
}


#' @import dply `%>%`


#' GSEA to pathways
#' @title GSE
#' @name GSEA
#' @param data clusterprofiler GSEA data
#' @param org tai_id 10090 is mouse
#' @return GSEA 
#' @export
#'
#' @examples
#'
##' \dontrun{ 
##'   system.file("data", "Tes.RData", package = "pathways")
##' n@data$a[1,]%>%c%>%unlist%>%abs->data2
##' gsea<-GESA(data=data2,org="10090")
##' # or
##' }



GESA<-  function(data,org="10090",exponent = 1, minGSSize = 2, maxGSSize = 500,
                                      eps  = 1e-10,
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH",
                                      verbose = TRUE,
                                      seed = FALSE,
                                      by = 'fgsea'){
  get_org(org)->n
  data%>%names->name
  if(name[1]%>%str_detect("cid:\\d")){
      cat("doing  Compound analyst \n")
                 return( get_GSEA(data,n, exponent = exponent,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps  = eps,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      verbose = verbose,
                                      seed = seed,
                                      by = by))
    }
    #par
   # pathways->kegg_pathways
    if(name[1]%>%str_detect("[a-zA-Z]")){
            cat("doing  Protacxns analyst \n")
            colnames(n)<- colnames(n)[c(7,5)]
                       return( get_GSEA(data,n, exponent = exponent,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps  = eps,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      verbose = verbose,
                                      seed = seed,
                                      by = by))
    }
     if(name[1]%>%str_detect("\\.")){
       cat("doing  Ecs analyst \n")
                   colnames(n)<- colnames(n)[c(8,5)]
                              return( get_GSEA(data,n, exponent = exponent,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps  = eps,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      verbose = verbose,
                                      seed = seed,
                                      by = by))
     }
     if(name[1]%>%str_remove_all("[0-9]")%>%{.==""}){
      cat("doing  genes analyst \n")
                         colnames(n)<- colnames(n)[c(6,5)]
                                    return( get_GSEA(data,n, exponent = exponent,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      eps  = eps,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      verbose = verbose,
                                      seed = seed,
                                      by = by))
     }
}


