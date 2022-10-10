#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate
#' project work
#' @importMethodsFrom DOSE enrichResult compareClusterResult
#' @import purrr
#' @import dplyr
#' @import data.table
#' @importFrom stringr str_split str_locate_all str_sub str_remove str_match
#' @importFrom methods setClass news
#' @importFrom data.table  fread
#' @importFrom ggplot2 `+` ggplot aes geom_point theme_bw

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

#import_model=c("betweenness","degree")
.get_import <- function(enrichKEGG,import_model="degree",kegg_org) {
  # n@kegg_analyst$enrichKEGG$a@result->d2
  if(is.null(enrichKEGG)){return(NULL)}
  enrichKEGG@result  ->d2
  d2%>%dplyr::filter(Count>1)->d2
  d2[,colnames(d2)!="import" ]->d2
  # ad<-""
  #.x<-"KEGG"
  ad<- names( table(d2$org)[table(d2$org)>1])%>%unique%>%purrr::map(function(.x){
    pathw<- tryCatch(  {d2[which(d2$org==.x),]->d23

      d23%>%dplyr::group_by(ID)%>%dplyr::summarise(geneID %>% str_split("/")%>%unlist())%>%
        setNames(c("ID", "Compound_ID" ))%>%
        ungroup()%>%as.data.table()->da

      if(sum(table(da$Compound_ID)>1)){
        d23%>%ungroup()%>%dplyr::select(ID,geneID)->d
        da[ , if (.N > 1)
          transpose(combn(ID, 2L, simplify = FALSE)), by = Compound_ID
        ][ , .(Sum = .N), keyby = .(Group.1 = V1, Group.2 = V2)]->l
        colnames(d)<-c("Group.1","n")

        full_join(l,d, by = "Group.1")%>%na.omit()%>%full_join(d,by=c("Group.2"="Group.1"))%>%na.omit()->l2
        apply(l2[,4:5],1,function(.x){.x%>% str_split("/")%>%unlist%>%unique%>% length })->l2$everysum
        data.frame(l2$Group.1,l2$Group.2,l2$Sum/l2$everysum)->l
        colnames(l)<-c("from","to","weight")
        l%>%dplyr::filter(weight>0.25)->l
        net_pc<-igraph::graph_from_data_frame(
          d=l[,1:2],
          directed=F)
        igraph::E(net_pc)$weight <- c(l[,3]);
        if(import_model[1]=="betweenness"){
          igraph::betweenness(net_pc,weights=l[,3],nobigint = F,directed=F,normalized=T)%>%
            as.data.frame()->import
        }
        if(import_model[1]=='degree'){
          a<-rbind(l[,c(1,3)],setNames(l[,2:3],names(l[,c(1,3)])))%>%ungroup%>%group_by(from)%>%dplyr::summarise(import=sum(weight/(length(unique(c(l[,1:2]%>%unlist)))-1)))%>%as.data.frame()
          a[,-1]%>%as.data.frame()->import
          rownames(import)<- a$from ;
        }
        if(import_model[1]=="cloness"){
          igraph::closeness(net_pc,normalized=F)%>%as.data.frame()->import
        }
        colnames(import)[1]<-"import"
        import$ID <-import%>%rownames()
        enrichKEGG@keytype<-"kegg"
        d23%>%full_join(import,by = "ID")->d23;
        d23[d23$ID!="1",]->d23
        return(d23%>%.[,c(which(colnames(d23)!="org"),which(colnames(d23)=="org"))])
      }
    }, error = function(e) return({
      d2[which(d2$org==.x),]->d23
      d23$import<-0
      return(d23)
    }))
    return(pathw)
  }
  )%>%purrr::reduce(rbind)
  ad%>%data.frame ->enrichKEGG@result
  # enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
  rownames(enrichKEGG@result)<-enrichKEGG@result$ID

  return(enrichKEGG)
}
get_all_p <- function(g){
  g%>%dplyr::group_by(Cluster)%>%dplyr::summarise(w=BgRatio%>%str_split("/")%>%unlist()%>%.[2])%>%ungroup()->n
  g%>%dplyr::select(Cluster,ID ,pvalue)%>%dplyr::group_split(Cluster)%>%purrr::map(~.x%>%.[,-1])%>%purrr::reduce(full_join,by="ID")->k
  k[is.na(k)]<-1
  .get_p_combine(k[1,-1],method = 'fisher',w =n[,2])$p_value
  map_dfr(k[,-1]%>%t%>%as.data.frame,~.get_p_combine(.x,method = 'fisher',w =n[,2])$p_value)%>%t%>%as.data.frame()->p
  #rownames(p)<-k$ID
  p.adjust(p$V1,method = "BH")->p.adjust
  p.adjust(p$V1,method = "fdr")->qvalue
  data.frame( pvalue=p,p.adjust=p.adjust,qvalue=qvalue)->p.data
  p.data$ID<-k$ID
  g%>%dplyr::group_by(ID,Description)%>%dplyr::summarise(geneID=geneID%>%
                                                           paste(sep = "/",collapse = "/"),Count=sum(Count))->data



  full_join(data,p.data)->data
  colnames(data)[5]<-"pvalue"
  return(data)
}
.get_p_combine <- function(p, method = c("fisher", "SL", "MG", "tippett"), w = NULL) {
  method <- match.arg(method, choices = c("fisher", "SL", "MG", "tippett"))
  p <- p[!is.na(p)]
  n <- length(p)
  if (max(p) - .Machine$double.ep > 1) {
    stop("invalid input, p > 1")
  }

  if (n == 0) {
    warning("vector of p-values is empty")
    if (method == "fisher") {
      return(list(statistic = NA, p_value = NA, method = "Fisher (1932)",
                  statistic_name = "Xsq"))
    } else if (method == "SL") {
      return(list(
        statistic = NA, p_value = NA,
        method = "Stouffer (1949), Liptak (1958)",
        statistic_name = "Z"
      ))
    } else if (method == "MG") {
      return(list(
        statistic = NA, p_value = NA,
        method = "Mudholkar and George (1979)",
        statistic_name = "L"
      ))
    } else if (method == "tippett") {
      return(list(
        statistic = NA, p_value = NA,
        method = "Tippett (1931)",
        statistic_name = "p_min"
      ))
    }
  } else {
    if (method == "fisher") {
      Xsq <- -2 * sum(log(p))
      p_val <- stats::pchisq(Xsq, df = 2 * n, lower.tail = FALSE)
      return(list(
        statistic = Xsq, p_value = p_val, method = "Fisher (1932)",
        statistic_name = "Xsq"
      ))
    } else if (method == "SL") {
      if (is.null(w)) {
        w <- rep(1, n) / n
      } else {
        if (length(w) != n) {
          stop("length of p and w must be equal")
        }
      }
      Zi <- stats::qnorm(1 - p)
      Z <- sum(w * Zi) / sqrt(sum(w^2))
      p_value <- 1 - stats::pnorm(Z)
      return(list(
        statistic = Z, p_value = p_value,
        method = "Stouffer (1949), Liptak (1958)",
        statistic_name = "Z"
      ))
    } else if (method == "MG") {
      L <- sum((-1) * log(p / (1 - p)))
      p_value <- stats::pt(L * sqrt((15 * n + 12) / (pi^2 * n * (5 * n + 2))),
                           df = 5 * n + 4, lower.tail = FALSE
      )
      return(list(
        statistic = L, p_value = p_value,
        method = "Mudholkar and George (1979)",
        statistic_name = "L"
      ))
    } else if (method == "tippett") {
      return(list(
        statistic = min(p), p_value = 1 - (1 - min(p))^n,
        method = "Tippett (1931)",
        statistic_name = "p_min"
      ))
    }
  }
}
get_org<-function(org=org) {
  ls<-ls(envir=.GlobalEnv)
  get_class <- function(.x) {
    tryCatch(
      return(eval(parse(text = .x))%>%class%>%{"metaProfiler"%in%.})  ,
      error = function(e){return(F)
      },finally = {
      })
  }
  ls[purrr::map(ls,~get_class(.x))%>%unlist()]->l

  get_org_1 <- function(l) {
    l%>%parse(text = .)%>%eval%>%.@org_organism%>%return()
  }
  if(length(l)>1){purrr::map(l,get_org_1)%>%unlist->k
    if(org%in%k){
      return(parse(text = l[k==org]%>%.[1])%>%eval%>%.@org)
    }
  }
  if(length(l)==1){
    l%>%parse(text = .)%>%eval%>%.@org_organism->k
    if(org%in%k){
      return(parse(text = l[k==org])%>%eval%>%.@org)
    }
  }


  return(.get_all_org(org))

}
#.get_kegg_org("9606")->data
#org=9606
.get_all_org<-function(org){
  options(dplyr.summarise.inform = FALSE)
  y<- readLines("https://www.genome.jp/kegg-bin/download_htext?htext=br08610&format=htext&filedir=")
  kegg_org<-y[y%>% stringr::str_detect(paste("TAX:",org,sep=""))%>% which()%>%{.+1}]%>%str_sub(start = 2, end = -1L) %>% str_split("\\W")%>%.[[1]]%>%.[.!=""]%>%.[1]
  #options(download.file.method="libcurl")
  url =(paste("https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22where%22:{%22ands%22:[{%22taxid%22:%22",org,"%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22TaxID_",org,"_pathway%22}",sep=""))
 cat("download by" ,url)
  d = read_csv(url)

  datad<-d%>%dplyr::summarise(ID=pwacc,name,group=category,org=source ,
                              Compound_ID=cids%>%str_split("\\|")%>%lapply(function(x)x%>%paste("cid:",.,sep="")%>%unlist()%>%paste(sep = "\\|",collapse = "|")),
                              Gene_id=geneids,Protacxns_id=protacxns,Ecs_id=ecs)
  datad$Compound_ID[str_detect(datad$Compound_ID,"NULL")]<-"NULL"
  kegg_pathways<-.get_kegg_org(kegg_org)
  kegg_pathways%>%rbind(datad,.)->pathways
}
.get_kegg_org<- function(org=org) {
  name_rmends <- function(.x) {
    .x %>%str_locate_all(" - ")%>%.[[1]]->jk
    jk[nrow(jk),]%>%.[1]%>%{.-1}%>%str_sub(.x, end =.)%>%return
  }# getname
  x<-readLines("https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=")
  x%>%str_detect("!" )%>%which()->s
  x[(s[1]+1):(s[2]-1)]%>%str_replace("<b>","  ")%>%str_remove("</b>")->a
  x[(s[1]+1):(s[2]-1)]%>%str_replace("<b>","  ")%>%str_remove("</b>")%>%as.data.frame()->das
  das$A<-NA
  das$B<-NA
  das$C<-NA
  a[str_detect(a,"^A")%>%which()]%>%str_remove("A  ")->das$A[str_detect(a,"^A")%>%which()]
  a[str_detect(a,"^B")%>%which()]%>%str_remove("B  ")->das$B[str_detect(a,"^B")%>%which()]
  a[str_detect(a,"^C")%>%which()]%>%str_remove("C  ")%>%str_remove("  [0-9][0-9][0-9][0-9][0-9]  ")->das$C[str_detect(a,"^C")%>%which()]
  das %>%
    fill(A, .direction = "up") %>%
    fill(B, .direction = "up") %>%na.omit()->s
  x1<-fread("http://rest.kegg.jp/list/pathway",head=F)
  x2<-fread(paste("http://rest.kegg.jp/list/pathway/",org,sep=""),head=F)

  x2$V2%>%purrr::map(~name_rmends(.x))%>%unlist()%>%cbind(x2[,1])%>% setNames(c("V2","name"))  ->pathways# getname

  y<-fread(paste("http://rest.kegg.jp/link/",org,"/pathway",sep=""),head=F)%>%setNames(c("name","V2"))
  x<-fread(paste("http://rest.kegg.jp/conv/",org,"/uniprot",sep=""),head=F)
  par=full_join(pathways,x1,by="V2")%>%na.omit%>%full_join( full_join(x,y,by="V2")%>%na.omit,by="name")%>%
    .[,c(2,1,4)]%>%setNames(c("ID","name","Protacxns_id"))%>% dplyr::group_by(ID ,name )%>%
    dplyr::summarise(Protacxns_id=Protacxns_id%>% str_remove("^up:")%>% paste( sep = "|", collapse = "|" )  )%>%
    ungroup%>% dplyr::mutate(ID= ID %>% str_remove("^path:") )%>%inner_join(s[,c(2,4)],by=c("name"="C"))%>%
    dplyr::summarise( ID = ID,name=name,group=A,org="KEGG",Protacxns_id=Protacxns_id)

  gene=full_join(pathways,x1,by="V2")%>%na.omit%>%full_join( y,by="name")%>%
    .[,c(2,1,4)]%>%setNames(c("ID","name","Gene_id"))%>% dplyr::group_by(ID ,name )%>%
    dplyr::summarise(Gene_id=Gene_id%>% str_remove(.,"^(.+?)\\:") %>% paste( sep = "|", collapse = "|" ))%>%
    ungroup%>% dplyr::mutate(ID= ID %>% str_remove("^path:") )%>%inner_join(s[,c(2,4)],by=c("name"="C"))%>%
    dplyr::summarise( ID = ID,name=name,group=A,org="KEGG",Gene_id=Gene_id)


  k<-fread("http://rest.kegg.jp/link/pathway/cpd",header = F)%>% setNames(c("Compound_ID","V1"))
  b=fread("http://rest.kegg.jp/list/pathway",head=F)%>%inner_join(k,by = "V1")%>% inner_join(pathways,by="V2")
  y<-readxl::read_xlsx(system.file("data", "PubChem_and_kegg.xlsx", package = "pathways"))
  cdp<-b%>% dplyr::mutate( ID=name%>%str_remove("path:"),name=V2,
                           KEGGID=Compound_ID%>%str_remove( "cpd:"))%>%inner_join(y,by="KEGGID") %>%inner_join(s,by=c("name"="C"))%>%
    dplyr::summarise( ID = ID,name=name,group=A,org="KEGG",Compound_ID=cid)%>% tibble %>%
    group_by(ID,name,group,org) %>% dplyr::summarise(Compound_ID=paste(Compound_ID,sep="|",collapse = "|"))%>%
    ungroup

  k<-fread("http://rest.kegg.jp/link/pathway/enzyme",header = F)%>% setNames(c("Ecs_id","V1"))%>%.[str_detect(.$V1,"^path:map"), ]
  Ecs=fread("http://rest.kegg.jp/list/pathway",head=F)%>%inner_join(k,by = "V1")%>% inner_join(pathways,by="V2")%>%
    dplyr::mutate( ID=name%>%str_remove("path:"),name=V2,
                   KEGGID=Ecs_id%>%str_remove( "ec:"))%>%inner_join(s,by=c("name"="C"))%>%
    dplyr::summarise( ID = ID,name=name,group=A,org="KEGG",Ecs_id=KEGGID)%>% tibble %>%
    group_by(ID,name,group,org) %>% dplyr::summarise( Ecs_id=paste(Ecs_id,sep="|",collapse = "|"))%>%
    ungroup

  fuu<-full_join(cdp,gene,by = c("ID", "name", "group", "org"))%>%
    full_join(par,by = c("ID", "name", "group", "org"))%>%
    full_join(Ecs,by = c("ID", "name", "group", "org"))%>%
    dplyr::mutate(name, group,org,Compound_ID,Gene_id,  Protacxns_id,Ecs_id)
  fuu[is.na(fuu)]<-"NULL"

  fuu %>%return

  #  ,model="ec"
  #  fread(paste("http://rest.kegg.jp/list/pathway/",org,sep = ""),header = F)->k
  # name_rmends <- function(.x) {
  #  .x %>%str_locate_all(" - ")%>%.[[1]]->jk
  # jk[nrow(jk),]%>%.[1]%>%{.-1}%>%str_sub(.x, end =.)%>%return
  #}# getname
  #k$V2%>%purrr::map(~name_rmends(.x))%>%unlist()%>%cbind(k[,1])->pathways# getname

  #fread("http://rest.kegg.jp/list/pathway/map",header = F)->enknowpathways
  # colnames(pathways)[1]<-"path"
  # colnames(enknowpathways)[2]<-"path"
  # dplyr::left_join(pathways,enknowpathways,by=colnames(enknowpathways)[2])->k2
  # colnames(k2)[2:3]<-c("ko","pathway")
  #  k2$pathway%>%str_remove("path:")%>%.[1]
  # fread(paste("http://rest.kegg.jp/link/",model, "/map00010",sep = ""),header = F,verbose=F,showProgress=F)
  # pb <- dplyr::progress_estimated(length(k2$pathway))
  # get_kegg_comp <- function(.x) {
  #  pb$tick()$print()

  # tryCatch({
  #     return(fread(paste("http://rest.kegg.jp/link/",model,"/",
  #                        .x,sep = ""),header = F,verbose=F,showProgress=F))
  #   },
  #  error = function(e){return(NULL)
  #  },finally = {
  #  })
  # }

  #   k2$pathway%>%str_remove("path:")%>%purrr::map(~get_kegg_comp(.x))%>%
  #   purrr::reduce(bind_rows)%>%setNames(c("pathway","Compound_ID"))->kegg

  #   inner_join(k2,kegg)->b

  #   inner_join(b,s[,c(2,4)],by=c("path"="C"))->b
  #  if(model=="cpd"){
  #    y<-readxl::read_xlsx(system.file("data", "PubChem_and_kegg.xlsx", package = "pathways"))
  #   b %>% dplyr::mutate( ID=ko%>%str_remove("path:"),name=path,
  #     KEGGID=Compound_ID%>%str_remove( paste(model,":",sep="")))%>%inner_join(y,by="KEGGID")%>%
  #     dplyr::summarise( ID = ID,name=path,group=A,org="KEGG",Compound_ID=cid)%>%as.data.frame%>%
  #     group_by(ID) %>% dplyr::summarise( ID = ID,name,group,org,Compound_ID=paste(Compound_ID,sep="|",collapse = "|"))%>%
  #     ungroup%>%return
  # }
  #   if(model=="ec"){
  #   b %>% dplyr::mutate( ID=ko%>%str_remove("path:"),name=path,
  #     KEGGID=Compound_ID%>%str_remove( paste(model,":",sep="")))%>%
  #     dplyr::summarise( ID = ID,name=path,group=A,org="KEGG",Compound_ID=KEGGID)%>%as.data.frame%>%
  #     group_by(ID) %>% dplyr::summarise( ID = ID,name,group,org,Ecs_id=paste(Compound_ID,sep="|",collapse = "|"))%>%
  #     ungroup%>%return
  # }
}

get_anong<-function(org=org) {
  ls<-ls(envir=.GlobalEnv)
  get_class <- function(.x) {
    tryCatch(
      return(eval(parse(text = .x))%>%class%>%{"metaProfiler"%in%.})  ,
      error = function(e){return(F)
      },finally = {
      })
  }
  ls[purrr::map(ls,~get_class(.x))%>%unlist()]->l

  get_org_1 <- function(l) {
    l%>%parse(text = .)%>%eval%>%.@org_organism%>%return()
  }
  if(length(l)>1){purrr::map(l,get_org_1)%>%unlist->k
    if(org%in%k){
      return(parse(text = l[k==org]%>%.[1])%>%eval%>%.@anong) }}
  if(length(l)==1){
    l%>%parse(text = .)%>%eval%>%.@org_organism->k
    if(org%in%k){
      return(parse(text = l[k==org])%>%eval%>%.@anong)
    }
  }
  return(NULL)
}

.get_kegg_p<- function(value=data,kegg_pathways=kegg_pathways,p_model=p_model,
                       p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods) {
  value%>%unique()->data
  #kegg_pathways%>%as.data.frame()->match.df
  #orgline<-"KEGG"
  ananan<-purrr::map(kegg_pathways$org%>%unique,function(orgline){
    kegg_pathways[kegg_pathways$org==orgline,]->l
    l[!duplicated(l$name),]->l
    l%>%dplyr::filter(Compound_ID!="NULL")%>%.$Compound_ID%>%str_split("\\|")%>%unlist%>%unique%>%length->M
    sum( data%in%c(l$Compound_ID%>%str_split("\\|")%>%unlist))->m
    l%>%group_by(ID,org)%>%dplyr::filter(Compound_ID!="NULL")%>%
      summarise(n=sum( data%in%c(Compound_ID%>%str_split("\\|")%>%unlist)),N=Compound_ID%>%str_split("\\|")%>%unlist%>%unique%>%length,
                m=m,M=M,Compound_ID=paste(unique(data[data%in%c(Compound_ID%>%str_split("\\|")%>%unlist)]),collapse = "/"))%>%ungroup()%>%
      dplyr::filter(n>0)->n
    fisher<-function(data,model=c("phyper","fisher")){
      data -> d
      gene.not.interest=data.frame(c(d[3]-d[1]),c(d[4]-d[2]-d[3]+d[1])) %>% t;
      gene.in.interest=data.frame(d[1],d[3]) %>% t
      d <- data.frame(gene.not.interest,gene.in.interest)
      if(model[1]=="fisher"){return(fisher.test(d, alternative ="greater" )[["p.value"]])}
      # k -》 富集到通路，n 总的差异基因，N 所有基因，M 这条通路基因
      data -> d
      if(model[1]=="phyper"){return(phyper(d[1]-1,d[3], d[4]-d[3], d[2], lower.tail=F))}
      cat("using model is ","phyper or","fisher")
    }
    apply(n[,c(3:6)],1, fisher,model=p_model)->n$pvalue
    p.adjust(n$pvalue,method =p.adjust.methods)->n$p.adjust
    p.adjust(n$pvalue,method =qvalue.methods)->n$qvalue
    inner_join(kegg_pathways[,1:4]%>%unique,n,by = c("ID", "org"))->data
    data$GeneRatio<-paste(data$n,data$m,sep = "/");data$BgRatio<-paste(data$N,data$M,sep = "/");data$RichFactor<-paste(data$n,data$N,sep = "/")
    data%>%dplyr::select(ID=ID,Description=name,GeneRatio,BgRatio,RichFactor,
                         pvalue,p.adjust,qvalue,geneID=Compound_ID,Count=n,org,group)%>%as.data.frame()->data
    rownames(data)<-data$ID
    data%>%return }) %>% do.call(bind_rows,.)%>% return

}


#name<-data=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1")
.get_kegg_analyst<-function(name="list Compound_ID or gean or ec or par", org=org,
                            p_model=p_model,kegg_pathways=kegg_pathways,
                            p.adjust.methods=p.adjust.methods,
                            qvalue.methods="fdr",enrichKEGG=kk,import_model=import_model,everyorg=everyorg){
  if(name[1]%>%str_detect("cid:\\d")){
    cat("doing  Compound analyst \n")
    name%>%unique()->f
    c(kegg_pathways$Compound_ID%>% str_split("\\|")%>%unlist )->bg
    if(!sum(f%in%bg)){return(NULL)}
    .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                p_model=p_model,
                p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
    rownames(y)<-y$ID
    enrichKEGG@result<-y
    if(is.null( (enrichKEGG))) return(enrichKEGG)
    if( ! nrow(enrichKEGG@result)>1) return(enrichKEGG)
    enrichKEGG@gene<-name
    # enrichKEGG@universe<-kegg_pathways$Compound_ID%>%unique()
    enrichKEGG@organism<-org
    enrichKEGG@pvalueCutoff<-1
    enrichKEGG@qvalueCutoff<-1
    enrichKEGG@pAdjustMethod<-p.adjust.methods
    enrichKEGG@universe<-""
    #n@kegg_analyst$enrichKEGG$b->enrichKEGG
    .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
    enrichKEGG@keytype<- "Compound"
    enrichKEGG@ontology<- "Compound"
    enrichKEGG@readable<- FALSE
    enrichKEGG@result[order(enrichKEGG@result$pvalue),]->enrichKEGG@result
    return(enrichKEGG)
  }
  #par
  # pathways->kegg_pathways
  if(name[1]%>%str_detect("[a-zA-Z]")){
    cat("doing  Protacxns analyst \n")

    colnames(kegg_pathways)[c(5,7)]<-colnames(kegg_pathways)[c(7,5)]
    name%>%unique()->f
    c(kegg_pathways$Compound_ID%>% str_split("\\|")%>%unlist )->bg
    if(!sum(f%in%bg)){return(NULL)}
    .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                p_model=p_model,
                p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
    rownames(y)<-y$ID
    enrichKEGG@result<-y
    if(is.null( (enrichKEGG))) return(enrichKEGG)
    if( ! nrow(enrichKEGG@result)>1) return(enrichKEGG)
    enrichKEGG@gene<-name
    # enrichKEGG@universe<-kegg_pathways$Compound_ID%>%unique()
    enrichKEGG@organism<-org
    enrichKEGG@pvalueCutoff<-1
    enrichKEGG@qvalueCutoff<-1
    enrichKEGG@pAdjustMethod<-p.adjust.methods
    enrichKEGG@universe<-""
    #n@kegg_analyst$enrichKEGG$b->enrichKEGG
    .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
    enrichKEGG@keytype<- "par"
    enrichKEGG@ontology<- "par"
    enrichKEGG@readable<- FALSE
    enrichKEGG@result[order(enrichKEGG@result$pvalue),]->enrichKEGG@result
    return(enrichKEGG)
  }
  if(name[1]%>%str_detect("\\.")){
    cat("doing  Ecs analyst \n")

    colnames(kegg_pathways)[c(5,8)]<-colnames(kegg_pathways)[c(8,5)]
    .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                p_model=p_model,
                p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
    name%>%unique()->f
    c(kegg_pathways$Compound_ID%>% str_split("\\|")%>%unlist )->bg
    if(!sum(f%in%bg)){return(NULL)}
    .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                p_model=p_model,
                p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
    rownames(y)<-y$ID
    enrichKEGG@result<-y
    if(is.null( (enrichKEGG))) return(enrichKEGG)
    if( ! nrow(enrichKEGG@result)>1) return(enrichKEGG)
    enrichKEGG@gene<-name
    # enrichKEGG@universe<-kegg_pathways$Compound_ID%>%unique()
    enrichKEGG@organism<-org
    enrichKEGG@pvalueCutoff<-1
    enrichKEGG@qvalueCutoff<-1
    enrichKEGG@pAdjustMethod<-p.adjust.methods
    enrichKEGG@universe<-""
    #n@kegg_analyst$enrichKEGG$b->enrichKEGG
    .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
    enrichKEGG@keytype<- "Ecs"
    enrichKEGG@ontology<- "Ecs"
    enrichKEGG@readable<- FALSE
    .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0

    enrichKEGG@result[order(enrichKEGG@result$pvalue),]->enrichKEGG@result
    return(enrichKEGG)
  }
  if(name[1]%>%str_remove_all("[0-9]")%>%{.==""}){
    cat("doing  genes analyst \n")

    colnames(kegg_pathways)[c(5,6)]<-colnames(kegg_pathways)[c(6,5)]
    name%>%unique()->f
    c(kegg_pathways$Compound_ID%>% str_split("\\|")%>%unlist )->bg
    if(!sum(f%in%bg)){return(NULL)}
    .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                p_model=p_model,
                p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
    rownames(y)<-y$ID
    enrichKEGG@result<-y
    if(is.null( (enrichKEGG))) return(enrichKEGG)
    if( ! nrow(enrichKEGG@result)>1) return(enrichKEGG)
    enrichKEGG@gene<-name
    # enrichKEGG@universe<-kegg_pathways$Compound_ID%>%unique()
    enrichKEGG@organism<-org
    enrichKEGG@pvalueCutoff<-1
    enrichKEGG@qvalueCutoff<-1
    enrichKEGG@pAdjustMethod<-p.adjust.methods
    enrichKEGG@universe<-""
    #n@kegg_analyst$enrichKEGG$b->enrichKEGG
    .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
    enrichKEGG@keytype<- "Gene"
    enrichKEGG@ontology<- "Gene"
    enrichKEGG@readable<- FALSE
    enrichKEGG@result[order(enrichKEGG@result$pvalue),]->enrichKEGG@result
    return(enrichKEGG)
  }

}


# kegg_pathway1(data=c("cid:5997","cid:65094","cid:5280335"))->da
# kegg_pathway1(data=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1"))->da
# kegg_pathway1(data=c("3.4.22.15","3.4.22.69","3.6.4.12","3.6.4.13"))->da
#kegg_pathway1(c("6921","6923","8453","8883","9039","9978","79699"))->da

kegg_pathway1 <- function(data=c("list or data"),org="9606",p_model=c("phyper","fisher"),
                          p.adjust.methods="holm",import_model=c("degree","betweenness")) {
  options(dplyr.summarise.inform = FALSE)
  pathways<-get_org(org)
  pathways$org<- factor(pathways$org, level=c("KEGG",pathways$org%>%unique%>%.[.!="KEGG"]%>% as.character()) )
  pathways<- pathways[order(pathways$org),]%>% .[!duplicated(pathways$name) ,]
  s <- new("metaProfiler",org_organism =org)
  s@org<-pathways%>%as.data.frame()
  # get_anong(org)%>%.[[1]]->bb
  # bb->>bb
  #s@universe<-kegg_pathways$Compound_ID%>%unique()
  kk <- new("enrichResult",organism=org,ontology="ALL")
  if(is.list(data)&(length(data)==1)){
    .get_kegg_analyst(name = data[[1]],org = org,kegg_pathways=pathways,import_model=import_model,
                      p_model=p_model,p.adjust.methods=p.adjust.methods,everyorg = pathways,
                      enrichKEGG=kk)%>%list->s@kegg_analyst
    names(s@kegg_analyst)<-"enrichKEGG"
    return(s)
  }
  if(!is.list(data)){
    .get_kegg_analyst(name = data,org = org,kegg_pathways=pathways,
                      p_model=p_model,everyorg = pathways,
                      p.adjust.methods=p.adjust.methods,import_model=import_model,
                      enrichKEGG=kk)%>%list->s@kegg_analyst
    names(s@kegg_analyst)<-"enrichKEGG"
    return(s)
  }
  #data[[1]]->name
  # data%>% lapply(function(.x){.get_kegg_analyst(name = .x,org = org,kegg_pathways=kegg_pathways,import_model=import_model,enrichKEGG=kk,everyorg = kegg_pathways,
  #                              p_model=p_model,p.adjust.methods=p.adjust.methods)})->y
  y <- data %>% purrr::map(function(.x) {
    tryCatch( return( .get_kegg_analyst(name = .x, org = org, kegg_pathways = pathways,
                                        import_model = import_model, enrichKEGG = kk, everyorg = pathways,
                                        p_model = p_model, p.adjust.methods = p.adjust.methods)), error = function(e) return(NA), finally = print("OK"))
  })

  num <- y %>% lapply(is.null) %>% unlist %>% {
    !.
  } %>% which
  y <- num %>% purrr::map(function(.x) {
    y[[.x]]
  })
  x <- new("compareClusterResult")
  Result <- function(.x, .y) {
    if (is.null(.x)) {
      return(NULL)
    }
    if (is.na(.x)) {
      return(NULL)
    }
    if (ncol(.x@result) == 9) {
      .x@result$import <- ""
    }
    data.frame(Cluster = .y, .x@result)
  }
  g <- purrr::map2(y, names(y), ~Result(.x, .y)) %>% do.call(rbind,
                                                             .)
  x@compareClusterResult <- g
  gene <- function(.x) {
    if (is.null(.x)) {
      return(NULL)
    }
    if (is.na(.x)) {
      return(NULL)
    }
    .x@gene
  }
  x@geneClusters <- purrr::map(y, ~gene(.x))
  b <- factor(x@compareClusterResult$Cluster)
  x@keytype <- "id"
  x@compareClusterResult$Cluster <- factor(x@compareClusterResult$Cluster)
  x@compareClusterResult <- x@compareClusterResult %>% as.data.frame()
  x@fun = "enrichDO"
  s@kegg_analyst <- list(enrichKEGG = y, compareClusterResult = x,
                         every = get_all_p(g = x@compareClusterResult))
  return(s)
}

comp_dotplot<- function(data) {
  data@result%>%ggplot(aes(y=-log10(qvalue),
                           x=import,col=import,size=-log10(qvalue),label=Description))+
    geom_point()+ ggsci::scale_color_material("orange")+
    ggrepel::geom_text_repel()+theme_bw()+ggplot2::scale_size(range = c(4,8) )+
    ggplot2::ylim(0,3)
}

#kegg_pathway1(data)
#.libPaths("D:/wdb/R-4.1.2/library")
