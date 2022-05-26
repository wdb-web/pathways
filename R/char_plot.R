#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @import circlize
#' @import dplyr
#' @importFrom ComplexHeatmap Legend
#' @import purrr
#' @importFrom randomcoloR distinctColorPalette 


plot_add <- function(f,begin=begin,scale_height=scale_height,scale_color=scale_color,
                     metaProfiler=metaProfiler,group2=group2) {
  metaProfiler->n
  scale_height<-scale_height
  grt <- function(data,begin,scale_height=0.4,scale_color=structure(c(-log(0.05),-log(0.01)),names =c("red","blue"))) {
    circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                 y = c( (begin-1)*scale_height, (begin-1)*scale_height),
                 sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                 y = c( (begin-1)*scale_height, (begin-1)*scale_height)+0.5*scale_height*0.95,
                 sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                 y = c( (begin-1)*scale_height, (begin-1)*scale_height)+1*scale_height*0.95,
                 sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.text(-0.5,track.index = 1,y=(begin-1)*scale_height+scale_height/2,labels =data$Cluster[1] ,
                sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1]
    )
    circos.text(x=rep(1.5,3),track.index = 1,y=c((begin-1)*scale_height+scale_height/2,(begin-1)*scale_height+scale_height*0.8,
                                                 (begin-1)*scale_height+0.05 ),labels = c(max(data$Count)/2,max(data$Count),0),
                sector.index = get.all.sector.index()[get.all.sector.index()%in% names(group2[group2=="pathways"])][sum(get.all.sector.index()%in%names(group2[group2=="pathways"]))]
    )
    add.do(x= data$ID ,h = data$Count,scale.col = scale_color,
           col= -log(data$pvalue), begin = (begin-1)*scale_height ,scale_height =scale_height*0.95 )%>%return()
    
  }
  
  f%>%dplyr::group_by(Cluster) %>% dplyr::group_split()%>%setNames(unique(f$Cluster))%>%
    purrr::map2(seq_along(unique(f$Cluster)),~grt(data =  .x,begin = .y))%>%.[[1]] -> length_pathways
  get.vip<-function(.x,.y) {
    tryCatch(
      .x$value%>%as.data.frame()->.x,
      error = function(e){
        return(NULL)
      },
      finally = {
      }
      
    )
    if(.y!="comp"){
      rownames(.x)<-paste(.y,"_",rownames(.x),sep = "");return(.x) }
    return(NULL)
  }
  mixOmics::selectVar(n@data_mixOmics_analyst[[1]])%>%
    purrr::map2((names(mixOmics::selectVar(n@data_mixOmics_analyst[[1]]))),~get.vip(.x,.y))->c
  c%>%purrr::reduce(rbind)->data
  
  get_add.names <- function(.x,.y) {
    colnames(.x)%>%paste(.y,.,sep = "_")->colnames(.x)
    .x[colnames(.x)%in%get.all.sector.index()]->.x
    return(.x)
  }
  #n@data%>%purrr::map2(n@data%>%names,get_add.names)%>%purrr::reduce(cbind)->data
  #scale_height=0.6
  .get_theme <- function(x,f=f,scale_height=0.6,begin=0,data=data) {
    #x="a"
    scale_height<-scale_height
    
    f[f$Cluster==x,]->k
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                 y = c((begin)*scale_height+scale_height, (begin)*scale_height+scale_height),
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                 y = c( 0, 0),
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                 y = c(((begin)*scale_height+scale_height)*abs(min(data[rownames(data)%in%
                                                                          get.all.sector.index(),]))/
                         (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                            max(max(data[rownames(data)%in%get.all.sector.index(),]))),
                       ((begin)*scale_height+scale_height)*abs(min(data[rownames(data)%in%get.all.sector.index(),]))/
                         (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                            max(max(data[rownames(data)%in%get.all.sector.index(),])))),
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                 col = "grey",track.index = 1,lty=2
    )
    
    circos.text(-1,track.index = 1,y=c((begin)*scale_height+scale_height)*0.9,labels =
                  max(data[rownames(data)%in%get.all.sector.index(),])%>%round(.,2) ,
                sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
    )
    
    circos.text(-1,track.index = 1,y=c((begin)*scale_height+scale_height)*0,labels =
                  min(data[rownames(data)%in%get.all.sector.index(),])%>%round(.,2) ,
                sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
    )
    
    circos.text(-1,track.index = 1,y = c(((begin)*scale_height+scale_height)*max(data[rownames(data)%in%get.all.sector.index(),])/
                                           (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                                              max(max(data[rownames(data)%in%get.all.sector.index(),]))),
                                         ((begin)*scale_height+scale_height)*max(data[rownames(data)%in%get.all.sector.index(),])/
                                           (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                                              max(max(data[rownames(data)%in%get.all.sector.index(),])))),labels =0 ,
                sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
    )
    
    circos.text(-0.5,track.index = 1,y=-0.1,labels =
                  k$Cluster[1] ,
                sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
    )
  }
  purrr::map(f$Cluster%>%unique,~.get_theme(x=.x,f=f,scale_height=scale_height,begin=0,data=data))
  #f$Cluster%>%unique%>%.[1]->.x
  #scale_height=0.6
  add.do(x=rownames(data)[rownames(data)%in%get.all.sector.index()],h=data[rownames(data)%in%get.all.sector.index(),],
         begin = scale_height/2,
         scale_height =scale_height,#scale.col = scale.col,
         col =data[rownames(data)%in%get.all.sector.index(),]
  )->length_hub
  
  
  n@data_mean%>%names%>%purrr::map(function(.x){n@data_mean%>%.[[.x]]->x;
    colnames(x)%>%str_remove("X")->colnames(x)
    colnames(x)[-1]<-paste(.x,colnames(x)[-1],sep = "_");return(x%>%as.data.frame())})%>%
    purrr::reduce(cbind)->mean
  mean[,c("group",get.all.sector.index()[get.all.sector.index()%in%f$line])]->mean
  mean->mean2
  (mean[,-1]-min(mean[,-1]))/max(mean[,-1]-min(mean[,-1]))->mean[,-1]
  
  #chord(f[,c(7,1)],directional = 0,group=group2,annotationTrack = NULL,
  #  annotationTrackHeight = c(0.01, 0.01))
  
  for (i in seq_along(table(f$Cluster))) {
    f[f$Cluster%in%(names(table(f$Cluster))[i]),]->l
    
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                 y = c( scale_height,scale_height)+0.05,
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                 col = "grey",track.index = 1,lty=2
    )
    
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                 y = c( scale_height+max(mean[,-1])*scale_height/2,scale_height+max(mean[,-1])*scale_height/2)+0.05,
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                 y = c( scale_height+max(mean[,-1])*scale_height/2/2,scale_height+max(mean[,-1])*scale_height/2/2)+0.05,
                 sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                 col = "grey",track.index = 1,lty=2
    )
    circos.text(-0.5,track.index = 1,y=scale_height+max(mean[,-1])*scale_height/2,labels =round(max(mean2[,-1]),2) ,
                sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
    )
    circos.text(-0.5,track.index = 1,y=scale_height+0.05,labels =round(min(mean2[,-1]),2),
                sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
    )
    circos.text(-0.5,track.index = 1,y=scale_height+max(mean[,-1])*scale_height/2/2+0.05,labels =
                  round(max(mean2[,-1])/2+min(mean2[,-1])/2,2),
                sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
    )
    col=c(ggsci::pal_aaas(palette = "default",alpha = 0.9)(10))
    pch<-16:25
    for (t in seq_len(mean%>%nrow())){
      {
        circos.lines(x=seq_along(get.all.sector.index()[get.all.sector.index()%in%l$line])-0.5,
                     y=c(mean[t,get.all.sector.index()[get.all.sector.index()%in%l$line]]/2/2+0.6)%>%as.numeric()+0.05,
                     col=col[t],track.index = 1,pt.col = col[t],
                     lwd=3,cex =1.5,pch = pch[t],
                     type="o",
                     sector.index=get.all.sector.index()[get.all.sector.index()%in%l$line][1])
      }
    }
    
    
  }
  length_mean =Legend(labels = c( table(f$Cluster)%>%names()), title = "length_mean", type = "points", pch = pch[seq_along(table(f$Cluster)%>%names())],
                      legend_gp = gpar(col =(c(ggsci::pal_aaas(palette = "default",
                                                               alpha = 0.9)(10))[seq_along(table(f$Cluster)%>%names())] )))
  for (i in get.all.sector.index()) {
    if(i%in%f$line){add.name(begin =scale_height+max(mean[,-1])*scale_height/2/2+0.05,
                             name =f$geneID[f$line==i]%>%unique(),x = i)
      
    }else{add.name(begin =scale_height+max(mean[,-1])*scale_height/2/2+0.05,
                   name =i,x = i)}
    
    
  }
  return(list(length_mean,length_pathways,length_hub))
}









#' pathways_analy to pathways and analy
#'
#' @param data daat@kegg_analyst$enrichKEGG$d@result
#' @param col data col like
#' @return ggplot2 
#' @export
#'
#' @examples
#'
##' \dontrun{ 
##' library(tidyverse)
##' library(data.table)
##' library(DOSE)
##' kegg_pathway1(data=list(a=c("cid:784","cid:962","cid:1004","cid:5862","cid:5886"),b=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1"),c=c("3.4.22.15","3.4.22.69","3.6.4.12","3.6.4.13"),
##' d=c("6921","6923","8453","8883","9039","9978","79699")))->da
##' save(da,file="pathways.RData")
##'  load("pathways.RData")
##' comp_dotplot(da@kegg_analyst$enrichKEGG$a)
##' da->lllll
##' lllll->da
##'  da@kegg_analyst$enrichKEGG$a@result%>% dplyr::filter(org=="KEGG") ->y
##'  da@kegg_analyst$enrichKEGG$d@result->d23
##' plot_funmap(d23)
##' plot_funmap(d23,col=c("none"="red"))
##' # or
##' }


#' @importFrom  igraph E V
#' @import ggraph
#' @import dplyr

plot_funmap<-function(result_data,col="auto"){
  result_data->d23
  d23%>%dplyr::group_by(ID)%>%dplyr::summarise(geneID %>% str_split("/")%>%unlist())%>%setNames(c("ID", "Compound_ID" ))%>%
      ungroup()%>%as.data.table()->da
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
    g<-igraph::graph_from_data_frame(
      d=l[,1:2],
      directed=F)
  E(g)$weight<- sqrt(l$weight * 5) * 0.2
  V(g)$group <- d23$group[match(V(g)$name,d23$ID)]
  V(g)$size <- d23$import[match(V(g)$name,d23$ID)]
  
  bb <- graphlayouts::layout_as_backbone(g,keep=0.2,  )

  data.frame(bb$xy,group=d23$group[match(V(g)$name,d23$ID)],ID=d23$ID[match(V(g)$name,d23$ID)])%>%group_by(group)%>% group_split%>% lapply( function(x){
  if(nrow(x)==1){return(  data.frame( X1 =0.5, X2 =0.5,x[,4])
  )}
   data.frame(x=(x[,1]-min(x[,1]))/(max(x[,1])-min(x[,1])),y=(x[,2]),x[,4])
   } )%>% reduce (rbind)->data
  data[match(V(g)$name,data$ID),]->bb$xy
  E(g)$w <- l[,3]
  V(g)$size <- d23$import[match(V(g)$name,d23$ID)]
  V(g)$grp <- d23$group[match(V(g)$name,d23$ID)]

  l[,1:2] %>% apply( 1,function(.x){
     g1=d23$group[d23$ID==.x[1]]
     g2=d23$group[d23$ID==.x[2]]
  if(g1==g2){ return(g1) }
  return("none")
  } )->E(g)$g

mypal = c(ggsci::pal_npg("nrc", alpha = 0.8)(9),ggsci::pal_d3("category20c", alpha = 0.8)(20))
if (col=="auto"){
col<-c("darkgrey",mypal[seq_along(E(g)$g%>% unique()%>%unlist%>%.[.!="none"])])
names(col)<-c("none",E(g)$g%>% unique()%>%unlist%>%.[.!="none"] )
} else if ("none"%in% names(col)) {
col<-c(mypal[seq_along(E(g)$g%>% unique()%>%unlist%>%.[! .%in%names(col)])])
names(col)<-c(E(g)$g%>% unique()%>%unlist%>%.[! .%in%names(col)] )
} else {
col<-c("darkgrey",mypal[seq_along(E(g)$g%>% unique()%>%unlist%>%.[.!="none"]%>%.[! .%in%names(col)])])
names(col)<-c("none",E(g)$g%>% unique()%>%unlist%>%.[.!="none"]%>%.[! .%in%names(col)] )
}
o<-ggraph(g,layout="manual",x=bb$xy[,1] + table(d23$group)%>% seq_along()%>%
.[match(d23$group[match(V(g)$name,d23$ID)],table(d23$group)%>% names())],y=bb$xy[,2])+
  geom_edge_link0(aes(width=w,col=g),alpha=.5)+
  scale_color_manual(values = col)+
  scale_edge_color_manual(  values = col)+

  geom_node_point(aes(col=grp,size=size))+
  theme_graph()+
  theme(legend.position = "none")+
  ggforce::geom_mark_hull(aes(x, y, group = grp, fill = grp, label=grp),
  concavity = 4,
    #expand = unit(2, "mm"),
  alpha = 0.1
  )+
  scale_fill_manual(  values = col)+
  scale_edge_width(range = c(0.2, 0.8))# control size
 return(o)
}











#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @import circlize  
#' @import ComplexHeatmap
#' @import randomcoloR
#' @return NULL 
#' @export
#'
#' @examples
#'
##' \dontrun{
##' library(tidyverse)
##' library(data.table)
##' library(DOSE)
##' kegg_pathway1(data=list(a=c("cid:784","cid:962","cid:1004","cid:5862","cid:5886"),b=c("P0DTD3","Q9BYF1","Q9NRS4","Q9NYK1"),c=c("3.4.22.15","3.4.22.69","3.6.4.12","3.6.4.13"),
##' d=c("6921","6923","8453","8883","9039","9978","79699")))->da
##' save(da,file="pathways.RData")
##'  load("pathways.RData")
##' comp_dotplot(da@kegg_analyst$enrichKEGG$a)
##' da->lllll
##' lllll->da
##'  da@kegg_analyst$enrichKEGG$a@result%>% dplyr::filter(org=="KEGG") ->y
##' library(randomcoloR)
##' library(circlize)
##' library(ComplexHeatmap)
##' FC<-c("ada"="mm")
##' FC<-c(0.5,1.2,0.4,0.1,2)
##' names(FC)<-y$geneID%>%unlist%>% str_split("/")%>%unlist %>%unique
##' plot_circos(y,FC)
##' }


plot_circos<-function(result_data,FC){
result_data->y
FC<-data.frame(FC=c(FC),geneID =names(FC))
circle_size = unit(1, 'snpc')
circos.par(gap.degree = 0.5, start.degree = 90)
ko_color <- c(rep('#F7CC13',10), rep('#954572',10), rep('#0796E0',9), rep('green', 6)) #各二级分类的颜色和数目
library(tidyverse)
y%>%dplyr::filter(pvalue<0.05)->y
col_chr = data.frame(group = c(y$group %>% unique()),
                     col = c(distinctColorPalette(length(y$group %>% unique()))))
full_join(y,col_chr)->y2
color_assign <- c(y2$col)
circos.genomicInitialize(data.frame(y$ID,0,1), plotType = NULL, major.by = 1)
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = color_assign,
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

plot_data <- tibble(y$ID,0,log2(y$Count/max(y$Count)+1),`-log10p`=-log10(y$pvalue))
label_data <- y$ID
p_max <- round(max(plot_data$`-log10p`)) + 1
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))
color_assign <- colorRamp2(breaks = 1:p_max, col = colorsChoice(p_max))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    ylim = get.cell.meta.data('ycenter')
    xlim = 1/2
    sector.name = get.cell.meta.data('sector.index')
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

y%>%group_by(ID)%>%summarise(geneID = geneID%>%str_split("/")%>%unlist() )%>%ungroup->da
left_join(da,FC)->FCdata
FCdata%>%group_by(ID)%>%summarise(0,FCda = sum(FC>1) / n(),1,num=n(),up = sum(FC > 1),down = sum(FC <= 1))->da
plot_data_up<-data.frame(id=da$ID,0,da$FCda,1)
plot_data_down<-data.frame(id=da$ID,da$FCda,1,2)
data.frame(da$up,da$down,da$FCda,1-da$FCda,row.names = da$ID)->label_data
colnames(label_data) <- c('up', 'down', 'up.regulated', 'down.regulated')
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))
colnames(plot_data_up)<-c("id","sta","end","group")
colnames(plot_data_down)<-c("id","sta","end","group")

rbind(plot_data_down,plot_data_up)->plot_data
circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = round(label_data[get.cell.meta.data('sector.index'),3],2)
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = round(label_data[get.cell.meta.data('sector.index'),4],2)
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

col_chr = data.frame(group = c( y$group  %>% unique()),
                     col = c(distinctColorPalette(length(y$group %>% unique()))))
full_join(y,col_chr)->y2
color_assign <- c(y2$col)#各二级分类的名称和颜色
names(color_assign)<-y2$ID

plot_data <- y2[,c(1,2,5)]
plot_data$Description<-0
plot_data$RichFactor%>%purrr::map(~parse(text =.x)%>%eval)%>%unlist()->plot_data$B
plot_data$RichFactor<-1

colnames(plot_data) <- c('id', 'gene_num.min', 'gene_num.max', 'value')

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,
  panel.fun = function(region, value, ...) {
    print( get.cell.meta.data('sector.index')  )
    sector.name = get.cell.meta.data('sector.index')
   circos.genomicRect(region, value, col = color_assign[get.cell.meta.data('sector.index')], border = NA, ytop.column = 1, ybottom = 0, ...)
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)
    sector.name = round(value%>%unlist()%>%.[[1]],2)
    circos.text(x = 0.5,y = value%>%unlist()%>%.[[1]] ,sector.name, cex = 1, niceFacing = FALSE)
} )
#library( ComplexHeatmap )
category_legend <- Legend(
  labels = c(col_chr$group),#各二级分类的名称
  type = 'points', pch = NA, background = c(col_chr$col), #各二级分类的颜色
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'),
  type = 'points', pch = NA, background = c('red', 'blue'),
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0),
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8),
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')
lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
grid.draw(lgd_list_vertical)
}



#' plot_chor is plot metaProfiler data 
#'
#' @param f using dplyr to filter data
#' @param metaProfiler  metaProfiler data
#' @param pathways_scale_color  pathways bar color
#' @param directional https://jokergoo.github.io/circlize_book/book/
#' @param annotationTrack https://jokergoo.github.io/circlize_book/book/
#' @param line_scale_col  line color
#' @param ... https://jokergoo.github.io/circlize_book/book/
#'
#' @return a list of length ,you can using packLegend and draw to add length
#' @export
#'
#' @examples
##' \dontrun{
##' load(system.file("data", "Tes.Rdata",package = "pathways"))
##'  n@kegg_analyst$compareClusterResult@compareClusterResult%>%dplyr::filter(
##' Description%in%c((n@kegg_analyst$compareClusterResult@compareClusterResult%>%group_by(Description)%>%
                    ##'                     summarise(n=n()>1))%>%dplyr::filter(n==T)%>%.$Description)|qvalue     <0.05
##' )%>%as.data.frame()->f
##' circos.clear()
##' circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,start.degree=90)
##' plot_chor(f,metaProfiler = n)->v
##' #looking for ComplexHeatmap[https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html]
##' pd = packLegend(list =  v)
##' draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
##' }



plot_chor<-function (f = f, metaProfiler, scale_height = 0.6, begin = 0, 
    pathways_scale_color = structure(c(-log(0.05), -log(0.01)),
        names = c("red", "blue")), directional = 0, annotationTrack = NULL,
    line_scale_col = NULL, ...){
    f <- f[, c("ID", "geneID", "Cluster", "Count", "pvalue",
        "GeneRatio")] %>% dplyr::group_by(ID, Cluster, Count,
        pvalue, GeneRatio) %>% dplyr::summarise(geneID = geneID %>%
        str_split("/") %>% unlist())
    f$line <- f$Cluster %>% paste(f$geneID, sep = "_")
    b <- f[, c(7, 1)] %>% unlist() %>% c
    group2 <- structure(c(f[, c(2, 7)] %>% unique %>% .$Cluster,
        rep("pathways", (!b %in% (f[, c(2, 7)] %>% unique %>%
            .$line)) %>% sum)), names = c(f[, c(2, 7)] %>% unique %>%
        .$line, b[!b %in% (f[, c(2, 7)] %>% unique %>% .$line)]))
    length_chord <- chord(data = f[, c(7, 1)], directional = directional,
        group = group2, target.prop.height = mm_h(1), diffHeight = mm_h(1),
        preAllocateTracks = 1, annotationTrack = annotationTrack,
        col = line_scale_col, small.gap = 0, annotationTrackHeight = c(0.01,
            0.01))
    k <- f
    {
        n <- metaProfiler
        scale_height <- scale_height
        grt <- function(data, begin, scale_height = 0.4, scale_color = structure(c(-log(0.05),
            -log(0.01)), names = c("red", "blue"))) {
            circos.lines(c(0, length(get.all.sector.index()[get.all.sector.index() %in%
                names(group2[group2 == "pathways"])]) + 0.5),
                y = c((begin - 1) * scale_height, (begin - 1) *
                  scale_height), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  names(group2[group2 == "pathways"])][1], col = "grey",
                track.index = 1, lty = 2)
            circos.lines(c(0, length(get.all.sector.index()[get.all.sector.index() %in%
                names(group2[group2 == "pathways"])]) + 0.5),
                y = c((begin - 1) * scale_height, (begin - 1) *
                  scale_height) + 0.5 * scale_height * 0.95,
                sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  names(group2[group2 == "pathways"])][1], col = "grey",
                track.index = 1, lty = 2)
            circos.lines(c(0, length(get.all.sector.index()[get.all.sector.index() %in%
                names(group2[group2 == "pathways"])]) + 0.5),
                y = c((begin - 1) * scale_height, (begin - 1) *
                  scale_height) + 1 * scale_height * 0.95, sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  names(group2[group2 == "pathways"])][1], col = "grey",
                track.index = 1, lty = 2)
            circos.text(-0.5, track.index = 1, y = (begin - 1) *
                scale_height + scale_height/2, labels = data$Cluster[1],
                sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  names(group2[group2 == "pathways"])][1])
            circos.text(x = rep(1.5, 3), track.index = 1, y = c((begin -
                1) * scale_height + scale_height/2, (begin -
                1) * scale_height + scale_height * 0.8, (begin -
                1) * scale_height + 0.05), labels = c(max(data$Count)/2,
                max(data$Count), 0), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                names(group2[group2 == "pathways"])][sum(get.all.sector.index() %in%
                names(group2[group2 == "pathways"]))])
            add.do(x = data$ID, h = data$Count, scale.col = scale_color,
                col = -log(data$pvalue), begin = (begin - 1) *
                  scale_height, scale_height = scale_height *
                  0.95) %>% return()
        }
        length_pathways <- f %>% dplyr::group_by(Cluster) %>%
            dplyr::group_split() %>% setNames(unique(f$Cluster)) %>%
            purrr::map2(seq_along(unique(f$Cluster)), ~grt(data = .x,
                begin = .y)) %>% .[[1]]
        get.vip <- function(.x, .y) {
            tryCatch(.x <- .x$value %>% as.data.frame(), error = function(e) {
                return(NULL)
            }, finally = {
            })
            if (.y != "comp") {
                rownames(.x) <- paste(.y, "_", rownames(.x),
                  sep = "")
                return(.x)
            }
            return(NULL)
        }
        c <- mixOmics::selectVar(n@data_mixOmics_analyst[[1]]) %>%
            purrr::map2((names(mixOmics::selectVar(n@data_mixOmics_analyst[[1]]))),
                ~get.vip(.x, .y))
        data <- c %>% purrr::reduce(rbind)
        get_add.names <- function(.x, .y) {
            colnames(.x) <- colnames(.x) %>% paste(.y, ., sep = "_")
            .x <- .x[colnames(.x) %in% get.all.sector.index()]
            return(.x)
        }
        .get_theme <- function(x, f = f, scale_height = 0.6,
            begin = 0, data = data) {
            scale_height <- scale_height
            k <- f[f$Cluster == x, ]
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                (k$line)] %>% length() + 0.5), y = c((begin) *
                scale_height + scale_height, (begin) * scale_height +
                scale_height), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1], col = "grey", track.index = 1,
                lty = 2)
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                (k$line)] %>% length() + 0.5), y = c(0, 0), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1], col = "grey", track.index = 1,
                lty = 2)
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                (k$line)] %>% length() + 0.5), y = c(((begin) *
                scale_height + scale_height) * abs(min(data[rownames(data) %in%
                get.all.sector.index(), ]))/(abs(min(data[rownames(data) %in%
                get.all.sector.index(), ])) + max(max(data[rownames(data) %in%
                get.all.sector.index(), ]))), ((begin) * scale_height +
                scale_height) * abs(min(data[rownames(data) %in%
                get.all.sector.index(), ]))/(abs(min(data[rownames(data) %in%
                get.all.sector.index(), ])) + max(max(data[rownames(data) %in%
                get.all.sector.index(), ])))), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1], col = "grey", track.index = 1,
                lty = 2)
            circos.text(-1, track.index = 1, y = c((begin) *
                scale_height + scale_height) * 0.9, labels = max(data[rownames(data) %in%
                get.all.sector.index(), ]) %>% round(., 2), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1])
            circos.text(-1, track.index = 1, y = c((begin) *
                scale_height + scale_height) * 0, labels = min(data[rownames(data) %in%
                get.all.sector.index(), ]) %>% round(., 2), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1])
            circos.text(-1, track.index = 1, y = c(((begin) *
                scale_height + scale_height) * max(data[rownames(data) %in%
                get.all.sector.index(), ])/(abs(min(data[rownames(data) %in%
                get.all.sector.index(), ])) + max(max(data[rownames(data) %in%
                get.all.sector.index(), ]))), ((begin) * scale_height +
                scale_height) * max(data[rownames(data) %in%
                get.all.sector.index(), ])/(abs(min(data[rownames(data) %in%
                get.all.sector.index(), ])) + max(max(data[rownames(data) %in%
                get.all.sector.index(), ])))), labels = 0, sector.index = get.all.sector.index()[get.all.sector.index() %in%
                (k$line)][1])
            circos.text(-0.5, track.index = 1, y = -0.1, labels = k$Cluster[1],
                sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  (k$line)][1])
        }
       # .x<- f$f$Cluster[1]
        purrr::map(f$Cluster %>% unique, ~.get_theme(x = .x,
            f = f, scale_height = scale_height, begin = 0, data = data))
        length_hub <- add.do(x = rownames(data)[rownames(data) %in%
            get.all.sector.index()], h = data[rownames(data) %in%
            get.all.sector.index(), ], begin = scale_height/2,
            scale_height = scale_height, col = data[rownames(data) %in%
                get.all.sector.index(), ])
        mean <- n@data_mean %>% names %>% purrr::map(function(.x) {
            x <- n@data_mean %>% .[[.x]]
            colnames(x) <- colnames(x) %>% str_remove("X")
            colnames(x)[-1] <- paste(.x, colnames(x)[-1], sep = "_")
            return(x %>% as.data.frame())
        }) %>% purrr::reduce(cbind)
        colnames(mean)%>% stringr::str_replace("cid\\.","cid:")  ->colnames(mean)
        mean <- mean[, c("group", get.all.sector.index()[get.all.sector.index() %in%
            f$line])]
        mean2 <- mean
        mean[, -1] <- (mean[, -1] - min(mean[, -1]))/max(mean[,
            -1] - min(mean[, -1]))
        for (i in seq_along(table(f$Cluster))) {
            l <- f[f$Cluster %in% (names(table(f$Cluster))[i]),
                ]
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                l$line] %>% length() + 0.5), y = c(scale_height,
                scale_height) + 0.05, sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1], col = "grey", track.index = 1, lty = 2)
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                l$line] %>% length() + 0.5), y = c(scale_height +
                max(mean[, -1]) * scale_height/2, scale_height +
                max(mean[, -1]) * scale_height/2) + 0.05, sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1], col = "grey", track.index = 1, lty = 2)
            circos.lines(c(0, get.all.sector.index()[get.all.sector.index() %in%
                l$line] %>% length() + 0.5), y = c(scale_height +
                max(mean[, -1]) * scale_height/2/2, scale_height +
                max(mean[, -1]) * scale_height/2/2) + 0.05, sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1], col = "grey", track.index = 1, lty = 2)
            circos.text(-0.5, track.index = 1, y = scale_height +
                max(mean[, -1]) * scale_height/2, labels = round(max(mean2[,
                -1]), 2), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1])
            circos.text(-0.5, track.index = 1, y = scale_height +
                0.05, labels = round(min(mean2[, -1]), 2), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1])
            circos.text(-0.5, track.index = 1, y = scale_height +
                max(mean[, -1]) * scale_height/2/2 + 0.05, labels = round(max(mean2[,
                -1])/2 + min(mean2[, -1])/2, 2), sector.index = get.all.sector.index()[get.all.sector.index() %in%
                l$line][1])
            col = c((ggsci::pal_aaas(palette = "default", alpha = 0.9))(10))
            pch <- 16:25
            for (t in seq_len(mean %>% nrow())) {
                {
                  circos.lines(x = seq_along(get.all.sector.index()[get.all.sector.index() %in%
                    l$line]) - 0.5, y = c(mean[t, get.all.sector.index()[get.all.sector.index() %in%
                    l$line]]/2/2 + 0.6) %>% as.numeric() + 0.05,
                    col = col[t], track.index = 1, pt.col = col[t],
                    lwd = 3, cex = 1.5, pch = pch[t], type = "o",
                    sector.index = get.all.sector.index()[get.all.sector.index() %in%
                      l$line][1])
                }
            }
        }
        length_mean = Legend(labels = c(table(f$Cluster) %>%
            names()), title = "length_mean", type = "points",
            pch = pch[seq_along(table(f$Cluster) %>% names())],
            legend_gp = gpar(col = (c((ggsci::pal_aaas(palette = "default",
                alpha = 0.9))(10))[seq_along(table(f$Cluster) %>%
                names())])))
        for (i in get.all.sector.index()) {
            circos.lines(c(0, 0), y = c(-0.1, scale_height),
                sector.index = get.all.sector.index()[get.all.sector.index() %in%
                  (i)][1], col = "grey", track.index = 1, lty = 2)
        }
        add.name <- function(x = x, name = name, begin = 0.8,
            col = "black", sector.index = NA, font = 1, cex = 1.2,
            facing = "clockwise", adj = c(0, (0)), niceFacing = TRUE,
            ...) {
            x <- x %>% as.data.frame() %>% unlist() %>% c
            name <- name %>% as.data.frame() %>% unlist() %>%
                c
            for (i in seq_along(x)) {
                circos.text(x = get.cell.meta.data("xlim", sector.index = x,
                  track.index = 1)[2]/2, y = begin, name, col = col,
                  sector.index = x, font = font, cex = cex, facing = facing,
                  adj = adj, niceFacing = niceFacing)
            }
        }
        for (i in get.all.sector.index()) {
            if (i %in% f$line) {
                add.name(begin = scale_height + max(mean[, -1]) *
                  scale_height/2/2 + 0.05, name = f$geneID[f$line ==
                  i] %>% unique(), x = i)
            }
            else {
                add.name(begin = scale_height + max(mean[, -1]) *
                  scale_height/2/2 + 0.05, name = i, x = i)
            }
        }
        length_col <- list(length_mean, length_pathways, length_hub)
    }
    return(list(chord = length_chord, mean = length_col[[1]],
        pathways = length_col[[2]], hub = length_col[[3]]))
}
