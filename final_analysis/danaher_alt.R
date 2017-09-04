### This is alternative code of Danaher's scoring for processing data from pbmc68K.rds and used to produce scores and corresponding graph in the thesis
### some of the analysis and projection(tsne,kmean,dispersion rank) require data from 10X workflow which is stored in nor_exp.rds and nor_10X_expr.rds
library(Matrix)
library(magrittr)
library(dplyr)
library(reshape2)
library(plyr)
library(cellrangerRkit)
source("danaher_def.R")
### function for specific genes visualization ### 
visualize_me<- function(gbm,gene_probes,projection,limits=c(0,10),marker_size=0.1,title=NULL) {
  moreone=FALSE
  gene_values <-gbm
  if(ncol(gene_values) >1){
    colnames(gene_values) <- gene_probes
    moreone = TRUE}
  projection_names <-  colnames(projection)
  colnames(projection) <- c('Component.1', 'Component.2')
  proj_gene <- data.frame(cbind(projection,gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars=c("Component.1", "Component.2"))
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour=value),size=marker_size) 
  if(moreone){p <-p+facet_wrap(~variable)}
  p<- p + scale_colour_gradient(low="grey",high="red",name = "val")  +
    labs(x=projection_names[1],y=projection_names[2]) +
    if (!is.null(title)) {  p <- p + ggtitle(title) }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(p)
}
nor_10X <-readRDS("nor_10X_expr.rds")
nor_exp <-readRDS("nor_exp.rds")
pbmc68k<-readRDS("pbmc68k_data.rds")
gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols
ngn <- sapply(nor_10X$use_gene, function(x){gen[x]})
mean_score <- list()
### mean score calcluation  
for(i in 1:length(cell_list)){
 type= cell_list[i]
 type_gene_select <- match(unlist(cell_list[i]),ngn)
 type_gene_select <- type_gene_select[!is.na(type_gene_select)]

 type_expr <-nor_10X$m[,type_gene_select]
 type_mean <-colMeans(t(type_expr))
 mean_score <- cbind(unlist(mean_score), type_mean)
}
colnames(mean_score) <- names(cell_list)
cluster_mean <- data.frame(Cluster=nor_exp$k$cluster,mean_score)

score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)

cluster_type <-apply(score_by_cluster[,-which(names(cluster_mean) %in% c("Cluster"))], 1,function(x) which(x == max(x)))

cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_tsne <-merge( x=data.frame(cbind(cluster=nor_exp$k$cluster,nor_exp$tsne)) , y=cluster_type, by.x="cluster", by.y="Cluster")

# mean expression for top 1000 vaiable genes 
expression_cluster <-data.frame( as.matrix(cbind(exp=rowMeans(nor_exp$exp),clu=nor_exp$k$cluster)))
reference_expression <- expression_cluster %>% group_by(clu) %>% summarise_all(funs(mean))
score_diff_p<- sapply(1:length(cluster_type), function(x){t.test(mean_score[which(nor_exp$k$cluster==x),cluster_type[x]],mu=reference_expression[[x,2]])$p.value})
all_gene_var_pos <- list()
## danaher marker dispersion: disperion rank data comes from 10X analysis 
for(i in 1:length(cell_list)){
type_gene_select <- match(unlist(cell_list[i]),nor_exp$gene$use_genes_n_id)
all_gene_var_pos[[length(all_gene_var_pos)+1]] <- type_gene_select
}
#This produces raw data for Table 2 
names(all_gene_var_pos) <- names(cell_list)
for(x in 1:length(cell_list)){names(all_gene_var_pos[[x]]) <- cell_list[[x]]}

### cluster labeling on tSNE plot ###                      
tsne_center <- list()
for(i in 1:length(unique(cluster_tsne$cluster))){
c.x<-mean(cluster_tsne$X1[cluster_tsne$cluster==i])
c.y<-mean(cluster_tsne$X2[cluster_tsne$cluster==i])
tsne_center <- rbind(unlist(tsne_center), c(c.x,c.y))
}
## Viusulaizing CD4 and CD8 genes -produce fig 1A and 1B ##
for(cell_type in c("CD4",c("CD8A","CD8B"))){
  type_gene_select <- match(cell_type,ngn)
  #this is for scale adjustment for the tSNE projection 
  if(length(cell_type) >1) {x <- log2(as.numeric((summary(nor_10X$m[which(nor_10X$m[,type_gene_select[1]]>0),type_gene_select[1]])))+1)}
  else{x <- log2(as.numeric(summary(nor_10X$m[which(nor_10X$m[,type_gene_select]>0),type_gene_select]))+1)}
  gens <- as.matrix(log2(nor_10X$m[,type_gene_select]+1))
  visualize_me(gens,cell_type, nor_exp$tsne[c("X1","X2")])+scale_colour_gradientn(colours = c("#e6e6e6","#ff3333","#e60000","#e60000","#4d0000","#4d0000"),values = rescale(c(0,x[2:5])),guide = "colorbar", limits=c(0, x[6])) + labs(x="TSNE.1",y="TSNE.2")
}
  
## This produces Fig 2a ### 
visualize_clusters(cluster_tsne$cluster,cluster_tsne[c("X1","X2")])+guides(col = guide_legend(title="Cluster", override.aes = list(size=3)))+annotate("text", x = tsne_center[, 1], y = tsne_center[, 2], size = 5,label=c(1:length(tsne_center[,1])))+labs( x= "TSNE1" , y="TSNE2")
### This produces Fig 2b ###
visualize_clusters(cluster_tsne$cell_type,cluster_tsne[c("X1","X2")],title="Danaher cell-type labels",legend_anno= sort(unique(cluster_tsne[,"name_type"])))+scale_color_manual(values = colorRampPalette(c("blue", "red","yellow"))( length(sort(unique(cluster_tsne[,"name_type"])))))+  annotate("text", x = tsne_center[, 1], y = tsne_center[, 2], size = 5,label=c(1:length(tsne_center[,1])))

### This produces Supp. Fig 1 ### 
sapply(1:length(unlist(score_by_cluster[,1])),function(x) { 
  plot_x <- plot(unlist(score_by_cluster[x,-1]),ylab="Cluster Cell Score",xaxt='n', xlab="Cell Type",xlim=c(0,length(score_by_cluster[1,]+1)), main=c("Cluster",x),ylim=c(-0.1,max(score_by_cluster[x,-1],reference_expression[x,2])*1.15))
  plot_x <- text(unlist(score_by_cluster[x,-1]), labels=(names(cell_list)), cex=0.5, pos=4,offset = 1 ,col="red")
  abline(h=reference_expression[x,2],col="green",lty = 2) 
  
})
