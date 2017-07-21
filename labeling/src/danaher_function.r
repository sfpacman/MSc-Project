library(cellrangerRkit)
library(magrittr)
library(dplyr)
library(reshape2)
library(Matrix)

get_sig_profile <- function(gbm,analysis_results,sig_def){
  cluster_result <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]
  tsne_proj <- analysis_results$tsne
  gbm <- gbm[,analysis_results$tsne$Barcode]
  gen= fData(gbm1)
  ##Selecting Genes marker from the expression data 
  mean_score <-list()
  type_score_table <-list()
  all_type_expr <-list()
  var_score <- list()
  for(i in 1:length(cell_list)){
    type= cell_list[i]
    type_gene_select <- gen[match(unlist(type),gen$symbol),]$id
    type_gene_select <-type_gene_select[!is.na(type_gene_select)]
    type_expr <- exprs(gbm1[type_gene_select,])
    ## computing mean score and variance  for each cells 
    ## Paper use log2 Transformation- subjected to change 
    #type_mean = log2(colMeans(type_expr))
    type_mean <- colMeans(type_expr)
    var_mean <- sd(type_expr)
    mean_score <- cbind(unlist(mean_score), type_mean)
    var_score <- cbind(unlist(var_score), type_mean)
  }
  colnames(mean_score) <- names(cell_list)
  colnames(var_score) <- names(cell_list)
  cluster_mean <- merge(x=cluster_result,y=mean_score,by.y="row.names" ,by.x="Barcode")
  cluster_var <-merge(x=cluster_result,y=var_score^2,by.y="row.names" ,by.x="Barcode")
  score_by_cluster <- round(cluster_mean[,-which(names(cluster_mean) %in% c("Barcode"))] %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
  score_by_cluster_sum <- round(cluster_mean [,-which(names(cluster_mean) %in% c("Barcode"))]%>% group_by(Cluster) %>% summarise_all(funs(sum)),3)
  score_by_cluster_var <- cluster_var[,-which(names(cluster_var) %in% c("Barcode"))]%>% group_by(Cluster) %>% summarise_all(funs(sum))
  score_by_cluster_var <- round(sqrt(score_by_cluster_var[,-which(names(score_by_cluster_var) %in% c("Cluster"))]),3)
  score_by_cluster_var <- cbind(Cluster=rownames(score_by_cluster_var),score_by_cluster_var)
  all_type_expr_table <- merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster,"cell_mean_expression"), by=c("Cluster","Signature") )
  all_type_expr_table <-  merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster_sum,"cell_mean_expression_sum"), by=c("Cluster","Signature") )
  all_type_expr_table <-  merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster_var,"cell_mean_expression_sd"), by=c("Cluster","Signature") )
  score_by_cluster <- score_by_cluster[, -which(names(score_by_cluster) %in% c("Cluster"))]
  cluster_type <-apply(score_by_cluster, 1,function(x) which(x == max(x)))
  cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
  cluster_assignment <- merge(x=cluster_result,y=cluster_type,by="Cluster")
  cluster_tsne <-merge(x=cluster_assignment, y=tsne_proj,by="Barcode")
  
}