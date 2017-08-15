library(cellrangerRkit)
library(magrittr)
library(dplyr)
library(reshape2)
library(Matrix)
### for rds file
load_purified_pbmc_types<-function(pure_select_file) {
  pure_select_id <- pure_select_file$pure_id   # from pure_select_file
  pure_select_avg <- pure_select_file$pure_avg # from pure_select_file
  pure_use_genes <- pure_select_file$pure_use_genes # from pure_select_file
  pure_use_genes_ens<-pure_select_file$pure_use_gene_name
  avg<-data.frame(t(pure_select_avg))
  rownames(avg)<-pure_select_file$pure_use_gene_name
  names(avg)<-pure_select_id
  return(avg)
}

.do_propack <- function(x,n) {
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  m <- sweep(m, 2, apply(m, 2, sd), '/')
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(ppk=ppk,pca=pca, m=m,use_genes=use_genes)
}

.compare_by_cor<-function(m_filt,use_gene_ids,dmap_data) {
  
  sig_genes <- intersect(use_gene_ids, rownames(dmap_data))
  m_forsig <- as.matrix(m_filt[,which(use_gene_ids %in% sig_genes)])
  sig_data_filt <- dmap_data[match(use_gene_ids[which(use_gene_ids %in% sig_genes)], rownames(dmap_data)),]
  
  z <- lapply(1:ncol(sig_data_filt), function(j) sapply(1:nrow(m_forsig), function(i) cor(m_forsig[i,], sig_data_filt[,j], method='spearman')))
  z <- do.call(cbind, z)
  colnames(z) <- colnames(sig_data_filt)
  z
}
                                                        
.reassign_pbmc_11<-function(z) {
  unlist(lapply(1:nrow(z),function(i) {
    best<-which.max(z[i,])
    x<-z[i,]
    nextbest<-which.max(x[x!=max(x)])
    # if best is CD4+ T helper, and the next best is cd4+/cd25+, or cd4+/cd45ro+ or cd4+/cd45ra+/cd25-, use the next best assignment
    if (best==9 & (nextbest==3 || nextbest==4 || nextbest==6)) {
      best=nextbest
    }
    # if best is CD8+, and the next best is CD8+/CD45RA+, use next best assignment
    if (best==7 & nextbest==5) {
      best=5
    }
    best
  }))
} 
                                                        
get_cor_assign(mat, ref_mat, use_genes_n_ens){
m <-.compare_by_cor(mat,use_genes_n_ens[1:1000],ref_mat)
test<-.reassign_pbmc_11(m)
cls_id<-factor(colnames(m)[test])
return(cls_id) 
}
                                                        
get_danaher_assign(mat,ngn,sig_def, cluster,tsne){
for(i in 1:length(sig_def)){
 type= cell_list[i]
 type_gene_select <- match(unlist(cell_list[i]),ngn)
 type_gene_select <- type_gene_select[!is.na(type_gene_select)]
 type_expr <- mat[,type_gene_select]
 type_mean <-colMeans(t(type_expr))
 mean_score <- cbind(unlist(mean_score), type_mean)
 }
colnames(mean_score) <- names(cell_list)
cluster_mean <- data.frame(cbind(cluster,mean_score))

score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
cluster_type <-apply(score_by_cluster[,-which(names(cluster_mean) %in% c("Cluster"))], 1,function(x) which(x == max(x)))
cluster_type <-data.frame(cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_tsne <-merge( x=data.frame(cbind(cluster,tsne)) , y=cluster_type, by.x="cluster", by.y="cluster")
return(cluster_tsne)
}


### For cell ranger object 
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

load_purified_pbmc_types<-function(pure_select_file,pure_gene_ens) {
  pure_select_id <- pure_select_file$pure_id   # from pure_select_file
  pure_select_avg <- pure_select_file$pure_avg # from pure_select_file
  pure_use_genes <- pure_select_file$pure_use_genes # from pure_select_file
  pure_use_genes_ens<-pure_gene_ens[pure_use_genes]
  avg<-data.frame(t(pure_select_avg))
  rownames(avg)<-pure_use_genes_ens
  names(avg)<-pure_select_id
  return(avg)
}
pure_11 <- readRDS('~/bioinfo/Project/labeling/data/all_pure_select_11types.rds')
pure_pbmcs <- readRDS('~/bioinfo/Project/labeling/data/all_pure_pbmc_data.rds')


