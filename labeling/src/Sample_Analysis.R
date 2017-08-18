library(magrittr)
library(dplyr)
library(reshape2)
library(Matrix)
library(cellrangerRkit)
source("danaher_def.R")
source("danaher_function.R")
sample <- readRDS("sample.rds")

## modtified:not enough unquie breaking points- add unique functions on quantile on cut 
.get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  #df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,unique(quantile(mean,seq(0.1,1,0.05))),Inf)))
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
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
.normalize_by_umi <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=x_use_genes)
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
sample <- readRDS("~/bioinfo/Project/labeling/data/sample.rds")
pbmc68k <- readRDS("~/bioinfo/Project/labeling/data/pbmc68k.rds")
### UMI Normalization ###
m <-sample$exp 
l<-.normalize_by_umi(m)
m_n<-l$m
sample_analysis$X <- list(m_n= l$m, use_genes= l$use_genes)
### computesumfactor nomrmalization###
use_genes<- which(colMeans(m) > 0)
M <- m[,use_genes]
M <- t(as.matrix(M))
high.ab <- which(rowMeans(M) >1)
clusters <- quickCluster(M, subset.row=high.ab)
nor_fac <- computeSumFactors(M,subset.row=high.ab,cluster=clusters)
M_n <- nor_fac*t(M)
sample_analysis$S <- list(m_n= M_n, use_genes= use_genes)

### Corelation ###
df<-.get_variable_gene(m_n) 
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
use_genes_n<-order(-df$dispersion_norm)
#use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-pbmc68k$all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11) 
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])

### Danaher ###
pca_n_1000<-.do_propack(m_n_1000,50)
k_n_1000<-kmeans(pca_n_1000$pca,10,iter.max=150,algorithm="MacQueen")
gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols

mean_score <- list()
for(i in 1:length(cell_list)){
 type= cell_list[i]
 type_gene_select <- match(unlist(cell_list[i]),ngn)
 type_gene_select <- type_gene_select[!is.na(type_gene_select)]
 type_expr <- m_n[,type_gene_select]
 type_mean <-colMeans(t(type_expr))
 mean_score <- cbind(unlist(mean_score), type_mean)
 }

colnames(mean_score) <- names(cell_list)
cluster_mean <- data.frame(cbind(Cluster=k_n_1000$cluster,mean_score))

score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
score_by_cluster <- score_by_cluster[,-which(names(score_by_cluster) %in% c("Cluster"))]NA
cluster_type <- list()
for( i in 1:nrow(score_by_cluster)){
  x <- as.numeric(score_by_cluster[i,])
  if( mean(x) == 0 ) {cluster_type = c(unlist(cluster_type),0) }
  else {cluster_type = c(unlist(cluster_type), which(x== max(x)))} } 

rcell_assign <-sapply(cluster_mean$Cluster,function(x){cluster_type[x]})
act <- unlist(sample$summary[,2])
result <- cbind(pre=rcell_assign,act =act)