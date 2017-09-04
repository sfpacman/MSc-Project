#######################################################
## Processing PBMC68k data following 10X workflow    ##
#######################################################
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)

.get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
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
### loading data 
pbmc68k <- readRDS('pbmc68k_data.rds')
gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols
m<- pbmc68k$all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)
## PBMC68K expression matrix normalized by UMI count- used for Danaher methods ## 
saveRDS(l,"nor_10X_expr.rds")
### Obtained top 1000 variable genes for PCA and k-mean used for Danaher methods ##
m_n<-l$m
df<-.get_variable_gene(m_n) 
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)
k_n_1000<-kmeans(pca_n_1000$pca,10,iter.max=150,algorithm="MacQueen")
tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
tdf_n_1000<-data.frame(tsne_n_1000$Y)
nor_exp <- list(exp=m_n_1000,k=k_n_1000,tsne= tdf_n_1000)
saveRDS(nor_exp,"nor_exp.rds")


