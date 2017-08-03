library(cellrangerRkit)
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)
library(data.table)

# Functions
get_pure_mtx_indx <- function(id,pure_type_pca, sav_rds=FALSE) {  
  known_ids <- c("CD19+ B","CD14+ Monocyte","Dendritic","CD56+ NK","CD34+","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",                 "CD4+/CD45RO+ Memory","CD4+ T Helper2","CD8+/CD45RA+ Naive Cytotoxic","CD8+ Cytotoxic T") 
  if (id %in% known_ids)  {    k <- 1; sel <- c(1) # default is to select all populations  } 
  else {    stop('Unknown cell type:',id,'\n')  }  
  if (id=="CD19+ B") { check_genes <- c("CD74","CD27")}   
  if (id=="CD14+ Monocyte") { check_genes <- c("FTL");      k <- 5; sel <- c(2,4) }  # select based on FTL  
  if (id=="Dendritic") { check_genes <- c("CLEC9A","CD1C"); k <- 5; sel <- c(5) }      # select based on CLEC9A  
  if (id=="CD34+")  { check_genes <- "CD34";                k <- 8; sel <- c(1,2,3,7) }# select based on CD34  
  if (id=="CD56+ NK") {check_genes <- c("CD3D","NKG7") }  
  if (id=="CD4+/CD25 T Reg") { check_genes <- "NKG7" }  
  if (id=="CD4+/CD45RA+/CD25- Naive T") { check_genes <- "NKG7" }  
  if (id=="CD4+/CD45RO+ Memory") { check_genes <- "NKG7" }  
  if (id=="CD4+ T Helper2") { check_genes <- "NKG7" }  
  if (id=="CD8+/CD45RA+ Naive Cytotoxic") { check_genes <- "NKG7"}   
  if (id=="CD8+ Cytotoxic T") { check_genes <- c("CD8A","CD3D")}   
  k_means<-kmeans(pure_type_pca$pca,k)  
  selection <- k_means$cluster %in% sel 
  # return logical object for subsampled_purified_mats selection (dgCMatrix) cell  
  # such as as.matrix(subsampled_purified_mats[[1]][selection,])  
  return(selection)}
  
 sampling_mtx <- function ( mtx, id ,size ){
    sample_index= sample(1:nrow(x), size)
    mtx <- mtx[sample_index,]
    rownames(mtx) <- lapply(sample_index,function(x){paste(id,x,sep="")}) 
 return(mtx) 
 }
 
 get_sample <- function(pure_all_mtx, pure_all_pca,id , size = 1000,sum_list = FALSE, default = TRUE){
 all_mtx <- list()
 sum_list <- list()
 one_more <- FALSE
 for(i in length(pure_all_mtx){
 type = id[i]
  sel_indx <-get_pure_mtx_indx(pure_all_pca[i],type)
  sel_mtx <- pure_all_mtx[[i]][sel_idx,] 
  s_mtx <- sampling_mtx(sel_mtx,size)
  all_mtx <- c(unlist(all_mtx), s_mtx)
  if(sum_list){
  sum_list <- rbind(sum_list,cbind(rownames(s_mtx), rep(type, length(s_mtx))))
 }
 if(i == 10 & default) {
  i <- i - 1 
  id[10] <- id[11]
 }
 }
 result$exp <- sum_list
 result$summary <- sum_list 
 return(result)
 }
 
 
