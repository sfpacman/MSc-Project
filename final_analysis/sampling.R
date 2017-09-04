library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)
library(data.table)

# Function (get_pure_mtx_indx) adpoted and modtified from:
# https://github.com/10XGenomics/single-cell-3prime-paper/blob/master/pbmc68k_analysis/main_process_pure_pbmc.R
get_pure_mtx_indx <- function(id,pure_type_pca) {
  known_ids <- c("CD19+ B","CD14+ Monocyte"
,"Dendritic","CD56+ NK","CD34+"
,"CD4+/CD25 T Reg"
,"CD4+/CD45RA+/CD25- Naive T"
, "CD4+/CD45RO+ Memory"
,"CD4+ T Helper2","CD8+/CD45RA+ Naive Cytotoxic"
,"CD8+ Cytotoxic T")
  if (id %in% known_ids) { k <- 1; sel <- c(1) } #default is to select all populations
  else { stop('Unknown cell type:',id,'\n') }
  if (id=="CD19+ B") { check_genes <- c("CD74","CD27")}
  if (id=="CD14+ Monocyte") { check_genes <- c("FTL"); k <- 5; sel <- c(2,4) } # select based on FTL
  if (id=="Dendritic") { check_genes <- c("CLEC9A","CD1C"); k <- 5; sel <- c(5) } #select based on CLEC9A
  if (id=="CD34+") { check_genes <- "CD34"; k <- 8; sel <- c(1,2,3,7) } # select based on CD34
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
    sample_index= sample(1:nrow(mtx), size)
    mtx <- mtx[sample_index,]
    rownames(mtx) <- lapply(sample_index,function(x){paste(id,x,sep=":")})
 return(mtx)
 }

 get_sample <- function(pure_all_mtx, pure_all_pca,id , size = 1000,get_sum_list = FALSE, default = TRUE){
 addone <- FALSE
 sum_list <- list()
 i <- 1
 while(i <= length(pure_all_mtx)){
  if(addone){type= id[i+1]}
  else{type = id[i]}
  cat(type)
  sel_indx <-get_pure_mtx_indx(type,pure_all_pca[[i]])
  cat(paste("index return",length(sel_indx),"\n"))
  sel_mtx <- pure_all_mtx[[i]][sel_indx,]
  cat(paste("Selecting:",dim(sel_mtx),"\n"))
  s_mtx <- sampling_mtx(sel_mtx,type,size)
  cat("Sampling \n")
  if(!exists("all_mtx")){
    all_mtx <- s_mtx
   }
   else{
  all_mtx <- rbind(all_mtx, s_mtx)}
  cat("Combining \n")
  if(get_sum_list){
  sum_list <- rbind(unlist(sum_list),cbind(rownames(s_mtx), rep(type, nrow(s_mtx))))
 }
  if(default & i == length(pure_all_mtx)) {
    addone=TRUE
    default = FALSE
  }
  else{ i <- i + 1} }
 
 result <-list()
 result$exp <- all_mtx
 result$summary <- sum_list 
 colnames(result$summary) <-c("id","result")
 return(result)
 }
 
 ######




