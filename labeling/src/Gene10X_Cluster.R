# ---
#   Reading Gene10X filtered for CD34 cell data 
# ---
### This is a trial run for combining two Gene10X data and clustering 
library(cellrangerRkit)
# Modtifited Functions
sparse_pca_c <- function(x, n_pcs, mu=NULL, s=NULL, center_scale=TRUE) {
  if (is.null(mu) && center_scale) mu <- colMeans(x)
  if (is.null(s) && center_scale) s <- apply(x, 2, sd, na.rm=TRUE)
  
  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba::irlba(x, n_pcs, center=mu, scale=s,fastpath=FALSE)
  } else {
    svd_res <- irlba::irlba(x, n_pcs,fastpath=FALSE)
  }
  
  # compute explained variance
  n <- dim(x)[1]
  variance_sum <- sum(apply(x,2,var,na.rm=TRUE)/(s^2)) # sample variance sum
  var_pcs <- svd_res$d^2/(n-1)/variance_sum
  
  return(list(x=svd_res$u %*% diag(svd_res$d), rotation=svd_res$v, sdev=svd_res$d/sqrt(n-1),
              tot_var=variance_sum, var_pcs=var_pcs))
}

run_pca_c <- function(gbm, n_pcs=10, logscale=FALSE, ...) {
  if (logscale) {
    gbm_log <- gbm
  } else {
    use_genes <- get_nonzero_genes(gbm)
    gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
    gbm_log <- log_gene_bc_matrix(gbm_bcnorm, ...)
  }
  
  pca <- sparse_pca_c(t(exprs(gbm_log)), n_pcs)
  
  pc_names <- sprintf("PC%d", 1:ncol(pca$x))
  rownames(pca$x) <- colnames(gbm)
  colnames(pca$x) <- pc_names
  
  rownames(pca$rotation) <- rownames(gbm)[use_genes]
  colnames(pca$rotation) <- pc_names
  var_explained <- sum(pca$var_pcs)
  cat('Variance explained by PCs:',var_explained,'\n')
  return(list(x=pca$x,
              rotation=pca$rotation,
              sdev=pca$sdev,
              tot_var=pca$tot_var,
              var_pcs=pca$var_pcs,
              use_genes=use_genes,
              normalized_mat=exprs(gbm_log),
              var_explained=var_explained))
}

#####
gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data')
gbm2 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/')
#c= exprs(gene_bc_matrix)
#dim(c)
# has 32738 cell profile among 9232 gene
# Getting supplied analysis result 
#analysis_results <- load_cellranger_analysis_results("~/bioinfo/Project/labeling/data")

 use_genes <- get_nonzero_genes(merged_gbm)
 gbm_bcnorm <- normalize_barcode_sums_to_median(merged_gbm[use_genes,])
 gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
 set.seed(0)
 gbm_list <- list(gbm1, gbm2)
 gbm_list <- lapply(gbm_list,load_molecule_info) # load sample molecule information
 gbm_list_equalized <- equalize_gbms(gbm_list) # equalize the gene-barcode matrices
 merged_gbm <- concatenate_gene_bc_matrices(gbm_list_equalized)
 merged_ID <- unlist(lapply(1:length(gbm_list), function(x) rep(x,dim(gbm_list[[x]])[2])))
 set.seed(0)
 n_clust <- 5
 pca_result <- run_pca_c(merged_gbm)
 tsne_result <- run_tsne(pca_result)
 kmeans_result <- run_kmeans_clustering(pca_result,k=n_clust)
 # included re-named cell IDs, t-SNE and k-means result in a merged data frame
 merged_tsne_clust <- data.frame(Barcode=as.factor(1:tsne_result$N),
                                TSNE.1=tsne_result$Y[,1],TSNE.2=tsne_result$Y[,2])
# dim(gbm_log) use 19K genes 
# 19259  9232
#tsne_proj <- analysis_results$tsne
gen[match('CD34',gen$symbol),]$id #fData(gbm1) 
colnames(merged_tsne_clust) <- c('TSNE.1', 'TSNE.2')
visualize_gene_markers(gbm_log,"CD34",merged_tsne_clust[c("TSNE.1","TSNE.2")],limits=c(0,1.5))

#DATA_DIR ="~/bioinfo/Project/labeling/data"
#pure_11 = readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))
#n = length(pure_11[[1]])
#DF = structure(pure_11, row.names = c(NA, -n), class = "data.frame")
# purl_all= readRDS(file.path(DATA_DIR,'all_pure_pbmc_data.rds'))
#Test by corrleation 
# Get data from  pure_11
#for example, such as CD34
#CD34 <- pure_11$pure_avg[1,]
# select genes that are found in donor population
#CD34_select = CD34[,which(dummy[,2] %in% CD34[,2])]
# find the highest correlaiton 

#CD34+ <- pure_11$pure_avg[1,]
#note to self: need to find out what line 29 on util.R file 