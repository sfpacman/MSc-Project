# ---
#   Reading Gene10X filtered for CD34 cell data 
# ---
### This is a trial run for combining two Gene10X data and clustering 
library(cellrangerRkit)
# Modtifited Functions
.compare_by_cor<-function(m_filt,use_gene_ids,dmap_data) {
  
  sig_genes <- intersect(use_gene_ids, rownames(dmap_data))
  #expression matrix of gene_id X cell(barcode) 
  
  m_forsig <- as.matrix(m_filt[which(use_gene_ids %in% sig_genes),])
  #expression matrix gene_id X pure_11 pop(barcode)
  # msig <- exp_gbm1[which(use_genes_n_ens %in% sig_genes),]
  #sig_data_filt <- ppt[which(rownames(ppt) %in% sig_genes),]
  sig_data_filt <- dmap_data[match(use_gene_ids[which(use_gene_ids %in% sig_genes)], rownames(dmap_data)),]
  #co
  z <- lapply(1:ncol(sig_data_filt), function(j) sapply(1:ncol(m_forsig), function(i) cor(m_forsig[,i], sig_data_filt[,j], method='spearman')))
  z <- do.call(cbind, z)
  colnames(z) <- colnames(sig_data_filt)
  z
}
.downsample_gene_bc_mtx <- function(json, orig_matrix_data, mol_info, tgt_rpc, target_type='raw_reads', transcriptome='hg19') {
  if (!(target_type %in% c('raw_reads', 'conf_mapped_reads'))) {
    stop(sprintf('Unsupported target_type: %s', target_type))
  }
  if (length(orig_matrix_data) > 2) {
    warning('Multiple transcriptomes not yet implemented')
  }
  cat("Filtering mol info\n")
  mol_info <- mol_info[reads > 0]
  cat("Sorting mol info\n")
  setkey(mol_info, barcode, gene, umi)
  cat("Aggregating mol info\n")
  bc_gene_umi <- mol_info[, j=list(reads=sum(reads)), by=c('barcode', 'gene', 'umi')]
  n_cells <- json[[sprintf("%s_filtered_bcs", transcriptome)]]
  tot_reads <- json$total_reads
  candidate_reads <- sum(bc_gene_umi$reads)
  candidate_read_frac <- candidate_reads / json$total_reads
  orig_mat_barcodes <- sub('-.*$', '', orig_matrix_data[[transcriptome]]$barcodes)
  subsampled_mats <- lapply(tgt_rpc, function(tgt_rpc_i) {
    cat('.')
    if (target_type == 'raw_reads') {
      tgt_candidate_reads <- tgt_rpc_i * n_cells * candidate_read_frac
      candidate_rpc <- tgt_candidate_reads / n_cells
      raw_reads_per_cell <- tgt_rpc_i
    } else if (target_type == 'conf_mapped_reads') {
      tgt_candidate_reads <- tgt_rpc_i * n_cells
      candidate_rpc <- tgt_candidate_reads / n_cells
      raw_reads_per_cell <- candidate_rpc / candidate_read_frac
    }
    if (tgt_candidate_reads > candidate_reads) {
      return(NA)
    }
    subsample_rate <- tgt_candidate_reads / candidate_reads
    cat("Subsampling\n")
    bc_gene_umi_subsampled <- bc_gene_umi %>% mutate(reads=rbinom(length(reads), reads, subsample_rate))
    cat("Sorting\n")
    setkey(bc_gene_umi_subsampled, barcode, gene)
    cat("Re-aggregating\n")
    bc_gene_counts <- bc_gene_umi_subsampled[barcode %in% orig_mat_barcodes, j=list(count=sum(reads > 0)), by=c('barcode', 'gene')]
    cat("Building matrix\n")
    with(bc_gene_counts, sparseMatrix(i = match(barcode, orig_mat_barcodes),
                                      j = 1 + gene, x = count, dims=dim(orig_matrix_data[[transcriptome]]$mat)))
  } )
  return(subsampled_mats)
}
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

.normalize_by_umi <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=x_use_genes)
}
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
#####
gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex')
analysis_results <- load_cellranger_analysis_results("~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex")
pure_11 <- readRDS('~/bioinfo/Project/labeling/data/all_pure_select_11types.rds')
tdf<- analysis_results$tsne
gen<- fData(gbm1)
#getting gene id
pure_11$pure_use_gene_name <-sapply(pure_11$pure_use_genes,function(x){gen[x,]$id})
#need to invert this for cellRangaer function 
#Normalize the umi count
use_genes <- get_nonzero_genes(gbm1)

#gbm_bcnorm <- normalize_barcode_sums_to_median(gbm1[use_genes,tdf$Barcode])
#Error: cannot allocate vector of size 10.4 Gb
gbm_log <- log_gene_bc_matrix(gbm1,base=10)
gbm_e_log<- exprs(gbm_log)
#get variable gene 
df<-.get_variable_gene(t(gbm_e_log))
## This is used for PCA and K-mean steps 
#disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
#df$used<-df$dispersion_norm >= disp_cut_off
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_ens <- unlist(gbm_e_log@Dimnames[1])[use_genes_n]
z_1000_11<-get_correletion_scores(gbm_log,t(ppt),use_genes_n_ens[1:1000])
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])
tdf$cls_id<-cls_id
#note to self: need to find out what line 29 on util.R file 
