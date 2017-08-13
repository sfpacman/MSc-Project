### This is alternative code for processing data from pbmc64k.rds 

gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols
ngn <- sapply(nor_10x$use_gene, function(x){gen[x]})
for(i in 1:length(cell_list)){
   type= cell_list[i]
 type_gene_select <- match(unlist(cell_list[i]),ngn)
 type_gene_select <- type_gene_select[!is.na(type_gene_select)]
 type_expr <-nor_10x$m[,type_gene_select]
 type_mean <-colMeans(t(type_expr))
 mean_score <- cbind(unlist(mean_score), type_mean)
 }
colnames(mean_score) <- names(cell_list)
cluster_mean <- data.frame(cbind(Cluster=nor_exp$k$cluster,mean_score))
score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
cluster_type <-apply(score_by_cluster[,-which(names(cluster_mean) %in% c("Cluster"))], 1,function(x) which(x == max(x)))
cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_assignment <-merge( x=data.frame(cbind(cluster=nor_exp$k$cluster,nor_exp$tsne) , y=cluster_type, by.x="k", by.y="Cluster")
sig_plot <- visualize_clusters(cluster_tsne$cell_type,cluster_tsne[c("TSNE.1","TSNE.2")],title="Danaher cell-type labels",legend_anno= sort(unique(cluster_tsne[,"name_type"])))+scale_color_manual(values = colorRampPalette(c("blue", "red","yellow"))( length(sort(unique(cluster_tsne[,"name_type"])))))+  annotate("text", x = tsne_center[, 1], y = tsne_center[, 2], size = 5,label=c(1:length(tsne_center[,1])))

all_gene_var_pos <- list()
## danaher marker dispersion 
for(i in 1:length(cell_list)){
type_gene_select <- match(unlist(cell_list[i]),nor_exp$gene$use_genes_n_id)
all_gene_var_pos[[length(all_gene_var_pos)+1]] <- type_gene_select
}
names(all_gene_var_pos) <- cell_list
