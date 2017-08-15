### This is alternative code for processing data from pbmc68K.rds  
### some of the analysis and projection(tsne,kmean,dispersion rank) require data from 10X workflow which is stored in nor_exp 

gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols
ngn <- sapply(nor_10X$use_gene, function(x){gen[x]})
mean_score <- list()
### Two tasks are perform in this for loop: mean score calcluation, cluster analysis table  
for(i in 1:length(cell_list)){
 type= cell_list[i]
 type_gene_select <- match(unlist(cell_list[i]),ngn)
 type_gene_select <- type_gene_select[!is.na(type_gene_select)]
 type_expr <-nor_10X$m[,type_gene_select]
 type_mean <-colMeans(t(type_expr))
 mean_score <- cbind(unlist(mean_score), type_mean)
 
 #computing expression summary for each cluster 
 # type_sum <- summary(type_expr)
 #type_expr_list <-cbind(Gene  = type_sum$i, Barcode = type_sum$j, Expression = type_sum$x,Signature= names(cell_list[i]))
 # type_expr_frame <-type_expr_list 
 # type_expr_frame <- merge(x=type_expr_frame, y=nor_exp$k$Cluster, by)
 # all_type_expr <- rbind(all_type_expr,type_expr_frame)   
 }
colnames(mean_score) <- names(cell_list)
cluster_mean <- data.frame(cbind(Cluster=nor_exp$k$cluster,mean_score))

score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
cluster_type <-apply(score_by_cluster[,-which(names(cluster_mean) %in% c("Cluster"))], 1,function(x) which(x == max(x)))

cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_tsne <-merge( x=data.frame(cbind(cluster=nor_exp$k$cluster,nor_exp$tsne)) , y=cluster_type, by.x="cluster", by.y="Cluster")
# reassign cluster 
#pnorm(x, mean= 1.043, sd = sd(x[,cluster_type$cell_type[6]]))
 
# mean expression for top 1000 vaiable genes
expression_cluster <-data.frame( exp=nor_exp$exp,clu=nor_exp$k$cluster)
reference_expression <- expression_cluster %>% group_by(clu) %>% summarise_all(funs(mean))

all_gene_var_pos <- list()
## danaher marker dispersion: disperion rank data comes from 10X analysis 
for(i in 1:length(cell_list)){
type_gene_select <- match(unlist(cell_list[i]),nor_exp$gene$use_genes_n_id)
all_gene_var_pos[[length(all_gene_var_pos)+1]] <- type_gene_select
}
names(all_gene_var_pos) <- names(cell_list)
for(x in 1:length(cell_list)){names(all_gene_var_pos[[x]]) <- cell_list[[x]]}

### cluster labeling                     
tsne_center <- list()
for(i in 1:length(unique(cluster_tsne$cluster))){
c.x<-mean(cluster_tsne$X1[cluster_tsne$cluster==i])
c.y<-mean(cluster_tsne$X2[cluster_tsne$cluster==i])
tsne_center <- rbind(unlist(tsne_center), c(c.x,c.y))
}
sig_plot <- visualize_clusters(cluster_tsne$cell_type,cluster_tsne[c("X1","X2")],title="Danaher cell-type labels",legend_anno= sort(unique(cluster_tsne[,"name_type"])))+scale_color_manual(values = colorRampPalette(c("blue", "red","yellow"))( length(sort(unique(cluster_tsne[,"name_type"])))))+  annotate("text", x = tsne_center[, 1], y = tsne_center[, 2], size = 5,label=c(1:length(tsne_center[,1])))

