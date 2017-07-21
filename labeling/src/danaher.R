# This method can only indicates the relative level of each cell type
#Importing gene markers per cell type for Danaher et al. 
library(cellrangerRkit)
library(magrittr)
library(dplyr)
library(reshape2)
library(Matrix)
B_cells <- c( "BLK"
              ,"CD19"
              ,"MS4A1"
              ,"TNFRSF17"
              ,"FCRL2"
              ,"KIAA0125"
              ,"PNOC"
              ,"SPIB"
              ,"TCL1A"
)
CD45 <- c( "PTPRC"
)
CD8_T_cells <- c( "CD8A"
                  ,"CD8B"
)
Cytotoxic_cells <- c( "CTSW"
                      ,"GNLY"
                      ,"GZMA"
                      ,"GZMB"
                      ,"GZMH"
                      ,"KLRB1"
                      ,"KLRD1"
                      ,"KLRK1"
                      ,"PRF1"
                      ,"NKG7"
)
DC <- c( "CCL13"
         ,"CD209"
         ,"HSD11B1"
)
Exhausted_CD8 <- c( "CD244"
                    ,"EOMES"
                    ,"LAG3"
                    ,"PTGER4"
)
Macrophages <- c( "CD163"
                  ,"CD68"
                  ,"CD84"
                  ,"MS4A4A"
)
Mast_cells <- c( "MS4A2"
                 ,"TPSAB1"
                 ,"CPA3"
                 ,"HDC"
                 ,"TPSB2"
)
Neutrophils <- c( "CSF3R"
                  ,"S100A12"
                  ,"CEACAM3"
                  ,"FCAR"
                  ,"FCGR3B"
                  ,"FPR1"
                  ,"SIGLEC5"
)
NK_CD56dim_cells <- c( "IL21R"
                       ,"KIR2DL3"
                       ,"KIR3DL1"
                       ,"KIR3DL2"
)
NK_cells <- c( "NCR1"
               ,"XCL2"
               ,"XCL1"
)
T_cells <- c( "CD3D"
              ,"CD3E"
              ,"CD3G"
              ,"CD6"
              ,"SH2D1A" 
              ,"TRAT1"
)
Th1_cells <- c( "TBX21"
)
Treg <- c( "FOXP3"
)
#Function 
visualize_me<- function(gbm,gene_probes,projection,limits=c(0,10),marker_size=0.1,title=NULL) {
  gbm[gbm<limits[1]] <- limits[1]
  gbm[gbm>limits[2]] <- limits[2]
  gene_values <-gbm
  #colnames(gene_values) <- gene_probes
  projection_names <-  colnames(projection)
  colnames(projection) <- c('Component.1', 'Component.2')
  proj_gene <- data.frame(cbind(projection,gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars=c("Component.1", "Component.2"))
  p<- ggplot(proj_gene_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour=value),size=marker_size) + facet_wrap(~variable) +
    scale_colour_gradient(low="grey",high="red",name = "val") +
    labs(x=projection_names[1],y=projection_names[2])
  if (!is.null(title)) {  p <- p + ggtitle(title) }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

get_signature_matrix <- function(test,c,i){
  all_list_gene <-list()
  for(j in unique(test[test$Cluster==c,]$Barcode)){
    ds <-test[test$Cluster == c & test$Barcode == j & test$Signature == i ,]
    if(sum(as.numeric(ds$Expression))>0){
      Gene_set <-unique(test[test$Cluster==c & test$Signature == i,]$Gene)  
      Gene_name <- ds$Gene
      list_gene <- as.numeric(ds$Expression)
      no_expr_gene <- Gene_set[!(Gene_set %in% Gene_name)] 
      no_expr_gene_list <- rep(0,length(no_expr_gene))
      names(no_expr_gene_list) <- no_expr_gene
      names(list_gene) <- Gene_name
      list_gene <-cbind(t(list_gene),  t(no_expr_gene_list))
      all_list_gene <-rbind(unlist(all_list_gene), unlist(list_gene))
    }}
  return (all_list_gene)}
### To covert cluster type table back to list for list combination 
deconstruct_summary_table <- function(df,n){
  final.list = list()
  tabl = t(df[,-which(names(df) %in% c("Cluster"))])
  tabl =setNames(split(tabl,seq(nrow(tabl))),rownames(tabl))
  for(x in 1: length(tabl)){
    entry = cbind(  as.numeric(unlist(tabl[[x]])),  rep(names(tabl[x]),length(tabl[[x]])),  c(1:length(tabl[[x]])))
  final.list = rbind(final.list, entry) }
  colnames(final.list) = c(n,"Signature","Cluster")
  return(final.list)
}
cell_list <- list(B_cells,CD45,CD8_T_cells,Cytotoxic_cells,DC,Exhausted_CD8,Macrophages,Mast_cells,Neutrophils,NK_CD56dim_cells,NK_cells,T_cells,Th1_cells,Treg)
names(cell_list) <- c("B_cells","CD45","CD8_T_cells","Cytotoxic_cells","DC","Exhausted_CD8","Macrophages","Mast_cells","Neutrophils","NK_CD56dim_cells","NK_cells","T_cells","Th1_cells","Treg"
)
###loading the data set 
gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex')
analysis_results <- load_cellranger_analysis_results("~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex")
#gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/cytotoxic_t_filtered_gene_bc_matrices.mex')
#analysis_results <- load_cellranger_analysis_results('~/bioinfo/Project/labeling/data/cytotoxic_t_filtered_gene_bc_matrices.mex')
cluster_result <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]
ar_cluster <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]$Cluster
tsne_proj <- analysis_results$tsne
gbm1 <- gbm1[,analysis_results$tsne$Barcode]
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
  #computing expression summary for each cluster 
  type_sum <- summary(type_expr)
  type_expr_list <-cbind(Gene  = rownames(type_expr)[type_sum$i], Barcode = colnames(type_expr)[type_sum$j], Expression = type_sum$x,Signature= names(cell_list[i]))
  type_expr_frame <-type_expr_list 
  type_expr_frame <- merge(x=type_expr_frame, y=cluster_result,by="Barcode")
  all_type_expr <- rbind(all_type_expr,type_expr_frame)
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
###
cluster_count <- cluster_result %>% group_by(Cluster) %>% summarise(cell_total_count=sum(n_distinct(Barcode)))
all_type_expr_table <- all_type_expr %>% group_by(Cluster,Signature) %>% summarise(Barcode_count=sum(n_distinct(Barcode)), total_all_gene_Exp= sum(as.numeric(Expression)), avg_non_zero=mean(as.numeric(Expression)) , SD= sd(as.numeric(Expression)) )
all_type_expr_table<-merge(all_type_expr_table,cluster_count,by="Cluster")
all_type_expr_table<- cbind(all_type_expr_table,precent_count =all_type_expr_table$Barcode_count/all_type_expr_table$cell_total_count *100)

all_list_gene <-list()


cluster_mean <- merge(x=cluster_result,y=mean_score,by.y="row.names" ,by.x="Barcode")
cluster_var <-merge(x=cluster_result,y=var_score^2,by.y="row.names" ,by.x="Barcode")
score_by_cluster <- cluster_mean[,-which(names(cluster_mean) %in% c("Barcode"))] %>% group_by(Cluster) %>% summarise_all(funs(mean))
score_by_cluster_sum <- cluster_mean [,-which(names(cluster_mean) %in% c("Barcode"))]%>% group_by(Cluster) %>% summarise_all(funs(sum))
score_by_cluster_var <- cluster_var[,-which(names(cluster_var) %in% c("Barcode"))]%>% group_by(Cluster) %>% summarise_all(funs(sum))
score_by_cluster_var <- sqrt(score_by_cluster_var[,-which(names(score_by_cluster_var) %in% c("Cluster"))])
score_by_cluster_var <- cbind(Cluster=rownames(score_by_cluster_var),score_by_cluster_var)
all_type_expr_table <- merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster,"cell_mean_expression"), by=c("Cluster","Signature") )
all_type_expr_table <-  merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster_sum,"cell_mean_expression_sum"), by=c("Cluster","Signature") )
all_type_expr_table <-  merge(x=all_type_expr_table,y=deconstruct_summary_table(score_by_cluster_var,"cell_mean_expression_sd"), by=c("Cluster","Signature") )
score_by_cluster <- score_by_cluster[, -which(names(score_by_cluster) %in% c("Cluster"))]
cluster_type <-apply(score_by_cluster, 1,function(x) which(x == max(x)))
cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_assignment <- merge(x=cluster_result,y=cluster_type,by="Cluster")
cluster_tsne <-merge(x=cluster_assignment, y=tsne_proj,by="Barcode")

#Cluster centers for TSNE projections
tsne_center <- list()
for(i in 1:length(unique(cluster_tsne$Cluster))){
c.x<-mean(cluster_tsne$TSNE.1[cluster_tsne$Cluster==i])
c.y<-mean(cluster_tsne$TSNE.2[cluster_tsne$Cluster==i])
tsne_center <- rbind(unlist(tsne_center), c(c.x,c.y))
}
