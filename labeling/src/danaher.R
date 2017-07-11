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
cell_list <- list(B_cells,CD45,CD8_T_cells,Cytotoxic_cells,DC,Exhausted_CD8,Macrophages,Mast_cells,Neutrophils,NK_CD56dim_cells,NK_cells,T_cells,Th1_cells,Treg)
names(cell_list) <- c("B_cells","CD45","CD8_T_cells","Cytotoxic_cells","DC","Exhausted_CD8","Macrophages","Mast_cells","Neutrophils","NK_CD56dim_cells","NK_cells","T_cells","Th1_cells","Treg"
)
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
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
###loading the data set 
gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex')
analysis_results <- load_cellranger_analysis_results("~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex")
cluster_result <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]
ar_cluster <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]$Cluster
tsne_proj <- analysis_results$tsne
gbm1 <- gbm1[,analysis_results$tsne$Barcode]
gen= fData(gbm1)
##cell_expr =exprs(gbm1)
##Selecting Genes marker from the expression data 
mean_score <-list()
type_score_table <-list()
all_type_expr <-list()
for(i in 1:length(cell_list)){
  type= cell_list[i]
  type_gene_select <- gen[match(unlist(type),gen$symbol),]$id
  type_gene_select <-type_gene_select[!is.na(type_gene_select)]
  type_expr <- exprs(gbm1[type_gene_select,])
  type_sum <- summary(type_expr)
  type_expr_list <-cbind(Gene  = rownames(type_expr)[type_sum$i], Barcode = colnames(type_expr)[type_sum$j], Expression = type_sum$x,Signature= names(cell_list[i]))
  type_expr_frame <-type_expr_list 
  #type_expr_frame <- data.frame(Gene  = rownames(type_expr)[type_sum$i], Barcode = colnames(type_expr)[type_sum$j], Expression = type_sum$x )
  type_expr_frame <- merge(x=type_expr_frame, y=cluster_result,by="Barcode")
  all_type_expr <- rbind(all_type_expr,type_expr_frame)
  #all_type_expr <- c(all_type_expr,list(type_expr_frame))
  
  ### Paper use log2 Transformation- subjected to change 
  #type_mean = log2(colMeans(type_expr))
  type_mean <- colMeans(type_expr)
  mean_score <- cbind(unlist(mean_score), type_mean)
}
###
cluster_count <- cluster_result %>% group_by(Cluster) %>% summarise(cell_total_count=sum(n_distinct(Barcode)))
all_type_expr_table <- all_type_expr %>% group_by(Cluster,Signature) %>% summarise(Barcode_count=sum(n_distinct(Barcode)), total_Exp= sum(as.numeric(Expression)), avg_non_zero=mean(as.numeric(Expression)) , SD= sd(as.numeric(Expression)) )
all_type_expr_table<-merge(all_type_expr_table,cluster_count,by="Cluster")
all_type_expr_table<- cbind(all_type_expr_table,precent_count =all_type_expr_table$Barcode_count/all_type_expr_table$cell_total_count *100)
all_list_gene <-list()


#Indivdual cell assignment
colnames(mean_score) <- names(cell_list)
# ar_cell_assign=list()
# for( i in 1:length(mean_score[,1])) {
#   if(as.numeric(max(mean_score[i,])) == 0 ){type= 999}
#   else{
#     type = unlist(as.numeric(names(which(mean_score[i,]==max(mean_score[i,])))))
#     if(length(type) > 1 ){
#       type = paste(as.character(sort(type)),collapse = '+')
#     }
#   }
#   ar_cell_assign <- c(unlist(ar_cell_assign), type)
# }
# for( i in 1:length(mean_score[1,])){
# visualize_me(mean_score[,i],cell_list[i],analysis_results$tsne[c("TSNE.1","TSNE.2")],title=names(cell_list)[i] )
# }
# 
# ##Fisrt parma is just a list of number then cbind with the data frame of tsne for labelling
# # visualize_clusters(ar_cell_assign,analysis_results$tsne[c("TSNE.1","TSNE.2")],title="Danaher cell-type labels")

# # cluster_counts <- cbind(cluster_result, ar_cell_assign)
# ## determin cell type population in clusters
#colnames(mean_score) <-names(cell_list)
cluster_mean <- data.frame(ar_cluster,mean_score)
score_by_cluster <- cluster_mean %>% group_by(ar_cluster) %>% summarise_all(funs(mean))
score_by_cluster <- score_by_cluster[, -which(names(score_by_cluster) %in% c("ar_cluster"))]
cluster_type <-apply(score_by_cluster, 1,function(x) which(x == max(x)))
cluster_type <-data.frame(Cluster=1:length(cluster_type),cell_type=cluster_type, name_type=names(cell_list)[cluster_type])
cluster_assignment <- merge(x=cluster_result,y=cluster_type,by="Cluster")
cluster_tsne <-merge(x=cluster_assignment, y=tsne_proj,by="Barcode")

# #bar chart for cell distrubtuion per cluster
# count_t <- table(cluster_counts$ar_cell_assign,cluster_counts$Cluster)
# barplot(count_t,	legend = rownames(count_t))
# ## Calculate average TIL score/cell assign type
# total_score_cluster=apply(mean_score>0,1,mean)
# # total_score_cluster= cbind(total_score_cluster, cluster_result)
# result <- aggregate(total_score_cluster[,1],by= list(total_score_cluster[,2]), FUN= sum)
# 
# CD45 <- lapply(mean_score,ar_cell_assign)
# plot(total_score)