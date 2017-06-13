# This method can only indicates the relative level of each cell type
#Importing gene markers per cell type for Danaher et al. 
library(cellrangerRkit)
library(magrittr)
library(dplyr)
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
cell_list <- list(B_cells,CD45,CD8_T_cells,Cytotoxic_cells,DC,Exhausted_CD8,Macrophages,Mast_cells,Neutrophils,NK_CD56dim_cells,NK_cells,T_cells,Th1_cells,Treg)
names(cell_list) <- c("B_cells","CD45","CD8_T_cells","Cytotoxic_cells","DC","Exhausted_CD8","Macrophages","Mast_cells","Neutrophils","NK_CD56dim_cells","NK_cells","T_cells","Th1_cells","Treg"
)
###loading the data set 
gbm1 <- load_cellranger_matrix('~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex')
analysis_results <- load_cellranger_analysis_results("~/bioinfo/Project/labeling/data/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.mex")
ar_cluster <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]$Cluster
gbm1 <- gbm1[,analysis_results$tsne$Barcode]
gen= fData(gbm1)
##cell_expr =exprs(gbm1)
##Selecting Genes marker from the expression data 
mean_score <-list()
for(type in cell_list){
  type_gene_select <- gen[match(type,gen$symbol),]$id
  type_gene_select <-type_gene_select[!is.na(type_gene_select)]
  type_expr = gbm1[type_gene_select,]
  ### Paper use log2 Transformation- subjected to change 
  type_mean = log2(colMeans(type_expr))
  mean_score <- cbind(unlist(mean_score), type_mean)
}
colnames(mean_score) <-(1:length(mean_score[1,]))
# tsne_proj <- analysis_results$tsne
# proj_gene <- data.frame(cbind(projection,cell_class))
# proj_gene_melt <- melt(proj_gene, id.vars=c("Component.1", "Component.2"))
ar_cell_assign=list() 
for( i in 1:length(mean_score[,1])) {
  if(as.numeric(max(mean_score[i,])) == -Inf ){type= 999}
  else{ 
    type = unlist(as.numeric(names(which(mean_score[i,]==max(mean_score[i,])))))
    if(length(type) > 1 ){
      type = paste(as.character(sort(type)),collapse = '+')
    }
  }
  ar_cell_assign <- c(unlist(ar_cell_assign), type)
}

##Fisrt parma is just a list of number then cbind with the data frame of tsne for labelling
visualize_clusters(ar_cell_assign,analysis_results$tsne[c("TSNE.1","TSNE.2")],title="Danaher cell-type labels")
cluster_result <- analysis_results$kmeans[[paste(10,"clusters",sep="_")]]
cluster_counts <- cbind(cluster_result, ar_cell_assign)
## determin cell type population in clusters
colnames(mean_score) <-names(cell_list)
clu <- cluster_result[,2]
cluster_mean <- data.frame(clu,mean_score)
cluster_mean %>% group_by(clu) %>% summarise_each(funs(mean))
#bar chart for cell distrubtuion per cluster
count_t <- table(cluster_counts$ar_cell_assign,cluster_counts$Cluster)
barplot(count_t,	legend = rownames(count_t))
## Calculate average TIL score/cell assign type
# total_score_cluster=apply(mean_score>-Inf,1,mean)
# total_score_cluster= cbind(total_score_cluster, cluster_result)
# result <- aggregate(total_score_cluster[,1],by= list(total_score_cluster[,2]), FUN= sum)

CD45 <- lapply(mean_score,ar_cell_assign)
plot(total_score)