
library(Matrix)
library(ggplot2)
library(dplyr)
source("danaher_def.R")
# Five combination of normalization and clustering are applied and analyzed by the following codes:correlation, UMI+kmean, computeSumFactors+kmeans, UMI+kmeans, UMI+SNN-cliq, computeSumFactors+ SNN-cliq 
# computeSumFactors+ SNN-cliq have two set of data (k = 11 and k = 5)
### This is the example of correlation (Figure. 3A and Supp.table 1) 
sample_analysis <-readRDS("sample_analysis.rds")
sample <- readRDS("sample.rds")
rcell_assign <- sample_analysis$X$cls_id
act <- unlist(sample$summary[,2])
result <- data.frame(pre=rcell_assign,act =act)
x<- result %>% group_by(act,pre) %>% summarise (n = n()) %>% mutate(Precentage = n / sum(n)*100)
# This produces supplmentary table 1   
kable(x,"latex",col.names=c("Actual Identity", "Prediction","cell count","% cell count" ))
freq <- as.numeric(table(result$pre)) 
# This produces figure 3A 
pie(freq, labels = paste(cluster_type,signif(clu$S_s_11[,2]/sum(clu$S_s_11[,2])*100,2),"%",sep=" "))
### This is the example UMI+kmean (Figure 3B and Supp. table 3)
### The rest of Danher's marker results are produced in this methods by change m_n (expression matrix) and cluster_assign (cluster assignment) 
### necessary for getting the gene name for danaher's markers 
pbmc68k <- readRDS("pbmc68k.rds")
gen <- pbmc68k$all_data[[1]]$hg19$gene_symbols
ngn <- sapply(use_gene, function(x){gen[x]})
m_n <- sample_analysis$X$$m_n
### mean score computing ### 
mean_score <- list()
for(i in 1:length(cell_list)){
  type= cell_list[i]
  type_gene_select <- match(unlist(cell_list[i]),ngn)
  type_gene_select <- type_gene_select[!is.na(type_gene_select)]
  type_expr <- m_n[,type_gene_select]
  type_mean <-colMeans(t(type_expr))
  mean_score <- cbind(unlist(mean_score), type_mean)
}

colnames(mean_score) <- names(cell_list)
cluster_assign <- sample_analysis$X$k$cluster
cluster_mean <- data.frame(Cluster= cluster_assign,mean_score)

score_by_cluster <- round(cluster_mean %>% group_by(Cluster) %>% summarise_all(funs(mean)),3)
score_by_cluster <- score_by_cluster[,-which(names(score_by_cluster) %in% c("Cluster"))]
### selecting cell identity for each cluster ###
cluster_type <- list()
for( i in 1:nrow(score_by_cluster)){
  x <- as.numeric(score_by_cluster[i,])
  if( mean(x) == 0 ) {cluster_type = c(unlist(cluster_type),0) }
  else {cluster_type = c(unlist(cluster_type), which(x== max(x)))} } 

rcell_assign <-sapply(cluster_mean$Cluster,function(x){cluster_type[x]})
rcell_assign <- factor(rcell_assign,levels=c(0:length(cell_list)),labels= c("unknown",names(cell_list)))
act <- unlist(sample$summary[,2])
result <- data.frame(pre=rcell_assign,act =act)

x<- result %>% group_by(act,pre) %>% summarise (n = n()) %>% mutate(Precentage = n / sum(n)*100)
# This produces supplmentary table 3
kable(x,"latex",col.names=c("Actual Identity", "Prediction","cell count","% cell count" ))
#This produces figure 3B
freq <- as.numeric(table(result$pre)) 
pie(freq, labels = paste(cluster_type,signif(clu$S_s_11[,2]/sum(clu$S_s_11[,2])*100,2),"%",sep=" "))
