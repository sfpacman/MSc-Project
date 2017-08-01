# known_ids <- c("CD19+ B","CD14+ Monocyte","Dendritic","CD56+ NK","CD34+","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
#                "CD4+/CD45RO+ Memory","CD4+ T Helper2","CD8+/CD45RA+ Naive Cytotoxic","CD8+ Cytotoxic T")
# cell_id <-list()
# for(id in known_ids){
# if (id=="CD19+ B") { check_genes <- c("CD74","CD27")} 
# if (id=="CD14+ Monocyte") { check_genes <- c("FTL");       }  # select based on FTL
# if (id=="Dendritic") { check_genes <- c("CLEC9A","CD1C");  }      # select based on CLEC9A
# if (id=="CD34+")  { check_genes <- "CD34";  }# select based on CD34
# if (id=="CD56+ NK") {check_genes <- c("CD3D","NKG7") }
# if (id=="CD4+/CD25 T Reg") { check_genes <- "NKG7" }
# if (id=="CD4+/CD45RA+/CD25- Naive T") { check_genes <- "NKG7" }
# if (id=="CD4+/CD45RO+ Memory") { check_genes <- "NKG7" }
# if (id=="CD4+ T Helper2") { check_genes <- "NKG7" }
# if (id=="CD8+/CD45RA+ Naive Cytotoxic") { check_genes <- "NKG7"} 
# if (id=="CD8+ Cytotoxic T") { check_genes <- c("CD8A","CD3D")}
#   cell_id[id] <- check_genes 
# }
# Marker Def for 10X 
myeloid_cells<- c( S100A8, S100A9)
B_cells<- c( CD79A)
NK_cells<- c( NKG7)
dendritic_cells <- c( FCER1A)
T_cells<- c( CD3D)
memory_T_cells<- c( CCR10)
regulatory_T_cells<-c(TNFRSF18)
CD8 <- c(CD8A)
CD4 <- c(CD4)
naive_T-cell<- c( ID3)
Activated_cytotoxic_T_cells<- c( NKG7)


