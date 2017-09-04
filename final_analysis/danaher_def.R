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
