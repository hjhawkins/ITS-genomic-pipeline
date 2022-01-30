library(gplots)
library(RcolorBrewer)
matrix_envi <- data.matrix(Env)

heatmap.2(matrix_envi, cellnote = matrix_envi)

EnvBsimple <- dplyr::select(EnvB, -Sample.ID, -Vfsand_., -SMBC_ugC_gdw, -V_C_., V_N_., -V_dN15_14, -Csand_., -Vcsand_., -Fsand_.,-Msand_., -S_dN15_14, -V_dC13_12, -Fsand_., -Bdfine_g_cm3, -Bdhybrid_g_cm3, -Residmoist_.airdw, -LOI_.dw, -L550.1000_.dw, -L550_.dw, -S_dC13_12, -Msand_., -totAl_.,-totSi_., -totP_.,-totS_., -totCl_., -totK_., -totCa_., -totTi_., -totV_., -totCr_., -totMn_., -totFe_., -totNi_., -totCu_., -totZn_., -totGa_., -totBr_., -totRb_., -totSr_., -totY_., -totZr_., totNb_., -totBa_., -totHf_., -totTa_., -totPb_., -totTh_.)

B_mat <- data.matrix(EnvBsimple[,12:ncol(EnvBsimple)])

heatmap.2(B_mat, cellnote = B_mat)

EnvBmean<- aggregate(EnvBsimple[,12:33], list(EnvBsimple$Treatm), mean)

row.names(EnvBmean) <- EnvBmean$Group.1
B <- data.matrix(EnvBmean)
B <- decostand(B, "stand") # 
heatmap.2(B, cellnote = B)

heatmap.2(t(B),symkey=FALSE, density.info="none", trace="none" )

#col = c("#FF0000", "#FF0000", "#FF0000")
par(mar=c(7,4,4,2)+0.1) 
#png(filename='test.png', width=800, height=750)
#heatmap.2(t(B), col=redgreen(75), scale="row", ColSideColors=col, key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)

heatmap.2(t(B), scale="row", key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none")




