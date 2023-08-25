library(scales); library(ggplot2); library(gplots); library(reshape2); library(tidyr); library(ComplexHeatmap)

full <- read.csv("/Users/mfrenkel/Desktop/MyCellRanger/FF_motif_significances_full_df_more.csv")
full[is.na(full)] <- 0

motif_sum <- NULL
motif_sum <- as.data.frame(motif_sum)
motif_list <- unique(full$TF)
for(motif in motif_list){
  if(max(full[which(full$TF == motif),]$normalized_p) > 0.9){
    motif_sum <- rbind(motif_sum, full[which(full$TF == motif),])
  }
}

motif_sum2 <- motif_sum[,c(2,5,6)]

formatted_df <- tidyr::pivot_wider(motif_sum2, values_from = normalized_p, id_cols = Fusion, names_from = TF)
formatted_df <- as.data.frame(formatted_df)
rownames(formatted_df) <- formatted_df$Fusion
formatted_df <- formatted_df[,-1]
formatted_df <- data.matrix(formatted_df)
formatted_df2 <- formatted_df[-which(row.names(test3) == "NEK7-SMYD3" | row.names(test3) == "KIF26B-SMYD3"),]

Colors=c("white","#74698C")
Colors=colorRampPalette(Colors)(100)

pdf("/Users/mfrenkel/Desktop/MyCellRanger/FF_Marker_Motif_Heatmap_Dendrogram.pdf", width = 14, height = 8)
Heatmap(formatted_df2, col = Colors, show_row_dend = T, 
        heatmap_legend_param = list(title = ""),
        row_names_side = "left", column_names_rot = 45)
dev.off()


