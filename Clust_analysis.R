cluster_df <- read.csv("/Users/mfrenkel/Desktop/MyCellRanger/cluster_identity_df.csv") #This came from ArchR and just has cluster identity for each nuclei
cluster_df <- cluster_df[,-1]
colnames(cluster_df) <- c("Cell","Fusion", "Cluster")

total <- cluster_df %>% group_by(Cluster) %>% summarize(clust_total_cells = n())

grouped <- cluster_df %>% group_by(Cluster,Fusion) %>% summarize(number = n())

merge_df <- merge(total, grouped, by.x = "Cluster", by.y = "Cluster", all.y = TRUE)

merge_df$PercentFusionPerCluster <- merge_df$number / merge_df$clust_total_cells

colourCount = length(unique(merge_df$Fusion))
getPalette = colorRampPalette(brewer.pal(10, "Paired"))(colourCount)

fusion_dist <- ggplot(merge_df, aes(fill=Fusion, y=PercentFusionPerCluster, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                    legend.position = "none") +
  scale_fill_manual(values = getPalette) + ylab("Percentage of nuclei assigned to each variant")
ggsave(fusion_dist, filename = "/Users/mfrenkel/Desktop/FusionsInClusters.pdf", 
       width = 14, height = 8, units ="in")

total2 <- cluster_df %>% group_by(Fusion) %>% summarize(fusion_total_cells = n())
merge_df2 <- merge(total2, grouped, by.x = "Fusion", by.y = "Fusion", all.y = TRUE)
merge_df2$PercentofClusterInFusion <- merge_df2$number / merge_df2$fusion_total_cells

colourCount2 = length(unique(merge_df$Cluster))
getPalette2 = colorRampPalette(brewer.pal(10, "Paired"))(colourCount2)
getPalette3 = c("#A6CEE3", "#559AC6", "#6A3D9A", "#94CA92", "#7FC564", "#33A02C",
                "#AB9C6D", "#F6807F", "#3C8CAB", "#ED5B3D", "#FDBF6F", "#FE982C",
                "#F4892A", "#D4A7AB", "#A383BD", "#E73335")

merge_df2_temp <- merge_df2 %>% group_by(Fusion) %>% summarise(number = sum(PercentofClusterInFusion[Cluster %in% c("C2", "C3")]))
merge_df2$Fusion <- factor(merge_df2$Fusion, levels = merge_df2_temp$Fusion[order(merge_df2_temp$number)])
ggplot(merge_df2, aes(fill=Cluster, y=PercentofClusterInFusion, x=Fusion)) + 
  geom_bar(position="fill", stat="identity", width = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = getPalette3) + ylab("Percentage of nuclei assigned to each cluster")

clust_dist <- ggplot(merge_df2, aes(fill=Cluster, y=PercentofClusterInFusion, x=Fusion)) + 
  geom_bar(position="fill", stat="identity", width = 0.98) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = getPalette3) + ylab("Percentage of nuclei assigned to each cluster") +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))
ggsave(clust_dist, filename = "/Users/mfrenkel/Desktop/MyCellRanger/ClustersPerFusions_Wide.pdf", 
       width = 30, height = 10, units ="in")

for(i in 1:nrow(cluster_df)){
  cluster_df$name[i] <-paste0(cluster_df$Cell[i], "-1")
}

clust_merged <- merge(cluster_df, df_UMAP_known, by.x = "name", by.y = "X")

getPalette3 = c("#A6CEE3", "#559AC6", "#6A3D9A", "#94CA92", "#7FC564", "#33A02C",
                "#AB9C6D", "#F6807F", "#3C8CAB", "#ED5B3D", "#FDBF6F", "#FE982C",
                "#F4892A", "#D4A7AB", "#A383BD", "#E73335")

clust_plt <- ggplot(clust_merged, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, col=factor(Cluster))) + geom_point(size =1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.position = "none") +
  xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") +
  scale_color_manual(values = c("C1" = "#A6CEE3",
                                "C10"= "#559AC6",
                                "C11" = "#6A3D9A",
                                "C12" = "#94CA92",
                                "C13" = "#7FC564",
                                "C14" = "#33A02C",
                                "C15" = "#AB9C6D",
                                "C16" = "#F6807F",
                                "C2" = "#3C8CAB",
                                "C3" = "#ED5B3D",
                                "C4" = "#FDBF6F",
                                "C5" = "#FE982C",
                                "C6" = "#F4892A",
                                "C7" = "#D4A7AB",
                                "C8" = "#A383BD",
                                "C9" = "#E73335")) 
ggsave(clust_plt, filename = "/Users/mfrenkel/Desktop/MyCellRanger/UMAP_ColoredByCluster.pdf", width = 9, height = 8)
