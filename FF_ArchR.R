library(ArchR)
library(hexbin)
library(Biostrings)
library(stringr)
library(stringi)
library(Biostrings)
library(parallel)
set.seed(1)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

tsv_files <- c("fullfusion_fragments-1.tsv.gz", "fullfusion_fragments-2.tsv.gz", "fullfusion_fragments-3.tsv.gz", "fullfusion_fragments-4.tsv.gz", "fullfusion_fragments-5.tsv.gz", "fullfusion_fragments-6.tsv.gz", "fullfusion_fragments-7.tsv.gz", "fullfusion_fragments-8.tsv.gz", "fullfusion_fragments-9.tsv.gz", "fullfusion_fragments-10.tsv.gz", "fullfusion_fragments-11.tsv.gz", "fullfusion_fragments-12.tsv.gz")

names(tsv_files) <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8", "sample9", "sample10", "sample11", "sample12")


FF_ArrowFiles <- createArrowFiles(
  inputFiles = tsv_files,
  sampleNames = names(tsv_files),
  minTSS = 7, 
  minFrags = 30000,
  threads = floor(detectCores()*(5/6)), 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE 
)

#Note I need to create the doublet scores otherwise subsequent ArchR calls will complain
#but I will _not_ ever filter on them or use them
FF_doubScores <- addDoubletScores(
  input = FF_ArrowFiles,
  k = 10,
  threads = floor(detectCores()*(5/6)), 
  knnMethod = "UMAP",
  LSIMethod = 1
)


FF_proj_AllSamples <- ArchRProject(
  ArrowFiles = FF_ArrowFiles,
  threads = floor(detectCores()*(5/6)), 
  outputDirectory = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/",
  copyArrows = FALSE
)

FF_proj_AllSamples_Filt <- FF_proj_AllSamples[which(FF_proj_AllSamples$TSSEnrichment >= 7 & FF_proj_AllSamples$nFrags >= 40000),]


FF_proj_AllSamples_Filt <- addIterativeLSI(ArchRProj = FF_proj_AllSamples_Filt, useMatrix = "TileMatrix", name = "IterativeLSI", threads = floor(detectCores()*(5/6)))

FF_proj_AllSamples_Filt <- addClusters(input = FF_proj_AllSamples_Filt, reducedDims = "IterativeLSI")

FF_proj_AllSamples_Filt <- addUMAP(ArchRProj = FF_proj_AllSamples_Filt, reducedDims = "IterativeLSI")

FF_AllSample_plt <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", labelMeans = FALSE)

plotPDF(FF_AllSample_plt, name = "FF-UMAP-Clusters-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt, addDOC = FALSE, width = 5, height = 5)


sample1_ass <- read.csv("Sample1_merged_associations.csv", header = TRUE)
sample2_ass <- read.csv("Sample2_merged_associations.csv", header = TRUE)
sample3_ass <- read.csv("Sample3_merged_associations.csv", header = TRUE)
sample4_ass <- read.csv("Sample4_merged_associations.csv", header = TRUE)
sample5_ass <- read.csv("Sample5_merged_associations.csv", header = TRUE)
sample6_ass <- read.csv("Sample6_merged_associations.csv", header = TRUE)
sample7_ass <- read.csv("Sample7_merged_associations.csv", header = TRUE)
sample8_ass <- read.csv("Sample8_merged_associations.csv", header = TRUE)
sample9_ass <- read.csv("Sample9_merged_associations.csv", header = TRUE)
sample10_ass <- read.csv("Sample10_merged_associations.csv", header = TRUE)
sample11_ass <- read.csv("Sample11_merged_associations.csv", header = TRUE)
sample12_ass <- read.csv("Sample12_merged_associations.csv", header = TRUE)

FF_proj_AllSamples_Filt@cellColData$cellBC_id = stringr::str_split_fixed(rownames(FF_proj_AllSamples_Filt@cellColData),"-",n=2)[,1]

FF_proj_AllSamples_Filt@cellColData$cellBC = stringr::str_split_fixed(FF_proj_AllSamples_Filt@cellColData$cellBC_id,"#",n=2)[,2]

FF_proj_AllSamples_Filt@cellColData$condition = stringr::str_split_fixed(FF_proj_AllSamples_Filt@cellColData$cellBC_id,"#",n=2)[,1]

newDF = data.frame(cellBC = FF_proj_AllSamples_Filt@cellColData$cellBC, condition = FF_proj_AllSamples_Filt@cellColData$condition, order= 1:length(FF_proj_AllSamples_Filt@cellColData$cellBC))

sample1_ass$revComp <- sapply(sample1_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample2_ass$revComp <- sapply(sample2_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample3_ass$revComp <- sapply(sample3_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample4_ass$revComp <- sapply(sample4_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample5_ass$revComp <- sapply(sample5_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample6_ass$revComp <- sapply(sample6_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample7_ass$revComp <- sapply(sample7_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample8_ass$revComp <- sapply(sample8_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample9_ass$revComp <- sapply(sample9_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample10_ass$revComp <- sapply(sample10_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample11_ass$revComp <- sapply(sample11_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample12_ass$revComp <- sapply(sample12_ass$CBC, function(x) as.character(reverseComplement(DNAString(x))))

sample1_ass_SUB = sample1_ass[sample1_ass $revComp %in% newDF$cellBC,]

sample2_ass_SUB = sample2_ass[sample2_ass $revComp %in% newDF$cellBC,]

sample3_ass_SUB = sample3_ass[sample3_ass $revComp %in% newDF$cellBC,]

sample4_ass_SUB = sample4_ass[sample4_ass $revComp %in% newDF$cellBC,]

sample5_ass_SUB = sample5_ass[sample5_ass $revComp %in% newDF$cellBC,]

sample6_ass_SUB = sample6_ass[sample6_ass $revComp %in% newDF$cellBC,]

sample7_ass_SUB = sample7_ass[sample7_ass $revComp %in% newDF$cellBC,]

sample8_ass_SUB = sample8_ass[sample8_ass $revComp %in% newDF$cellBC,]

sample9_ass_SUB = sample9_ass[sample9_ass $revComp %in% newDF$cellBC,]

sample10_ass_SUB = sample10_ass[sample10_ass $revComp %in% newDF$cellBC,]

sample11_ass_SUB = sample11_ass[sample11_ass $revComp %in% newDF$cellBC,]

sample12_ass_SUB = sample12_ass[sample12_ass $revComp %in% newDF$cellBC,]


sample1_ass_SUB_test = sample1_ass[sample1_ass $revComp %in% newDF$cellBC[newDF$condition == "sample1"],]

sample2_ass_SUB_test = sample2_ass[sample2_ass $revComp %in% newDF$cellBC[newDF$condition == "sample2"],]

sample3_ass_SUB_test = sample3_ass[sample3_ass $revComp %in% newDF$cellBC[newDF$condition == "sample3"],]

sample4_ass_SUB_test = sample4_ass[sample4_ass $revComp %in% newDF$cellBC[newDF$condition == "sample4"],]

sample5_ass_SUB_test = sample5_ass[sample5_ass $revComp %in% newDF$cellBC[newDF$condition == "sample5"],]

sample6_ass_SUB_test = sample6_ass[sample6_ass $revComp %in% newDF$cellBC[newDF$condition == "sample6"],]

sample7_ass_SUB_test = sample7_ass[sample7_ass $revComp %in% newDF$cellBC[newDF$condition == "sample7"],]

sample8_ass_SUB_test = sample8_ass[sample8_ass $revComp %in% newDF$cellBC[newDF$condition == "sample8"],]

sample9_ass_SUB_test = sample9_ass[sample9_ass $revComp %in% newDF$cellBC[newDF$condition == "sample9"],]

sample10_ass_SUB_test = sample10_ass[sample10_ass $revComp %in% newDF$cellBC[newDF$condition == "sample10"],]

sample11_ass_SUB_test = sample11_ass[sample11_ass $revComp %in% newDF$cellBC[newDF$condition == "sample11"],]

sample12_ass_SUB_test = sample12_ass[sample12_ass $revComp %in% newDF$cellBC[newDF$condition == "sample12"],]


newDF$Fusion = "None"

for(i in 1:nrow(newDF)){
if(newDF$condition[i] == "sample1"){
sub = sample1_ass_SUB$Fusion[sample1_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample2"){
sub = sample2_ass_SUB$Fusion[sample2_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample3"){
sub = sample3_ass_SUB$Fusion[sample3_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample4"){
sub = sample4_ass_SUB$Fusion[sample4_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample5"){
sub = sample5_ass_SUB$Fusion[sample5_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample6"){
sub = sample6_ass_SUB$Fusion[sample6_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample7"){
sub = sample7_ass_SUB$Fusion[sample7_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample8"){
sub = sample8_ass_SUB$Fusion[sample8_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample9"){
sub = sample9_ass_SUB$Fusion[sample9_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample10"){
sub = sample10_ass_SUB$Fusion[sample10_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample11"){
sub = sample11_ass_SUB$Fusion[sample11_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
if(newDF$condition[i] == "sample12"){
sub = sample12_ass_SUB$Fusion[sample12_ass_SUB$revComp == newDF$cellBC[i]]
if(length(sub) == 1){ newDF$Fusion[i] = sub }
}
}

newDF$Fusion_modified = newDF$Fusion
newDF$Fusion_modified[is.na(newDF$Fusion_modified) | newDF$Fusion_modified == "None"] = "Unknown"
FF_proj_AllSamples_Filt@cellColData$Fusion = newDF$Fusion_modified


plt_fusion_color <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt, colorBy = "cellColData", name = "Fusion", embedding = "UMAP", labelMeans = FALSE)

plotPDF(plt_fusion_color, name = "Plot-UMAP-FusionColors-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt, addDOC = FALSE, width = 5, height = 5)

# This project now contains only the nuclei that were able to be genotyped and it has the variant assignments for each genotype
                               
FF_proj_AllSamples_Filt_known <- FF_proj_AllSamples_Filt[which(FF_proj_AllSamples_Filt@cellColData$Fusion != "Unknown"),]

plt_fusion_color_known <- plotEmbedding(ArchRProj = FF_proj_AllSamples_Filt_known, colorBy = "cellColData", name = "Fusion", embedding = "UMAP", labelMeans = FALSE)

plotPDF(plt_fusion_color_known, name = "Plot-UMAP-FusionColors-KnownOnly-MyCellRanger.pdf",
        ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE, width = 5, height = 5)

df_forUMAP <- FF_proj_AllSamples_Filt@embeddings$UMAP$df

write.csv(df_forUMAP, file ="/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/FF_df_forUMAP_MyCellRanger.csv")


pathToMacs2 <- findMacs2()

FF_proj_AllSamples_Filt_known <- addGroupCoverages(ArchRProj = FF_proj_AllSamples_Filt_known, groupBy = "Fusion", threads = 38)

FF_proj_AllSamples_Filt_known <- addReproduciblePeakSet(
    ArchRProj = FF_proj_AllSamples_Filt_known, 
    groupBy = "Fusion", 
    pathToMacs2 = pathToMacs2,
    threads = floor(detectCores()*(5/6))
)

FF_proj_AllSamples_Filt_known <- addPeakMatrix(ArchRProj = FF_proj_AllSamples_Filt_known, threads = floor(detectCores()*(5/6)))

FF_proj_AllSamples_Filt_known <- FF_proj_AllSamples_Filt_known[which(FF_proj_AllSamples_Filt_known@cellColData$Fusion != "NUP98-KDM5A" & FF_proj_AllSamples_Filt_known@cellColData$Fusion != "SET-NUP214" & FF_proj_AllSamples_Filt_known@cellColData$Fusion != "EWSR1-SMARCA5" & FF_proj_AllSamples_Filt_known@cellColData$Fusion != "ASPSCR1-TFE3" & FF_proj_AllSamples_Filt_known@cellColData$Fusion != "BCR-JAK2" & FF_proj_AllSamples_Filt_known@cellColData$Fusion != "NEO1-ZFAT"),]

markersPeaks <- getMarkerFeatures(
    ArchRProj = FF_proj_AllSamples_Filt_known, 
    useMatrix = "PeakMatrix", 
    groupBy = "Fusion",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads = floor(detectCores()*(5/6))
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

test_heatmap_mat <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 0,
  binaryClusterRows = TRUE,
  returnMatrix = TRUE
)

write.csv(test_heatmap_mat, "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/Peak_Heatmap_Matrix_MyCellRanger.csv", row.names=T)

test_heatmap <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  nLabel = 0,
  binaryClusterRows = TRUE
)

plotPDF(test_heatmap, name = "Peak-Marker-Heatmap-FF-MyCellRanger", width = 8, height = 6, ArchRProj = FF_proj_AllSamples_Filt_known, addDOC = FALSE)

# Making all pairwise comparisons 

devtools::install_github("GreenleafLab/chromVARmotifs")
FF_proj_AllSamples_Filt_known <- addMotifAnnotations(ArchRProj = FF_proj_AllSamples_Filt_known, motifSet = "homer", name = "Motif")
                               
fusion_list <- c(unique(FF_proj_AllSamples_Filt_known@cellColData$Fusion))

comp_list <- vector(mode='list', length=length(fusion_list))
pma_list <- vector(mode='list', length=length(fusion_list))
pv_list <- vector(mode='list', length=length(fusion_list))
motif_up_list <- vector(mode='list', length=length(fusion_list))
motif_down_list <- vector(mode='list', length=length(fusion_list))
df_down_combo <- NULL
df_up_combo <- NULL

for(i in 1:length(fusion_list)){
comp_list[[i]] <- getMarkerFeatures(
  ArchRProj = FF_proj_AllSamples_Filt_known, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = fusion_list[[i]],
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

pma_list[[i]] <- plotMarkers(seMarker = comp_list[[i]], name = fusion_list[i], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")

pv_list[[i]] <- plotMarkers(seMarker = comp_list[[i]], name = fusion_list[i], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

motifsUp <- peakAnnoEnrichment(
    seMarker = comp_list[[i]],
    ArchRProj = FF_proj_AllSamples_Filt_known,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

temp_df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
temp_df <- temp_df[order(temp_df$mlog10Padj, decreasing = TRUE),]
temp_df$rank <- seq_len(nrow(temp_df))
temp_df$Fusion <- fusion_list[[i]]
df_up_combo <- rbind(df_up_combo, temp_df)

rm(temp_df)

motifsDo <- peakAnnoEnrichment(
    seMarker = comp_list[[i]],
    ArchRProj = FF_proj_AllSamples_Filt_known,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
temp_df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
temp_df <- temp_df[order(temp_df$mlog10Padj, decreasing = TRUE),]
temp_df$rank <- seq_len(nrow(temp_df))
temp_df$Fusion <- fusion_list[[i]]

df_down_combo <- rbind(df_down_combo, temp_df)
rm(temp_df)

}

write.csv(df_down_combo, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/AllComparisons_DownMotifs_Dataframe.csv")

write.csv(df_up_combo, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/AllComparisons_UpMotifs_Dataframe.csv")

pma_df_combo <- NULL
pma_df_combo <- as.data.frame(pma_df_combo)
pma_df_temp <- NULL

for(i in 1:length(pma_list)){
pma_df_temp$Log2FC <- comp_list[[i]]@assays@data$Log2FC
pma_df_temp$FDR <- comp_list[[i]]@assays@data$FDR
pma_df_temp$Mean <- comp_list[[i]]@assays@data$Mean
pma_df_temp$Fusion <- fusion_list[[i]]
pma_df_temp <- as.data.frame(pma_df_temp)
colnames(pma_df_temp) <- c("Log2FC", "FDR", "Mean", "Fusion")
pma_df_combo <- rbind(pma_df_combo, pma_df_temp)
rm(pma_df_temp)
pma_df_temp <- NULL
}

write.csv(pma_df_combo, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/AllComparisons_PMAandPV_Dataframe.csv")
