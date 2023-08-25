#This assumes you start with an ArchR project and then downsamples a few variants of interest 
library(ArchR); library(hexbin); library(Biostrings); library(stringr)
library(stringi); library(Biostrings); library(ggplot2); library(parallel)
library(reshape2); library(viridis)
set.seed(1)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

proj <- loadArchRProject(path = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/", force = FALSE, showLogo = TRUE)

EWSR1_YY1 <- proj[which(proj@cellColData$Fusion == "EWSR1-YY1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(EWSR1_YY1@cellColData$Fusion == "EWSR1-YY1"))
n_empty <- length(which(EWSR1_YY1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(EWSR1_YY1@cellColData$Fusion == "EWSR1-YY1")
empty_rows <- which(EWSR1_YY1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-YY1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- EWSR1_YY1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-YY1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/EWSR1YY1_downsample_df.csv")
rm(full_comp)
rm(EWSR1_YY1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

EWSR1_ATF1 <- proj[which(proj@cellColData$Fusion == "EWSR1-ATF1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(EWSR1_ATF1@cellColData$Fusion == "EWSR1-ATF1"))
n_empty <- length(which(EWSR1_ATF1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(EWSR1_ATF1@cellColData$Fusion == "EWSR1-ATF1")
empty_rows <- which(EWSR1_ATF1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-ATF1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- EWSR1_ATF1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-ATF1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/EWSR1ATF1_downsample_df.csv")
rm(full_comp)
rm(EWSR1_ATF1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

CCDC6_RET <- proj[which(proj@cellColData$Fusion == "CCDC6-RET" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(CCDC6_RET@cellColData$Fusion == "CCDC6-RET"))
n_empty <- length(which(CCDC6_RET@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(CCDC6_RET@cellColData$Fusion == "CCDC6-RET")
empty_rows <- which(CCDC6_RET@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CCDC6-RET",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- CCDC6_RET[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CCDC6-RET",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/CCDC6RET_downsample_df.csv")
rm(full_comp)
rm(CCDC6_RET)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

ARFGEF2_HNF4A <- proj[which(proj@cellColData$Fusion == "ARFGEF2-HNF4A" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(ARFGEF2_HNF4A@cellColData$Fusion == "ARFGEF2-HNF4A"))
n_empty <- length(which(ARFGEF2_HNF4A@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(ARFGEF2_HNF4A@cellColData$Fusion == "ARFGEF2-HNF4A")
empty_rows <- which(ARFGEF2_HNF4A@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ARFGEF2-HNF4A",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- ARFGEF2_HNF4A[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ARFGEF2-HNF4A",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/ARFGEF2HNF4A_downsample_df.csv")
rm(full_comp)
rm(ARFGEF2_HNF4A)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

IRF2BP2_CDX1 <- proj[which(proj@cellColData$Fusion == "IRF2BP2-CDX1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(IRF2BP2_CDX1@cellColData$Fusion == "IRF2BP2-CDX1"))
n_empty <- length(which(IRF2BP2_CDX1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(IRF2BP2_CDX1@cellColData$Fusion == "IRF2BP2-CDX1")
empty_rows <- which(IRF2BP2_CDX1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "IRF2BP2-CDX1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- IRF2BP2_CDX1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "IRF2BP2-CDX1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/IRF2BP2CDX1_downsample_df.csv")
rm(full_comp)
rm(IRF2BP2_CDX1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

PAX7_FOXO1 <- proj[which(proj@cellColData$Fusion == "PAX7-FOXO1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(PAX7_FOXO1@cellColData$Fusion == "PAX7-FOXO1"))
n_empty <- length(which(PAX7_FOXO1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(PAX7_FOXO1@cellColData$Fusion == "PAX7-FOXO1")
empty_rows <- which(PAX7_FOXO1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PAX7-FOXO1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- PAX7_FOXO1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PAX7-FOXO1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/PAX7FOXO1_downsample_df.csv")
rm(full_comp)
rm(PAX7_FOXO1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

ACTB_GLI1 <- proj[which(proj@cellColData$Fusion == "ACTB-GLI1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(ACTB_GLI1@cellColData$Fusion == "ACTB-GLI1"))
n_empty <- length(which(ACTB_GLI1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(ACTB_GLI1@cellColData$Fusion == "ACTB-GLI1")
empty_rows <- which(ACTB_GLI1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ACTB-GLI1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- ACTB_GLI1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ACTB-GLI1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/ACTBGLI1_downsample_df.csv")
rm(full_comp)
rm(ACTB_GLI1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Next fusion

EWSR1_FLI1 <- proj[which(proj@cellColData$Fusion == "EWSR1-FLI1" | proj@cellColData$Fusion == "Empty"),]

n_fusion <- length(which(EWSR1_FLI1@cellColData$Fusion == "EWSR1-FLI1"))
n_empty <- length(which(EWSR1_FLI1@cellColData$Fusion == "Empty"))

fus_sample_list <- seq(10, n_fusion, 20)
empty_sample_list <- seq(10, 1520, 100)
fus_rows <- which(EWSR1_FLI1@cellColData$Fusion == "EWSR1-FLI1")
empty_rows <- which(EWSR1_FLI1@cellColData$Fusion == "Empty")

full_comp <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-FLI1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df <- NULL
downsample_df_temp <- NULL

for(fus in fus_sample_list){
for(emp in empty_sample_list){
fus_sample <- sample(fus_rows, fus)
empty_sample <- sample(empty_rows, emp)
temp <- EWSR1_FLI1[c(fus_sample, empty_sample),]
temp_comp <- getMarkerFeatures(
  ArchRProj = temp, 
  useMatrix = "PeakMatrix",
  groupBy = "Fusion",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "EWSR1-FLI1",
  bgdGroups = "Empty",
  threads = floor(detectCores()*(5/6))
)

downsample_df_temp$n_up <- length(which(temp_comp@assays@data$Log2FC >= 1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_down <- length(which(temp_comp@assays@data$Log2FC <= -1 & temp_comp@assays@data$FDR <= 0.1))
downsample_df_temp$n_total <- downsample_df_temp$n_up + downsample_df_temp$n_down

downsample_df_temp$pearson <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "pearson")
downsample_df_temp$spearman <- cor(temp_comp@assays@data$Log2FC, full_comp@assays@data$Log2FC, method = "spearman")

downsample_df_temp$cell_fus <- fus
downsample_df_temp$cell_emp <- emp

downsample_df_temp <- data.frame(downsample_df_temp)
colnames(downsample_df_temp) <- c("n_up", "n_down", "n_total", "pearson", "spearman", "Fusion_cells", "Empty_cells")


downsample_df <- rbind(downsample_df, downsample_df_temp)

rm(downsample_df_temp)
downsample_df_temp <- NULL
rm(temp)
rm(temp_comp)
rm(fus_sample)
rm(empty_sample)
}
}


write.csv(downsample_df, file = "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/EWSR1FLI1_downsample_df.csv")
rm(full_comp)
rm(EWSR1_FLI1)
rm(fus_sample)
rm(empty_sample)
rm(n_fusion)
rm(n_empty)
rm(fus_rows)
rm(empty_rows)

#Only making plots for one example but obviously this can be extended to any number

ewsr1yy1_downsample <- read.csv("/Users/mfrenkel/Desktop/FusionDownsample/EWSR1YY1_downsample_df.csv")
ewsr1yy1_downsample <- ewsr1yy1_downsample[,-1]
for(i in 1:nrow(ewsr1yy1_downsample)){
  ewsr1yy1_downsample$frac_total[i] <- ewsr1yy1_downsample$n_total[i] / max(ewsr1yy1_downsample$n_total)
}
ewsr1yy1_pearson <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='pearson')
ewsr1yy1_spearman <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='spearman')
ewsr1yy1_total <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='n_total')
ewsr1yy1_up <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='n_up')
ewsr1yy1_down <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='n_down')
ewsr1yy1_frac <- acast(ewsr1yy1_downsample, Fusion_cells~Empty_cells, value.var='frac_total')

Heatmap(ewsr1yy1_pearson, cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(title = "Pearson correlation 
of Log2(FC)", direction = "horizontal", title_position = "topcenter"), col=mako(100))

Heatmap(ewsr1yy1_total, cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(title = "Number of differentially
accessible peaks", direction = "horizontal", title_position = "topcenter"), col=mako(100))

Heatmap(ewsr1yy1_frac, cluster_rows = FALSE, cluster_columns = FALSE,
                     heatmap_legend_param = list(title = "Number of differentially
accessible peaks", direction = "horizontal", title_position = "topcenter", legend_width = unit(4, "cm")), col=mako(100))

#Of course can combine this dataframe with others of all the other fusions to make comparisons etc. 




