#Takes all of ChIP-Atlas Enrichment output to make a heatmap

library(tidyverse); library(ggdendro)
library(cowplot); library(ggtree)

AG <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/ACTBGLI1.txt")
AG <- AG[,-c(2,3)]
colnames(AG) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                         "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                         "Log_q_val", "Fold_enrichment")
AG$Fusion_shared_peak_count <- as.numeric(str_split_fixed(AG$Overlaps_with_fusion, pattern = "/", n =2)[,1])
AG$Fusion_input_peak_count <- as.numeric(str_split_fixed(AG$Overlaps_with_fusion, pattern = "/", n =2)[,2])
AG$Control_shared_peak_count <- as.numeric(str_split_fixed(AG$Overlaps_with_control, pattern = "/", n =2)[,1])
AG$Frac_fusionpeaks_in_cell <- AG$Fusion_shared_peak_count/AG$Fusion_input_peak_count
AG$ComboID <- paste0(AG$ID, AG$Cell_class, AG$Cell_type)
AG$Fusion <- "ACTB-GLI1"

AH <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/ARFGEF2HNF4A.txt")
AH <- AH[,-c(2,3)]
colnames(AH) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
AH$Fusion_shared_peak_count <- as.numeric(str_split_fixed(AH$Overlaps_with_fusion, pattern = "/", n =2)[,1])
AH$Fusion_input_peak_count <- as.numeric(str_split_fixed(AH$Overlaps_with_fusion, pattern = "/", n =2)[,2])
AH$Control_shared_peak_count <- as.numeric(str_split_fixed(AH$Overlaps_with_control, pattern = "/", n =2)[,1])
AH$Frac_fusionpeaks_in_cell <- AH$Fusion_shared_peak_count/AH$Fusion_input_peak_count
AH$ComboID <- paste0(AH$ID, AH$Cell_class, AH$Cell_type)
AH$Fusion <- "ARFGEF2-HNF4A"

CR <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/CCDC6RET.txt")
CR <- CR[,-c(2,3)]
colnames(CR) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
CR$Fusion_shared_peak_count <- as.numeric(str_split_fixed(CR$Overlaps_with_fusion, pattern = "/", n =2)[,1])
CR$Fusion_input_peak_count <- as.numeric(str_split_fixed(CR$Overlaps_with_fusion, pattern = "/", n =2)[,2])
CR$Control_shared_peak_count <- as.numeric(str_split_fixed(CR$Overlaps_with_control, pattern = "/", n =2)[,1])
CR$Frac_fusionpeaks_in_cell <- CR$Fusion_shared_peak_count/CR$Fusion_input_peak_count
CR$ComboID <- paste0(CR$ID, CR$Cell_class, CR$Cell_type)
CR$Fusion <- "CCDC6-RET"

EP <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EPC1PHF1.txt")
EP <- EP[,-c(2,3)]
colnames(EP) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EP$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EP$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EP$Fusion_input_peak_count <- as.numeric(str_split_fixed(EP$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EP$Control_shared_peak_count <- as.numeric(str_split_fixed(EP$Overlaps_with_control, pattern = "/", n =2)[,1])
EP$Frac_fusionpeaks_in_cell <- EP$Fusion_shared_peak_count/EP$Fusion_input_peak_count
EP$ComboID <- paste0(EP$ID, EP$Cell_class, EP$Cell_type)
EP$Fusion <- "EPC1-PHF1"

EN <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/ETV6NTRK3.txt")
EN <- EN[,-c(2,3)]
colnames(EN) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EN$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EN$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EN$Fusion_input_peak_count <- as.numeric(str_split_fixed(EN$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EN$Control_shared_peak_count <- as.numeric(str_split_fixed(EN$Overlaps_with_control, pattern = "/", n =2)[,1])
EN$Frac_fusionpeaks_in_cell <- EN$Fusion_shared_peak_count/EN$Fusion_input_peak_count
EN$ComboID <- paste0(EN$ID, EN$Cell_class, EN$Cell_type)
EN$Fusion <- "ETV6-NTRK3"

EA <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1ATF1.txt")
EA <- EA[,-c(2,3)]
colnames(EA) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EA$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EA$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EA$Fusion_input_peak_count <- as.numeric(str_split_fixed(EA$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EA$Control_shared_peak_count <- as.numeric(str_split_fixed(EA$Overlaps_with_control, pattern = "/", n =2)[,1])
EA$Frac_fusionpeaks_in_cell <- EA$Fusion_shared_peak_count/EA$Fusion_input_peak_count
EA$ComboID <- paste0(EA$ID, EA$Cell_class, EA$Cell_type)
EA$Fusion <- "EWSR1-ATF1"

EC <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1CREB1.txt")
EC <- EC[,-c(2,3)]
colnames(EC) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EC$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EC$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EC$Fusion_input_peak_count <- as.numeric(str_split_fixed(EC$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EC$Control_shared_peak_count <- as.numeric(str_split_fixed(EC$Overlaps_with_control, pattern = "/", n =2)[,1])
EC$Frac_fusionpeaks_in_cell <- EC$Fusion_shared_peak_count/EC$Fusion_input_peak_count
EC$ComboID <- paste0(EC$ID, EC$Cell_class, EC$Cell_type)
EC$Fusion <- "EWSR1-CREB1"

ED <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1DDIT3.txt")
ED <- ED[,-c(2,3)]
colnames(ED) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
ED$Fusion_shared_peak_count <- as.numeric(str_split_fixed(ED$Overlaps_with_fusion, pattern = "/", n =2)[,1])
ED$Fusion_input_peak_count <- as.numeric(str_split_fixed(ED$Overlaps_with_fusion, pattern = "/", n =2)[,2])
ED$Control_shared_peak_count <- as.numeric(str_split_fixed(ED$Overlaps_with_control, pattern = "/", n =2)[,1])
ED$Frac_fusionpeaks_in_cell <- ED$Fusion_shared_peak_count/ED$Fusion_input_peak_count
ED$ComboID <- paste0(ED$ID, ED$Cell_class, ED$Cell_type)
ED$Fusion <- "EWSR1-DDIT3"

EV1 <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1ETV1.txt")
EV1 <- EV1[,-c(2,3)]
colnames(EV1) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EV1$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EV1$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EV1$Fusion_input_peak_count <- as.numeric(str_split_fixed(EV1$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EV1$Control_shared_peak_count <- as.numeric(str_split_fixed(EV1$Overlaps_with_control, pattern = "/", n =2)[,1])
EV1$Frac_fusionpeaks_in_cell <- EV1$Fusion_shared_peak_count/EV1$Fusion_input_peak_count
EV1$ComboID <- paste0(EV1$ID, EV1$Cell_class, EV1$Cell_type)
EV1$Fusion <- "EWSR1-ETV1"

EV4 <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1ETV4.txt")
EV4 <- EV4[,-c(2,3)]
colnames(EV4) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EV4$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EV4$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EV4$Fusion_input_peak_count <- as.numeric(str_split_fixed(EV4$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EV4$Control_shared_peak_count <- as.numeric(str_split_fixed(EV4$Overlaps_with_control, pattern = "/", n =2)[,1])
EV4$Frac_fusionpeaks_in_cell <- EV4$Fusion_shared_peak_count/EV4$Fusion_input_peak_count
EV4$ComboID <- paste0(EV4$ID, EV4$Cell_class, EV4$Cell_type)
EV4$Fusion <- "EWSR1-ETV4"

EFEV <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1FEV.txt")
EFEV <- EFEV[,-c(2,3)]
colnames(EFEV) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EFEV$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EFEV$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EFEV$Fusion_input_peak_count <- as.numeric(str_split_fixed(EFEV$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EFEV$Control_shared_peak_count <- as.numeric(str_split_fixed(EFEV$Overlaps_with_control, pattern = "/", n =2)[,1])
EFEV$Frac_fusionpeaks_in_cell <- EFEV$Fusion_shared_peak_count/EFEV$Fusion_input_peak_count
EFEV$ComboID <- paste0(EFEV$ID, EFEV$Cell_class, EFEV$Cell_type)
EFEV$Fusion <- "EWSR1-FEV"

EFLI <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1FLI1.txt")
EFLI <- EFLI[,-c(2,3)]
colnames(EFLI) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EFLI$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EFLI$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EFLI$Fusion_input_peak_count <- as.numeric(str_split_fixed(EFLI$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EFLI$Control_shared_peak_count <- as.numeric(str_split_fixed(EFLI$Overlaps_with_control, pattern = "/", n =2)[,1])
EFLI$Frac_fusionpeaks_in_cell <- EFLI$Fusion_shared_peak_count/EFLI$Fusion_input_peak_count
EFLI$ComboID <- paste0(EFLI$ID, EFLI$Cell_class, EFLI$Cell_type)
EFLI$Fusion <- "EWSR1-FLI1"

ENFAT <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1NFATC2.txt")
ENFAT <- ENFAT[,-c(2,3)]
colnames(ENFAT) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
ENFAT$Fusion_shared_peak_count <- as.numeric(str_split_fixed(ENFAT$Overlaps_with_fusion, pattern = "/", n =2)[,1])
ENFAT$Fusion_input_peak_count <- as.numeric(str_split_fixed(ENFAT$Overlaps_with_fusion, pattern = "/", n =2)[,2])
ENFAT$Control_shared_peak_count <- as.numeric(str_split_fixed(ENFAT$Overlaps_with_control, pattern = "/", n =2)[,1])
ENFAT$Frac_fusionpeaks_in_cell <- ENFAT$Fusion_shared_peak_count/ENFAT$Fusion_input_peak_count
ENFAT$ComboID <- paste0(ENFAT$ID, ENFAT$Cell_class, ENFAT$Cell_type)
ENFAT$Fusion <- "EWSR1-NFATC2"

ENR4 <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1NR4A3.txt")
ENR4 <- ENR4[,-c(2,3)]
colnames(ENR4) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
ENR4$Fusion_shared_peak_count <- as.numeric(str_split_fixed(ENR4$Overlaps_with_fusion, pattern = "/", n =2)[,1])
ENR4$Fusion_input_peak_count <- as.numeric(str_split_fixed(ENR4$Overlaps_with_fusion, pattern = "/", n =2)[,2])
ENR4$Control_shared_peak_count <- as.numeric(str_split_fixed(ENR4$Overlaps_with_control, pattern = "/", n =2)[,1])
ENR4$Frac_fusionpeaks_in_cell <- ENR4$Fusion_shared_peak_count/ENR4$Fusion_input_peak_count
ENR4$ComboID <- paste0(ENR4$ID, ENR4$Cell_class, ENR4$Cell_type)
ENR4$Fusion <- "EWSR1-NR4A3"

EPB <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1PBX1.txt")
EPB <- EPB[,-c(2,3)]
colnames(EPB) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EPB$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EPB$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EPB$Fusion_input_peak_count <- as.numeric(str_split_fixed(EPB$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EPB$Control_shared_peak_count <- as.numeric(str_split_fixed(EPB$Overlaps_with_control, pattern = "/", n =2)[,1])
EPB$Frac_fusionpeaks_in_cell <- EPB$Fusion_shared_peak_count/EPB$Fusion_input_peak_count
EPB$ComboID <- paste0(EPB$ID, EPB$Cell_class, EPB$Cell_type)
EPB$Fusion <- "EWSR1-PBX1"

EPOU <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1POU5F1.txt")
EPOU <- EPOU[,-c(2,3)]
colnames(EPOU) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EPOU$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EPOU$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EPOU$Fusion_input_peak_count <- as.numeric(str_split_fixed(EPOU$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EPOU$Control_shared_peak_count <- as.numeric(str_split_fixed(EPOU$Overlaps_with_control, pattern = "/", n =2)[,1])
EPOU$Frac_fusionpeaks_in_cell <- EPOU$Fusion_shared_peak_count/EPOU$Fusion_input_peak_count
EPOU$ComboID <- paste0(EPOU$ID, EPOU$Cell_class, EPOU$Cell_type)
EPOU$Fusion <- "EWSR1-POU5F1"

ES <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1SP3.txt")
ES <- ES[,-c(2,3)]
colnames(ES) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
ES$Fusion_shared_peak_count <- as.numeric(str_split_fixed(ES$Overlaps_with_fusion, pattern = "/", n =2)[,1])
ES$Fusion_input_peak_count <- as.numeric(str_split_fixed(ES$Overlaps_with_fusion, pattern = "/", n =2)[,2])
ES$Control_shared_peak_count <- as.numeric(str_split_fixed(ES$Overlaps_with_control, pattern = "/", n =2)[,1])
ES$Frac_fusionpeaks_in_cell <- ES$Fusion_shared_peak_count/ES$Fusion_input_peak_count
ES$ComboID <- paste0(ES$ID, ES$Cell_class, ES$Cell_type)
ES$Fusion <- "EWSR1-SP3"

EY <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/EWSR1YY1.txt")
EY <- EY[,-c(2,3)]
colnames(EY) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
EY$Fusion_shared_peak_count <- as.numeric(str_split_fixed(EY$Overlaps_with_fusion, pattern = "/", n =2)[,1])
EY$Fusion_input_peak_count <- as.numeric(str_split_fixed(EY$Overlaps_with_fusion, pattern = "/", n =2)[,2])
EY$Control_shared_peak_count <- as.numeric(str_split_fixed(EY$Overlaps_with_control, pattern = "/", n =2)[,1])
EY$Frac_fusionpeaks_in_cell <- EY$Fusion_shared_peak_count/EY$Fusion_input_peak_count
EY$ComboID <- paste0(EY$ID, EY$Cell_class, EY$Cell_type)
EY$Fusion <- "EWSR1-YY1"

FT <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/FGFR3TACC3.txt")
FT <- FT[,-c(2,3)]
colnames(FT) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
FT$Fusion_shared_peak_count <- as.numeric(str_split_fixed(FT$Overlaps_with_fusion, pattern = "/", n =2)[,1])
FT$Fusion_input_peak_count <- as.numeric(str_split_fixed(FT$Overlaps_with_fusion, pattern = "/", n =2)[,2])
FT$Control_shared_peak_count <- as.numeric(str_split_fixed(FT$Overlaps_with_control, pattern = "/", n =2)[,1])
FT$Frac_fusionpeaks_in_cell <- FT$Fusion_shared_peak_count/FT$Fusion_input_peak_count
FT$ComboID <- paste0(FT$ID, FT$Cell_class, FT$Cell_type)
FT$Fusion <- "FGFR3-TACC3"

FA <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/FUSATF1.txt")
FA <- FA[,-c(2,3)]
colnames(FA) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
FA$Fusion_shared_peak_count <- as.numeric(str_split_fixed(FA$Overlaps_with_fusion, pattern = "/", n =2)[,1])
FA$Fusion_input_peak_count <- as.numeric(str_split_fixed(FA$Overlaps_with_fusion, pattern = "/", n =2)[,2])
FA$Control_shared_peak_count <- as.numeric(str_split_fixed(FA$Overlaps_with_control, pattern = "/", n =2)[,1])
FA$Frac_fusionpeaks_in_cell <- FA$Fusion_shared_peak_count/FA$Fusion_input_peak_count
FA$ComboID <- paste0(FA$ID, FA$Cell_class, FA$Cell_type)
FA$Fusion <- "FUS-ATF1"

FD <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/FUSDDIT3.txt")
FD <- FD[,-c(2,3)]
colnames(FD) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
FD$Fusion_shared_peak_count <- as.numeric(str_split_fixed(FD$Overlaps_with_fusion, pattern = "/", n =2)[,1])
FD$Fusion_input_peak_count <- as.numeric(str_split_fixed(FD$Overlaps_with_fusion, pattern = "/", n =2)[,2])
FD$Control_shared_peak_count <- as.numeric(str_split_fixed(FD$Overlaps_with_control, pattern = "/", n =2)[,1])
FD$Frac_fusionpeaks_in_cell <- FD$Fusion_shared_peak_count/FD$Fusion_input_peak_count
FD$ComboID <- paste0(FD$ID, FD$Cell_class, FD$Cell_type)
FD$Fusion <- "FUS-DDIT3"

FF <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/FUSFEV.txt")
FF <- FF[,-c(2,3)]
colnames(FF) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
FF$Fusion_shared_peak_count <- as.numeric(str_split_fixed(FF$Overlaps_with_fusion, pattern = "/", n =2)[,1])
FF$Fusion_input_peak_count <- as.numeric(str_split_fixed(FF$Overlaps_with_fusion, pattern = "/", n =2)[,2])
FF$Control_shared_peak_count <- as.numeric(str_split_fixed(FF$Overlaps_with_control, pattern = "/", n =2)[,1])
FF$Frac_fusionpeaks_in_cell <- FF$Fusion_shared_peak_count/FF$Fusion_input_peak_count
FF$ComboID <- paste0(FF$ID, FF$Cell_class, FF$Cell_type)
FF$Fusion <- "FUS-FEV"

HN <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/HEY1NCOA2.txt")
HN <- HN[,-c(2,3)]
colnames(HN) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
HN$Fusion_shared_peak_count <- as.numeric(str_split_fixed(HN$Overlaps_with_fusion, pattern = "/", n =2)[,1])
HN$Fusion_input_peak_count <- as.numeric(str_split_fixed(HN$Overlaps_with_fusion, pattern = "/", n =2)[,2])
HN$Control_shared_peak_count <- as.numeric(str_split_fixed(HN$Overlaps_with_control, pattern = "/", n =2)[,1])
HN$Frac_fusionpeaks_in_cell <- HN$Fusion_shared_peak_count/HN$Fusion_input_peak_count
HN$ComboID <- paste0(HN$ID, HN$Cell_class, HN$Cell_type)
HN$Fusion <- "HEY1-NCOA2"

IC <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/IRF2BP2CDX1.txt")
IC <- IC[,-c(2,3)]
colnames(IC) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
IC$Fusion_shared_peak_count <- as.numeric(str_split_fixed(IC$Overlaps_with_fusion, pattern = "/", n =2)[,1])
IC$Fusion_input_peak_count <- as.numeric(str_split_fixed(IC$Overlaps_with_fusion, pattern = "/", n =2)[,2])
IC$Control_shared_peak_count <- as.numeric(str_split_fixed(IC$Overlaps_with_control, pattern = "/", n =2)[,1])
IC$Frac_fusionpeaks_in_cell <- IC$Fusion_shared_peak_count/IC$Fusion_input_peak_count
IC$ComboID <- paste0(IC$ID, IC$Cell_class, IC$Cell_type)
IC$Fusion <- "IRF2BP2-CDX1"

MP <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/MEAF6PHF1.txt")
MP <- MP[,-c(2,3)]
colnames(MP) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
MP$Fusion_shared_peak_count <- as.numeric(str_split_fixed(MP$Overlaps_with_fusion, pattern = "/", n =2)[,1])
MP$Fusion_input_peak_count <- as.numeric(str_split_fixed(MP$Overlaps_with_fusion, pattern = "/", n =2)[,2])
MP$Control_shared_peak_count <- as.numeric(str_split_fixed(MP$Overlaps_with_control, pattern = "/", n =2)[,1])
MP$Frac_fusionpeaks_in_cell <- MP$Fusion_shared_peak_count/MP$Fusion_input_peak_count
MP$ComboID <- paste0(MP$ID, MP$Cell_class, MP$Cell_type)
MP$Fusion <- "MEAF6-PHF1"

PF <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/PAX7FOXO1.txt")
PF <- PF[,-c(2,3)]
colnames(PF) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
PF$Fusion_shared_peak_count <- as.numeric(str_split_fixed(PF$Overlaps_with_fusion, pattern = "/", n =2)[,1])
PF$Fusion_input_peak_count <- as.numeric(str_split_fixed(PF$Overlaps_with_fusion, pattern = "/", n =2)[,2])
PF$Control_shared_peak_count <- as.numeric(str_split_fixed(PF$Overlaps_with_control, pattern = "/", n =2)[,1])
PF$Frac_fusionpeaks_in_cell <- PF$Fusion_shared_peak_count/PF$Fusion_input_peak_count
PF$ComboID <- paste0(PF$ID, PF$Cell_class, PF$Cell_type)
PF$Fusion <- "PAX7-FOXO1"

SS <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/SS18L1SSX1.txt")
SS <- SS[,-c(2,3)]
colnames(SS) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
SS$Fusion_shared_peak_count <- as.numeric(str_split_fixed(SS$Overlaps_with_fusion, pattern = "/", n =2)[,1])
SS$Fusion_input_peak_count <- as.numeric(str_split_fixed(SS$Overlaps_with_fusion, pattern = "/", n =2)[,2])
SS$Control_shared_peak_count <- as.numeric(str_split_fixed(SS$Overlaps_with_control, pattern = "/", n =2)[,1])
SS$Frac_fusionpeaks_in_cell <- SS$Fusion_shared_peak_count/SS$Fusion_input_peak_count
SS$ComboID <- paste0(SS$ID, SS$Cell_class, SS$Cell_type)
SS$Fusion <- "SS18L1-SSX1"

TAFN <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/TAF15NR4A3.txt")
TAFN <- TAFN[,-c(2,3)]
colnames(TAFN) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
TAFN$Fusion_shared_peak_count <- as.numeric(str_split_fixed(TAFN$Overlaps_with_fusion, pattern = "/", n =2)[,1])
TAFN$Fusion_input_peak_count <- as.numeric(str_split_fixed(TAFN$Overlaps_with_fusion, pattern = "/", n =2)[,2])
TAFN$Control_shared_peak_count <- as.numeric(str_split_fixed(TAFN$Overlaps_with_control, pattern = "/", n =2)[,1])
TAFN$Frac_fusionpeaks_in_cell <- TAFN$Fusion_shared_peak_count/TAFN$Fusion_input_peak_count
TAFN$ComboID <- paste0(TAFN$ID, TAFN$Cell_class, TAFN$Cell_type)
TAFN$Fusion <- "TAF15-NR4A3"

TCN <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/TCF12NR4A3.txt")
TCN <- TCN[,-c(2,3)]
colnames(TCN) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
TCN$Fusion_shared_peak_count <- as.numeric(str_split_fixed(TCN$Overlaps_with_fusion, pattern = "/", n =2)[,1])
TCN$Fusion_input_peak_count <- as.numeric(str_split_fixed(TCN$Overlaps_with_fusion, pattern = "/", n =2)[,2])
TCN$Control_shared_peak_count <- as.numeric(str_split_fixed(TCN$Overlaps_with_control, pattern = "/", n =2)[,1])
TCN$Frac_fusionpeaks_in_cell <- TCN$Fusion_shared_peak_count/TCN$Fusion_input_peak_count
TCN$ComboID <- paste0(TCN$ID, TCN$Cell_class, TCN$Cell_type)
TCN$Fusion <- "TCF12-NR4A3"

TFG <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/TFGNR4A3.txt")
TFG <- TFG[,-c(2,3)]
colnames(TFG) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                  "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                  "Log_q_val", "Fold_enrichment")
TFG$Fusion_shared_peak_count <- as.numeric(str_split_fixed(TFG$Overlaps_with_fusion, pattern = "/", n =2)[,1])
TFG$Fusion_input_peak_count <- as.numeric(str_split_fixed(TFG$Overlaps_with_fusion, pattern = "/", n =2)[,2])
TFG$Control_shared_peak_count <- as.numeric(str_split_fixed(TFG$Overlaps_with_control, pattern = "/", n =2)[,1])
TFG$Frac_fusionpeaks_in_cell <- TFG$Fusion_shared_peak_count/TFG$Fusion_input_peak_count
TFG$ComboID <- paste0(TFG$ID, TFG$Cell_class, TFG$Cell_type)
TFG$Fusion <- "TFG-NR4A3"

TE <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/TMPRSS2ERG.txt")
TE <- TE[,-c(2,3)]
colnames(TE) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                   "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                   "Log_q_val", "Fold_enrichment")
TE$Fusion_shared_peak_count <- as.numeric(str_split_fixed(TE$Overlaps_with_fusion, pattern = "/", n =2)[,1])
TE$Fusion_input_peak_count <- as.numeric(str_split_fixed(TE$Overlaps_with_fusion, pattern = "/", n =2)[,2])
TE$Control_shared_peak_count <- as.numeric(str_split_fixed(TE$Overlaps_with_control, pattern = "/", n =2)[,1])
TE$Frac_fusionpeaks_in_cell <- TE$Fusion_shared_peak_count/TE$Fusion_input_peak_count
TE$ComboID <- paste0(TE$ID, TE$Cell_class, TE$Cell_type)
TE$Fusion <- "TMPRSS2-ERG"

TNRK1 <- fread("/Users/mfrenkel/Desktop/FunctionalAnalysis/ChipAtlas_Sig200/TPRNTRK1.txt")
TNRK1 <- TNRK1[,-c(2,3)]
colnames(TNRK1) <- c("ID", "Cell_class", "Cell_type", "Total_peaks_in_cell_type",
                   "Overlaps_with_fusion", "Overlaps_with_control", "Log_p_val",
                   "Log_q_val", "Fold_enrichment")
TNRK1$Fusion_shared_peak_count <- as.numeric(str_split_fixed(TNRK1$Overlaps_with_fusion, pattern = "/", n =2)[,1])
TNRK1$Fusion_input_peak_count <- as.numeric(str_split_fixed(TNRK1$Overlaps_with_fusion, pattern = "/", n =2)[,2])
TNRK1$Control_shared_peak_count <- as.numeric(str_split_fixed(TNRK1$Overlaps_with_control, pattern = "/", n =2)[,1])
TNRK1$Frac_fusionpeaks_in_cell <- TNRK1$Fusion_shared_peak_count/TNRK1$Fusion_input_peak_count
TNRK1$ComboID <- paste0(TNRK1$ID, TNRK1$Cell_class, TNRK1$Cell_type)
TNRK1$Fusion <- "TPR-NTRK1"

combo_chipatlas_sample <- rbind(AG, AH)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, CR)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EP)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EN)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EA)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EC)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, ED)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EV1)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EV4)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EFEV)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EFLI)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, ENFAT)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, ENR4)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EPB)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EPOU)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, ES)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, EY)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, FT)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, FA)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, FD)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, FF)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, HN)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, IC)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, MP)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, PF)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, SS)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, TAFN)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, TCN)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, TFG)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, TE)
combo_chipatlas_sample <- rbind(combo_chipatlas_sample, TNRK1)


combo_filt2 <- combo_chipatlas_sample[-which(combo_chipatlas_sample$Fusion == "ETV6-NTRK3" |
                                               combo_chipatlas_sample$Fusion == "CCDC6-RET" |
                                               combo_chipatlas_sample$Fusion == "FGFR3-TACC3" |
                                               combo_chipatlas_sample$Fusion == "TPR-NTRK1"),]

combo_filt2 <- combo_filt2 %>% filter(Frac_fusionpeaks_in_cell > 0.2 &
                                                Log_p_val < -3)

IDs <- unique(combo_filt2$ComboID)
seqs <- seq(1:length(unique(combo_filt2$ComboID)))
ID_df <- data.frame(IDs, seqs)

combo_filt_merge <- merge(combo_filt2, ID_df, by.x = "ComboID", by.y = "IDs", all.x = T, all.y = T)
combo_filt_merge$ID_label <- paste0(combo_filt_merge$Cell_class, "-", combo_filt_merge$Cell_type, "-",
                                    combo_filt_merge$seqs)

ggplot(combo_filt_merge, aes(x=ID_label, y = Fusion, color = -Log_p_val, size = Frac_fusionpeaks_in_cell)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_color_viridis_c(name = 'p value')

#Can also rearrange rows/columns based on clustering: 

mat <- combo_filt_merge %>% 
  select(ID_label, Fusion, Frac_fusionpeaks_in_cell) %>%  
  pivot_wider(names_from = ID_label, values_from = Frac_fusionpeaks_in_cell) %>% 
  data.frame() 
row.names(mat) <- mat$Fusion 
mat <- mat[,-1]
mat[is.na(mat)] <- 0
clust <- hclust(dist(mat %>% as.matrix()))

mat2 <- combo_filt_merge %>% 
  select(ID_label, Fusion, Frac_fusionpeaks_in_cell) %>%  
  pivot_wider(names_from = Fusion, values_from = Frac_fusionpeaks_in_cell) %>% 
  data.frame() 
row.names(mat2) <- mat2$ID_label  
mat2 <- mat2[,-1] 
mat2[is.na(mat2)] <- 0
clust2 <- hclust(dist(mat2 %>% as.matrix())) 


combo_filt_merge2 <- combo_filt_merge %>% mutate(Fusion = factor(Fusion, levels = clust$labels[clust$order]),
                                                 ID_label = factor(ID_label, levels = clust2$labels[clust2$order]))

ggplot(combo_filt_merge2, aes(x=ID_label, y = Fusion, color = -Log_p_val, size = Frac_fusionpeaks_in_cell)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_viridis_c(name = 'p value')+
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size =6),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.y = element_text(size =16))

ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot
ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()
ggtree_plot_col






