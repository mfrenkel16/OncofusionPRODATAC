#This takes in sequencing data of the dail-out samples (which tells us which random barcode was associated with which nuclei) and 
#tests to see if that random barcode was also seen when directly sequencing the cell library's genomic DNA. If so and there is a match
#Then we known which 10X Genomics cell/nuclei barcode goes whith which random barcode that has a matched variant

library(Biostrings); library(ShortRead); library(stringdist)
library(ggplot2); library(dplyr); library(gridExtra); library(grid)
library(lattice); library(stringi); library(stringr); library(data.table)

FF1 <- readFastq("FFDialout1_S2_R1_001.fastq.gz")
FF2 <- readFastq("FFDialout2_S3_R1_001.fastq.gz")
FF3 <- readFastq("FFDialout3_S4_R1_001.fastq.gz")
FF4 <- readFastq("FFDialout4_S5_R1_001.fastq.gz")
FF5 <- readFastq("FFDialout5_S6_R1_001.fastq.gz")
FF6 <- readFastq("FFDialout6_S7_R1_001.fastq.gz")
FF7 <- readFastq("FFDialout7_S8_R1_001.fastq.gz")
FF8 <- readFastq("FFDialout8_S9_R1_001.fastq.gz")
FF9 <- readFastq("FFDialout9_S10_R1_001.fastq.gz")
FF10 <- readFastq("FFDialout10_S11_R1_001.fastq.gz")
FF11 <- readFastq("FFDialout11_S12_R1_001.fastq.gz")
FF12 <- readFastq("FFDialout12_S13_R1_001.fastq.gz")

ass1 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF1))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF1))$x)

#Note I'm cutting back to 10 bp here to try to preserve more reads
ass1$is_string_present = str_detect(ass1$VBC_string,"ATACTGATTC") &
  str_detect(ass1$VBC_string,"CAAGGTGGTT") 

ass1_filt = ass1[ass1$is_string_present == 1,]

ass1_filt$VBC_formatted = str_split_fixed(ass1_filt$VBC_string,
                                                            pattern="ATACTGATTC",n=2)[,2]

ass1_filt$VBC_formatted = str_split_fixed(ass1_filt$VBC_formatted,
                                                            pattern="CAAGGTGGTT",n=2)[,1]

#Repeat for all of them

ass2 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF2))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF2))$x)
ass2$is_string_present = str_detect(ass2$VBC_string,"ATACTGATTC") &
str_detect(ass2$VBC_string,"CAAGGTGGTT") 
ass2_filt = ass2[ass2$is_string_present == 1,]
ass2_filt$VBC_formatted = str_split_fixed(ass2_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass2_filt$VBC_formatted = str_split_fixed(ass2_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass2_filt)

ass3 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF3))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF3))$x)
ass3$is_string_present = str_detect(ass3$VBC_string,"ATACTGATTC") &
str_detect(ass3$VBC_string,"CAAGGTGGTT") 
ass3_filt = ass3[ass3$is_string_present == 1,]
ass3_filt$VBC_formatted = str_split_fixed(ass3_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass3_filt$VBC_formatted = str_split_fixed(ass3_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass3_filt)

ass4 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF4))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF4))$x)
ass4$is_string_present = str_detect(ass4$VBC_string,"ATACTGATTC") &
str_detect(ass4$VBC_string,"CAAGGTGGTT") 
ass4_filt = ass4[ass4$is_string_present == 1,]
ass4_filt$VBC_formatted = str_split_fixed(ass4_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass4_filt$VBC_formatted = str_split_fixed(ass4_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass4_filt)

ass5 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF5))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF5))$x)
ass5$is_string_present = str_detect(ass5$VBC_string,"ATACTGATTC") &
str_detect(ass5$VBC_string,"CAAGGTGGTT") 
ass5_filt = ass5[ass5$is_string_present == 1,]
ass5_filt$VBC_formatted = str_split_fixed(ass5_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass5_filt$VBC_formatted = str_split_fixed(ass5_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass5_filt)

ass6 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF6))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF6))$x)
ass6$is_string_present = str_detect(ass6$VBC_string,"ATACTGATTC") &
str_detect(ass6$VBC_string,"CAAGGTGGTT") 
ass6_filt = ass6[ass6$is_string_present == 1,]
ass6_filt$VBC_formatted = str_split_fixed(ass6_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass6_filt$VBC_formatted = str_split_fixed(ass6_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass6_filt)

ass7 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF7))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF7))$x)
ass7$is_string_present = str_detect(ass7$VBC_string,"ATACTGATTC") &
str_detect(ass7$VBC_string,"CAAGGTGGTT") 
ass7_filt = ass7[ass7$is_string_present == 1,]
ass7_filt$VBC_formatted = str_split_fixed(ass7_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass7_filt$VBC_formatted = str_split_fixed(ass7_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass7_filt)

ass8 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF8))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF8))$x)
ass8$is_string_present = str_detect(ass8$VBC_string,"ATACTGATTC") &
str_detect(ass8$VBC_string,"CAAGGTGGTT") 
ass8_filt = ass8[ass8$is_string_present == 1,]
ass8_filt$VBC_formatted = str_split_fixed(ass8_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass8_filt$VBC_formatted = str_split_fixed(ass8_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass8_filt)

ass9 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF9))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF9))$x)
ass9$is_string_present = str_detect(ass9$VBC_string,"ATACTGATTC") &
str_detect(ass9$VBC_string,"CAAGGTGGTT") 
ass9_filt = ass9[ass9$is_string_present == 1,]
ass9_filt$VBC_formatted = str_split_fixed(ass9_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass9_filt$VBC_formatted = str_split_fixed(ass9_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass9_filt)

ass10 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF10))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF10))$x)
ass10$is_string_present = str_detect(ass10$VBC_string,"ATACTGATTC") &
str_detect(ass10$VBC_string,"CAAGGTGGTT") 
ass10_filt = ass10[ass10$is_string_present == 1,]
ass10_filt$VBC_formatted = str_split_fixed(ass10_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass10_filt$VBC_formatted = str_split_fixed(ass10_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass10_filt)

ass11 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF11))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF11))$x)
ass11$is_string_present = str_detect(ass11$VBC_string,"ATACTGATTC") &
str_detect(ass11$VBC_string,"CAAGGTGGTT") 
ass11_filt = ass11[ass11$is_string_present == 1,]
ass11_filt$VBC_formatted = str_split_fixed(ass11_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass11_filt$VBC_formatted = str_split_fixed(ass11_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass11_filt)

ass12 = data.frame(CBC = substr((str_split_fixed(as.data.frame(ShortRead::id(FF12))$x, pattern = ":r", n=2)[,2]),1,16), 
                             VBC_string = as.data.frame(ShortRead::sread(FF12))$x)
ass12$is_string_present = str_detect(ass12$VBC_string,"ATACTGATTC") &
str_detect(ass12$VBC_string,"CAAGGTGGTT") 
ass12_filt = ass12[ass12$is_string_present == 1,]
ass12_filt$VBC_formatted = str_split_fixed(ass12_filt$VBC_string,
  pattern="ATACTGATTC",n=2)[,2]
ass12_filt$VBC_formatted = str_split_fixed(ass12_filt$VBC_formatted,
   pattern="CAAGGTGGTT",n=2)[,1]
nrow(ass12_filt)
write.csv(ass1_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample1_dialout_associations.csv")
write.csv(ass2_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample2_dialout_associations.csv")
write.csv(ass3_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample3_dialout_associations.csv")
write.csv(ass4_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample4_dialout_associations.csv")
write.csv(ass5_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample5_dialout_associations.csv")
write.csv(ass6_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample6_dialout_associations.csv")
write.csv(ass7_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample7_dialout_associations.csv")
write.csv(ass8_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample8_dialout_associations.csv")
write.csv(ass9_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample9_dialout_associations.csv")
write.csv(ass10_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample10_dialout_associations.csv")
write.csv(ass11_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample11_dialout_associations.csv")
write.csv(ass12_filt, file = "/staging/mfrenkel/10X/FullFusionDialout/sample12_dialout_associations.csv")

df1 <- ass1_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df1_filt <- df1[which(df1 $TotalReadsForCBC >= 3),]
df1_filt2 <- df1_filt[which(df1_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df1_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample1_association_RC3Frac9.csv")
nrow(df1)
nrow(df1_filt)
nrow(df1_filt2)

df2 <- ass2_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df2_filt <- df2[which(df2 $TotalReadsForCBC >= 3),]
df2_filt2 <- df2_filt[which(df2_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df2_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample2_association_RC3Frac9.csv")
nrow(df2)
nrow(df2_filt)
nrow(df2_filt2)

df3 <- ass3_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df3_filt <- df3[which(df3 $TotalReadsForCBC >= 3),]
df3_filt2 <- df3_filt[which(df3_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df3_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample3_association_RC3Frac9.csv")
nrow(df3)
nrow(df3_filt)
nrow(df3_filt2)

df4 <- ass4_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df4_filt <- df4[which(df4 $TotalReadsForCBC >= 3),]
df4_filt2 <- df4_filt[which(df4_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df4_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample4_association_RC3Frac9.csv")
nrow(df4)
nrow(df4_filt)
nrow(df4_filt2)

df5 <- ass5_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df5_filt <- df5[which(df5 $TotalReadsForCBC >= 3),]
df5_filt2 <- df5_filt[which(df5_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df5_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample5_association_RC3Frac9.csv")
nrow(df5)
nrow(df5_filt)
nrow(df5_filt2)

df6 <- ass6_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df6_filt <- df6[which(df6 $TotalReadsForCBC >= 3),]
df6_filt2 <- df6_filt[which(df6_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df6_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample6_association_RC3Frac9.csv")
nrow(df6)
nrow(df6_filt)
nrow(df6_filt2)

df7 <- ass7_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df7_filt <- df7[which(df7 $TotalReadsForCBC >= 3),]
df7_filt2 <- df7_filt[which(df7_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df7_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample7_association_RC3Frac9.csv")
nrow(df7)
nrow(df7_filt)
nrow(df7_filt2)

df8 <- ass8_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df8_filt <- df8[which(df8 $TotalReadsForCBC >= 3),]
df8_filt2 <- df8_filt[which(df8_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df8_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample8_association_RC3Frac9.csv")
nrow(df8)
nrow(df8_filt)
nrow(df8_filt2)
rm(df8)
rm(df8_filt)
rm(df8_filt2)

df9 <- ass9_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df9_filt <- df9[which(df9 $TotalReadsForCBC >= 3),]
df9_filt2 <- df9_filt[which(df9_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df9_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample9_association_RC3Frac9.csv")
nrow(df9)
nrow(df9_filt)
nrow(df9_filt2)

df10 <- ass10_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df10_filt <- df10[which(df10 $TotalReadsForCBC >= 3),]
df10_filt2 <- df10_filt[which(df10_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df10_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample10_association_RC3Frac9.csv")
nrow(df10)
nrow(df10_filt)
nrow(df10_filt2)

df11 <- ass11_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df11_filt <- df11[which(df11 $TotalReadsForCBC >= 3),]
df11_filt2 <- df11_filt[which(df11_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df11_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample11_association_RC3Frac9.csv")
nrow(df11)
nrow(df11_filt)
nrow(df11_filt2)

df12 <- ass12_filt %>% group_by(CBC, VBC_formatted) %>% summarize(n=n()) %>% group_by(CBC) %>% mutate(TotalReadsForCBC = sum(n), NumberRandomBCs = length(unique(VBC_formatted)), FractionReadsForGivenRandomBC = n/TotalReadsForCBC)
df12_filt <- df12[which(df12 $TotalReadsForCBC >= 3),]
df12_filt2 <- df12_filt[which(df12_filt $FractionReadsForGivenRandomBC >= 0.9),]
write.csv(df12_filt2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample12_association_RC3Frac9.csv")
nrow(df12)
nrow(df12_filt)
nrow(df12_filt2)

#Now associate with known random barcodes that were directly sequenced off of gDNA before (Full_Fusion_Barcode_Mapping_gDNA.csv)

gDNA_map <- read.csv("/staging/mfrenkel/10X/FullFusionDialout/Full_Fusion_Barcode_Mapping_gDNA.csv")
gDNA_map <- gDNA_map[,-c(1,2,4,5)]
gDNA_map$Random_BC_RevComp <- sapply(gDNA_map$Random_BC, function(x) as.character(reverseComplement(DNAString(x))))
gDNA_map <- gDNA_map[,-1]

#Do a sanity check and then merge
mean(df1_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df1 <- merge(df1_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df1$Fusion))
write.csv(merged_df1, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample1_merged_associations.csv")

mean(df2_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df2 <- merge(df2_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df2$Fusion))
write.csv(merged_df2, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample2_merged_associations.csv")

mean(df3_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df3 <- merge(df3_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df3$Fusion))
write.csv(merged_df3, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample3_merged_associations.csv")

mean(df4_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df4 <- merge(df4_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df4$Fusion))
write.csv(merged_df4, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample4_merged_associations.csv")

mean(df5_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df5 <- merge(df5_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df5$Fusion))
write.csv(merged_df5, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample5_merged_associations.csv")

mean(df6_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df6 <- merge(df6_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df6$Fusion))
write.csv(merged_df6, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample6_merged_associations.csv")

mean(df7_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df7 <- merge(df7_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df7$Fusion))
write.csv(merged_df7, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample7_merged_associations.csv")

mean(df8_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df8 <- merge(df8_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df8$Fusion))
write.csv(merged_df8, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample8_merged_associations.csv")

mean(df9_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df9<- merge(df9_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df9$Fusion))
write.csv(merged_df9, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample9_merged_associations.csv")

mean(df10_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df10 <- merge(df10_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df10$Fusion))
write.csv(merged_df10, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample10_merged_associations.csv")

mean(df11_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df11 <- merge(df11_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df11$Fusion))
write.csv(merged_df11, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample11_merged_associations.csv")

mean(df12_filt2$VBC_formatted %in% gDNA_map$Random_BC_RevComp)
merged_df12 <- merge(df12_filt2, gDNA_map, by.x ="VBC_formatted", by.y = "Random_BC_RevComp")
length(unique(merged_df12$Fusion))
write.csv(merged_df12, file = "/staging/mfrenkel/10X/FullFusionDialout/Sample12_merged_associations.csv")

#These 12 files are the final association lists

