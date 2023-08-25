#This takes in sequencing data of the barcodes off of genomic DNA
#These are first filtered in bash

gzcat Sample_S1_L001_R1_001.fastq.gz | grep "AAATCCAAGC" > FullFusion_Filtered_Amplicons2
grep "CCAGAGCATG" FullFusion_Filtered_Amplicons > FullFusion_Filtered_Amplicons2
grep "CAAGGTGGTT" FullFusion_Filtered_Amplicons2 > FullFusion_Filtered_Amplicons3
grep "ATACTGATTC" FullFusion_Filtered_Amplicons3 > FullFusion_Filtered_Amplicons4_gDNA

#Subsequently, they are analyzed in R

#The known barcode list is known from the synthetic oligo library ordered from Twist
Known_BCs <- read.csv("/Volumes/Max/Full_Fusion_Library/Barcode_Mapping/FullFusion_TwistBC_List.csv", header = TRUE)
Known_BCs$Rev_Comp <- sapply(Known_BCs$Twist_BC, function(x) as.character(reverseComplement(DNAString(x))))

#These are the reads that were filtered above
Seqs <- fread("/Volumes/Max/Full_Fusion_Library/gDNA_Distribution/FullFusion_Filtered_Amplicons4_gDNA")
colnames(Seqs) <- "Full_Seq"
Seqs$revComp <- sapply(Seqs$Full_Seq, function(x) as.character(reverseComplement(DNAString(x))))
Seqs$Twist_BC <- str_split_fixed(str_split_fixed(Seqs$revComp, pattern = "GCTTGGATTT",  n= 2)[,2], pattern = "CATGCTCTGG", n = 2)[,1]
Seqs$Random_BC <- str_split_fixed(str_split_fixed(Seqs$revComp, pattern = "AACCACCTTG",  n = 2)[,2], pattern = "GAATCAGTAT", n = 2)[,1]

Twist_BC_df <- Seqs %>% dplyr::count(Twist_BC, name = "Twist_BC_Count")
Random_BC_df <- Seqs %>% dplyr::count(Random_BC, name = "Random_BC_Count")

Seqs$Matched_BC <- sapply(Seqs$Twist_BC, function(x) Known_BCs$Twist_BC[amatch(x,Known_BCs$Twist_BC,maxDist=3)])
Seqs$Match_found = "yes"; Seqs$Match_found[is.na(Seqs$Matched_BC)] = "no"
Seqs$Matched_BC[is.na(Seqs$Matched_BC)] = Seqs$Twist_BC[is.na(Seqs$Matched_BC)]

Merged_BCs <- merge(Seqs, Known_BCs, by.x = "Matched_BC", by.y = "Twist_BC")
Merged_BCs2 <- Merged_BCs[,-c(2,3,8)]
Merged_BCs2 <- Merged_BCs2[,-c(2,4)]

test_sum <- Merged_BCs2 %>% group_by(Matched_BC) %>% summarize(reads = n())
test_sum2 <- merge(test_sum, Known_BCs, by.x = "Matched_BC", by.y = "Twist_BC")

Summarized_BCs <- Merged_BCs2 %>% group_by(Matched_BC,Random_BC) %>% summarize(number_pairs = n()) %>% 
  group_by(Random_BC) %>% mutate(total_number_RandomBC = sum(number_pairs)) 


Summarized_BCs2 <- Summarized_BCs[which(Summarized_BCs$number_pairs >1),]

#Two different fractional cutoffs were tested (90% and 95% but all subsequent analysis in the paper was with 95%)
Summarized_BCs3 <- Summarized_BCs2[-(which(Summarized_BCs2$number_pairs/Summarized_BCs2$total_number_RandomBC <0.9)),]
ordered <- Summarized_BCs3[order(-Summarized_BCs3$total_number_RandomBC),]

Summarized_BCs4 <- Summarized_BCs2[-(which(Summarized_BCs2$number_pairs/Summarized_BCs2$total_number_RandomBC <0.95)),]

Final_BC_Mapping <- merge(Summarized_BCs3, Known_BCs, by.x = "Matched_BC", by.y = "Twist_BC" )
Final_BC_Mapping <- Final_BC_Mapping[,-6]

Final_BC_Mapping2 <- merge(Summarized_BCs4, Known_BCs, by.x = "Matched_BC", by.y = "Twist_BC" )
Final_BC_Mapping2 <- Final_BC_Mapping2[,-6]

#The second of these was the final list of associations
write.csv(Final_BC_Mapping, file = "FullFusion_Barcode_Mapping_gDNA.csv")
write.csv(Final_BC_Mapping2, file = "FullFusion_Barcode_Mapping2_gDNA.csv")





