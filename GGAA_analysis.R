#I'm assuming you're starting with an ArchR project here called "proj" and I am doing this only for 
#two example variants below but this can be extended to any variants of interest obviously

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Fusion",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

EWSR1_FLI1_df <- markerList$'EWSR1-FLI1'
EWSR1_FLI1_df <- EWSR1_FLI1_df[,c(1,3,4)]
EWSR1_FLI1_df <- EWSR1_FLI1_df[order(EWSR1_FLI1_df$start),]
EWSR1_FLI1_df$seqnames = factor(EWSR1_FLI1_df$seqnames, levels= c(paste0("chr",seq(1,22)),"chrX"))
EWSR1_FLI1_df = EWSR1_FLI1_df[order(EWSR1_FLI1_df$seqnames),] 
write.table(EWSR1_FLI1_df, col.names = F, row.names =F,
            quote = F, "EWSR1_FLI1_markers.bed", sep = "\t")
write.table(EWSR1_FLI1_df, col.names = F, row.names =F,
            quote = F, "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/GGAA/EWSR1_FLI1_markers.bed", sep = "\t")

EWSR1_ATF1_df <- markerList$'EWSR1-ATF1'
EWSR1_ATF1_df <- EWSR1_ATF1_df[,c(1,3,4)]
EWSR1_ATF1_df <- EWSR1_ATF1_df[order(EWSR1_ATF1_df$start),]
EWSR1_ATF1_df$seqnames = factor(EWSR1_ATF1_df$seqnames, levels= c(paste0("chr",seq(1,22)),"chrX"))
EWSR1_ATF1_df = EWSR1_ATF1_df[order(EWSR1_ATF1_df$seqnames),] 
write.table(EWSR1_ATF1_df, col.names = F, row.names =F,
            quote = F, "EWSR1_ATF1_markers.bed", sep = "\t")
write.table(EWSR1_ATF1_df, col.names = F, row.names =F,
            quote = F, "/staging/mfrenkel/10X/FullFusionATACLib/MyCellRanger/GGAA/EWSR1_ATF1_markers.bed", sep = "\t")

#Get fasta files from those BED files in bash separately. Then back in R:
EWSR1_FLI1_Seqs <- readFasta("OrthEtAlEWSR1FLI1_ChIP_Sequences.fa")
EWSR1_ATF1_Seqs <- readFasta("Moller_common_ChIPEWSR1ATF1_Sequences.fa")

EWSR1_FLI1_df <- data.frame(as.character(id(EWSR1_FLI1_Seqs)), as.character(sread(EWSR1_FLI1_Seqs)))
colnames(EWSR1_FLI1_df) <- c("id", "seq")
EWSR1_ATF1_df <- data.frame(as.character(id(EWSR1_ATF1_Seqs)), as.character(sread(EWSR1_ATF1_Seqs)))
colnames(EWSR1_ATF1_df) <- c("id", "seq")

for(i in 1:nrow(EWSR1FLI1)){
  EWSR1FLI1$GGAA_ct[i] <- str_count(EWSR1FLI1$seq[i], pattern = "GGAA")
  EWSR1FLI1$TTCC_ct[i] <- str_count(EWSR1FLI1$seq[i], pattern = "TTCC")
  EWSR1FLI1$both_ct[i] <- EWSR1FLI1$GGAA_ct[i] + EWSR1FLI1$TTCC_ct[i]
}

EWSR1FLI1$max_GGAA <- 0
for(i in 1:nrow(EWSR1FLI1)){
  for(j in 1:100){
    search_str <- paste0(rep("GGAA",j), collapse = "")
    if(str_detect(EWSR1FLI1$seq[i], pattern = search_str)){
      EWSR1FLI1$max_GGAA[i] <- j
    }
  }
}

EWSR1FLI1$max_TTCC <- 0
for(i in 1:nrow(EWSR1FLI1)){
  for(j in 1:100){
    search_str <- paste0(rep("TTCC",j), collapse = "")
    if(str_detect(EWSR1FLI1$seq[i], pattern = search_str)){
      EWSR1FLI1$max_TTCC[i] <- j
    }
  }
}

for(i in 1:nrow(EWSR1FLI1)){
  EWSR1FLI1$max_combo[i] <- max(EWSR1FLI1$max_GGAA[i], EWSR1FLI1$max_TTCC[i])
}

summary(EWSR1FLI1$max_combo)

for(i in 1:nrow(EWSR1FLI1)){
  EWSR1FLI1$sample[i] <- "EWSR1-FLI1"
}

for(i in 1:nrow(EWSR1ATF1)){
  EWSR1ATF1$GGAA_ct[i] <- str_count(EWSR1ATF1$seq[i], pattern = "GGAA")
  EWSR1ATF1$TTCC_ct[i] <- str_count(EWSR1ATF1$seq[i], pattern = "TTCC")
  EWSR1ATF1$both_ct[i] <- EWSR1ATF1$GGAA_ct[i] + EWSR1ATF1$TTCC_ct[i]
}

EWSR1ATF1$max_GGAA <- 0
for(i in 1:nrow(EWSR1ATF1)){
  for(j in 1:100){
    search_str <- paste0(rep("GGAA",j), collapse = "")
    if(str_detect(EWSR1ATF1$seq[i], pattern = search_str)){
      EWSR1ATF1$max_GGAA[i] <- j
    }
  }
}

EWSR1ATF1$max_TTCC <- 0
for(i in 1:nrow(EWSR1ATF1)){
  for(j in 1:100){
    search_str <- paste0(rep("TTCC",j), collapse = "")
    if(str_detect(EWSR1ATF1$seq[i], pattern = search_str)){
      EWSR1ATF1$max_TTCC[i] <- j
    }
  }
}

for(i in 1:nrow(EWSR1ATF1)){
  EWSR1ATF1$max_combo[i] <- max(EWSR1ATF1$max_GGAA[i], EWSR1ATF1$max_TTCC[i])
}

summary(EWSR1ATF1$max_combo)

for(i in 1:nrow(EWSR1ATF1)){
  EWSR1ATF1$sample[i] <- "EWSR1-ATF1"
}

markerscombo <- rbind(EWSR1FLI1,EWSR1ATF1)
ggplot(markerscombo, aes(x = sample, y = max_combo)) + geom_violin() 
ggplot(markerscombo, aes(x = sample, y = both_ct)) + geom_violin() 

#Obviously, this can be manipulated in any way to count any sort of summary statistic that you think is relevant for GGAA repeat content (i.e. max length,
#mean length, etc etc) which I don't show here


