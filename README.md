# Oncofusion PROD-ATAC

This workflow assumes you start with the output of Cellranger-atac count for each sample (i.e. lane of 10X chip) separately (i.e. not aggregated by Cellranger). 

BC_Mapping.R: This takes in NGS short read sequencing data off the recombined library genomic DNA and associates hardcoded barcodes with random barcodes such that when random barcodes are later seen in the Spear-ATAC dial-out library they can be assigned to variants in the library.

Dialout_Analysis.R: This takes in the sequence of each dial-out sample individually (i.e. each lane of a 10X chip has its own scATAC library and its own dial-out library) and produces the data frame that says, for each sample, which nuclei (cell) barcode is associated with which random barcode (random barcodes are then associated to variants using the output of BC_Mapping.R). Note this produces several different mapping outputs for each sample with different filter cutoffs. The one used in all of this work is a read count cutoff of 3 and a fraction cutoff of 0.9. 

FF_ArchR.R: works with Cellranger's tsv files to create ArchR Arrow Files and an ArchR project. Nuclei are quality filtered and genotyping information is used here to identify perturbations in nuclei, a UMAP is created (Fig 1d, 2a), a csv for making the peak heatmap is created (Fig 3), all pairwise comparisons in pseudobulk are made between variants and empty vector (Fig 2c) and other relevant comparisons can be made the same way as in Fig 5c, and motif discovery is performed (Fig 4c,d). 

Motif.R: takes the data frames produced earlier and creates Fig 4c,d 

Clust_analysis.R: takes the cluster information generated by ArchR/Seurat and creates Fig 2a,b

GGAA_analysis.R: First uses ArchR to create BED files of marker peaks for variants of interest and then counts the GGAA/TTCC content at those sites

Downsampling.R: Down samples the number of fusion and empty vector nuclei and calculates metrics with the original dataset like fraction of peaks recovered, pearson/spearman correlation of log2FC values etc for Fig 6 b,c.

ChipAtlas.R: Takes the output of ChIP-Atlas Enrichment and creates Fig 4e
