# OncofusionPerturbSeq

This workflow assumes you start with the output of Cellranger-atac count for each sample (i.e. lane of 10X chip) separately (i.e. not aggregated by Cellranger). 

FF_ArchR.R works with Cellranger's tsv files to create ArchR Arrow Files and an ArchR project. Nuclei are quality filtered and genotyping information is used here to identify perturbations in nuclei, a UMAP is created (Fig 1d, 2a), a csv for making the peak heatmap is created (Fig 3), all pairwise comparisons in pseudobulk are made between variants and empty vector (Fig 2b), and the motif discovery (Fig 4c,d)
