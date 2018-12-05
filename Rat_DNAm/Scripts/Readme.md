
# Scripts for the rat DNA methylation manuscript

#### • Preprocessing pipeline
    1. Filter bam files using samtools 
        • samtools_filter_newdata.py
        • Removes duplicate reads, unpaired reads, and reads with a minimum quality score below 10.

    2. Merge bam files for each sample 
        • samtools_merge_newdata.py
        • Merges the data from both sequencing runs into one for each sample 
            • 60 samples were split into 3 pools of 20 samples, each run on two lanes (6 lanes total)
    
    3. Call MeDIP peaks using MACS2 for each sample
        • https://github.com/taoliu/MACS
        • macs2 callpeaks –f BAMPE –m 5 50 –bw 300 –g 2.9e9 –q 0.05
        • Split into three scripts to parallelize (one for each sequencing pool)

#### • Analysis - Developmental profile of the hypothalamus
    1. Generate peaksets for analysis in R
        • Peaksets_developmental.Rmd
            • Creates a peakset for all hypothalamus samples across 4 ages (P1,P8,P15,P22)

    2. Normalize and correct data
        • Normalization_developmental.Rmd
            • Data filtering
            • RPKM normalization
            • ComBat correction
    
    3. Analyze data and plot results
        • Analysis_developmental.Rmd
            • Final glm used in manuscript
            • Genomic enrichment analysis
            • TFBS analysis
            • ermineJ table generation
            • Publication figures

#### • Analysis - Tissue-concordant analysis of hypothalamus and WBC
    1. Generate peaksets for analysis in R
        • Peaksets_tissue_concordant.Rmd
            • Creates a peakset for P22 samples (WBC and hypothalamus)

    2. Normalize and correct data
        • Normalization_tissue_concordant.Rmd
            • Data filtering
            • RPKM normalization
            • ComBat correction
    
    3. Analyze data and plot results
        • Analysis_tissue_concordant.Rmd
            • Final glm used in manuscript
            • Genomic enrichment analysis
            • TFBS analysis
            • ermineJ table generation
            • Publication figures

#### • Other
    1. annotation.R
        • Custom annotations the developmental and tissue-concordant datasets (P22). 
        • Genomic features, DNA sequences, and CG content.
    2. cell_type_correction.Rmd
        • Moderate cell type correction using RNA-seq data from Cahoy et al., 2008
        • Used in the developmental analysis

### • Coming soon:
    1. Analysis - Typical development of the hypothalamus
    2. Unused scripts
        • Analysis of repeat sequences
        • read.bedgraph
        • Plotting bedgraph data in R
