# Scripts for the tissue-concordant analysis

### 1. Generate peaksets for analysis in R
    • Peaksets_tissue_concordant.Rmd
        • Creates a peakset for P22 samples (WBC and hypothalamus)

### 2. Normalize and correct data
    • Normalization_tissue_concordant.Rmd
        • Data filtering
        • RPKM normalization
        • ComBat correction
        • Please note that this script might be missing some library() commands (but hopefully not!)
    
### 3. Analyze data and plot results
    • Analysis_tissue_concordant.Rmd
        • Final glm used in manuscript
        • Genomic enrichment analysis
        • TFBS analysis
        • ermineJ table generation
        • Publication figures
        • Please note that this script might be missing some library() commands (but hopefully not!)
