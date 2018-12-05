# Scripts for the developmental profile analysis

### 1. Generate peaksets for analysis in R
    • Peaksets_developmental.Rmd
        • Creates a peakset for all hypothalamus samples across 4 ages (P1,P8,P15,P22)

### 2. Normalize and correct data
    • Normalization_developmental.Rmd
        • Data filtering
        • RPKM normalization
        • ComBat correction
        • Please note that this script might be missing some library() commands (but hopefully not!)
    
### 3. Analyze data and plot results
    • Analysis_developmental.Rmd
        • Final glm used in manuscript
        • Genomic enrichment analysis
        • TFBS analysis
        • ermineJ table generation
        • Publication figures
        • Please note that this script might be missing some library() commands (but hopefully not!)
