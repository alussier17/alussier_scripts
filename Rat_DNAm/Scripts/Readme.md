
# Scripts

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

#### • Analysis - Tissue-concordant analysis of hypothalamus and WBC

#### • Analysis - Typical development of the hypothalamus

#### • annotation.R
    • Custom annotations the developmental and tissue-concordant datasets (P22). 
    • Genomic features, DNA sequences, and CG content.

#### • Other
