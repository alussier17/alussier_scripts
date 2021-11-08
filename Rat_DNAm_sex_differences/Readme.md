# Prenatal Adversity Alters the Epigenetic Profile of the Pre-frontal Cortex: Sexually Dimorphic Effects of Prenatal Alcohol Exposure and Food-related Stress  

## A. Manuscript
  • Citation: 
  
  • doi: 
  
  • GEO accession number: [GSE186840](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186840)
  
  • Source: 

  • [Supplementary materials](/Rat_DNAm_sex_differences/Supplementary_materials)

## B. [Scripts](/Rat_DNAm/Scripts)
  • 1.Preprocessing
        - 1.filtering.sh:       Filtering bam files
        - 2.peaks.sh:           Calling meDIP-seq peaks
        - 3.peakset.diff.R:     Generating the "diff" peakset for downstream analyses
        - 4.peakset.processing: Processing the peakset file for input into edgeR and downstream analyses
    
  • 2.Analysis  
        - 1.diff.linear.models.Rmd: Script for edgeR analyses of sex-concordant and sex-specific differentially methylated regions
        - 2.diff.analysis.Rmd:      Analyses of results from linear models (main figures, overlap with [initial study](/Rat_DNAm/), gene ontology)
        - 3. ASD.overlap.Rmd:       Checking overlap with results of studies of autism spectrum disorder, summarized in ASD_genes_summary_2021-08-12.xlsx
  
  • 3.Annotation
        - 1.diff.annotation.Rmd:   Script to generate annotation specific to the "diff" peakset, shown in diff.peaks.info.bed
        - Bed files were downloaded from the UCSC genome browser using the settings shown in annotation.info.png
           
## C. Project summary
**Background**: Prenatal adversity or stress can have long-term consequences on developmental trajectories and health outcomes. Although the biological mechanisms underlying these effects are poorly under-stood, epigenetic modifications, such as DNA methylation, have the potential to link early-life en-vironments to alterations in physiological systems, with long-term functional implications. We investigated the consequences of two prenatal insults, prenatal alcohol exposure (PAE) and food-related stress, on DNA methylation profiles of the rat brain during early development, with a focus on sex-specific and sex-concordant effects. 

**Methods**: We used methylated DNA immunoprecipitated and next-generation sequencing (meDIP-seq)to analyze epigenome-wide DNA methylation patterns in prefrontal cortex, a key brain region involved in cognition, execu-tive function, and behavior, of both males and females. 

**Results**: We found sex-dependent and sex-concordant influences of these insults on epigenetic patterns. These alterations occurred in genes and pathways related to brain development and immune function, suggesting that PAE and food-related stress may reprogram neurobiological/physiological systems partly through central epigenetic changes, and may do so in a sex-dependent manner. We also found several overlapping genes with autism spectrum disorder, highlighting potential shared etiologies between fetal alcohol spectrum disorder and autism.

**Conclusions**: These epigenetic changes may reflect the sex-specific effects of prenatal insults on long-term functional and health outcomes and have important implications for understanding possible mechanisms underlying fetal alcohol spectrum disorder and other neurodevelopmental disorders.
