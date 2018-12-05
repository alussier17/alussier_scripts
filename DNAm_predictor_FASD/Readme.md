# DNA methylation as a predictor of fetal alcohol spectrum disorder.

## A. Manuscript
  • Citation: Lussier AA, Morin AM, MacIsaac JL, Salmon J, Weinberg J, Reynolds JN, Pavlidis P, Chudley AE, Kobor MS. DNA methylation as a predictor of fetal alcohol spectrum disorder. Clin Epigenetics. 2018 Jan 12;10:5. 
  
  • doi: 10.1186/s13148-018-0439-6.
  
  • GEO accession number: [GSE10904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109042)
  
  • Source: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-018-0439-6
  
  • Predictor website: https://fasdpredictor.shinyapps.io/fasdpredictorapp/
  
  • Supplementary materials

## [B. Scripts](/DNAm_predictor_FASD/Scripts)
   • Preprocessing and normalization
   
   • Replication analysis
   
   • Predictor generation
   
   • Website source code

## C. Summary of project
**Background**: Fetal alcohol spectrum disorder (FASD) is a developmental disorder that manifests through a range of cognitive, adaptive, physiological, and neurobiological deficits resulting from prenatal alcohol exposure. Although the North American prevalence is currently estimated at 2-5%, FASD has proven difficult to identify in the absence of the overt physical features characteristic of fetal alcohol syndrome. As interventions may have the greatest impact at an early age, accurate biomarkers are needed to identify children at risk for FASD. Building on our previous work identifying distinct DNA methylation patterns in children and adolescents with FASD, we have attempted to validate these associations in a different clinical cohort and to use our DNA methylation signature to develop a possible epigenetic predictor of FASD.

**Methods**: Genome-wide DNA methylation patterns were analyzed using the Illumina HumanMethylation450 array in the buccal epithelial cells of a cohort of 48 individuals aged 3.5-18 (24 FASD cases, 24 controls). The DNA methylation predictor of FASD was built using a stochastic gradient boosting model on our previously published dataset FASD cases and controls (GSE80261). The predictor was tested on the current dataset and an independent dataset of 48 autism spectrum disorder cases and 48 controls (GSE50759).

**Results**: We validated findings from our previous study that identified a DNA methylation signature of FASD, replicating the altered DNA methylation levels of 161/648 CpGs in this independent cohort, which may represent a robust signature of FASD in the epigenome. We also generated a predictive model of FASD using machine learning in a subset of our previously published cohort of 179 samples (83 FASD cases, 96 controls), which was tested in this novel cohort of 48 samples and resulted in a moderately accurate predictor of FASD status. Upon testing the algorithm in an independent cohort of individuals with autism spectrum disorder, we did not detect any bias towards autism, sex, age, or ethnicity.

**Conclusions**: These findings further support the association of FASD with distinct DNA methylation patterns, while providing a possible entry point towards the development of epigenetic biomarkers of FASD.
